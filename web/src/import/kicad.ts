/**
 * KiCad board import (build step 6): parses a .kicad_pcb (KiCad 6–8 S-expression
 * format) into a solver-neutral board model.
 *
 * Coordinates: KiCad uses millimetres with y pointing down. The model uses
 * metres with y pointing up (y_model = -y_kicad), z up from the bottom copper.
 */

export type SExpr = string | SExpr[];

/** Minimal S-expression reader: lists, bare atoms and "quoted strings". */
export function parseSExpr(text: string): SExpr {
  let i = 0;
  const n = text.length;
  const skip = () => { while (i < n && /\s/.test(text[i])) i++; };
  const read = (): SExpr => {
    skip();
    if (text[i] === '(') {
      i++;
      const out: SExpr[] = [];
      for (;;) {
        skip();
        if (i >= n) throw new Error('unbalanced parentheses');
        if (text[i] === ')') { i++; return out; }
        out.push(read());
      }
    }
    if (text[i] === '"') {
      i++;
      let s = '';
      while (i < n && text[i] !== '"') {
        if (text[i] === '\\' && i + 1 < n) { s += text[i + 1]; i += 2; } else s += text[i++];
      }
      i++;
      return s;
    }
    const start = i;
    while (i < n && !/[\s()]/.test(text[i])) i++;
    return text.slice(start, i);
  };
  return read();
}

const isList = (e: SExpr): e is SExpr[] => Array.isArray(e);
const head = (e: SExpr) => (isList(e) && typeof e[0] === 'string' ? e[0] : '');
const kids = (e: SExpr, name: string) => (isList(e) ? e.filter((c) => head(c) === name) : []) as SExpr[][];
const kid = (e: SExpr, name: string) => kids(e, name)[0];
const num = (e: SExpr | undefined, k = 1) => (e && isList(e) ? Number(e[k]) : NaN);
const str = (e: SExpr | undefined, k = 1) => (e && isList(e) ? String(e[k]) : '');

export interface StackLayer {
  name: string;
  type: 'copper' | 'dielectric' | 'other';
  thickness: number;       // m
  epsR?: number;
  lossTangent?: number;
  material?: string;
  roughness?: number;      // copper: rms roughness Rq (m), Groisse model
  mask?: boolean;          // top solder mask: a dielectric slab over the top copper and its gaps
  z0: number; z1: number;  // filled in bottom-up
}

export type Pt = [number, number];
export interface Track { a: Pt; b: Pt; width: number; layer: string; net: number }
export interface Via { at: Pt; drill: number; size: number; layers: [string, string]; net: number }
export interface CopperPoly { pts: Pt[]; layer: string; net: number }
export interface Board {
  layers: StackLayer[];               // bottom to top
  copperLayers: string[];             // names, bottom to top
  tracks: Track[];
  vias: Via[];
  polys: CopperPoly[];                // zones (filled) and pads, as polygons
  outline: { min: Pt; max: Pt };      // bounding box of Edge.Cuts
  nets: Map<number, string>;
  thickness: number;
}

const toM = (v: number) => v * 1e-3;
const P = (x: number, y: number): Pt => [toM(x), -toM(y)];

function rotate(p: Pt, deg: number): Pt {
  const r = (deg * Math.PI) / 180, c = Math.cos(r), s = Math.sin(r);
  return [p[0] * c - p[1] * s, p[0] * s + p[1] * c];
}

export function importKicad(text: string): Board {
  const root = parseSExpr(text);
  if (head(root) !== 'kicad_pcb') throw new Error('not a kicad_pcb file');
  const nets = new Map<number, string>();
  for (const n of kids(root, 'net')) nets.set(Number(n[1]), String(n[2] ?? ''));

  // ---- stackup (setup (stackup (layer ...))) ----
  const setup = kid(root, 'setup');
  const stack = setup ? kid(setup, 'stackup') : undefined;
  const layersTopDown: StackLayer[] = [];
  if (stack) {
    for (const l of kids(stack, 'layer')) {
      const name = String(l[1]);
      const type = str(kid(l, 'type'));
      const t = toM(num(kid(l, 'thickness')));
      if (type === 'copper') layersTopDown.push({ name, type: 'copper', thickness: t || 35e-6, z0: 0, z1: 0 });
      else if (/solder mask/i.test(type) && !layersTopDown.some((q) => q.type === 'copper')) {
        // top solder mask only (the bottom one is below the ground floor)
        layersTopDown.push({ name, type: 'dielectric', mask: true, thickness: t || 20e-6, epsR: num(kid(l, 'epsilon_r')) || 3.5,
          lossTangent: num(kid(l, 'loss_tangent')) || 0.02, z0: 0, z1: 0 });
      }
      else if (/core|prepreg/.test(type)) {
        layersTopDown.push({
          name, type: 'dielectric', thickness: t, epsR: num(kid(l, 'epsilon_r')) || 4.5,
          lossTangent: num(kid(l, 'loss_tangent')) || 0, material: str(kid(l, 'material')), z0: 0, z1: 0,
        });
      }
    }
  }
  if (!layersTopDown.length) {
    // no stackup section: two copper layers on a 1.6 mm FR-4 core
    layersTopDown.push({ name: 'F.Cu', type: 'copper', thickness: 35e-6, z0: 0, z1: 0 },
      { name: 'core', type: 'dielectric', thickness: 1.6e-3 - 70e-6, epsR: 4.5, lossTangent: 0.02, z0: 0, z1: 0 },
      { name: 'B.Cu', type: 'copper', thickness: 35e-6, z0: 0, z1: 0 });
  }
  // bottom-up z: copper is modelled as a zero-thickness sheet at the dielectric
  // surface; its thickness is kept for the thick-copper option
  const layers = [...layersTopDown].reverse();
  let z = 0;
  for (const l of layers) {
    if (l.type === 'dielectric') { l.z0 = z; z += l.thickness; l.z1 = z; } else { l.z0 = z; l.z1 = z; }
  }

  const tracks: Track[] = [];
  for (const s of kids(root, 'segment')) {
    const a = kid(s, 'start'), b = kid(s, 'end');
    tracks.push({ a: P(num(a, 1), num(a, 2)), b: P(num(b, 1), num(b, 2)), width: toM(num(kid(s, 'width'))),
      layer: str(kid(s, 'layer')), net: num(kid(s, 'net')) || 0 });
  }
  // arc tracks (KiCad 7+): start / mid / end, split into short straight segments
  for (const a of kids(root, 'arc')) {
    const p0 = P(num(kid(a, 'start'), 1), num(kid(a, 'start'), 2)), pm = P(num(kid(a, 'mid'), 1), num(kid(a, 'mid'), 2));
    const p1 = P(num(kid(a, 'end'), 1), num(kid(a, 'end'), 2));
    const pts = arcThrough(p0, pm, p1);
    for (let k = 0; k + 1 < pts.length; k++) {
      tracks.push({ a: pts[k], b: pts[k + 1], width: toM(num(kid(a, 'width'))), layer: str(kid(a, 'layer')), net: num(kid(a, 'net')) || 0 });
    }
  }
  const vias: Via[] = [];
  for (const v of kids(root, 'via')) {
    const at = kid(v, 'at'), ly = kid(v, 'layers');
    vias.push({ at: P(num(at, 1), num(at, 2)), size: toM(num(kid(v, 'size'))), drill: toM(num(kid(v, 'drill'))),
      layers: [String(ly?.[1] ?? 'F.Cu'), String(ly?.[2] ?? 'B.Cu')], net: num(kid(v, 'net')) || 0 });
  }
  const polys: CopperPoly[] = [];
  for (const zn of kids(root, 'zone')) {
    const net = num(kid(zn, 'net')) || 0;
    for (const fp of kids(zn, 'filled_polygon')) {
      const pts = kid(fp, 'pts');
      polys.push({ pts: kids(pts, 'xy').map((q) => P(Number(q[1]), Number(q[2]))), layer: str(kid(fp, 'layer')), net });
    }
  }
  // footprints: rectangular / round-rect / circular SMD and through-hole pads
  for (const fpn of kids(root, 'footprint')) {
    const at = kid(fpn, 'at');
    const fx = num(at, 1), fy = num(at, 2), frot = num(at, 3) || 0;
    for (const pad of kids(fpn, 'pad')) {
      const pat = kid(pad, 'at'), sz = kid(pad, 'size');
      const shape = String(pad[3] ?? 'rect');
      const prot = num(pat, 3) || 0;
      const w = num(sz, 1), h = num(sz, 2);
      const net = num(kid(pad, 'net')) || 0;
      // pad centre in board coordinates (KiCad: footprint rotation applies to pad offsets, y down)
      const off = rotate([num(pat, 1), -num(pat, 2)], frot);
      const cx = fx + off[0], cy = fy - off[1];
      const local = padShape(pad, shape, w, h);
      const pts = local.map((q) => { const r = rotate([q[0], -q[1]], prot); return P(cx + r[0], cy - r[1]); });
      const ly = kid(pad, 'layers');
      const names = ly ? ly.slice(1).map(String) : [];
      for (const lname of names) {
        if (lname === '*.Cu') for (const c of layers.filter((l) => l.type === 'copper')) polys.push({ pts, layer: c.name, net });
        else if (lname.endsWith('.Cu')) polys.push({ pts, layer: lname, net });
      }
      // plated through-hole pads get a barrel, like a via
      const drill = kid(pad, 'drill');
      if (String(pad[2]) === 'thru_hole' && drill) {
        const d = drill[1] === 'oval' ? Math.min(Number(drill[2]), Number(drill[3] ?? drill[2])) : Number(drill[1]);
        const cu = layers.filter((l) => l.type === 'copper');
        vias.push({ at: P(cx, cy), drill: toM(d), size: toM(Math.min(w, h)), layers: [cu[cu.length - 1].name, cu[0].name], net });
      }
    }
  }
  // outline: bounding box of Edge.Cuts graphics
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const grow = (p: Pt) => { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); };
  for (const kind of ['gr_line', 'gr_rect']) {
    for (const g of kids(root, kind)) {
      if (str(kid(g, 'layer')) !== 'Edge.Cuts') continue;
      const a = kid(g, 'start'), b = kid(g, 'end');
      grow(P(num(a, 1), num(a, 2))); grow(P(num(b, 1), num(b, 2)));
    }
  }
  for (const g of kids(root, 'gr_arc')) {
    if (str(kid(g, 'layer')) !== 'Edge.Cuts' || !kid(g, 'mid')) continue;
    arcThrough(P(num(kid(g, 'start'), 1), num(kid(g, 'start'), 2)), P(num(kid(g, 'mid'), 1), num(kid(g, 'mid'), 2)),
      P(num(kid(g, 'end'), 1), num(kid(g, 'end'), 2))).forEach(grow);
  }
  for (const g of kids(root, 'gr_circle')) {
    if (str(kid(g, 'layer')) !== 'Edge.Cuts') continue;
    const c = P(num(kid(g, 'center'), 1), num(kid(g, 'center'), 2)), e = P(num(kid(g, 'end'), 1), num(kid(g, 'end'), 2));
    const r = Math.hypot(e[0] - c[0], e[1] - c[1]);
    grow([c[0] - r, c[1] - r]); grow([c[0] + r, c[1] + r]);
  }
  for (const g of kids(root, 'gr_poly')) {
    if (str(kid(g, 'layer')) !== 'Edge.Cuts') continue;
    for (const q of kids(kid(g, 'pts'), 'xy')) grow(P(Number(q[1]), Number(q[2])));
  }
  if (!isFinite(x0)) throw new Error('board outline (Edge.Cuts) not found');
  return {
    layers, copperLayers: layers.filter((l) => l.type === 'copper').map((l) => l.name), tracks, vias, polys,
    outline: { min: [x0, y0], max: [x1, y1] }, nets, thickness: z,
  };
}

/** Points along the circular arc through p0, pm, p1 (p0 and p1 included), at most 10 degrees apart. */
export function arcThrough(p0: Pt, pm: Pt, p1: Pt): Pt[] {
  const [ax, ay] = p0, [bx, by] = pm, [cx, cy] = p1;
  const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  if (Math.abs(d) < 1e-30) return [p0, p1]; // collinear
  const ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
  const uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
  const r = Math.hypot(ax - ux, ay - uy);
  const t0 = Math.atan2(ay - uy, ax - ux), tm = Math.atan2(by - uy, bx - ux), t1 = Math.atan2(cy - uy, cx - ux);
  const norm = (t: number) => { while (t < 0) t += 2 * Math.PI; return t; };
  // sweep from t0 through tm to t1, in whichever direction passes the midpoint
  let sw = norm(t1 - t0);
  if (norm(tm - t0) > sw) sw -= 2 * Math.PI;
  const n = Math.max(1, Math.ceil(Math.abs(sw) / (Math.PI / 18)));
  const out: Pt[] = [p0];
  for (let k = 1; k < n; k++) { const t = t0 + (sw * k) / n; out.push([ux + r * Math.cos(t), uy + r * Math.sin(t)]); }
  out.push(p1);
  return out;
}

/** Pad outline in pad-local KiCad units (mm, y down), before rotation. */
function padShape(pad: SExpr[], shape: string, w: number, h: number): Pt[] {
  const circ = (cx: number, cy: number, r: number, a0: number, a1: number, n: number): Pt[] =>
    Array.from({ length: n + 1 }, (_, k) => [cx + r * Math.cos(a0 + ((a1 - a0) * k) / n), cy + r * Math.sin(a0 + ((a1 - a0) * k) / n)] as Pt);
  const rounded = (r: number): Pt[] => {
    if (r <= 0) return [[-w / 2, -h / 2], [w / 2, -h / 2], [w / 2, h / 2], [-w / 2, h / 2]];
    const ax = w / 2 - r, ay = h / 2 - r, q = Math.PI / 2;
    return [...circ(ax, ay, r, 0, q, 6), ...circ(-ax, ay, r, q, 2 * q, 6), ...circ(-ax, -ay, r, 2 * q, 3 * q, 6), ...circ(ax, -ay, r, 3 * q, 4 * q, 6)];
  };
  switch (shape) {
    case 'circle': return circ(0, 0, w / 2, 0, 2 * Math.PI, 24).slice(0, 24);
    case 'oval': return rounded(Math.min(w, h) / 2); // a stadium
    case 'roundrect': return rounded((num(kid(pad, 'roundrect_rratio')) || 0.25) * Math.min(w, h));
    case 'custom': {
      // anchor shape plus gr_poly primitives (others are ignored)
      const prims = kid(pad, 'primitives');
      const polys = prims ? kids(prims, 'gr_poly').map((g) => kids(kid(g, 'pts'), 'xy').map((q) => [Number(q[1]), Number(q[2])] as Pt)) : [];
      const anchor = rounded(0);
      // a single outline is needed here: the convex hull of anchor and primitives covers the usual custom pads
      return hullPts([...anchor, ...polys.flat()]);
    }
    default: return rounded(0); // rect, trapezoid and chamfered rects as rectangles
  }
}

function hullPts(pts: Pt[]): Pt[] {
  const p = [...pts].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cross = (o: Pt, a: Pt, b: Pt) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo: Pt[] = [], hi: Pt[] = [];
  for (const q of p) { while (lo.length >= 2 && cross(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop(); lo.push(q); }
  for (let i = p.length - 1; i >= 0; i--) { const q = p[i]; while (hi.length >= 2 && cross(hi[hi.length - 2], hi[hi.length - 1], q) <= 0) hi.pop(); hi.push(q); }
  return lo.slice(0, -1).concat(hi.slice(0, -1));
}

// ---- geometry helpers used by the rasteriser and the mesher ----

/** Is point p inside polygon pts (even-odd rule)? */
export function inPoly(p: Pt, pts: Pt[]): boolean {
  let inside = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    const [xi, yi] = pts[i], [xj, yj] = pts[j];
    if ((yi > p[1]) !== (yj > p[1]) && p[0] < ((xj - xi) * (p[1] - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/** Is point p on a track (stadium of half-width w/2 around segment ab)? */
export function onTrack(p: Pt, t: Track): boolean {
  const [ax, ay] = t.a, [bx, by] = t.b;
  const dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy;
  let u = L2 > 0 ? ((p[0] - ax) * dx + (p[1] - ay) * dy) / L2 : 0;
  u = Math.max(0, Math.min(1, u));
  const qx = ax + u * dx - p[0], qy = ay + u * dy - p[1];
  return qx * qx + qy * qy <= (t.width / 2) ** 2 * (1 + 1e-9);
}

/** Copper test on one layer. */
export function copperAt(b: Board, layer: string, p: Pt): boolean {
  for (const t of b.tracks) if (t.layer === layer && onTrack(p, t)) return true;
  for (const q of b.polys) if (q.layer === layer && inPoly(p, q.pts)) return true;
  return false;
}

/**
 * Copper along the row y = const on one layer, as intervals in x, giving exactly
 * the same answers as copperAt() on that row: polygons use the even-odd rule with
 * the same crossing formula (half-open [lo, hi)), tracks the same stadium test
 * (closed). Rasterising row by row costs O(rows x objects) instead of
 * O(samples x objects).
 */
export interface RowCopper { lo: Float64Array; hi: Float64Array; closed: Uint8Array }

export function copperRow(b: Board, layer: string, y: number): RowCopper {
  const lo: number[] = [], hi: number[] = [], closed: number[] = [];
  for (const t of b.tracks) {
    if (t.layer !== layer) continue;
    const h2 = (t.width / 2) ** 2 * (1 + 1e-9);
    const [ax, ay] = t.a, [bx, by] = t.b;
    let l = Infinity, r = -Infinity;
    // end caps
    for (const [cx, cy] of [t.a, t.b]) {
      const d = h2 - (y - cy) ** 2;
      if (d >= 0) { const s = Math.sqrt(d); l = Math.min(l, cx - s); r = Math.max(r, cx + s); }
    }
    // body: 0 <= u <= 1 and |perpendicular distance| <= w/2, both linear in x along the row
    const dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy;
    if (L2 > 0) {
      const L = Math.sqrt(L2), hw = Math.sqrt(h2);
      // u(x) = ((x-ax)dx + (y-ay)dy)/L2 in [0,1];  n(x) = ((x-ax)(-dy) + (y-ay)dx)/L in [-hw, hw]
      let a = -Infinity, c = Infinity;
      const clip = (k: number, m: number, mn: number, mx: number) => { // mn <= k*x + m <= mx
        if (Math.abs(k) < 1e-300) { if (m < mn || m > mx) { a = Infinity; c = -Infinity; } return; }
        const x1 = (mn - m) / k, x2 = (mx - m) / k;
        a = Math.max(a, Math.min(x1, x2)); c = Math.min(c, Math.max(x1, x2));
      };
      clip(dx / L2, (-ax * dx + (y - ay) * dy) / L2, 0, 1);
      clip(-dy / L, (ax * dy + (y - ay) * dx) / L, -hw, hw);
      if (a <= c) { l = Math.min(l, a); r = Math.max(r, c); }
    }
    if (l <= r) { lo.push(l); hi.push(r); closed.push(1); }
  }
  for (const q of b.polys) {
    if (q.layer !== layer) continue;
    const xs: number[] = [];
    const pts = q.pts;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const [xi, yi] = pts[i], [xj, yj] = pts[j];
      if ((yi > y) !== (yj > y)) xs.push(((xj - xi) * (y - yi)) / (yj - yi) + xi);
    }
    xs.sort((m, n) => m - n);
    for (let k = 0; k + 1 < xs.length; k += 2) { lo.push(xs[k]); hi.push(xs[k + 1]); closed.push(0); }
  }
  return { lo: Float64Array.from(lo), hi: Float64Array.from(hi), closed: Uint8Array.from(closed) };
}

export function inRow(r: RowCopper, x: number): boolean {
  for (let k = 0; k < r.lo.length; k++) {
    if (x >= r.lo[k] && (x < r.hi[k] || (r.closed[k] === 1 && x <= r.hi[k]))) return true;
  }
  return false;
}

/** Axis-aligned copper boundaries on a layer, with the metal side, for the mesher. */
export function copperEdges(b: Board, layer: string, axis: 0 | 1): { x: number; metal: 'lo' | 'hi' }[] {
  const out = new Map<string, { x: number; metal: 'lo' | 'hi' }>();
  const add = (x: number, metal: 'lo' | 'hi') => out.set(`${x.toFixed(9)}${metal}`, { x, metal });
  for (const t of b.tracks) {
    if (t.layer !== layer) continue;
    const along = Math.abs(t.b[0] - t.a[0]) > Math.abs(t.b[1] - t.a[1]) ? 0 : 1;
    if (along === axis) continue; // edges across the track's length are its sides only when perpendicular
    const c = axis === 0 ? t.a[0] : t.a[1];
    add(c - t.width / 2, 'hi'); add(c + t.width / 2, 'lo');
  }
  for (const q of b.polys) {
    if (q.layer !== layer) continue;
    for (let i = 0; i < q.pts.length; i++) {
      const p0 = q.pts[i], p1 = q.pts[(i + 1) % q.pts.length];
      const o = 1 - axis; // edge is a line of constant coordinate `axis` if it runs along the other axis
      if (Math.abs(p0[axis] - p1[axis]) > 1e-12 || Math.abs(p0[o] - p1[o]) < 1e-12) continue;
      const x = p0[axis];
      const mid: Pt = [0, 0];
      mid[o] = (p0[o] + p1[o]) / 2;
      mid[axis] = x + 1e-7;
      add(x, inPoly(mid, q.pts) ? 'hi' : 'lo');
    }
  }
  return [...out.values()].sort((a, c) => a.x - c.x);
}
