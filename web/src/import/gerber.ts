/**
 * Gerber (RS-274X, with X2 attributes when present) + Excellon drill import
 * (build step 6, second half). Produces the same solver-neutral Board model as
 * the KiCad importer; the stackup, which Gerbers do not carry, comes from the
 * app's stackup editor (or a two-layer FR-4 default).
 *
 * Supported: %FS/%MO, standard apertures C/R/O/P, aperture macros (primitives
 * 1 circle, 4 outline, 5 polygon, 20 vector line, 21 centre line; arithmetic in
 * parameters), D01/D02/D03, G01/G02/G03 with G74/G75, regions G36/G37, %LP
 * (clear polarity is reported, not applied), %TF.FileFunction, %TO.N net names.
 * Not supported: step-and-repeat (%SR), block apertures (%AB), mirroring and
 * rotation of apertures (%LM/%LR/%LS) — each raises a warning.
 *
 * Coordinates: Gerber y points up, metres in the model. KiCad's Gerber export
 * writes y_gerber = -y_pcbnew, the same flip the KiCad importer applies, so a
 * board imported both ways lands at the same coordinates.
 */
import type { Board, CopperPoly, Pt, StackLayer, Track, Via } from './kicad.ts';
import { inPoly } from './kicad.ts';

// ------------------------------------------------------------------ Gerber

export interface GerberImage {
  tracks: { a: Pt; b: Pt; width: number; net?: string }[];
  polys: { pts: Pt[]; net?: string; flash?: { at: Pt; dia: number } }[];
  fileFunction?: string[];      // e.g. ['Copper', 'L1', 'Top']
  warnings: string[];
}

interface Aperture { kind: 'C' | 'R' | 'O' | 'P' | 'M'; p: number[]; macro?: string }

const CIRCLE_N = 24;

function circlePts(c: Pt, r: number, n = CIRCLE_N): Pt[] {
  return Array.from({ length: n }, (_, k) => [c[0] + r * Math.cos((2 * Math.PI * k) / n), c[1] + r * Math.sin((2 * Math.PI * k) / n)] as Pt);
}

function rot(p: Pt, deg: number): Pt {
  if (!deg) return p;
  const r = (deg * Math.PI) / 180, c = Math.cos(r), s = Math.sin(r);
  return [p[0] * c - p[1] * s, p[0] * s + p[1] * c];
}

/** Convex hull (monotone chain), counter-clockwise. */
function hull(pts: Pt[]): Pt[] {
  const p = [...pts].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cross = (o: Pt, a: Pt, b: Pt) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo: Pt[] = [], hi: Pt[] = [];
  for (const q of p) { while (lo.length >= 2 && cross(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop(); lo.push(q); }
  for (let i = p.length - 1; i >= 0; i--) { const q = p[i]; while (hi.length >= 2 && cross(hi[hi.length - 2], hi[hi.length - 1], q) <= 0) hi.pop(); hi.push(q); }
  return lo.slice(0, -1).concat(hi.slice(0, -1));
}

/** Arithmetic in macro parameters: $n variables, + - x X / and parentheses. */
function evalExpr(src: string, vars: Map<number, number>): number {
  const s = src.replace(/\s+/g, '');
  let i = 0;
  const atom = (): number => {
    if (s[i] === '(') { i++; const v = sum(); i++; return v; }
    if (s[i] === '-') { i++; return -atom(); }
    if (s[i] === '+') { i++; return atom(); }
    if (s[i] === '$') { i++; const m = /^\d+/.exec(s.slice(i))!; i += m[0].length; return vars.get(Number(m[0])) ?? 0; }
    const m = /^(\d+\.?\d*|\.\d+)/.exec(s.slice(i));
    if (!m) throw new Error(`bad macro expression "${src}"`);
    i += m[0].length;
    return Number(m[0]);
  };
  const prod = (): number => {
    let v = atom();
    while (i < s.length && /[xX/]/.test(s[i])) { const op = s[i++]; const w = atom(); v = op === '/' ? v / w : v * w; }
    return v;
  };
  const sum = (): number => {
    let v = prod();
    while (i < s.length && /[+-]/.test(s[i])) { const op = s[i++]; const w = prod(); v = op === '+' ? v + w : v - w; }
    return v;
  };
  return sum();
}

/** Flash an aperture at `at`: list of polygons (model units) and, for round shapes, the diameter. */
function flashShape(ap: Aperture, macros: Map<string, string[]>, at: Pt, scale: number, warn: (m: string) => void): { polys: Pt[][]; dia?: number } {
  const S = (v: number) => v * scale;
  const tr = (q: Pt, deg = 0): Pt => { const r = rot(q, deg); return [at[0] + S(r[0]), at[1] + S(r[1])]; };
  switch (ap.kind) {
    case 'C': return { polys: [circlePts(at, S(ap.p[0]) / 2)], dia: S(ap.p[0]) };
    case 'R': { const [w, h] = ap.p; return { polys: [[tr([-w / 2, -h / 2]), tr([w / 2, -h / 2]), tr([w / 2, h / 2]), tr([-w / 2, h / 2])]] }; }
    case 'O': {
      const [w, h] = ap.p, r = Math.min(w, h) / 2, pts: Pt[] = [];
      const c1: Pt = w > h ? [-(w / 2 - r), 0] : [0, -(h / 2 - r)], c2: Pt = [-c1[0], -c1[1]];
      for (const c of [c1, c2]) for (const q of circlePts(c, r)) pts.push(q);
      return { polys: [hull(pts).map((q) => tr(q))], dia: w === h ? S(w) : undefined };
    }
    case 'P': {
      const [d, n, a = 0] = ap.p;
      return { polys: [Array.from({ length: n }, (_, k) => tr(rot([d / 2, 0], a + (360 * k) / n)))] };
    }
    case 'M': {
      const body = macros.get(ap.macro!);
      if (!body) { warn(`undefined aperture macro ${ap.macro}`); return { polys: [] }; }
      const vars = new Map<number, number>();
      ap.p.forEach((v, k) => vars.set(k + 1, v));
      const polys: Pt[][] = [];
      for (const stmt of body) {
        if (!stmt || stmt.startsWith('0')) continue; // comment
        const asg = /^\$(\d+)=(.*)$/.exec(stmt);
        if (asg) { vars.set(Number(asg[1]), evalExpr(asg[2], vars)); continue; }
        const f = stmt.split(',').map((t) => evalExpr(t, vars));
        const code = f[0], exposure = f[1];
        if (exposure === 0) { warn('aperture macro exposure-off primitive ignored'); continue; }
        if (code === 1) { const [, , d, x, y, a = 0] = f; const c = rot([x, y], a); polys.push(circlePts(tr(c), S(d) / 2)); }
        else if (code === 4) {
          const n = f[2], pts: Pt[] = [];
          for (let k = 0; k < n; k++) pts.push([f[3 + 2 * k], f[4 + 2 * k]]);
          const a = f[3 + 2 * (n + 1)] ?? 0;
          polys.push(pts.map((q) => tr(q, a)));
        } else if (code === 5) { const [, , n, x, y, d, a = 0] = f; polys.push(Array.from({ length: n }, (_, k) => tr(rot([x + (d / 2) * Math.cos((2 * Math.PI * k) / n), y + (d / 2) * Math.sin((2 * Math.PI * k) / n)], a)))); }
        else if (code === 20) {
          const [, , w, x1, y1, x2, y2, a = 0] = f, L = Math.hypot(x2 - x1, y2 - y1) || 1;
          const nx = (-(y2 - y1) / L) * w / 2, ny = ((x2 - x1) / L) * w / 2;
          polys.push(([[x1 + nx, y1 + ny], [x2 + nx, y2 + ny], [x2 - nx, y2 - ny], [x1 - nx, y1 - ny]] as Pt[]).map((q) => tr(q, a)));
        } else if (code === 21) {
          const [, , w, h, x, y, a = 0] = f;
          polys.push(([[x - w / 2, y - h / 2], [x + w / 2, y - h / 2], [x + w / 2, y + h / 2], [x - w / 2, y + h / 2]] as Pt[]).map((q) => tr(q, a)));
        } else warn(`aperture macro primitive ${code} not supported`);
      }
      return { polys };
    }
  }
}

/** Points along an arc from p0 to p1 about centre c (excluding p0), ccw or cw. */
function arcPts(p0: Pt, p1: Pt, c: Pt, ccw: boolean, full: boolean): Pt[] {
  const r0 = Math.hypot(p0[0] - c[0], p0[1] - c[1]), r1 = Math.hypot(p1[0] - c[0], p1[1] - c[1]);
  const a0 = Math.atan2(p0[1] - c[1], p0[0] - c[0]);
  let a1 = Math.atan2(p1[1] - c[1], p1[0] - c[0]);
  let sweep = a1 - a0;
  if (ccw) { while (sweep <= 0) sweep += 2 * Math.PI; } else { while (sweep >= 0) sweep -= 2 * Math.PI; }
  if (!full && Math.abs(Math.abs(sweep) - 2 * Math.PI) < 1e-9) sweep = 0; // coincident endpoints without a full-circle intent
  a1 = a0 + sweep;
  const n = Math.max(1, Math.ceil(Math.abs(sweep) / (Math.PI / 24)));
  const out: Pt[] = [];
  for (let k = 1; k <= n; k++) {
    const t = k / n, a = a0 + sweep * t, r = r0 + (r1 - r0) * t;
    out.push(k === n ? [p1[0], p1[1]] : [c[0] + r * Math.cos(a), c[1] + r * Math.sin(a)]);
  }
  return out;
}

export function parseGerber(text: string): GerberImage {
  const img: GerberImage = { tracks: [], polys: [], warnings: [] };
  const warned = new Set<string>();
  const warn = (m: string) => { if (!warned.has(m)) { warned.add(m); img.warnings.push(m); } };
  let intD = 2, decD = 4, trailing = false, unit = 25.4e-3; // default inches until %MO
  let unitSet = false;
  const aps = new Map<number, Aperture>();
  const macros = new Map<string, string[]>();
  let cur: Aperture | undefined;
  let x = 0, y = 0, interp: 'G01' | 'G02' | 'G03' = 'G01', multi = true, region = false, dark = true, lastD = 2;
  let regionPts: Pt[] = [];
  let net: string | undefined;

  const coord = (s: string): number => {
    if (s.includes('.')) return Number(s);
    const neg = s.startsWith('-');
    let d = s.replace(/^[+-]/, '');
    if (trailing) d = d.padEnd(intD + decD, '0');
    const v = Number(d) / 10 ** decD;
    return neg ? -v : v;
  };
  const M = (v: number) => v * unit;

  const flushRegion = () => {
    if (regionPts.length >= 3) {
      if (dark) img.polys.push({ pts: regionPts, net });
      else warn('clear-polarity (%LPC) objects are not subtracted');
    }
    regionPts = [];
  };

  // split into extended (%...%) blocks and word commands (...*)
  let i = 0;
  const n = text.length;
  while (i < n) {
    const ch = text[i];
    if (/\s/.test(ch)) { i++; continue; }
    if (ch === '%') {
      const end = text.indexOf('%', i + 1);
      if (end < 0) throw new Error('unterminated % block');
      const blk = text.slice(i + 1, end).replace(/[\r\n]/g, '');
      i = end + 1;
      const cmds = blk.split('*').filter((c) => c.length);
      if (!cmds.length) continue;
      const c0 = cmds[0];
      if (c0.startsWith('FS')) {
        const m = /^FS([LTD])?A?X(\d)(\d)Y(\d)(\d)/.exec(c0);
        if (!m) throw new Error(`bad format spec ${c0}`);
        trailing = m[1] === 'T'; intD = Number(m[2]); decD = Number(m[3]);
      } else if (c0.startsWith('MO')) { unit = c0 === 'MOMM' ? 1e-3 : 25.4e-3; unitSet = true; }
      else if (c0.startsWith('AD')) {
        const m = /^ADD(\d+)([^,]+),?(.*)$/.exec(c0);
        if (!m) { warn(`bad aperture ${c0}`); continue; }
        const p = m[3] ? m[3].split('X').map(Number) : [];
        const name = m[2];
        const kind = (['C', 'R', 'O', 'P'].includes(name) ? name : 'M') as Aperture['kind'];
        aps.set(Number(m[1]), { kind, p, macro: kind === 'M' ? name : undefined });
      } else if (c0.startsWith('AM')) {
        macros.set(c0.slice(2), cmds.slice(1).map((c) => c.trim()));
      } else if (c0.startsWith('LP')) { if (region) flushRegion(); dark = c0 === 'LPD'; if (!dark) warn('clear-polarity (%LPC) objects are not subtracted'); }
      else if (c0.startsWith('TF.FileFunction')) img.fileFunction = c0.split(',').slice(1);
      else if (c0.startsWith('TO.N,')) net = c0.slice(5).split(',')[0] || undefined;
      else if (c0 === 'TD' || c0.startsWith('TD.N')) net = undefined;
      else if (/^(SR|AB|LM|LR|LS)/.test(c0) && !/^SR$/.test(c0) && !/^SRX1Y1/.test(c0)) warn(`%${c0.slice(0, 2)} is not supported`);
      continue;
    }
    // word command
    const end = text.indexOf('*', i);
    if (end < 0) break;
    let w = text.slice(i, end).replace(/[\r\n\s]/g, '');
    i = end + 1;
    if (!w) continue;
    if (w.startsWith('G04')) continue;
    if (w === 'M02' || w === 'M00' || w === 'M01') break;
    if (w === 'G36') { region = true; regionPts = []; continue; }
    if (w === 'G37') { flushRegion(); region = false; continue; }
    if (w === 'G74') { multi = false; continue; }
    if (w === 'G75') { multi = true; continue; }
    if (w === 'G70') { unit = 25.4e-3; unitSet = true; continue; }
    if (w === 'G71') { unit = 1e-3; unitSet = true; continue; }
    const g = /^G0?([123])/.exec(w);
    if (g) { interp = (`G0${g[1]}`) as typeof interp; w = w.slice(g[0].length); }
    w = w.replace(/^G54/, '').replace(/^G55/, '');
    const sel = /^D(\d+)$/.exec(w);
    if (sel && Number(sel[1]) >= 10) { cur = aps.get(Number(sel[1])); if (!cur) warn(`undefined aperture D${sel[1]}`); continue; }
    if (!w) continue;
    const m = /^(?:X([+-]?[\d.]+))?(?:Y([+-]?[\d.]+))?(?:I([+-]?[\d.]+))?(?:J([+-]?[\d.]+))?(?:D0?([123]))?$/.exec(w);
    if (!m) { warn(`unrecognised command ${w.slice(0, 20)}`); continue; }
    if (!unitSet) { warn('no %MO units; assuming inches'); unitSet = true; }
    const nx = m[1] !== undefined ? M(coord(m[1])) : x, ny = m[2] !== undefined ? M(coord(m[2])) : y;
    const I = m[3] !== undefined ? M(coord(m[3])) : 0, J = m[4] !== undefined ? M(coord(m[4])) : 0;
    const d = m[5] ? Number(m[5]) : lastD;
    lastD = d;
    const p0: Pt = [x, y], p1: Pt = [nx, ny];
    if (d === 2) { if (region) { flushRegion(); } x = nx; y = ny; continue; }
    if (d === 3) {
      x = nx; y = ny;
      if (!cur) { warn('flash without aperture'); continue; }
      if (!dark) { warn('clear-polarity (%LPC) objects are not subtracted'); continue; }
      const sh = flashShape(cur, macros, p1, unit, warn);
      for (const pts of sh.polys) img.polys.push({ pts, net, flash: { at: p1, dia: sh.dia ?? 0 } });
      continue;
    }
    // D01: draw
    let path: Pt[];
    if (interp === 'G01') path = [p1];
    else {
      let c: Pt;
      if (multi) c = [x + I, y + J];
      else {
        // single-quadrant: pick the centre sign that gives a sweep <= 90 deg
        const cands: Pt[] = [[x + I, y + J], [x - I, y + J], [x + I, y - J], [x - I, y - J]];
        c = cands.reduce((best, q) => {
          const e = Math.abs(Math.hypot(p1[0] - q[0], p1[1] - q[1]) - Math.hypot(p0[0] - q[0], p0[1] - q[1]));
          const eb = Math.abs(Math.hypot(p1[0] - best[0], p1[1] - best[1]) - Math.hypot(p0[0] - best[0], p0[1] - best[1]));
          return e < eb ? q : best;
        });
      }
      const full = multi && Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) < 1e-12;
      path = arcPts(p0, p1, c, interp === 'G03', full);
    }
    if (region) {
      if (!regionPts.length) regionPts.push(p0);
      regionPts.push(...path);
    } else if (!cur) warn('draw without aperture');
    else if (!dark) warn('clear-polarity (%LPC) objects are not subtracted');
    else {
      let a = p0;
      for (const b of path) {
        if (cur.kind === 'C') img.tracks.push({ a, b, width: cur.p[0] * unit, net });
        else {
          const s0 = flashShape(cur, macros, a, unit, warn).polys.flat(), s1 = flashShape(cur, macros, b, unit, warn).polys.flat();
          img.polys.push({ pts: hull([...s0, ...s1]), net });
          if (cur.kind !== 'R') warn('draws with non-circular apertures use the convex hull');
        }
        a = b;
      }
    }
    x = nx; y = ny;
  }
  return img;
}

// ------------------------------------------------------------------ Excellon

export interface DrillHit { at: Pt; dia: number; plated: boolean }

export function parseExcellon(text: string, platedDefault = true): DrillHit[] {
  let unit = 25.4e-3, lz = true, intD = 2, decD = 4, plated = platedDefault;
  const tools = new Map<number, number>();
  let tool = 0, x = 0, y = 0;
  const hits: DrillHit[] = [];
  const coord = (s: string) => {
    if (s.includes('.')) return Number(s) * unit;
    const neg = s.startsWith('-');
    let d = s.replace(/^[+-]/, '');
    // LZ: leading zeros kept, trailing omitted -> pad right; TZ: the reverse
    d = lz ? d.padEnd(intD + decD, '0') : d.padStart(intD + decD, '0');
    const v = (Number(d) / 10 ** decD) * unit;
    return neg ? -v : v;
  };
  for (const raw of text.split(/\r?\n/)) {
    const l = raw.trim();
    if (!l) continue;
    if (l.startsWith(';')) {
      if (/TYPE=NON_PLATED|NPTH/i.test(l)) plated = false;
      else if (/TYPE=PLATED/i.test(l)) plated = true;
      continue;
    }
    if (/^(METRIC|M71)/.test(l)) { unit = 1e-3; intD = 3; decD = 3; }
    else if (/^(INCH|M72)/.test(l)) { unit = 25.4e-3; intD = 2; decD = 4; }
    if (/^(METRIC|INCH)/.test(l)) {
      if (/,TZ/.test(l)) lz = false; else if (/,LZ/.test(l)) lz = true;
      const f = /,(0+)\.(0+)/.exec(l);
      if (f) { intD = f[1].length; decD = f[2].length; }
      continue;
    }
    if (/^%|^M48|^M30|^G90|^G05|^FMAT|^ICI|^VER|^M95|^T0$/.test(l)) continue;
    const td = /^T(\d+)(?:F\d+)?(?:S\d+)?C([\d.]+)/.exec(l);
    if (td) { tools.set(Number(td[1]), Number(td[2]) * unit); continue; }
    const ts = /^T(\d+)$/.exec(l);
    if (ts) { tool = Number(ts[1]); continue; }
    const xy = /^(?:X([+-]?[\d.]+))?(?:Y([+-]?[\d.]+))?/.exec(l);
    if (xy && (xy[1] || xy[2])) {
      if (xy[1]) x = coord(xy[1]);
      if (xy[2]) y = coord(xy[2]);
      if (/G85/.test(l)) continue; // slots are not modelled
      hits.push({ at: [x, y], dia: tools.get(tool) ?? 0, plated });
    }
  }
  return hits;
}

// ------------------------------------------------------------------ board assembly

export type GerberRole = { kind: 'copper'; index: number; side: 'top' | 'inner' | 'bottom' } | { kind: 'profile' } |
  { kind: 'drill'; plated: boolean } | { kind: 'other' };

/** Guess what a file is from its X2 attributes, falling back to KiCad / Altium / Eagle naming. */
export function classify(name: string, text: string): GerberRole {
  const lower = name.toLowerCase();
  if (/\.(drl|xln|exc|drd|nc)$/.test(lower) || /^M48/m.test(text.slice(0, 200))) {
    return { kind: 'drill', plated: !/npth|non.?plated/i.test(lower + text.slice(0, 400)) };
  }
  const ff = /%TF\.FileFunction,([^*]*)\*%/.exec(text);
  if (ff) {
    const f = ff[1].split(',');
    if (f[0] === 'Copper') {
      const idx = Number((f[1] ?? 'L1').replace(/^L/, ''));
      const side = f[2] === 'Top' ? 'top' : f[2] === 'Bot' ? 'bottom' : 'inner';
      return { kind: 'copper', index: idx, side };
    }
    if (f[0] === 'Profile') return { kind: 'profile' };
    if (f[0] === 'Plated' || f[0] === 'NonPlated') return { kind: 'drill', plated: f[0] === 'Plated' };
    return { kind: 'other' };
  }
  if (/[-_.]f[._]cu\b|\.gtl$|\.cmp$|top.?copper|\.top$/.test(lower)) return { kind: 'copper', index: 1, side: 'top' };
  if (/[-_.]b[._]cu\b|\.gbl$|\.sol$|bot(tom)?.?copper|\.bot$/.test(lower)) return { kind: 'copper', index: 99, side: 'bottom' };
  const inner = /[-_.]in(\d+)[._]cu\b/.exec(lower) ?? /\.g(\d+)$/.exec(lower) ?? /\.gp?(\d+)$/.exec(lower);
  if (inner) return { kind: 'copper', index: Number(inner[1]) + 1, side: 'inner' };
  if (/edge[._]cuts|\.gko$|\.gm1$|\.gml$|profile|outline/.test(lower)) return { kind: 'profile' };
  return { kind: 'other' };
}

export interface StackupRow { name: string; type: 'copper' | 'dielectric'; thickness: number; epsR?: number; lossTangent?: number; roughness?: number; mask?: boolean }

/** Two-layer FR-4 default, or an evenly split stack for n copper layers (top to bottom). */
export function defaultStackup(nCopper: number, names?: string[]): StackupRow[] {
  const rows: StackupRow[] = [];
  const nd = Math.max(1, nCopper - 1), t = (1.6e-3 - nCopper * 35e-6) / nd;
  for (let k = 0; k < nCopper; k++) {
    rows.push({ name: names?.[k] ?? (k === 0 ? 'Top' : k === nCopper - 1 ? 'Bottom' : `In${k}`), type: 'copper', thickness: 35e-6 });
    if (k < nCopper - 1) rows.push({ name: `diel${k + 1}`, type: 'dielectric', thickness: t, epsR: 4.4, lossTangent: 0.02 });
  }
  return rows;
}

export interface GerberFile { name: string; text: string }

export interface GerberImport {
  board: Board;
  copperFiles: string[];        // top to bottom, matching the stackup's copper rows
  warnings: string[];
}

/**
 * Assemble a Board from Gerber copper layers, an optional profile, and drill
 * files. `stackup` lists rows top to bottom; its copper rows are matched to the
 * copper files in top-to-bottom order (defaultStackup() when omitted).
 */
export function importGerber(files: GerberFile[], stackup?: StackupRow[]): GerberImport {
  const warnings: string[] = [];
  const copper: { file: GerberFile; role: Extract<GerberRole, { kind: 'copper' }>; img: GerberImage }[] = [];
  let profile: GerberImage | undefined;
  const hits: DrillHit[] = [];
  for (const f of files) {
    const role = classify(f.name, f.text);
    if (role.kind === 'copper') copper.push({ file: f, role, img: parseGerber(f.text) });
    else if (role.kind === 'profile') profile = parseGerber(f.text);
    else if (role.kind === 'drill') hits.push(...parseExcellon(f.text, role.plated));
  }
  if (!copper.length) throw new Error('no copper layers found (expected X2 FileFunction or KiCad/Altium file names)');
  copper.sort((a, b) => (a.role.side === 'top' ? -1e9 : a.role.side === 'bottom' ? 1e9 : a.role.index) -
    (b.role.side === 'top' ? -1e9 : b.role.side === 'bottom' ? 1e9 : b.role.index));
  for (const c of copper) for (const w of c.img.warnings) warnings.push(`${c.file.name}: ${w}`);

  const rows = stackup ?? defaultStackup(copper.length);
  const cuRows = rows.filter((r) => r.type === 'copper');
  if (cuRows.length !== copper.length) throw new Error(`stackup has ${cuRows.length} copper layers but ${copper.length} copper files were found`);

  // layers bottom-up with z, as the KiCad importer builds them
  const layers: StackLayer[] = [...rows].reverse().map((r) => ({ name: r.name, type: r.type, thickness: r.thickness,
    epsR: r.epsR, lossTangent: r.lossTangent, roughness: r.roughness, mask: r.mask, z0: 0, z1: 0 }));
  let z = 0;
  for (const l of layers) { if (l.type === 'dielectric') { l.z0 = z; z += l.thickness; l.z1 = z; } else { l.z0 = z; l.z1 = z; } }
  const layerName = (k: number) => cuRows[k].name; // k: top-down copper index

  // ---- vias: plated hits; pad flashes centred on them are the via lands ----
  const plated = hits.filter((h) => h.plated && h.dia > 0);
  const near = (a: Pt, b: Pt, tol: number) => Math.abs(a[0] - b[0]) < tol && Math.abs(a[1] - b[1]) < tol;
  const tracks: Track[] = [];
  const polys: CopperPoly[] = [];
  const trackNames: (string | undefined)[] = [], polyNames: (string | undefined)[] = []; // X2 names if present
  const vias: Via[] = [];
  const viaPad = new Map<number, number>(); // via index -> land diameter
  const viaLayers = new Map<number, Set<number>>();
  copper.forEach((c, k) => {
    for (const t of c.img.tracks) { tracks.push({ a: t.a, b: t.b, width: t.width, layer: layerName(k), net: 0 }); trackNames.push(t.net); }
    for (const p of c.img.polys) {
      if (p.flash && p.flash.dia > 0) {
        const v = plated.findIndex((h) => near(h.at, p.flash!.at, 1e-6));
        if (v >= 0) {
          viaPad.set(v, Math.max(viaPad.get(v) ?? 0, p.flash.dia));
          if (!viaLayers.has(v)) viaLayers.set(v, new Set());
          viaLayers.get(v)!.add(k);
          continue;
        }
      }
      polys.push({ pts: p.pts, layer: layerName(k), net: 0 });
      polyNames.push(p.net);
    }
  });
  plated.forEach((h, v) => {
    const ks = [...(viaLayers.get(v) ?? new Set([0, copper.length - 1]))].sort((a, b) => a - b);
    // a barrel runs between the outermost lands it has (through-hole when only drilled)
    const top = ks[0], bot = ks[ks.length - 1];
    vias.push({ at: h.at, drill: h.dia, size: viaPad.get(v) ?? h.dia + 0.2e-3,
      layers: [layerName(top), layerName(bot === top ? copper.length - 1 : bot)], net: 0 });
  });

  // ---- outline ----
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  const grow = (p: Pt) => { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); };
  if (profile && (profile.tracks.length || profile.polys.length)) {
    // profile lines are drawn with a thin aperture; the board edge is the line centre
    for (const t of profile.tracks) { grow(t.a); grow(t.b); }
    for (const p of profile.polys) p.pts.forEach(grow);
  } else {
    warnings.push('no board profile; using the copper bounding box as the outline');
    for (const t of tracks) { grow(t.a); grow(t.b); }
    for (const p of polys) p.pts.forEach(grow);
  }

  // ---- nets: X2 names where present, otherwise connectivity ----
  const nets = assignNets(tracks, polys, vias, layers, [...trackNames, ...polyNames]);
  return {
    board: { layers, copperLayers: layers.filter((l) => l.type === 'copper').map((l) => l.name), tracks, vias, polys,
      outline: { min: [x0, y0], max: [x1, y1] }, nets, thickness: z },
    copperFiles: copper.map((c) => c.file.name), warnings,
  };
}

// ---- connectivity ----

type Obj = { layer: string; box: [number, number, number, number]; track?: Track; poly?: CopperPoly };

function segDist(p: Pt, a: Pt, b: Pt): number {
  const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy;
  const u = L2 > 0 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) : 0;
  return Math.hypot(a[0] + u * dx - p[0], a[1] + u * dy - p[1]);
}
function segSeg(a: Pt, b: Pt, c: Pt, d: Pt): number {
  const o = (p: Pt, q: Pt, r: Pt) => Math.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]));
  if (o(a, b, c) * o(a, b, d) < 0 && o(c, d, a) * o(c, d, b) < 0) return 0;
  return Math.min(segDist(a, c, d), segDist(b, c, d), segDist(c, a, b), segDist(d, a, b));
}
function polyEdges(pts: Pt[]): [Pt, Pt][] { return pts.map((p, i) => [p, pts[(i + 1) % pts.length]]); }

function touches(A: Obj, B: Obj): boolean {
  const tol = 1e-9;
  if (A.track && B.track) return segSeg(A.track.a, A.track.b, B.track.a, B.track.b) <= (A.track.width + B.track.width) / 2 + tol;
  const t = A.track ?? B.track, p = (A.poly ?? B.poly)!;
  if (t) {
    if (inPoly(t.a, p.pts) || inPoly(t.b, p.pts)) return true;
    return polyEdges(p.pts).some(([c, d]) => segSeg(t.a, t.b, c, d) <= t.width / 2 + tol);
  }
  const P1 = A.poly!, P2 = B.poly!;
  if (inPoly(P1.pts[0], P2.pts) || inPoly(P2.pts[0], P1.pts)) return true;
  const e2 = polyEdges(P2.pts);
  return polyEdges(P1.pts).some(([a, b]) => e2.some(([c, d]) => segSeg(a, b, c, d) <= tol));
}

function assignNets(tracks: Track[], polys: CopperPoly[], vias: Via[], layers: StackLayer[], names: (string | undefined)[]): Map<number, string> {
  const objs: Obj[] = [];
  for (const t of tracks) {
    const h = t.width / 2;
    objs.push({ layer: t.layer, track: t, box: [Math.min(t.a[0], t.b[0]) - h, Math.min(t.a[1], t.b[1]) - h, Math.max(t.a[0], t.b[0]) + h, Math.max(t.a[1], t.b[1]) + h] });
  }
  for (const p of polys) {
    const xs = p.pts.map((q) => q[0]), ys = p.pts.map((q) => q[1]);
    objs.push({ layer: p.layer, poly: p, box: [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)] });
  }
  const parent = objs.map((_, i) => i);
  const find = (i: number): number => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
  const union = (a: number, b: number) => { const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; };
  const overlap = (a: Obj['box'], b: Obj['box']) => a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3];
  // sweep by x to keep this near-linear on real boards
  const order = objs.map((_, i) => i).sort((a, b) => objs[a].box[0] - objs[b].box[0]);
  for (let s = 0; s < order.length; s++) {
    const A = objs[order[s]];
    for (let r = s + 1; r < order.length && objs[order[r]].box[0] <= A.box[2]; r++) {
      const B = objs[order[r]];
      if (A.layer === B.layer && overlap(A.box, B.box) && touches(A, B)) union(order[s], order[r]);
    }
  }
  // X2 net names join objects that share a name (e.g. split pours of one net)
  const byName = new Map<string, number>();
  names.forEach((nm, i) => { if (!nm) return; const j = byName.get(nm); if (j === undefined) byName.set(nm, i); else union(i, j); });
  // vias join every layer between their lands
  const zOf = new Map(layers.map((l) => [l.name, l.z0]));
  const viaRoot: number[] = [];
  for (const v of vias) {
    const za = zOf.get(v.layers[0])!, zb = zOf.get(v.layers[1])!;
    const lo = Math.min(za, zb), hi = Math.max(za, zb);
    let first = -1;
    objs.forEach((o, i) => {
      const zl = zOf.get(o.layer)!;
      if (zl < lo - 1e-12 || zl > hi + 1e-12) return;
      const hit = o.track ? segDist(v.at, o.track.a, o.track.b) <= o.track.width / 2 + v.size / 2 : inPoly(v.at, o.poly!.pts);
      if (!hit) return;
      if (first < 0) first = i; else union(first, i);
    });
    viaRoot.push(first);
  }
  // number the components; name them from X2 attributes when available
  const id = new Map<number, number>();
  const nets = new Map<number, string>();
  nets.set(0, '');
  objs.forEach((o, i) => {
    const r = find(i);
    if (!id.has(r)) id.set(r, id.size + 1);
    const n = id.get(r)!;
    if (o.track) o.track.net = n; else o.poly!.net = n;
    if (names[i] && !nets.get(n)) nets.set(n, names[i]!);
  });
  for (const [, n] of id) if (!nets.has(n)) nets.set(n, `Net-${n}`);
  vias.forEach((v, k) => { if (viaRoot[k] >= 0) v.net = id.get(find(viaRoot[k]))!; });
  return nets;
}
