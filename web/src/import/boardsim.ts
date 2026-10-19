/**
 * Board model -> Courant simulation: rule-based mesh, rasterised copper, via
 * barrels, dielectric layers, and mode ports where tracks leave the board.
 *
 * Milestone-1 scope: planar boards whose port lines run along x and leave the
 * board through its left/right edges. The board is extruded past those edges
 * by a feed length so the modal source and monitor sit on a uniform line, and
 * the line then runs into the x CPML.
 */
import { Axis, Boundary, Face, type Courant, type Sim } from '../solver/courant.ts';
import { lambdaOver20, meshAxis, withPml } from '../mesh/mesher.ts';
import { copperAt, copperEdges, copperRow, inRow, type Board, type Track } from './kicad.ts';

export interface PortSpec {
  net: string | number;       // net name or number of the track that leaves the board
  edge: 'xlo' | 'xhi';
  /**
   * Port plane, measured into the board from its edge (m), default 0. The board
   * is cut there (connector launches beyond it are dropped), the line's cross-
   * section at the plane is extruded outwards into the feed, and S-parameters
   * are referred to the plane — e.g. the plane of a TRL-calibrated measurement.
   */
  offset?: number;
}

export interface BoardSimOptions {
  fmax: number;
  df: number;
  hmin: number;             // min cell at copper edges (m)
  feed?: number;            // extrusion past the board edge for each port (m), default 4 mm
  monitor?: number;         // source-to-monitor distance (m), default 2.5 mm
  airAbove?: number;        // default 1.5 mm
  airSide?: number;         // air beyond the board's y edges, default 1 mm
  npml?: number;            // default 8
  lossy?: boolean;          // DS dielectrics from the stackup + copper loss
  /**
   * Model copper with its stackup thickness (PEC bars extruded upwards from each
   * layer's surface) instead of zero-thickness sheets. Matters for Z0: 35 um copper
   * lowers the test GCPW from 49.9 to 46.8 ohm. Side edges keep the thirds rule
   * (the metal ends h/3 inside the drawn edge), which cuts the Z0 error at 20 um
   * cells from -1.9% (lines on the faces) to +0.5%; the top face gets a grid line.
   */
  thickCopper?: boolean;
  /**
   * Simulate only a band of the board around the port net: y within this distance
   * (m) of the port tracks. Copper and vias outside are dropped and the band's
   * edges become air + CPML. 0 / undefined = whole board. Makes large boards
   * tractable when only one line matters.
   */
  yWindow?: number;
  ports: PortSpec[];
  excite?: number;          // which port to excite (default 0)
}

export interface BoardSim {
  sim: Sim;
  freqs: number[];
  ports: number[];
  cells: [number, number, number];
  tStop: number;
  mesh: { x: number[]; y: number[]; z: number[] };
}

const nearest = (a: number[], v: number) => { let b = 0; for (let i = 1; i < a.length; i++) if (Math.abs(a[i] - v) < Math.abs(a[b] - v)) b = i; return b; };

/** The board restricted to a y band around the port tracks (see BoardSimOptions.yWindow). */
function cropToNet(b: Board, o: BoardSimOptions): Board {
  const W = o.yWindow!;
  const ys: number[] = [];
  for (const ps of o.ports) {
    const plane = ps.edge === 'xlo' ? b.outline.min[0] + (ps.offset ?? 0) : b.outline.max[0] - (ps.offset ?? 0);
    const t = trackAcross(b, netNumber(b, ps.net), plane);
    if (t) ys.push(t.a[1] - t.width / 2, t.a[1] + t.width / 2);
  }
  if (!ys.length) return b;
  const lo = Math.max(b.outline.min[1], Math.min(...ys) - W), hi = Math.min(b.outline.max[1], Math.max(...ys) + W);
  const yIn = (y0: number, y1: number) => y1 >= lo && y0 <= hi;
  return {
    ...b,
    tracks: b.tracks.filter((t) => yIn(Math.min(t.a[1], t.b[1]) - t.width / 2, Math.max(t.a[1], t.b[1]) + t.width / 2)),
    polys: b.polys.filter((q) => yIn(Math.min(...q.pts.map((p) => p[1])), Math.max(...q.pts.map((p) => p[1])))),
    vias: b.vias.filter((v) => v.at[1] >= lo && v.at[1] <= hi),
    outline: { min: [b.outline.min[0], lo], max: [b.outline.max[0], hi] },
  };
}

/** Effective index guess for a strip of width w at height z0 above the floor. */
function neffGuess(b: Board, z0: number, w: number): number {
  const below = b.layers.filter((l) => l.type === 'dielectric' && !l.mask && l.z1 <= z0 + 1e-12);
  const h = Math.max(1e-6, below.reduce((a, l) => a + l.thickness, 0));
  const er = below.length ? below.reduce((a, l) => a + (l.epsR ?? 4.5) * l.thickness, 0) / h : 1;
  const ee = (er + 1) / 2 + ((er - 1) / 2) / Math.sqrt(1 + 12 * h / w);
  return Math.sqrt(ee);
}

/** An x-directed track of `net` that reaches or crosses the plane x = xc. */
export function trackAcross(b: Board, net: number, xc: number): Track | undefined {
  return b.tracks.find((t) => t.net === net && Math.abs(t.a[1] - t.b[1]) < 1e-9 &&
    Math.min(t.a[0], t.b[0]) <= xc + 1e-6 && Math.max(t.a[0], t.b[0]) >= xc - 1e-6);
}

function netNumber(b: Board, n: string | number): number {
  if (typeof n === 'number') return n;
  for (const [k, v] of b.nets) if (v === n) return k;
  throw new Error(`net ${n} not found`);
}

/** Mesh and geometry decisions for a board run, without building the simulation
 * (the app uses it for the pre-run estimate). */
export function planBoard(board: Board, o: BoardSimOptions) {
  const b = o.yWindow ? cropToNet(board, o) : board;
  const feed = o.feed ?? 4e-3, mon = o.monitor ?? 2.5e-3, airAbove = o.airAbove ?? 1.5e-3, airSide = o.airSide ?? 1e-3;
  const npml = o.npml ?? 8;
  const [bx0, by0] = b.outline.min, [bx1, by1] = b.outline.max;
  // where the band cut through the board, the cross-section continues into the y CPML
  // (no air margin): parallel-plate waves between the grounds leave the band instead of
  // ringing in an open-edged, nearly lossless cavity
  const cutLo = by0 > board.outline.min[1] + 1e-9, cutHi = by1 < board.outline.max[1] - 1e-9;
  const yLo = cutLo ? -Infinity : by0, yHi = cutHi ? Infinity : by1;
  const clampY = (v: number) => Math.min(by1, Math.max(by0, v));
  const dielectrics = b.layers.filter((l) => l.type === 'dielectric');
  const erMax = Math.max(1, ...dielectrics.map((l) => l.epsR ?? 1));
  const hAir = Math.min(lambdaOver20(o.fmax, 1), 200e-6);
  const hDiel = lambdaOver20(o.fmax, erMax);
  const copper = b.layers.filter((l) => l.type === 'copper');
  const bottom = copper[0];
  // the bottom copper is the domain's PEC floor if it covers the board
  const floorPlane = (() => {
    let hits = 0, tot = 0;
    for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) {
      tot++;
      if (copperAt(b, bottom.name, [bx0 + (bx1 - bx0) * (i + 0.5) / 20, by0 + (by1 - by0) * (j + 0.5) / 20])) hits++;
    }
    return hits / tot > 0.95;
  })();
  if (!floorPlane) throw new Error('milestone 1 needs a solid bottom ground plane');
  const hasLo = o.ports.some((p) => p.edge === 'xlo'), hasHi = o.ports.some((p) => p.edge === 'xhi');
  // cut planes: the board between them is simulated; outside, the cross-section at the plane is extruded
  const offLo = o.ports.find((p) => p.edge === 'xlo')?.offset ?? 0, offHi = o.ports.find((p) => p.edge === 'xhi')?.offset ?? 0;
  const cx0 = bx0 + (hasLo ? offLo : 0), cx1 = bx1 - (hasHi ? offHi : 0);
  if (!(cx1 - cx0 > 1e-4)) throw new Error('port planes leave less than 0.1 mm of board');
  const X0 = cx0 - (hasLo ? feed : 0), X1 = cx1 + (hasHi ? feed : 0);
  const clampX = (x: number) => Math.min(cx1, Math.max(cx0, x)); // extrusion past the port planes

  // ---- mesh ----
  const top = copper[copper.length - 1];
  const zTop = top.z1 + (o.thickCopper ? top.thickness : 0);
  const viaRefX: number[] = [], viaRefY: number[] = [];
  const inside = b.vias.filter((v) => v.at[0] > cx0 && v.at[0] < cx1); // vias beyond the port planes are cut away
  for (const v of inside) { viaRefX.push(v.at[0]); viaRefY.push(v.at[1]); }
  const sheetY = copper.flatMap((c) => copperEdges(b, c.name, 1)).filter((e) => e.x > by0 - 1e-9 && e.x < by1 + 1e-9);
  const sheetX = copper.flatMap((c) => copperEdges(b, c.name, 0)).filter((e) => e.x > cx0 + 1e-9 && e.x < cx1 - 1e-9);
  // merge duplicates of the same edge from several layers
  const uniq = (arr: { x: number; metal: 'lo' | 'hi' }[]) => {
    const m = new Map<string, { x: number; metal: 'lo' | 'hi' }>();
    for (const e of arr) m.set(`${e.x.toFixed(9)}${e.metal}`, e);
    return [...m.values()];
  };
  const viaR = inside.length ? Math.min(...inside.map((v) => v.drill / 2)) : 0;
  const hVia = viaR ? Math.max(o.hmin, viaR / 2) : Infinity;
  const thick = !!o.thickCopper;
  const cuTop = (l: { z0: number; thickness: number }) => l.z0 + (thick ? l.thickness : 0);
  const ymesh = meshAxis({
    lo: cutLo ? by0 : by0 - airSide, hi: cutHi ? by1 : by1 + airSide,
    sheetEdges: uniq(sheetY), planes: [by0, by1],
    refine: [...uniq(sheetY).map((e) => e.x), ...viaRefY],
    maxCell: [{ a: by0, b: by1, h: hDiel }, ...viaRefY.map((y) => ({ a: y - viaR, b: y + viaR, h: hVia }))],
  }, { hmin: o.hmin, hmax: hAir });
  const xmesh = meshAxis({
    lo: X0, hi: X1,
    sheetEdges: uniq(sheetX), planes: [cx0, cx1],
    refine: [...uniq(sheetX).map((e) => e.x), ...viaRefX],
    maxCell: [{ a: X0, b: X1, h: Math.min(hDiel, 4 * o.hmin) }, ...viaRefX.map((x) => ({ a: x - viaR, b: x + viaR, h: hVia }))],
  }, { hmin: o.hmin, hmax: hAir });
  const zplanes = [...new Set([...b.layers.flatMap((l) => [l.z0, l.z1]), ...copper.filter((c) => c !== bottom).map(cuTop)])];
  const zmesh = meshAxis({
    lo: 0, hi: zTop + airAbove, planes: zplanes,
    spans: dielectrics.map((l) => ({ a: l.z0, b: l.z1, n: l.mask ? 1 : 4 })),
    maxCell: dielectrics.map((l) => ({ a: l.z0, b: l.z1, h: hDiel })),
    refine: [...copper.map((c) => c.z0), ...(thick ? copper.filter((c) => c !== bottom).map(cuTop) : [])],
  }, { hmin: o.hmin, hmax: hAir });
  const hp = hAir / 2;
  const x = withPml(xmesh.nodes, npml, hp, hasLo, hasHi);
  const y = withPml(ymesh.nodes, npml, hp);
  const z = withPml(zmesh.nodes, npml, hp, false, true);
  const nx = x.length - 1, ny = y.length - 1, nz = z.length - 1;
  return { b, feed, mon, npml, bx0, bx1, by0, by1, dielectrics, erMax, copper, bottom, hasLo, hasHi, cx0, cx1, X0, X1, clampX, yLo, yHi, clampY, top, zTop, inside, thick, cuTop, hp, x, y, z, nx, ny, nz };
}

export function buildBoardSim(cem: Courant, board: Board, o: BoardSimOptions): BoardSim {
  const { b, feed, mon, npml, bx0, bx1, by0, by1, dielectrics, erMax, copper, bottom, hasLo, hasHi, cx0, cx1, X0, X1, clampX, yLo, yHi, clampY, top, zTop, inside, thick, cuTop, hp, x, y, z, nx, ny, nz } = planBoard(board, o);

  const sim = cem.createSim(x, y, z);
  if (hasLo) sim.setBoundary(Face.XLo, Boundary.CPML, npml);
  if (hasHi) sim.setBoundary(Face.XHi, Boundary.CPML, npml);
  sim.setBoundary(Face.YLo, Boundary.CPML, npml);
  sim.setBoundary(Face.YHi, Boundary.CPML, npml);
  sim.setBoundary(Face.ZHi, Boundary.CPML, npml);
  const big = 1;

  // ---- dielectrics (extruded in x into the CPML) ----
  if (o.lossy) {
    const fits = dielectrics.map((l) => cem.fitDjordjevicSarkar(l.epsR ?? 4.5, l.lossTangent ?? 0, 10e9, o.df, o.fmax, 8));
    // all materials share one pole set: refit is not needed because the pole
    // frequencies depend only on the band
    sim.setDebyePoles(fits[0].fp);
    dielectrics.forEach((l, i) => sim.materialBoxDebye(-big, big, Math.max(-big, yLo), Math.min(big, yHi), l.z0, l.z1, fits[i], 0));
    // roughness of the port tracks' copper (SIBC is one model for the whole board)
    const rq = Math.max(0, ...copper.map((l) => l.roughness ?? 0));
    sim.setConductorLoss(5.8e7, o.df, o.fmax, rq > 0 ? 1 : 0, rq);
  } else {
    for (const l of dielectrics) sim.materialBox(-big, big, Math.max(-big, yLo), Math.min(big, yHi), l.z0, l.z1, l.epsR ?? 4.5, 0);
  }

  // ---- copper sheets (rasterised edge by edge; bottom layer is the PEC floor) ----
  for (const c of copper) {
    if (c === bottom) continue;
    const k = nearest(z, c.z0);
    const mu = new Uint8Array(nx * (ny + 1)), mv = new Uint8Array((nx + 1) * ny);
    // copper sampled row by row (same answers as copperAt, far faster on real boards)
    const xc = x.map(clampX), xm = x.slice(0, -1).map((v, i) => clampX((v + x[i + 1]) / 2));
    const node: Uint8Array[] = [];
    for (let j = 0; j <= ny; j++) {
      const r = y[j] < yLo - 1e-12 || y[j] > yHi + 1e-12 ? null : copperRow(b, c.name, clampY(y[j]));
      node.push(r ? Uint8Array.from(xc, (v) => +inRow(r, v)) : new Uint8Array(nx + 1));
      // edges along x (u = x for a z-normal plane), fully inside copper
      if (r) for (let i = 0; i < nx; i++) if (node[j][i] && node[j][i + 1] && inRow(r, xm[i])) mu[i * (ny + 1) + j] = 1;
    }
    for (let j = 0; j < ny; j++) {
      if (y[j] < yLo - 1e-12 || y[j + 1] > yHi + 1e-12) continue;
      const r = copperRow(b, c.name, clampY((y[j] + y[j + 1]) / 2));
      for (let i = 0; i <= nx; i++) if (node[j][i] && node[j + 1][i] && inRow(r, xc[i])) mv[i * ny + j] = 1;
    }
    if (!thick) { cemEdges(cem, sim, k, mu, mv); continue; }
    // thick copper: cells between the layer surface and its top face whose centre is on copper
    const kt = nearest(z, cuTop(c));
    const nk = kt - k;
    if (nk <= 0) { cemEdges(cem, sim, k, mu, mv); continue; }
    const cellsM = new Uint8Array(nx * ny * nk);
    const xcen = x.slice(0, -1).map((v, i) => clampX((v + x[i + 1]) / 2));
    for (let j = 0; j < ny; j++) {
      const yc = (y[j] + y[j + 1]) / 2;
      if (yc < yLo || yc > yHi) continue;
      const r = copperRow(b, c.name, clampY(yc));
      for (let i = 0; i < nx; i++) if (inRow(r, xcen[i])) cellsM.fill(1, (i * ny + j) * nk, (i * ny + j + 1) * nk);
    }
    cemCells(cem, sim, 0, 0, k, nx, ny, nk, cellsM);
  }

  // ---- via barrels: solid PEC cells inside the drill radius, between their layers ----
  for (const v of inside) {
    const zl = Math.min(...v.layers.map((n) => b.layers.find((l) => l.name === n)?.z0 ?? 0));
    const zh = Math.max(...v.layers.map((n) => b.layers.find((l) => l.name === n)?.z0 ?? zTop));
    const r = v.drill / 2;
    const i0 = Math.max(0, nearest(x, v.at[0] - r) - 1), i1 = Math.min(nx - 1, nearest(x, v.at[0] + r) + 1);
    const j0 = Math.max(0, nearest(y, v.at[1] - r) - 1), j1 = Math.min(ny - 1, nearest(y, v.at[1] + r) + 1);
    const k0 = nearest(z, zl), k1 = nearest(z, zh);
    const ni = i1 - i0 + 1, nj = j1 - j0 + 1, nk = k1 - k0;
    if (nk <= 0) continue;
    const m = new Uint8Array(ni * nj * nk);
    for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
      const cx = (x[i0 + i] + x[i0 + i + 1]) / 2, cy = (y[j0 + j] + y[j0 + j + 1]) / 2;
      if ((cx - v.at[0]) ** 2 + (cy - v.at[1]) ** 2 <= r * r) for (let k = 0; k < nk; k++) m[(i * nj + j) * nk + k] = 1;
    }
    cemCells(cem, sim, i0, j0, k0, ni, nj, nk, m);
  }

  // ---- feed extensions: via rows become continuous walls ----
  // The board is extruded past its port edges so the modal source and monitor sit
  // on a uniform line. Discrete vias can't be extruded, so each via's y extent is
  // carried through the feed as a solid wall between its layers (the idealised
  // fence of the milestone-1 GCPW); otherwise the feed's top grounds would float.
  // Only ground vias (the bottom plane's net) near the port plane are carried through:
  // a via elsewhere, or of another net, must not become a wall across the feed.
  const floorNet = (() => {
    const count = new Map<number, number>();
    for (const q of b.polys) if (q.layer === bottom.name) count.set(q.net, (count.get(q.net) ?? 0) + q.pts.length);
    return [...count.entries()].sort((p1, p2) => p2[1] - p1[1])[0]?.[0];
  })();
  const nearPlane = 2e-3;
  const feedWalls = (xa: number, xb: number, plane: number) => {
    const i0 = nearest(x, xa), i1 = nearest(x, xb);
    if (i1 <= i0) return;
    const seen = new Set<string>();
    for (const v of inside) {
      if (Math.abs(v.at[0] - plane) > nearPlane || (floorNet !== undefined && v.net !== floorNet && v.net !== 0)) continue;
      const key = v.at[1].toFixed(7);
      if (seen.has(key)) continue;
      seen.add(key);
      const r = v.drill / 2;
      const zl = Math.min(...v.layers.map((n) => b.layers.find((l) => l.name === n)?.z0 ?? 0));
      const zh = Math.max(...v.layers.map((n) => b.layers.find((l) => l.name === n)?.z0 ?? zTop));
      sim.pecBox(x[i0], x[i1], v.at[1] - r, v.at[1] + r, zl, zh);
    }
  };
  if (hasLo) feedWalls(X0 - npml * hp, cx0, cx0);
  if (hasHi) feedWalls(cx1, X1 + npml * hp, cx1);

  // ---- ports ----
  const tau = Math.sqrt(6.9) / (Math.PI * 1.25 * o.fmax);
  const portIds: number[] = [];
  o.ports.forEach((ps) => {
    const net = netNumber(b, ps.net);
    const edgeX = ps.edge === 'xlo' ? cx0 : cx1;
    const t: Track | undefined = trackAcross(b, net, edgeX);
    if (!t) throw new Error(`no x-directed track of net ${ps.net} crosses the ${ps.edge} port plane`);
    const layer = b.layers.find((l) => l.name === t.layer)!;
    const yc = t.a[1], kh = nearest(z, layer.z0);
    const jl = nearest(y, yc - t.width / 2), jr = nearest(y, yc + t.width / 2);
    // current loop: dual nodes one cell outside each strip edge
    const ja = Math.max(0, (y[jl] <= yc - t.width / 2 ? jl : jl - 1) - 1);
    const jb = Math.min(ny - 1, (y[jr] >= yc + t.width / 2 ? jr : jr + 1));
    const dir = ps.edge === 'xlo' ? 1 : -1;
    const iSrc = ps.edge === 'xlo' ? npml + 2 : nx - npml - 2;
    const iMon = ps.edge === 'xlo' ? nearest(x, x[iSrc] + mon) : nearest(x, x[iSrc] - mon) - 1;
    const id = sim.addModePort({
      dir, iSrc, iMon, j0: 0, j1: ny, k0: 0, k1: nz, vj: nearest(y, yc), vk0: 0, vk1: kh, vsign: -1,
      ja, jb, ka: kh - 1, kb: thick ? Math.max(kh, nearest(z, cuTop(layer))) : kh, xRef: edgeX, f0: 0.5 * o.fmax,
      // quasi-static microstrip estimate (Hammerstad) from the strip's own substrate: a shift
      // close to the line mode keeps inverse iteration away from the slab's surface-wave mode
      neff0: neffGuess(b, layer.z0, t.width),
    });
    portIds.push(id);
  });
  sim.setModeExcitation(portIds[o.excite ?? 0], 4.5 * tau, tau);
  const nf = Math.round(o.fmax / o.df);
  const freqs = Array.from({ length: nf }, (_, i) => (i + 1) * o.df);
  sim.setFrequencies(freqs);
  sim.finalize(0.99);
  return { sim, freqs, ports: portIds, cells: [nx, ny, nz], tStop: 10.5 * tau, mesh: { x, y, z } };
}

// thin helpers around the raw exports for mask uploads
function cemEdges(cem: Courant, sim: Sim, k: number, mu: Uint8Array, mv: Uint8Array): void {
  const ex = cem.ex as unknown as { malloc(n: number): number; free(p: number): void;
    cem_pec_edges_plane(s: number, axis: number, m: number, mu: number, mv: number): number; memory: WebAssembly.Memory };
  const pu = ex.malloc(mu.length), pv = ex.malloc(mv.length);
  new Uint8Array(ex.memory.buffer, pu, mu.length).set(mu);
  new Uint8Array(ex.memory.buffer, pv, mv.length).set(mv);
  const rc = ex.cem_pec_edges_plane(sim.handleForImport(), Axis.Z, k, pu, pv);
  ex.free(pu); ex.free(pv);
  if (rc) throw new Error(`pec_edges_plane: ${cem.lastError()}`);
}

function cemCells(cem: Courant, sim: Sim, i0: number, j0: number, k0: number, ni: number, nj: number, nk: number,
                  m: Uint8Array): void {
  const ex = cem.ex as unknown as { malloc(n: number): number; free(p: number): void; memory: WebAssembly.Memory;
    cem_pec_cells(s: number, i0: number, j0: number, k0: number, ni: number, nj: number, nk: number, m: number): number };
  const p = ex.malloc(m.length);
  new Uint8Array(ex.memory.buffer, p, m.length).set(m);
  const rc = ex.cem_pec_cells(sim.handleForImport(), i0, j0, k0, ni, nj, nk, p);
  ex.free(p);
  if (rc) throw new Error(`pec_cells: ${cem.lastError()}`);
}
