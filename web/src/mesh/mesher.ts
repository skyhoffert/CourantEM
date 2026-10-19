/**
 * Rule-based mesher (build step 5). Produces the non-uniform node lines of one
 * axis from geometric features, following the milestone-1 rules:
 *
 *  1. Metal edges: zero-thickness copper uses the thirds rule (a line h/3 inside
 *     the metal and 2h/3 outside, docs/design/0002); thick copper faces and
 *     dielectric interfaces get a line exactly on them.
 *  2. A minimum number of cells across every gap and every dielectric layer.
 *  3. Maximum cell lambda/20 at fmax in the densest material (per region).
 *  4. Neighbouring cells change by at most `ratio` (default 1.3).
 *  5. Manual overrides: extra fixed lines and per-interval maximum cell sizes.
 */

export type MetalSide = 'lo' | 'hi';

export interface AxisFeatures {
  /** domain extent (the CPML is added outside by the caller) */
  lo: number;
  hi: number;
  /** zero-thickness metal edges: coordinate and which side the metal is on */
  sheetEdges?: { x: number; metal: MetalSide }[];
  /** lines that must lie exactly on the grid (thick copper faces, layer interfaces) */
  planes?: number[];
  /** intervals that need at least n cells across (CPW gaps, dielectric layers) */
  spans?: { a: number; b: number; n: number }[];
  /** regions with their own maximum cell (e.g. lambda/20 in the substrate) */
  maxCell?: { a: number; b: number; h: number }[];
  /** points where cells should be smallest (defaults to every edge and plane) */
  refine?: number[];
  /** manual extra fixed lines */
  extra?: number[];
}

export interface MeshOptions {
  hmin: number;    // target cell at metal edges (m)
  hmax: number;    // global maximum cell (m)
  ratio?: number;  // max neighbour ratio (default 1.3)
}

export interface AxisMesh {
  nodes: number[];
  report: { cells: number; minCell: number; maxCell: number; maxRatio: number; warnings: string[] };
}

const EPS = 1e-12;

/** lambda/20 at fmax in a material of relative permittivity er. */
export function lambdaOver20(fmax: number, er: number): number {
  return 299792458.0 / (fmax * Math.sqrt(er)) / 20;
}

export function meshAxis(feat: AxisFeatures, opt: MeshOptions): AxisMesh {
  const r = opt.ratio ?? 1.3;
  const h0 = opt.hmin;
  const warnings: string[] = [];
  // ---- 1: fixed lines ----
  const fixed: { x: number; pri: number }[] = [{ x: feat.lo, pri: 3 }, { x: feat.hi, pri: 3 }];
  const refine: number[] = [...(feat.refine ?? [])];
  for (const p of feat.planes ?? []) { fixed.push({ x: p, pri: 2 }); if (!feat.refine) refine.push(p); }
  for (const e of feat.sheetEdges ?? []) {
    const s = e.metal === 'lo' ? -1 : 1; // direction into the metal
    fixed.push({ x: e.x + s * h0 / 3, pri: 1 }, { x: e.x - s * 2 * h0 / 3, pri: 1 });
    if (!feat.refine) refine.push(e.x);
  }
  for (const x of feat.extra ?? []) fixed.push({ x, pri: 2 });
  fixed.sort((a, b) => a.x - b.x);
  // merge lines closer than hmin/4, keeping the higher-priority one
  const lines: number[] = [];
  const pri: number[] = [];
  for (const f of fixed) {
    if (f.x < feat.lo - EPS || f.x > feat.hi + EPS) continue;
    const n = lines.length;
    if (n && f.x - lines[n - 1] < h0 / 4) {
      if (f.pri > pri[n - 1]) { lines[n - 1] = f.x; pri[n - 1] = f.pri; }
      else if (f.pri === pri[n - 1] && f.pri >= 2 && f.x - lines[n - 1] > EPS) {
        warnings.push(`lines ${lines[n - 1]} and ${f.x} closer than hmin/4; merged`);
      }
      continue;
    }
    lines.push(f.x); pri.push(f.pri);
  }
  // ---- 2/3: size field ----
  // Every constraint is a cap that relaxes away from where it applies with slope
  // (ratio - 1), so the field is graded by construction and the fill below never
  // needs to split cells afterwards.
  const g = 0.8 * (r - 1); // margin for the discrete steps and the per-interval rounding
  const dist = (x: number, a: number, b: number) => (x < a ? a - x : x > b ? x - b : 0);
  const size = (x: number): number => {
    let s = opt.hmax;
    for (const m of feat.maxCell ?? []) s = Math.min(s, m.h + g * dist(x, m.a, m.b));
    for (const sp of feat.spans ?? []) s = Math.min(s, (sp.b - sp.a) / sp.n + g * dist(x, sp.a, sp.b));
    for (const q of refine) s = Math.min(s, h0 + g * Math.abs(x - q));
    return Math.max(s, h0 * 0.5);
  };
  // ---- fill each interval ----
  // n = ceil(integral of dx / s(x)); nodes at equal steps of that integral, so
  // every cell is at most s(x) and follows the field's grading.
  const nodes: number[] = [lines[0]];
  for (let q = 0; q + 1 < lines.length; q++) {
    const a = lines[q], b = lines[q + 1];
    const M = Math.max(64, Math.ceil((b - a) / (h0 / 8)));
    const cum: number[] = [0];
    for (let m = 0; m < M; m++) {
      const x0 = a + ((b - a) * m) / M, x1 = a + ((b - a) * (m + 1)) / M;
      cum.push(cum[m] + (x1 - x0) * 0.5 * (1 / size(x0) + 1 / size(x1)));
    }
    const I = cum[M];
    const n = Math.max(1, Math.ceil(I - 1e-6));
    let m = 0;
    for (let k = 1; k < n; k++) {
      const target = (I * k) / n;
      while (cum[m + 1] < target) m++;
      const t = (target - cum[m]) / (cum[m + 1] - cum[m]);
      nodes.push(a + ((b - a) * (m + t)) / M);
    }
    nodes.push(b);
  }
  // ---- 4: the size field keeps neighbour ratios near `ratio`; the rescaling of each
  // interval to end on its fixed line can push a pair slightly over, which is reported
  // ---- 2 again: make sure every span got its cell count ----
  for (const sp of feat.spans ?? []) {
    const inside = nodes.filter((v) => v > sp.a + EPS && v < sp.b - EPS).length + 1;
    if (inside < sp.n) warnings.push(`span [${sp.a}, ${sp.b}] has ${inside} cells, wanted ${sp.n}`);
  }
  let minCell = Infinity, maxCell = 0, maxRatio = 1;
  for (let i = 0; i + 1 < nodes.length; i++) {
    const c = nodes[i + 1] - nodes[i];
    minCell = Math.min(minCell, c); maxCell = Math.max(maxCell, c);
    if (i > 0) { const p = nodes[i] - nodes[i - 1]; maxRatio = Math.max(maxRatio, c / p, p / c); }
  }
  return { nodes, report: { cells: nodes.length - 1, minCell, maxCell, maxRatio, warnings } };
}

/** Append n CPML cells of size h on both ends. */
export function withPml(nodes: number[], n: number, h: number, lo = true, hi = true): number[] {
  const out = [...nodes];
  if (lo) for (let i = 1; i <= n; i++) out.unshift(nodes[0] - i * h);
  if (hi) for (let i = 1; i <= n; i++) out.push(nodes[nodes.length - 1] + i * h);
  return out;
}

/** Pre-run estimate: cells, memory, time step, steps and work for a run of tSim seconds. */
export function estimate(x: number[], y: number[], z: number[], opts: { tSim: number; backend: 'cpu' | 'gpu'; poles?: number; npml?: number }) {
  const n = [x.length - 1, y.length - 1, z.length - 1];
  const N = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
  const min = (a: number[]) => { let m = Infinity; for (let i = 0; i + 1 < a.length; i++) m = Math.min(m, a[i + 1] - a[i]); return m; };
  const dt = 0.99 / (299792458.0 * Math.sqrt(1 / min(x) ** 2 + 1 / min(y) ** 2 + 1 / min(z) ** 2));
  const steps = Math.ceil(opts.tSim / dt);
  // per point: E, H, ca, cb (12 values) (+ cell eps, sigma and a PEC flag in the C path);
  // CPML auxiliaries only on PML planes: 4 per axis on 2·npml planes
  const np = opts.npml ?? 8;
  const per = [(n[1] + 1) * (n[2] + 1), (n[0] + 1) * (n[2] + 1), (n[0] + 1) * (n[1] + 1)];
  const psi = per.reduce((acc, p) => acc + 4 * 2 * np * p, 0);
  const bytes = opts.backend === 'cpu' ? N * 8 * (12 + 2) + N + 8 * psi : 4 * (N * 12 + psi);
  return { cells: n[0] * n[1] * n[2], points: N, dt, steps, bytes, cellSteps: n[0] * n[1] * n[2] * steps };
}
