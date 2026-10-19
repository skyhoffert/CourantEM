/**
 * Build step 1 scenario: hardcoded GCPW line between two 50-ohm lumped ports.
 * Mirrors solver/apps/gcpw_line.c exactly (same mesh arithmetic), so a browser
 * run can be compared number-for-number with the native app.
 */
import { Axis, Boundary, Face, sParams, type Complex, type Courant, type Sim } from '../solver/courant.ts';

export interface GcpwParams {
  hmin: number;   // m
  hmax: number;   // m
  dx: number;     // m
  len: number;    // m
  er: number;
  fmax: number;   // Hz
  df: number;     // Hz
  decay: number;
  ports: 'mode' | 'lumped';
  /** RO4350B-like Djordjevic–Sarkar substrate (Df 0.0037 at 10 GHz) and smooth copper sheets */
  lossy: boolean;
}

export const defaultGcpw: GcpwParams = {
  hmin: 20e-6, hmax: 100e-6, dx: 100e-6, len: 10e-3, er: 3.66, fmax: 40e9, df: 0.5e9, decay: 1e-5, ports: 'mode', lossy: false,
};

/** Graded 1D mesh, same algorithm and operation order as graded() in gcpw_line.c */
export function graded(fix: number[], ref: number[], hmin: number, hmax: number, ratio: number): number[] {
  const pts = [...fix].sort((a, b) => a - b);
  const out: number[] = [pts[0]];
  const sizeAt = (x: number): number => {
    let d = 1e300;
    for (const r of ref) if (Math.abs(x - r) < d) d = Math.abs(x - r);
    return Math.min(hmax, hmin + (ratio - 1) * d);
  };
  for (let s = 0; s + 1 < pts.length; s++) {
    const a = pts[s], b = pts[s + 1];
    if (b - a < 1e-15) continue;
    const seg: number[] = [];
    let x = a;
    for (;;) {
      let h = sizeAt(x);
      h = Math.min(h, sizeAt(x + h));
      if (x + h >= b - 0.5 * h) break;
      x += h;
      seg.push(x);
    }
    if (seg.length) {
      const last = seg[seg.length - 1];
      const prev = seg.length > 1 ? last - seg[seg.length - 2] : last - a;
      const k = (b - a) / ((last - a) + prev);
      for (const v of seg) out.push(a + (v - a) * k);
    }
    out.push(b);
  }
  return out;
}

function nearest(a: number[], x: number): number {
  let b = 0;
  for (let i = 1; i < a.length; i++) if (Math.abs(a[i] - x) < Math.abs(a[b] - x)) b = i;
  return b;
}

export interface GcpwSetup {
  mesh: { x: number[]; y: number[]; z: number[]; kh: number };
  sim: Sim;
  freqs: number[];
  ports: [number, number];
  portType: 'mode' | 'lumped';
  cells: [number, number, number];
  tStop: number;
}

/** Mesh lines of the built-in GCPW (no solver needed; used by the pre-run estimator). */
export function gcpwMesh(p: GcpwParams) {
  const h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6;
  const a = w / 2, b = w / 2 + s, yf = b + fence;
  const yext = yf + 1.0e-3, ztop = h + 1.5e-3, ext = 0.5e-3, npml = 8;
  const t1 = p.hmin / 3, t2 = 2 * p.hmin / 3;
  const Y = graded([-yext, -yf, -b - t1, -b + t2, -a - t2, -a + t1, 0, a - t1, a + t2, b - t2, b + t1, yf, yext],
    [-b, -a, a, b], p.hmin, p.hmax, 1.3);
  const Z = graded([0, h, ztop], [h], p.hmin, p.hmax, 1.3);
  const Y2: number[] = [];
  for (let i = npml; i >= 1; i--) Y2.push(Y[0] - i * p.hmax);
  Y2.push(...Y);
  for (let i = 1; i <= npml; i++) Y2.push(Y[Y.length - 1] + i * p.hmax);
  for (let i = 1; i <= npml; i++) Z.push(ztop + i * p.hmax);
  const X: number[] = [];
  const nlen = Math.round(p.len / p.dx), next = Math.round((p.ports === 'mode' ? 4.0e-3 : ext) / p.dx);
  for (let i = -next - npml; i <= nlen + next + npml; i++) X.push(i * p.dx);
  return { x: X, y: Y2, z: Z };
}

export function buildGcpw(cem: Courant, p: GcpwParams): GcpwSetup {
  const h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6;
  const a = w / 2, b = w / 2 + s, yf = b + fence;
  const yext = yf + 1.0e-3, ztop = h + 1.5e-3, ext = 0.5e-3;
  const npml = 8;
  // the C app uses lround(); the arguments here are never at .5 boundaries
  const lround = (v: number) => Math.round(v);

  // thirds rule at the metal edges (see gcpw_line.c / docs/design/0002-edge-meshing.md)
  const t1 = p.hmin / 3, t2 = 2 * p.hmin / 3;
  const Y = graded([-yext, -yf, -b - t1, -b + t2, -a - t2, -a + t1, 0, a - t1, a + t2, b - t2, b + t1, yf, yext],
    [-b, -a, a, b], p.hmin, p.hmax, 1.3);
  const Z = graded([0, h, ztop], [h], p.hmin, p.hmax, 1.3);
  const Y2: number[] = [];
  for (let i = npml; i >= 1; i--) Y2.push(Y[0] - i * p.hmax);
  Y2.push(...Y);
  for (let i = 1; i <= npml; i++) Y2.push(Y[Y.length - 1] + i * p.hmax);
  for (let i = 1; i <= npml; i++) Z.push(ztop + i * p.hmax);
  const X: number[] = [];
  const mode = p.ports === 'mode';
  // mode ports: 4 mm feed lines between the CPML and the reference planes,
  // monitor 2.5 mm from the source (same defaults as gcpw_line --modeports)
  const nlen = lround(p.len / p.dx), next = lround((mode ? 4.0e-3 : ext) / p.dx);
  for (let i = -next - npml; i <= nlen + next + npml; i++) X.push(i * p.dx);
  const nx = X.length - 1, ny = Y2.length - 1, nz = Z.length - 1;
  const ip1 = next + npml, ip2 = ip1 + nlen;
  const j0 = nearest(Y2, 0), kh = nearest(Z, h);

  const sim = cem.createSim(X, Y2, Z);
  const xlo = X[0] - 1, xhi = X[nx] + 1, ylo = Y2[0] - 1, yhi = Y2[ny] + 1;
  sim.setBoundary(Face.XLo, Boundary.CPML, npml);
  sim.setBoundary(Face.XHi, Boundary.CPML, npml);
  sim.setBoundary(Face.YLo, Boundary.CPML, npml);
  sim.setBoundary(Face.YHi, Boundary.CPML, npml);
  sim.setBoundary(Face.ZHi, Boundary.CPML, npml);
  if (p.lossy) {
    const m = cem.fitDjordjevicSarkar(p.er, 0.0037, 10e9, p.df, p.fmax, 8);
    sim.setDebyePoles(m.fp);
    sim.materialBoxDebye(xlo, xhi, ylo, yhi, 0, h, m, 0);
    sim.setConductorLoss(5.8e7, p.df, p.fmax);
  } else {
    sim.materialBox(xlo, xhi, ylo, yhi, 0, h, p.er, 0);
  }
  sim.pecBox(xlo, xhi, ylo, -yf, 0, h);
  sim.pecBox(xlo, xhi, yf, yhi, 0, h);
  sim.pecSheet(Axis.Z, h, xlo, xhi, ylo, -b);
  sim.pecSheet(Axis.Z, h, xlo, xhi, b, yhi);
  const tau = Math.sqrt(6.9) / (Math.PI * 1.25 * p.fmax);
  let p1: number, p2: number;
  if (!mode) {
    sim.pecSheet(Axis.Z, h, X[ip1], X[ip2], -a, a);
    p1 = sim.addLumpedPort(Axis.Z, ip1, j0, 0, kh, -1, 50);
    p2 = sim.addLumpedPort(Axis.Z, ip2, j0, 0, kh, -1, 50);
    sim.setGaussian(p1, 1.0, 4.5 * tau, tau);
  } else {
    sim.pecSheet(Axis.Z, h, xlo, xhi, -a, a);
    const jl = nearest(Y2, -(a + b) / 2), jr = nearest(Y2, (a + b) / 2);
    const mon = lround(2.5e-3 / p.dx);
    const common = { j0: 0, j1: ny, k0: 0, k1: nz, vj: j0, vk0: 0, vk1: kh, vsign: -1 as const,
      ja: jl, jb: jr - 1, ka: kh - 1, kb: kh, f0: 0.5 * p.fmax, neff0: 1.65 };
    const s1 = npml + 2, s2 = nx - npml - 2;
    p1 = sim.addModePort({ ...common, dir: 1, iSrc: s1, iMon: s1 + mon, xRef: X[ip1] });
    p2 = sim.addModePort({ ...common, dir: -1, iSrc: s2, iMon: s2 - mon - 1, xRef: X[ip2] });
    sim.setModeExcitation(p1, 4.5 * tau, tau);
  }
  const nf = lround(p.fmax / p.df);
  const freqs: number[] = [];
  for (let i = 0; i < nf; i++) freqs.push((i + 1) * p.df);
  sim.setFrequencies(freqs);
  sim.finalize(0.99);
  return { sim, freqs, ports: [p1, p2], portType: p.ports, cells: [nx, ny, nz], tStop: 4.5 * tau + 6 * tau,
    mesh: { x: X, y: Y2, z: Z, kh } };
}

export interface RunProgress { steps: number; time: number; energyRatio: number }

/** Same stopping rule as cem_run(): after the source, stop when energy < decay * peak. */
export function runToDecay(setup: GcpwSetup, decay: number, maxTime: number, checkEvery: number,
                           onProgress?: (p: RunProgress) => void): void {
  const sim = setup.sim;
  const maxSteps = Math.floor(maxTime / sim.dt);
  let peak = 0;
  while (sim.steps < maxSteps) {
    const chunk = Math.min(checkEvery, maxSteps - sim.steps);
    const before = sim.steps;
    sim.step(chunk);
    if (sim.steps === before) throw new Error('solver did not advance (not finalized?)');
    const W = sim.energy();
    if (W > peak) peak = W;
    onProgress?.({ steps: sim.steps, time: sim.steps * sim.dt, energyRatio: peak > 0 ? W / peak : 1 });
    if (sim.steps * sim.dt > setup.tStop && peak > 0 && W < decay * peak) break;
  }
}

export function gcpwResults(setup: GcpwSetup, cem?: Courant):
  { f: number[]; s11: Complex[]; s21: Complex[]; z0?: number[]; s11_50?: Complex[]; s21_50?: Complex[] } {
  const [p1, p2] = setup.ports;
  const sim = setup.sim;
  if (setup.portType === 'lumped') {
    const { s11, s21 } = sParams(sim.portV(p1), sim.portI(p1), sim.portV(p2), sim.portI(p2), 50);
    return { f: setup.freqs, s11, s21 };
  }
  // mode ports: S normalised to the line's own Z0(f) (renormalisation to 50 ohm is a later step)
  const w1 = sim.modePortWaves(p1), w2 = sim.modePortWaves(p2);
  const div = (x: Complex, y: Complex): Complex => {
    const d = y.re * y.re + y.im * y.im;
    return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d };
  };
  const s11 = w1.b.map((b, i) => div(b, w1.a[i]));
  const s21 = w2.b.map((b, i) => div(b, w1.a[i]));
  const out: ReturnType<typeof gcpwResults> = { f: setup.freqs, s11, s21, z0: Array.from(w1.z0) };
  if (cem) {
    // 50-ohm view; the line is symmetric, so S22 = S11 and S12 = S21
    out.s11_50 = []; out.s21_50 = [];
    s11.forEach((a, i) => {
      const r = cem.renormalize2([a, s21[i], s21[i], a], w1.z0[i], w1.z0[i], 50);
      out.s11_50!.push(r[0]); out.s21_50!.push(r[2]);
    });
  }
  return out;
}
