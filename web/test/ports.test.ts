/**
 * Port planes inside the board: the microstrip fixture gets a wide pad near each
 * end, 0.2–0.8 mm from the edge (a step discontinuity, like a connector launch). With the port planes
 * moved 1 mm into the board, past the pads, the pads are cut away and the result
 * must match the plain line with the same planes: matched, and the same S21.
 * Without the planes moved, the pads must show up as a mismatch.
 * Run: node web/test/ports.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant, type Complex } from '../src/solver/courant.ts';
import { importKicad, type Board, type Pt } from '../src/import/kicad.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { portNets } from '../src/import/load.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));

const plain = importKicad(readFileSync(join(here, 'fixtures', 'msl.kicad_pcb'), 'utf8'));
const t = plain.tracks[0], yc = t.a[1];
const [x0] = plain.outline.min, [x1] = plain.outline.max;
const pad = (xa: number, xb: number): Pt[] => [[xa, yc - 0.75e-3], [xb, yc - 0.75e-3], [xb, yc + 0.75e-3], [xa, yc + 0.75e-3]];
const padded: Board = { ...plain, polys: [...plain.polys,
  { pts: pad(x0 + 0.2e-3, x0 + 0.8e-3), layer: t.layer, net: t.net }, { pts: pad(x1 - 0.8e-3, x1 - 0.2e-3), layer: t.layer, net: t.net }] };
check(portNets(padded, 1e-3, 1e-3).join() === 'RF', 'port net found at planes 1 mm inside the board');

function run(b: Board, off: number, excite = 0) {
  const bs = buildBoardSim(cem, b, { fmax: 40e9, df: 1e9, hmin: 40e-6, excite,
    ports: [{ net: 'RF', edge: 'xlo', offset: off }, { net: 'RF', edge: 'xhi', offset: off }] });
  let peak = 0;
  for (;;) {
    bs.sim.step(500);
    const W = bs.sim.energy(); peak = Math.max(peak, W);
    if ((bs.sim.steps * bs.sim.dt > bs.tStop && W < 1e-5 * peak) || bs.sim.steps * bs.sim.dt > 3e-9) break;
  }
  const w1 = bs.sim.modePortWaves(bs.ports[0]), w2 = bs.sim.modePortWaves(bs.ports[1]);
  const div = (x: Complex, y: Complex) => { const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d }; };
  const a = excite ? w2.a : w1.a;
  const out = { f: bs.freqs, s11: w1.b.map((v, i) => div(v, a[i])), s21: w2.b.map((v, i) => div(v, a[i])), cells: bs.cells,
    w: { a1: w1.a, b1: w1.b, a2: w2.a, b2: w2.b } };
  bs.sim.destroy();
  return out;
}
const db = (c: Complex) => 20 * Math.log10(Math.hypot(c.re, c.im));
const t0 = Date.now();
const A = run(plain, 1e-3), B = run(padded, 1e-3), C = run(padded, 0);
console.log(`3 runs in ${((Date.now() - t0) / 1000).toFixed(0)} s; meshes ${A.cells.join('x')}, ${B.cells.join('x')}, ${C.cells.join('x')}`);
let s11B = -999, s11C = -999, dmag = 0, dph = 0;
A.f.forEach((f, i) => {
  if (f < 1e9) return;
  s11B = Math.max(s11B, db(B.s11[i])); s11C = Math.max(s11C, db(C.s11[i]));
  const a = A.s21[i], b = B.s21[i];
  dmag = Math.max(dmag, Math.abs(db(a) - db(b)));
  let d = Math.atan2(b.im, b.re) - Math.atan2(a.im, a.re);
  d = Math.atan2(Math.sin(d), Math.cos(d));
  dph = Math.max(dph, Math.abs(d) * 180 / Math.PI);
});
console.log(`pads cut away: |S11| max ${s11B.toFixed(1)} dB, S21 vs plain line ${dmag.toFixed(4)} dB, ${dph.toFixed(3)} deg; pads kept: |S11| max ${s11C.toFixed(1)} dB`);
check(s11B < -35, 'planes past the pads: matched');
check(dmag < 0.02 && dph < 0.5, 'planes past the pads: same S21 as the plain line with the same planes');
check(s11C > -25, 'planes at the board edges: the launch pads show up as a mismatch');
// reciprocity on an asymmetric board: one pad, near the left edge only
const oneSided: Board = { ...plain, polys: [...plain.polys, { pts: pad(x0 + 0.2e-3, x0 + 0.8e-3), layer: t.layer, net: t.net }] };
const P1 = run(oneSided, 0, 0), P2 = run(oneSided, 0, 1);
// Two runs give S = B A^-1 with A = [a1 a1'; a2 a2'], B = [b1 b1'; b2 b2'], which also
// removes the small waves reflected back by the far port's CPML termination.
const div = (x: Complex, y: Complex) => { const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d }; };
const mul = (x: Complex, y: Complex) => ({ re: x.re * y.re - x.im * y.im, im: x.re * y.im + x.im * y.re });
const sub = (x: Complex, y: Complex) => ({ re: x.re - y.re, im: x.im - y.im });
let recipRatio = 0, recipFull = 0, asym = 0;
P1.f.forEach((f, i) => {
  if (f < 1e9) return;
  // ratio extraction: S21 = b2/a1 (run 1), S12 = b1'/a2' (run 2)
  recipRatio = Math.max(recipRatio, Math.hypot(P1.s21[i].re - P2.s11[i].re, P1.s21[i].im - P2.s11[i].im));
  const A = [[P1.w.a1[i], P2.w.a1[i]], [P1.w.a2[i], P2.w.a2[i]]], B = [[P1.w.b1[i], P2.w.b1[i]], [P1.w.b2[i], P2.w.b2[i]]];
  const det = sub(mul(A[0][0], A[1][1]), mul(A[0][1], A[1][0]));
  const inv = [[div(A[1][1], det), div({ re: -A[0][1].re, im: -A[0][1].im }, det)], [div({ re: -A[1][0].re, im: -A[1][0].im }, det), div(A[0][0], det)]];
  const S = [0, 1].map((r) => [0, 1].map((c) => { const p = mul(B[r][0], inv[0][c]), q = mul(B[r][1], inv[1][c]); return { re: p.re + q.re, im: p.im + q.im }; }));
  recipFull = Math.max(recipFull, Math.hypot(S[0][1].re - S[1][0].re, S[0][1].im - S[1][0].im));
  let dph = Math.atan2(S[0][0].im, S[0][0].re) - Math.atan2(S[1][1].im, S[1][1].re);
  dph = Math.abs(Math.atan2(Math.sin(dph), Math.cos(dph))) * 180 / Math.PI;
  if (db(S[0][0]) > -20) asym = Math.max(asym, dph);
});
console.log(`asymmetric board: |S21 - S12| max ${recipRatio.toExponential(2)} from wave ratios, ${recipFull.toExponential(2)} from S = B A^-1 (so not the far-end termination); ` +
  `S11 and S22 phases differ by up to ${asym.toFixed(0)} deg`);
// Not exact: the pad radiates (open microstrip), and space and surface waves reach the
// monitor planes; the error grows with frequency (1e-4 at 2 GHz, ~6e-3 at 40 GHz, where
// |S21| is 0.83) and shrinks with longer feeds. A shielded line would be reciprocal to
// the mode-port accuracy.
check(recipRatio < 1e-2, 'reciprocity: S12 = S21 on an asymmetric, radiating board (two independent runs), within 1e-2');
check(asym > 10, 'S11 and S22 differ in phase on the asymmetric board');
console.log(fail ? 'ports test FAILED' : 'ports test passed');
process.exit(fail);
