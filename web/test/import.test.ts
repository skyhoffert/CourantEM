/**
 * End-to-end KiCad import test: parse test/fixtures/gcpw.kicad_pcb (GCPW with
 * via fences), mesh it, rasterise copper and vias, put mode ports where the RF
 * track leaves the board, and solve with the f64 WebAssembly core.
 * A uniform line must come out matched (|S11| small) with |S21| ~ 1 and a port
 * Z0 close to the hand-built GCPW (via fences instead of walls change it slightly).
 * Run: node web/test/import.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant } from '../src/solver/courant.ts';
import { copperEdges, importKicad } from '../src/import/kicad.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };

const board = importKicad(readFileSync(join(here, 'fixtures', 'gcpw.kicad_pcb'), 'utf8'));
check(board.copperLayers.join(',') === 'B.Cu,F.Cu', `copper layers bottom-up: ${board.copperLayers}`);
const core = board.layers.find((l) => l.type === 'dielectric')!;
check(Math.abs(core.thickness - 254e-6) < 1e-12 && core.epsR === 3.66, 'core 254 um, er 3.66');
check(board.vias.length === 24 && board.tracks.length === 1 && board.polys.length === 3, 'vias, tracks, zones parsed');
const ye = copperEdges(board, 'F.Cu', 1).map((e) => `${(e.x * 1e3).toFixed(3)}${e.metal}`);
check(ye.join(' ') === '-6.000hi -3.450lo -3.250hi -2.750lo -2.550hi 0.000lo', `F.Cu y-edges with metal sides: ${ye.join(' ')}`);

const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));

/** Hammerstad-Jensen static microstrip impedance (zero thickness). */
function msl(w: number, h: number, er: number) {
  const u = w / h;
  const a = 1 + Math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + Math.log(1 + (u / 18.1) ** 3) / 18.7;
  const b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053;
  const ee = (er + 1) / 2 + ((er - 1) / 2) * (1 + 10 / u) ** (-a * b);
  const f = 6 + (2 * Math.PI - 6) * Math.exp(-((30.666 / u) ** 0.7528));
  const z = (60 / Math.sqrt(ee)) * Math.log(f / u + Math.sqrt(1 + (2 / u) ** 2));
  return { z, ee };
}

for (const fx of ['msl', 'gcpw'] as const) {
const bd = fx === 'msl' ? importKicad(readFileSync(join(here, 'fixtures', 'msl.kicad_pcb'), 'utf8')) : board;
const t0 = Date.now();
const bs = buildBoardSim(cem, bd, {
  fmax: 40e9, df: 1e9, hmin: 40e-6, ports: [{ net: 'RF', edge: 'xlo' }, { net: 'RF', edge: 'xhi' }],
});
console.log(`== ${fx}`);
console.log(`mesh ${bs.cells.join(' x ')} = ${(bs.cells[0] * bs.cells[1] * bs.cells[2] / 1e6).toFixed(2)} M cells, built in ${Date.now() - t0} ms`);
const dt = bs.sim.dt;
let peak = 0;
for (;;) {
  bs.sim.step(500);
  const W = bs.sim.energy();
  peak = Math.max(peak, W);
  if (bs.sim.steps * dt > bs.tStop && W < 1e-5 * peak) break;
  if (bs.sim.steps * dt > 3e-9) break;
}
console.log(`ran ${bs.sim.steps} steps in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
const w1 = bs.sim.modePortWaves(bs.ports[0]), w2 = bs.sim.modePortWaves(bs.ports[1]);
const mag = (c: { re: number; im: number }) => Math.hypot(c.re, c.im);
let s11max = 0, s21dev = 0;
bs.freqs.forEach((f, i) => {
  if (f < 1e9) return;
  const a = mag(w1.a[i]);
  s11max = Math.max(s11max, mag(w1.b[i]) / a);
  s21dev = Math.max(s21dev, Math.abs(20 * Math.log10(mag(w2.b[i]) / a)));
});
const z0 = w1.z0[Math.floor(w1.z0.length / 4)];
console.log(`|S11| max ${(20 * Math.log10(s11max)).toFixed(1)} dB, ||S21|-1| max ${s21dev.toFixed(4)} dB, Z0(10 GHz) ${z0.toFixed(2)} ohm`);
if (fx === 'msl') {
  const ref = msl(0.55e-3, 254e-6, 3.66);
  const zlow = w1.z0[1];
  console.log(`Hammerstad-Jensen Z0 ${ref.z.toFixed(2)} ohm; port Z0 at 2 GHz ${zlow.toFixed(2)} ohm`);
  check(20 * Math.log10(s11max) < -35, 'imported microstrip is matched (|S11| < -35 dB)');
  check(s21dev < 0.1, '|S21| within 0.1 dB of 0 dB (lossless; radiation from the open line is small)');
  check(Math.abs(zlow - ref.z) / ref.z < 0.03, 'port Z0 within 3% of Hammerstad-Jensen');
} else {
  // a real via fence (not a wall): small mismatch to the wall-fenced feeds and some leakage
  // into the parallel-plate region outside the fence are physical
  check(20 * Math.log10(s11max) < -20, 'GCPW with via fences: |S11| < -20 dB');
  check(s21dev < 0.3, 'GCPW with via fences: |S21| within 0.3 dB');
  check(z0 > 47 && z0 < 53, 'port Z0 near 50 ohm');
}
bs.sim.destroy();
}
console.log(fail ? 'import test FAILED' : 'import test passed');
process.exit(fail);
