/**
 * The app's built-in GCPW (scenarios/gcpwboard.ts) through the board path: via
 * fences and 1 oz copper are in the model, port Z0 matches the open Laplace
 * reference of tools/gcpw_ref.py with 35 um copper (46.761 ohm) and with
 * zero-thickness sheets (49.92 ohm), and a full solve gives a matched line.
 * Run: node web/test/builtin.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant } from '../src/solver/courant.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { GCPW_NET, gcpwBoard } from '../src/scenarios/gcpwboard.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));
const ports = [{ net: GCPW_NET, edge: 'xlo' as const }, { net: GCPW_NET, edge: 'xhi' as const }];

const b = gcpwBoard({ len: 10e-3, er: 3.66 });
const rows = [...new Set(b.vias.map((v) => v.at[1].toFixed(6)))];
check(b.vias.length === 24 && rows.length === 2, `two via rows along the line: ${b.vias.length} vias, rows at y = ${rows.join(', ')} m`);
check(b.vias.every((v) => Math.abs(v.at[1]) - v.drill / 2 - 450e-6 > -1e-12), 'via barrels start 300 um beyond the gaps');
check(b.layers.filter((l) => l.type === 'copper').every((l) => l.thickness === 35e-6), 'copper layers are 1 oz (35 um)');

for (const [thick, ref] of [[false, 49.92], [true, 46.761]] as const) {
  const bs = buildBoardSim(cem, b, { fmax: 40e9, df: 1e9, hmin: 20e-6, thickCopper: thick, ports });
  const z0 = bs.sim.modePortWaves(bs.ports[0]).z0[1];
  const err = (z0 - ref) / ref;
  check(Math.abs(err) < 0.01, `${thick ? '35 um copper' : 'sheets'}: Z0 ${z0.toFixed(3)} vs Laplace ${ref} (${(err * 100).toFixed(2)}%)`);
  bs.sim.destroy();
}

// full solve, coarse mesh as in import.test.ts, shorter line
const t0 = Date.now();
const bs = buildBoardSim(cem, gcpwBoard({ len: 6e-3, er: 3.66 }), { fmax: 40e9, df: 1e9, hmin: 40e-6, thickCopper: true, ports });
let peak = 0;
for (;;) {
  bs.sim.step(500);
  const W = bs.sim.energy();
  peak = Math.max(peak, W);
  if (bs.sim.steps * bs.sim.dt > bs.tStop && W < 1e-5 * peak) break;
  if (bs.sim.steps * bs.sim.dt > 3e-9) break;
}
const w1 = bs.sim.modePortWaves(bs.ports[0]), w2 = bs.sim.modePortWaves(bs.ports[1]);
const mag = (c: { re: number; im: number }) => Math.hypot(c.re, c.im);
let s11max = 0, s21dev = 0;
bs.freqs.forEach((f, i) => {
  if (f < 1e9) return;
  const a = mag(w1.a[i]);
  s11max = Math.max(s11max, mag(w1.b[i]) / a);
  s21dev = Math.max(s21dev, Math.abs(20 * Math.log10(mag(w2.b[i]) / a)));
});
console.log(`mesh ${bs.cells.join(' x ')}, ${bs.sim.steps} steps in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
check(20 * Math.log10(s11max) < -20, `matched: |S11| max ${(20 * Math.log10(s11max)).toFixed(1)} dB < -20 dB`);
check(s21dev < 0.3, `||S21| - 1| max ${s21dev.toFixed(4)} dB < 0.3 dB`);
bs.sim.destroy();
console.log(fail ? 'builtin test FAILED' : 'builtin test passed');
process.exit(fail);
