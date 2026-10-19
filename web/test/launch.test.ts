/**
 * A realistic coupon (test/fixtures/launch.kicad_pcb: edge launches, a jog with
 * 45-degree segments and an arc, top pours, stitching vias, solder mask, rounded
 * outline) through the whole board path: port planes 3.5 mm in, past the launches,
 * a +-3 mm band around the line, copper thickness and mask from the stackup.
 * Checks the run decays (no open-edged parallel-plate ringing at the band edges),
 * the port impedance, and a sensible S21.
 * Run: node web/test/launch.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant, type Complex } from '../src/solver/courant.ts';
import { importKicad } from '../src/import/kicad.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { portNets } from '../src/import/load.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const b = importKicad(readFileSync(join(here, 'fixtures', 'launch.kicad_pcb'), 'utf8'));
check(b.layers.some((l) => l.mask) && b.vias.length > 90, `mask layer and ${b.vias.length} vias (incl. 8 TH pins)`);
check(portNets(b).length === 0 && portNets(b, 3.5e-3, 3.5e-3).join() === 'RF', 'RF found only once the planes pass the launch pads');
const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));
const t0 = Date.now();
const bs = buildBoardSim(cem, b, { fmax: 40e9, df: 1e9, hmin: 60e-6, thickCopper: true, yWindow: 3e-3,
  ports: [{ net: 'RF', edge: 'xlo', offset: 3.5e-3 }, { net: 'RF', edge: 'xhi', offset: 3.5e-3 }] });
let peak = 0, W = 0;
for (;;) {
  bs.sim.step(500); W = bs.sim.energy(); peak = Math.max(peak, W);
  if ((bs.sim.steps * bs.sim.dt > bs.tStop && W < 1e-5 * peak) || bs.sim.steps * bs.sim.dt > 2e-9) break;
}
const w1 = bs.sim.modePortWaves(bs.ports[0]), w2 = bs.sim.modePortWaves(bs.ports[1]);
const div = (x: Complex, y: Complex) => { const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d }; };
const db = (c: Complex) => 20 * Math.log10(Math.hypot(c.re, c.im));
const s11 = w1.b.map((v, i) => db(div(v, w1.a[i]))), s21 = w2.b.map((v, i) => db(div(v, w1.a[i])));
console.log(`mesh ${bs.cells.join(' x ')}, ${bs.sim.steps} steps (${(bs.sim.steps * bs.sim.dt * 1e12).toFixed(0)} ps) in ${((Date.now() - t0) / 1000).toFixed(0)} s, ` +
  `energy ${(10 * Math.log10(W / peak)).toFixed(0)} dB; Z0 ${w1.z0[1].toFixed(2)} ohm; |S11| max ${Math.max(...s11).toFixed(1)} dB; ` +
  `|S21| at 5/20/40 GHz ${[4, 19, 39].map((i) => s21[i].toFixed(2)).join(' / ')} dB`);
check(W < 1e-5 * peak, 'field energy decays below -50 dB before the 2 ns limit');
check(w1.z0[1] > 44 && w1.z0[1] < 49, 'port Z0 of the GCPW with 35 um copper and 10 um mask in 44-49 ohm');
check(s21.slice(0, 20).every((v) => v > -1.5), '|S21| above -1.5 dB to 20 GHz (jog and stepped-back pour are mild discontinuities)');
console.log(fail ? 'launch test FAILED' : 'launch test passed');
process.exit(fail);
