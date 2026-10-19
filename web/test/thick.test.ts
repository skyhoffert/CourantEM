/**
 * Thick copper on the imported GCPW fixture: port Z0 (2D eigenmode on the FDTD
 * mesh, 2 GHz) vs the open Laplace reference of tools/gcpw_ref.py, with 35 um
 * copper (46.761 ohm) and with zero-thickness sheets (49.92 ohm), then with a
 * 20 um, er 3.5 solder-mask slab added in the stackup (45.625 / 48.626 ohm).
 * Run: node web/test/thick.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant } from '../src/solver/courant.ts';
import { importKicad } from '../src/import/kicad.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { applyStackup, stackupOf } from '../src/import/load.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));
const b = importKicad(readFileSync(join(here, 'fixtures', 'gcpw.kicad_pcb'), 'utf8'));
const masked = applyStackup(b, [{ name: 'F.Mask', type: 'dielectric', mask: true, thickness: 20e-6, epsR: 3.5, lossTangent: 0.02 }, ...stackupOf(b)]);
for (const [thick, ref, bd] of [[false, 49.92, b], [true, 46.761, b], [false, 48.626, masked], [true, 45.625, masked]] as const) {
  const bs = buildBoardSim(cem, bd, { fmax: 40e9, df: 1e9, hmin: 20e-6, thickCopper: thick,
    ports: [{ net: 'RF', edge: 'xlo' }, { net: 'RF', edge: 'xhi' }] });
  const z0 = bs.sim.modePortWaves(bs.ports[0]).z0[1];
  const err = (z0 - ref) / ref;
  check(Math.abs(err) < 0.01, `${thick ? '35 um copper' : 'sheets'}${bd === masked ? ' + mask' : ''}: Z0 ${z0.toFixed(3)} vs Laplace ${ref} (${(err * 100).toFixed(2)}%)`);
  bs.sim.destroy();
}
console.log(fail ? 'thick test FAILED' : 'thick test passed');
process.exit(fail);
