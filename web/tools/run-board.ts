/**
 * Command-line board run (Node 22+, no browser): import a KiCad board or a Gerber
 * set, simulate the two-port net between the left and right port planes with the
 * f64 WebAssembly core, and write a Touchstone file stamped with the version.
 *
 *   node web/tools/run-board.ts board.kicad_pcb --net RF --out board.s2p
 *   node web/tools/run-board.ts top.gbr bottom.gbr edge.gm1 holes.drl --stackup stack.json --net Net-1 \
 *        --planes 1.5,1.5 --hmin 20 --fmax 40 --df 0.25 --lossy --full2 --out board.s2p
 *
 * Options: --net NAME (default: the first net crossing both planes), --planes LO,HI
 * (mm into the board), --hmin UM (default 20), --fmax GHZ (40), --df GHZ (0.5),
 * --lossy, --sheets (zero-thickness copper instead of stackup thickness), --full2
 * (excite port 2 too), --stackup FILE (JSON array of stackup rows, top to bottom,
 * SI units, as the app's stackup editor holds them), --decay (1e-5), --tmax NS (5), --window MM (simulate only
 * a band of +-MM around the port net).
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, basename } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant, type Complex } from '../src/solver/courant.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { loadBoard, type BoardInput } from '../src/import/load.ts';
import type { StackupRow } from '../src/import/gerber.ts';
import { toTouchstone } from '../src/app/touchstone.ts';

const here = dirname(fileURLToPath(import.meta.url));
const args = process.argv.slice(2);
const files: string[] = [];
const opt: Record<string, string | boolean> = {};
for (let i = 0; i < args.length; i++) {
  const a = args[i];
  if (a.startsWith('--')) {
    const k = a.slice(2);
    if (['lossy', 'sheets', 'full2'].includes(k)) opt[k] = true; else opt[k] = args[++i];
  } else files.push(a);
}
if (!files.length) { console.error('usage: node web/tools/run-board.ts <board files> [options] (see the header)'); process.exit(2); }

const input: BoardInput = files.length === 1 && files[0].endsWith('.kicad_pcb')
  ? { kind: 'kicad', name: basename(files[0]), text: readFileSync(files[0], 'utf8') }
  : { kind: 'gerber', name: basename(files[0]), files: files.map((f) => ({ name: basename(f), text: readFileSync(f, 'utf8') })) };
const stackup = opt.stackup ? JSON.parse(readFileSync(String(opt.stackup), 'utf8')) as StackupRow[] : undefined;
const planes = (opt.planes ? String(opt.planes).split(',').map((v) => Number(v) * 1e-3) : [0, 0]) as [number, number];
const lb = loadBoard(input, stackup, planes);
for (const w of lb.warnings) console.warn('warning:', w);
const net = opt.net ? String(opt.net) : lb.portNets[0];
if (!net) { console.error('no net crosses both port planes; nets crossing them:', lb.portNets); process.exit(1); }
const fmax = Number(opt.fmax ?? 40) * 1e9, df = Number(opt.df ?? 0.5) * 1e9, hmin = Number(opt.hmin ?? 20) * 1e-6;
const decay = Number(opt.decay ?? 1e-5), tmax = Number(opt.tmax ?? 5) * 1e-9;
const version = readFileSync(join(here, '..', '..', 'VERSION'), 'utf8').trim();
(globalThis as unknown as { __COURANT_VERSION__: string }).__COURANT_VERSION__ = version;

const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));
const div = (x: Complex, y: Complex) => { const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d }; };
const t0 = Date.now();
const waves = [];
let freqs: number[] = [], steps = 0;
for (let run = 0; run < (opt.full2 ? 2 : 1); run++) {
  const bs = buildBoardSim(cem, lb.board, { fmax, df, hmin, lossy: !!opt.lossy, thickCopper: !opt.sheets, excite: run, yWindow: opt.window ? Number(opt.window) * 1e-3 : undefined,
    ports: [{ net, edge: 'xlo', offset: planes[0] }, { net, edge: 'xhi', offset: planes[1] }] });
  console.log(`run ${run + 1}: mesh ${bs.cells.join(' x ')} = ${(bs.cells[0] * bs.cells[1] * bs.cells[2] / 1e6).toFixed(2)} M cells`);
  let peak = 0;
  for (;;) {
    bs.sim.step(500);
    const W = bs.sim.energy(); peak = Math.max(peak, W);
    if (bs.sim.steps % 5000 === 0) console.log(`  step ${bs.sim.steps}, energy ${(10 * Math.log10(W / peak)).toFixed(1)} dB`);
    if ((bs.sim.steps * bs.sim.dt > bs.tStop && W < decay * peak) || bs.sim.steps * bs.sim.dt > tmax) {
      if (bs.sim.steps * bs.sim.dt > tmax) console.warn(`warning: stopped at the ${tmax * 1e9} ns limit with energy ${(10 * Math.log10(W / peak)).toFixed(0)} dB below peak (slowly decaying resonances; see --tmax, vias, losses, --window)`);
      break;
    }
  }
  waves.push({ w1: bs.sim.modePortWaves(bs.ports[0]), w2: bs.sim.modePortWaves(bs.ports[1]) });
  freqs = bs.freqs; steps += bs.sim.steps;
  bs.sim.destroy();
}
const { w1, w2 } = waves[0];
const s11 = w1.b.map((b, i) => div(b, w1.a[i])), s21 = w2.b.map((b, i) => div(b, w1.a[i]));
const s12 = waves[1]?.w1.b.map((b, i) => div(b, waves[1].w2.a[i])), s22 = waves[1]?.w2.b.map((b, i) => div(b, waves[1].w2.a[i]));
const r50 = s11.map((a, i) => cem.renormalize2([a, s12?.[i] ?? s21[i], s21[i], s22?.[i] ?? a], w1.z0[i], w2.z0[i], 50));
const result = {
  f: freqs, s11, s21, s12, s22, z0: Array.from(w1.z0), steps, wall: (Date.now() - t0) / 1000,
  s11_50: r50.map((r) => r[0]), s12_50: s12 ? r50.map((r) => r[1]) : undefined, s21_50: r50.map((r) => r[2]), s22_50: s22 ? r50.map((r) => r[3]) : undefined,
};
const out = String(opt.out ?? `${input.name.replace(/\.[^.]*$/, '')}.s2p`);
writeFileSync(out, toTouchstone(result, `${input.name}, net ${net}, planes ${planes.map((p) => p * 1e3).join('/')} mm`));
console.log(`wrote ${out}: ${steps} steps in ${result.wall.toFixed(0)} s; Z0 ${w1.z0[0].toFixed(2)}-${w1.z0[w1.z0.length - 1].toFixed(2)} ohm`);
