/**
 * Gerber + Excellon import tests.
 *  1. Parser units: formats, units, standard apertures, KiCad's RoundRect macro,
 *     arcs, regions, Excellon LZ/TZ.
 *  2. Round trip: the KiCad GCPW fixture is written out as KiCad-style Gerbers
 *     (copper as tracks, regions and via-land flashes; Edge.Cuts profile; PTH
 *     drill file), imported back with the same stackup, and must give the same
 *     board model: copper on a dense grid, mesh edges, vias, the port net, and
 *     a bit-identical first 300 FDTD steps.
 * Run: node web/test/gerber.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Courant } from '../src/solver/courant.ts';
import { copperAt, copperEdges, importKicad, type Board } from '../src/import/kicad.ts';
import { classify, importGerber, parseExcellon, parseGerber, type StackupRow } from '../src/import/gerber.ts';
import { buildBoardSim } from '../src/import/boardsim.ts';
import { toGerbers } from './fixtures/gerberwriter.ts';

const here = dirname(fileURLToPath(import.meta.url));
let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const close = (a: number, b: number, tol = 1e-12) => Math.abs(a - b) <= tol;

// ---------------------------------------------------------------- 1. parser units
{
  const g = parseGerber(`%FSLAX46Y46*%
%MOMM*%
%TF.FileFunction,Copper,L1,Top*%
%ADD10C,0.250000*%
%ADD11R,1.000000X0.500000*%
%ADD12O,2.000000X1.000000*%
%AMRoundRect*
0 Rectangle with rounded corners*
0 $1 Rounding radius*
4,1,4,$2,$3,$4,$5,$6,$7,$8,$9,$2,$3,0*
1,1,$1+$1,$2,$3*
1,1,$1+$1,$4,$5*
1,1,$1+$1,$6,$7*
1,1,$1+$1,$8,$9*
20,1,$1+$1,$2,$3,$4,$5,0*
20,1,$1+$1,$4,$5,$6,$7,0*
20,1,$1+$1,$6,$7,$8,$9,0*
20,1,$1+$1,$8,$9,$2,$3,0*%
%ADD13RoundRect,0.100000X-0.400000X-0.300000X0.400000X-0.300000X0.400000X0.300000X-0.400000X0.300000X0*%
G04 comment*
%TO.N,SIG*%
D10*
X0Y0D02*
X5000000Y0D01*
%TD*%
D11*
X10000000Y2000000D03*
D12*
X20000000Y0D03*
D13*
X30000000Y0D03*
G75*
G36*
X0Y10000000D02*
G01X2000000Y10000000D01*
G03X0Y12000000I-2000000J0D01*
G01X0Y10000000D01*
G37*
M02*
`);
  check(g.fileFunction?.join(',') === 'Copper,L1,Top', 'X2 file function');
  check(g.tracks.length === 1 && close(g.tracks[0].b[0], 5e-3) && close(g.tracks[0].width, 0.25e-3) && g.tracks[0].net === 'SIG',
    'D01 draw with C aperture, mm, leading-zero omission, TO.N net');
  const rect = g.polys.find((p) => p.flash && close(p.flash.at[0], 10e-3))!;
  const xs = rect.pts.map((q) => q[0]), ys = rect.pts.map((q) => q[1]);
  check(close(Math.min(...xs), 9.5e-3, 1e-12) && close(Math.max(...ys), 2.25e-3, 1e-12), 'R aperture flash 1.0 x 0.5 mm');
  const obl = g.polys.find((p) => p.flash && close(p.flash.at[0], 20e-3))!;
  check(close(Math.max(...obl.pts.map((q) => q[0])), 21e-3, 1e-9) && close(Math.max(...obl.pts.map((q) => q[1])), 0.5e-3, 2e-6), 'O aperture flash');
  const rr = g.polys.filter((p) => p.flash && close(p.flash.at[0], 30e-3));
  const inRR = (x: number, y: number) => rr.some((p) => copperAt({ tracks: [], polys: [{ pts: p.pts, layer: 'L', net: 0 }] } as unknown as Board, 'L', [x, y]));
  check(rr.length === 9 && inRR(30e-3 + 0.49e-3, 0) && inRR(30e-3, 0.39e-3) && !inRR(30e-3 + 0.495e-3, 0.395e-3) && !inRR(30.6e-3, 0),
    `RoundRect macro: 9 primitives, 1.0 x 0.8 mm body, corners rounded (r 0.1 mm)`);
  const reg = g.polys.find((p) => !p.flash)!;
  const quarterArea = Math.PI * 4 / 4; // mm^2: quarter disc of radius 2 mm
  let area = 0;
  reg.pts.forEach((p, i) => { const q = reg.pts[(i + 1) % reg.pts.length]; area += (p[0] * q[1] - q[0] * p[1]) / 2; });
  check(close(Math.abs(area) * 1e6, quarterArea, 0.01), `region with G03 arc: quarter disc area ${(Math.abs(area) * 1e6).toFixed(4)} mm^2`);

  const inch = parseGerber('%FSLAX24Y24*%%MOIN*%%ADD10C,0.010*%D10*X10000Y0D02*X20000Y0D01*M02*');
  check(close(inch.tracks[0].a[0], 25.4e-3, 1e-12) && close(inch.tracks[0].width, 0.254e-3, 1e-15), 'inch, 2.4 format');
  const tz = parseGerber('%FSTAX24Y24*%%MOIN*%%ADD10C,0.010*%D10*X01Y0D02*X02Y0D01*M02*');
  check(close(tz.tracks[0].a[0], 25.4e-3, 1e-12) && close(tz.tracks[0].b[0], 50.8e-3, 1e-12), 'trailing-zero omission');

  const d1 = parseExcellon('M48\nMETRIC,TZ\nT1C0.300\n%\nT1\nX12.5Y-3.0\nM30\n');
  check(d1.length === 1 && close(d1[0].at[0], 12.5e-3) && close(d1[0].dia, 0.3e-3), 'Excellon metric with decimals');
  const d2 = parseExcellon('M48\nINCH,LZ\nT01C0.0120\n%\nT01\nX01Y-005\nM30\n');
  check(close(d2[0].at[0], 25.4e-3, 1e-12) && close(d2[0].at[1], -0.5 * 25.4e-3, 1e-12) && close(d2[0].dia, 0.3048e-3, 1e-12), 'Excellon inch LZ (00.0000)');
  const d3 = parseExcellon('M48\nINCH,TZ\nT1C0.012\n%\nT1\nX1Y-5\nM30\n');
  check(close(d3[0].at[0], 0.0001 * 25.4e-3, 1e-15) && close(d3[0].at[1], -0.0005 * 25.4e-3, 1e-15), 'Excellon inch TZ');
  check(classify('board-F_Cu.gbr', '').kind === 'copper' && classify('board.GBL', '').kind === 'copper' &&
    classify('board-Edge_Cuts.gm1', '').kind === 'profile' && classify('board-PTH.drl', 'M48').kind === 'drill' &&
    classify('x.gbr', '%TF.FileFunction,Copper,L2,Bot*%').kind === 'copper', 'file classification (X2, KiCad, Altium)');
}

// ---------------------------------------------------------------- 2. round trip
const kb = importKicad(readFileSync(join(here, 'fixtures', 'gcpw.kicad_pcb'), 'utf8'));
const stack: StackupRow[] = [...kb.layers].reverse().map((l) => ({ name: l.name, type: l.type as 'copper' | 'dielectric', thickness: l.thickness, epsR: l.epsR, lossTangent: l.lossTangent }));
const cem = await Courant.load(readFileSync(join(here, '..', '..', 'solver', 'build', 'courant.wasm')));

for (const x2 of [true, false]) {
  const gi = importGerber(toGerbers(kb, x2), stack);
  const gb = gi.board;
  console.log(`== round trip, ${x2 ? 'with' : 'without'} X2 net attributes; warnings: ${gi.warnings.length ? gi.warnings.join('; ') : 'none'}`);
  check(gb.copperLayers.join(',') === kb.copperLayers.join(','), `copper layers ${gb.copperLayers}`);
  check(close(gb.outline.min[0], kb.outline.min[0], 1e-9) && close(gb.outline.max[1], kb.outline.max[1], 1e-9), 'outline from profile');
  check(gb.vias.length === kb.vias.length && gb.vias.every((v, i) => close(v.drill, kb.vias[i].drill, 1e-9) && close(v.size, kb.vias[i].size, 1e-9)),
    `${gb.vias.length} vias with drill and land size`);
  let mism = 0, tot = 0;
  for (const layer of kb.copperLayers) for (let i = 0; i <= 200; i++) for (let j = 0; j <= 120; j++) {
    const p: [number, number] = [kb.outline.min[0] + (kb.outline.max[0] - kb.outline.min[0]) * i / 200 + 1.3e-7,
      kb.outline.min[1] + (kb.outline.max[1] - kb.outline.min[1]) * j / 120 + 1.7e-7];
    tot++; if (copperAt(kb, layer, p) !== copperAt(gb, layer, p)) mism++;
  }
  check(mism === 0, `copper raster identical (${tot} samples, ${mism} differ)`);
  for (const layer of kb.copperLayers) for (const ax of [0, 1] as const) {
    const a = copperEdges(kb, layer, ax).map((e) => e.x.toFixed(9) + e.metal).join(' ');
    const c = copperEdges(gb, layer, ax).map((e) => e.x.toFixed(9) + e.metal).join(' ');
    check(a === c, `${layer} mesh edges along ${'xy'[ax]}`);
  }
  // the RF track must be its own net running edge to edge; the grounds one net
  const rf = gb.tracks[0].net;
  const rfName = gb.nets.get(rf)!;
  check(x2 ? rfName === 'RF' : /^Net-/.test(rfName), `port net name: ${rfName}`);
  const gnd = new Set([...gb.polys.map((p) => p.net), ...gb.vias.map((v) => v.net)]);
  check(gnd.size === 1 && !gnd.has(rf), `grounds connected through vias into one net (${gnd.size})`);
  // identical simulation: mesh, coefficients and the first 300 steps
  const opts = { fmax: 40e9, df: 1e9, hmin: 40e-6, ports: [{ net: kb.nets.get(kb.tracks[0].net)!, edge: 'xlo' as const }, { net: kb.nets.get(kb.tracks[0].net)!, edge: 'xhi' as const }] };
  const A = buildBoardSim(cem, kb, opts);
  const B = buildBoardSim(cem, gb, { ...opts, ports: [{ net: rfName, edge: 'xlo' }, { net: rfName, edge: 'xhi' }] });
  check(A.cells.join() === B.cells.join() && A.sim.dt === B.sim.dt, `same mesh ${A.cells.join(' x ')} and dt`);
  A.sim.step(300); B.sim.step(300);
  let diff = 0;
  for (let c = 0; c < 6; c++) { const fa = A.sim.field(c), fb = B.sim.field(c); for (let i = 0; i < fa.length; i++) diff = Math.max(diff, Math.abs(fa[i] - fb[i])); }
  check(diff === 0, `first 300 FDTD steps bit-identical (max diff ${diff})`);
  A.sim.destroy(); B.sim.destroy();
}

console.log(fail ? 'gerber test FAILED' : 'gerber test passed');
process.exit(fail);
