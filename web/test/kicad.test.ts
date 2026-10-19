/**
 * KiCad importer on constructs real boards use: arc tracks, rounded-rect / oval /
 * custom / through-hole pads in a rotated footprint, and an outline made of
 * arcs, a circle-free rounded rectangle. Geometry only, no simulation.
 * Run: node web/test/kicad.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { copperAt, copperRow, importKicad, inRow, type Board } from '../src/import/kicad.ts';

let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };
const mm = (x: number, y: number): [number, number] => [x * 1e-3, -y * 1e-3]; // KiCad (y down) -> model

const b = importKicad(`(kicad_pcb (version 20240108) (generator "pcbnew") (generator_version "8.0")
  (general (thickness 0.324))
  (layers (0 "F.Cu" signal) (31 "B.Cu" signal) (44 "Edge.Cuts" user))
  (setup (stackup
    (layer "F.Cu" (type "copper") (thickness 0.035))
    (layer "dielectric 1" (type "core") (thickness 0.254) (material "RO4350B") (epsilon_r 3.66) (loss_tangent 0.0037))
    (layer "B.Cu" (type "copper") (thickness 0.035))))
  (net 0 "") (net 1 "RF") (net 2 "GND")
  (footprint "Test:Launch" (layer "F.Cu") (at 10 10 90)
    (pad "1" smd roundrect (at 0 0 90) (size 2 1) (layers "F.Cu" "F.Mask") (roundrect_rratio 0.25) (net 1 "RF"))
    (pad "2" thru_hole circle (at 3 0 90) (size 1.2 1.2) (drill 0.6) (layers "*.Cu" "*.Mask") (net 2 "GND"))
    (pad "3" smd oval (at -3 0 90) (size 2 1) (layers "F.Cu") (net 2 "GND"))
    (pad "4" smd custom (at 0 4 90) (size 0.4 0.4) (layers "F.Cu") (net 2 "GND")
      (primitives (gr_poly (pts (xy -1 -0.2) (xy 1 -0.2) (xy 1 0.2) (xy -1 0.2)) (width 0) (fill yes)))))
  (arc (start 20 10) (mid 25 15) (end 20 20) (width 0.4) (layer "F.Cu") (net 1))
  (gr_line (start 1 0) (end 29 0) (layer "Edge.Cuts"))
  (gr_arc (start 29 0) (mid 29.7071 0.2929) (end 30 1) (layer "Edge.Cuts"))
  (gr_line (start 30 1) (end 30 29) (layer "Edge.Cuts"))
  (gr_line (start 0 1) (end 0 30) (layer "Edge.Cuts"))
)`);

check(Math.abs(b.outline.min[0]) < 1e-12 && Math.abs(b.outline.max[0] - 30e-3) < 1e-9 && Math.abs(b.outline.max[1] - 0) < 1e-9 &&
  Math.abs(b.outline.min[1] + 30e-3) < 1e-9, 'outline with an arc corner: 30 x 30 mm');
// arc track: centre (20, 15), radius 5, through (25, 15) on the right
const arcTracks = b.tracks.filter((t) => t.net === 1);
check(arcTracks.length >= 18, `arc track split into ${arcTracks.length} segments`);
check(copperAt(b, 'F.Cu', mm(25, 15)) && copperAt(b, 'F.Cu', mm(25.15, 15)) && !copperAt(b, 'F.Cu', mm(25.3, 15)) && !copperAt(b, 'F.Cu', mm(15, 15)),
  'arc track passes (25, 15) with 0.4 mm width, not the far side');
// footprint at (10,10) rotated 90: pad offsets (x, y) -> (y, -x) in KiCad's y-down frame
// roundrect 2 x 1 rotated 90 at (10,10): 1 wide in x, 2 tall in y, corner radius 0.25
check(copperAt(b, 'F.Cu', mm(10, 10.95)) && copperAt(b, 'F.Cu', mm(10.45, 10)) && !copperAt(b, 'F.Cu', mm(10.47, 10.97)) && !copperAt(b, 'F.Cu', mm(10.6, 10)),
  'rotated round-rect pad: 1 x 2 mm with rounded corners');
// through-hole pad 3 mm along local x -> (10, 7): a barrel and copper on both layers
const th = b.vias.find((v) => Math.abs(v.at[0] - 10e-3) < 1e-9 && Math.abs(v.at[1] + 7e-3) < 1e-9);
check(!!th && Math.abs(th.drill - 0.6e-3) < 1e-12 && th.net === 2, 'plated through-hole pad becomes a barrel (0.6 mm)');
check(copperAt(b, 'F.Cu', mm(10.5, 7)) && copperAt(b, 'B.Cu', mm(10.5, 7)), 'through-hole pad copper on both layers');
// oval 2 x 1 at local (-3, 0) -> (10, 13), rotated: 1 wide in x, 2 tall in y, a stadium
check(copperAt(b, 'F.Cu', mm(10, 13.95)) && !copperAt(b, 'F.Cu', mm(10.45, 13.95)), 'oval pad is a stadium');
// custom pad at local (0, 4) -> (14, 10): primitive 2 x 0.4 rotated -> 0.4 wide in x, 2 tall in y
check(copperAt(b, 'F.Cu', mm(14, 10.9)) && !copperAt(b, 'F.Cu', mm(14.3, 10.9)), 'custom pad primitive');
// the row rasteriser must agree with copperAt everywhere, including exactly on copper edges
const here = dirname(fileURLToPath(import.meta.url));
for (const [name, bd] of [['features', b], ['gcpw fixture', importKicad(readFileSync(join(here, 'fixtures', 'gcpw.kicad_pcb'), 'utf8'))]] as [string, Board][]) {
  let bad = 0, n = 0;
  const xs: number[] = [], ys: number[] = [];
  const [x0, y0] = bd.outline.min, [x1, y1] = bd.outline.max;
  for (let i = 0; i <= 300; i++) xs.push(x0 + ((x1 - x0) * i) / 300);
  for (let j = 0; j <= 200; j++) ys.push(y0 + ((y1 - y0) * j) / 200);
  // and every copper vertex coordinate, to hit edges exactly
  for (const q of bd.polys) for (const p of q.pts) { xs.push(p[0]); ys.push(p[1]); }
  for (const t of bd.tracks) { xs.push(t.a[0] - t.width / 2, t.a[0] + t.width / 2); ys.push(t.a[1] - t.width / 2, t.a[1] + t.width / 2); }
  for (const layer of bd.copperLayers) for (const y of ys) {
    const r = copperRow(bd, layer, y);
    for (const x of xs) { n++; if (inRow(r, x) !== copperAt(bd, layer, [x, y])) bad++; }
  }
  check(bad === 0, `row rasteriser = copperAt on the ${name} (${n} samples, ${bad} differ)`);
}
console.log(fail ? 'kicad test FAILED' : 'kicad test passed');
process.exit(fail);
