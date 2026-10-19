/**
 * Writes a KiCad 8 board that looks like a real RF test coupon: 30 x 16 mm, a
 * 50-ohm GCPW on 10-mil RO4350B with an edge-launch connector footprint at each
 * end (round-rect centre pad, four plated ground pins, rotated footprint), a
 * jog in the middle made of 45-degree segments and arcs, top ground pour with a
 * fractured clearance, a solid bottom plane, via fences along the line, a grid
 * of stitching vias, and a rounded board outline. Used by test/launch.test.ts.
 * Run: node web/test/fixtures/make_launch_board.ts > web/test/fixtures/launch.kicad_pcb
 */
const L = 30, W = 16;
const yc = 8, w = 0.5, s = 0.2, jog = 1.0;    // line centre, width, gap, jog height
const out: string[] = [];
out.push(`(kicad_pcb (version 20240108) (generator "pcbnew") (generator_version "8.0")
  (general (thickness 0.344))
  (paper "A4")
  (layers (0 "F.Cu" signal) (31 "B.Cu" signal) (38 "B.Mask" user) (39 "F.Mask" user) (44 "Edge.Cuts" user))
  (setup
    (stackup
      (layer "F.SilkS" (type "Top Silk Screen"))
      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "core") (thickness 0.254) (material "RO4350B") (epsilon_r 3.66) (loss_tangent 0.0037))
      (layer "B.Cu" (type "copper") (thickness 0.035))
      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))
    )
  )
  (net 0 "")
  (net 1 "RF")
  (net 2 "GND")`);
// outline with 1 mm corner arcs
const r = 1;
out.push(`  (gr_line (start ${r} 0) (end ${L - r} 0) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_arc (start ${L - r} 0) (mid ${L - r + r * Math.SQRT1_2} ${r - r * Math.SQRT1_2}) (end ${L} ${r}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_line (start ${L} ${r}) (end ${L} ${W - r}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_arc (start ${L} ${W - r}) (mid ${L - r + r * Math.SQRT1_2} ${W - r + r * Math.SQRT1_2}) (end ${L - r} ${W}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_line (start ${L - r} ${W}) (end ${r} ${W}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_arc (start ${r} ${W}) (mid ${r - r * Math.SQRT1_2} ${W - r + r * Math.SQRT1_2}) (end 0 ${W - r}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_line (start 0 ${W - r}) (end 0 ${r}) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);
out.push(`  (gr_arc (start 0 ${r}) (mid ${r - r * Math.SQRT1_2} ${r - r * Math.SQRT1_2}) (end ${r} 0) (layer "Edge.Cuts") (stroke (width 0.1) (type solid)))`);

// edge-launch footprints: centre pad at the edge, four ground pins (rotated 90 / 270)
const launch = (x: number, rot: number) => {
  out.push(`  (footprint "Connector_Coaxial:EdgeLaunch" (layer "F.Cu") (at ${x} ${yc} ${rot})
    (pad "1" smd roundrect (at 0 0 ${rot}) (size 0.8 2) (layers "F.Cu" "F.Mask") (roundrect_rratio 0.2) (net 1 "RF"))
    (pad "2" thru_hole circle (at -2.5 1.5 ${rot}) (size 1.6 1.6) (drill 1.0) (layers "*.Cu" "*.Mask") (net 2 "GND"))
    (pad "3" thru_hole circle (at 2.5 1.5 ${rot}) (size 1.6 1.6) (drill 1.0) (layers "*.Cu" "*.Mask") (net 2 "GND"))
    (pad "4" smd rect (at -2 0.2 ${rot}) (size 1.5 2.4) (layers "F.Cu") (net 2 "GND"))
    (pad "5" smd rect (at 2 0.2 ${rot}) (size 1.5 2.4) (layers "F.Cu") (net 2 "GND"))
  )`);
};
launch(1.0, 90);
launch(L - 1.0, 270);

// the line: straight from the launches, a jog of 45-degree segments with arcs in the middle
const xa = 12, xb = 18;
const seg = (x0: number, y0: number, x1: number, y1: number) =>
  out.push(`  (segment (start ${x0} ${y0}) (end ${x1} ${y1}) (width ${w}) (layer "F.Cu") (net 1))`);
seg(1.0, yc, xa, yc);
seg(xa, yc, xa + jog, yc - jog);
seg(xa + jog, yc - jog, xb - jog, yc - jog);
// a gentle arc (sagitta 0.05 mm) for the down-going diagonal, inside the clearance band
out.push(`  (arc (start ${xb - jog} ${yc - jog}) (mid ${(xb - jog / 2 + 0.05 * Math.SQRT1_2).toFixed(4)} ${(yc - jog / 2 + 0.05 * Math.SQRT1_2).toFixed(4)}) (end ${xb} ${yc}) (width ${w}) (layer "F.Cu") (net 1))`);
seg(xb, yc, L - 1.0, yc);

// top ground pour outside a clearance band around the line (two polygons; above the jog the
// pour steps back to clear it, leaving a wider-gapped section there), and the solid bottom plane
const c = w / 2 + s;
const zone = (layer: string, pts: number[][]) =>
  `  (zone (net 2) (net_name "GND") (layer "${layer}") (uuid "00000000-0000-0000-0000-00000000000${layer === 'F.Cu' ? 1 : 2}") (hatch edge 0.5)
    (connect_pads (clearance 0.2)) (min_thickness 0.25) (fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5))
    (filled_polygon (layer "${layer}") (pts ${pts.map((p) => `(xy ${p[0].toFixed(4)} ${p[1].toFixed(4)})`).join(' ')}))
  )`;
const lower = [[0.2, yc + c], [L - 0.2, yc + c], [L - 0.2, W - 0.2], [0.2, W - 0.2]];
const upper = [[0.2, 0.2], [L - 0.2, 0.2], [L - 0.2, yc - c], [xb + c, yc - c], [xb + c, yc - jog - c], [xa - c, yc - jog - c], [xa - c, yc - c], [0.2, yc - c]];
out.push(zone('F.Cu', upper));
out.push(zone('F.Cu', lower));
out.push(zone('B.Cu', [[0.2, 0.2], [L - 0.2, 0.2], [L - 0.2, W - 0.2], [0.2, W - 0.2]]));

// via fences along the straight parts, stitching vias on a 2 mm grid elsewhere
const via = (x: number, y: number) => out.push(`  (via (at ${x.toFixed(3)} ${y.toFixed(3)}) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net 2))`);
const fence = c + 0.45;
for (let x = 4; x <= xa - 0.5; x += 0.8) { via(x, yc - fence); via(x, yc + fence); }
for (let x = xb + 0.5; x <= L - 4; x += 0.8) { via(x, yc - fence); via(x, yc + fence); }
for (let x = 2; x < L - 1; x += 2) for (let y = 2; y < W - 1; y += 2) if (Math.abs(y - yc) > 3) via(x, y);
out.push(')');
console.log(out.join('\n'));

export {};
