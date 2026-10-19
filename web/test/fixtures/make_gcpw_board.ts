/**
 * Writes a KiCad 7 board with a 50-ohm GCPW line on 10-mil RO4350B, top ground
 * pours, a solid bottom plane and via fences. Used by the import tests.
 * Run: node web/test/fixtures/make_gcpw_board.ts > web/test/fixtures/gcpw.kicad_pcb
 */
const L = 10, W = 6;                 // board size, mm
const yc = W / 2, w = 0.5, s = 0.2;   // strip centre, width, gap
const drill = 0.3, pitch = 0.8;
const fenceOff = 0.3 + drill / 2;     // via centre from gap edge
const out: string[] = [];
out.push(`(kicad_pcb (version 20221018) (generator pcbnew)
  (general (thickness 0.324))
  (paper "A4")
  (layers (0 "F.Cu" signal) (31 "B.Cu" signal) (44 "Edge.Cuts" user))
  (setup
    (stackup
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "core") (thickness 0.254) (material "RO4350B") (epsilon_r 3.66) (loss_tangent 0.0037))
      (layer "B.Cu" (type "copper") (thickness 0.035))
    )
  )
  (net 0 "")
  (net 1 "RF")
  (net 2 "GND")`);
for (const [x0, y0, x1, y1] of [[0, 0, L, 0], [L, 0, L, W], [L, W, 0, W], [0, W, 0, 0]]) {
  out.push(`  (gr_line (start ${x0} ${y0}) (end ${x1} ${y1}) (layer "Edge.Cuts") (width 0.1))`);
}
out.push(`  (segment (start 0 ${yc}) (end ${L} ${yc}) (width ${w}) (layer "F.Cu") (net 1))`);
const zone = (layer: string, pts: number[][]) =>
  `  (zone (net 2) (net_name "GND") (layer "${layer}") (hatch edge 0.5)
    (filled_polygon (layer "${layer}") (pts ${pts.map((p) => `(xy ${p[0]} ${p[1]})`).join(' ')}))
  )`;
const g1 = yc - w / 2 - s, g2 = yc + w / 2 + s;
out.push(zone('F.Cu', [[0, 0], [L, 0], [L, g1], [0, g1]]));
out.push(zone('F.Cu', [[0, g2], [L, g2], [L, W], [0, W]]));
out.push(zone('B.Cu', [[0, 0], [L, 0], [L, W], [0, W]]));
for (let x = pitch / 2; x < L; x += pitch) {
  for (const y of [g1 - fenceOff + drill / 2 - drill / 2, g2 + fenceOff - drill / 2 + drill / 2]) {
    out.push(`  (via (at ${x.toFixed(3)} ${y.toFixed(3)}) (size 0.6) (drill ${drill}) (layers "F.Cu" "B.Cu") (net 2))`);
  }
}
out.push(')');
console.log(out.join('\n'));

export {};
