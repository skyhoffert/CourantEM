/**
 * The app's built-in GCPW line as a board model: 1 oz (35 um) copper, a solid
 * bottom plane and a row of stitching vias on each side tying the top grounds to
 * it. It runs through the board path (boardsim.ts), so what the 3D view draws is
 * what is simulated. Same line as test/fixtures/gcpw.kicad_pcb: strip 500 um,
 * gaps 200 um on 254 um RO4350B-like core, 0.3 mm drills at 0.8 mm pitch whose
 * barrels start 300 um from the gap edges (where scenarios/gcpw.ts, the native
 * reference, puts its idealised solid fence walls). The copper stops 0.1 mm short
 * of the board edge: copper edges on the outline put mesh lines a third of a cell
 * apart there, which would halve the time step for no physical gain.
 */
import type { Board, Pt, Via } from '../import/kicad.ts';

export interface GcpwBoardParams {
  len: number;   // m
  er: number;
}

export const GCPW_NET = 'RF';

export function gcpwBoard(p: GcpwBoardParams): Board {
  const h = 254e-6, w = 500e-6, s = 200e-6, tcu = 35e-6;
  const drill = 0.3e-3, pad = 0.6e-3, pitch = 0.8e-3, fence = 300e-6;
  const a = w / 2, b = a + s, yv = b + fence + drill / 2, W = yv + 1.5e-3, edge = W + 0.1e-3;
  const rect = (x0: number, y0: number, x1: number, y1: number): Pt[] => [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
  const vias: Via[] = [];
  for (let i = 0; (i + 0.5) * pitch < p.len - 1e-9; i++) {
    for (const y of [-yv, yv]) vias.push({ at: [(i + 0.5) * pitch, y], drill, size: pad, layers: ['F.Cu', 'B.Cu'], net: 2 });
  }
  return {
    layers: [
      { name: 'B.Cu', type: 'copper', thickness: tcu, z0: 0, z1: 0 },
      { name: 'core', type: 'dielectric', thickness: h, epsR: p.er, lossTangent: 0.0037, z0: 0, z1: h },
      { name: 'F.Cu', type: 'copper', thickness: tcu, z0: h, z1: h },
    ],
    copperLayers: ['B.Cu', 'F.Cu'],
    tracks: [{ a: [0, 0], b: [p.len, 0], width: w, layer: 'F.Cu', net: 1 }],
    vias,
    polys: [
      { pts: rect(0, -W, p.len, W), layer: 'B.Cu', net: 2 },
      { pts: rect(0, -W, p.len, -b), layer: 'F.Cu', net: 2 },
      { pts: rect(0, b, p.len, W), layer: 'F.Cu', net: 2 },
    ],
    outline: { min: [0, -edge], max: [p.len, edge] },
    nets: new Map([[1, GCPW_NET], [2, 'GND']]),
    thickness: h,
  };
}
