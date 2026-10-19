/**
 * One entry point for board inputs, shared by the app and the solver worker:
 * a KiCad file or a Gerber + drill set, plus the stackup the user edited.
 */
import { importKicad, type Board, type StackLayer } from './kicad.ts';
import { defaultStackup, importGerber, type GerberFile, type StackupRow } from './gerber.ts';
import { trackAcross } from './boardsim.ts';

export type BoardInput = { kind: 'kicad'; name: string; text: string } | { kind: 'gerber'; name: string; files: GerberFile[] };

export interface LoadedBoard { board: Board; stackup: StackupRow[]; warnings: string[]; portNets: string[] }

/** Stackup rows (top to bottom) as the editor shows them. */
export function stackupOf(b: Board): StackupRow[] {
  return [...b.layers].reverse().filter((l) => l.type !== 'other').map((l) => ({ name: l.name, type: l.type as StackupRow['type'],
    thickness: l.thickness, epsR: l.epsR, lossTangent: l.lossTangent, roughness: l.roughness, mask: l.mask }));
}

/** Replace a board's layer stack by edited rows (same copper layers, top to bottom). */
export function applyStackup(b: Board, rows: StackupRow[]): Board {
  const cu = rows.filter((r) => r.type === 'copper').map((r) => r.name);
  const have = [...b.copperLayers].reverse();
  if (cu.join() !== have.join()) throw new Error(`stackup copper layers ${cu.join(', ')} do not match the board's ${have.join(', ')}`);
  const layers: StackLayer[] = [...rows].reverse().map((r) => ({ name: r.name, type: r.type, thickness: r.thickness, epsR: r.epsR,
    lossTangent: r.lossTangent, roughness: r.roughness, mask: r.mask, z0: 0, z1: 0 }));
  let z = 0;
  for (const l of layers) { if (l.type === 'dielectric') { l.z0 = z; z += l.thickness; l.z1 = z; } else { l.z0 = z; l.z1 = z; } }
  return { ...b, layers, thickness: z };
}

/** Nets with an x-directed track crossing both port planes (offsets measured into the board from its edges). */
export function portNets(b: Board, offLo = 0, offHi = 0): string[] {
  const x0 = b.outline.min[0] + offLo, x1 = b.outline.max[0] - offHi, out: string[] = [];
  for (const [n, name] of b.nets) {
    if (n && trackAcross(b, n, x0) && trackAcross(b, n, x1)) out.push(name || String(n));
  }
  return out;
}

export function loadBoard(input: BoardInput, stackup?: StackupRow[], offsets: [number, number] = [0, 0]): LoadedBoard {
  if (input.kind === 'kicad') {
    let board = importKicad(input.text);
    if (stackup) board = applyStackup(board, stackup);
    return { board, stackup: stackupOf(board), warnings: [], portNets: portNets(board, ...offsets) };
  }
  const g = importGerber(input.files, stackup);
  return { board: g.board, stackup: stackup ?? stackupOf(g.board), warnings: g.warnings, portNets: portNets(g.board, ...offsets) };
}

/** Classify a set of picked files: one .kicad_pcb, or Gerbers + drills. */
export async function inputFromFiles(files: File[]): Promise<BoardInput> {
  const kicad = files.find((f) => /\.kicad_pcb$/i.test(f.name));
  if (kicad) return { kind: 'kicad', name: kicad.name, text: await kicad.text() };
  const gf = await Promise.all(files.map(async (f) => ({ name: f.name, text: await f.text() })));
  // common prefix of the file names, e.g. "board-F_Cu.gbr", "board-B_Cu.gbr" -> "board"
  let stem = files[0]?.name ?? 'board';
  for (const f of files) while (!f.name.startsWith(stem)) stem = stem.slice(0, -1);
  stem = stem.replace(/[-_.\s]+$/, '') || 'board';
  return { kind: 'gerber', name: `${stem} (Gerber)`, files: gf };
}

export { defaultStackup };
