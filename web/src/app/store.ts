/**
 * Central typed state store with one-way data flow. Panels read state in
 * update(state) and change it only through dispatch(action).
 */
import type { Board } from '../import/kicad.ts';
import type { StackupRow } from '../import/gerber.ts';
import type { BoardInput } from '../import/load.ts';
import type { Complex } from '../solver/courant.ts';

export interface RunParams {
  hmin: number; dx: number; len: number; er: number;
  fmax: number; df: number; decay: number;
  ports: 'mode' | 'lumped'; lossy: boolean; thick: boolean; full2: boolean; ywin: number; tmax: number;
}

export interface Result {
  f: number[]; s11: Complex[]; s21: Complex[];
  s11_50?: Complex[]; s21_50?: Complex[]; z0?: number[];
  /** present when port 2 was excited too (full 2-port) */
  s12?: Complex[]; s22?: Complex[]; s12_50?: Complex[]; s22_50?: Complex[];
  steps: number; wall: number;
  /** set when the run hit its time limit before the field energy decayed */
  capped?: { energyDb: number; tMax: number };
  snapshot?: { nx: number; ny: number; x: number[]; y: number[]; z: number; data: Float32Array };
}

export interface State {
  source: 'gcpw' | 'board';
  params: RunParams;
  backend: 'cpu' | 'gpu';
  threads: number;           // 0 = auto
  threadsUsed?: number;
  isolated?: boolean;
  gpuAvailable: boolean;
  ref50: boolean;
  board?: { name: string; input: BoardInput; stackup: StackupRow[]; offsets: [number, number]; model: Board; nets: string[]; net: string; rev: number };
  run: { phase: 'idle' | 'running' | 'done' | 'error'; message: string; progress: number; cells?: number[] };
  result?: Result;
}

export type Action =
  | { type: 'params'; patch: Partial<RunParams> }
  | { type: 'source'; source: State['source'] }
  | { type: 'backend'; backend: State['backend'] }
  | { type: 'ref50'; on: boolean }
  | { type: 'threads'; threads: number }
  | { type: 'threads-used'; threads: number; isolated: boolean }
  | { type: 'board'; board: Omit<NonNullable<State['board']>, 'rev'> }
  | { type: 'net'; net: string }
  | { type: 'run-start' }
  | { type: 'run-setup'; cells: number[]; message: string }
  | { type: 'run-progress'; progress: number; message: string }
  | { type: 'run-done'; result: Result; message: string }
  | { type: 'run-error'; message: string }
  | { type: 'status'; message: string };

export function reduce(s: State, a: Action): State {
  switch (a.type) {
    case 'params': return { ...s, params: { ...s.params, ...a.patch } };
    case 'source': return { ...s, source: a.source };
    case 'backend': return { ...s, backend: a.backend };
    case 'ref50': return { ...s, ref50: a.on };
    case 'threads': return { ...s, threads: a.threads };
    case 'threads-used': return { ...s, threadsUsed: a.threads, isolated: a.isolated };
    case 'board': return { ...s, board: { ...a.board, rev: (s.board?.rev ?? 0) + 1 }, source: 'board', result: undefined };
    case 'net': return s.board ? { ...s, board: { ...s.board, net: a.net } } : s;
    case 'run-start': return { ...s, run: { phase: 'running', message: 'Building mesh…', progress: 0 }, result: undefined };
    case 'run-setup': return { ...s, run: { ...s.run, cells: a.cells, message: a.message } };
    case 'run-progress': return { ...s, run: { ...s.run, progress: a.progress, message: a.message } };
    case 'run-done': return { ...s, run: { ...s.run, phase: 'done', progress: 1, message: a.message }, result: a.result };
    case 'run-error': return { ...s, run: { ...s.run, phase: 'error', message: a.message } };
    case 'status': return { ...s, run: { ...s.run, message: a.message } };
  }
}

export interface Panel { update(s: State, prev: State | undefined): void }

export class Store {
  private state: State;
  private panels: Panel[] = [];
  constructor(initial: State) { this.state = initial; }
  get(): State { return this.state; }
  dispatch(a: Action): void {
    const prev = this.state;
    this.state = reduce(prev, a);
    for (const p of this.panels) p.update(this.state, prev);
  }
  mount(p: Panel): void { this.panels.push(p); p.update(this.state, undefined); }
}
