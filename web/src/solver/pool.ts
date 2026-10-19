/**
 * Thread pool for the multithreaded CPU build (courant-mt.wasm): helper workers
 * instantiate the same module on the same shared memory, each with its own
 * shadow stack, and run the bulk H and E updates (cem_mt_H / cem_mt_E) on
 * i-slabs while the calling worker runs its own slab and the serial phases.
 * The operation order per cell is unchanged, so results are bit-identical to
 * the single-threaded build. Requires cross-origin isolation (SharedArrayBuffer)
 * and must run in a worker (Atomics.wait).
 *
 * Control block (Int32Array on a SharedArrayBuffer):
 *   [0] generation  [1] command (1 H, 2 E, 3 exit)  [2] sim handle  [3] done count
 *   [4 + 2k], [5 + 2k]  slab [i0, i1) for helper k
 */
export interface MtExports {
  cem_mt_pre(s: number): void; cem_mt_H(s: number, i0: number, i1: number): void; cem_mt_mid(s: number): void;
  cem_mt_E(s: number, i0: number, i1: number): void; cem_mt_post(s: number): void;
  malloc(n: number): number;
}

export const STACK_BYTES = 1 << 20;

export class ThreadPool {
  private readonly ctrl: Int32Array;
  private readonly workers: Worker[] = [];
  private readonly ex: MtExports;
  readonly threads: number;
  private constructor(ex: MtExports, threads: number) {
    this.ex = ex; this.threads = threads;
    this.ctrl = new Int32Array(new SharedArrayBuffer(4 * (4 + 2 * threads)));
  }

  /** Start threads-1 helpers. `workerUrl` is the bundled mtworker.js. */
  static async create(module: WebAssembly.Module, memory: WebAssembly.Memory, ex: MtExports, threads: number,
                      workerUrl: string | URL): Promise<ThreadPool> {
    const pool = new ThreadPool(ex, threads);
    const ready: Promise<void>[] = [];
    for (let k = 0; k < threads - 1; k++) {
      const stack = ex.malloc(STACK_BYTES);
      if (!stack) throw new Error('out of memory for thread stacks');
      const w = new Worker(workerUrl, { type: 'module' });
      ready.push(new Promise((res, rej) => {
        w.onmessage = (e) => (e.data === 'ready' ? res() : rej(new Error(String(e.data))));
        w.onerror = (e) => rej(new Error(e.message));
      }));
      // 16-byte aligned top of this helper's stack
      w.postMessage({ module, memory, ctrl: pool.ctrl, id: k, stackTop: (stack + STACK_BYTES) & ~15 });
      pool.workers.push(w);
    }
    await Promise.all(ready);
    return pool;
  }

  private parallel(cmd: 1 | 2, sim: number, lo: number, hi: number): void {
    const T = this.threads, c = this.ctrl, n = hi - lo;
    const cut = (k: number) => lo + Math.floor((n * k) / T);
    for (let k = 1; k < T; k++) { c[4 + 2 * (k - 1)] = cut(k); c[5 + 2 * (k - 1)] = cut(k + 1); }
    c[1] = cmd; c[2] = sim;
    Atomics.store(c, 3, 0);
    Atomics.add(c, 0, 1);
    Atomics.notify(c, 0);
    if (cmd === 1) this.ex.cem_mt_H(sim, lo, cut(1)); else this.ex.cem_mt_E(sim, lo, cut(1));
    for (;;) {
      const d = Atomics.load(c, 3);
      if (d >= T - 1) break;
      Atomics.wait(c, 3, d);
    }
  }

  /** n time steps of the simulation `sim` with nx cells along x. */
  step(sim: number, nx: number, n: number): void {
    const ex = this.ex;
    for (let it = 0; it < n; it++) {
      ex.cem_mt_pre(sim);
      this.parallel(1, sim, 0, nx + 1);
      ex.cem_mt_mid(sim);
      this.parallel(2, sim, 0, nx);
      ex.cem_mt_post(sim);
    }
  }

  terminate(): void {
    this.ctrl[1] = 3;
    Atomics.add(this.ctrl, 0, 1);
    Atomics.notify(this.ctrl, 0);
    for (const w of this.workers) w.terminate();
  }
}
