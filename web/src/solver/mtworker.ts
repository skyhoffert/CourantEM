/// <reference lib="webworker" />
/** Helper thread for ThreadPool (see pool.ts). */
declare const self: DedicatedWorkerGlobalScope;
export {};

self.onmessage = async (ev: MessageEvent) => {
  const { module, memory, ctrl, id, stackTop } = ev.data as { module: WebAssembly.Module; memory: WebAssembly.Memory;
    ctrl: Int32Array; id: number; stackTop: number };
  try {
    const inst = await WebAssembly.instantiate(module, { env: { memory }, wasi_snapshot_preview1: { sched_yield: () => 0 } });
    const ex = inst.exports as unknown as { __stack_pointer: WebAssembly.Global; cem_mt_H(s: number, a: number, b: number): void;
      cem_mt_E(s: number, a: number, b: number): void };
    ex.__stack_pointer.value = stackTop; // own shadow stack; never _initialize (the main instance did)
    self.postMessage('ready');
    let gen = Atomics.load(ctrl, 0);
    for (;;) {
      Atomics.wait(ctrl, 0, gen);
      gen = Atomics.load(ctrl, 0);
      const cmd = ctrl[1];
      if (cmd === 3) return;
      const s = ctrl[2], i0 = ctrl[4 + 2 * id], i1 = ctrl[5 + 2 * id];
      if (i1 > i0) { if (cmd === 1) ex.cem_mt_H(s, i0, i1); else ex.cem_mt_E(s, i0, i1); }
      Atomics.add(ctrl, 3, 1);
      Atomics.notify(ctrl, 3);
    }
  } catch (e) {
    self.postMessage(e instanceof Error ? e.message : String(e));
  }
};
