# 0004 — WebAssembly toolchain

Status: clang + wasi-libc in use; Emscripten deferred. **This changes a plan item; Sky to choose.**

## What happened

The plan says C → WebAssembly with Emscripten, pthreads and SIMD. Emscripten could not be installed in the overnight build sandbox (its SDK downloads come from a blocked host). The solver core only needs `malloc`, `libm` and `<complex.h>`, so it was built instead with:

```
clang --target=wasm32-wasip1 --sysroot=<wasi-sysroot> -O3 -msimd128 -mexec-model=reactor \
      -Wl,--export-dynamic -Wl,--export=malloc -Wl,--export=free -Wl,--max-memory=4294967296
```

The core avoids stdio entirely, so the resulting module has **zero imports**. It instantiates with `WebAssembly.instantiate(bytes, {})`: no JS glue, no WASI shim. Size about 200 KB.

Public API functions carry `CEM_API` (`used, visibility("default")` on wasm) so `--export-dynamic` exports exactly the C API.

## Parity

`make -C solver parity` builds `apps/parity.c` both natively and as a WASI command module (run with Node's `node:wasi`) and compares every printed number: ≤ 1e-14 relative (differences are libm rounding). `web/test/parity.test.ts` rebuilds the same scenario through the TypeScript wrapper. `web/test/browser.test.mjs` runs the GCPW demo in headless Chromium and compares S-parameters with the native app.

## Threads (open)

Not used yet. Options:

1. **wasi-sdk `wasm32-wasip1-threads` + a small Worker pool**: keeps the zero-glue approach; needs shared memory, so COOP/COEP via coi-serviceworker as planned.
2. **Emscripten pthreads**: the plan's original choice; more glue but well-trodden.

The OpenMP pragmas in `sim.c` already mark the parallel loops (i-slabs), so either option partitions the same way.

## Other notes

- `long` is 32-bit on wasm32: step counts are capped at 2³¹ (fine; a 40 GHz run is ~10⁵ steps).
- `M_PI` is not defined under `-std=c11`; every source defines a fallback.
- No npm registry access in the sandbox either: the web build uses esbuild directly (`web/build.mjs`) instead of Vite. Swapping to Vite later is a small change; esbuild is also what Vite uses underneath.


## Multithreaded build

`make wasm` also builds `courant-mt.wasm` for `wasm32-wasip1-threads` (atomics, bulk memory, imported shared memory, `__stack_pointer` exported). There are no pthreads in the C code: `cem_step` is split into phases (`cem_mt_pre`, `cem_mt_H(i0, i1)`, `cem_mt_mid`, `cem_mt_E(i0, i1)`, `cem_mt_post`), where the H and E phases — bulk Yee updates plus the CPML corrections, the bulk of the work — run over independent i-slabs. `web/src/solver/pool.ts` starts helper workers that instantiate the same module on the same shared memory, each with its own 1 MB shadow stack set through the exported `__stack_pointer` (only the main instance runs `_initialize`; the threads build's memory initialisation is guarded so it runs once). Phases are handed out through a small `SharedArrayBuffer` control block with `Atomics.wait/notify`, two barriers per time step.

Every cell sees the same operations in the same order as the single-threaded `cem_step`, so results are bit-identical: `web/test/threads.test.mjs` runs the GCPW with 1 and 2 threads and requires zero difference (2 threads: ×1.53 in the 2-core sandbox, browser included). The single-threaded `courant.wasm` stays the default whenever the page is not cross-origin isolated.

Cross-origin isolation on static hosting comes from the vendored `coi-serviceworker` (MIT, v0.1.7, `web/public/`): it re-serves the site with COOP/COEP headers after one reload on first visit, which the same test checks with `COI=sw` (no headers from the server).
