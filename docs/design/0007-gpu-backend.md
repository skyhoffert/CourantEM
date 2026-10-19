# 0007 — WebGPU backend (f32)

Status: implemented (`web/src/gpu/`). Every feature of the C core runs on the GPU: Yee updates, CPML, lumped and modal sources, probes, Debye dispersion and surface-impedance copper.

## Split of responsibilities

The C core (WebAssembly, f64) builds and finalizes the simulation: mesh, materials, edge coefficients, CPML profiles, ports, mode solves, probe definitions. The GPU engine:

1. reads the finalized plan through `cem_plan_*` (coefficient arrays, CPML slabs, source and probe lists, ADE and SIBC tables),
2. uploads it as f32 buffers,
3. time-steps on the GPU, recording probe sums each step into a small series buffer,
4. every chunk (256–500 steps) reads the series back and hands it to `cem_push_samples`, which performs exactly the DFTs `cem_step` would.

So geometry, port de-embedding, renormalisation and Touchstone output are a single code path, and "same mesh, same S-parameters" is by construction.

## Kernels

`fdtd.wgsl.ts`: `updateH`, `updateE` (one thread per grid point, three components each, same update ranges as C), `pmlH{x,y,z}`, `pmlE{x,y,z}` (one thread per point of the low and high slabs), `inject`, `probes`, `tick` (step counter on the GPU so a whole chunk is encoded in one pass), `dispSave/Correct/Update`, `sibc`. Dispatch order per step mirrors `cem_step`.

Storage: E, H, ca, cb as component-major arrays of 3N f32; 1D spacing and CPML profiles packed into one buffer with offsets in the uniform block; 12N CPML auxiliaries. This keeps every entry point within the default 8 storage bindings. One buffer per quantity caps problems at ~10 M points under the default 128 MB binding limit; the multi-buffer split from the plan is still to do.

## Verification (headless Chromium, SwiftShader)

| Case | GPU f32 vs CPU f64 |
|---|---|
| Parity scene (non-uniform grid, CPML, lossy dielectric, PEC box and sheet, lumped ports, probes), fields after 600 steps | 1.3e-6 of max |
| GCPW with mode ports, S21 | 5e-6 dB, 4e-5° |
| Same with Debye substrate and copper SIBC | 4e-6 dB, 3e-5° |

Plan tolerance is ±0.01 dB and ±0.1°. Speed was not measured: SwiftShader is a CPU emulation (about 7× slower than the wasm path here). Run `node web/test/gpu.test.mjs` on a machine with a GPU for real numbers.

## Run control

The CPU path stops on total field energy. Reading the fields back every chunk would cost more than the step, so the GPU path stops when the squared peak probe amplitude of the last chunk falls below `decay` times its maximum.


## Buffer limits

The engine requests the adapter's maximum storage-binding and buffer sizes. With the compact CPML layout (note 0001) the largest binding is a 3N-float field array, so the point limit is maxStorageBufferBindingSize / 12 bytes (~11 M points at the 128 MB default, far more on desktop GPUs, which usually allow 1–4 GB). A mesh over the limit is refused up front with a message suggesting the CPU backend; splitting E/H per component would buy another 3× but needs more storage bindings per stage than the WebGPU default of 8.
