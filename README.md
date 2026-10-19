# Courant EM

A browser-based, full-wave 3D electromagnetic solver that turns PCB designs into validated S-parameters. No install, no backend: the solver runs in the visitor's browser.

**Status: milestone 1 in progress (CalVer `20261017.0`, pre-release).** Built so far: the double-precision FDTD core (WebAssembly, multithreaded when the page is cross-origin isolated), the single-precision WebGPU backend with full feature parity, mode ports, dielectric and copper loss models, a rule-based mesher, KiCad and Gerber + drill import with a stackup editor, checkpoint/resume, and an in-browser app (3D model view, field snapshot, S-parameter plots, Touchstone export). Validated on uniform GCPW and microstrip lines against Laplace references (with copper thickness and solder mask). Not done yet: openEMS/HFSS/measured comparisons (`validation/PROTOCOL.md`); conductor loss on edge-dominated lines runs ~25% low (issue 9). See `docs/` for design notes and `validation/` for results.

## What works today

| Piece | Where | State |
|---|---|---|
| FDTD core, f64: non-uniform Yee grid, CPML, PEC cells and zero-thickness sheets, lossy dielectrics, lumped ports, DFT V/I probes | `solver/src/sim.c` | tested |
| 2D full-wave eigenmode solver on the FDTD's own cross-section (complex, PML-stretched) | `solver/src/mode2d.c` | matches Python prototype to 6+ digits |
| Mode ports (modal source + V/I wave separation with discrete-dispersion corrections), 50 Ω renormalisation | `solver/src/modeport.c` | uniform-line S11 ≤ −62 dB, ∠S21 error ≤ 0.05° to 40 GHz |
| Djordjevic–Sarkar dielectrics as fitted Debye poles (ADE) | `solver/src/dispersive.c` | lossy line within 0.0044 dB of the 2D mode |
| Copper loss: surface-impedance boundary, Groisse/Huray roughness | `solver/src/sibc.c` | cavity Q within 1% |
| WebAssembly build (clang + wasi-libc, zero imports) | `solver/Makefile` → `make wasm` | parity with native ≤ 1e-14 |
| WebGPU f32 backend: all kernels incl. Debye and SIBC | `web/src/gpu/` | GCPW S21 within 5e-6 dB / 4e-5° of the f64 core |
| Rule-based mesher + pre-run estimator | `web/src/mesh/mesher.ts` | rules checked in `test/mesher.test.ts` |
| KiCad import (stackup incl. solder mask, tracks and arcs, zones, pads, vias) and Gerber + Excellon import with a stackup editor | `web/src/import/` | microstrip Z0 within 1% of Hammerstad–Jensen; Gerber and KiCad imports of one board simulate identically |
| Board → simulation: copper thickness (default) or sheets, port planes inside the board, mode ports, full 2-port | `web/src/import/boardsim.ts` | GCPW with 35 µm copper within 0.5% of Laplace |
| Checkpoint/resume (IndexedDB), CPU threads (wasm threads + worker pool), coi-serviceworker | `web/src/solver/` | bit-identical resume and threading |
| App: store + panels, three.js model view with field snapshot, uPlot plots, Web Worker, CPU/GPU, Touchstone export | `web/src/app/`, `web/src/main.ts` | browser tests in headless Chromium |
| Step-9 kit: Touchstone comparison with the plan's tolerances, openEMS script, CLI board runner | `tools/compare_s2p.py`, `tools/openems_gcpw.py`, `web/tools/run-board.ts` | `validation/PROTOCOL.md` |
| Reference models: GCPW conformal mapping, 2D Laplace FV with Richardson extrapolation | `tools/gcpw_ref.py` | |
| Python eigenmode prototype | `tools/modesolver2d.py` | |

## Layout

```
solver/            C solver core (f64 reference; OpenMP natively, wasm threads in the browser)
  include/courant.h  public C API (flat, WebAssembly-friendly)
  src/               sim.c (FDTD), mode2d.c (2D eigenmodes), modeport.c (mode ports)
  tests/             cavity resonance, CPML reflection, native/wasm parity
  apps/              gcpw_line (build-step-1 scenario), parity
web/               TypeScript app
  src/solver/        wasm wrapper (courant.ts), Web Worker, thread pool, checkpoints
  src/gpu/           WebGPU engine and WGSL kernels
  src/mesh/          rule-based mesher, estimator
  src/import/        KiCad and Gerber/Excellon readers, stackup, board -> simulation
  tools/             run-board.ts (command-line board runs)
  src/app/           store, panels, 3D view, plots
  test/              parity, mesher, import, browser, GPU tests; fixtures
tools/             Python references, analysis and plotting
docs/design/       design notes: equations, choices, rejected alternatives
validation/        validation notes and figures
```

## Building

Native solver and tests (needs gcc or clang, OpenMP optional):

```sh
make -C solver           # library, tests, apps
make -C solver test      # cavity + CPML tests
```

WebAssembly (needs clang ≥ 17 with the wasm32 target, plus the wasi-sysroot and builtins from a [wasi-sdk release](https://github.com/WebAssembly/wasi-sdk/releases)):

```sh
make -C solver wasm WASI_SYSROOT=/path/to/wasi-sysroot WASI_RT=/path/to/libclang_rt.builtins-wasm32-wasi
make -C solver parity    # native vs wasm on a scenario that touches every kernel feature
```

Web app (Node ≥ 22.18 for type stripping; esbuild, TypeScript, three, uplot from npm):

```sh
cd web && npm install && npm run build   # writes web/dist
python3 -m http.server -d web/dist 8000  # then open http://localhost:8000
npm test                                 # parity, mesher, KiCad, Gerber, thickness, import, port planes (Node + wasm)
npm run test:browser                     # headless Chromium: UI, GPU vs CPU, board uploads, resume, threads
```

`web/test/browser.test.mjs` optionally takes a native run directory (from `gcpw_line`) to compare against; `gpu.test.mjs` needs WebGPU (`--enable-unsafe-webgpu`; SwiftShader works for correctness). Multithreading needs a cross-origin-isolated page: `python3 -m http.server` doesn't send COOP/COEP, so the bundled service worker sets them after one reload.

Boards from the command line (writes a version-stamped .s2p):

```sh
node web/tools/run-board.ts board.kicad_pcb --net RF --planes 1.5,1.5 --lossy --full2 --out board.s2p
python3 tools/compare_s2p.py board.s2p reference.s2p --rule hfss --plot compare.png
```

## Running the milestone-1 line

```sh
./solver/build/gcpw_line --modeports --hmin 20 --dx 50 --df 0.25 --out runs/gcpw
python3 tools/eval_modeports.py runs/gcpw
```

Options: `--hmin` (µm), `--dx` (µm), `--len` (mm), `--er`, `--fmax` (GHz), `--df` (GHz), `--modeports` or lumped (default), `--closed` (PEC box instead of CPML), `--snap` (edge on a grid line instead of the thirds rule), `--modes-only`, `--debye FILE` (from `tools/dielectric.py`), `--cu SIGMA` (copper loss), `--rq UM` (Groisse roughness), `--tcu UM` (copper thickness), `--yair/--zair MM` (air margins).

## License

MIT. See `LICENSE`.
