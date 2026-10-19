# 0009 — Checkpoints and resume

Status: implemented (`solver/src/export.c`, `web/src/solver/checkpoints.ts`, worker).

## What a checkpoint holds

`cem_checkpoint_save` writes one flat `double` buffer: an 8-value header (magic `20261017`, format version 1, cell count N, frequency count, step, body size, a port/probe signature, dispersive + surface-impedance edge count) and the body — E and H, CPML ψ, Debye currents and previous E, SIBC filter states, lumped-port V/I accumulators, probe DFT accumulators, the source spectrum, and the DFT phasors. Geometry, materials and coefficients are *not* saved: a checkpoint is only valid for the simulation that produced it, rebuilt from the same inputs, and `cem_checkpoint_load` refuses one whose header does not match.

Resuming is bit-exact on the CPU: `test_checkpoint` runs 2·n steps straight and n + save + load into a fresh sim + n, and compares every field and spectrum exactly.

## In the app

The worker hashes the run request (SHA-256 of the JSON message, minus the wasm URL) to a key and stores checkpoints in IndexedDB under it, every ~20 s of wall time. Starting the same run again — after a closed tab, a crash or a reload — finds the key, loads the checkpoint (plus the run-control energy peak, appended as one extra value) and continues; the status line says "Resumed from a checkpoint at step n". The checkpoint is deleted when the run finishes. A checkpoint that no longer matches (e.g. after a solver update changed the layout) is ignored and the run starts over.

On the GPU the state lives in GPU buffers. Before saving, `GpuEngine.syncToSim()` copies E, H, ψ, Debye and SIBC states back into the C sim (f32 → f64); the DFT accumulators are already in C because the GPU pushes probe samples every chunk. A GPU engine created on a restored sim uploads that state and continues the source waveforms from the restored step, so a checkpoint can be resumed on either backend.

Bug found on the way: `cem_push_samples` set the DFT phasors directly for each sample but left them at the last sample's time, whereas `cem_step` leaves them advanced to the next step. A checkpoint taken after GPU steps, resumed on the CPU, therefore summed every later sample with a one-step phase lag (0.56 dB / 2.9° error in S21). `cem_push_samples` now ends by setting the phasors for the next step.

## Validation

| Test | Result |
|---|---|
| CPU, lossy GCPW (Debye + SIBC), n + checkpoint + n vs 2n | fields and S21 bit-identical |
| GPU n steps → `syncToSim` → checkpoint → fresh CPU sim → n steps vs CPU 2n, lossy GCPW | S21 within 2e-6 dB, 1e-5° (was 0.56 dB / 2.9° before the phasor fix) |
