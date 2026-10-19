# 0003 — Mode ports

Status: implemented for x-directed lines (`solver/src/mode2d.c`, `solver/src/modeport.c`). Validated on a uniform GCPW.

## Goal

Inject and measure the line's own guided mode so S-parameters have no port parasitics and reference planes sit where the TRL-corrected measurements will be.

## 2D eigenmode solver

Fields `∝ exp(−jβx)`, transverse plane (y, z), placed exactly where the 3D Yee grid puts them:

```
E_u = Ey(j+½, k)   E_v = Ez(j, k+½)   E_w = Ex(j, k)
H_u = Hy(j, k+½)   H_v = Hz(j+½, k)   H_w = Hx(j+½, k+½)
```

Eliminating the longitudinal components:

```
β E_u =  ωμ H_v + (1/ω) ∂_u [ (∂_u H_v − ∂_v H_u) / ε_w ]
β E_v = −ωμ H_u + (1/ω) ∂_v [ (∂_u H_v − ∂_v H_u) / ε_w ]
β H_u = −ω ε_v E_v − (1/ωμ) ∂_u (∂_u E_v − ∂_v E_u)
β H_v =  ω ε_u E_u − (1/ωμ) ∂_v (∂_u E_v − ∂_v E_u)
```

i.e. `β e = P h`, `β h = Q e`, so `β² e = P Q e`. Because the cross-section uses the FDTD's own mesh and material averaging, the discrete mode is the one the 3D grid supports.

- Complex arithmetic: lossy media (`ε − jσ/ω`) and stretched-coordinate PML (`s = κ + σ/(α + jωε0)`, same profiles as the 3D CPML) on every derivative.
- `A = P Q` is applied matrix-free and assembled into band storage by probing with coloured unit vectors (columns spaced `2·bw+1` apart). Unknown ordering `(j outer, k inner, [E_u, E_v])` gives bandwidth ≈ `2(nz+2)`.
- Banded LU with partial pivoting, then inverse iteration with shift `σ = (k0·n_guess)²`. Convergence requires both λ and the phase-aligned eigenvector to settle.
- Frequency sweeps walk outward from `f0`, seeding each solve with the previous mode so the same mode is tracked.

Verified: on a rectangular waveguide the Python prototype reproduces the discrete dispersion exactly; the C solver matches the Python prototype to all printed digits; against a 2D Laplace solution in the same PEC box the quasi-static limit agrees within 0.1%.

**Why the PML matters**: a PEC lid 1.5 mm above the GCPW lowers ε_eff by 1.2%. With the cross-section closed by PEC while the 3D domain is open, the ports reflected at −49 dB and S21 phase drifted 1.2°. With the CPML carried into the 2D solve: −62 dB and 0.05°.

## Excitation

A soft E-field source with the f0 mode profile on a node plane two cells in front of the x CPML. It launches the mode both ways; the backward half goes into the CPML, along which the line continues. Time function: differentiated Gaussian (zero DC, so no static charge is left behind).

## Wave separation

At a monitor plane between source and device:

- `V` at node planes `m` and `m+1`, `I` (loop around the conductor) at the half plane between.
- Discrete consistency: the 2D solve runs at the pre-warped frequency `(2/Δt) sin(ωΔt/2)/2π`; its eigenvalue κ maps to the 3D constant by `κ = (2/Δx) sin(β Δx/2)`; V is moved to the half plane with `V = (V_m + V_{m+1}) / (2 cosh(γΔx/2))`.
- `Z0(f) = V/I` of the 2D mode (same path and loop), so `a = (V + Z0 I)/2√Z0`, `b = (V − Z0 I)/2√Z0`.
- Shift to the reference plane: `a·e^{−γd}`, `b·e^{+γd}`.

S-parameters come out normalised to each port's `Z0(f)`; renormalisation to 50 Ω is a post-processing step (not yet implemented).

The monitor sits 2.5 mm from the source. At 1 mm, the evanescent channel mode near its ~52 GHz cutoff still reached the monitor and S11 was only −42 dB at 40 GHz.

## Results (uniform GCPW, 10 mm, hmin 20 µm, open)

|S11| ≤ −62 dB, ||S21| − 1| ≤ 0.004 dB, ∠S21 error vs `−βL` ≤ 0.05° from 1 to 40 GHz. See `validation/gcpw-m1.md`.

## Known limits

- Platform-to-platform agreement of mode-port results is ~1e-6 (vs 1e-14 for the bare kernel): at low frequency `k0² ≪ 1/Δ²` and the eigenproblem is poorly conditioned. Physically irrelevant; it sets the wasm-vs-native test tolerance.
- V/I separation assumes a single quasi-TEM mode and needs the voltage path and current loop from the port definition. Overlap-integral separation would remove that need and reject other modes; planned when ports come from KiCad nets.
- Only x-directed ports so far.

## Rejected alternatives

- **TF/SF or Huygens (J + M) unidirectional source**: not needed when the line continues into the CPML; the bidirectional soft source is simpler and the wasted half costs nothing.
- **Mode profile per output frequency in the source**: impossible in one time-domain run; a single f0 profile plus per-frequency decomposition is the standard approach.
- **Z0 and β from multi-plane probe fits on the feed line**: works (it is how `analyze_gcpw.py` extracts the line), but it needs long feed lines and is noisier at low frequency.
