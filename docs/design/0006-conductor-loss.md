# 0006 — Copper loss: surface-impedance boundary and roughness

Status: implemented (`solver/src/sibc.c`); validated on cavity Q; GCPW agreement with the perturbation estimate is 4–14% and open (issue 9).

## Why SIBC

Skin depth in copper is 0.66 µm at 10 GHz, far below any cell. Instead of meshing the conductor, every metal edge that faces a field region carries

```
E_t = Z_s(ω) · J_s,   J_s = n × H   (H taken half a cell outside the surface)
Z_s = (1 + j) √(ωμ0 / 2σ) · K(ω)
```

## Roughness

- Groisse: `K = 1 + exp(−(δ / 2R_q)^1.6)`
- Huray, single sphere size: `K = 1 + (N·4πa²) / (1 + δ/a + δ²/2a²)` with N spheres per unit flat area.

K multiplies the whole Z_s; the fitted model below makes the result causal and passive regardless.

## Time domain

Z_s is fitted over the band by a Foster RL form, passive by construction:

```
Z(s) ≈ Σ_k c_k · s/(s + a_k),   c_k ≥ 0,   a_k log-spaced over [f_min/10, 10 f_max]
```

(non-negative least squares, 12 terms; fit error 0.02% smooth, 2.2% with Groisse at R_q = 1 µm). Each term is a first-order recursion: `y_k ← e^{−a_k Δt} y_k + (1 − e^{−a_k Δt}) J`, `E = Σ c_k (J − y_k)`, applied right after the E update.

## Which edges

An edge is metal if its update coefficient is zero (PEC cells, sheets) or it lies on a PEC outer face. It is a surface edge if a parallel neighbour edge is not metal; each such neighbour direction is a face, and J is the average over faces. Edges inside CPML layers and lumped-port edges stay lossless.

Zero-thickness sheets get two faces; a single E forces equal currents on top and bottom, which under-counts loss when the currents differ (microstrip). Use thick copper (`--tcu`) when loss matters.

## Port model

Mode ports add the conductor attenuation by perturbation of the lossless 2D mode: `α_c = (R_s/2) ∮|n × H|² dl / 2P`, over metal surfaces that face a field region and are outside the PML. On-metal H samples (normal components in this formulation) are excluded.

## Validation

| Test | Result |
|---|---|
| Copper cavity TE101 Q, smooth (analytic, Pozar 6.46) | +0.36% |
| Same with Groisse R_q = 1 µm | +0.98% |
| WR-28 TE10 α_c, 2D perturbation vs closed form | −0.05% |
| GCPW, 35 µm copper, FDTD α vs 2D perturbation, 2–40 GHz | 20 µm cells: −1% to +19%; 10 µm: −4% to +14% |

## Issue 9: what the reference says

`solver/apps/condloss2d.c` solves the same GCPW cross-section in 2D with every metal region a finite-conductivity *volume*, meshed down to 0.6 µm so the skin depth (5 µm at σ = 1e6 S/m, 10 GHz) and the current crowding at the copper corners are resolved. Converged attenuation: 4.70, 4.51, 4.46 Np/m at 2.5, 1.2, 0.6 µm surface cells → about 4.43 Np/m.

| Method (σ = 1e6 S/m, 10 GHz, 35 µm copper) | α (Np/m) | vs volumetric |
|---|---|---|
| 2D volumetric, skin depth resolved (reference) | ~4.43 | — |
| 3D FDTD + SIBC, 20 µm cells | 3.30 | −26% |
| 3D FDTD + SIBC, 10 µm cells | 3.50 | −21% |
| 2D perturbation, same 10 µm mesh | 3.16 | −29% |
| 2D perturbation, 0.6 µm mesh | 4.05 | −9% |

So the SIBC itself is right (cavity Q within 0.4%), and FDTD was the *closer* of the two estimates: both under-count loss because 10–20 µm cells can't see the current crowding within a few skin depths of the copper corners. Expect conductor loss on edge-dominated lines to come out ~20–25% low at practical cell sizes.

Options, in order of preference:

1. **Per-cross-section calibration**: for each port line, compute the ratio (2D volumetric, fine) / (2D perturbation, FDTD mesh) and scale R_s on that line's copper by it. Cheap (seconds in 2D) and exact for uniform lines; approximate for discontinuities.
2. **Corner refinement** to ~δ near copper corners: accurate but cuts Δt by 3–5× in 3D.
3. **Edge-singularity correction** of the surface-impedance coefficients near corners, analogous to the thirds rule.

Not done yet; recorded as issue 9.

### Follow-up: can the reference reach real copper?

Calibration needs a trustworthy fine reference at σ = 5.8e7 S/m (δ = 0.66 µm at 10 GHz). Sweeping `condloss2d` (α in Np/m at 10 GHz; "pert" is the PEC perturbation on the same fine mesh, scaled by that σ's R_s):

| σ (S/m), δ | min cell 1.2 µm | 0.6 µm | 0.3 µm | 0.15 µm | 0.075 µm |
|---|---|---|---|---|---|
| 1e6, 5.0 µm — volumetric | 4.506 | 4.462 | 4.450 | | |
| 1e6 — pert | 3.900 | 4.045 | 4.174 | | |
| 1e7, 1.6 µm — volumetric | 1.570 | 1.463 | 1.426 | | |
| 5.8e7, 0.66 µm — volumetric | 0.593 | 0.669 | 0.621 | 0.604 | 0.643 |
| 5.8e7 — pert | 0.512 | 0.531 | 0.548 | 0.692 | (30.4, breaks down) |

- The perturbation estimate climbs slowly with refinement (the corner field singularity, |H|² ~ r^(−2/3), converges like h^(1/3)), so even 0.3 µm cells leave it ~10% short.
- At real copper the volumetric solve does not settle (±5% scatter between 0.6 and 0.075 µm) and the perturbation breaks down below 0.15 µm: this mesher's grading (×1.25 up to 25 µm) and the direct eigen-solve are not well conditioned there. The best copper estimate is α ≈ 0.60–0.67 Np/m ≈ (23–26)·R_s, against FDTD's ≈ 17.6·R_s at 10 µm cells, i.e. FDTD **25–33% low** on this edge-dominated line.

Conclusion: per-line calibration (option 1) is the right shape, but its reference has to be a dedicated quasi-static conductor-loss solve rather than this eigen-solver. Candidates: Wheeler's incremental-inductance rule on the Laplace solver of `tools/gcpw_ref.py` with corner rounding at r ≈ δ, or a surface-integral (BEM) current-distribution solve. Needs Sky's call on how much accuracy milestone 1 needs before the measured comparisons.
