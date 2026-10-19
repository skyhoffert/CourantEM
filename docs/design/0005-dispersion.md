# 0005 — Dispersive dielectrics

Status: implemented (`solver/src/dispersive.c`, C fit in `cem_fit_djordjevic_sarkar`, Python twin in `tools/dielectric.py`).

## Model

Laminates are specified the way datasheets give them: Dk and Df at one frequency. The Djordjevic–Sarkar (wideband Debye) model turns that into a causal ε(ω):

```
ε(ω) = ε∞ + Δε/(m2 − m1) · log10((ω2 + jω)/(ω1 + jω)),   ω1 = 2π·10^m1, ω2 = 2π·10^m2
```

with m1 = 3, m2 = 12 (1 kHz to 1 THz, the common default). Δε and ε∞ follow from Dk(1 − j·Df) at f_ref.

## Debye poles

The DS model is exactly a log-uniform distribution of Debye relaxations between ω1 and ω2. FDTD needs a finite sum, so N poles are placed log-uniformly over [f_min/10, 10·f_max] and their weights (plus ε∞) are fitted by non-negative least squares over the band, with the imaginary part weighted up because loss is the small, important quantity. Relaxations far outside the band are absorbed into ε∞ (above) or dropped (below, they only act below the band).

Fit quality for RO4350B-like Dk 3.66 / Df 0.0037 @ 10 GHz, 0.1–40 GHz:

| poles | max ε′ error | max tan δ error |
|---|---|---|
| 4 | 0.10% | 28.5% |
| 6 | 0.02% | 4.4% |
| 8 | 0.03% | 0.66% |
| 10 | 0.006% | 0.10% |

Default: 8 poles. All dispersive materials share one pole set (frequencies); each material has its own weights, so interface edges average weights exactly like ε and σ.

## Time stepping (ADE)

Per pole, polarisation current `τ dJ/dt + J = ε0 Δε dE/dt`, trapezoidal at n+½:

```
J^{n+1} = k J^n + β (E^{n+1} − E^n)/Δt,   k = (1 − Δt/2τ)/(1 + Δt/2τ),   β = ε0 Δε Δt / (τ (1 + Δt/2τ))
```

Ampère with J averaged over the step keeps the ordinary update if ε → ε0ε∞ + Σβ/2 in `ca/cb`, followed by `E −= cb Σ (1 + k)/2 · J^n`, then the J update. Storage: one J per pole per dispersive edge. Lumped-port edges stay non-dispersive.

## 2D solver

The mode solver evaluates the same pole sum in the frequency domain, so port Z0(f), β(f) and α(f) match the 3D medium.

## Validation

Lossy GCPW with mode ports: |S21| agrees with exp(−αL) from the 2D mode within 0.0044 dB over 1–40 GHz (0.2 dB total loss at 40 GHz). α_d is within ±5% of the quasi-static filling-factor formula, which itself is approximate.
