# 0001 — FDTD kernel (double-precision reference)

Status: implemented (`solver/src/sim.c`), 2026-10-17.

## Scheme

Standard Yee leapfrog on a non-uniform Cartesian grid. For a cell `(i, j, k)` spanning `[x_i, x_{i+1}] × [y_j, y_{j+1}] × [z_k, z_{k+1}]`:

```
Ex(i+½, j, k)   Ey(i, j+½, k)   Ez(i, j, k+½)
Hx(i, j+½, k+½) Hy(i+½, j, k+½) Hz(i+½, j+½, k)
```

All six components share one `(nx+1)(ny+1)(nz+1)` index space (k fastest), which keeps index arithmetic identical for every component and maps directly to a GPU buffer layout later.

- H update uses primary spacings `d[i] = x_{i+1} − x_i`; E update uses dual spacings `dd[i] = (d[i−1] + d[i]) / 2`.
- Per-edge coefficients `ca = (1 − h)/(1 + h)`, `cb = (Δt/ε)/(1 + h)`, `h = σΔt/(2ε)`. ε and σ are area-weighted averages of the four cells that share the edge.
- Time step `Δt = CFL / (c √(Δx_min⁻² + Δy_min⁻² + Δz_min⁻²))`, CFL 0.99.
- Tangential E on the outer faces is never updated (PEC behind everything).

**Why non-uniform Yee and not a subgridded or conformal scheme first?** It is second-order accurate when neighbouring cells change by ≤ 1.3–1.5×, it's simple to port to WGSL, and PCB geometry is mostly Manhattan. Conformal (Dey–Mittra) cells are planned for round vias and angled traces.

## Materials and metal

- `cem_material_box`: cells whose centre lies in a box.
- `cem_pec_box`: solid metal; every edge of a PEC cell is PEC.
- `cem_pec_sheet`: zero-thickness metal on a node plane; every tangential edge lying fully inside the rectangle is PEC. With the thirds-rule mesh (note 0002) this puts the effective metal edge where the field singularity wants it.

## CPML

Roden–Gedney CPML with CFS parameters, polynomial grading m = 3:

```
σ(ρ) = σ_max ρ^m,  κ(ρ) = 1 + (κ_max − 1) ρ^m,  α(ρ) = α_max (1 − ρ)
σ_max = 0.8 (m+1) / (η0 Δ),  κ_max = 5,  α_max = 0.05 S/m
b = exp(−(σ/κ + α) Δt/ε0),  c = σ (b − 1) / (σκ + κ²α)
```

The layer cells are part of the grid. Main updates use `1/(κ Δ)`; ψ auxiliaries are updated only inside the slabs. Measured reflection for a 10-cell layer: −61 dB (test_cpml).

The 2D mode solver reuses the same σ, κ, α as a frequency-domain stretch `s = κ + σ/(α + jωε0)`, so port cross-sections are terminated the same way as the 3D domain.

## Lumped ports

A column of edges along one axis. Each edge carries a share `R·Δℓ/L` of the port resistance and `Vs·Δℓ/L` of the source, updated semi-implicitly:

```
E' = dir·E,   ε dE'/dt = dir·(∇×H) − (E'Δℓ − Vs_e)/(R_e A)
β = Δt Δℓ / (2 R_e A ε)
E'^{n+1} = (1−β)/(1+β) E'^n + (Δt/ε)/(1+β) dir·(∇×H) + Δt Vs_e /(R_e A ε (1+β))
```

Port voltage and current are sampled at `t = (n+½)Δt`; `a = (V + Z I)/2√Z`, `b = (V − Z I)/2√Z`. Lumped ports are kept for component-like locations; for line ports they add 0.3 mm of vertical-column inductance (S11 up to −6 dB at 35 GHz on the test line), which is why mode ports exist (note 0003).

## Probes and DFT

Voltage probes integrate E along an edge column; current probes take the circulation of H around a dual-grid rectangle. All spectra are running DFTs evaluated at the exact sample times (E at `nΔt`, H at `(n+½)Δt`), with phasors advanced by complex rotation and re-synchronised with `cos/sin` every 1024 steps.

## Run control

`cem_run` stops after the source has finished when total field energy falls below `decay × peak` (checked every N steps).

## Verification

- `test_cavity`: PEC cavity TE101/TE201/TE102. Uniform grid matches the Yee discrete dispersion relation to ≤ 1e-5; a ±16% non-uniform grid is within 0.1% of the analytic frequency.
- `test_cpml`: dipole near a 10-cell CPML vs a 128³ reference box: −61 dB.
- `parity`: native vs WebAssembly ≤ 1e-14 relative.

## Rejected alternatives

- **Single precision in the reference path**: the C path is the f64 reference by plan; f32 is the GPU's job.
- **Per-slab ψ storage**: memory-optimal but more index bookkeeping; full-size ψ arrays are allocated for now and only slabs are visited. Revisit for the GPU layout.
- **Separate arrays per component size**: cleaner memory use but every kernel would need its own strides; one index space wins for portability.
