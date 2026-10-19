# 0002 — Meshing zero-thickness metal edges: the thirds rule

Status: adopted for the GCPW scenario (`gcpw_line.c`, `web/src/scenarios/gcpw.ts`); to be built into the rule-based mesher (build step 5). **This changes a plan item; Sky to confirm.**

## Problem

The plan's mesh rule 1 says grid lines snap to every copper edge. With zero-thickness sheets, the field near the edge is singular (∝ r^−½). FDTD samples tangential E on edges, and a sheet whose last PEC edge sits on the physical edge behaves as if it were wider by roughly a third of a cell. Capacitance goes up and Z0 goes down, and the error shrinks only slowly with refinement (~h^0.5–1).

Measured on the default GCPW (2D eigenmode on the FDTD's own cross-section, quasi-static, against a converged Laplace reference in the same box):

| min cell | snap to edge | thirds rule |
|---|---|---|
| 20 µm | −2.46% | −0.18% |
| 10 µm | −1.37% | −0.23% |
| 5 µm | −0.67% | −0.09% |
| 2.5 µm | −0.34% | −0.06% |

The 3D FDTD showed the same pattern: with snapping, Z0 at 2–5 GHz was −3.3%, −2.2%, −1.3% at 20, 10, 5 µm.

## Rule

At each zero-thickness metal edge, place one grid line `h/3` inside the metal and one `2h/3` outside it, where `h` is the local minimum cell. The metal in the FDTD then ends at the inner line. (This is the "thirds rule" known from openEMS practice.) Grading continues from these lines with ratio ≤ 1.3.

No staircasing is introduced: rectangular features still align with grid lines; only the line positions near edges move.

## Consequences

- The mesher needs to know which side of each edge is metal. For Gerber/KiCad copper polygons this comes from polygon orientation.
- Two lines per edge instead of one: ~3% more cells across a CPW.
- Thick copper: the edges are 270° corners with a weaker singularity (r^−⅓), and snapping lines to both faces is the usual choice. Measured on the imported GCPW with 35 µm copper (Laplace reference with the same bars: 46.761 Ω, `tools/gcpw_ref.py --t 35e-6`): lines on the side faces give −1.9% / −0.9% at 20 / 10 µm cells; the thirds rule on the side edges (bar ends h/3 inside the drawn edge) gives +0.48% / +0.25%. So the thirds rule is kept for thick copper's side edges too; the top face gets a plain grid line. Note that copper thickness itself moves Z0 by −6% on this line (49.92 → 46.76 Ω), so boards are modelled with their stackup copper thickness by default.

## Rejected alternatives

- **Finer cells only**: 4× more cells per halving and ~8× more steps; too slow in a browser.
- **Edge-singularity-corrected update coefficients** (modified ε/μ near the edge): more accurate in principle, but harder to keep consistent between the C, WGSL and 2D mode-solver paths. Possible later.
