# Step 9 protocol: openEMS, HFSS and measured comparisons

What to run, with which settings, and how to judge it against the plan's tolerances. Everything here can be done on Sky's machine; nothing in it has run against openEMS, HFSS or a VNA yet.

## Tolerances (plan, "Validation")

| Comparison | Range | Tolerance | `compare_s2p.py --rule` |
|---|---|---|---|
| GPU (f32) vs CPU (f64) | to 40 GHz | \|S21\| ±0.01 dB, phase ±0.1° | `gpu` |
| Analytical, uniform GCPW | to 40 GHz | Z0 and ε_eff ±1%, attenuation ±10% | (in `validation/gcpw-m1.md`) |
| openEMS, same geometry | to 40 GHz | \|S21\| ±0.05 dB, phase ±1°; S11 resonances ±1%, magnitude ±1 dB above −30 dB | `openems` |
| HFSS | to 40 GHz | \|S21\| ±0.15 dB, phase ±3°; same S11 rule | `hfss` |
| Measured, TRL-corrected | to 20 GHz | \|S21\| ±0.3 dB | `measured` |

Known gap before starting: conductor attenuation on edge-dominated lines comes out ~25–33% low (issue 9), which fails the analytical ±10% attenuation rule and will show in lossy |S21| against HFSS and measurement on long lines. Lossless comparisons and Z0/phase are not affected.

## 1. Courant EM results

Either the app (load the board, set the stackup, move the port planes past the launches, tick "Full 2-port" for anything not mirror-symmetric, Run, Save .s2p) or the command line, which gives the same numbers:

```sh
make -C solver && make -C solver wasm
node web/tools/run-board.ts board.kicad_pcb --net RF --planes 1.5,1.5 --hmin 20 --fmax 40 --df 0.25 --lossy --full2 --out courant.s2p
# Gerbers: list the copper, outline and drill files and pass the stackup as JSON
node web/tools/run-board.ts *.gbr *.gm1 *.drl --stackup stack.json --net Net-1 --planes 1.5,1.5 --out courant.s2p
```

Every .s2p is stamped with the solver version (`! Courant EM 20261017.0 …`). Keep the file, the board files and the command line together in `validation/`.

Mesh convergence: run at `--hmin 20` and `--hmin 10`; the comparison is meaningful when the two agree well inside the tolerance being tested.

## 2. openEMS

`tools/openems_gcpw.py` builds the test GCPW in openEMS with the same geometry and mesh rules (thirds rule, ratio 1.3, λ/20) and MSL-type ports, and writes `gcpw.s2p`. It was written without openEMS available: check the port setup on first use.

```sh
python3 tools/openems_gcpw.py --out runs/openems --len 10 --t 0      # sheets
python3 tools/openems_gcpw.py --out runs/openems-t35 --len 10 --t 35 # 35 um copper
./solver/build/gcpw_line --modeports --hmin 20 --dx 50 --df 0.5 --out runs/courant   # writes line_50ohm.s2p
python3 tools/compare_s2p.py runs/courant/line_50ohm.s2p runs/openems/gcpw.s2p --rule openems --plot runs/openems/compare.png
```

Compare lossless first (geometry and ports), then with the same Debye/conductor models if openEMS's are set up equivalently.

## 3. HFSS

Model the same board (import the KiCad/ODB++/Gerbers or draw the line), with:

- Wave ports at the same reference planes as Courant's port planes, 50 Ω renormalisation, de-embedding to the planes.
- Substrate: Djordjevic–Sarkar with Dk/Df at 10 GHz as in the stackup editor.
- Copper: finite conductivity 5.8e7 S/m, thickness from the stackup, the same roughness model (Groisse, Rq) if roughness is on.
- Solder mask if present: Courant models it as a slab filling the gaps (not conformal over copper); either match that in HFSS or note the difference.
- Radiation boundary at least λ/4 away at the lowest frequency; interpolating sweep 0.25–40 GHz.

Export Touchstone (RI or MA, 50 Ω) and run `compare_s2p.py --rule hfss`.

## 4. Measured, TRL-corrected

Lines of several lengths with SMA (or 2.92 mm) edge launches and TRL standards on each panel; de-embed to on-board reference planes. In Courant, move the port planes to the TRL reference planes (the "port plane, mm in" inputs or `--planes`), so the launches are cut away exactly as the calibration removes them. Run `compare_s2p.py --rule measured` (to 20 GHz) and plot against the measurement uncertainty.

## 5. Recording

For each comparison, add a row to `validation/gcpw-m1.md` (or a new note per board) and to the tracker's validation table: date, solver version, files, rule, worst numbers, pass/fail. Raw and corrected measurements, board files and the exact command lines go in the repo.
