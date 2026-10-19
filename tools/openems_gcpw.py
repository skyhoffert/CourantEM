#!/usr/bin/env python3
"""openEMS model of the milestone-1 test GCPW, for the step-9 comparison.

    python3 tools/openems_gcpw.py --out runs/openems [--len 10] [--t 0] [--hmin 20] [--fmax 40e9]
    python3 tools/compare_s2p.py <courant>.s2p runs/openems/gcpw.s2p --rule openems

NOTE: written without openEMS available (it can't be installed in Claude's
sandbox), so this script has not been run. It uses the documented Python API
(CSXCAD / openEMS 0.0.36): check the port setup on first use.

Geometry, matching the built-in Courant EM line (src/scenarios/gcpw.ts) and the
KiCad fixture: εr 3.66 (lossless), h 254 µm, strip 500 µm, gaps 200 µm, via
fences idealised as PEC walls 300 µm beyond each gap edge, bottom ground =
PEC boundary, copper zero-thickness (--t 0) or bars of thickness t. Both ends
run into PML; ports are openEMS MSL ports (voltage strip-to-ground, current
around the strip), which for this quasi-TEM line are equivalent to Courant's
mode ports; S-parameters are referred to 50 Ω and to planes `len` apart.

Mesh: lines on every copper edge plus the thirds rule (h/3 inside, 2h/3
outside), smoothed with ratio 1.3, λ/20 maximum — the same rules as Courant's
mesher, so differences are the solvers', not the meshes'.
"""
import argparse
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/openems")
    ap.add_argument("--len", type=float, default=10.0, help="line length between reference planes (mm)")
    ap.add_argument("--t", type=float, default=0.0, help="copper thickness (um); 0 = sheets")
    ap.add_argument("--hmin", type=float, default=20.0, help="cell at copper edges (um)")
    ap.add_argument("--fmax", type=float, default=40e9)
    ap.add_argument("--er", type=float, default=3.66)
    a = ap.parse_args()

    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    from openEMS.physical_constants import C0

    um = 1e-6
    h, w, s, fence = 254.0, 500.0, 200.0, 300.0           # um
    L = a.len * 1e3                                        # um between reference planes
    feed = 4000.0                                          # line beyond each reference plane
    yb = w / 2 + s                                         # gap outer edge
    yf = yb + fence                                        # inner face of the fence wall
    ywall = yf + 200.0
    air = 1500.0
    t = a.t
    hmin = a.hmin

    FDTD = openEMS(EndCriteria=1e-5)
    f0 = a.fmax / 2
    FDTD.SetGaussExcite(f0, f0)
    # x: PML both ends; y: PML; z: PEC floor (bottom ground), PML on top
    FDTD.SetBoundaryCond(["PML_8", "PML_8", "PML_8", "PML_8", "PEC", "PML_8"])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(um)

    x0, x1 = -L / 2 - feed, L / 2 + feed
    ymax = ywall + 1000.0
    sub = CSX.AddMaterial("substrate", epsilon=a.er)
    sub.AddBox([x0, -ymax, 0], [x1, ymax, h], priority=0)

    pec = CSX.AddMetal("copper")
    ztop = h + t
    # top grounds (both sides) and the fence walls
    for sgn in (1, -1):
        pec.AddBox([x0, sgn * yb, h], [x1, sgn * ymax, ztop], priority=10)
        pec.AddBox([x0, sgn * yf, 0], [x1, sgn * ywall, h], priority=10)

    # mesh: copper edges with the thirds rule, dielectric and copper planes
    ylines = [0.0, ymax, -ymax]
    # thirds rule: one line h/3 inside the metal, one 2h/3 outside (the metal then ends at the inner line)
    for sgn in (1, -1):
        ylines += [sgn * (w / 2 - hmin / 3), sgn * (w / 2 + 2 * hmin / 3)]   # strip edges (metal toward y = 0)
        ylines += [sgn * (yb + hmin / 3), sgn * (yb - 2 * hmin / 3)]         # ground edges (metal outward)
    ylines += [yf, -yf, ywall, -ywall]
    mesh.AddLine("y", ylines)
    zlines = [0.0, h, h + air] + ([ztop] if t > 0 else [])
    zlines += list(np.linspace(0, h, 5))
    mesh.AddLine("z", zlines)
    mesh.AddLine("x", [x0, -L / 2, L / 2, x1])
    lam20 = C0 / a.fmax / np.sqrt(a.er) / 20 / um
    mesh.SmoothMeshLines("y", min(lam20, 200.0), 1.3)
    mesh.SmoothMeshLines("z", min(lam20, 200.0), 1.3)
    mesh.SmoothMeshLines("x", min(lam20, 4 * hmin), 1.3)

    # MSL ports: the strip is the port's metal; propagation along x, excitation along z.
    # The excitation sits 1.5 mm in from the domain end, outside the 8-cell PML.
    ports = []
    ports.append(FDTD.AddMSLPort(1, pec, [x0, -w / 2, ztop], [-L / 2, w / 2, 0], "x", "z",
                                 excite=-1, FeedShift=1500.0, MeasPlaneShift=feed, priority=20))
    ports.append(FDTD.AddMSLPort(2, pec, [x1, -w / 2, ztop], [L / 2, w / 2, 0], "x", "z",
                                 MeasPlaneShift=feed, priority=20))
    # the strip over the whole length (a bar when t > 0; the ports' own strips are sheets at its top face)
    pec.AddBox([x0, -w / 2, h], [x1, w / 2, ztop], priority=20)

    os.makedirs(a.out, exist_ok=True)
    sim_path = os.path.abspath(os.path.join(a.out, "sim"))
    FDTD.Run(sim_path, cleanup=True)

    f = np.linspace(0.5e9, a.fmax, 80)
    for p in ports:
        p.CalcPort(sim_path, f, ref_impedance=50)
    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc
    out = os.path.join(a.out, "gcpw.s2p")
    with open(out, "w") as fh:
        fh.write(f"! openEMS GCPW, len {a.len} mm, t {t} um, hmin {hmin} um; port 1 excited, S22 = S11, S12 = S21 assumed\n")
        fh.write("# HZ S RI R 50\n")
        for k in range(len(f)):
            v = [s11[k], s21[k], s21[k], s11[k]]
            fh.write(" ".join([f"{f[k]:.6e}"] + [f"{x.real:.9e} {x.imag:.9e}" for x in v]) + "\n")
    print("wrote", out)
    print("port 1 Z0 at", f[1] / 1e9, "GHz:", np.real(ports[0].Z_ref) if hasattr(ports[0], "Z_ref") else "n/a")


if __name__ == "__main__":
    main()
