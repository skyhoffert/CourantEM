#!/usr/bin/env python3
"""Figure for the milestone-1 GCPW validation note.
Usage: plot_m1_validation.py MODEPORT_RUN LUMPED_RUN OUT.png"""
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "tools")
from gcpw_ref import gcpw_laplace, gcpw_conformal  # noqa: E402

mp, lp, out = sys.argv[1], sys.argv[2], sys.argv[3]
d = np.genfromtxt(f"{mp}/waves.csv", delimiter=",", names=True)
f = d["f_hz"] / 1e9
a1 = d["a1_re"] + 1j * d["a1_im"]; b1 = d["b1_re"] + 1j * d["b1_im"]; b2 = d["b2_re"] + 1j * d["b2_im"]
S11, S21 = b1 / a1, b2 / a1
L = 10e-3
alpha = d["alpha"] if "alpha" in d.dtype.names else 0 * d["beta"]
ideal = np.exp(-(alpha + 1j * d["beta"]) * L)
eps = (d["beta"] * 299792458.0 / (2 * np.pi * d["f_hz"])) ** 2
ls = np.genfromtxt(f"{lp}/sparams.csv", delimiter=",", names=True)
L11 = ls["s11_re"] + 1j * ls["s11_im"]
meta = json.load(open(f"{mp}/mesh.json"))["geometry"]
ref = gcpw_laplace(meta["er"], meta["h"], meta["w"], meta["s"], meta["fence"])["extrapolated"]

plt.rcParams.update({"font.size": 9})
fg, ax = plt.subplots(2, 2, figsize=(10, 6.5))
ax[0, 0].plot(f, d["z0"], color="#2a7de1", label="mode port Z0(f), 2D eigenmode on the FDTD mesh")
ax[0, 0].axhline(ref["z0"], color="k", ls="--", lw=1, label=f"2D Laplace, open, converged: {ref['z0']:.2f} Ω")
ax[0, 0].set_ylabel("Z0 (Ω)")
ax[0, 1].plot(f, eps, color="#2a7de1", label="mode port ε_eff(f)")
ax[0, 1].axhline(ref["eps_eff"], color="k", ls="--", lw=1, label=f"2D Laplace (static): {ref['eps_eff']:.4f}")
ax[0, 1].set_ylabel("ε_eff")
ax[1, 0].plot(f, 20 * np.log10(np.abs(L11)), color="#e8590c", label="|S11|, 50 Ω lumped ports")
ax[1, 0].plot(f, 20 * np.log10(np.abs(S11)), color="#2a7de1", label="|S11|, mode ports")
ax[1, 0].set_ylabel("dB"); ax[1, 0].set_ylim(-90, 0)
ax[1, 1].plot(f, np.angle(S21 / ideal, deg=True), color="#2a7de1", label="∠S21 + βL (mode ports)")
ax[1, 1].plot(f, 20 * np.log10(np.abs(S21)) * 10, color="#7a7a7a", label="|S21| (dB) × 10")
ax[1, 1].set_ylabel("deg / (dB×10)"); ax[1, 1].set_ylim(-0.3, 0.3)
for x in ax.flat:
    x.set_xlabel("GHz"); x.grid(alpha=0.3); x.legend(fontsize=7.5)
fg.suptitle("Milestone 1, uniform GCPW (εr 3.66, h 254 µm, w 500 µm, s 200 µm), 10 mm, hmin 20 µm, open (CPML)")
fg.tight_layout()
fg.savefig(out, dpi=120)
print("saved", out)
