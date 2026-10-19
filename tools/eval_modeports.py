#!/usr/bin/env python3
"""Evaluate a mode-port run of a uniform line: S11 should vanish and S21 should
equal exp(-j beta L). Usage: eval_modeports.py RUN_DIR [--len m]"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--len", type=float, default=10e-3)
ap.add_argument("--plot")
a = ap.parse_args()
d = np.genfromtxt(f"{a.run}/waves.csv", delimiter=",", names=True)
f = d["f_hz"]
a1 = d["a1_re"] + 1j * d["a1_im"]; b1 = d["b1_re"] + 1j * d["b1_im"]
a2 = d["a2_re"] + 1j * d["a2_im"]; b2 = d["b2_re"] + 1j * d["b2_im"]
S11, S21 = b1 / a1, b2 / a1
alpha = d["alpha"] if "alpha" in d.dtype.names else 0 * d["beta"]
ideal = np.exp(-(alpha + 1j * d["beta"]) * a.len)
band = (f >= 1e9) & (f <= 40e9)
for fq in [1, 2, 5, 10, 20, 30, 40]:
    q = np.argmin(abs(f - fq * 1e9))
    print(f"{f[q]/1e9:5.1f} GHz |S11| {20*np.log10(abs(S11[q])):7.1f} dB  |S21| {20*np.log10(abs(S21[q])):+.4f} dB  "
          f"phase err {np.angle(S21[q]/ideal[q], deg=True):+.3f} deg  |a2/a1| {20*np.log10(abs(a2[q]/a1[q])):6.1f} dB  "
          f"Z0 {d['z0'][q]:.3f}")
print("|S21| vs exp(-alpha L) max dev (dB):", f"{np.max(np.abs(20*np.log10(np.abs(S21[band]/ideal[band])))):.4f}")
print(f"max |S11| 1-40 GHz: {20*np.log10(np.max(abs(S11[band]))):.1f} dB; "
      f"max ||S21|-1|: {np.max(abs(20*np.log10(abs(S21[band])))):.4f} dB; "
      f"max phase err: {np.max(abs(np.angle(S21[band]/ideal[band], deg=True))):.3f} deg")
if a.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fg, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(f / 1e9, 20 * np.log10(abs(S11)), label="|S11|")
    ax[0].plot(f / 1e9, 20 * np.log10(abs(S21)), label="|S21|")
    ax[0].set_ylim(-90, 5); ax[0].set_ylabel("dB"); ax[0].legend()
    ax[1].plot(f / 1e9, np.angle(S21 / ideal, deg=True))
    ax[1].set_ylabel("∠S21 − (−βL) (deg)")
    for x in ax:
        x.set_xlabel("GHz"); x.grid(alpha=.3)
    fg.tight_layout(); fg.savefig(a.plot, dpi=110)
