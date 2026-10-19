#!/usr/bin/env python3
"""Post-process a gcpw_line run: extract Z0 and eps_eff from the line probes and
compare S-parameters with an ideal transmission-line model.

gamma from three equally spaced voltage probes (valid for any mix of forward and
backward waves on a uniform line):
    cosh(gamma d) = (V(x-d) + V(x+d)) / (2 V(x))
Z0 from two adjacent voltage probes and the current probe midway between them:
    V(x+dx/2) - V(x-dx/2) = -2 Z0 sinh(gamma dx/2) I(x)

Usage: analyze_gcpw.py RUN_DIR [--plot out.png] [--json out.json]
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from gcpw_ref import gcpw_conformal, gcpw_laplace  # noqa: E402

C0 = 299792458.0


def load(run):
    with open(os.path.join(run, "meta.json")) as fh:
        meta = json.load(fh)
    with open(os.path.join(run, "probes.csv")) as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        data = np.array([[float(v) for v in row] for row in rd])
    f = data[:, 0]
    V, I = {}, {}
    for c in range(1, len(hdr), 2):
        name = hdr[c][:-3]
        val = data[:, c] + 1j * data[:, c + 1]
        if name.startswith("V"):
            V[int(name[1:])] = val
        else:
            I[int(name[1:])] = val
    with open(os.path.join(run, "sparams.csv")) as fh:
        rd = csv.reader(fh)
        next(rd)
        sd = np.array([[float(v) for v in row] for row in rd])
    S11 = sd[:, 1] + 1j * sd[:, 2]
    S21 = sd[:, 3] + 1j * sd[:, 4]
    return meta, f, V, I, S11, S21


def extract(meta, f, V, I):
    """Fit V(x) = A e^{-gx} + B e^{gx} over all voltage probes (least squares in
    gamma, linear in A, B), then I(x) = (A e^{-gx} - B e^{gx}) / Z0 over all
    current probes. The three-probe cosh formula seeds the fit."""
    from scipy.optimize import least_squares
    xs = {int(k): v for k, v in meta["x_nodes_for_probes"].items()}
    base = sorted(I.keys())
    dx = meta["dx"]
    vnodes = sorted(V.keys())
    xv = np.array([xs[m] for m in vnodes])
    xv = xv - xv.mean()
    xi = np.array([0.5 * (xs[m] + xs[m + 1]) for m in base]) - (np.array([xs[m] for m in vnodes]).mean())
    Vm = np.array([V[m] for m in vnodes])        # (nprobe, nf)
    Im = np.array([I[m] for m in base])
    # seed from the widest equally spaced triplet
    m0, m1, m2 = base[0], base[len(base) // 2], base[2 * (len(base) // 2)]
    d = xs[m1] - xs[m0]
    g0 = np.arccosh(((V[m0] + V[m2]) / (2 * V[m1])).astype(complex)) / d
    g0 = np.where(g0.imag < 0, -g0, g0)
    nf = len(f)
    gamma = np.zeros(nf, complex)
    Z0 = np.zeros(nf, complex)
    for q in range(nf):
        v = Vm[:, q]

        def basis(g, x):
            return np.stack([np.exp(-g * x), np.exp(g * x)], axis=1)

        def resid(p):
            g = p[0] + 1j * p[1]
            M = basis(g, xv)
            c, *_ = np.linalg.lstsq(M, v, rcond=None)
            r = M @ c - v
            return np.r_[r.real, r.imag] / (np.abs(v).max() + 1e-300)

        beta_est = 2 * math.pi * f[q] * math.sqrt(2.7) / C0
        start = g0[q] if np.isfinite(g0[q]) and g0[q].imag > 0 else 1j * beta_est
        sol = least_squares(resid, [start.real, start.imag], x_scale=[beta_est, beta_est])
        g = sol.x[0] + 1j * sol.x[1]
        if g.imag < 0:
            g = -g
        M = basis(g, xv)
        c, *_ = np.linalg.lstsq(M, v, rcond=None)
        A_, B_ = c
        u = A_ * np.exp(-g * xi) - B_ * np.exp(g * xi)   # = Z0 * I
        y0 = np.vdot(u, Im[:, q]) / np.vdot(u, u)        # least-squares 1/Z0
        gamma[q] = g
        Z0[q] = 1.0 / y0
    beta = gamma.imag
    eps_eff = (beta * C0 / (2 * math.pi * f)) ** 2
    return gamma, eps_eff, Z0, None


def tl_sparams(f, z0, eps_eff, L, zr=50.0):
    g = 1j * 2 * math.pi * f * np.sqrt(eps_eff) / C0
    D = 2 * z0 * zr * np.cosh(g * L) + (z0 ** 2 + zr ** 2) * np.sinh(g * L)
    return (z0 ** 2 - zr ** 2) * np.sinh(g * L) / D, 2 * z0 * zr / D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--plot")
    ap.add_argument("--json")
    ap.add_argument("--no-laplace", action="store_true")
    ap.add_argument("--eigen", action="store_true",
                    help="also solve the 2D mode on the same mesh (closed box) and compare")
    a = ap.parse_args()
    meta, f, V, I, S11, S21 = load(a.run)
    gamma, eps_eff, Z0, zall = extract(meta, f, V, I)
    er, h, w, s, fence = meta["er"], meta["h"], meta["w"], meta["s"], meta["fence"]
    cm = gcpw_conformal(er, h, w, s)
    if a.no_laplace:
        ref = cm
    else:
        ref = gcpw_laplace(er, h, w, s, fence)["extrapolated"]
    lowband = (f >= 2e9) & (f <= 5e9)
    z0_lo = float(np.mean(Z0[lowband].real))
    ee_lo = float(np.mean(eps_eff[lowband]))
    band = (f >= 2e9) & (f <= 40e9)
    out = {
        "run": a.run,
        "ref_laplace": ref,
        "ref_conformal": cm,
        "z0_2to5GHz": z0_lo,
        "eps_eff_2to5GHz": ee_lo,
        "z0_err_pct": 100 * (z0_lo - ref["z0"]) / ref["z0"],
        "eps_eff_err_pct": 100 * (ee_lo - ref["eps_eff"]) / ref["eps_eff"],
        "z0_max_dev_to40_pct": float(100 * np.max(np.abs(Z0[band].real - ref["z0"])) / ref["z0"]),
        "eps_eff_max_dev_to40_pct": float(100 * np.max(np.abs(eps_eff[band] - ref["eps_eff"])) / ref["eps_eff"]),
        "atten_max_db_per_m": float(np.max(np.abs(gamma.real[band])) * 8.686),
    }
    L = meta["x_port2"] - meta["x_port1"]
    s11m, s21m = tl_sparams(f, Z0.real, eps_eff, L)
    s11r, s21r = tl_sparams(f, ref["z0"], ref["eps_eff"], L)
    out["s21_vs_tl_extracted_max_db"] = float(np.max(np.abs(20 * np.log10(np.abs(S21[band])) - 20 * np.log10(np.abs(s21m[band])))))
    out["s21_phase_vs_ref_tl_max_deg"] = float(np.max(np.abs(np.angle(S21[band] / s21r[band], deg=True))))
    eig = None
    if a.eigen:
        from modesolver2d import gcpw_modes
        with open(os.path.join(a.run, "mesh.json")) as fh:
            mesh = json.load(fh)
        fe = f[(f >= 2e9) & (f <= 40e9)][::4]
        res, _ = gcpw_modes(mesh, list(fe))
        e_eps = np.array([r["eps_eff"] for r in res])
        e_z = np.array([r["z_vi"].real for r in res])
        idx = [int(np.argmin(np.abs(f - x))) for x in fe]
        out["eigen_eps_eff_max_dev_pct"] = float(np.max(np.abs(eps_eff[idx] - e_eps) / e_eps) * 100)
        out["eigen_z0_max_dev_pct"] = float(np.max(np.abs(Z0[idx].real - e_z) / e_z) * 100)
        eig = (fe, e_eps, e_z)
    print(json.dumps(out, indent=2))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fg, ax = plt.subplots(2, 2, figsize=(11, 7.5))
        fG = f / 1e9
        ax[0, 0].plot(fG, Z0.real, label="FDTD (probe extraction)")
        ax[0, 0].axhline(ref["z0"], color="k", ls="--", label=f"2D Laplace {ref['z0']:.2f} Ω")
        ax[0, 0].axhline(cm["z0"], color="gray", ls=":", label=f"conformal (no fence) {cm['z0']:.2f} Ω")
        ax[0, 0].set_ylabel("Z0 (Ω)"); ax[0, 0].legend(fontsize=8)
        ax[0, 1].plot(fG, eps_eff, label="FDTD")
        ax[0, 1].axhline(ref["eps_eff"], color="k", ls="--", label=f"2D Laplace {ref['eps_eff']:.4f}")
        ax[0, 1].axhline(cm["eps_eff"], color="gray", ls=":", label=f"conformal {cm['eps_eff']:.4f}")
        if eig is not None:
            ax[0, 0].plot(eig[0] / 1e9, eig[2], "o", ms=3, mfc="none", label="2D eigenmode, same mesh")
            ax[0, 1].plot(eig[0] / 1e9, eig[1], "o", ms=3, mfc="none", label="2D eigenmode, same mesh")
            ax[0, 0].legend(fontsize=8)
        ax[0, 1].set_ylabel("ε_eff"); ax[0, 1].legend(fontsize=8)
        ax[1, 0].plot(fG, 20 * np.log10(np.abs(S21)), label="|S21| FDTD")
        ax[1, 0].plot(fG, 20 * np.log10(np.abs(S11)), label="|S11| FDTD")
        ax[1, 0].plot(fG, 20 * np.log10(np.abs(s11r) + 1e-12), "k--", lw=0.8, label="|S11| ideal line (ref)")
        ax[1, 0].set_ylabel("dB"); ax[1, 0].set_ylim(-60, 2); ax[1, 0].legend(fontsize=8)
        ax[1, 1].plot(fG, np.unwrap(np.angle(S21)) * 180 / math.pi, label="∠S21 FDTD")
        ax[1, 1].plot(fG, np.unwrap(np.angle(s21r)) * 180 / math.pi, "k--", lw=0.8, label="∠S21 ideal line (ref)")
        ax[1, 1].set_ylabel("deg"); ax[1, 1].legend(fontsize=8)
        for x in ax.flat:
            x.set_xlabel("GHz"); x.grid(alpha=0.3)
        fg.suptitle(f"GCPW line, L = {L*1e3:.1f} mm, cells {meta['cells']}, hmin {meta['hmin']*1e6:.0f} µm")
        fg.tight_layout()
        fg.savefig(a.plot, dpi=110)


if __name__ == "__main__":
    main()
