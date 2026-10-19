#!/usr/bin/env python3
"""Wideband dielectric models for Courant EM.

Djordjevic–Sarkar (wideband Debye) model, set from Dk/Df at one frequency:
    eps(w) = eps_inf + d_eps / (m2 - m1) * log10((w2 + j w) / (w1 + j w)),
    w1 = 2 pi 10^m1, w2 = 2 pi 10^m2   (defaults: 1 kHz .. 1 THz)

FDTD needs a sum of Debye poles. The DS model is exactly a log-uniform
distribution of Debye relaxations between w1 and w2, so poles are placed
log-uniformly across the band of interest and their weights are fitted by
non-negative least squares; relaxations far above the band fold into eps_inf
and those far below are dropped (they only matter below the band).

Usage:
    python3 tools/dielectric.py --dk 3.66 --df 0.0037 --fref 10e9 --fmin 0.1e9 --fmax 40e9 --poles 6 [--out poles.txt]
Output file format (read by gcpw_line --debye): first line eps_inf, then one
line per pole: "f_p d_eps" (f_p in Hz, tau = 1/(2 pi f_p)).
"""
import argparse
import math

import numpy as np
from scipy.optimize import nnls


def ds_eps(f, dk, df, fref, m1=3.0, m2=12.0):
    """Complex relative permittivity (eps' - j eps'') of the DS model."""
    w1, w2 = 2 * math.pi * 10 ** m1, 2 * math.pi * 10 ** m2

    def L(fq):
        w = 2 * math.pi * fq
        return np.log10((w2 + 1j * w) / (w1 + 1j * w)) / (m2 - m1)

    Lr = L(fref)
    target = dk * (1 - 1j * df)
    d_eps = target.imag / Lr.imag
    eps_inf = dk - d_eps * Lr.real
    return eps_inf + d_eps * L(np.asarray(f, float)), eps_inf, d_eps


def fit_debye(dk, df, fref, fmin, fmax, npoles, m1=3.0, m2=12.0):
    f = np.geomspace(fmin, fmax, 400)
    eps, _, _ = ds_eps(f, dk, df, fref, m1, m2)
    # poles spread one decade beyond each end of the band
    fp = np.geomspace(fmin / 10, fmax * 10, npoles)
    w = 2 * math.pi * f
    M = np.zeros((2 * len(f), npoles + 1))
    for p, fq in enumerate(fp):
        r = 1 / (1 + 1j * f / fq)
        M[:len(f), p] = r.real
        M[len(f):, p] = r.imag
    M[:len(f), -1] = 1.0   # eps_inf (real only)
    # weight the imaginary part up: loss is the small, important quantity
    wi = 1.0 / max(dk * df, 1e-6)
    wr = 1.0 / dk
    A = np.vstack([M[:len(f)] * wr, M[len(f):] * wi])
    b = np.r_[eps.real * wr, eps.imag * wi]
    x, _ = nnls(A, b)
    fit = x[-1] + sum(x[p] / (1 + 1j * f / fp[p]) for p in range(npoles))
    err_re = np.max(np.abs(fit.real - eps.real) / eps.real)
    err_tan = np.max(np.abs((-fit.imag / fit.real) - (-eps.imag / eps.real)) / (-eps.imag / eps.real))
    return x[-1], fp, x[:-1], err_re, err_tan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dk", type=float, default=3.66)
    ap.add_argument("--df", type=float, default=0.0037)
    ap.add_argument("--fref", type=float, default=10e9)
    ap.add_argument("--fmin", type=float, default=0.1e9)
    ap.add_argument("--fmax", type=float, default=40e9)
    ap.add_argument("--poles", type=int, default=6)
    ap.add_argument("--out")
    a = ap.parse_args()
    einf, fp, de, er, et = fit_debye(a.dk, a.df, a.fref, a.fmin, a.fmax, a.poles)
    print(f"eps_inf {einf:.6f}; max error over band: eps' {er*100:.4f}%, tan d {et*100:.2f}%")
    for f, d in zip(fp, de):
        print(f"  pole {f/1e9:10.4f} GHz  d_eps {d:.6f}")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(f"{einf:.12g}\n")
            for f, d in zip(fp, de):
                fh.write(f"{f:.12g} {d:.12g}\n")


if __name__ == "__main__":
    main()
