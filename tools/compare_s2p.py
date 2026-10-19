#!/usr/bin/env python3
"""Compare two Touchstone 2-port files against the milestone-1 validation rules.

    python3 tools/compare_s2p.py courant.s2p reference.s2p --rule openems [--plot out.png] [--fmax 40e9]

Rules (plan, "Validation"):
  openems   |S21| +-0.05 dB, phase +-1 deg; S11 resonances +-1 %, |S11| +-1 dB where above -30 dB
  hfss      |S21| +-0.15 dB, phase +-3 deg; same S11 rule
  measured  |S21| +-0.3 dB (to 20 GHz); S11 reported only
  gpu       |S21| +-0.01 dB, phase +-0.1 deg

Both files are read in any Touchstone v1 format (RI/MA/DB, Hz..GHz, any R) and
renormalised to a common 50 ohm reference if their R differs. The comparison
runs on the first file's frequencies inside the overlap, with the second file
interpolated (magnitude in dB and unwrapped phase, linearly). S11 "resonances"
are the dips of |S11| below -10 dB; each dip of the first file is matched to the
nearest dip of the second.
"""
import argparse
import json
import math
import sys

import numpy as np

RULES = {
    "openems": {"s21_db": 0.05, "s21_deg": 1.0, "s11_res": 0.01, "s11_db": 1.0, "fmax": 40e9},
    "hfss": {"s21_db": 0.15, "s21_deg": 3.0, "s11_res": 0.01, "s11_db": 1.0, "fmax": 40e9},
    "measured": {"s21_db": 0.3, "s21_deg": None, "s11_res": None, "s11_db": None, "fmax": 20e9},
    "gpu": {"s21_db": 0.01, "s21_deg": 0.1, "s11_res": None, "s11_db": None, "fmax": 40e9},
}
UNITS = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}


def read_s2p(path):
    """Return f (Hz), S (n, 2, 2) complex, reference R."""
    unit, fmt, R = 1e9, "ma", 50.0
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.split("!", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                tok = line[1:].lower().split()
                for i, t in enumerate(tok):
                    if t in UNITS:
                        unit = UNITS[t]
                    elif t in ("ri", "ma", "db"):
                        fmt = t
                    elif t == "r" and i + 1 < len(tok):
                        R = float(tok[i + 1])
                continue
            rows.extend(float(v) for v in line.split())
    data = np.array(rows).reshape(-1, 9)
    f = data[:, 0] * unit
    pairs = data[:, 1:].reshape(-1, 4, 2)
    if fmt == "ri":
        s = pairs[..., 0] + 1j * pairs[..., 1]
    else:
        mag = 10 ** (pairs[..., 0] / 20) if fmt == "db" else pairs[..., 0]
        s = mag * np.exp(1j * np.radians(pairs[..., 1]))
    # Touchstone 2-port order: S11 S21 S12 S22
    S = np.empty((len(f), 2, 2), complex)
    S[:, 0, 0], S[:, 1, 0], S[:, 0, 1], S[:, 1, 1] = s[:, 0], s[:, 1], s[:, 2], s[:, 3]
    return f, S, R


def renorm(S, R, R0=50.0):
    if abs(R - R0) < 1e-12:
        return S
    I = np.eye(2)
    out = np.empty_like(S)
    for k in range(len(S)):
        Z = R * (I + S[k]) @ np.linalg.inv(I - S[k])
        out[k] = (Z - R0 * I) @ np.linalg.inv(Z + R0 * I)
    return out


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-30))


def interp(f_to, f, x):
    """Interpolate a complex response: dB magnitude and unwrapped phase, linearly."""
    m = np.interp(f_to, f, db(x))
    p = np.interp(f_to, f, np.unwrap(np.angle(x)))
    return 10 ** (m / 20) * np.exp(1j * p)


def dips(f, s11, floor=-10.0):
    d = db(s11)
    return [f[i] for i in range(1, len(f) - 1) if d[i] < floor and d[i] <= d[i - 1] and d[i] <= d[i + 1]]


def compare(a, b, rule, fmax=None):
    fa, Sa, Ra = a
    fb, Sb, Rb = b
    Sa, Sb = renorm(Sa, Ra), renorm(Sb, Rb)
    lim = RULES[rule]
    top = min(fa.max(), fb.max(), fmax or lim["fmax"])
    sel = (fa >= fb.min()) & (fa <= top)
    f = fa[sel]
    out = {"rule": rule, "f_min": float(f.min()), "f_max": float(f.max()), "points": int(len(f)), "checks": []}

    def add(name, value, tol, unit):
        ok = None if tol is None else bool(value <= tol)
        out["checks"].append({"name": name, "value": value, "tol": tol, "unit": unit, "pass": ok})

    for name, (r, c) in (("S21", (1, 0)), ("S12", (0, 1))):
        xa, xb = Sa[sel, r, c], interp(f, fb, Sb[:, r, c])
        dmag = float(np.max(np.abs(db(xa) - db(xb))))
        dph = np.degrees(np.angle(xa / xb))
        add(f"|{name}| max difference", dmag, lim["s21_db"], "dB")
        add(f"phase {name} max difference", float(np.max(np.abs(dph))), lim["s21_deg"], "deg")
    for name, k in (("S11", 0), ("S22", 1)):
        xa, xb = Sa[sel, k, k], interp(f, fb, Sb[:, k, k])
        above = (db(xa) > -30) | (db(xb) > -30)
        dmag = float(np.max(np.abs(db(xa) - db(xb))[above])) if above.any() else 0.0
        add(f"|{name}| max difference where above -30 dB", dmag, lim["s11_db"], "dB")
        da, dbb = dips(f, xa), dips(f, xb)
        if da and dbb:
            worst = max(min(abs(x - y) for y in dbb) / x for x in da)
            add(f"{name} dip frequencies (nearest match)", float(worst), lim["s11_res"], "relative")
    out["pass"] = all(c["pass"] is not False for c in out["checks"])
    return out, (f, Sa[sel], np.array([[interp(f, fb, Sb[:, r, c]) for c in range(2)] for r in range(2)]).transpose(2, 0, 1))


def plot(path, f, Sa, Sb, la, lb):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    g = f / 1e9
    for (r, c), a in zip(((1, 0), (0, 0)), (ax[0, 0], ax[0, 1])):
        a.plot(g, db(Sa[:, r, c]), label=la)
        a.plot(g, db(Sb[:, r, c]), "--", label=lb)
        a.set_title(f"|S{r + 1}{c + 1}| (dB)")
    ax[1, 0].plot(g, db(Sa[:, 1, 0]) - db(Sb[:, 1, 0]))
    ax[1, 0].set_title("|S21| difference (dB)")
    ax[1, 1].plot(g, np.degrees(np.angle(Sa[:, 1, 0] / Sb[:, 1, 0])))
    ax[1, 1].set_title("∠S21 difference (deg)")
    for a in ax.flat:
        a.grid(alpha=0.3)
        a.set_xlabel("GHz")
    ax[0, 0].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ours")
    ap.add_argument("reference")
    ap.add_argument("--rule", choices=sorted(RULES), default="openems")
    ap.add_argument("--fmax", type=float)
    ap.add_argument("--plot")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    res, (f, Sa, Sb) = compare(read_s2p(a.ours), read_s2p(a.reference), a.rule, a.fmax)
    if a.plot:
        plot(a.plot, f, Sa, Sb, "Courant EM", a.reference)
    if a.json:
        print(json.dumps(res, indent=2))
    else:
        print(f"{a.ours} vs {a.reference}: rule {a.rule}, {res['points']} points, "
              f"{res['f_min'] / 1e9:.3g}-{res['f_max'] / 1e9:.3g} GHz")
        for c in res["checks"]:
            tol = "" if c["tol"] is None else f" (tol {c['tol']:g})"
            mark = {True: "pass", False: "FAIL", None: "info"}[c["pass"]]
            print(f"  {mark:4}  {c['name']}: {c['value']:.4g} {c['unit']}{tol}")
        print("PASS" if res["pass"] else "FAIL")
    sys.exit(0 if res["pass"] else 1)


if __name__ == "__main__":
    main()
