#!/usr/bin/env python3
"""Reference models for a grounded coplanar waveguide (GCPW / conductor-backed CPW).

Two independent references, both quasi-static, zero-thickness, lossless:

1. Conformal mapping (Ghione & Naldi / Simons), infinite lateral ground planes, no
   via fences.  Closed form, good to ~1% for typical geometries.

2. A 2D finite-volume Laplace solver on a graded tensor grid that can include the
   via-fence walls and finite domain used by the FDTD model.  Capacitance is taken
   from the discrete field energy (C = phi^T A phi with phi = 1 on the strip), which
   is variational, then Richardson-extrapolated over a sequence of refinements.

Z0 = 1 / (c * sqrt(C * C0)),  eps_eff = C / C0,
where C is per-unit-length capacitance with the dielectric and C0 with air only.

Usage:
    python3 tools/gcpw_ref.py --er 3.66 --h 254e-6 --w 500e-6 --s 200e-6 --fence 300e-6
"""
import argparse
import json
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import ellipk

C0 = 299792458.0
EPS0 = 8.8541878128e-12


def _kratio(k):
    """K(k)/K(k') using scipy's parameter m = k^2."""
    return ellipk(k * k) / ellipk(1.0 - k * k)


def gcpw_conformal(er, h, w, s):
    """Conductor-backed CPW, infinite grounds, zero thickness."""
    a = w / 2.0
    b = w / 2.0 + s
    k1 = a / b
    k3 = math.tanh(math.pi * a / (2 * h)) / math.tanh(math.pi * b / (2 * h))
    r1 = _kratio(k1)      # K(k1)/K(k1')
    r3 = _kratio(k3)      # K(k3)/K(k3')
    q = r3 / r1           # (K(k3)/K(k3')) * (K(k1')/K(k1))
    eps_eff = (1.0 + er * q) / (1.0 + q)
    z0 = 60.0 * math.pi / math.sqrt(eps_eff) / (r1 + r3)
    return {"z0": z0, "eps_eff": eps_eff}


# ---------------------------------------------------------------------------
# 2D finite-volume Laplace solver
# ---------------------------------------------------------------------------

def graded_line(x0, x1, points, hmin, hmax, ratio=1.15):
    """Grid lines on [x0, x1] that include `points`, with spacing hmin at each point
    growing geometrically by `ratio` up to hmax."""
    refine = [p for p in points if x0 <= p <= x1]
    pts = sorted(set([x0, x1] + refine))

    def size_at(x):
        # spacing grows linearly with distance from the nearest refinement point,
        # which is what a geometric progression with factor `ratio` produces
        d = min(abs(x - p) for p in refine) if refine else float("inf")
        return min(hmax, hmin + (ratio - 1.0) * d)

    lines = [x0]
    for a, b in zip(pts[:-1], pts[1:]):
        x = a
        seg = []
        while True:
            hs = size_at(x)
            hs = min(hs, size_at(x + hs))
            if x + hs >= b - 0.5 * hs:
                break
            x += hs
            seg.append(x)
        # rescale the segment so it ends exactly at b with smooth spacing
        if seg:
            edges = np.array([a] + seg + [b])
            gaps = np.diff(edges)
            gaps[-1] = max(gaps[-1], 1e-30)
            gaps *= (b - a) / gaps.sum()
            seg = list(a + np.cumsum(gaps)[:-1])
        lines.extend(seg)
        lines.append(b)
    return np.array(lines)


def laplace_capacitance(er, h, w, s, fence, ztop, hmin, hmax, ymax=None, outer_pec=False, t=0.0, tm=0.0, erm=1.0):
    """Per-unit-length capacitance (F/m) of the strip in the half-symmetric GCPW model.

    Geometry (half domain, y >= 0, magnetic wall at y = 0):
      z = 0         : ground (Dirichlet 0)
      z = h, y<=w/2 : strip (Dirichlet 1)
      z = h, y>=w/2+s : top ground (Dirichlet 0)
      y >= yf, z<=h : solid metal (via fence idealized as a wall), yf = w/2+s+fence
      dielectric er for 0<z<h, air above; Neumann at y = ymax and z = ztop.
    With t > 0 the strip and top grounds are copper bars from z = h to h + t.
    With tm > 0 a solder-mask slab of permittivity erm fills h < z < h + tm (the
    app's approximation: gaps filled, copper bars override).
    """
    a = w / 2.0
    b = w / 2.0 + s
    yf = b + fence
    if ymax is None:
        ymax = yf + 4 * h
    ys = graded_line(0.0, ymax, [a, b, yf], hmin, hmax)
    zs = graded_line(0.0, ztop, sorted(set([h] + ([h + t] if t > 0 else []) + ([h + tm] if tm > 0 else []))), hmin, hmax)
    ny, nz = len(ys), len(zs)
    dy = np.diff(ys)
    dz = np.diff(zs)
    # cell permittivity (cells between lines)
    zc = 0.5 * (zs[:-1] + zs[1:])
    yc = 0.5 * (ys[:-1] + ys[1:])
    epsc = np.where(zc[None, :] < h, er, np.where(zc[None, :] < h + tm, erm, 1.0)) * np.ones((ny - 1, 1))

    idx = lambda j, k: j * nz + k
    # Dirichlet values: NaN = unknown
    fixed = np.full((ny, nz), np.nan)
    tol = 1e-12
    fixed[:, 0] = 0.0
    kh = int(np.argmin(np.abs(zs - h)))
    kt = int(np.argmin(np.abs(zs - (h + t)))) if t > 0 else kh
    for j in range(ny):
        if ys[j] <= a + tol:
            fixed[j, kh:kt + 1] = 1.0
        elif ys[j] >= b - tol:
            fixed[j, kh:kt + 1] = 0.0
        if ys[j] >= yf - tol:
            fixed[j, :kh + 1] = 0.0
    if outer_pec:  # shielded box instead of open (Neumann) outer walls
        fixed[-1, :] = 0.0
        fixed[:, -1] = 0.0

    # Assemble the FV stiffness matrix: edge conductances
    rows, cols, vals = [], [], []
    edges = []  # (n1, n2, g)
    # y-directed edges between (j,k) and (j+1,k): length dy[j], dual width in z
    for j in range(ny - 1):
        for k in range(nz):
            g = 0.0
            if k > 0:
                g += epsc[j, k - 1] * dz[k - 1] / 2
            if k < nz - 1:
                g += epsc[j, k] * dz[k] / 2
            edges.append((idx(j, k), idx(j + 1, k), g / dy[j]))
    for j in range(ny):
        for k in range(nz - 1):
            g = 0.0
            if j > 0:
                g += epsc[j - 1, k] * dy[j - 1] / 2
            if j < ny - 1:
                g += epsc[j, k] * dy[j] / 2
            edges.append((idx(j, k), idx(j, k + 1), g / dz[k]))
    e = np.array(edges)
    n1 = e[:, 0].astype(int)
    n2 = e[:, 1].astype(int)
    g = e[:, 2]
    N = ny * nz
    A = sp.coo_matrix((np.concatenate([g, g, -g, -g]),
                       (np.concatenate([n1, n2, n1, n2]),
                        np.concatenate([n1, n2, n2, n1]))), shape=(N, N)).tocsr()
    fx = fixed.ravel()
    known = ~np.isnan(fx)
    unk = ~known
    phi = np.zeros(N)
    phi[known] = fx[known]
    Auu = A[unk][:, unk]
    rhs = -A[unk][:, known] @ phi[known]
    phi[unk] = spla.spsolve(Auu.tocsc(), rhs)
    energy2 = phi @ (A @ phi)          # = C/eps0 for V = 1 (half domain)
    return 2.0 * EPS0 * energy2, N   # x2 for the mirrored half


def gcpw_laplace(er, h, w, s, fence, ztop=None, levels=(1e-6, 0.5e-6, 0.25e-6), hmax=None, ymax=None, t=0.0, tm=0.0, erm=1.0):
    """Open-structure reference: the box is large enough (10 mm of air above, 6 mm
    beyond the fences) that Neumann and PEC outer walls agree to ~0.05%."""
    if ztop is None:
        ztop = h + 10e-3
    if ymax is None:
        ymax = w / 2 + s + fence + 6e-3
    if hmax is None:
        hmax = 100e-6
    res = []
    for hmin in levels:
        C, N = laplace_capacitance(er, h, w, s, fence, ztop, hmin, hmax, ymax=ymax, t=t, tm=tm, erm=erm)
        Ca, _ = laplace_capacitance(1.0, h, w, s, fence, ztop, hmin, hmax, ymax=ymax, t=t, tm=tm, erm=1.0)
        res.append((hmin, C, Ca, N))
    out = {"levels": []}
    for hmin, C, Ca, N in res:
        out["levels"].append({"hmin": hmin, "unknowns": N,
                              "z0": 1.0 / (C0 * math.sqrt(C * Ca)), "eps_eff": float(C / Ca)})
    # Richardson on C and Ca assuming error ~ hmin^p, estimate p from three levels
    if len(res) >= 3:
        def extrap(vals):
            v1, v2, v3 = vals[-3:]
            r = res[-3][0] / res[-2][0]
            denom = (v2 - v3)
            if abs(denom) < 1e-30 or (v1 - v2) / denom <= 0:
                return v3, float("nan")
            p = math.log((v1 - v2) / denom) / math.log(r)
            return v3 + (v3 - v2) / (r ** p - 1), p
        Cx, pc = extrap([r[1] for r in res])
        Cax, pa = extrap([r[2] for r in res])
        out["extrapolated"] = {"z0": 1.0 / (C0 * math.sqrt(Cx * Cax)), "eps_eff": float(Cx / Cax),
                               "order_C": pc, "order_C0": pa}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--er", type=float, default=3.66)
    ap.add_argument("--h", type=float, default=254e-6)
    ap.add_argument("--w", type=float, default=500e-6)
    ap.add_argument("--s", type=float, default=200e-6)
    ap.add_argument("--fence", type=float, default=300e-6,
                    help="distance from gap edge to via-fence wall")
    ap.add_argument("--t", type=float, default=0.0, help="copper thickness (m); 0 = zero-thickness sheets")
    ap.add_argument("--mask-t", dest="tm", type=float, default=0.0, help="solder-mask slab thickness (m)")
    ap.add_argument("--mask-er", dest="erm", type=float, default=3.5)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    cm = gcpw_conformal(a.er, a.h, a.w, a.s)
    lp = gcpw_laplace(a.er, a.h, a.w, a.s, a.fence, t=a.t, tm=a.tm, erm=a.erm)
    out = {"geometry": vars(a), "conformal": cm, "laplace": lp}
    if a.json:
        print(json.dumps(out, indent=2))
    else:
        print(f"conformal (no fence): Z0 = {cm['z0']:.3f} ohm, eps_eff = {cm['eps_eff']:.4f}")
        for lv in lp["levels"]:
            print(f"laplace hmin={lv['hmin']*1e6:4.1f}um N={lv['unknowns']:6d}: "
                  f"Z0 = {lv['z0']:.3f}, eps_eff = {lv['eps_eff']:.4f}")
        if "extrapolated" in lp:
            x = lp["extrapolated"]
            print(f"laplace extrapolated: Z0 = {x['z0']:.3f}, eps_eff = {x['eps_eff']:.4f} "
                  f"(orders {x['order_C']:.2f}, {x['order_C0']:.2f})")


if __name__ == "__main__":
    main()
