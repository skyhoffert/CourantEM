#!/usr/bin/env python3
"""Full-wave 2D eigenmode solver for a waveguide cross-section on a Yee grid.

Prototype of the mode-port solver (build step 3) and a frequency-dependent
reference for the FDTD line extraction. The cross-section is the (y, z) plane of
the 3D grid, propagation along +x with fields ~ exp(-j beta x). Fields live where
the 3D Yee scheme puts them, so on the FDTD's own mesh the discrete mode is the
one the 3D solver supports (up to x/time discretisation):

    E_u = Ey (j+1/2, k)   E_v = Ez (j, k+1/2)   E_w = Ex (j, k)      [nodes]
    H_u = Hy (j, k+1/2)   H_v = Hz (j+1/2, k)   H_w = Hx (j+1/2, k+1/2) [cells]

Eliminating the longitudinal components gives  beta e = P h,  beta h = Q e,
so  beta^2 e = P Q e  is solved for the transverse E field with shift-invert.

Impedance definitions reported: Z_vi = V/I, Z_pv = |V|^2/(2P), Z_pi = 2P/|I|^2.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6


class CrossSection:
    """Cross-section geometry on a tensor grid (y nodes, z nodes)."""

    def __init__(self, y, z):
        self.y = np.asarray(y, float)
        self.z = np.asarray(z, float)
        self.ny = len(y) - 1
        self.nz = len(z) - 1
        self.eps = np.ones((self.ny, self.nz))
        self.pec_cell = np.zeros((self.ny, self.nz), bool)
        self.pec_eu = np.zeros((self.ny, self.nz + 1), bool)   # Ey edges
        self.pec_ev = np.zeros((self.ny + 1, self.nz), bool)   # Ez edges
        self.pec_node = np.zeros((self.ny + 1, self.nz + 1), bool)
        # outer boundary: PEC walls
        self.pec_eu[:, 0] = self.pec_eu[:, -1] = True
        self.pec_ev[0, :] = self.pec_ev[-1, :] = True
        self.pec_node[0, :] = self.pec_node[-1, :] = True
        self.pec_node[:, 0] = self.pec_node[:, -1] = True

    def _cells(self, y0, y1, z0, z1):
        yc = 0.5 * (self.y[:-1] + self.y[1:])
        zc = 0.5 * (self.z[:-1] + self.z[1:])
        return ((yc >= y0) & (yc <= y1))[:, None] & ((zc >= z0) & (zc <= z1))[None, :]

    def material(self, y0, y1, z0, z1, er):
        self.eps[self._cells(y0, y1, z0, z1)] = er

    def pec_box(self, y0, y1, z0, z1):
        m = self._cells(y0, y1, z0, z1)
        self.pec_cell |= m
        for j, k in zip(*np.nonzero(m)):
            self.pec_eu[j, k] = self.pec_eu[j, k + 1] = True
            self.pec_ev[j, k] = self.pec_ev[j + 1, k] = True
            self.pec_node[j:j + 2, k:k + 2] = True

    def pec_sheet_z(self, zpos, y0, y1):
        """Zero-thickness horizontal sheet at node plane nearest zpos."""
        k = int(np.argmin(np.abs(self.z - zpos)))
        tol = 1e-12
        for j in range(self.ny):
            if self.y[j] >= y0 - tol and self.y[j + 1] <= y1 + tol:
                self.pec_eu[j, k] = True
        for j in range(self.ny + 1):
            if y0 - tol <= self.y[j] <= y1 + tol:
                self.pec_node[j, k] = True
        return k


def gcpw_section(y, z, er, h, w, s, fence):
    cs = CrossSection(y, z)
    a, b = w / 2, w / 2 + s
    yf = b + fence
    big = 1.0
    cs.material(-big, big, 0, h, er)
    cs.pec_box(-big, -yf, 0, h)
    cs.pec_box(yf, big, 0, h)
    cs.pec_sheet_z(h, -big, -b)
    cs.pec_sheet_z(h, b, big)
    kh = cs.pec_sheet_z(h, -a, a)
    return cs, kh


def _d(n):
    return n


class ModeSolver:
    def __init__(self, cs):
        self.cs = cs
        ny, nz = cs.ny, cs.nz
        y, z = cs.y, cs.z
        dy, dz = np.diff(y), np.diff(z)
        ddy = np.zeros(ny + 1); ddy[:-1] += dy / 2; ddy[1:] += dy / 2
        ddz = np.zeros(nz + 1); ddz[:-1] += dz / 2; ddz[1:] += dz / 2
        self.dy, self.dz, self.ddy, self.ddz = dy, dz, ddy, ddz
        # index maps
        self.iu = np.arange(ny * (nz + 1)).reshape(ny, nz + 1)          # E_u / H_v
        self.iv = np.arange((ny + 1) * nz).reshape(ny + 1, nz)          # E_v / H_u
        self.inode = np.arange((ny + 1) * (nz + 1)).reshape(ny + 1, nz + 1)
        self.icell = np.arange(ny * nz).reshape(ny, nz)
        Nu, Nv, Nn, Nc = ny * (nz + 1), (ny + 1) * nz, (ny + 1) * (nz + 1), ny * nz
        self.N = (Nu, Nv, Nn, Nc)

        def mat(rows, cols, vals, shape):
            return sp.csr_matrix((vals, (rows, cols)), shape=shape)

        # Du_e: nodes -> E_u ; Dv_e: nodes -> E_v
        J, K = np.meshgrid(np.arange(ny), np.arange(nz + 1), indexing="ij")
        r = self.iu[J, K].ravel()
        Du_e = mat(np.r_[r, r], np.r_[self.inode[J + 1, K].ravel(), self.inode[J, K].ravel()],
                   np.r_[1 / dy[J].ravel(), -1 / dy[J].ravel()], (Nu, Nn))
        J, K = np.meshgrid(np.arange(ny + 1), np.arange(nz), indexing="ij")
        r = self.iv[J, K].ravel()
        Dv_e = mat(np.r_[r, r], np.r_[self.inode[J, K + 1].ravel(), self.inode[J, K].ravel()],
                   np.r_[1 / dz[K].ravel(), -1 / dz[K].ravel()], (Nv, Nn))
        # Du_h: H_v (at E_u loc) -> nodes ; Dv_h: H_u (at E_v loc) -> nodes
        rows, cols, vals = [], [], []
        for j in range(ny + 1):
            for k in range(nz + 1):
                n = self.inode[j, k]
                if j < ny:
                    rows.append(n); cols.append(self.iu[j, k]); vals.append(1 / ddy[j])
                if j > 0:
                    rows.append(n); cols.append(self.iu[j - 1, k]); vals.append(-1 / ddy[j])
        Du_h = mat(rows, cols, vals, (Nn, Nu))
        rows, cols, vals = [], [], []
        for j in range(ny + 1):
            for k in range(nz + 1):
                n = self.inode[j, k]
                if k < nz:
                    rows.append(n); cols.append(self.iv[j, k]); vals.append(1 / ddz[k])
                if k > 0:
                    rows.append(n); cols.append(self.iv[j, k - 1]); vals.append(-1 / ddz[k])
        Dv_h = mat(rows, cols, vals, (Nn, Nv))
        # cells: Du_c E_v -> cells, Dv_c E_u -> cells
        J, K = np.meshgrid(np.arange(ny), np.arange(nz), indexing="ij")
        r = self.icell[J, K].ravel()
        Du_c = mat(np.r_[r, r], np.r_[self.iv[J + 1, K].ravel(), self.iv[J, K].ravel()],
                   np.r_[1 / dy[J].ravel(), -1 / dy[J].ravel()], (Nc, Nv))
        Dv_c = mat(np.r_[r, r], np.r_[self.iu[J, K + 1].ravel(), self.iu[J, K].ravel()],
                   np.r_[1 / dz[K].ravel(), -1 / dz[K].ravel()], (Nc, Nu))
        # cells -> H_u loc (j, k+1/2): diff in y (dual) ; cells -> H_v loc (j+1/2, k): diff in z
        rows, cols, vals = [], [], []
        for j in range(ny + 1):
            for k in range(nz):
                n = self.iv[j, k]
                if j < ny:
                    rows.append(n); cols.append(self.icell[j, k]); vals.append(1 / ddy[j])
                if j > 0:
                    rows.append(n); cols.append(self.icell[j - 1, k]); vals.append(-1 / ddy[j])
        Du_h2 = mat(rows, cols, vals, (Nv, Nc))
        rows, cols, vals = [], [], []
        for j in range(ny):
            for k in range(nz + 1):
                n = self.iu[j, k]
                if k < nz:
                    rows.append(n); cols.append(self.icell[j, k]); vals.append(1 / ddz[k])
                if k > 0:
                    rows.append(n); cols.append(self.icell[j, k - 1]); vals.append(-1 / ddz[k])
        Dv_h2 = mat(rows, cols, vals, (Nu, Nc))
        self.ops = (Du_e, Dv_e, Du_h, Dv_h, Du_c, Dv_c, Du_h2, Dv_h2)

        # edge / node permittivities (area weighted, same as the 3D solver)
        eps = cs.eps
        eu = np.zeros((ny, nz + 1)); wu = np.zeros((ny, nz + 1))
        eu[:, :-1] += eps * dz[None, :]; wu[:, :-1] += dz[None, :]
        eu[:, 1:] += eps * dz[None, :]; wu[:, 1:] += dz[None, :]
        ev = np.zeros((ny + 1, nz)); wv = np.zeros((ny + 1, nz))
        ev[:-1, :] += eps * dy[:, None]; wv[:-1, :] += dy[:, None]
        ev[1:, :] += eps * dy[:, None]; wv[1:, :] += dy[:, None]
        en = np.zeros((ny + 1, nz + 1)); wn = np.zeros((ny + 1, nz + 1))
        A = dy[:, None] * dz[None, :]
        for dj in (0, 1):
            for dk in (0, 1):
                en[dj:dj + ny, dk:dk + nz] += eps * A
                wn[dj:dj + ny, dk:dk + nz] += A
        self.eps_u = (eu / wu).ravel() * EPS0
        self.eps_v = (ev / wv).ravel() * EPS0
        inv_n = (wn / en).ravel() / EPS0
        inv_n[cs.pec_node.ravel()] = 0.0
        self.inv_eps_w = inv_n
        self.free = np.r_[~cs.pec_eu.ravel(), ~cs.pec_ev.ravel()]

    def operators(self, f):
        w = 2 * math.pi * f
        Du_e, Dv_e, Du_h, Dv_h, Du_c, Dv_c, Du_h2, Dv_h2 = self.ops
        Nu, Nv, Nn, Nc = self.N
        Ie = sp.diags(self.inv_eps_w)
        Cw = sp.hstack([-Dv_h, Du_h])                  # on [H_u; H_v] -> nodes
        G = sp.vstack([Du_e, Dv_e])                     # nodes -> [E_u; E_v]
        # P maps [H_u; H_v] -> [E_u; E_v]
        rot = sp.bmat([[None, sp.identity(Nu)], [-sp.identity(Nv), None]])  # E_u <- H_v, E_v <- -H_u
        P = w * MU0 * rot + (1 / w) * G @ Ie @ Cw
        # Q maps [E_u; E_v] -> [H_u; H_v]
        Cwe = sp.hstack([-Dv_c, Du_c])                  # on [E_u; E_v] -> cells
        Gh = sp.vstack([Du_h2, Dv_h2])                  # cells -> [H_u; H_v]
        rot2 = sp.bmat([[None, -w * sp.diags(self.eps_v)], [w * sp.diags(self.eps_u), None]])
        Q = rot2 - (1 / (w * MU0)) * Gh @ Cwe
        return P.tocsr(), Q.tocsr()

    def solve(self, f, neff_guess, nmodes=4):
        P, Q = self.operators(f)
        fr = np.nonzero(self.free)[0]
        A = (P[fr, :] @ Q[:, fr]).tocsc()
        k0 = 2 * math.pi * f / C0
        sigma = (k0 * neff_guess) ** 2
        # explicit LU with full partial pivoting: the default threshold pivoting
        # loses the guided mode at low frequency where |A| >> beta^2
        lu = spla.splu((A - sigma * sp.identity(A.shape[0], format="csc")).tocsc(),
                       diag_pivot_thresh=1.0, permc_spec="COLAMD")
        op = spla.LinearOperator(A.shape, matvec=lu.solve, dtype=float)
        vals, vecs = spla.eigs(A, k=nmodes, sigma=sigma, OPinv=op)
        order = np.argsort(np.abs(vals - sigma))
        modes = []
        for idx in order:
            b2 = vals[idx]
            beta = np.sqrt(b2)
            if beta.real < 0:
                beta = -beta
            e = np.zeros(self.free.size, complex)
            e[fr] = vecs[:, idx]
            h = (Q @ e) / beta
            modes.append({"beta": beta, "neff": beta / k0, "e": e, "h": h})
        return modes

    def split(self, v):
        Nu, Nv, _, _ = self.N
        ny, nz = self.cs.ny, self.cs.nz
        return v[:Nu].reshape(ny, nz + 1), v[Nu:].reshape(ny + 1, nz)

    def power(self, mode):
        eu, ev = self.split(mode["e"])
        Nu, Nv, _, _ = self.N
        ny, nz = self.cs.ny, self.cs.nz
        Hu = mode["h"][:Nv].reshape(ny + 1, nz)
        Hv = mode["h"][Nv:].reshape(ny, nz + 1)
        Au = self.dy[:, None] * self.ddz[None, :]
        Av = self.ddy[:, None] * self.dz[None, :]
        return 0.5 * np.real(np.sum(eu * np.conj(Hv) * Au) - np.sum(ev * np.conj(Hu) * Av)), Hu, Hv

    def vi(self, mode, j0, k0, k1, loop):
        """V = -sum Ez dz along node column j0 from k0 to k1 (ground to strip), and
        I = circulation of H around the dual rectangle loop=(ja, jb, ka, kb) (current in +x)."""
        eu, ev = self.split(mode["e"])
        P, Hu, Hv = self.power(mode)
        V = -np.sum(ev[j0, k0:k1] * self.dz[k0:k1])
        ja, jb, ka, kb = loop
        I = 0
        for j in range(ja + 1, jb + 1):
            I += Hu[j, ka] * self.ddy[j] - Hu[j, kb] * self.ddy[j]
        for k in range(ka + 1, kb + 1):
            I += Hv[jb, k] * self.ddz[k] - Hv[ja, k] * self.ddz[k]
        return V, I, P


def gcpw_modes(mesh, freqs, crop_pml=True):
    """Solve the GCPW mode on the FDTD mesh (mesh.json from gcpw_line)."""
    g = mesh["geometry"]
    y = np.array(mesh["y"]); z = np.array(mesh["z"])
    npml = mesh["npml"]
    if crop_pml:
        y = y[npml:len(y) - npml]
        z = z[:len(z) - npml]
    cs, kh = gcpw_section(y, z, g["er"], g["h"], g["w"], g["s"], g["fence"])
    ms = ModeSolver(cs)
    j0 = int(np.argmin(np.abs(y)))
    a, b = g["w"] / 2, g["w"] / 2 + g["s"]
    jl = int(np.argmin(np.abs(y + (a + b) / 2))); jr = int(np.argmin(np.abs(y - (a + b) / 2)))
    loop = (jl, jr - 1, kh - 1, kh)
    out = []
    guess = None
    from gcpw_ref import gcpw_conformal
    qs = gcpw_conformal(g["er"], g["h"], g["w"], g["s"])
    for f in freqs:
        ng = math.sqrt(qs["eps_eff"]) if guess is None else guess
        modes = ms.solve(f, ng, nmodes=3)
        # pick the mode with the largest |V| relative to sqrt(P): the strip mode
        best = None
        for m in modes:
            V, I, P = ms.vi(m, j0, 0, kh, loop)
            if P <= 0:
                continue
            score = abs(V) ** 2 / P
            if best is None or score > best[0]:
                best = (score, m, V, I, P)
        _, m, V, I, P = best
        guess = m["neff"].real
        out.append({"f": f, "neff": complex(m["neff"]), "eps_eff": float((m["neff"] ** 2).real),
                    "z_vi": complex(V / I), "z_pv": float(abs(V) ** 2 / (2 * P)), "z_pi": float(2 * P / abs(I) ** 2)})
    return out, ms


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", help="mesh.json written by gcpw_line")
    ap.add_argument("--freqs", default="0.5,5,10,20,30,40", help="GHz list")
    ap.add_argument("--json")
    a = ap.parse_args()
    with open(a.mesh) as fh:
        mesh = json.load(fh)
    freqs = [float(v) * 1e9 for v in a.freqs.split(",")]
    res, _ = gcpw_modes(mesh, freqs)
    for r in res:
        print(f"{r['f']/1e9:6.2f} GHz  eps_eff {r['eps_eff']:.5f}  Zvi {r['z_vi'].real:8.3f}{r['z_vi'].imag:+.3f}j  "
              f"Zpv {r['z_pv']:8.3f}  Zpi {r['z_pi']:8.3f}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump([{k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in r.items()} for r in res], fh, indent=1)


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(__file__))
    main()
