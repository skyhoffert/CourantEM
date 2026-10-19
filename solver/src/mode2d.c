/*
 * Courant EM - 2D full-wave eigenmode solver on a Yee cross-section.
 *
 * Same formulation as tools/modesolver2d.py (the verified Python prototype):
 * propagation along +x, fields ~ exp(-j beta x), transverse plane (y, z) with
 *   E_u = Ey (j+1/2, k)   E_v = Ez (j, k+1/2)   E_w = Ex (j, k)
 *   H_u = Hy (j, k+1/2)   H_v = Hz (j+1/2, k)   H_w = Hx (j+1/2, k+1/2)
 * and  beta e = P h,  beta h = Q e  =>  beta^2 e = P Q e.
 *
 * Complex arithmetic so the cross-section can carry the same stretched-coordinate
 * PML as the 3D CPML (s = kappa + sigma / (alpha + j w eps0)) and lossy media
 * (eps - j sigma / w). Every derivative uses 1/(s * spacing) evaluated where the
 * derivative lands: cells for differences of E, nodes for differences of H.
 *
 * The operator A = P Q is applied matrix-free, assembled into band storage by
 * probing with coloured unit vectors, factorised once with a banded LU (partial
 * pivoting), and the mode nearest sigma = (k0 n_guess)^2 is found by inverse
 * iteration. Unknowns are ordered (j outer, k inner, [E_u, E_v] interleaved) so
 * the bandwidth is about 2 (nz + 2).
 */
#include "courant.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double complex cplx;

struct cem_mode2d {
    int ny, nz;
    double *y, *z, *dy, *dz, *ddy, *ddz;
    double *er_u, *er_v, *er_w;      /* relative permittivity per unknown */
    double *sg_u, *sg_v, *sg_w;      /* conductivity (S/m) per unknown */
    unsigned char *pec_u, *pec_v, *pec_w;
    /* Debye dispersion */
    int npole;
    double pole_f[16];
    double *dp_u, *dp_v, *dp_w;      /* npole blocks */
    /* PML profiles: sigma, kappa, alpha at y cells / y nodes / z cells / z nodes */
    double *ps[4], *pk[4], *pa[4];
    /* last solution */
    double freq;
    cplx beta;
    cplx *eu, *ev, *hu, *hv;
};

#define IU(j, k) ((j) * (nz + 1) + (k))   /* E_u / H_v: j < ny, k <= nz */
#define IV(j, k) ((j) * (nz) + (k))       /* E_v / H_u: j <= ny, k < nz */
#define IN(j, k) ((j) * (nz + 1) + (k))   /* nodes */
#define IC(j, k) ((j) * (nz) + (k))       /* cells */

cem_mode2d *cem_mode2d_create(int ny, int nz, const double *y, const double *z) {
    cem_mode2d *m = calloc(1, sizeof *m);
    m->ny = ny; m->nz = nz;
    m->y = malloc((ny + 1) * sizeof(double)); memcpy(m->y, y, (ny + 1) * sizeof(double));
    m->z = malloc((nz + 1) * sizeof(double)); memcpy(m->z, z, (nz + 1) * sizeof(double));
    m->dy = malloc(ny * sizeof(double)); m->dz = malloc(nz * sizeof(double));
    m->ddy = calloc(ny + 1, sizeof(double)); m->ddz = calloc(nz + 1, sizeof(double));
    for (int j = 0; j < ny; j++) { m->dy[j] = y[j + 1] - y[j]; m->ddy[j] += m->dy[j] / 2; m->ddy[j + 1] += m->dy[j] / 2; }
    for (int k = 0; k < nz; k++) { m->dz[k] = z[k + 1] - z[k]; m->ddz[k] += m->dz[k] / 2; m->ddz[k + 1] += m->dz[k] / 2; }
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    m->er_u = malloc(nu * sizeof(double)); m->er_v = malloc(nv * sizeof(double)); m->er_w = malloc(nn * sizeof(double));
    m->sg_u = calloc(nu, sizeof(double)); m->sg_v = calloc(nv, sizeof(double)); m->sg_w = calloc(nn, sizeof(double));
    m->pec_u = calloc(nu, 1); m->pec_v = calloc(nv, 1); m->pec_w = calloc(nn, 1);
    for (int q = 0; q < nu; q++) m->er_u[q] = 1;
    for (int q = 0; q < nv; q++) m->er_v[q] = 1;
    for (int q = 0; q < nn; q++) m->er_w[q] = 1;
    int len[4] = {ny, ny + 1, nz, nz + 1};
    for (int a = 0; a < 4; a++) {
        m->ps[a] = calloc(len[a], sizeof(double));
        m->pk[a] = malloc(len[a] * sizeof(double));
        m->pa[a] = calloc(len[a], sizeof(double));
        for (int i = 0; i < len[a]; i++) m->pk[a][i] = 1;
    }
    /* outer boundary is PEC */
    for (int j = 0; j < ny; j++) { m->pec_u[IU(j, 0)] = m->pec_u[IU(j, nz)] = 1; }
    for (int k = 0; k < nz; k++) { m->pec_v[IV(0, k)] = m->pec_v[IV(ny, k)] = 1; }
    for (int j = 0; j <= ny; j++) { m->pec_w[IN(j, 0)] = m->pec_w[IN(j, nz)] = 1; }
    for (int k = 0; k <= nz; k++) { m->pec_w[IN(0, k)] = m->pec_w[IN(ny, k)] = 1; }
    m->eu = calloc(nu, sizeof(cplx)); m->ev = calloc(nv, sizeof(cplx));
    m->hu = calloc(nv, sizeof(cplx)); m->hv = calloc(nu, sizeof(cplx));
    return m;
}

void cem_mode2d_destroy(cem_mode2d *m) {
    if (!m) return;
    free(m->y); free(m->z); free(m->dy); free(m->dz); free(m->ddy); free(m->ddz);
    free(m->er_u); free(m->er_v); free(m->er_w); free(m->sg_u); free(m->sg_v); free(m->sg_w);
    free(m->pec_u); free(m->pec_v); free(m->pec_w);
    for (int a = 0; a < 4; a++) { free(m->ps[a]); free(m->pk[a]); free(m->pa[a]); }
    free(m->eu); free(m->ev); free(m->hu); free(m->hv);
    free(m->dp_u); free(m->dp_v); free(m->dp_w);
    free(m);
}

void cem_mode2d_set_debye(cem_mode2d *m, int npole, const double *f_p, const double *deps_u, const double *deps_v,
                          const double *deps_w) {
    int ny = m->ny, nz = m->nz;
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    if (npole > 16) npole = 16;
    m->npole = npole;
    memcpy(m->pole_f, f_p, npole * sizeof(double));
    free(m->dp_u); free(m->dp_v); free(m->dp_w);
    m->dp_u = malloc((size_t)npole * nu * sizeof(double)); memcpy(m->dp_u, deps_u, (size_t)npole * nu * sizeof(double));
    m->dp_v = malloc((size_t)npole * nv * sizeof(double)); memcpy(m->dp_v, deps_v, (size_t)npole * nv * sizeof(double));
    m->dp_w = malloc((size_t)npole * nn * sizeof(double)); memcpy(m->dp_w, deps_w, (size_t)npole * nn * sizeof(double));
}

void cem_mode2d_set_media(cem_mode2d *m, const double *epsr_u, const double *epsr_v, const double *epsr_w,
                          const unsigned char *pec_u, const unsigned char *pec_v, const unsigned char *pec_w) {
    int ny = m->ny, nz = m->nz;
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    for (int q = 0; q < nu; q++) { m->er_u[q] = epsr_u[q]; m->pec_u[q] |= pec_u[q]; }
    for (int q = 0; q < nv; q++) { m->er_v[q] = epsr_v[q]; m->pec_v[q] |= pec_v[q]; }
    for (int q = 0; q < nn; q++) { m->er_w[q] = epsr_w[q]; m->pec_w[q] |= pec_w[q]; }
}

void cem_mode2d_set_loss(cem_mode2d *m, const double *sig_u, const double *sig_v, const double *sig_w) {
    int ny = m->ny, nz = m->nz;
    memcpy(m->sg_u, sig_u, ny * (nz + 1) * sizeof(double));
    memcpy(m->sg_v, sig_v, (ny + 1) * nz * sizeof(double));
    memcpy(m->sg_w, sig_w, (ny + 1) * (nz + 1) * sizeof(double));
}

void cem_mode2d_set_pml(cem_mode2d *m, int which, const double *sigma, const double *kappa, const double *alpha) {
    int len[4] = {m->ny, m->ny + 1, m->nz, m->nz + 1};
    if (which < 0 || which > 3) return;
    memcpy(m->ps[which], sigma, len[which] * sizeof(double));
    memcpy(m->pk[which], kappa, len[which] * sizeof(double));
    memcpy(m->pa[which], alpha, len[which] * sizeof(double));
}

/* per-solve coefficient set */
typedef struct {
    cplx *iyc, *iyn, *izc, *izn;   /* 1/(s*d) at y cells, y nodes, z cells, z nodes */
    cplx *eu, *ev, *iw;            /* complex permittivity (absolute) and 1/eps at E_w nodes (0 on PEC) */
} coef_t;

static void make_coef(const cem_mode2d *m, double w, coef_t *c) {
    int ny = m->ny, nz = m->nz;
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    c->iyc = malloc(ny * sizeof(cplx)); c->iyn = malloc((ny + 1) * sizeof(cplx));
    c->izc = malloc(nz * sizeof(cplx)); c->izn = malloc((nz + 1) * sizeof(cplx));
    cplx jw = I * w * CEM_EPS0;
#define STRETCH(a, i) (m->pk[a][i] + m->ps[a][i] / (m->pa[a][i] + jw))
    for (int j = 0; j < ny; j++) c->iyc[j] = 1.0 / (STRETCH(0, j) * m->dy[j]);
    for (int j = 0; j <= ny; j++) c->iyn[j] = 1.0 / (STRETCH(1, j) * m->ddy[j]);
    for (int k = 0; k < nz; k++) c->izc[k] = 1.0 / (STRETCH(2, k) * m->dz[k]);
    for (int k = 0; k <= nz; k++) c->izn[k] = 1.0 / (STRETCH(3, k) * m->ddz[k]);
#undef STRETCH
    c->eu = malloc(nu * sizeof(cplx)); c->ev = malloc(nv * sizeof(cplx)); c->iw = malloc(nn * sizeof(cplx));
    cplx dfac[16];
    for (int p = 0; p < m->npole; p++) dfac[p] = 1.0 / (1.0 + I * w / (2 * M_PI * m->pole_f[p]));
    for (int q = 0; q < nu; q++) {
        cplx e = m->er_u[q];
        for (int p = 0; p < m->npole; p++) e += m->dp_u[(size_t)p * nu + q] * dfac[p];
        c->eu[q] = e * CEM_EPS0 - I * m->sg_u[q] / w;
    }
    for (int q = 0; q < nv; q++) {
        cplx e = m->er_v[q];
        for (int p = 0; p < m->npole; p++) e += m->dp_v[(size_t)p * nv + q] * dfac[p];
        c->ev[q] = e * CEM_EPS0 - I * m->sg_v[q] / w;
    }
    for (int q = 0; q < nn; q++) {
        cplx e = m->er_w[q];
        for (int p = 0; p < m->npole; p++) e += m->dp_w[(size_t)p * nn + q] * dfac[p];
        c->iw[q] = m->pec_w[q] ? 0 : 1.0 / (e * CEM_EPS0 - I * m->sg_w[q] / w);
    }
}

static void free_coef(coef_t *c) {
    free(c->iyc); free(c->iyn); free(c->izc); free(c->izn); free(c->eu); free(c->ev); free(c->iw);
}

/* h = Q e ; hw: cells work array */
static void apply_Q(const cem_mode2d *m, const coef_t *c, double w, const cplx *eu, const cplx *ev, cplx *hu,
                    cplx *hv, cplx *hw) {
    int ny = m->ny, nz = m->nz;
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            hw[IC(j, k)] = (ev[IV(j + 1, k)] - ev[IV(j, k)]) * c->iyc[j] - (eu[IU(j, k + 1)] - eu[IU(j, k)]) * c->izc[k];
    double cc = 1.0 / (w * CEM_MU0);
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k < nz; k++) {
            cplx r = (j < ny ? hw[IC(j, k)] : 0.0) - (j > 0 ? hw[IC(j - 1, k)] : 0.0);
            hu[IV(j, k)] = -w * c->ev[IV(j, k)] * ev[IV(j, k)] - cc * r * c->iyn[j];
        }
    for (int j = 0; j < ny; j++)
        for (int k = 0; k <= nz; k++) {
            cplx r = (k < nz ? hw[IC(j, k)] : 0.0) - (k > 0 ? hw[IC(j, k - 1)] : 0.0);
            hv[IU(j, k)] = w * c->eu[IU(j, k)] * eu[IU(j, k)] - cc * r * c->izn[k];
        }
}

/* e = P h ; ew: nodes work array */
static void apply_P(const cem_mode2d *m, const coef_t *c, double w, const cplx *hu, const cplx *hv, cplx *eu,
                    cplx *ev, cplx *ew) {
    int ny = m->ny, nz = m->nz;
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k <= nz; k++) {
            cplx du = ((j < ny ? hv[IU(j, k)] : 0.0) - (j > 0 ? hv[IU(j - 1, k)] : 0.0)) * c->iyn[j];
            cplx dv = ((k < nz ? hu[IV(j, k)] : 0.0) - (k > 0 ? hu[IV(j, k - 1)] : 0.0)) * c->izn[k];
            ew[IN(j, k)] = (du - dv) * c->iw[IN(j, k)];
        }
    double cc = 1.0 / w;
    for (int j = 0; j < ny; j++)
        for (int k = 0; k <= nz; k++)
            eu[IU(j, k)] = w * CEM_MU0 * hv[IU(j, k)] + cc * (ew[IN(j + 1, k)] - ew[IN(j, k)]) * c->iyc[j];
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k < nz; k++)
            ev[IV(j, k)] = -w * CEM_MU0 * hu[IV(j, k)] + cc * (ew[IN(j, k + 1)] - ew[IN(j, k)]) * c->izc[k];
}

/* ---- ordering ------------------------------------------------------------ */
typedef struct { int n; int *kind; int *idx; } order_t;

static void make_order(const cem_mode2d *m, order_t *o) {
    int ny = m->ny, nz = m->nz;
    int nmax = ny * (nz + 1) + (ny + 1) * nz;
    o->kind = malloc(nmax * sizeof(int));
    o->idx = malloc(nmax * sizeof(int));
    int n = 0;
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k <= nz; k++) {
            if (j < ny) { o->kind[n] = 0; o->idx[n] = IU(j, k); n++; }
            if (k < nz) { o->kind[n] = 1; o->idx[n] = IV(j, k); n++; }
        }
    o->n = n;
}

static int is_pec(const cem_mode2d *m, const order_t *o, int g) {
    return o->kind[g] == 0 ? m->pec_u[o->idx[g]] : m->pec_v[o->idx[g]];
}

/* ---- complex band LU (LAPACK zgbtrf-like layout) -------------------------- */
typedef struct { int n, kl, ku, ldab; cplx *ab; int *ipiv; } band_t;

static inline cplx *bel(band_t *b, int i, int j) { return &b->ab[(size_t)j * b->ldab + (b->kl + b->ku + i - j)]; }

static int band_lu(band_t *b) {
    int n = b->n, kl = b->kl, kv = b->ku + b->kl;
    for (int j = 0; j < n; j++) {
        int last = j + kl < n - 1 ? j + kl : n - 1;
        int p = j;
        double best = cabs(*bel(b, j, j));
        for (int i = j + 1; i <= last; i++) {
            double v = cabs(*bel(b, i, j));
            if (v > best) { best = v; p = i; }
        }
        b->ipiv[j] = p;
        if (best == 0) return -1;
        int jmax = j + kv < n - 1 ? j + kv : n - 1;
        if (p != j)
            for (int c = j; c <= jmax; c++) {
                cplx *a = bel(b, j, c), *q = bel(b, p, c), t = *a;
                *a = *q; *q = t;
            }
        cplx piv = *bel(b, j, j);
        for (int i = j + 1; i <= last; i++) {
            cplx *lij = bel(b, i, j);
            *lij /= piv;
            cplx l = *lij;
            if (l == 0) continue;
            for (int c = j + 1; c <= jmax; c++) *bel(b, i, c) -= l * *bel(b, j, c);
        }
    }
    return 0;
}

static void band_solve(band_t *b, cplx *x) {
    int n = b->n, kl = b->kl, kv = b->ku + b->kl;
    for (int j = 0; j < n; j++) {
        int p = b->ipiv[j];
        if (p != j) { cplx t = x[j]; x[j] = x[p]; x[p] = t; }
        int last = j + kl < n - 1 ? j + kl : n - 1;
        for (int i = j + 1; i <= last; i++) x[i] -= *bel(b, i, j) * x[j];
    }
    for (int j = n - 1; j >= 0; j--) {
        x[j] /= *bel(b, j, j);
        int first = j - kv > 0 ? j - kv : 0;
        for (int i = first; i < j; i++) x[i] -= *bel(b, i, j) * x[j];
    }
}

typedef struct { cplx *eu, *ev, *hu, *hv, *wk; } work_t;

static void apply_A(const cem_mode2d *m, const coef_t *c, const order_t *o, double w, const cplx *x, cplx *y,
                    work_t *W) {
    int ny = m->ny, nz = m->nz;
    memset(W->eu, 0, ny * (nz + 1) * sizeof(cplx));
    memset(W->ev, 0, (ny + 1) * nz * sizeof(cplx));
    for (int g = 0; g < o->n; g++) {
        if (is_pec(m, o, g)) continue;
        if (o->kind[g] == 0) W->eu[o->idx[g]] = x[g]; else W->ev[o->idx[g]] = x[g];
    }
    apply_Q(m, c, w, W->eu, W->ev, W->hu, W->hv, W->wk);
    apply_P(m, c, w, W->hu, W->hv, W->eu, W->ev, W->wk);
    for (int g = 0; g < o->n; g++)
        y[g] = is_pec(m, o, g) ? 0 : (o->kind[g] == 0 ? W->eu[o->idx[g]] : W->ev[o->idx[g]]);
}

int cem_mode2d_solve(cem_mode2d *m, double f, double neff_guess, const double *start, int maxit, double tol) {
    int ny = m->ny, nz = m->nz;
    double w = 2 * M_PI * f, k0 = w / CEM_C0;
    cplx sigma = k0 * k0 * neff_guess * neff_guess;
    coef_t c;
    make_coef(m, w, &c);
    order_t o;
    make_order(m, &o);
    int n = o.n, bw = 2 * (nz + 2) + 2;
    band_t B = {n, bw, bw, 3 * bw + 1, NULL, NULL};
    B.ab = calloc((size_t)B.ldab * n, sizeof(cplx));
    B.ipiv = malloc(n * sizeof(int));
    cplx *x = calloc(n, sizeof(cplx)), *y = calloc(n, sizeof(cplx));
    int nn = (ny + 1) * (nz + 1);
    work_t W = {calloc(ny * (nz + 1), sizeof(cplx)), calloc((ny + 1) * nz, sizeof(cplx)),
                calloc((ny + 1) * nz, sizeof(cplx)), calloc(ny * (nz + 1), sizeof(cplx)), calloc(nn, sizeof(cplx))};
    int rc = -1;
    if (!B.ab || !B.ipiv || !x || !y) goto done;
    int Wc = 2 * bw + 1;
    for (int col = 0; col < Wc; col++) {
        memset(x, 0, n * sizeof(cplx));
        for (int g = col; g < n; g += Wc) x[g] = 1.0;
        apply_A(m, &c, &o, w, x, y, &W);
        for (int g = col; g < n; g += Wc) {
            int lo = g - bw < 0 ? 0 : g - bw, hi = g + bw > n - 1 ? n - 1 : g + bw;
            for (int r = lo; r <= hi; r++) *bel(&B, r, g) = y[r];
        }
    }
    double big = -1e6 * (cabs(sigma) + 1.0);
    for (int g = 0; g < n; g++) {
        if (is_pec(m, &o, g)) *bel(&B, g, g) = big;
        *bel(&B, g, g) -= sigma;
    }
    if (band_lu(&B)) goto done;
    if (start) {
        const cplx *st = (const cplx *)start;
        for (int g = 0; g < n; g++) x[g] = o.kind[g] == 0 ? st[o.idx[g]] : st[ny * (nz + 1) + o.idx[g]];
    } else {
        unsigned s = 12345u;
        for (int g = 0; g < n; g++) { s = s * 1103515245u + 12345u; x[g] = ((s >> 8) & 0xffff) / 65536.0 - 0.5; }
    }
    cplx lam = sigma, lam_old = 0;
    cplx *xo = calloc(n, sizeof(cplx));
    rc = 1;
    for (int it = 0; it < maxit; it++) {
        double nx2 = 0;
        for (int g = 0; g < n; g++) nx2 += creal(x[g] * conj(x[g]));
        double inv = 1 / sqrt(nx2);
        for (int g = 0; g < n; g++) x[g] *= inv;
        /* convergence of the vector itself (phase-aligned), not only of lambda:
         * the eigenvector converges more slowly and sets the source profile */
        cplx ph = 0;
        for (int g = 0; g < n; g++) ph += conj(xo[g]) * x[g];
        double dv = 0;
        if (cabs(ph) > 0) {
            ph = conj(ph) / cabs(ph);
            for (int g = 0; g < n; g++) { cplx d = x[g] * ph - xo[g]; dv += creal(d * conj(d)); }
        } else dv = 1;
        memcpy(xo, x, n * sizeof(cplx));
        memcpy(y, x, n * sizeof(cplx));
        band_solve(&B, y);
        cplx mu = 0;
        for (int g = 0; g < n; g++) mu += conj(x[g]) * y[g];
        lam = sigma + 1 / mu;
        memcpy(x, y, n * sizeof(cplx));
        if (it > 2 && cabs(lam - lam_old) < tol * cabs(lam) && sqrt(dv) < 1e-10) { rc = 0; break; }
        lam_old = lam;
    }
    free(xo);
    {
        /* normalise; fix the global phase so the largest component is real */
        double nx2 = 0, amax = 0;
        cplx ph = 1;
        for (int g = 0; g < n; g++) {
            nx2 += creal(x[g] * conj(x[g]));
            if (cabs(x[g]) > amax) { amax = cabs(x[g]); ph = conj(x[g]) / amax; }
        }
        double inv = 1 / sqrt(nx2);
        memset(m->eu, 0, ny * (nz + 1) * sizeof(cplx));
        memset(m->ev, 0, (ny + 1) * nz * sizeof(cplx));
        for (int g = 0; g < n; g++) {
            if (is_pec(m, &o, g)) continue;
            if (o.kind[g] == 0) m->eu[o.idx[g]] = x[g] * inv * ph;
            else m->ev[o.idx[g]] = x[g] * inv * ph;
        }
        cplx bt = csqrt(lam);
        if (creal(bt) < 0) bt = -bt;
        m->beta = bt;
        m->freq = f;
        apply_Q(m, &c, w, m->eu, m->ev, m->hu, m->hv, W.wk);
        for (int q = 0; q < (ny + 1) * nz; q++) m->hu[q] /= bt;
        for (int q = 0; q < ny * (nz + 1); q++) m->hv[q] /= bt;
    }
done:
    free(B.ab); free(B.ipiv); free(x); free(y);
    free(W.eu); free(W.ev); free(W.hu); free(W.hv); free(W.wk);
    free(o.kind); free(o.idx);
    free_coef(&c);
    return rc;
}

double cem_mode2d_beta(const cem_mode2d *m) { return creal(m->beta); }
double cem_mode2d_alpha(const cem_mode2d *m) { return -cimag(m->beta); }
double cem_mode_freq_hint(const cem_mode2d *m) { return m->freq; }
const double *cem_mode2d_eu(const cem_mode2d *m) { return (const double *)m->eu; }
const double *cem_mode2d_ev(const cem_mode2d *m) { return (const double *)m->ev; }
const double *cem_mode2d_hu(const cem_mode2d *m) { return (const double *)m->hu; }
const double *cem_mode2d_hv(const cem_mode2d *m) { return (const double *)m->hv; }

/* Complex power 1/2 sum (Eu Hv* - Ev Hu*) dA, real part returned. */
double cem_mode2d_power(const cem_mode2d *m) {
    int ny = m->ny, nz = m->nz;
    cplx P = 0;
    for (int j = 0; j < ny; j++)
        for (int k = 0; k <= nz; k++) P += m->eu[IU(j, k)] * conj(m->hv[IU(j, k)]) * m->dy[j] * m->ddz[k];
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k < nz; k++) P -= m->ev[IV(j, k)] * conj(m->hu[IV(j, k)]) * m->ddy[j] * m->dz[k];
    return 0.5 * creal(P);
}

/* V = -sum Ez dz on node column j0 from k0 to k1 ; I = circulation around the dual
 * rectangle (ja+1/2..jb+1/2) x (ka+1/2..kb+1/2), current in +x. Complex outputs. */
void cem_mode2d_vi(const cem_mode2d *m, int j0, int k0, int k1, int ja, int jb, int ka, int kb, double *V,
                   double *Iout) {
    int nz = m->nz;
    cplx v = 0, c = 0;
    for (int k = k0; k < k1; k++) v -= m->ev[IV(j0, k)] * m->dz[k];
    for (int j = ja + 1; j <= jb; j++) c += (m->hu[IV(j, ka)] - m->hu[IV(j, kb)]) * m->ddy[j];
    for (int k = ka + 1; k <= kb; k++) c += (m->hv[IU(jb, k)] - m->hv[IU(ja, k)]) * m->ddz[k];
    V[0] = creal(v); V[1] = cimag(v);
    Iout[0] = creal(c); Iout[1] = cimag(c);
}

/* Conductor attenuation (Np/m) by perturbation: P_loss = Rs/2 * sum |n x H|^2 dl over
 * every metal surface facing a field region, alpha_c = P_loss / (2 P). Uses the last
 * solved (lossless-conductor) mode. Longitudinal currents come from the transverse H
 * next to the surface, transverse currents from H_w. */
double cem_mode2d_conductor_alpha(const cem_mode2d *m, double Rs) {
    int ny = m->ny, nz = m->nz;
#define SOLID(j, k) ((j) >= 0 && (j) < ny && (k) >= 0 && (k) < nz && m->pec_u[IU(j, k)] && m->pec_u[IU(j, (k) + 1)] && \
                     m->pec_v[IV(j, k)] && m->pec_v[IV((j) + 1, k)])
    /* H_w at cells: H_w = -curl_w(E) / (j w mu) -> |H_w| = |curl_w E| / (w mu) */
    double w = 2 * M_PI * m->freq;
    double loss = 0;
    /* horizontal surfaces: E_u segments (j, k) */
    for (int j = 0; j < ny; j++)
        for (int k = 0; k <= nz; k++) {
            if (!m->pec_u[IU(j, k)]) continue;
            /* metal inside the PML is not part of the device (and its stretched
             * fields are not physical); the 3D SIBC skips it too */
            if (m->ps[0][j] > 0 || m->ps[3][k] > 0) continue;
            for (int side = 0; side < 2; side++) {
                int kc = side == 0 ? k : k - 1;          /* cell on the field side */
                if (kc < 0 || kc >= nz || SOLID(j, kc)) continue;
                /* transverse current |H_w| on the segment */
                cplx cw = (m->ev[IV(j + 1, kc)] - m->ev[IV(j, kc)]) / m->dy[j] - (m->eu[IU(j, kc + 1)] - m->eu[IU(j, kc)]) / m->dz[kc];
                double hw = cabs(cw) / (w * CEM_MU0);
                loss += 0.5 * Rs * hw * hw * m->dy[j];
                /* longitudinal current from H_u at the two end nodes, half segment each */
                for (int e = 0; e < 2; e++) {
                    int jn = j + e;
                    /* H_u shares its location with E_v; on metal it is a normal
                     * component with no physical value in this formulation */
                    if (m->pec_v[IV(jn, kc)]) continue;
                    double hu = cabs(m->hu[IV(jn, kc)]);
                    loss += 0.5 * Rs * hu * hu * 0.5 * m->dy[j];
                }
            }
        }
    /* vertical surfaces: E_v segments (j, k) */
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k < nz; k++) {
            if (!m->pec_v[IV(j, k)]) continue;
            if (m->ps[1][j] > 0 || m->ps[2][k] > 0) continue;
            for (int side = 0; side < 2; side++) {
                int jc = side == 0 ? j : j - 1;
                if (jc < 0 || jc >= ny || SOLID(jc, k)) continue;
                cplx cw = (m->ev[IV(jc + 1, k)] - m->ev[IV(jc, k)]) / m->dy[jc] - (m->eu[IU(jc, k + 1)] - m->eu[IU(jc, k)]) / m->dz[k];
                double hw = cabs(cw) / (w * CEM_MU0);
                loss += 0.5 * Rs * hw * hw * m->dz[k];
                for (int e = 0; e < 2; e++) {
                    int kn = k + e;
                    if (m->pec_u[IU(jc, kn)]) continue;
                    double hv = cabs(m->hv[IU(jc, kn)]);
                    loss += 0.5 * Rs * hv * hv * 0.5 * m->dz[k];
                }
            }
        }
#undef SOLID
    double P = cem_mode2d_power(m);
    return P > 0 ? loss / (2 * P) : 0;
}
