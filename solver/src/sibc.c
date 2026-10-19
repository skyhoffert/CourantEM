/*
 * Courant EM - surface-impedance boundary condition (SIBC) for copper.
 *
 * Every PEC edge with field on at least one side becomes a lossy surface edge:
 *     E_t = Z_s(w) J_s,   J_s = n x H   (H taken half a cell outside the surface)
 * with Z_s(w) = (1 + j) sqrt(w mu0 / (2 sigma)) * K(w) for a good conductor and a
 * roughness factor K (none, Groisse or Huray snowball).
 *
 * Time domain: Z_s is fitted over the band by a passive Foster form
 *     Z(s) ~ sum_k c_k s / (s + a_k),  c_k >= 0, a_k log-spaced,
 * (non-negative least squares), and each term is advanced by recursive
 * convolution:  y_k' = a_k (J - y_k),  term_k = c_k (J - y_k).
 * The skin depth (0.5 um at 20 GHz) is far below any cell, so the conductor
 * interior is never meshed. Design notes: docs/design/0006-conductor-loss.md
 */
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NZP 16

/* ---- target impedance and NNLS fit -------------------------------------- */

static void zs_target(double f, double sigma, int rough, double p1, double p2, double *re, double *im) {
    double w = 2 * M_PI * f;
    double rs = sqrt(w * CEM_MU0 / (2 * sigma));
    double K = 1;
    double delta = sqrt(2 / (w * CEM_MU0 * sigma));
    if (rough == 1) {
        /* Groisse: K = 1 + exp(-(delta / (2 Rq))^1.6), p1 = Rq (rms roughness, m) */
        K = 1 + exp(-pow(delta / (2 * p1), 1.6));
    } else if (rough == 2) {
        /* Huray single-size snowball: K = 1 + (p2 * 4 pi a^2) / (1 + delta/a + delta^2/(2 a^2)),
         * p1 = sphere radius a (m), p2 = spheres per unit flat area (1/m^2) */
        double a = p1;
        K = 1 + p2 * 4 * M_PI * a * a / (1 + delta / a + delta * delta / (2 * a * a));
    }
    *re = rs * K;
    *im = rs * K;
}

/* Lawson-Hanson NNLS for small dense problems (m x n, column-major A). */
void cem_nnls(int m, int n, const double *A, const double *b, double *x) {
    int *P = calloc(n, sizeof(int));
    double *w = malloc(n * sizeof(double)), *z = malloc(n * sizeof(double)), *r = malloc(m * sizeof(double));
    double *M = malloc((size_t)n * n * sizeof(double)), *g = malloc(n * sizeof(double));
    for (int j = 0; j < n; j++) x[j] = 0;
    for (int outer = 0; outer < 3 * n; outer++) {
        for (int i = 0; i < m; i++) {
            r[i] = b[i];
            for (int j = 0; j < n; j++) r[i] -= A[(size_t)j * m + i] * x[j];
        }
        int jmax = -1;
        double wmax = 1e-14;
        for (int j = 0; j < n; j++) {
            w[j] = 0;
            for (int i = 0; i < m; i++) w[j] += A[(size_t)j * m + i] * r[i];
            if (!P[j] && w[j] > wmax) { wmax = w[j]; jmax = j; }
        }
        if (jmax < 0) break;
        P[jmax] = 1;
        for (int inner = 0; inner < 3 * n; inner++) {
            /* least squares on the passive set via normal equations + Cholesky-free Gauss */
            int idx[64], np = 0;
            for (int j = 0; j < n; j++) if (P[j]) idx[np++] = j;
            for (int a = 0; a < np; a++) {
                g[a] = 0;
                for (int i = 0; i < m; i++) g[a] += A[(size_t)idx[a] * m + i] * b[i];
                for (int c = 0; c < np; c++) {
                    double s = 0;
                    for (int i = 0; i < m; i++) s += A[(size_t)idx[a] * m + i] * A[(size_t)idx[c] * m + i];
                    M[a * np + c] = s;
                }
            }
            for (int a = 0; a < np; a++) { /* Gaussian elimination with pivoting */
                int p = a;
                for (int q = a + 1; q < np; q++) if (fabs(M[q * np + a]) > fabs(M[p * np + a])) p = q;
                if (p != a) {
                    for (int c = 0; c < np; c++) { double t = M[a * np + c]; M[a * np + c] = M[p * np + c]; M[p * np + c] = t; }
                    double t = g[a]; g[a] = g[p]; g[p] = t;
                }
                for (int q = a + 1; q < np; q++) {
                    double l = M[q * np + a] / M[a * np + a];
                    for (int c = a; c < np; c++) M[q * np + c] -= l * M[a * np + c];
                    g[q] -= l * g[a];
                }
            }
            for (int a = np - 1; a >= 0; a--) {
                double s = g[a];
                for (int c = a + 1; c < np; c++) s -= M[a * np + c] * z[idx[c]];
                z[idx[a]] = s / M[a * np + a];
            }
            int ok = 1;
            for (int a = 0; a < np; a++) if (z[idx[a]] <= 0) ok = 0;
            if (ok) { for (int a = 0; a < np; a++) x[idx[a]] = z[idx[a]]; break; }
            double alpha = 1;
            for (int a = 0; a < np; a++) {
                int j = idx[a];
                if (z[j] <= 0) { double t = x[j] / (x[j] - z[j]); if (t < alpha) alpha = t; }
            }
            for (int a = 0; a < np; a++) {
                int j = idx[a];
                x[j] += alpha * (z[j] - x[j]);
                if (x[j] <= 1e-300) { x[j] = 0; P[j] = 0; }
            }
        }
    }
    free(P); free(w); free(z); free(r); free(M); free(g);
}

/* ---- API ----------------------------------------------------------------- */

int cem_set_conductor_loss(cem_sim *s, double sigma, int rough, double p1, double p2, double fmin, double fmax) {
    if (s->finalized) { cem_set_err("already finalized"); return -1; }
    if (!(sigma > 0) || !(fmax > fmin) || !(fmin > 0)) { cem_set_err("bad conductor-loss parameters"); return -1; }
    s->sibc_on = 1;
    s->sibc_sigma = sigma;
    s->sibc_rough = rough; s->sibc_p1 = p1; s->sibc_p2 = p2;
    /* fit: poles one decade beyond the band each side */
    int K = 12, M = 200;
    double *A = malloc((size_t)2 * M * K * sizeof(double)), *b = malloc(2 * M * sizeof(double));
    for (int k = 0; k < K; k++) s->sibc_a[k] = 2 * M_PI * fmin / 10 * pow(100 * fmax / fmin, (double)k / (K - 1));
    for (int i = 0; i < M; i++) {
        double f = fmin * pow(fmax / fmin, (double)i / (M - 1)), w = 2 * M_PI * f, re, im;
        zs_target(f, sigma, rough, p1, p2, &re, &im);
        double scale = 1 / hypot(re, im); /* relative error */
        b[i] = re * scale; b[M + i] = im * scale;
        for (int k = 0; k < K; k++) {
            double a = s->sibc_a[k];
            /* s/(s+a) at s = jw: (w^2 + j w a) / (w^2 + a^2) */
            double d = w * w + a * a;
            A[(size_t)k * 2 * M + i] = w * w / d * scale;
            A[(size_t)k * 2 * M + M + i] = w * a / d * scale;
        }
    }
    cem_nnls(2 * M, K, A, b, s->sibc_c);
    s->sibc_npole = K;
    /* fit quality */
    double worst = 0;
    for (int i = 0; i < M; i++) {
        double f = fmin * pow(fmax / fmin, (double)i / (M - 1)), w = 2 * M_PI * f, re, im, fr = 0, fi = 0;
        zs_target(f, sigma, rough, p1, p2, &re, &im);
        for (int k = 0; k < K; k++) {
            double a = s->sibc_a[k], d = w * w + a * a;
            fr += s->sibc_c[k] * w * w / d;
            fi += s->sibc_c[k] * w * a / d;
        }
        double e = hypot(fr - re, fi - im) / hypot(re, im);
        if (e > worst) worst = e;
    }
    s->sibc_fit_err = worst;
    free(A); free(b);
    return 0;
}

double cem_conductor_fit_error(const cem_sim *s) { return s->sibc_fit_err; }

/* Surface impedance of the fitted model at frequency f (re, im). */
void cem_conductor_zs(const cem_sim *s, double f, double *re, double *im) {
    double w = 2 * M_PI * f;
    *re = *im = 0;
    for (int k = 0; k < s->sibc_npole; k++) {
        double a = s->sibc_a[k], d = w * w + a * a;
        *re += s->sibc_c[k] * w * w / d;
        *im += s->sibc_c[k] * w * a / d;
    }
}

/* Is the E edge along c at ijk metal (PEC)? Outer-boundary tangential edges
 * count as metal when that face is a PEC boundary. Out-of-range -> metal. */
static int edge_is_metal(const cem_sim *s, int c, const int *ijk) {
    for (int a = 0; a < 3; a++) {
        int lim = a == c ? s->n[a] - 1 : s->n[a];
        if (ijk[a] < 0 || ijk[a] > lim) return 1;
    }
    for (int a = 0; a < 3; a++) {
        if (a == c) continue;
        if (ijk[a] == 0 && s->bc[2 * a] == CEM_BC_PEC) return 1;
        if (ijk[a] == s->n[a] && s->bc[2 * a + 1] == CEM_BC_PEC) return 1;
    }
    return s->cb[c][IDX(s, ijk[0], ijk[1], ijk[2])] == 0;
}

int cem_sibc_finalize(cem_sim *s, const unsigned char *is_port) {
    if (!s->sibc_on) return 0;
    int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    for (int pass = 0; pass < 2; pass++) {
        int ne = 0;
        for (int c = 0; c < 3; c++) {
            int ua = (c + 1) % 3, va = (c + 2) % 3;
            for (int i = 0; i <= nx; i++)
                for (int j = 0; j <= ny; j++)
                    for (int k = 0; k <= nz; k++) {
                        int ijk[3] = {i, j, k};
                        if (ijk[c] >= s->n[c]) continue;
                        size_t q = IDX(s, i, j, k);
                        if (!edge_is_metal(s, c, ijk)) continue;
                        if (is_port && is_port[c * s->N + q]) continue;
                        /* edges lying in a CPML slab stay PEC */
                        int in_pml = 0;
                        for (int a = 0; a < 3; a++) {
                            if (a == c) { if (ijk[a] < s->npml[2 * a] || ijk[a] >= s->n[a] - s->npml[2 * a + 1]) in_pml = 1; }
                            else if (ijk[a] < s->npml[2 * a] || ijk[a] > s->n[a] - s->npml[2 * a + 1]) in_pml = 1;
                        }
                        if (in_pml) continue;
                        /* faces: a non-metal neighbour edge parallel to this one */
                        int nf = 0, hc[4];
                        size_t hi[4];
                        double sg[4];
                        int nb[3];
                        /* +v: J_c = -H_u(v+1/2) */
                        memcpy(nb, ijk, sizeof nb); nb[va] += 1;
                        if (!edge_is_metal(s, c, nb)) { hc[nf] = ua; hi[nf] = q; sg[nf] = -1; nf++; }
                        /* -v: J_c = +H_u(v-1/2) */
                        memcpy(nb, ijk, sizeof nb); nb[va] -= 1;
                        if (!edge_is_metal(s, c, nb)) { int h3[3] = {ijk[0], ijk[1], ijk[2]}; h3[va] -= 1; hc[nf] = ua; hi[nf] = IDX(s, h3[0], h3[1], h3[2]); sg[nf] = 1; nf++; }
                        /* +u: J_c = +H_v(u+1/2) */
                        memcpy(nb, ijk, sizeof nb); nb[ua] += 1;
                        if (!edge_is_metal(s, c, nb)) { hc[nf] = va; hi[nf] = q; sg[nf] = 1; nf++; }
                        /* -u: J_c = -H_v(u-1/2) */
                        memcpy(nb, ijk, sizeof nb); nb[ua] -= 1;
                        if (!edge_is_metal(s, c, nb)) { int h3[3] = {ijk[0], ijk[1], ijk[2]}; h3[ua] -= 1; hc[nf] = va; hi[nf] = IDX(s, h3[0], h3[1], h3[2]); sg[nf] = -1; nf++; }
                        if (!nf) continue;
                        if (pass == 1) {
                            s->sb_comp[ne] = c;
                            s->sb_idx[ne] = q;
                            s->sb_nf[ne] = nf;
                            for (int f = 0; f < nf; f++) {
                                s->sb_hc[4 * ne + f] = hc[f];
                                s->sb_hi[4 * ne + f] = hi[f];
                                s->sb_sg[4 * ne + f] = sg[f] / nf; /* average over faces */
                            }
                        }
                        ne++;
                    }
        }
        if (pass == 0) {
            s->nsb = ne;
            s->sb_comp = calloc(ne ? ne : 1, sizeof(int));
            s->sb_idx = calloc(ne ? ne : 1, sizeof(size_t));
            s->sb_nf = calloc(ne ? ne : 1, sizeof(int));
            s->sb_hc = calloc(4 * (ne ? ne : 1), sizeof(int));
            s->sb_hi = calloc(4 * (ne ? ne : 1), sizeof(size_t));
            s->sb_sg = calloc(4 * (ne ? ne : 1), sizeof(double));
            s->sb_y = calloc((size_t)(ne ? ne : 1) * s->sibc_npole, sizeof(double));
            if (!s->sb_y) { cem_set_err("out of memory"); return -1; }
        }
    }
    for (int k = 0; k < s->sibc_npole; k++) s->sibc_e[k] = exp(-s->sibc_a[k] * s->dt);
    return 0;
}

/* After the E update: E = sum_k c_k (J - y_k), with y_k advanced to n+1 from J at n+1/2. */
void cem_sibc_apply(cem_sim *s) {
    int K = s->sibc_npole, ne = s->nsb;
    for (int e = 0; e < ne; e++) {
        double J = 0;
        for (int f = 0; f < s->sb_nf[e]; f++) J += s->sb_sg[4 * e + f] * s->H[s->sb_hc[4 * e + f]][s->sb_hi[4 * e + f]];
        double E = 0;
        double *y = &s->sb_y[(size_t)e * K];
        for (int k = 0; k < K; k++) {
            y[k] = s->sibc_e[k] * y[k] + (1 - s->sibc_e[k]) * J;
            E += s->sibc_c[k] * (J - y[k]);
        }
        s->E[s->sb_comp[e]][s->sb_idx[e]] = E;
    }
}

void cem_sibc_free(cem_sim *s) {
    free(s->sb_comp); free(s->sb_idx); free(s->sb_nf); free(s->sb_hc); free(s->sb_hi); free(s->sb_sg); free(s->sb_y);
}
