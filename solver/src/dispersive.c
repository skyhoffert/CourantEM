/*
 * Courant EM - Debye dispersion by the auxiliary differential equation method.
 *
 * Each pole p contributes a polarisation current J_p with
 *     tau_p dJ_p/dt + J_p = eps0 deps_p dE/dt,
 * discretised at n+1/2 (trapezoidal):
 *     J_p^{n+1} = k_p J_p^n + beta_p (E^{n+1} - E^n) / dt,
 *     k_p = (1 - dt/2tau) / (1 + dt/2tau),  beta_p = eps0 deps_p dt / (tau (1 + dt/2tau)).
 * Ampere's law with J averaged over the step gives the usual update with
 *     eps -> eps0 eps_inf + sum_p beta_p / 2
 * in ca/cb, minus cb * sum_p (1 + k_p)/2 J_p^n. Poles share frequencies across
 * materials; edge weights are area averages of the surrounding cells, the same
 * rule used for eps and sigma. Design notes: docs/design/0005-dispersion.md
 */
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static double fmax_d(double a, double b) { return a > b ? a : b; }

int cem_set_debye_poles(cem_sim *s, int npole, const double *f_p) {
    if (s->finalized) { cem_set_err("already finalized"); return -1; }
    if (s->npole) { cem_set_err("Debye poles already set"); return -1; }
    if (npole < 1 || npole > MAXPOLES) { cem_set_err("bad pole count"); return -1; }
    s->npole = npole;
    for (int p = 0; p < npole; p++) {
        if (!(f_p[p] > 0)) { cem_set_err("pole frequencies must be positive"); return -1; }
        s->pole_f[p] = f_p[p];
        s->cell_deps[p] = calloc(s->N, sizeof(double));
        if (!s->cell_deps[p]) { cem_set_err("out of memory"); return -1; }
    }
    return 0;
}

int cem_material_box_debye(cem_sim *s, double x0, double x1, double y0, double y1, double z0, double z1,
                           double eps_inf, double sigma, const double *deps) {
    if (!s->npole) { cem_set_err("set Debye poles first"); return -1; }
    if (cem_material_box(s, x0, x1, y0, y1, z0, z1, eps_inf, sigma)) return -1;
    double lo[3] = {x0, y0, z0}, hi[3] = {x1, y1, z1};
    int r[3][2];
    for (int a = 0; a < 3; a++) {
        r[a][0] = s->n[a]; r[a][1] = -1;
        for (int i = 0; i < s->n[a]; i++) {
            double c = 0.5 * (s->node[a][i] + s->node[a][i + 1]);
            if (c >= lo[a] && c <= hi[a]) {
                if (i < r[a][0]) r[a][0] = i;
                if (i > r[a][1]) r[a][1] = i;
            }
        }
        if (r[a][1] < r[a][0]) return 0;
    }
    for (int i = r[0][0]; i <= r[0][1]; i++)
        for (int j = r[1][0]; j <= r[1][1]; j++)
            for (int k = r[2][0]; k <= r[2][1]; k++) {
                size_t q = IDX(s, i, j, k);
                for (int p = 0; p < s->npole; p++) s->cell_deps[p][q] = deps[p];
            }
    return 0;
}

/* Area-weighted average of the per-pole weights of the cells around an edge
 * along axis c. Returns 1 if any weight is non-zero. */
int cem_edge_debye(const cem_sim *s, int c, int i, int j, int k, double *deps) {
    if (!s->npole) return 0;
    int ua = (c + 1) % 3, va = (c + 2) % 3;
    int ijk[3] = {i, j, k};
    double wsum = 0;
    for (int p = 0; p < s->npole; p++) deps[p] = 0;
    for (int du = -1; du <= 0; du++)
        for (int dv = -1; dv <= 0; dv++) {
            int cc[3] = {ijk[0], ijk[1], ijk[2]};
            cc[ua] += du; cc[va] += dv;
            if (cc[ua] < 0 || cc[ua] >= s->n[ua] || cc[va] < 0 || cc[va] >= s->n[va]) continue;
            size_t q = IDX(s, cc[0], cc[1], cc[2]);
            double w = s->d[ua][cc[ua]] * s->d[va][cc[va]];
            for (int p = 0; p < s->npole; p++) deps[p] += w * s->cell_deps[p][q];
            wsum += w;
        }
    int any = 0;
    for (int p = 0; p < s->npole; p++) {
        deps[p] = wsum > 0 ? deps[p] / wsum : 0;
        if (deps[p] != 0) any = 1;
    }
    return any;
}

void cem_disp_save(cem_sim *s) {
    for (int e = 0; e < s->ndisp; e++) s->d_Eold[e] = s->E[s->d_comp[e]][s->d_idx[e]];
}

void cem_disp_correct(cem_sim *s) {
    int np = s->npole, nd = s->ndisp;
    for (int e = 0; e < nd; e++) {
        double acc = 0;
        for (int p = 0; p < np; p++) acc += 0.5 * (1 + s->d_k[p]) * s->d_J[(size_t)p * nd + e];
        s->E[s->d_comp[e]][s->d_idx[e]] -= s->d_cb[e] * acc;
    }
}

void cem_disp_update(cem_sim *s) {
    int np = s->npole, nd = s->ndisp;
    double idt = 1.0 / s->dt;
    for (int e = 0; e < nd; e++) {
        double dE = (s->E[s->d_comp[e]][s->d_idx[e]] - s->d_Eold[e]) * idt;
        for (int p = 0; p < np; p++) {
            size_t q = (size_t)p * nd + e;
            s->d_J[q] = s->d_k[p] * s->d_J[q] + s->d_beta[q] * dE;
        }
    }
}

void cem_disp_free(cem_sim *s) {
    for (int p = 0; p < s->npole; p++) free(s->cell_deps[p]);
    free(s->d_comp); free(s->d_idx); free(s->d_beta); free(s->d_J); free(s->d_Eold); free(s->d_cb);
}

/* Djordjevic-Sarkar model from Dk/Df at fref (1 kHz .. 1 THz relaxation band),
 * fitted by npole Debye poles log-spaced over [fmin/10, 10 fmax] with
 * non-negative weights (same method as tools/dielectric.py). Returns the
 * maximum relative error in the loss tangent over [fmin, fmax]. */
double cem_fit_djordjevic_sarkar(double dk, double df, double fref, double fmin, double fmax, int npole,
                                 double *eps_inf, double *f_p, double *deps) {
    const double m1 = 3, m2 = 12;
    double w1 = 2 * M_PI * pow(10, m1), w2 = 2 * M_PI * pow(10, m2);
#define LOGTERM(f, re, im) do { double w_ = 2 * M_PI * (f); \
        double nr = w2 * w2 + w_ * w_, dr = w1 * w1 + w_ * w_; \
        double ar = atan2(w_, w2) - atan2(w_, w1); \
        re = 0.5 * log(nr / dr) / log(10.0) / (m2 - m1); im = ar / log(10.0) / (m2 - m1); } while (0)
    double Lr, Li;
    LOGTERM(fref, Lr, Li);
    double d_eps = -dk * df / Li;          /* eps'' = -Im: target imag = -dk df */
    double e_inf_ds = dk - d_eps * Lr;
    int M = 200, K = npole + 1;
    double *A = malloc((size_t)2 * M * K * sizeof(double)), *b = malloc(2 * M * sizeof(double)), *x = malloc(K * sizeof(double));
    for (int p = 0; p < npole; p++) f_p[p] = fmin / 10 * pow(100 * fmax / fmin, npole > 1 ? (double)p / (npole - 1) : 0.5);
    double wr = 1 / dk, wi = 1 / fmax_d(dk * df, 1e-6);
    for (int i = 0; i < M; i++) {
        double f = fmin * pow(fmax / fmin, (double)i / (M - 1)), re, im;
        LOGTERM(f, re, im);
        double er = e_inf_ds + d_eps * re, ei = d_eps * im;
        b[i] = er * wr; b[M + i] = ei * wi;
        for (int p = 0; p < npole; p++) {
            double r = f / f_p[p], d = 1 + r * r;
            A[(size_t)p * 2 * M + i] = (1 / d) * wr;
            A[(size_t)p * 2 * M + M + i] = (-r / d) * wi;
        }
        A[(size_t)npole * 2 * M + i] = wr;
        A[(size_t)npole * 2 * M + M + i] = 0;
    }
    cem_nnls(2 * M, K, A, b, x);
    *eps_inf = x[npole];
    for (int p = 0; p < npole; p++) deps[p] = x[p];
    double worst = 0;
    for (int i = 0; i < M; i++) {
        double f = fmin * pow(fmax / fmin, (double)i / (M - 1)), re, im;
        LOGTERM(f, re, im);
        double er = e_inf_ds + d_eps * re, ei = d_eps * im;
        double fr = x[npole], fi = 0;
        for (int p = 0; p < npole; p++) {
            double r = f / f_p[p], d = 1 + r * r;
            fr += x[p] / d; fi += -x[p] * r / d;
        }
        double t0 = -ei / er, t1 = -fi / fr, e = fabs(t1 - t0) / t0;
        if (e > worst) worst = e;
    }
#undef LOGTERM
    free(A); free(b); free(x);
    return worst;
}
