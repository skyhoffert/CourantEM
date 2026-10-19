/* Conductor-loss reference for issue 9: the GCPW cross-section with every metal
 * region a finite-conductivity volume, meshed finely enough to resolve the skin
 * depth, solved by the 2D eigenmode solver. Its attenuation is compared with the
 * perturbation estimate on the same geometry with the metal as PEC.
 *
 * Conductivity is lowered (default 1e6 S/m, delta = 5 um at 10 GHz) so the skin
 * depth can be resolved in 2D; the surface-impedance model in 3D is run with the
 * same conductivity for the comparison.
 *
 * Usage: condloss2d [--sigma S/m] [--f GHz] [--hs um] [--tcu um]
 */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double *v; int n, cap; } vec;
static void push(vec *a, double x) {
    if (a->n == a->cap) { a->cap = a->cap ? 2 * a->cap : 64; a->v = realloc(a->v, a->cap * sizeof(double)); }
    a->v[a->n++] = x;
}
static int cmpd(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }

/* graded lines: fixed points, refinement points with spacing hs growing by `ratio` up to hmax */
static vec graded(double *fix, int nfix, const double *ref, int nref, double hs, double hmax, double ratio) {
    qsort(fix, nfix, sizeof(double), cmpd);
    vec out = {0};
    push(&out, fix[0]);
    for (int s = 0; s + 1 < nfix; s++) {
        double a = fix[s], b = fix[s + 1];
        if (b - a < 1e-15) continue;
        vec seg = {0};
        double x = a;
        for (;;) {
            double d = 1e300, d2 = 1e300;
            for (int r = 0; r < nref; r++) if (fabs(x - ref[r]) < d) d = fabs(x - ref[r]);
            double h = fmin(hmax, hs + (ratio - 1) * d);
            for (int r = 0; r < nref; r++) if (fabs(x + h - ref[r]) < d2) d2 = fabs(x + h - ref[r]);
            h = fmin(h, fmin(hmax, hs + (ratio - 1) * d2));
            if (x + h >= b - 0.5 * h) break;
            x += h;
            push(&seg, x);
        }
        if (seg.n) {
            double last = seg.v[seg.n - 1], prev = seg.n > 1 ? last - seg.v[seg.n - 2] : last - a;
            double k = (b - a) / ((last - a) + prev);
            for (int i = 0; i < seg.n; i++) push(&out, a + (seg.v[i] - a) * k);
        }
        push(&out, b);
        free(seg.v);
    }
    return out;
}

int main(int argc, char **argv) {
    double sigma = 1e6, f = 10e9, hs = 1.2e-6, tcu = 35e-6;
    for (int i = 1; i + 1 < argc; i += 2) {
        double v = atof(argv[i + 1]);
        if (!strcmp(argv[i], "--sigma")) sigma = v;
        else if (!strcmp(argv[i], "--f")) f = v * 1e9;
        else if (!strcmp(argv[i], "--hs")) hs = v * 1e-6;
        else if (!strcmp(argv[i], "--tcu")) tcu = v * 1e-6;
    }
    const double h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6, er = 3.66;
    const double a = w / 2, b = a + s, yf = b + fence, ymax = yf + 0.5e-3, zbot = -12e-6, ztop = h + tcu + 1.0e-3;
    double w0 = 2 * M_PI * f, delta = sqrt(2 / (w0 * CEM_MU0 * sigma)), Rs = 1 / (sigma * delta);
    double yfix[] = {-ymax, -yf, -b, -a, 0, a, b, yf, ymax};
    double yref[] = {-yf, -b, -a, a, b, yf};
    double zfix[] = {zbot, 0, h, h + tcu, ztop};
    double zref[] = {0, h, h + tcu};
    vec Y = graded(yfix, 9, yref, 6, hs, 25e-6, 1.25);
    vec Z = graded(zfix, 5, zref, 3, hs, 25e-6, 1.25);
    int ny = Y.n - 1, nz = Z.n - 1;
    printf("sigma %.3g S/m, f %.2f GHz, skin depth %.2f um, Rs %.4f ohm; mesh %d x %d (min %.2f um)\n", sigma, f / 1e9,
           delta * 1e6, Rs, ny, nz, hs * 1e6);
    /* cells: material and metal flag */
    double *ce = malloc(ny * nz * sizeof(double));
    unsigned char *cm = calloc(ny * nz, 1);
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            double yc = 0.5 * (Y.v[j] + Y.v[j + 1]), zc = 0.5 * (Z.v[k] + Z.v[k + 1]), ay = fabs(yc);
            int metal = zc < 0 || (zc > h && zc < h + tcu && (ay < a || ay > b)) || (zc > 0 && zc < h && ay > yf);
            cm[j * nz + k] = (unsigned char)metal;
            ce[j * nz + k] = (zc > 0 && zc < h) ? er : 1.0;
        }
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    double *eu = calloc(nu, sizeof(double)), *ev = calloc(nv, sizeof(double)), *ew = calloc(nn, sizeof(double));
    double *su = calloc(nu, sizeof(double)), *sv = calloc(nv, sizeof(double)), *sw = calloc(nn, sizeof(double));
    unsigned char *pu = calloc(nu, 1), *pv = calloc(nv, 1), *pw = calloc(nn, 1);
    unsigned char *qu = calloc(nu, 1), *qv = calloc(nv, 1), *qw = calloc(nn, 1); /* PEC version */
    for (int pass = 0; pass < 3; pass++) {
        int NJ = pass == 0 ? ny : ny + 1, NK = pass == 1 ? nz : nz + 1;
        for (int j = 0; j < NJ; j++)
            for (int k = 0; k < NK; k++) {
                /* adjacent cells: E_u (j, k-1..k); E_v (j-1..j, k); E_w (j-1..j, k-1..k) */
                double es = 0, ss = 0, wsum = 0;
                int anym = 0;
                for (int dj = (pass == 0 ? 0 : -1); dj <= 0; dj++)
                    for (int dk = (pass == 1 ? 0 : -1); dk <= 0; dk++) {
                        int jj = j + dj, kk = k + dk;
                        if (jj < 0 || jj >= ny || kk < 0 || kk >= nz) continue;
                        double wgt = (pass == 0 ? 1 : (Y.v[jj + 1] - Y.v[jj])) * (pass == 1 ? 1 : (Z.v[kk + 1] - Z.v[kk]));
                        es += wgt * ce[jj * nz + kk];
                        ss += wgt * (cm[jj * nz + kk] ? sigma : 0);
                        wsum += wgt;
                        if (cm[jj * nz + kk]) anym = 1;
                    }
                double e = wsum > 0 ? es / wsum : 1, sg = wsum > 0 ? ss / wsum : 0;
                if (pass == 0) { int o = j * (nz + 1) + k; eu[o] = e; su[o] = sg; qu[o] = anym; }
                else if (pass == 1) { int o = j * nz + k; ev[o] = e; sv[o] = sg; qv[o] = anym; }
                else { int o = j * (nz + 1) + k; ew[o] = e; sw[o] = sg; qw[o] = anym; }
            }
    }
    /* lossy-volume solve */
    cem_mode2d *M = cem_mode2d_create(ny, nz, Y.v, Z.v);
    cem_mode2d_set_media(M, eu, ev, ew, pu, pv, pw);
    cem_mode2d_set_loss(M, su, sv, sw);
    int rc = cem_mode2d_solve(M, f, 1.64, NULL, 400, 1e-12);
    double alpha_vol = cem_mode2d_alpha(M), beta_vol = cem_mode2d_beta(M);
    /* PEC + perturbation on the same mesh */
    cem_mode2d *Q = cem_mode2d_create(ny, nz, Y.v, Z.v);
    cem_mode2d_set_media(Q, eu, ev, ew, qu, qv, qw);
    int rc2 = cem_mode2d_solve(Q, f, beta_vol / (w0 / CEM_C0), NULL, 400, 1e-12);
    double alpha_pert = cem_mode2d_conductor_alpha(Q, Rs);
    printf("volumetric conductors: rc %d, beta %.4f rad/m, alpha %.5f Np/m\n", rc, beta_vol, alpha_vol);
    printf("PEC + perturbation:    rc %d, beta %.4f rad/m, alpha_c %.5f Np/m (ratio pert/vol %.3f)\n", rc2,
           cem_mode2d_beta(Q), alpha_pert, alpha_pert / alpha_vol);
    cem_mode2d_destroy(M);
    cem_mode2d_destroy(Q);
    return 0;
}
