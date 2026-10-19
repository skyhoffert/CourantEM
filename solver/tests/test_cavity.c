/* PEC cavity resonance test.
 *
 * Uniform grid: the Yee scheme's discrete dispersion relation predicts the
 * resonant frequency of each cavity mode exactly, so the measured peaks must
 * match it to within the DFT peak-location error (~1e-4 relative).
 * Non-uniform grid: compare against the continuous resonance within 1%. */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double peak_freq(cem_sim *s, int probe, const double *f, int nf) {
    const double *X = cem_vprobe(s, probe);
    int best = 0;
    double bm = 0;
    for (int i = 0; i < nf; i++) {
        double m = hypot(X[2 * i], X[2 * i + 1]);
        if (m > bm) { bm = m; best = i; }
    }
    if (best == 0 || best == nf - 1) return f[best];
    double a = hypot(X[2 * best - 2], X[2 * best - 1]), b = bm, c = hypot(X[2 * best + 2], X[2 * best + 3]);
    double off = 0.5 * (a - c) / (a - 2 * b + c);
    return f[best] + off * (f[1] - f[0]);
}

static int run(int nonuniform) {
    const double a = 10e-3, b = 5e-3, d = 8e-3;
    const int nx = 40, ny = 20, nz = 32;
    double x[nx + 1], y[ny + 1], z[nz + 1];
    for (int i = 0; i <= nx; i++) {
        double u = (double)i / nx;
        x[i] = a * (nonuniform ? u + 0.08 * sin(2 * M_PI * u) / (2 * M_PI) * 2.0 : u);
    }
    for (int j = 0; j <= ny; j++) y[j] = b * j / ny;
    for (int k = 0; k <= nz; k++) {
        double u = (double)k / nz;
        z[k] = d * (nonuniform ? u - 0.1 * sin(2 * M_PI * u) / (2 * M_PI) * 2.0 : u);
    }
    cem_sim *s = cem_create(nx, ny, nz, x, y, z);
    if (!s) { printf("create failed: %s\n", cem_last_error()); return 1; }

    /* modes with Ey only: TE_m0p, kx = m pi/a, kz = p pi/d */
    int modes[3][2] = {{1, 1}, {2, 1}, {1, 2}};
    double fc[3];
    for (int m = 0; m < 3; m++)
        fc[m] = CEM_C0 / 2 * sqrt(pow(modes[m][0] / a, 2) + pow(modes[m][1] / d, 2));

    const int nf = 4001;
    double *f = malloc(nf * sizeof(double));
    for (int i = 0; i < nf; i++) f[i] = 15e9 + 30e9 * i / (nf - 1);
    cem_set_frequencies(s, nf, f);
    int pr = cem_add_vprobe(s, CEM_Y, 13, 9, 7, 8, 1); /* Ey edge at x node 13, z node 9 (cyclic: z, x) */
    if (pr < 0) { printf("probe: %s\n", cem_last_error()); return 1; }
    if (cem_finalize(s, 0.99)) { printf("finalize: %s\n", cem_last_error()); return 1; }

    /* initial Ey: off-centre Gaussian bump (excites all Ey modes) */
    double *Ey = cem_field(s, CEM_Y);
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++) {
                double dx = x[i] - 0.31 * a, dz = z[k] - 0.37 * d;
                Ey[cem_index(s, i, j, k)] = exp(-(dx * dx + dz * dz) / (2 * 1.5e-3 * 1.5e-3));
            }
    double dt = cem_dt(s);
    long steps = (long)(20e-9 / dt);
    cem_step(s, steps);

    int fail = 0;
    for (int m = 0; m < 3; m++) {
        /* search window +-3% around the expected mode */
        int i0 = (int)((fc[m] * 0.97 - f[0]) / (f[1] - f[0])), i1 = (int)((fc[m] * 1.03 - f[0]) / (f[1] - f[0]));
        if (i0 < 0) i0 = 0;
        if (i1 > nf - 1) i1 = nf - 1;
        /* temporarily view a window by offsetting pointers */
        const double *X = cem_vprobe(s, pr);
        int best = i0;
        double bm = 0;
        for (int i = i0; i <= i1; i++) {
            double mm = hypot(X[2 * i], X[2 * i + 1]);
            if (mm > bm) { bm = mm; best = i; }
        }
        double A = hypot(X[2 * best - 2], X[2 * best - 1]), C = hypot(X[2 * best + 2], X[2 * best + 3]);
        double fpk = f[best] + 0.5 * (A - C) / (A - 2 * bm + C) * (f[1] - f[0]);
        double ref, tol;
        if (!nonuniform) {
            double dx = a / nx, dz = d / nz;
            double kx = modes[m][0] * M_PI / a, kz = modes[m][1] * M_PI / d;
            double S = CEM_C0 * dt * sqrt(pow(sin(kx * dx / 2) / dx, 2) + pow(sin(kz * dz / 2) / dz, 2));
            ref = 2 / dt * asin(S) / (2 * M_PI);
            tol = 2e-4;
        } else {
            ref = fc[m];
            tol = 1e-2;
        }
        double err = (fpk - ref) / ref;
        int ok = fabs(err) < tol;
        if (!ok) fail = 1;
        printf("  TE%d0%d: measured %.6f GHz, expected %.6f GHz (continuous %.6f), rel err %+.2e %s\n",
               modes[m][0], modes[m][1], fpk / 1e9, ref / 1e9, fc[m] / 1e9, err, ok ? "PASS" : "FAIL");
    }
    (void)peak_freq;
    free(f);
    cem_destroy(s);
    return fail;
}

int main(void) {
    int fail = 0;
    printf("uniform grid:\n");
    fail |= run(0);
    printf("non-uniform grid:\n");
    fail |= run(1);
    printf(fail ? "test_cavity FAILED\n" : "test_cavity passed\n");
    return fail;
}
