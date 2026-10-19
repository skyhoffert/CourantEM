/* Surface-impedance (copper) test: Q of a copper cavity.
 *
 * A 10 x 5 x 8 mm box is started in its TE101 mode. With every wall lossy
 * (sigma = 5.8e7 S/m), the stored energy decays as exp(-w t / Q). The measured
 * Q is compared with the closed form (Pozar, eq. 6.46):
 *   Q_c = (k a d)^3 b eta / (2 pi^2 Rs) / (2 a^3 b + 2 b d^3 + a^3 d + a d^3).
 * Require agreement within 3% (wall H sampled half a cell from the wall). */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int run(int rough, double rq) {
    const double a = 10e-3, b = 5e-3, d = 8e-3, sigma = 5.8e7;
    const int nx = 40, ny = 20, nz = 32;
    double x[nx + 1], y[ny + 1], z[nz + 1];
    for (int i = 0; i <= nx; i++) x[i] = a * i / nx;
    for (int j = 0; j <= ny; j++) y[j] = b * j / ny;
    for (int k = 0; k <= nz; k++) z[k] = d * k / nz;
    cem_sim *s = cem_create(nx, ny, nz, x, y, z);
    double f0 = CEM_C0 / 2 * sqrt(1 / (a * a) + 1 / (d * d));
    if (cem_set_conductor_loss(s, sigma, rough, rq, 0, f0 / 3, f0 * 3)) { printf("%s\n", cem_last_error()); return 1; }
    double fr = f0;
    cem_set_frequencies(s, 1, &fr);
    if (cem_finalize(s, 0.99)) { printf("finalize: %s\n", cem_last_error()); return 1; }
    printf("  Zs fit max error over band: %.3f%%\n", 100 * cem_conductor_fit_error(s));
    double *Ey = cem_field(s, CEM_Y);
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++)
                Ey[cem_index(s, i, j, k)] = sin(M_PI * x[i] / a) * sin(M_PI * z[k] / d);
    double dt = cem_dt(s);
    /* let the start-up transient pass, then fit ln W(t) over 20 ns */
    cem_step(s, (long)(1e-9 / dt));
    int nsamp = 40;
    long chunk = (long)(0.5e-9 / dt);
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int n = 0; n < nsamp; n++) {
        /* average energy over a period to remove the E/H exchange ripple */
        double W = 0;
        int nav = 16;
        long per = (long)(1 / f0 / dt / nav);
        for (int q = 0; q < nav; q++) { cem_step(s, per); W += cem_energy(s) / nav; }
        cem_step(s, chunk - per * nav);
        double t = cem_step_count(s) * dt, lw = log(W);
        sx += t; sy += lw; sxx += t * t; sxy += t * lw;
    }
    double slope = (nsamp * sxy - sx * sy) / (nsamp * sxx - sx * sx);
    double w0 = 2 * M_PI * f0;
    double Qm = -w0 / slope;
    double Rs = sqrt(w0 * CEM_MU0 / (2 * sigma)), k = w0 / CEM_C0, eta = 376.730313668;
    double Qa = pow(k * a * d, 3) * b * eta / (2 * M_PI * M_PI * Rs) /
                (2 * pow(a, 3) * b + 2 * b * pow(d, 3) + pow(a, 3) * d + a * pow(d, 3));
    if (rough == 1) { /* Groisse: Rs scales by K */
        double delta = sqrt(2 / (w0 * CEM_MU0 * sigma));
        Qa /= 1 + exp(-pow(delta / (2 * rq), 1.6));
    }
    double err = (Qm - Qa) / Qa;
    int ok = fabs(err) < 0.03;
    printf("  TE101 %.3f GHz %s: Q measured %.1f, analytic %.1f, error %+.2f%%\n", f0 / 1e9,
           rough ? "Groisse Rq 1 um" : "smooth", Qm, Qa, 100 * err);
    cem_destroy(s);
    return !ok;
}

int main(void) {
    int fail = run(0, 0) | run(1, 1e-6);
    printf(fail ? "test_sibc FAILED\n" : "test_sibc passed\n");
    return fail;
}
