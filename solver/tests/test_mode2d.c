/* 2D mode solver test: WR-28 rectangular waveguide, TE10 at 35 GHz.
 * beta must match the Yee discrete dispersion (kc from sin(pi dy / 2a)) and the
 * perturbational copper loss must match the closed form
 *   alpha_c = Rs / (a^3 b beta k eta) (2 b pi^2 + a^3 k^2). */
#include "courant.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    const double a = 7.112e-3, b = 3.556e-3;
    enum { NY = 48, NZ = 24 };
    double y[NY + 1], z[NZ + 1];
    for (int j = 0; j <= NY; j++) y[j] = a * j / NY;
    for (int k = 0; k <= NZ; k++) z[k] = b * k / NZ;
    cem_mode2d *m = cem_mode2d_create(NY, NZ, y, z);
    double f = 35e9, w = 2 * M_PI * f, k0 = w / CEM_C0;
    double dy = a / NY, kc_d = 2 / dy * sin(M_PI / a * dy / 2), kc = M_PI / a;
    double beta_d = sqrt(k0 * k0 - kc_d * kc_d), beta = sqrt(k0 * k0 - kc * kc);
    int rc = cem_mode2d_solve(m, f, beta / k0, NULL, 300, 1e-13);
    double Rs = sqrt(w * CEM_MU0 / (2 * 5.8e7)), eta = 376.730313668;
    double ac = cem_mode2d_conductor_alpha(m, Rs);
    double an = Rs / (a * a * a * b * beta * k0 * eta) * (2 * b * M_PI * M_PI + a * a * a * k0 * k0);
    double eb = (cem_mode2d_beta(m) - beta_d) / beta_d, ea = (ac - an) / an;
    printf("  TE10: beta %.4f vs discrete %.4f (%+.1e), alpha_c %.5f vs %.5f Np/m (%+.2f%%)\n", cem_mode2d_beta(m),
           beta_d, eb, ac, an, 100 * ea);
    int ok = rc == 0 && fabs(eb) < 1e-8 && fabs(ea) < 0.01;
    printf(ok ? "test_mode2d passed\n" : "test_mode2d FAILED\n");
    cem_mode2d_destroy(m);
    return !ok;
}
