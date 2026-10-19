/* CPML reflection test.
 *
 * A short lumped dipole radiates in free space. The field at a probe near the
 * CPML is compared with a reference run in a much larger PEC box, over a time
 * window that ends before the reference box's wall reflection reaches the probe.
 * The difference is the CPML's reflection; require < 1% (-40 dB) of the peak. */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define H 1e-3

static double *run(int interior, int npml, int nsteps, double *dt_out) {
    int n = interior + 2 * npml;
    double *c = malloc((n + 1) * sizeof(double));
    for (int i = 0; i <= n; i++) c[i] = (i - n / 2.0) * H;
    cem_sim *s = cem_create(n, n, n, c, c, c);
    if (npml)
        for (int f = 0; f < 6; f++) cem_set_boundary(s, f, CEM_BC_CPML, npml);
    int m = n / 2;
    int p = cem_add_lumped_port(s, CEM_Z, m, m, m, m + 1, 1, 50.0);
    double tau = 15e-12;
    cem_set_gaussian(s, p, 1.0, 4 * tau, tau);
    double fr = 1e9;
    cem_set_frequencies(s, 1, &fr);
    cem_finalize(s, 0.99);
    *dt_out = cem_dt(s);
    double *rec = malloc(nsteps * sizeof(double));
    double *Ez = cem_field(s, CEM_Z);
    /* probe 4 cells inside the CPML interface along x, diagonal offset in y */
    size_t q = cem_index(s, m + interior / 2 - 4, m + 3, m);
    for (int t = 0; t < nsteps; t++) {
        cem_step(s, 1);
        rec[t] = Ez[q];
    }
    cem_destroy(s);
    free(c);
    return rec;
}

int main(void) {
    const int interior = 30, npml = 10;
    double dt;
    /* window: reflection from the reference walls must not arrive. Reference box
     * half size 64 cells; path source->wall->probe >= 64 + (64-11) = 117 cells. */
    int nsteps = (int)(100 * H / CEM_C0 / (0.99 * H / (CEM_C0 * sqrt(3.0))));
    double *a = run(interior, npml, nsteps, &dt);
    /* reference: same cell indexing near the source, interior 128, no PML */
    int ref_int = 128;
    double dt2;
    double *b = NULL;
    {
        int n = ref_int;
        double *c = malloc((n + 1) * sizeof(double));
        for (int i = 0; i <= n; i++) c[i] = (i - n / 2.0) * H;
        cem_sim *s = cem_create(n, n, n, c, c, c);
        int m = n / 2;
        int p = cem_add_lumped_port(s, CEM_Z, m, m, m, m + 1, 1, 50.0);
        double tau = 15e-12;
        cem_set_gaussian(s, p, 1.0, 4 * tau, tau);
        double fr = 1e9;
        cem_set_frequencies(s, 1, &fr);
        cem_finalize(s, 0.99);
        dt2 = cem_dt(s);
        b = malloc(nsteps * sizeof(double));
        double *Ez = cem_field(s, CEM_Z);
        size_t q = cem_index(s, m + interior / 2 - 4, m + 3, m);
        for (int t = 0; t < nsteps; t++) {
            cem_step(s, 1);
            b[t] = Ez[q];
        }
        cem_destroy(s);
        free(c);
    }
    double peak = 0, err = 0;
    for (int t = 0; t < nsteps; t++) {
        if (fabs(b[t]) > peak) peak = fabs(b[t]);
        if (fabs(a[t] - b[t]) > err) err = fabs(a[t] - b[t]);
    }
    double rel = err / peak;
    printf("  steps %d, dt %.3e (ref %.3e), peak %.3e, max diff %.3e, reflection %.1f dB\n", nsteps, dt, dt2,
           peak, err, 20 * log10(rel));
    int ok = rel < 1e-2;
    printf(ok ? "test_cpml passed\n" : "test_cpml FAILED\n");
    free(a);
    free(b);
    return !ok;
}
