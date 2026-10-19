/* Checkpoint round trip: run A for 900 steps; run B for 400, save, rebuild C,
 * load, run C for 500 more. C's fields and spectra must equal A's bit for bit.
 * Uses every stateful feature: CPML, Debye poles, copper SIBC, lumped ports, probes. */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static cem_sim *scene(void) {
    enum { NX = 24, NY = 18, NZ = 14 };
    double x[NX + 1], y[NY + 1], z[NZ + 1];
    for (int i = 0; i <= NX; i++) x[i] = 1e-4 * i + 2e-6 * i * i;
    for (int j = 0; j <= NY; j++) y[j] = 1e-4 * (j - NY / 2);
    for (int k = 0; k <= NZ; k++) z[k] = 5e-5 * k + 1e-6 * k * k;
    cem_sim *s = cem_create(NX, NY, NZ, x, y, z);
    for (int f = 0; f < 6; f++) if (f != 4) cem_set_boundary(s, f, CEM_BC_CPML, 4);
    double fp[3] = {1e9, 10e9, 100e9}, de[3] = {0.02, 0.03, 0.02};
    cem_set_debye_poles(s, 3, fp);
    cem_material_box_debye(s, -1, 1, -1, 1, 0, z[5], 3.4, 0.0, de);
    cem_set_conductor_loss(s, 5.8e7, 0, 0, 0, 1e9, 40e9);
    cem_pec_sheet(s, CEM_Z, z[5], x[5], x[19], y[8], y[10]);
    int p1 = cem_add_lumped_port(s, CEM_Z, 5, 9, 0, 5, -1, 50);
    cem_add_lumped_port(s, CEM_Z, 19, 9, 0, 5, -1, 50);
    cem_set_gaussian(s, p1, 1.0, 60e-12, 15e-12);
    double f[3] = {5e9, 20e9, 40e9};
    cem_set_frequencies(s, 3, f);
    cem_add_vprobe(s, CEM_Z, 12, 9, 0, 5, -1);
    cem_add_iprobe(s, CEM_X, 12, 7, 11, 4, 5);
    if (cem_finalize(s, 0.95)) { printf("finalize: %s\n", cem_last_error()); exit(1); }
    return s;
}

int main(void) {
    cem_sim *A = scene(), *B = scene(), *C = scene();
    cem_step(A, 900);
    cem_step(B, 400);
    size_t n = cem_checkpoint_size(B);
    double *buf = malloc(n * sizeof(double));
    cem_checkpoint_save(B, buf);
    if (cem_checkpoint_load(C, buf)) { printf("load: %s\n", cem_last_error()); return 1; }
    cem_step(C, 500);
    int fail = 0;
    int nx, ny, nz;
    cem_dims(A, &nx, &ny, &nz);
    size_t N = (size_t)(nx + 1) * (ny + 1) * (nz + 1);
    for (int c = 0; c < 6; c++)
        if (memcmp(cem_field(A, c), cem_field(C, c), N * sizeof(double))) { fail = 1; printf("  field %d differs\n", c); }
    if (memcmp(cem_port_v(A, 1), cem_port_v(C, 1), 6 * sizeof(double))) { fail = 1; printf("  port spectra differ\n"); }
    if (memcmp(cem_vprobe(A, 0), cem_vprobe(C, 0), 6 * sizeof(double))) { fail = 1; printf("  probe spectra differ\n"); }
    printf("  checkpoint: %zu doubles (%.1f MB); resumed run identical: %s\n", n, n * 8 / 1e6, fail ? "no" : "yes");
    printf(fail ? "test_checkpoint FAILED\n" : "test_checkpoint passed\n");
    return fail;
}
