/* Deterministic small scenario that exercises every kernel feature (non-uniform
 * grid, CPML, materials with loss, PEC box and sheets, two lumped ports, V and I
 * probes). Built natively and as WebAssembly; the printed numbers must agree.
 * The same scenario is rebuilt through the TypeScript wrapper in
 * web/test/parity.test.ts. */
#include "courant.h"

#include <math.h>
#include <stdio.h>

int main(void) {
    enum { NX = 24, NY = 18, NZ = 14 };
    double x[NX + 1], y[NY + 1], z[NZ + 1];
    for (int i = 0; i <= NX; i++) x[i] = 1e-4 * i + 2e-6 * i * i;
    for (int j = 0; j <= NY; j++) y[j] = 1e-4 * (j - NY / 2) * (1 + 0.01 * fabs(j - NY / 2.0));
    for (int k = 0; k <= NZ; k++) z[k] = 5e-5 * k + 1e-6 * k * k;
    cem_sim *s = cem_create(NX, NY, NZ, x, y, z);
    for (int f = 0; f < 6; f++)
        if (f != 4) cem_set_boundary(s, f, CEM_BC_CPML, 4);
    cem_material_box(s, -1, 1, -1, 1, 0, z[5], 3.4, 0.02);
    cem_pec_box(s, -1, 1, y[NY - 6], 1, 0, z[5]);
    cem_pec_sheet(s, CEM_Z, z[5], x[5], x[19], y[8], y[10]);
    int p1 = cem_add_lumped_port(s, CEM_Z, 5, 9, 0, 5, -1, 50);
    int p2 = cem_add_lumped_port(s, CEM_Z, 19, 9, 0, 5, -1, 50);
    cem_set_gaussian(s, p1, 1.0, 60e-12, 15e-12);
    double f[3] = {5e9, 20e9, 40e9};
    cem_set_frequencies(s, 3, f);
    int v = cem_add_vprobe(s, CEM_Z, 12, 9, 0, 5, -1);
    int ip = cem_add_iprobe(s, CEM_X, 12, 7, 11, 4, 5);
    if (cem_finalize(s, 0.95)) { printf("finalize failed: %s\n", cem_last_error()); return 1; }
    cem_step(s, 1500);
    printf("dt %.17g\n", cem_dt(s));
    printf("energy %.17g\n", cem_energy(s));
    const double *arr[4] = {cem_port_v(s, p1), cem_port_v(s, p2), cem_vprobe(s, v), cem_iprobe(s, ip)};
    const char *nm[4] = {"V1", "V2", "Vp", "Ip"};
    for (int a = 0; a < 4; a++)
        for (int q = 0; q < 3; q++) printf("%s[%d] %.17g %.17g\n", nm[a], q, arr[a][2 * q], arr[a][2 * q + 1]);
    cem_destroy(s);
    return 0;
}
