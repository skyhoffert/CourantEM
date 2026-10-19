/* Build step 1: hardcoded GCPW line between two 50-ohm lumped ports.
 *
 * Geometry (defaults, RO4350B-like, lossless):
 *   substrate eps_r 3.66, h = 254 um; strip w = 500 um; gaps s = 200 um;
 *   via fences idealised as solid PEC walls 300 um beyond each gap edge;
 *   bottom ground = PEC domain boundary; top grounds extend to the CPML.
 *   Line length L between the port columns; board continues 0.5 mm past each
 *   port into the x CPML (strip ends at the ports).
 * Outputs (in --out dir):
 *   sparams.csv   f, S11 re/im, S21 re/im (port 1 excited)
 *   line.s2p      Touchstone, S22 = S11 and S12 = S21 by symmetry
 *   probes.csv    f, then V and I probe spectra (re/im) along the line
 *   meta.json     mesh, timing, probe positions
 *
 * Usage: gcpw_line [--hmin um] [--hmax um] [--dx um] [--len mm] [--er v]
 *                  [--fmax GHz] [--df GHz] [--decay v] [--out dir]
 */
#include "courant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double *v; int n, cap; } vec;
static void push(vec *a, double x) {
    if (a->n == a->cap) { a->cap = a->cap ? 2 * a->cap : 64; a->v = realloc(a->v, a->cap * sizeof(double)); }
    a->v[a->n++] = x;
}
static int cmpd(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* Graded 1D mesh: fixed lines at `fix` (always kept); `ref` are refinement points
 * where spacing is hmin, growing by `ratio` per cell up to hmax. */
static vec graded(const double *fix, int nfix, const double *ref, int nref, double hmin, double hmax, double ratio) {
    double *pts = malloc(nfix * sizeof(double));
    memcpy(pts, fix, nfix * sizeof(double));
    qsort(pts, nfix, sizeof(double), cmpd);
    vec out = {0};
    push(&out, pts[0]);
    for (int s = 0; s + 1 < nfix; s++) {
        double a = pts[s], b = pts[s + 1];
        if (b - a < 1e-15) continue;
        vec seg = {0};
        double x = a;
        for (;;) {
            double d = 1e300;
            for (int r = 0; r < nref; r++) if (fabs(x - ref[r]) < d) d = fabs(x - ref[r]);
            double h = fmin(hmax, hmin + (ratio - 1) * d);
            double d2 = 1e300;
            for (int r = 0; r < nref; r++) if (fabs(x + h - ref[r]) < d2) d2 = fabs(x + h - ref[r]);
            h = fmin(h, fmin(hmax, hmin + (ratio - 1) * d2));
            if (x + h >= b - 0.5 * h) break;
            x += h;
            push(&seg, x);
        }
        if (seg.n) {
            /* the remainder b - last lies in (0.5h, 1.5h]; scale the segment
             * uniformly so the final gap equals its neighbour */
            double last = seg.v[seg.n - 1];
            double prev = seg.n > 1 ? last - seg.v[seg.n - 2] : last - a;
            double k = (b - a) / ((last - a) + prev);
            for (int i = 0; i < seg.n; i++) push(&out, a + (seg.v[i] - a) * k);
        }
        push(&out, b);
        free(seg.v);
    }
    free(pts);
    return out;
}

static int find_node(const vec *a, double x) {
    int b = 0;
    for (int i = 1; i < a->n; i++) if (fabs(a->v[i] - x) < fabs(a->v[b] - x)) b = i;
    return b;
}

int main(int argc, char **argv) {
    double hmin = 10e-6, hmax = 100e-6, dxl = 50e-6, len = 10e-3, er = 3.66, fmax = 40e9, df = 0.1e9;
    double decay = 1e-6, mon_mm = 2.5, feed_mm = 4.0, yair_mm = 1.0, zair_mm = 1.5;
    double cu_sigma = 0, rq_um = 0, tcu_um = 0;
    int mesh_only = 0, closed = 0, snap = 0, modes_only = 0, modeports = 0;
    const char *out = "run_gcpw";
    const char *debye = NULL;
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "--mesh-only") || !strcmp(argv[i], "--closed") || !strcmp(argv[i], "--snap") ||
            !strcmp(argv[i], "--modes-only") || !strcmp(argv[i], "--modeports")) {
            if (!strcmp(argv[i], "--mesh-only")) mesh_only = 1;
            else if (!strcmp(argv[i], "--closed")) closed = 1;
            else if (!strcmp(argv[i], "--snap")) snap = 1;
            else if (!strcmp(argv[i], "--modeports")) modeports = 1;
            else modes_only = 1;
            for (int k = i; k + 1 < argc; k++) argv[k] = argv[k + 1];
            argc--; i--;
        }
    for (int i = 1; i + 1 < argc; i += 2) {
        double v = atof(argv[i + 1]);
        if (!strcmp(argv[i], "--hmin")) hmin = v * 1e-6;
        else if (!strcmp(argv[i], "--hmax")) hmax = v * 1e-6;
        else if (!strcmp(argv[i], "--dx")) dxl = v * 1e-6;
        else if (!strcmp(argv[i], "--len")) len = v * 1e-3;
        else if (!strcmp(argv[i], "--er")) er = v;
        else if (!strcmp(argv[i], "--fmax")) fmax = v * 1e9;
        else if (!strcmp(argv[i], "--df")) df = v * 1e9;
        else if (!strcmp(argv[i], "--decay")) decay = v;
        else if (!strcmp(argv[i], "--mon")) mon_mm = v;
        else if (!strcmp(argv[i], "--feed")) feed_mm = v;
        else if (!strcmp(argv[i], "--yair")) yair_mm = v;
        else if (!strcmp(argv[i], "--zair")) zair_mm = v;
        else if (!strcmp(argv[i], "--cu")) cu_sigma = v;
        else if (!strcmp(argv[i], "--rq")) rq_um = v;
        else if (!strcmp(argv[i], "--tcu")) tcu_um = v;
        else if (!strcmp(argv[i], "--out")) out = argv[i + 1];
        else if (!strcmp(argv[i], "--debye")) debye = argv[i + 1];
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    const double h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6;
    const double a = w / 2, b = w / 2 + s, yf = b + fence;
    const double yext = yf + yair_mm * 1e-3, ztop = h + zair_mm * 1e-3, ext = 0.5e-3;
    const int npml = 8;

    /* ---- mesh ---- */
    /* y: symmetric about 0; CPML cells appended at the ends with hmax spacing */
    /* Metal edges in y. Default is the "thirds rule": a grid line hmin/3 inside the
     * metal and one 2*hmin/3 outside, which cancels most of the FDTD error from the
     * field singularity at zero-thickness edges (2.5% -> 0.2% in Z0 at 20 um cells,
     * see docs/design/0002-edge-meshing.md). --snap puts a line on the edge. */
    double t1 = hmin / 3, t2 = 2 * hmin / 3;
    double yfix_thirds[] = {-yext, -yf, -b - t1, -b + t2, -a - t2, -a + t1, 0, a - t1, a + t2, b - t2, b + t1, yf, yext};
    double yfix_snap[] = {-yext, -yf, -b, -a, 0, a, b, yf, yext};
    double yref[] = {-b, -a, a, b};
    if (tcu_um > 0) snap = 1; /* thick copper: grid lines on the copper faces */
    vec Y = snap ? graded(yfix_snap, 9, yref, 4, hmin, hmax, 1.3) : graded(yfix_thirds, 13, yref, 4, hmin, hmax, 1.3);
    const double tcu = tcu_um * 1e-6;
    double zfix[] = {0, h, ztop, h + tcu};
    double zref[] = {h, h + tcu};
    vec Z = tcu > 0 ? graded(zfix, 4, zref, 2, hmin, hmax, 1.3) : graded(zfix, 3, zref, 1, hmin, hmax, 1.3);
    /* append CPML cells */
    vec Y2 = {0};
    for (int i = npml; i >= 1; i--) push(&Y2, Y.v[0] - i * hmax);
    for (int i = 0; i < Y.n; i++) push(&Y2, Y.v[i]);
    for (int i = 1; i <= npml; i++) push(&Y2, Y.v[Y.n - 1] + i * hmax);
    for (int i = 1; i <= npml; i++) push(&Z, ztop + i * hmax);
    /* x: uniform dxl; ports at 0 and len; board extends ext past each, then CPML */
    vec X = {0};
    int nlen = (int)lround(len / dxl), next = (int)lround((modeports ? feed_mm * 1e-3 : ext) / dxl);
    for (int i = -next - npml; i <= nlen + next + npml; i++) push(&X, i * dxl);
    int nx = X.n - 1, ny = Y2.n - 1, nz = Z.n - 1;
    int ip1 = next + npml, ip2 = ip1 + nlen;
    int j0 = find_node(&Y2, 0.0), kh = find_node(&Z, h);
    int kt = tcu_um > 0 ? find_node(&Z, h + tcu_um * 1e-6) : kh; /* top of the copper */
    printf("mesh: %d x %d x %d = %.2f M cells (hmin %.1f um, dx %.1f um)\n", nx, ny, nz, nx * (double)ny * nz / 1e6,
           hmin * 1e6, dxl * 1e6);

    {
        mkdir(out, 0755);
        char mp[512];
        snprintf(mp, sizeof mp, "%s/mesh.json", out);
        FILE *fm = fopen(mp, "w");
        fprintf(fm, "{\n  \"snap\": %d, \"npml\": %d, \"closed\": %d, \"ip1\": %d, \"ip2\": %d, \"j0\": %d, \"kh\": %d,\n", snap, npml, closed, ip1, ip2, j0, kh);
        fprintf(fm, "  \"geometry\": {\"er\": %.17g, \"h\": %.17g, \"w\": %.17g, \"s\": %.17g, \"fence\": %.17g},\n", er, h, w, s, fence);
        const vec *ax[3] = {&X, &Y2, &Z};
        const char *nm[3] = {"x", "y", "z"};
        for (int q = 0; q < 3; q++) {
            fprintf(fm, "  \"%s\": [", nm[q]);
            for (int i = 0; i < ax[q]->n; i++) fprintf(fm, "%s%.17g", i ? ", " : "", ax[q]->v[i]);
            fprintf(fm, "]%s\n", q < 2 ? "," : "");
        }
        fprintf(fm, "}\n");
        fclose(fm);
        if (mesh_only) { printf("wrote %s\n", mp); return 0; }
    }
    cem_sim *S = cem_create(nx, ny, nz, X.v, Y2.v, Z.v);
    if (!S) { fprintf(stderr, "create: %s\n", cem_last_error()); return 1; }
    double xlo = X.v[0] - 1, xhi = X.v[nx] + 1, ylo = Y2.v[0] - 1, yhi = Y2.v[ny] + 1;
    cem_set_boundary(S, CEM_XLO, CEM_BC_CPML, npml);
    cem_set_boundary(S, CEM_XHI, CEM_BC_CPML, npml);
    if (!closed) {
        cem_set_boundary(S, CEM_YLO, CEM_BC_CPML, npml);
        cem_set_boundary(S, CEM_YHI, CEM_BC_CPML, npml);
        cem_set_boundary(S, CEM_ZHI, CEM_BC_CPML, npml);
    } else {
        /* shielded box: PEC on the y sides and top at the inner edge of where the
         * CPML would be (the layer cells become ordinary air), for comparison with
         * the 2D mode solver's closed cross-section */
        cem_pec_box(S, xlo, xhi, ylo, Y2.v[npml], ylo, yhi);
        cem_pec_box(S, xlo, xhi, Y2.v[ny - npml], yhi, ylo, yhi);
        cem_pec_box(S, xlo, xhi, ylo, yhi, Z.v[nz - npml], yhi);
    }
    /* ZLO stays PEC: the bottom ground plane */

    if (debye) {
        /* dispersive substrate: file from tools/dielectric.py (eps_inf, then "f_p deps" lines) */
        FILE *fd = fopen(debye, "r");
        if (!fd) { fprintf(stderr, "cannot open %s\n", debye); return 1; }
        double einf, fp[16], dep[16];
        int np = 0;
        if (fscanf(fd, "%lf", &einf) != 1) { fprintf(stderr, "bad debye file\n"); return 1; }
        while (np < 16 && fscanf(fd, "%lf %lf", &fp[np], &dep[np]) == 2) np++;
        fclose(fd);
        if (cem_set_debye_poles(S, np, fp) ||
            cem_material_box_debye(S, xlo, xhi, ylo, yhi, 0, h, einf, 0, dep)) {
            fprintf(stderr, "debye: %s\n", cem_last_error()); return 1;
        }
        printf("dispersive substrate: eps_inf %.4f, %d Debye poles\n", einf, np);
    } else {
        cem_material_box(S, xlo, xhi, ylo, yhi, 0, h, er, 0);
    }
    /* via-fence walls: solid metal outside |y| > yf below the top ground */
    cem_pec_box(S, xlo, xhi, ylo, -yf, 0, h);
    cem_pec_box(S, xlo, xhi, yf, yhi, 0, h);
    /* top grounds (sheet at z = h, or copper boxes of thickness tcu) */
    double tc = tcu_um * 1e-6;
    if (tc > 0) {
        cem_pec_box(S, xlo, xhi, ylo, -b, h, h + tc);
        cem_pec_box(S, xlo, xhi, b, yhi, h, h + tc);
    } else {
        cem_pec_sheet(S, CEM_Z, h, xlo, xhi, ylo, -b);
        cem_pec_sheet(S, CEM_Z, h, xlo, xhi, b, yhi);
    }
    /* strip between the port columns */
    if (tc > 0) {
        if (modeports) cem_pec_box(S, xlo, xhi, -a, a, h, h + tc);
        else cem_pec_box(S, X.v[ip1], X.v[ip2], -a, a, h, h + tc);
    } else if (modeports) cem_pec_sheet(S, CEM_Z, h, xlo, xhi, -a, a);
    else cem_pec_sheet(S, CEM_Z, h, X.v[ip1], X.v[ip2], -a, a);
    if (cu_sigma > 0) {
        if (cem_set_conductor_loss(S, cu_sigma, rq_um > 0 ? 1 : 0, rq_um * 1e-6, 0, df, fmax)) {
            fprintf(stderr, "conductor: %s\n", cem_last_error()); return 1;
        }
    }

    /* ports: z-directed columns from ground (k=0) to strip (kh), V = strip - ground */
    /* pulse: spectrum down ~ -60 dB at 1.25 fmax */
    double tau = sqrt(6.9) / (M_PI * 1.25 * fmax);
    int p1 = -1, p2 = -1, m1 = -1, m2 = -1;
    int jl = find_node(&Y2, -(a + b) / 2), jr = find_node(&Y2, (a + b) / 2);
    int ja = jl - 1 < 0 ? 0 : jl, jb = jr - 1; /* dual nodes ja+1/2, jb+1/2 */
    if (!modeports) {
        p1 = cem_add_lumped_port(S, CEM_Z, ip1, j0, 0, kh, -1, 50.0);
        p2 = cem_add_lumped_port(S, CEM_Z, ip2, j0, 0, kh, -1, 50.0);
        if (p1 < 0 || p2 < 0) { fprintf(stderr, "port: %s\n", cem_last_error()); return 1; }
        cem_set_gaussian(S, p1, 1.0, 4.5 * tau, tau);
    } else {
        /* source 2 cells in front of the CPML, monitor 1 mm further in, reference
         * planes at the ends of the length-L section (ip1, ip2) */
        int mon = (int)lround(mon_mm * 1e-3 / dxl);
        cem_mode_port_def d = {0};
        if (closed) { d.j0 = npml; d.j1 = ny - npml; d.k0 = 0; d.k1 = nz - npml; }
        else { d.j0 = 0; d.j1 = ny; d.k0 = 0; d.k1 = nz; } /* window includes the CPML */
        d.vj = j0; d.vk0 = 0; d.vk1 = kh; d.vsign = -1;
        d.ja = ja; d.jb = jb; d.ka = kh - 1; d.kb = kt;
        d.f0 = 0.5 * fmax; d.neff0 = 1.65;
        d.dir = 1; d.i_src = npml + 2; d.i_mon = d.i_src + mon; d.x_ref = X.v[ip1];
        m1 = cem_add_mode_port(S, &d);
        d.dir = -1; d.i_src = nx - npml - 2; d.i_mon = d.i_src - mon - 1; d.x_ref = X.v[ip2];
        m2 = cem_add_mode_port(S, &d);
        if (m1 < 0 || m2 < 0) { fprintf(stderr, "mode port: %s\n", cem_last_error()); return 1; }
        cem_set_mode_excitation(S, m1, 4.5 * tau, tau);
    }

    int nf = (int)lround(fmax / df);
    double *f = malloc(nf * sizeof(double));
    for (int i = 0; i < nf; i++) f[i] = (i + 1) * df;
    cem_set_frequencies(S, nf, f);

    /* probes in the middle of the line: V at nodes, I at cells. Gap centre for the
     * current loop: dual node nearest y = -(a+b)/2 and +(a+b)/2 */
    int nprobe_sets = 0;
    int vpi[64], ipi[64], vnode[64], icell[64];
    int nv = 0, ni = 0;
    /* V probes every 1 mm from 2 mm to 8 mm of a 10 mm line (scaled) and adjacent nodes */
    int step = (int)lround(1e-3 / dxl);
    for (int m = ip1 + 2 * step; m <= ip2 - 2 * step && nv < 60; m += step) {
        vnode[nv] = m;
        vpi[nv++] = cem_add_vprobe(S, CEM_Z, m, j0, 0, kh, -1);
        vnode[nv] = m + 1;
        vpi[nv++] = cem_add_vprobe(S, CEM_Z, m + 1, j0, 0, kh, -1);
        icell[ni] = m;
        ipi[ni++] = cem_add_iprobe(S, CEM_X, m, ja, jb, kh - 1, kt);
        nprobe_sets++;
    }
    for (int i = 0; i < nv; i++) if (vpi[i] < 0) { fprintf(stderr, "vprobe: %s\n", cem_last_error()); return 1; }
    for (int i = 0; i < ni; i++) if (ipi[i] < 0) { fprintf(stderr, "iprobe: %s\n", cem_last_error()); return 1; }

    if (cem_finalize(S, 0.99)) { fprintf(stderr, "finalize: %s\n", cem_last_error()); return 1; }
    if (cu_sigma > 0) printf("copper sigma %.3g S/m, roughness Rq %.2f um, Zs fit error %.3f%%\n", cu_sigma, rq_um,
                             100 * cem_conductor_fit_error(S));
    double dt = cem_dt(S);
    if (modes_only) {
        /* 2D eigenmode of the line cross-section (PEC window at the CPML edge) */
        int im = (ip1 + ip2) / 2;
        cem_mode2d *M = closed ? cem_cross_section_x(S, im, npml, ny - npml, 0, nz - npml)
                               : cem_cross_section_x(S, im, 0, ny, 0, nz); /* include the CPML */
        if (!M) { fprintf(stderr, "cross-section: %s\n", cem_last_error()); return 1; }
        double fl[] = {0.5e9, 1e9, 2e9, 5e9, 10e9, 20e9, 30e9, 40e9};
        /* indices inside the window */
        int off = closed ? npml : 0;
        int jw0 = j0 - off, jla = ja - off, jrb = jb - off;
        double neff = 1.63;
        for (int q = 0; q < 8; q++) {
            int rc = cem_mode2d_solve(M, fl[q], neff, NULL, 200, 1e-13);
            double beta = cem_mode2d_beta(M), k0 = 2 * M_PI * fl[q] / CEM_C0, V[2], I[2];
            cem_mode2d_vi(M, jw0, 0, kh, jla, jrb, kh - 1, kt, V, I);
            double P = cem_mode2d_power(M);
            neff = beta / k0;
            double vv = V[0] * V[0] + V[1] * V[1], ii = I[0] * I[0] + I[1] * I[1];
            double zr = (V[0] * I[0] + V[1] * I[1]) / ii, zi = (V[1] * I[0] - V[0] * I[1]) / ii;
            double Rs = sqrt(M_PI * fl[q] * CEM_MU0 / 5.8e7);
            printf("%6.2f GHz rc %d eps_eff %.6f alpha %.3e Zvi %.4f%+.4fj Zpv %.4f Zpi %.4f alpha_c(Cu) %.4f\n", fl[q] / 1e9, rc,
                   neff * neff, cem_mode2d_alpha(M), zr, zi, vv / (2 * fabs(P)), 2 * fabs(P) / ii,
                   cem_mode2d_conductor_alpha(M, Rs));
        }
        cem_mode2d_destroy(M);
        return 0;
    }
    printf("dt = %.4e s, pulse tau = %.2f ps\n", dt, tau * 1e12);
    clock_t c0 = clock();
    time_t w0 = time(NULL);
    long steps = cem_run(S, (long)(5e-9 / dt), decay, 500);
    double wall = difftime(time(NULL), w0);
    printf("ran %ld steps (%.3f ns) in %.0f s wall, %.1f Mcell-steps/s\n", steps, steps * dt * 1e9, wall,
           nx * (double)ny * nz * steps / 1e6 / (wall > 0 ? wall : 1));
    (void)c0;

    mkdir(out, 0755);
    char path[512];
    if (modeports) {
        double *A1 = malloc(2 * nf * sizeof(double)), *B1 = malloc(2 * nf * sizeof(double));
        double *A2 = malloc(2 * nf * sizeof(double)), *B2 = malloc(2 * nf * sizeof(double));
        double *Z1 = malloc(nf * sizeof(double)), *BT = malloc(nf * sizeof(double)), *AL = malloc(nf * sizeof(double));
        double *AC = malloc(nf * sizeof(double));
        if (cem_mode_port_waves(S, m1, A1, B1, Z1, BT) || cem_mode_port_waves(S, m2, A2, B2, NULL, NULL) ||
            cem_mode_port_alpha(S, m1, AL) || cem_mode_port_alpha_c(S, m1, AC)) {
            fprintf(stderr, "waves: %s\n", cem_last_error()); return 1;
        }
        snprintf(path, sizeof path, "%s/waves.csv", out);
        FILE *fw = fopen(path, "w");
        fprintf(fw, "f_hz,z0,beta,alpha,alpha_c,a1_re,a1_im,b1_re,b1_im,a2_re,a2_im,b2_re,b2_im\n");
        for (int i = 0; i < nf; i++)
            fprintf(fw, "%.6e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", f[i], Z1[i], BT[i], AL[i], AC[i], A1[2 * i],
                    A1[2 * i + 1], B1[2 * i], B1[2 * i + 1], A2[2 * i], A2[2 * i + 1], B2[2 * i], B2[2 * i + 1]);
        fclose(fw);
        printf("wrote %s\n", path);
        /* Touchstone, renormalised to 50 ohm; symmetric line assumed (S22 = S11, S12 = S21) */
        snprintf(path, sizeof path, "%s/line_50ohm.s2p", out);
        FILE *ft2 = fopen(path, "w");
        fprintf(ft2, "! Courant EM %s gcpw_line, mode ports, renormalised from Z0(f) to 50 ohm\n", cem_version());
        fprintf(ft2, "! port 1 excited; S22 = S11 and S12 = S21 by symmetry\n# HZ S RI R 50\n");
        for (int i = 0; i < nf; i++) {
            double ar = A1[2 * i], ai = A1[2 * i + 1], d = ar * ar + ai * ai;
            double s11r = (B1[2 * i] * ar + B1[2 * i + 1] * ai) / d, s11i = (B1[2 * i + 1] * ar - B1[2 * i] * ai) / d;
            double s21r = (B2[2 * i] * ar + B2[2 * i + 1] * ai) / d, s21i = (B2[2 * i + 1] * ar - B2[2 * i] * ai) / d;
            double Sm[8] = {s11r, s11i, s21r, s21i, s21r, s21i, s11r, s11i};
            cem_s2_renormalize(Sm, Z1[i], Z1[i], 50.0);
            fprintf(ft2, "%.6e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n", f[i], Sm[0], Sm[1], Sm[4], Sm[5], Sm[2], Sm[3],
                    Sm[6], Sm[7]);
        }
        fclose(ft2);
        printf("wrote %s\n", path);
    }
    const double *V1 = modeports ? NULL : cem_port_v(S, p1), *I1 = modeports ? NULL : cem_port_i(S, p1);
    const double *V2 = modeports ? NULL : cem_port_v(S, p2), *I2 = modeports ? NULL : cem_port_i(S, p2);
    snprintf(path, sizeof path, "%s/sparams.csv", out);
    FILE *fp = fopen(path, "w");
    snprintf(path, sizeof path, "%s/line.s2p", out);
    FILE *ft = fopen(path, "w");
    fprintf(fp, "f_hz,s11_re,s11_im,s21_re,s21_im\n");
    fprintf(ft, "! Courant EM %s gcpw_line, port 1 excited; S22=S11, S12=S21 assumed by symmetry\n# HZ S RI R 50\n", cem_version());
    double Zr = 50.0, sq = sqrt(Zr);
    for (int i = 0; i < nf && !modeports; i++) {
        double v1r = V1[2 * i], v1i = V1[2 * i + 1], i1r = I1[2 * i], i1i = I1[2 * i + 1];
        double v2r = V2[2 * i], v2i = V2[2 * i + 1], i2r = I2[2 * i], i2i = I2[2 * i + 1];
        double a1r = (v1r + Zr * i1r) / (2 * sq), a1i = (v1i + Zr * i1i) / (2 * sq);
        double b1r = (v1r - Zr * i1r) / (2 * sq), b1i = (v1i - Zr * i1i) / (2 * sq);
        double b2r = (v2r - Zr * i2r) / (2 * sq), b2i = (v2i - Zr * i2i) / (2 * sq);
        double den = a1r * a1r + a1i * a1i;
        double s11r = (b1r * a1r + b1i * a1i) / den, s11i = (b1i * a1r - b1r * a1i) / den;
        double s21r = (b2r * a1r + b2i * a1i) / den, s21i = (b2i * a1r - b2r * a1i) / den;
        fprintf(fp, "%.6e,%.9e,%.9e,%.9e,%.9e\n", f[i], s11r, s11i, s21r, s21i);
        fprintf(ft, "%.6e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n", f[i], s11r, s11i, s21r, s21i, s21r, s21i,
                s11r, s11i);
    }
    fclose(fp);
    fclose(ft);

    snprintf(path, sizeof path, "%s/probes.csv", out);
    fp = fopen(path, "w");
    fprintf(fp, "f_hz");
    for (int i = 0; i < nv; i++) fprintf(fp, ",V%d_re,V%d_im", vnode[i], vnode[i]);
    for (int i = 0; i < ni; i++) fprintf(fp, ",I%d_re,I%d_im", icell[i], icell[i]);
    fprintf(fp, "\n");
    for (int q = 0; q < nf; q++) {
        fprintf(fp, "%.6e", f[q]);
        for (int i = 0; i < nv; i++) { const double *X_ = cem_vprobe(S, vpi[i]); fprintf(fp, ",%.9e,%.9e", X_[2 * q], X_[2 * q + 1]); }
        for (int i = 0; i < ni; i++) { const double *X_ = cem_iprobe(S, ipi[i]); fprintf(fp, ",%.9e,%.9e", X_[2 * q], X_[2 * q + 1]); }
        fprintf(fp, "\n");
    }
    fclose(fp);

    snprintf(path, sizeof path, "%s/meta.json", out);
    fp = fopen(path, "w");
    fprintf(fp, "{\n  \"cells\": [%d, %d, %d],\n  \"hmin\": %g, \"hmax\": %g, \"dx\": %g, \"len\": %g, \"er\": %g,\n", nx, ny, nz,
            hmin, hmax, dxl, len, er);
    fprintf(fp, "  \"h\": %g, \"w\": %g, \"s\": %g, \"fence\": %g,\n", h, w, s, fence);
    fprintf(fp, "  \"dt\": %.6e, \"steps\": %ld, \"wall_s\": %.0f,\n", dt, steps, wall);
    fprintf(fp, "  \"x_nodes_for_probes\": {");
    for (int i = 0; i < nv; i++) fprintf(fp, "%s\"%d\": %.9e", i ? ", " : "", vnode[i], X.v[vnode[i]]);
    fprintf(fp, "},\n  \"x_port1\": %.9e, \"x_port2\": %.9e,\n", X.v[ip1], X.v[ip2]);
    fprintf(fp, "  \"iprobe_y\": [%.9e, %.9e]\n}\n", 0.5 * (Y2.v[ja] + Y2.v[ja + 1]), 0.5 * (Y2.v[jb] + Y2.v[jb + 1]));
    fclose(fp);
    printf("wrote %s/\n", out);
    cem_destroy(S);
    return 0;
}
