/*
 * Courant EM - mode ports for x-directed transmission lines.
 *
 * A mode port launches the line's own discrete eigenmode (from the 2D solver on
 * the FDTD cross-section) with a soft E-field source on a plane in front of the
 * x CPML; the line continues into the CPML, which absorbs the backward-launched
 * half and anything reflected by the device. Waves are separated at a monitor
 * plane between source and device from the voltage and current there, using the
 * mode impedance Z0(f) and propagation constant from the 2D solver at each
 * output frequency, then shifted to the reference plane.
 *
 * Discrete consistency: the 2D solve uses the FDTD's own transverse grid; its
 * frequency is pre-warped to (2/dt) sin(w dt/2) and its eigenvalue kappa maps to
 * the 3D propagation constant by kappa = (2/dx) sin(beta dx/2). The voltage is
 * interpolated to the current's half plane with V = (V_m + V_m+1)/(2 cos(beta dx/2)).
 * Design notes: docs/design/0003-mode-ports.md
 */
#include "sim_internal.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int cem_add_mode_port(cem_sim *s, const cem_mode_port_def *d) {
    if (s->finalized) { cem_set_err("already finalized"); return -1; }
    if (s->nmport >= MAXMPORTS) { cem_set_err("too many mode ports"); return -1; }
    if (d->i_mon < 1 || d->i_mon + 1 > s->n[0] || d->i_src < 1 || d->i_src >= s->n[0] ||
        d->j0 < 0 || d->j1 > s->n[1] || d->k0 < 0 || d->k1 > s->n[2] || (d->dir != 1 && d->dir != -1)) {
        cem_set_err("bad mode port definition");
        return -1;
    }
    if ((d->dir == 1 && d->i_src >= d->i_mon) || (d->dir == -1 && d->i_src <= d->i_mon + 1)) {
        cem_set_err("mode port source must lie behind the monitor (away from the device)");
        return -1;
    }
    mport_t *p = &s->mport[s->nmport];
    memset(p, 0, sizeof *p);
    p->dir = d->dir; p->i_src = d->i_src; p->i_mon = d->i_mon;
    p->j0 = d->j0; p->j1 = d->j1; p->k0 = d->k0; p->k1 = d->k1;
    p->vj = d->vj; p->vk0 = d->vk0; p->vk1 = d->vk1; p->vsign = d->vsign;
    p->ja = d->ja; p->jb = d->jb; p->ka = d->ka; p->kb = d->kb;
    p->x_ref = d->x_ref; p->f0 = d->f0; p->neff0 = d->neff0 > 0 ? d->neff0 : 1.0;
    p->vp0 = cem_add_vprobe(s, CEM_Z, d->i_mon, d->vj, d->vk0, d->vk1, d->vsign);
    p->vp1 = cem_add_vprobe(s, CEM_Z, d->i_mon + 1, d->vj, d->vk0, d->vk1, d->vsign);
    p->ip = cem_add_iprobe(s, CEM_X, d->i_mon, d->ja, d->jb, d->ka, d->kb);
    if (p->vp0 < 0 || p->vp1 < 0 || p->ip < 0) return -1;
    return s->nmport++;
}

int cem_add_mode_port_args(cem_sim *s, int dir, int i_src, int i_mon, int j0, int j1, int k0, int k1, int vj,
                           int vk0, int vk1, int vsign, int ja, int jb, int ka, int kb, double x_ref, double f0,
                           double neff0) {
    cem_mode_port_def d = {dir, i_src, i_mon, j0, j1, k0, k1, vj, vk0, vk1, vsign, ja, jb, ka, kb, x_ref, f0, neff0};
    return cem_add_mode_port(s, &d);
}

int cem_set_mode_excitation(cem_sim *s, int port, double t0, double tau) {
    if (port < 0 || port >= s->nmport) { cem_set_err("bad mode port"); return -1; }
    s->msrc = port; s->msrc_t0 = t0; s->msrc_tau = tau;
    return 0;
}

static double warp(const cem_sim *s, double f) {
    /* frequency at which the continuous-time 2D problem matches the leapfrog scheme */
    return 2.0 / s->dt * sin(M_PI * f * s->dt) / (2 * M_PI);
}

/* Called at the end of cem_finalize: mode profile at f0 -> source weights. */
int cem_modeports_finalize(cem_sim *s) {
    for (int q = 0; q < s->nmport; q++) {
        mport_t *p = &s->mport[q];
        p->xs = cem_cross_section_x(s, p->i_mon, p->j0, p->j1, p->k0, p->k1);
        if (!p->xs) return -1;
        /* inverse iteration converges slowly when two modes straddle the shift
         * (e.g. a line mode and the substrate's surface wave): retry nearby */
        static const double tries[] = {1.0, 1.08, 0.93, 1.16, 0.87};
        int ok = 0;
        for (int t = 0; t < 5 && !ok; t++) {
            if (!cem_mode2d_solve(p->xs, warp(s, p->f0), p->neff0 * tries[t], NULL, 300, 1e-13)) { ok = 1; p->neff0 *= tries[t]; }
        }
        if (!ok) {
            cem_set_err("mode port: eigenmode did not converge");
            return -1;
        }
        int ny = p->j1 - p->j0, nz = p->k1 - p->k0;
        const double *eu = cem_mode2d_eu(p->xs), *ev = cem_mode2d_ev(p->xs);
        /* fields are complex (interleaved); the mode is phase-normalised so the
         * real part carries the profile for a lossless line */
        double vmax = 0;
        for (int i = 0; i < ny * (nz + 1); i++) if (fabs(eu[2 * i]) > vmax) vmax = fabs(eu[2 * i]);
        for (int i = 0; i < (ny + 1) * nz; i++) if (fabs(ev[2 * i]) > vmax) vmax = fabs(ev[2 * i]);
        int nmax = ny * (nz + 1) + (ny + 1) * nz, n = 0;
        p->src_comp = malloc(nmax * sizeof(int));
        p->src_idx = malloc(nmax * sizeof(size_t));
        p->src_w = malloc(nmax * sizeof(double));
        for (int j = 0; j < ny; j++)
            for (int k = 0; k <= nz; k++) {
                double w = eu[2 * (j * (nz + 1) + k)] / vmax;
                if (w == 0) continue;
                size_t idx = IDX(s, p->i_src, p->j0 + j, p->k0 + k);
                if (s->cb[1][idx] == 0) continue;
                p->src_comp[n] = 1; p->src_idx[n] = idx; p->src_w[n] = w; n++;
            }
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k < nz; k++) {
                double w = ev[2 * (j * nz + k)] / vmax;
                if (w == 0) continue;
                size_t idx = IDX(s, p->i_src, p->j0 + j, p->k0 + k);
                if (s->cb[2][idx] == 0) continue;
                p->src_comp[n] = 2; p->src_idx[n] = idx; p->src_w[n] = w; n++;
            }
        p->nsrc = n;
        /* orient the profile so the source voltage is positive */
        double V[2], Iv[2];
        cem_mode2d_vi(p->xs, p->vj - p->j0, p->vk0 - p->k0, p->vk1 - p->k0, p->ja - p->j0, p->jb - p->j0,
                      p->ka - p->k0, p->kb - p->k0, V, Iv);
        if (-V[0] * p->vsign < 0)
            for (int i = 0; i < n; i++) p->src_w[i] = -p->src_w[i];
    }
    return 0;
}

double cem_modeport_source(const cem_sim *s, double t) {
    if (s->msrc < 0) return 0;
    double x = (t - s->msrc_t0) / s->msrc_tau;
    return -2 * x * exp(-x * x); /* tau * d/dt of exp(-x^2): zero DC, peak ~0.86 */
}

void cem_modeports_inject(cem_sim *s, double t) {
    if (s->msrc < 0) return;
    mport_t *p = &s->mport[s->msrc];
    double v = cem_modeport_source(s, t);
    for (int i = 0; i < p->nsrc; i++) s->E[p->src_comp[i]][p->src_idx[i]] += p->src_w[i] * v;
}

void cem_modeports_free(cem_sim *s) {
    for (int q = 0; q < s->nmport; q++) {
        mport_t *p = &s->mport[q];
        free(p->src_comp); free(p->src_idx); free(p->src_w);
        free(p->z0); free(p->beta); free(p->kappa); free(p->alpha); free(p->alpha_c);
        cem_mode2d_destroy(p->xs);
    }
}

static int solve_port_modes(cem_sim *s, mport_t *p) {
    int nf = s->nf;
    p->z0 = malloc(2 * nf * sizeof(double));
    p->beta = malloc(nf * sizeof(double));
    p->alpha = malloc(nf * sizeof(double));
    p->kappa = malloc(2 * nf * sizeof(double));
    p->alpha_c = calloc(nf, sizeof(double));
    int ny = p->j1 - p->j0, nz = p->k1 - p->k0;
    int nst = ny * (nz + 1) + (ny + 1) * nz;
    double *start = malloc(2 * nst * sizeof(double));
    double dx = s->d[0][p->i_mon];
    /* walk outward from the frequency nearest f0, seeding each solve with the
     * previous mode so the same mode is tracked */
    int q0 = 0;
    for (int q = 0; q < nf; q++) if (fabs(s->freq[q] - p->f0) < fabs(s->freq[q0] - p->f0)) q0 = q;
    for (int pass = 0; pass < 2; pass++) {
        cem_mode2d_solve(p->xs, warp(s, s->freq[q0]), p->neff0, NULL, 300, 1e-13);
        int step = pass == 0 ? 1 : -1;
        for (int q = q0; q >= 0 && q < nf; q += step) {
            double fw = warp(s, s->freq[q]);
            double ng = cem_mode2d_beta(p->xs) / (2 * M_PI * cem_mode_freq_hint(p->xs) / CEM_C0);
            memcpy(start, cem_mode2d_eu(p->xs), 2 * ny * (nz + 1) * sizeof(double));
            memcpy(start + 2 * ny * (nz + 1), cem_mode2d_ev(p->xs), 2 * (ny + 1) * nz * sizeof(double));
            if (cem_mode2d_solve(p->xs, fw, ng, start, 300, 1e-13)) {
                cem_set_err("mode port: eigenmode did not converge at an output frequency");
                free(start);
                return -1;
            }
            /* kappa = (2/dx) sin(gamma_3d dx / (2j)) ... in terms of the complex
             * propagation constant g = alpha + j beta: kappa_c = beta2d - j alpha2d,
             * beta3d_c = (2/dx) asin(kappa_c dx / 2) */
            double complex kc = cem_mode2d_beta(p->xs) - I * cem_mode2d_alpha(p->xs);
            double complex b3 = 2 / dx * casin(kc * dx / 2);
            p->kappa[2 * q] = creal(kc); p->kappa[2 * q + 1] = cimag(kc);
            p->beta[q] = creal(b3);
            p->alpha[q] = -cimag(b3);
            if (s->sibc_on) {
                /* copper loss by perturbation of the lossless-conductor mode, with
                 * Rs = Re Zs of the fitted surface-impedance model */
                double zr, zi;
                cem_conductor_zs(s, s->freq[q], &zr, &zi);
                p->alpha_c[q] = cem_mode2d_conductor_alpha(p->xs, zr);
                p->alpha[q] += p->alpha_c[q];
            }
            double V[2], Iv[2];
            cem_mode2d_vi(p->xs, p->vj - p->j0, p->vk0 - p->k0, p->vk1 - p->k0, p->ja - p->j0, p->jb - p->j0,
                          p->ka - p->k0, p->kb - p->k0, V, Iv);
            /* cem_mode2d_vi returns -sum Ez dz; the port's voltage is vsign * sum Ez dz */
            double complex Z = -p->vsign * (V[0] + I * V[1]) / (Iv[0] + I * Iv[1]);
            p->z0[2 * q] = creal(Z); p->z0[2 * q + 1] = cimag(Z);
        }
    }
    free(start);
    p->solved = 1;
    return 0;
}

int cem_mode_port_waves(cem_sim *s, int port, double *a, double *b, double *z0, double *beta) {
    if (port < 0 || port >= s->nmport || !s->finalized) { cem_set_err("bad mode port"); return -1; }
    mport_t *p = &s->mport[port];
    if (!p->solved && solve_port_modes(s, p)) return -1;
    const double *V0 = cem_vprobe(s, p->vp0), *V1 = cem_vprobe(s, p->vp1), *Ix = cem_iprobe(s, p->ip);
    double dx = s->d[0][p->i_mon];
    double xmon = 0.5 * (s->node[0][p->i_mon] + s->node[0][p->i_mon + 1]);
    double d = p->dir * (p->x_ref - xmon); /* distance from monitor to reference, into the device */
    for (int q = 0; q < s->nf; q++) {
        double complex Z = p->z0[2 * q] + I * p->z0[2 * q + 1];
        double complex g = p->alpha[q] + I * p->beta[q];      /* gamma */
        double complex c = ccosh(g * dx / 2);                  /* = cos(beta dx/2) when lossless */
        double complex V = ((V0[2 * q] + V1[2 * q]) + I * (V0[2 * q + 1] + V1[2 * q + 1])) / (2 * c);
        double complex Iin = p->dir * (Ix[2 * q] + I * Ix[2 * q + 1]); /* current flowing into the device */
        double complex sq = 2 * csqrt(Z);
        double complex A = (V + Z * Iin) / sq, Bw = (V - Z * Iin) / sq;
        A *= cexp(-g * d);   /* inward wave travels to the reference plane */
        Bw *= cexp(g * d);   /* outward wave came from it */
        a[2 * q] = creal(A); a[2 * q + 1] = cimag(A);
        b[2 * q] = creal(Bw); b[2 * q + 1] = cimag(Bw);
        if (z0) z0[q] = creal(Z);
        if (beta) beta[q] = p->beta[q];
    }
    return 0;
}

int cem_mode_port_alpha(const cem_sim *s, int port, double *alpha) {
    if (port < 0 || port >= s->nmport || !s->mport[port].solved) { cem_set_err("mode port not solved"); return -1; }
    memcpy(alpha, s->mport[port].alpha, s->nf * sizeof(double));
    return 0;
}

int cem_mode_port_alpha_c(const cem_sim *s, int port, double *alpha_c) {
    if (port < 0 || port >= s->nmport || !s->mport[port].solved) { cem_set_err("mode port not solved"); return -1; }
    memcpy(alpha_c, s->mport[port].alpha_c, s->nf * sizeof(double));
    return 0;
}

/* Z = sqrt(Zd) (I + S) (I - S)^-1 sqrt(Zd); S' = (Z - zref I)(Z + zref I)^-1 */
void cem_s2_renormalize(double *Sv, double z1, double z2, double zref) {
    double complex S[2][2] = {{Sv[0] + I * Sv[1], Sv[2] + I * Sv[3]}, {Sv[4] + I * Sv[5], Sv[6] + I * Sv[7]}};
    double complex A[2][2] = {{1 - S[0][0], -S[0][1]}, {-S[1][0], 1 - S[1][1]}};
    double complex det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    double complex Ai[2][2] = {{A[1][1] / det, -A[0][1] / det}, {-A[1][0] / det, A[0][0] / det}};
    double complex B[2][2] = {{1 + S[0][0], S[0][1]}, {S[1][0], 1 + S[1][1]}};
    double complex Zn[2][2];
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) Zn[r][c] = B[r][0] * Ai[0][c] + B[r][1] * Ai[1][c];
    double sq[2] = {sqrt(z1), sqrt(z2)};
    double complex Z[2][2];
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) Z[r][c] = sq[r] * Zn[r][c] * sq[c];
    double complex P[2][2] = {{Z[0][0] + zref, Z[0][1]}, {Z[1][0], Z[1][1] + zref}};
    double complex M[2][2] = {{Z[0][0] - zref, Z[0][1]}, {Z[1][0], Z[1][1] - zref}};
    double complex dp = P[0][0] * P[1][1] - P[0][1] * P[1][0];
    double complex Pi[2][2] = {{P[1][1] / dp, -P[0][1] / dp}, {-P[1][0] / dp, P[0][0] / dp}};
    double complex R[2][2];
    for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++) R[r][c] = M[r][0] * Pi[0][c] + M[r][1] * Pi[1][c];
    Sv[0] = creal(R[0][0]); Sv[1] = cimag(R[0][0]); Sv[2] = creal(R[0][1]); Sv[3] = cimag(R[0][1]);
    Sv[4] = creal(R[1][0]); Sv[5] = cimag(R[1][0]); Sv[6] = creal(R[1][1]); Sv[7] = cimag(R[1][1]);
}
