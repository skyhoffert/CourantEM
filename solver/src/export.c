/*
 * Courant EM - access to the finalized update plan for external time steppers
 * (the WebGPU backend). The C core still builds the geometry, coefficients,
 * ports and probes; another engine runs the time loop and hands the probe
 * samples back with cem_push_samples(), so every post-processing path (DFTs,
 * mode-port waves, renormalisation) is shared with the reference solver.
 */
#include "sim_internal.h"

#include <math.h>
#include <string.h>

/* which: 0 ca, 1 cb, 2 ide (inverse dual spacing incl. kappa, n+1),
 * 3 idh (inverse primary spacing incl. kappa, n), 4 be, 5 ce (n+1), 6 bh, 7 ch (n),
 * 8 primary spacing d (n), 9 dual spacing dd (n+1) */
const double *cem_plan_array(const cem_sim *s, int which, int comp) {
    if (!s->finalized || comp < 0 || comp > 2) return NULL;
    switch (which) {
    case 0: return s->ca[comp];
    case 1: return s->cb[comp];
    case 2: return s->ide[comp];
    case 3: return s->idh[comp];
    case 4: return s->be[comp];
    case 5: return s->ce[comp];
    case 6: return s->bh[comp];
    case 7: return s->ch[comp];
    case 8: return s->d[comp];
    case 9: return s->dd[comp];
    default: return NULL;
    }
}

int cem_plan_npml(const cem_sim *s, int face) { return face >= 0 && face < 6 ? s->npml[face] : 0; }
int cem_plan_has_pml(const cem_sim *s) { return s->has_pml; }
int cem_plan_nport(const cem_sim *s) { return s->nport; }
int cem_plan_nvprobe(const cem_sim *s) { return s->nvp; }
int cem_plan_niprobe(const cem_sim *s) { return s->nip; }
int cem_plan_ndisp(const cem_sim *s) { return s->ndisp; }
int cem_plan_nsibc(const cem_sim *s) { return s->nsb; }

/* Lumped port p: number of edges; fills comp (axis), field index, dl*dir weights
 * for the voltage sum, and the per-unit-Vs source coefficient cs. */
int cem_plan_port(const cem_sim *s, int p, double *idx, double *vw, double *cs) {
    if (p < 0 || p >= s->nport) return -1;
    const port_t *pt = &s->port[p];
    for (int e = 0; e < pt->nedge; e++) {
        if (idx) idx[e] = (double)(pt->axis * s->N + pt->edge[e]);
        if (vw) vw[e] = pt->dl[e] * pt->dir;
        if (cs) cs[e] = pt->cs[e];
    }
    return pt->nedge;
}

/* Voltage probe p: global E index (comp*N + idx) and weights (signed dl). */
int cem_plan_vprobe(const cem_sim *s, int p, double *idx, double *w) {
    if (p < 0 || p >= s->nvp) return -1;
    const vprobe_t *v = &s->vp[p];
    for (int e = 0; e < v->n; e++) {
        if (idx) idx[e] = (double)(v->axis * s->N + v->edge[e]);
        if (w) w[e] = v->dl[e];
    }
    return v->n;
}

/* Current probe p: global H index ((comp-3)*N + idx) and weights. */
int cem_plan_iprobe(const cem_sim *s, int p, double *idx, double *w) {
    if (p < 0 || p >= s->nip) return -1;
    const iprobe_t *c = &s->ip[p];
    for (int e = 0; e < c->n; e++) {
        if (idx) idx[e] = (double)((c->comp[e] - 3) * s->N + c->idx[e]);
        if (w) w[e] = c->w[e];
    }
    return c->n;
}

/* Excited mode port's source entries: global E index and weight. Returns count (0 if none). */
int cem_plan_mode_source(const cem_sim *s, double *idx, double *w) {
    if (s->msrc < 0) return 0;
    const mport_t *p = &s->mport[s->msrc];
    for (int i = 0; i < p->nsrc; i++) {
        if (idx) idx[i] = (double)(p->src_comp[i] * s->N + p->src_idx[i]);
        if (w) w[i] = p->src_w[i];
    }
    return p->nsrc;
}

/* Source waveforms at H times t = (n + 1/2) dt for steps n0 .. n0+n-1 */
double cem_modeport_source(const cem_sim *s, double t);
void cem_plan_waveforms(const cem_sim *s, long n0, int n, double *lumped, double *mode) {
    for (int q = 0; q < n; q++) {
        double t = (n0 + q + 0.5) * s->dt;
        if (lumped) {
            double v = 0;
            if (s->src_port >= 0) { double x = (t - s->src_t0) / s->src_tau; v = s->src_amp * exp(-x * x); }
            lumped[q] = v;
        }
        if (mode) mode[q] = cem_modeport_source(s, t);
    }
}

int cem_plan_excited_lumped(const cem_sim *s) { return s->src_port; }

static void set_phasors_ext(cem_sim *s, double tE, double tH) {
    for (int f = 0; f < s->nf; f++) {
        double w = 2 * M_PI * s->freq[f];
        s->phE[2 * f] = cos(w * tE); s->phE[2 * f + 1] = -sin(w * tE);
        s->phH[2 * f] = cos(w * tH); s->phH[2 * f + 1] = -sin(w * tH);
    }
}

static inline void acc(double *a, const double *ph, int nf, double v) {
    for (int f = 0; f < nf; f++) { a[2 * f] += v * ph[2 * f]; a[2 * f + 1] += v * ph[2 * f + 1]; }
}

/* Feed n steps of samples recorded by an external engine (layout: step-major).
 *   ip:   n * nip current-probe values (H time)
 *   vold: n * nport port voltages before the E update (already times dir)
 *   vnew: n * nport port voltages after the E update incl. source (times dir)
 *   vp:   n * nvp voltage-probe values (E time)
 * Accumulates exactly the DFTs cem_step would, and advances the step counter. */
void cem_push_samples(cem_sim *s, int n, const double *ip, const double *vold, const double *vnew, const double *vp) {
    int nf = s->nf;
    for (int q = 0; q < n; q++) {
        long st = s->step;
        double tH = (st + 0.5) * s->dt, tE = (st + 1) * s->dt;
        set_phasors_ext(s, tE, tH);
        for (int p = 0; p < s->nip; p++) acc(s->ip[p].acc, s->phH, nf, ip[(size_t)q * s->nip + p] * s->dt);
        double Vs = 0;
        if (s->src_port >= 0) {
            double x = (tH - s->src_t0) / s->src_tau;
            Vs = s->src_amp * exp(-x * x);
            acc(s->src_acc, s->phH, nf, Vs * s->dt);
        }
        if (s->msrc >= 0) acc(s->src_acc, s->phH, nf, cem_modeport_source(s, tH) * s->dt);
        for (int p = 0; p < s->nport; p++) {
            port_t *pt = &s->port[p];
            double V = 0.5 * (vnew[(size_t)q * s->nport + p] + vold[(size_t)q * s->nport + p]);
            double vs = p == s->src_port ? Vs : 0.0;
            double I = (vs - V) / pt->R;
            acc(pt->V, s->phH, nf, V * s->dt);
            acc(pt->I, s->phH, nf, I * s->dt);
        }
        for (int p = 0; p < s->nvp; p++) acc(s->vp[p].acc, s->phE, nf, vp[(size_t)q * s->nvp + p] * s->dt);
        s->step++;
    }
    /* leave the phasors where cem_step expects them, so a checkpoint taken now
     * (or CPU steps after this) continue the DFTs correctly */
    set_phasors_ext(s, (s->step + 1) * s->dt, (s->step + 0.5) * s->dt);
}

/* Dispersion: number of poles; per-edge global E index, cb, and beta (pole-major, npole*ndisp); k per pole */
int cem_plan_npole(const cem_sim *s) { return s->npole; }
int cem_plan_disp(const cem_sim *s, double *idx, double *cb, double *beta, double *k) {
    for (int e = 0; e < s->ndisp; e++) {
        if (idx) idx[e] = (double)(s->d_comp[e] * s->N + s->d_idx[e]);
        if (cb) cb[e] = s->d_cb[e];
    }
    if (beta) memcpy(beta, s->d_beta, (size_t)s->npole * s->ndisp * sizeof(double));
    if (k) memcpy(k, s->d_k, s->npole * sizeof(double));
    return s->ndisp;
}

/* Surface impedance: per edge global E index and up to 4 (global H index, weight) faces;
 * per term c_k and e_k = exp(-a_k dt). Returns number of edges; *nterm gets the term count. */
int cem_plan_sibc(const cem_sim *s, double *idx, double *hidx, double *hw, double *c, double *e, int *nterm) {
    for (int q = 0; q < s->nsb; q++) {
        if (idx) idx[q] = (double)(s->sb_comp[q] * s->N + s->sb_idx[q]);
        for (int f = 0; f < 4; f++) {
            int has = f < s->sb_nf[q];
            if (hidx) hidx[4 * q + f] = has ? (double)(s->sb_hc[4 * q + f] * s->N + s->sb_hi[4 * q + f]) : 0;
            if (hw) hw[4 * q + f] = has ? s->sb_sg[4 * q + f] : 0;
        }
    }
    if (c) memcpy(c, s->sibc_c, s->sibc_npole * sizeof(double));
    if (e) memcpy(e, s->sibc_e, s->sibc_npole * sizeof(double));
    if (nterm) *nterm = s->sibc_npole;
    return s->nsb;
}

/* ---- checkpoints ------------------------------------------------------------
 * Everything that evolves in time: fields, CPML auxiliaries, Debye currents,
 * surface-impedance states, DFT accumulators and the step counter. The sim must
 * be rebuilt identically (same geometry, ports, probes, frequencies) and
 * finalized before a load; the header guards the sizes. Layout is doubles. */
#define CKPT_MAGIC 20261017.0

static size_t ckpt_count(const cem_sim *s, int write, double *buf, const double *in) {
    size_t o = 0;
#define BLOCK(ptr, n)                                                            \
    do {                                                                         \
        size_t n_ = (size_t)(n);                                                 \
        if (buf && (ptr)) memcpy(buf + o, (ptr), n_ * sizeof(double));           \
        if (in && (ptr)) memcpy((ptr), in + o, n_ * sizeof(double));             \
        o += n_;                                                                 \
    } while (0)
    (void)write;
    for (int c = 0; c < 3; c++) BLOCK(s->E[c], s->N);
    for (int c = 0; c < 3; c++) BLOCK(s->H[c], s->N);
    if (s->has_pml)
        for (int c = 0; c < 3; c++)
            for (int w = 0; w < 2; w++) { BLOCK(s->psiE[c][w], cem_psi_len(s, 1, c, w)); BLOCK(s->psiH[c][w], cem_psi_len(s, 0, c, w)); }
    if (s->ndisp) { BLOCK(s->d_J, (size_t)s->npole * s->ndisp); BLOCK(s->d_Eold, s->ndisp); }
    if (s->nsb) BLOCK(s->sb_y, (size_t)s->nsb * s->sibc_npole);
    for (int p = 0; p < s->nport; p++) { BLOCK(s->port[p].V, 2 * s->nf); BLOCK(s->port[p].I, 2 * s->nf); }
    for (int p = 0; p < s->nvp; p++) BLOCK(s->vp[p].acc, 2 * s->nf);
    for (int p = 0; p < s->nip; p++) BLOCK(s->ip[p].acc, 2 * s->nf);
    BLOCK(s->src_acc, 2 * s->nf);
    /* the DFT phasors are advanced by rotation, so they are state too */
    BLOCK(s->phE, 2 * s->nf);
    BLOCK(s->phH, 2 * s->nf);
#undef BLOCK
    return o;
}

size_t cem_checkpoint_size(const cem_sim *s) {
    if (!s->finalized) return 0;
    return 8 + ckpt_count(s, 0, NULL, NULL);
}

int cem_checkpoint_save(const cem_sim *s, double *buf) {
    if (!s->finalized) { cem_set_err("not finalized"); return -1; }
    size_t body = ckpt_count(s, 0, NULL, NULL);
    double hdr[8] = {CKPT_MAGIC, 1, (double)s->N, (double)s->nf, (double)s->step, (double)body,
                     (double)(s->nport + 100 * s->nvp + 10000 * s->nip), (double)(s->ndisp + s->nsb)};
    memcpy(buf, hdr, sizeof hdr);
    ckpt_count(s, 1, buf + 8, NULL);
    return 0;
}

int cem_checkpoint_load(cem_sim *s, const double *buf) {
    if (!s->finalized) { cem_set_err("not finalized"); return -1; }
    size_t body = ckpt_count(s, 0, NULL, NULL);
    if (buf[0] != CKPT_MAGIC || buf[1] != 1 || buf[2] != (double)s->N || buf[3] != (double)s->nf ||
        buf[5] != (double)body || buf[6] != (double)(s->nport + 100 * s->nvp + 10000 * s->nip) ||
        buf[7] != (double)(s->ndisp + s->nsb)) {
        cem_set_err("checkpoint does not match this simulation");
        return -1;
    }
    ckpt_count(s, 0, NULL, buf + 8);
    s->step = (long)buf[4];
    return 0;
}

/* Writable state arrays for engines that keep their own copy (GPU): which
 * 0: psiE[comp][w], 1: psiH[comp][w] (compact CPML layout; lengths from cem_state_len), 2: Debye J (npole*ndisp),
 * 3: Debye E_old (ndisp), 4: SIBC states (nsb*nterm). NULL if absent. */
double *cem_state_array(cem_sim *s, int which, int comp, int w) {
    switch (which) {
    case 0: return s->has_pml && comp >= 0 && comp < 3 && w >= 0 && w < 2 ? s->psiE[comp][w] : NULL;
    case 1: return s->has_pml && comp >= 0 && comp < 3 && w >= 0 && w < 2 ? s->psiH[comp][w] : NULL;
    case 2: return s->d_J;
    case 3: return s->d_Eold;
    case 4: return s->sb_y;
    default: return NULL;
    }
}

/* Length (doubles) of a state array from cem_state_array; 0 if absent. */
size_t cem_state_len(const cem_sim *s, int which, int comp, int w) {
    switch (which) {
    case 0: return cem_psi_len(s, 1, comp, w);
    case 1: return cem_psi_len(s, 0, comp, w);
    case 2: return s->d_J ? (size_t)s->npole * s->ndisp : 0;
    case 3: return s->d_Eold ? (size_t)s->ndisp : 0;
    case 4: return s->sb_y ? (size_t)s->nsb * s->sibc_npole : 0;
    default: return 0;
    }
}

/* Set the step counter (an external engine that ran k steps without pushing samples). */
void cem_set_step(cem_sim *s, long step) { s->step = step; }
