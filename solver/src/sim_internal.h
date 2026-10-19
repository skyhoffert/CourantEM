/* Internal definitions shared by the solver sources (not part of the API). */
#ifndef COURANT_SIM_INTERNAL_H
#define COURANT_SIM_INTERNAL_H

#include "courant.h"

#include <stddef.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXPORTS 16
#define MAXPROBES 256


typedef struct {
    int axis, dir, nedge;
    size_t *edge;    /* field indices of the port edges */
    double *dl;      /* edge lengths */
    double *cs;      /* source coefficient per edge */
    double R;
    double *V, *I;   /* DFT accumulators, 2*nf */
} port_t;

typedef struct {
    int axis, dir, n;
    size_t *edge;
    double *dl;
    double *acc;     /* 2*nf */
} vprobe_t;

typedef struct {
    int n;
    int *comp;       /* 3..5 */
    size_t *idx;
    double *w;       /* signed weights (dual lengths) */
    double *acc;
} iprobe_t;

#define MAXMPORTS 8
#define MAXPOLES 16

typedef struct {
    int dir;                    /* +1: DUT lies on the +x side of the port */
    int i_src, i_mon;           /* source node plane, monitor node plane (I at cell i_mon) */
    int j0, j1, k0, k1;         /* cross-section window (nodes) */
    int vj, vk0, vk1, vsign;    /* voltage path: vsign * sum Ez dz at (i, vj, vk0..vk1) */
    int ja, jb, ka, kb;         /* current loop (dual nodes), current in +x */
    double x_ref, f0, neff0;
    int vp0, vp1, ip;           /* probe indices */
    int nsrc;
    int *src_comp;
    size_t *src_idx;
    double *src_w;
    cem_mode2d *xs;             /* cross-section at the monitor plane */
    double *z0, *beta, *kappa;  /* per output frequency, filled lazily: z0 and kappa complex (2*nf) */
    double *alpha, *alpha_c;
    int solved;
} mport_t;

struct cem_sim {
    int n[3];
    double *node[3];
    double *d[3];   /* primary spacing, n[a] entries */
    double *dd[3];  /* dual spacing at nodes, n[a]+1 entries */
    size_t si, sj, N;

    double *E[3], *H[3];
    double *ca[3], *cb[3];
    uint8_t *edge_pec[3];

    /* cells: stored with the same (i,j,k) indexing; only i<nx, j<ny, k<nz used */
    double *eps_r, *sigma;
    uint8_t *cell_pec;

    int bc[6], npml[6];
    double *ide[3]; /* 1/(kappa*dual) at nodes */
    double *idh[3]; /* 1/(kappa*primary) at cells */
    double *be[3], *ce[3], *bh[3], *ch[3];
    double *psiE[3][2]; /* [component][which derivative]: Ex: (y,z) Ey: (z,x) Ez: (x,y) */
    double *pml_sn[3], *pml_kn[3], *pml_an[3]; /* CPML sigma/kappa/alpha at nodes */
    double *pml_sc[3], *pml_kc[3], *pml_ac[3]; /* ... and at cells */
    double *psiH[3][2];
    double vold[MAXPORTS];  /* port voltages before the E update (between step phases) */
    int has_pml;

    double dt;
    long step;
    int finalized;

    int nport;
    port_t port[MAXPORTS];
    int src_port;
    double src_amp, src_t0, src_tau;
    double *src_acc;

    int nf;
    double *freq;
    double *phE, *phH; /* current phasors e^{-j w t} at E and H times, 2*nf */
    double *rotE;      /* e^{-j w dt}, 2*nf */

    int nvp, nip;
    vprobe_t vp[MAXPROBES];
    iprobe_t ip[MAXPROBES];

    int nmport, msrc;
    mport_t mport[MAXMPORTS];
    double msrc_t0, msrc_tau;

    /* Debye dispersion: global pole set, per-cell weights, per-edge ADE state */
    int npole;
    double pole_f[MAXPOLES];
    double *cell_deps[MAXPOLES];
    int ndisp;
    int *d_comp;
    size_t *d_idx;
    double *d_beta, *d_J, *d_Eold, *d_cb; /* beta, J: npole*ndisp (pole-major) */
    double d_k[MAXPOLES];

    /* surface-impedance copper (sibc.c) */
    int sibc_on, sibc_rough, sibc_npole;
    double sibc_sigma, sibc_p1, sibc_p2, sibc_fit_err;
    double sibc_a[12], sibc_c[12], sibc_e[12];
    int nsb;
    int *sb_comp, *sb_nf, *sb_hc;
    size_t *sb_idx, *sb_hi;
    double *sb_sg, *sb_y;
};

void cem_set_err(const char *m);
int cem_modeports_finalize(cem_sim *s);
void cem_modeports_inject(cem_sim *s, double t);
void cem_modeports_free(cem_sim *s);
double cem_modeport_source(const cem_sim *s, double t);
double cem_mode_freq_hint(const cem_mode2d *m);
void cem_conductor_zs(const cem_sim *s, double f, double *re, double *im);
/* dispersion helpers (dispersive.c) */
int cem_edge_debye(const cem_sim *s, int c, int i, int j, int k, double *deps);
void cem_disp_save(cem_sim *s);
void cem_disp_correct(cem_sim *s);
void cem_disp_update(cem_sim *s);
void cem_disp_free(cem_sim *s);
int cem_sibc_finalize(cem_sim *s, const unsigned char *is_port);
void cem_sibc_apply(cem_sim *s);
void cem_sibc_free(cem_sim *s);
void cem_nnls(int m, int n, const double *A, const double *b, double *x);
void cem_mode2d_set_loss(cem_mode2d *m, const double *sig_u, const double *sig_v, const double *sig_w);

static inline size_t IDX(const cem_sim *s, int i, int j, int k) {
    return (size_t)i * s->si + (size_t)j * s->sj + (size_t)k;
}

/* length of the compact CPML auxiliary psiE (is_e=1) / psiH (is_e=0) [comp][w] */
size_t cem_psi_len(const cem_sim *s, int is_e, int comp, int w);

#endif
