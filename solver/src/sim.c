/*
 * Courant EM - FDTD core, double precision reference implementation.
 * See include/courant.h for the API and the Yee layout.
 * Design notes: docs/design/0001-fdtd-kernel.md
 */
#include "sim_internal.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static char g_err[256] = "";
static void set_err(const char *m) {
    /* no stdio in the core: keeps the wasm module free of WASI imports */
    size_t i = 0;
    for (; m[i] && i + 1 < sizeof g_err; i++) g_err[i] = m[i];
    g_err[i] = 0;
}
void cem_set_err(const char *m) { set_err(m); }
const char *cem_last_error(void) { return g_err; }

size_t cem_index(const cem_sim *s, int i, int j, int k) { return IDX(s, i, j, k); }
void cem_dims(const cem_sim *s, int *nx, int *ny, int *nz) {
    *nx = s->n[0]; *ny = s->n[1]; *nz = s->n[2];
}
double *cem_field(cem_sim *s, int comp) { return comp < 3 ? s->E[comp] : s->H[comp - 3]; }
double cem_dt(const cem_sim *s) { return s->dt; }
long cem_step_count(const cem_sim *s) { return s->step; }
int cem_nfreq(const cem_sim *s) { return s->nf; }

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n ? n : 1, sz);
    return p;
}

cem_sim *cem_create(int nx, int ny, int nz, const double *x, const double *y, const double *z) {
    if (nx < 2 || ny < 2 || nz < 2) { set_err("grid needs at least 2 cells per axis"); return NULL; }
    cem_sim *s = xcalloc(1, sizeof *s);
    if (!s) { set_err("out of memory"); return NULL; }
    s->n[0] = nx; s->n[1] = ny; s->n[2] = nz;
    const double *src[3] = {x, y, z};
    for (int a = 0; a < 3; a++) {
        int n = s->n[a];
        s->node[a] = xcalloc(n + 1, sizeof(double));
        s->d[a] = xcalloc(n, sizeof(double));
        s->dd[a] = xcalloc(n + 1, sizeof(double));
        memcpy(s->node[a], src[a], (n + 1) * sizeof(double));
        for (int i = 0; i < n; i++) {
            s->d[a][i] = src[a][i + 1] - src[a][i];
            if (!(s->d[a][i] > 0)) { set_err("node coordinates must be strictly increasing"); cem_destroy(s); return NULL; }
        }
        for (int i = 0; i <= n; i++) {
            double lo = i > 0 ? s->d[a][i - 1] : 0.0, hi = i < n ? s->d[a][i] : 0.0;
            s->dd[a][i] = 0.5 * (lo + hi);
        }
    }
    s->sj = (size_t)nz + 1;
    s->si = s->sj * ((size_t)ny + 1);
    s->N = s->si * ((size_t)nx + 1);
    for (int c = 0; c < 3; c++) {
        s->E[c] = xcalloc(s->N, sizeof(double));
        s->H[c] = xcalloc(s->N, sizeof(double));
        s->ca[c] = xcalloc(s->N, sizeof(double));
        s->cb[c] = xcalloc(s->N, sizeof(double));
        s->edge_pec[c] = xcalloc(s->N, 1);
        if (!s->E[c] || !s->H[c] || !s->ca[c] || !s->cb[c] || !s->edge_pec[c]) {
            set_err("out of memory"); cem_destroy(s); return NULL;
        }
    }
    s->eps_r = xcalloc(s->N, sizeof(double));
    s->sigma = xcalloc(s->N, sizeof(double));
    s->cell_pec = xcalloc(s->N, 1);
    if (!s->eps_r || !s->sigma || !s->cell_pec) { set_err("out of memory"); cem_destroy(s); return NULL; }
    for (size_t q = 0; q < s->N; q++) s->eps_r[q] = 1.0;
    for (int f = 0; f < 6; f++) { s->bc[f] = CEM_BC_PEC; s->npml[f] = 0; }
    s->src_port = -1;
    s->msrc = -1;
    return s;
}

void cem_destroy(cem_sim *s) {
    if (!s) return;
    for (int a = 0; a < 3; a++) {
        free(s->node[a]); free(s->d[a]); free(s->dd[a]);
        free(s->E[a]); free(s->H[a]); free(s->ca[a]); free(s->cb[a]); free(s->edge_pec[a]);
        free(s->ide[a]); free(s->idh[a]); free(s->be[a]); free(s->ce[a]); free(s->bh[a]); free(s->ch[a]);
        free(s->pml_sn[a]); free(s->pml_kn[a]); free(s->pml_an[a]);
        free(s->pml_sc[a]); free(s->pml_kc[a]); free(s->pml_ac[a]);
        for (int w = 0; w < 2; w++) { free(s->psiE[a][w]); free(s->psiH[a][w]); }
    }
    free(s->eps_r); free(s->sigma); free(s->cell_pec);
    for (int p = 0; p < s->nport; p++) {
        free(s->port[p].edge); free(s->port[p].dl); free(s->port[p].cs);
        free(s->port[p].V); free(s->port[p].I);
    }
    for (int p = 0; p < s->nvp; p++) { free(s->vp[p].edge); free(s->vp[p].dl); free(s->vp[p].acc); }
    for (int p = 0; p < s->nip; p++) { free(s->ip[p].comp); free(s->ip[p].idx); free(s->ip[p].w); free(s->ip[p].acc); }
    free(s->src_acc); free(s->freq); free(s->phE); free(s->phH); free(s->rotE);
    cem_modeports_free(s);
    cem_disp_free(s);
    cem_sibc_free(s);
    free(s);
}

int cem_set_boundary(cem_sim *s, int face, int type, int ncells) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (face < 0 || face > 5) { set_err("bad face"); return -1; }
    int a = face / 2;
    if (type == CEM_BC_CPML && (ncells < 1 || ncells * 2 >= s->n[a])) { set_err("bad CPML thickness"); return -1; }
    s->bc[face] = type;
    s->npml[face] = type == CEM_BC_CPML ? ncells : 0;
    return 0;
}

/* ---- geometry ------------------------------------------------------------ */

static int nearest_node(const cem_sim *s, int a, double pos) {
    int best = 0;
    double bd = fabs(s->node[a][0] - pos);
    for (int i = 1; i <= s->n[a]; i++) {
        double d = fabs(s->node[a][i] - pos);
        if (d < bd) { bd = d; best = i; }
    }
    return best;
}

static int cell_range(const cem_sim *s, int a, double lo, double hi, int *i0, int *i1) {
    /* cells whose centre lies in [lo, hi] */
    *i0 = s->n[a]; *i1 = -1;
    for (int i = 0; i < s->n[a]; i++) {
        double c = 0.5 * (s->node[a][i] + s->node[a][i + 1]);
        if (c >= lo && c <= hi) {
            if (i < *i0) *i0 = i;
            if (i > *i1) *i1 = i;
        }
    }
    return *i1 >= *i0;
}

int cem_material_box(cem_sim *s, double x0, double x1, double y0, double y1, double z0, double z1,
                     double eps_r, double sigma) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    int r[3][2];
    double lo[3] = {x0, y0, z0}, hi[3] = {x1, y1, z1};
    for (int a = 0; a < 3; a++)
        if (!cell_range(s, a, lo[a], hi[a], &r[a][0], &r[a][1])) return 0;
    for (int i = r[0][0]; i <= r[0][1]; i++)
        for (int j = r[1][0]; j <= r[1][1]; j++)
            for (int k = r[2][0]; k <= r[2][1]; k++) {
                size_t q = IDX(s, i, j, k);
                s->eps_r[q] = eps_r;
                s->sigma[q] = sigma;
            }
    return 0;
}

int cem_pec_box(cem_sim *s, double x0, double x1, double y0, double y1, double z0, double z1) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    int r[3][2];
    double lo[3] = {x0, y0, z0}, hi[3] = {x1, y1, z1};
    for (int a = 0; a < 3; a++)
        if (!cell_range(s, a, lo[a], hi[a], &r[a][0], &r[a][1])) return 0;
    for (int i = r[0][0]; i <= r[0][1]; i++)
        for (int j = r[1][0]; j <= r[1][1]; j++)
            for (int k = r[2][0]; k <= r[2][1]; k++) s->cell_pec[IDX(s, i, j, k)] = 1;
    return 0;
}

int cem_pec_sheet(cem_sim *s, int axis, double pos, double a0, double a1, double b0, double b1) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (axis < 0 || axis > 2) { set_err("bad axis"); return -1; }
    int m = nearest_node(s, axis, pos);
    int ua = (axis + 1) % 3, va = (axis + 2) % 3;
    double tol = 1e-9 * (s->node[ua][s->n[ua]] - s->node[ua][0] + s->node[va][s->n[va]] - s->node[va][0]);
    /* edges along ua: span [u, u+1] in ua, node v in va */
    for (int u = 0; u < s->n[ua]; u++) {
        if (s->node[ua][u] < a0 - tol || s->node[ua][u + 1] > a1 + tol) continue;
        for (int v = 0; v <= s->n[va]; v++) {
            if (s->node[va][v] < b0 - tol || s->node[va][v] > b1 + tol) continue;
            int ijk[3];
            ijk[axis] = m; ijk[ua] = u; ijk[va] = v;
            s->edge_pec[ua][IDX(s, ijk[0], ijk[1], ijk[2])] = 1;
        }
    }
    for (int u = 0; u <= s->n[ua]; u++) {
        if (s->node[ua][u] < a0 - tol || s->node[ua][u] > a1 + tol) continue;
        for (int v = 0; v < s->n[va]; v++) {
            if (s->node[va][v] < b0 - tol || s->node[va][v + 1] > b1 + tol) continue;
            int ijk[3];
            ijk[axis] = m; ijk[ua] = u; ijk[va] = v;
            s->edge_pec[va][IDX(s, ijk[0], ijk[1], ijk[2])] = 1;
        }
    }
    return 0;
}

/* Mark PEC edges on the node plane m normal to `axis` from masks (rasterised copper).
 * mask_u: edges along ua (n[ua] x (n[va]+1), ua index outer); mask_v: edges along
 * va ((n[ua]+1) x n[va]). ua, va are the cyclic remaining axes. */
int cem_pec_edges_plane(cem_sim *s, int axis, int m, const unsigned char *mask_u, const unsigned char *mask_v) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (axis < 0 || axis > 2 || m < 0 || m > s->n[axis]) { set_err("bad plane"); return -1; }
    int ua = (axis + 1) % 3, va = (axis + 2) % 3;
    int nu = s->n[ua], nv = s->n[va];
    int ijk[3];
    ijk[axis] = m;
    if (mask_u)
        for (int u = 0; u < nu; u++)
            for (int v = 0; v <= nv; v++)
                if (mask_u[(size_t)u * (nv + 1) + v]) {
                    ijk[ua] = u; ijk[va] = v;
                    s->edge_pec[ua][IDX(s, ijk[0], ijk[1], ijk[2])] = 1;
                }
    if (mask_v)
        for (int u = 0; u <= nu; u++)
            for (int v = 0; v < nv; v++)
                if (mask_v[(size_t)u * nv + v]) {
                    ijk[ua] = u; ijk[va] = v;
                    s->edge_pec[va][IDX(s, ijk[0], ijk[1], ijk[2])] = 1;
                }
    return 0;
}

/* Mark solid PEC cells from a mask over the cell box starting at (i0, j0, k0) of
 * size ni x nj x nk (k fastest). */
int cem_pec_cells(cem_sim *s, int i0, int j0, int k0, int ni, int nj, int nk, const unsigned char *mask) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (i0 < 0 || j0 < 0 || k0 < 0 || i0 + ni > s->n[0] || j0 + nj > s->n[1] || k0 + nk > s->n[2]) {
        set_err("bad cell box"); return -1;
    }
    for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
            for (int k = 0; k < nk; k++)
                if (mask[((size_t)i * nj + j) * nk + k]) s->cell_pec[IDX(s, i0 + i, j0 + j, k0 + k)] = 1;
    return 0;
}

/* ---- ports, probes ------------------------------------------------------- */

static int check_edge_column(const cem_sim *s, int axis, int n1, int n2, int m0, int m1) {
    int ua = (axis + 1) % 3, va = (axis + 2) % 3;
    if (axis < 0 || axis > 2) return 0;
    if (n1 < 0 || n1 > s->n[ua] || n2 < 0 || n2 > s->n[va]) return 0;
    if (m0 < 0 || m1 > s->n[axis] || m1 <= m0) return 0;
    return 1;
}

static size_t col_index(const cem_sim *s, int axis, int n1, int n2, int m) {
    int ijk[3];
    ijk[axis] = m; ijk[(axis + 1) % 3] = n1; ijk[(axis + 2) % 3] = n2;
    return IDX(s, ijk[0], ijk[1], ijk[2]);
}

int cem_add_lumped_port(cem_sim *s, int axis, int n1, int n2, int m0, int m1, int dir, double R) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (s->nport >= MAXPORTS) { set_err("too many ports"); return -1; }
    if (!check_edge_column(s, axis, n1, n2, m0, m1) || !(R > 0) || (dir != 1 && dir != -1)) {
        set_err("bad port definition"); return -1;
    }
    port_t *p = &s->port[s->nport];
    memset(p, 0, sizeof *p);
    p->axis = axis; p->dir = dir; p->R = R;
    p->nedge = m1 - m0;
    p->edge = xcalloc(p->nedge, sizeof(size_t));
    p->dl = xcalloc(p->nedge, sizeof(double));
    p->cs = xcalloc(p->nedge, sizeof(double));
    for (int e = 0; e < p->nedge; e++) {
        p->edge[e] = col_index(s, axis, n1, n2, m0 + e);
        p->dl[e] = s->d[axis][m0 + e];
    }
    return s->nport++;
}

int cem_set_gaussian(cem_sim *s, int port, double amp, double t0, double tau) {
    if (port < 0 || port >= s->nport) { set_err("bad port"); return -1; }
    s->src_port = port; s->src_amp = amp; s->src_t0 = t0; s->src_tau = tau;
    return 0;
}

int cem_set_frequencies(cem_sim *s, int nf, const double *f) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (nf < 1) { set_err("need at least one frequency"); return -1; }
    free(s->freq);
    s->nf = nf;
    s->freq = xcalloc(nf, sizeof(double));
    memcpy(s->freq, f, nf * sizeof(double));
    return 0;
}

int cem_add_vprobe(cem_sim *s, int axis, int n1, int n2, int m0, int m1, int dir) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (s->nvp >= MAXPROBES) { set_err("too many probes"); return -1; }
    if (!check_edge_column(s, axis, n1, n2, m0, m1)) { set_err("bad probe"); return -1; }
    vprobe_t *p = &s->vp[s->nvp];
    p->axis = axis; p->dir = dir; p->n = m1 - m0;
    p->edge = xcalloc(p->n, sizeof(size_t));
    p->dl = xcalloc(p->n, sizeof(double));
    for (int e = 0; e < p->n; e++) {
        p->edge[e] = col_index(s, axis, n1, n2, m0 + e);
        p->dl[e] = s->d[axis][m0 + e] * dir;
    }
    return s->nvp++;
}

int cem_add_iprobe(cem_sim *s, int axis, int m, int u0, int u1, int v0, int v1) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (s->nip >= MAXPROBES) { set_err("too many probes"); return -1; }
    int ua = (axis + 1) % 3, va = (axis + 2) % 3;
    if (m < 0 || m >= s->n[axis] || u0 < 0 || u1 >= s->n[ua] || u1 <= u0 || v0 < 0 ||
        v1 >= s->n[va] || v1 <= v0) {
        set_err("bad current probe"); return -1;
    }
    iprobe_t *p = &s->ip[s->nip];
    int nmax = 2 * (u1 - u0) + 2 * (v1 - v0);
    p->comp = xcalloc(nmax, sizeof(int));
    p->idx = xcalloc(nmax, sizeof(size_t));
    p->w = xcalloc(nmax, sizeof(double));
    int n = 0, ijk[3];
    ijk[axis] = m;
    /* bottom (v = v0+1/2) +u and top (v = v1+1/2) -u : H along ua at nodes u */
    for (int u = u0 + 1; u <= u1; u++) {
        ijk[ua] = u; ijk[va] = v0;
        p->comp[n] = 3 + ua; p->idx[n] = IDX(s, ijk[0], ijk[1], ijk[2]); p->w[n] = s->dd[ua][u]; n++;
        ijk[va] = v1;
        p->comp[n] = 3 + ua; p->idx[n] = IDX(s, ijk[0], ijk[1], ijk[2]); p->w[n] = -s->dd[ua][u]; n++;
    }
    /* right (u = u1+1/2) +v and left (u = u0+1/2) -v : H along va at nodes v */
    for (int v = v0 + 1; v <= v1; v++) {
        ijk[va] = v; ijk[ua] = u1;
        p->comp[n] = 3 + va; p->idx[n] = IDX(s, ijk[0], ijk[1], ijk[2]); p->w[n] = s->dd[va][v]; n++;
        ijk[ua] = u0;
        p->comp[n] = 3 + va; p->idx[n] = IDX(s, ijk[0], ijk[1], ijk[2]); p->w[n] = -s->dd[va][v]; n++;
    }
    p->n = n;
    return s->nip++;
}

/* ---- finalize ------------------------------------------------------------ */

/* Average cell material onto the edge along axis c at index (i,j,k): the four
 * cells sharing the edge, area weighted. */
static void edge_material(const cem_sim *s, int c, int i, int j, int k, double *eps, double *sig,
                          int *pec) {
    int ua = (c + 1) % 3, va = (c + 2) % 3;
    int ijk[3] = {i, j, k};
    double wsum = 0, es = 0, ss = 0;
    *pec = 0;
    for (int du = -1; du <= 0; du++)
        for (int dv = -1; dv <= 0; dv++) {
            int cc[3] = {ijk[0], ijk[1], ijk[2]};
            cc[ua] += du; cc[va] += dv;
            if (cc[ua] < 0 || cc[ua] >= s->n[ua] || cc[va] < 0 || cc[va] >= s->n[va]) continue;
            if (cc[c] >= s->n[c]) continue;
            size_t q = IDX(s, cc[0], cc[1], cc[2]);
            double w = s->d[ua][cc[ua]] * s->d[va][cc[va]];
            if (s->cell_pec[q]) *pec = 1;
            es += w * s->eps_r[q];
            ss += w * s->sigma[q];
            wsum += w;
        }
    *eps = wsum > 0 ? es / wsum : 1.0;
    *sig = wsum > 0 ? ss / wsum : 0.0;
}

static void setup_pml_axis(cem_sim *s, int a) {
    int n = s->n[a];
    s->ide[a] = xcalloc(n + 1, sizeof(double));
    s->idh[a] = xcalloc(n, sizeof(double));
    s->be[a] = xcalloc(n + 1, sizeof(double));
    s->ce[a] = xcalloc(n + 1, sizeof(double));
    s->bh[a] = xcalloc(n, sizeof(double));
    s->ch[a] = xcalloc(n, sizeof(double));
    s->pml_sn[a] = xcalloc(n + 1, sizeof(double)); s->pml_an[a] = xcalloc(n + 1, sizeof(double));
    s->pml_kn[a] = xcalloc(n + 1, sizeof(double));
    s->pml_sc[a] = xcalloc(n, sizeof(double)); s->pml_ac[a] = xcalloc(n, sizeof(double));
    s->pml_kc[a] = xcalloc(n, sizeof(double));
    for (int i = 0; i <= n; i++) s->pml_kn[a][i] = 1;
    for (int i = 0; i < n; i++) s->pml_kc[a][i] = 1;
    for (int i = 0; i <= n; i++) s->ide[a][i] = s->dd[a][i] > 0 ? 1.0 / s->dd[a][i] : 0.0;
    for (int i = 0; i < n; i++) s->idh[a][i] = 1.0 / s->d[a][i];

    const double m = 3.0, kmax = 5.0, amax = 0.05;
    for (int side = 0; side < 2; side++) {
        int np = s->npml[2 * a + side];
        if (!np) continue;
        double xi, L, dmean;
        if (side == 0) {
            xi = s->node[a][np];
            L = s->node[a][np] - s->node[a][0];
        } else {
            xi = s->node[a][n - np];
            L = s->node[a][n] - s->node[a][n - np];
        }
        dmean = L / np;
        double smax = 0.8 * (m + 1) / (376.730313668 * dmean);
        /* E nodes */
        for (int i = 0; i <= n; i++) {
            double depth = side == 0 ? (xi - s->node[a][i]) / L : (s->node[a][i] - xi) / L;
            if (depth <= 0) continue;
            if (depth > 1) depth = 1;
            double sg = smax * pow(depth, m), kp = 1 + (kmax - 1) * pow(depth, m), al = amax * (1 - depth);
            double b = exp(-(sg / kp + al) * s->dt / CEM_EPS0);
            s->be[a][i] = b;
            s->ce[a][i] = sg > 0 ? sg / (sg * kp + kp * kp * al) * (b - 1) : 0;
            s->ide[a][i] = 1.0 / (kp * s->dd[a][i]);
            s->pml_sn[a][i] = sg; s->pml_kn[a][i] = kp; s->pml_an[a][i] = al;
        }
        for (int i = 0; i < n; i++) {
            double cpos = 0.5 * (s->node[a][i] + s->node[a][i + 1]);
            double depth = side == 0 ? (xi - cpos) / L : (cpos - xi) / L;
            if (depth <= 0) continue;
            if (depth > 1) depth = 1;
            double sg = smax * pow(depth, m), kp = 1 + (kmax - 1) * pow(depth, m), al = amax * (1 - depth);
            double b = exp(-(sg / kp + al) * s->dt / CEM_EPS0);
            s->bh[a][i] = b;
            s->ch[a][i] = sg > 0 ? sg / (sg * kp + kp * kp * al) * (b - 1) : 0;
            s->idh[a][i] = 1.0 / (kp * s->d[a][i]);
            s->pml_sc[a][i] = sg; s->pml_kc[a][i] = kp; s->pml_ac[a][i] = al;
        }
    }
}

int cem_finalize(cem_sim *s, double cfl) {
    if (s->finalized) { set_err("already finalized"); return -1; }
    if (!s->nf) { double f = 1e9; cem_set_frequencies(s, 1, &f); }
    double dmin[3];
    for (int a = 0; a < 3; a++) {
        dmin[a] = 1e300;
        for (int i = 0; i < s->n[a]; i++) if (s->d[a][i] < dmin[a]) dmin[a] = s->d[a][i];
    }
    if (!(cfl > 0 && cfl <= 1)) cfl = 0.99;
    s->dt = cfl / (CEM_C0 * sqrt(1 / (dmin[0] * dmin[0]) + 1 / (dmin[1] * dmin[1]) + 1 / (dmin[2] * dmin[2])));

    for (int a = 0; a < 3; a++) setup_pml_axis(s, a);
    s->has_pml = 0;
    for (int f = 0; f < 6; f++) if (s->npml[f]) s->has_pml = 1;
    if (s->has_pml)
        for (int c = 0; c < 3; c++)
            for (int w = 0; w < 2; w++) {
                s->psiE[c][w] = xcalloc(cem_psi_len(s, 1, c, w) + 1, sizeof(double));
                s->psiH[c][w] = xcalloc(cem_psi_len(s, 0, c, w) + 1, sizeof(double));
                if (!s->psiE[c][w] || !s->psiH[c][w]) { set_err("out of memory"); return -1; }
            }

    /* E coefficients */
    int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    for (int c = 0; c < 3; c++)
        for (int i = 0; i <= nx; i++)
            for (int j = 0; j <= ny; j++)
                for (int k = 0; k <= nz; k++) {
                    int ijk[3] = {i, j, k};
                    if (ijk[c] >= s->n[c]) continue; /* no such edge */
                    size_t q = IDX(s, i, j, k);
                    double er, sg;
                    int pec;
                    edge_material(s, c, i, j, k, &er, &sg, &pec);
                    if (pec || s->edge_pec[c][q]) { s->ca[c][q] = 0; s->cb[c][q] = 0; continue; }
                    double eps = er * CEM_EPS0, h = sg * s->dt / (2 * eps);
                    s->ca[c][q] = (1 - h) / (1 + h);
                    s->cb[c][q] = (s->dt / eps) / (1 + h);
                }

    /* lumped ports: semi-implicit resistive voltage source on each edge */
    for (int p = 0; p < s->nport; p++) {
        port_t *pt = &s->port[p];
        double L = 0;
        for (int e = 0; e < pt->nedge; e++) L += pt->dl[e];
        for (int e = 0; e < pt->nedge; e++) {
            size_t q = pt->edge[e];
            int i = (int)(q / s->si), j = (int)((q % s->si) / s->sj), k = (int)(q % s->sj);
            int ijk[3] = {i, j, k};
            int ua = (pt->axis + 1) % 3, va = (pt->axis + 2) % 3;
            double A = s->dd[ua][ijk[ua]] * s->dd[va][ijk[va]];
            double er, sg;
            int pec;
            edge_material(s, pt->axis, i, j, k, &er, &sg, &pec);
            if (pec || s->edge_pec[pt->axis][q]) { set_err("port edge is PEC"); return -1; }
            double eps = er * CEM_EPS0;
            double Re = pt->R * pt->dl[e] / L;
            double beta = s->dt * pt->dl[e] / (2 * Re * A * eps) + sg * s->dt / (2 * eps);
            s->ca[pt->axis][q] = (1 - beta) / (1 + beta);
            s->cb[pt->axis][q] = (s->dt / eps) / (1 + beta);
            /* source term per unit total source voltage Vs: Vs_e = Vs dl/L */
            pt->cs[e] = s->dt / (Re * A * eps * (1 + beta)) * (pt->dl[e] / L) * pt->dir;
        }
        pt->V = xcalloc(2 * s->nf, sizeof(double));
        pt->I = xcalloc(2 * s->nf, sizeof(double));
    }
    /* Debye dispersion: modify ca/cb of dispersive edges and set up the ADE state.
     * Lumped-port edges keep their (non-dispersive, eps_inf) coefficients. */
    unsigned char *is_port = NULL;
    if (s->nport && (s->npole || s->sibc_on)) {
        is_port = xcalloc(s->N * 3, 1);
        for (int p = 0; p < s->nport; p++)
            for (int e = 0; e < s->port[p].nedge; e++) is_port[s->port[p].axis * s->N + s->port[p].edge[e]] = 1;
    }
    if (s->npole) {
        double dep[MAXPOLES];
        for (int pass = 0; pass < 2; pass++) {
            int nd = 0;
            for (int c = 0; c < 3; c++)
                for (int i = 0; i <= nx; i++)
                    for (int j = 0; j <= ny; j++)
                        for (int k = 0; k <= nz; k++) {
                            int ijk[3] = {i, j, k};
                            if (ijk[c] >= s->n[c]) continue;
                            size_t q = IDX(s, i, j, k);
                            if (s->cb[c][q] == 0 || (is_port && is_port[c * s->N + q])) continue;
                            if (!cem_edge_debye(s, c, i, j, k, dep)) continue;
                            if (pass == 1) {
                                double er, sg;
                                int pec;
                                edge_material(s, c, i, j, k, &er, &sg, &pec);
                                double einst = er * CEM_EPS0;
                                for (int p = 0; p < s->npole; p++) {
                                    double tau = 1.0 / (2 * M_PI * s->pole_f[p]);
                                    double bt = CEM_EPS0 * dep[p] * s->dt / (tau * (1 + s->dt / (2 * tau)));
                                    s->d_beta[(size_t)p * s->ndisp + nd] = bt;
                                    einst += 0.5 * bt;
                                }
                                double h = sg * s->dt / (2 * einst);
                                s->ca[c][q] = (1 - h) / (1 + h);
                                s->cb[c][q] = (s->dt / einst) / (1 + h);
                                s->d_comp[nd] = c;
                                s->d_idx[nd] = q;
                                s->d_cb[nd] = s->cb[c][q];
                            }
                            nd++;
                        }
            if (pass == 0) {
                s->ndisp = nd;
                s->d_comp = xcalloc(nd, sizeof(int));
                s->d_idx = xcalloc(nd, sizeof(size_t));
                s->d_cb = xcalloc(nd, sizeof(double));
                s->d_Eold = xcalloc(nd, sizeof(double));
                s->d_beta = xcalloc((size_t)nd * s->npole, sizeof(double));
                s->d_J = xcalloc((size_t)nd * s->npole, sizeof(double));
                if (!s->d_J || !s->d_beta) { set_err("out of memory"); free(is_port); return -1; }
            }
        }
        for (int p = 0; p < s->npole; p++) {
            double tau = 1.0 / (2 * M_PI * s->pole_f[p]);
            s->d_k[p] = (1 - s->dt / (2 * tau)) / (1 + s->dt / (2 * tau));
        }
    }
    if (cem_sibc_finalize(s, is_port)) { free(is_port); return -1; }
    free(is_port);
    for (int p = 0; p < s->nvp; p++) s->vp[p].acc = xcalloc(2 * s->nf, sizeof(double));
    for (int p = 0; p < s->nip; p++) s->ip[p].acc = xcalloc(2 * s->nf, sizeof(double));
    s->src_acc = xcalloc(2 * s->nf, sizeof(double));
    s->phE = xcalloc(2 * s->nf, sizeof(double));
    s->phH = xcalloc(2 * s->nf, sizeof(double));
    s->rotE = xcalloc(2 * s->nf, sizeof(double));
    for (int f = 0; f < s->nf; f++) {
        double w = 2 * M_PI * s->freq[f];
        s->rotE[2 * f] = cos(w * s->dt);
        s->rotE[2 * f + 1] = -sin(w * s->dt);
    }
    s->step = 0;
    s->finalized = 1;
    if (cem_modeports_finalize(s)) { s->finalized = 0; return -1; }
    return 0;
}

/* ---- updates ------------------------------------------------------------- */

static void pml_H(cem_sim *s, int i0, int i1);
static void pml_E(cem_sim *s, int i0, int i1);
static inline int imax(int a, int b) { return a > b ? a : b; }
static inline int imin(int a, int b) { return a < b ? a : b; }

/* Bulk + CPML updates over the slab i0 <= i < i1 (node index i). Slabs are
 * independent within one half step, so they can run on separate threads;
 * each cell sees the same operations in the same order as a whole-domain pass. */
void cem_mt_H(cem_sim *s, int i0, int i1) {
    const int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    const size_t si = s->si, sj = s->sj;
    const double dth = s->dt / CEM_MU0;
    double *Ex = s->E[0], *Ey = s->E[1], *Ez = s->E[2];
    double *Hx = s->H[0], *Hy = s->H[1], *Hz = s->H[2];
    const double *ix = s->idh[0], *iy = s->idh[1], *iz = s->idh[2];

    for (int i = imax(0, i0); i < imin(nx + 1, i1); i++)
        for (int j = 0; j < ny; j++) {
            size_t q = i * si + j * sj;
            for (int k = 0; k < nz; k++, q++)
                Hx[q] -= dth * ((Ez[q + sj] - Ez[q]) * iy[j] - (Ey[q + 1] - Ey[q]) * iz[k]);
        }
    for (int i = imax(0, i0); i < imin(nx, i1); i++)
        for (int j = 0; j <= ny; j++) {
            size_t q = i * si + j * sj;
            for (int k = 0; k < nz; k++, q++)
                Hy[q] -= dth * ((Ex[q + 1] - Ex[q]) * iz[k] - (Ez[q + si] - Ez[q]) * ix[i]);
        }
    for (int i = imax(0, i0); i < imin(nx, i1); i++)
        for (int j = 0; j < ny; j++) {
            size_t q = i * si + j * sj;
            for (int k = 0; k <= nz; k++, q++)
                Hz[q] -= dth * ((Ey[q + si] - Ey[q]) * ix[i] - (Ex[q + sj] - Ex[q]) * iy[j]);
        }
    if (s->has_pml) pml_H(s, i0, i1);
}

void cem_mt_E(cem_sim *s, int i0, int i1) {
    const int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    const size_t si = s->si, sj = s->sj;
    double *Ex = s->E[0], *Ey = s->E[1], *Ez = s->E[2];
    double *Hx = s->H[0], *Hy = s->H[1], *Hz = s->H[2];
    const double *ix = s->ide[0], *iy = s->ide[1], *iz = s->ide[2];

    for (int i = imax(0, i0); i < imin(nx, i1); i++)
        for (int j = 1; j < ny; j++) {
            size_t q = i * si + j * sj + 1;
            const double *ca = s->ca[0], *cb = s->cb[0];
            for (int k = 1; k < nz; k++, q++)
                Ex[q] = ca[q] * Ex[q] + cb[q] * ((Hz[q] - Hz[q - sj]) * iy[j] - (Hy[q] - Hy[q - 1]) * iz[k]);
        }
    for (int i = imax(1, i0); i < imin(nx, i1); i++)
        for (int j = 0; j < ny; j++) {
            size_t q = i * si + j * sj + 1;
            const double *ca = s->ca[1], *cb = s->cb[1];
            for (int k = 1; k < nz; k++, q++)
                Ey[q] = ca[q] * Ey[q] + cb[q] * ((Hx[q] - Hx[q - 1]) * iz[k] - (Hz[q] - Hz[q - si]) * ix[i]);
        }
    for (int i = imax(1, i0); i < imin(nx, i1); i++)
        for (int j = 1; j < ny; j++) {
            size_t q = i * si + j * sj;
            const double *ca = s->ca[2], *cb = s->cb[2];
            for (int k = 0; k < nz; k++, q++)
                Ez[q] = ca[q] * Ez[q] + cb[q] * ((Hy[q] - Hy[q - si]) * ix[i] - (Hx[q] - Hx[q - sj]) * iy[j]);
        }
    if (s->has_pml) pml_E(s, i0, i1);
}

/* CPML auxiliary updates. For each derivative direction a, only the slabs
 * inside the layers are visited. */
static void pml_ranges(const cem_sim *s, int a, int is_e, int r[2][2]) {
    int n = s->n[a];
    int lo = s->npml[2 * a], hi = s->npml[2 * a + 1];
    if (is_e) { /* nodes 1..n-1 */
        r[0][0] = 1; r[0][1] = lo;          /* inclusive end, exclusive-safe below */
        r[1][0] = n - hi; r[1][1] = n - 1;
        if (lo == 0) { r[0][0] = 1; r[0][1] = 0; }
        if (hi == 0) { r[1][0] = 1; r[1][1] = 0; }
    } else { /* cells 0..n-1 */
        r[0][0] = 0; r[0][1] = lo - 1;
        r[1][0] = n - hi; r[1][1] = n - 1;
    }
}

/* Compact CPML storage: the auxiliaries for derivative axis a exist only on
 * that axis's PML planes. Plane order is the low slab then the high slab; the
 * in-plane index runs over the other two axes, (j,k) / (i,k) / (i,j), with the
 * full node extent of each. psiX[c][w] belongs to axis (c + w + 1) % 3. */
static size_t pml_per(const cem_sim *s, int a) {
    const int *n = s->n;
    return a == 0 ? (size_t)(n[1] + 1) * (n[2] + 1) : a == 1 ? (size_t)(n[0] + 1) * (n[2] + 1) : (size_t)(n[0] + 1) * (n[1] + 1);
}
static size_t pml_planes(const cem_sim *s, int a, int is_e) {
    int r[2][2];
    pml_ranges(s, a, is_e, r);
    size_t n = 0;
    for (int sd = 0; sd < 2; sd++) if (r[sd][1] >= r[sd][0]) n += (size_t)(r[sd][1] - r[sd][0] + 1);
    return n;
}
static size_t pml_base(const cem_sim *s, int a, int is_e, int sd, int idx, const int r[2][2]) {
    size_t lo = r[0][1] >= r[0][0] ? (size_t)(r[0][1] - r[0][0] + 1) : 0;
    size_t plane = sd == 0 ? (size_t)(idx - r[0][0]) : lo + (size_t)(idx - r[1][0]);
    (void)is_e;
    return plane * pml_per(s, a);
}
size_t cem_psi_len(const cem_sim *s, int is_e, int comp, int w) {
    return s->has_pml ? pml_planes(s, (comp + w + 1) % 3, is_e) * pml_per(s, (comp + w + 1) % 3) : 0;
}

static void pml_H(cem_sim *s, int i0, int i1) {
    const int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    const size_t si = s->si, sj = s->sj;
    const double dth = s->dt / CEM_MU0;
    double *Ex = s->E[0], *Ey = s->E[1], *Ez = s->E[2];
    double *Hx = s->H[0], *Hy = s->H[1], *Hz = s->H[2];
    int r[2][2];

    /* derivative along x: Hy (dEz/dx, sign +), Hz (dEy/dx, sign -) ; cells i */
    pml_ranges(s, 0, 0, r);
    for (int sd = 0; sd < 2; sd++)
        for (int i = imax(r[sd][0], i0); i <= imin(r[sd][1], i1 - 1); i++) {
            double b = s->bh[0][i], c = s->ch[0][i], id = 1.0 / s->d[0][i];
            size_t base = pml_base(s, 0, 0, sd, i, (const int (*)[2])r);
            double *pHy = s->psiH[1][1], *pHz = s->psiH[2][0];
            for (int j = 0; j <= ny; j++)
                for (int k = 0; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)j * (nz + 1) + k;
                    pHy[p] = b * pHy[p] + c * (Ez[q + si] - Ez[q]) * id;
                    Hy[q] += dth * pHy[p];
                }
            for (int j = 0; j < ny; j++)
                for (int k = 0; k <= nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)j * (nz + 1) + k;
                    pHz[p] = b * pHz[p] + c * (Ey[q + si] - Ey[q]) * id;
                    Hz[q] -= dth * pHz[p];
                }
        }
    /* derivative along y: Hz (dEx/dy, +), Hx (dEz/dy, -) */
    pml_ranges(s, 1, 0, r);
    for (int sd = 0; sd < 2; sd++)
        for (int j = r[sd][0]; j <= r[sd][1]; j++) {
            double b = s->bh[1][j], c = s->ch[1][j], id = 1.0 / s->d[1][j];
            size_t base = pml_base(s, 1, 0, sd, j, (const int (*)[2])r);
            double *pHz = s->psiH[2][1], *pHx = s->psiH[0][0];
            for (int i = imax(0, i0); i < imin(nx, i1); i++)
                for (int k = 0; k <= nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (nz + 1) + k;
                    pHz[p] = b * pHz[p] + c * (Ex[q + sj] - Ex[q]) * id;
                    Hz[q] += dth * pHz[p];
                }
            for (int i = imax(0, i0); i < imin(nx + 1, i1); i++)
                for (int k = 0; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (nz + 1) + k;
                    pHx[p] = b * pHx[p] + c * (Ez[q + sj] - Ez[q]) * id;
                    Hx[q] -= dth * pHx[p];
                }
        }
    /* derivative along z: Hx (dEy/dz, +), Hy (dEx/dz, -) */
    pml_ranges(s, 2, 0, r);
    for (int sd = 0; sd < 2; sd++)
        for (int k = r[sd][0]; k <= r[sd][1]; k++) {
            double b = s->bh[2][k], c = s->ch[2][k], id = 1.0 / s->d[2][k];
            size_t base = pml_base(s, 2, 0, sd, k, (const int (*)[2])r);
            double *pHx = s->psiH[0][1], *pHy = s->psiH[1][0];
            for (int i = imax(0, i0); i < imin(nx + 1, i1); i++)
                for (int j = 0; j < ny; j++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (ny + 1) + j;
                    pHx[p] = b * pHx[p] + c * (Ey[q + 1] - Ey[q]) * id;
                    Hx[q] += dth * pHx[p];
                }
            for (int i = imax(0, i0); i < imin(nx, i1); i++)
                for (int j = 0; j <= ny; j++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (ny + 1) + j;
                    pHy[p] = b * pHy[p] + c * (Ex[q + 1] - Ex[q]) * id;
                    Hy[q] -= dth * pHy[p];
                }
        }
}

static void pml_E(cem_sim *s, int i0, int i1) {
    const int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    const size_t si = s->si, sj = s->sj;
    double *Ex = s->E[0], *Ey = s->E[1], *Ez = s->E[2];
    double *Hx = s->H[0], *Hy = s->H[1], *Hz = s->H[2];
    int r[2][2];

    /* along x: Ey (-dHz/dx), Ez (+dHy/dx); nodes i in 1..nx-1 */
    pml_ranges(s, 0, 1, r);
    for (int sd = 0; sd < 2; sd++)
        for (int i = imax(r[sd][0], i0); i <= imin(r[sd][1], i1 - 1); i++) {
            double b = s->be[0][i], c = s->ce[0][i], id = 1.0 / s->dd[0][i];
            size_t base = pml_base(s, 0, 1, sd, i, (const int (*)[2])r);
            double *pEy = s->psiE[1][1], *pEz = s->psiE[2][0];
            const double *cby = s->cb[1], *cbz = s->cb[2];
            for (int j = 0; j < ny; j++)
                for (int k = 1; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)j * (nz + 1) + k;
                    pEy[p] = b * pEy[p] + c * (Hz[q] - Hz[q - si]) * id;
                    Ey[q] -= cby[q] * pEy[p];
                }
            for (int j = 1; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)j * (nz + 1) + k;
                    pEz[p] = b * pEz[p] + c * (Hy[q] - Hy[q - si]) * id;
                    Ez[q] += cbz[q] * pEz[p];
                }
        }
    /* along y: Ez (-dHx/dy), Ex (+dHz/dy) */
    pml_ranges(s, 1, 1, r);
    for (int sd = 0; sd < 2; sd++)
        for (int j = r[sd][0]; j <= r[sd][1]; j++) {
            double b = s->be[1][j], c = s->ce[1][j], id = 1.0 / s->dd[1][j];
            size_t base = pml_base(s, 1, 1, sd, j, (const int (*)[2])r);
            double *pEz = s->psiE[2][1], *pEx = s->psiE[0][0];
            const double *cbz = s->cb[2], *cbx = s->cb[0];
            for (int i = imax(1, i0); i < imin(nx, i1); i++)
                for (int k = 0; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (nz + 1) + k;
                    pEz[p] = b * pEz[p] + c * (Hx[q] - Hx[q - sj]) * id;
                    Ez[q] -= cbz[q] * pEz[p];
                }
            for (int i = imax(0, i0); i < imin(nx, i1); i++)
                for (int k = 1; k < nz; k++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (nz + 1) + k;
                    pEx[p] = b * pEx[p] + c * (Hz[q] - Hz[q - sj]) * id;
                    Ex[q] += cbx[q] * pEx[p];
                }
        }
    /* along z: Ex (-dHy/dz), Ey (+dHx/dz) */
    pml_ranges(s, 2, 1, r);
    for (int sd = 0; sd < 2; sd++)
        for (int k = r[sd][0]; k <= r[sd][1]; k++) {
            double b = s->be[2][k], c = s->ce[2][k], id = 1.0 / s->dd[2][k];
            size_t base = pml_base(s, 2, 1, sd, k, (const int (*)[2])r);
            double *pEx = s->psiE[0][1], *pEy = s->psiE[1][0];
            const double *cbx = s->cb[0], *cby = s->cb[1];
            for (int i = imax(0, i0); i < imin(nx, i1); i++)
                for (int j = 1; j < ny; j++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (ny + 1) + j;
                    pEx[p] = b * pEx[p] + c * (Hy[q] - Hy[q - 1]) * id;
                    Ex[q] -= cbx[q] * pEx[p];
                }
            for (int i = imax(1, i0); i < imin(nx, i1); i++)
                for (int j = 0; j < ny; j++) {
                    size_t q = i * si + j * sj + k;
                    size_t p = base + (size_t)i * (ny + 1) + j;
                    pEy[p] = b * pEy[p] + c * (Hx[q] - Hx[q - 1]) * id;
                    Ey[q] += cby[q] * pEy[p];
                }
        }
}

static double source_value(const cem_sim *s, double t) {
    if (s->src_port < 0) return 0.0;
    double x = (t - s->src_t0) / s->src_tau;
    return s->src_amp * exp(-x * x);
}

static inline void dft_add(double *acc, const double *ph, int nf, double v) {
    for (int f = 0; f < nf; f++) {
        acc[2 * f] += v * ph[2 * f];
        acc[2 * f + 1] += v * ph[2 * f + 1];
    }
}

static void set_phasors(cem_sim *s, double tE, double tH) {
    for (int f = 0; f < s->nf; f++) {
        double w = 2 * M_PI * s->freq[f];
        s->phE[2 * f] = cos(w * tE); s->phE[2 * f + 1] = -sin(w * tE);
        s->phH[2 * f] = cos(w * tH); s->phH[2 * f + 1] = -sin(w * tH);
    }
}

static void rotate(double *ph, const double *rot, int nf) {
    for (int f = 0; f < nf; f++) {
        double a = ph[2 * f], b = ph[2 * f + 1], c = rot[2 * f], d = rot[2 * f + 1];
        ph[2 * f] = a * c - b * d;
        ph[2 * f + 1] = a * d + b * c;
    }
}

/* One time step in phases, so an outside thread pool can run the bulk updates
 * (cem_mt_H / cem_mt_E over i-slabs) in parallel:
 *   cem_mt_pre; cem_mt_H(0..nx+1); cem_mt_mid; cem_mt_E(0..nx); cem_mt_post.
 * cem_step runs exactly this sequence, so results are identical either way. */
void cem_mt_pre(cem_sim *s) {
    long n = s->step;
    if ((n & 1023) == 0) set_phasors(s, (n + 1) * s->dt, (n + 0.5) * s->dt);
}

void cem_mt_mid(cem_sim *s) {
    const int nf = s->nf;

    /* H-time probes */
    for (int p = 0; p < s->nip; p++) {
        iprobe_t *ip = &s->ip[p];
        double I = 0;
        for (int e = 0; e < ip->n; e++) I += s->H[ip->comp[e] - 3][ip->idx[e]] * ip->w[e];
        dft_add(ip->acc, s->phH, nf, I * s->dt);
    }

    /* port voltages before the E update (time n) */
    for (int p = 0; p < s->nport; p++) {
        port_t *pt = &s->port[p];
        double V = 0;
        for (int e = 0; e < pt->nedge; e++) V += s->E[pt->axis][pt->edge[e]] * pt->dl[e];
        s->vold[p] = V * pt->dir;
    }

    if (s->ndisp) cem_disp_save(s);
}

void cem_mt_post(cem_sim *s) {
    const int nf = s->nf;
    double tH = (s->step + 0.5) * s->dt;
    if (s->ndisp) cem_disp_correct(s);
    if (s->nsb) cem_sibc_apply(s);

    double Vs = source_value(s, tH);
    if (s->src_port >= 0) {
        port_t *pt = &s->port[s->src_port];
        for (int e = 0; e < pt->nedge; e++) s->E[pt->axis][pt->edge[e]] += pt->cs[e] * Vs;
        dft_add(s->src_acc, s->phH, nf, Vs * s->dt);
    }
    if (s->msrc >= 0) {
        cem_modeports_inject(s, tH);
        dft_add(s->src_acc, s->phH, nf, cem_modeport_source(s, tH) * s->dt);
    }
    if (s->ndisp) cem_disp_update(s);
    for (int p = 0; p < s->nport; p++) {
        port_t *pt = &s->port[p];
        double V = 0;
        for (int e = 0; e < pt->nedge; e++) V += s->E[pt->axis][pt->edge[e]] * pt->dl[e];
        V = 0.5 * (V * pt->dir + s->vold[p]);
        double vs = p == s->src_port ? Vs : 0.0;
        double I = (vs - V) / pt->R;
        dft_add(pt->V, s->phH, nf, V * s->dt);
        dft_add(pt->I, s->phH, nf, I * s->dt);
    }

    /* E-time probes */
    for (int p = 0; p < s->nvp; p++) {
        vprobe_t *vp = &s->vp[p];
        double V = 0;
        for (int e = 0; e < vp->n; e++) V += s->E[vp->axis][vp->edge[e]] * vp->dl[e];
        dft_add(vp->acc, s->phE, nf, V * s->dt);
    }

    rotate(s->phE, s->rotE, nf);
    rotate(s->phH, s->rotE, nf);
    s->step++;
}

/* Native builds with OpenMP split the slab phases across threads the same way
 * the browser's worker pool does; results are identical either way. */
static void slabs(cem_sim *s, int e) {
    int n = e ? s->n[0] : s->n[0] + 1;
#ifdef _OPENMP
#pragma omp parallel
    {
        int T = omp_get_num_threads(), id = omp_get_thread_num();
        int i0 = (int)((long)n * id / T), i1 = (int)((long)n * (id + 1) / T);
        if (e) cem_mt_E(s, i0, i1); else cem_mt_H(s, i0, i1);
    }
#else
    if (e) cem_mt_E(s, 0, n); else cem_mt_H(s, 0, n);
#endif
}

void cem_step(cem_sim *s, long nsteps) {
    if (!s->finalized) return;
    for (long it = 0; it < nsteps; it++) {
        cem_mt_pre(s);
        slabs(s, 0);
        cem_mt_mid(s);
        slabs(s, 1);
        cem_mt_post(s);
    }
}

const char *cem_version(void) {
#ifdef CEM_VERSION
    return CEM_VERSION;
#else
    return "dev";
#endif
}

double cem_energy(const cem_sim *s) {
    const int nx = s->n[0], ny = s->n[1], nz = s->n[2];
    double we = 0, wm = 0;
    for (int c = 0; c < 3; c++) {
        int ua = (c + 1) % 3, va = (c + 2) % 3;
        for (int i = 0; i <= nx; i++)
            for (int j = 0; j <= ny; j++)
                for (int k = 0; k <= nz; k++) {
                    int ijk[3] = {i, j, k};
                    size_t q = IDX(s, i, j, k);
                    /* E edge along c: length d[c], dual area */
                    if (ijk[c] < s->n[c] && s->cb[c][q] != 0) {
                        double vol = s->d[c][ijk[c]] * s->dd[ua][ijk[ua]] * s->dd[va][ijk[va]];
                        double eps = s->dt / s->cb[c][q]; /* approx (exact when lossless) */
                        we += 0.5 * eps * s->E[c][q] * s->E[c][q] * vol;
                    }
                    /* H face normal c: dual length along c, primary area */
                    if (ijk[ua] < s->n[ua] && ijk[va] < s->n[va]) {
                        double vol = s->dd[c][ijk[c]] * s->d[ua][ijk[ua]] * s->d[va][ijk[va]];
                        wm += 0.5 * CEM_MU0 * s->H[c][q] * s->H[c][q] * vol;
                    }
                }
    }
    return we + wm;
}

long cem_run(cem_sim *s, long max_steps, double decay, int check_every) {
    if (!s->finalized) return 0;
    if (check_every < 1) check_every = 200;
    long start = s->step;
    double t_src_end = s->src_port >= 0 ? s->src_t0 + 6 * s->src_tau : 0;
    if (s->msrc >= 0 && s->msrc_t0 + 6 * s->msrc_tau > t_src_end) t_src_end = s->msrc_t0 + 6 * s->msrc_tau;
    double peak = 0;
    while (s->step - start < max_steps) {
        long chunk = check_every;
        if (s->step - start + chunk > max_steps) chunk = max_steps - (s->step - start);
        cem_step(s, chunk);
        double W = cem_energy(s);
        if (W > peak) peak = W;
        if (s->step * s->dt > t_src_end && peak > 0 && W < decay * peak) break;
    }
    return s->step - start;
}

/* ---- results ------------------------------------------------------------- */

const double *cem_port_v(const cem_sim *s, int p) { return p >= 0 && p < s->nport ? s->port[p].V : NULL; }
const double *cem_port_i(const cem_sim *s, int p) { return p >= 0 && p < s->nport ? s->port[p].I : NULL; }
const double *cem_source_spectrum(const cem_sim *s) { return s->src_acc; }
const double *cem_vprobe(const cem_sim *s, int p) { return p >= 0 && p < s->nvp ? s->vp[p].acc : NULL; }
const double *cem_iprobe(const cem_sim *s, int p) { return p >= 0 && p < s->nip ? s->ip[p].acc : NULL; }
double cem_port_impedance(const cem_sim *s, int p) { return p >= 0 && p < s->nport ? s->port[p].R : 0; }
int cem_excited_port(const cem_sim *s) { return s->src_port; }

/* ---- cross-section extraction for the mode solver ------------------------ */

cem_mode2d *cem_cross_section_x(const cem_sim *s, int i, int j0, int j1, int k0, int k1) {
    if (!s->finalized || i < 0 || i >= s->n[0] || j0 < 0 || j1 > s->n[1] || k0 < 0 || k1 > s->n[2] ||
        j1 - j0 < 2 || k1 - k0 < 2) { set_err("bad cross-section window"); return NULL; }
    int ny = j1 - j0, nz = k1 - k0, np = s->npole;
    cem_mode2d *m = cem_mode2d_create(ny, nz, s->node[1] + j0, s->node[2] + k0);
    int nu = ny * (nz + 1), nv = (ny + 1) * nz, nn = (ny + 1) * (nz + 1);
    double *eu = malloc(nu * sizeof(double)), *ev = malloc(nv * sizeof(double)), *ew = malloc(nn * sizeof(double));
    double *su = calloc(nu, sizeof(double)), *sv = calloc(nv, sizeof(double)), *sw = calloc(nn, sizeof(double));
    unsigned char *pu = calloc(nu, 1), *pv = calloc(nv, 1), *pw = calloc(nn, 1);
    double *du = np ? calloc((size_t)np * nu, sizeof(double)) : NULL;
    double *dv = np ? calloc((size_t)np * nv, sizeof(double)) : NULL;
    double *dw = np ? calloc((size_t)np * nn, sizeof(double)) : NULL;
    double dep[MAXPOLES];
    int pec;
    /* materials from the cell averages (eps_inf, sigma, Debye weights); PEC from
     * the update coefficients, which include sheets and solid metal */
    for (int j = 0; j < ny; j++)
        for (int k = 0; k <= nz; k++) {
            size_t q = IDX(s, i, j0 + j, k0 + k);
            int o = j * (nz + 1) + k;
            if (s->cb[1][q] == 0) { pu[o] = 1; eu[o] = 1; continue; }
            edge_material(s, 1, i, j0 + j, k0 + k, &eu[o], &su[o], &pec);
            if (np && cem_edge_debye(s, 1, i, j0 + j, k0 + k, dep))
                for (int p = 0; p < np; p++) du[(size_t)p * nu + o] = dep[p];
        }
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k < nz; k++) {
            size_t q = IDX(s, i, j0 + j, k0 + k);
            int o = j * nz + k;
            if (s->cb[2][q] == 0) { pv[o] = 1; ev[o] = 1; continue; }
            edge_material(s, 2, i, j0 + j, k0 + k, &ev[o], &sv[o], &pec);
            if (np && cem_edge_debye(s, 2, i, j0 + j, k0 + k, dep))
                for (int p = 0; p < np; p++) dv[(size_t)p * nv + o] = dep[p];
        }
    for (int j = 0; j <= ny; j++)
        for (int k = 0; k <= nz; k++) {
            size_t q = IDX(s, i, j0 + j, k0 + k);
            int o = j * (nz + 1) + k;
            if (s->cb[0][q] == 0) { pw[o] = 1; ew[o] = 1; continue; }
            edge_material(s, 0, i, j0 + j, k0 + k, &ew[o], &sw[o], &pec);
            if (np && cem_edge_debye(s, 0, i, j0 + j, k0 + k, dep))
                for (int p = 0; p < np; p++) dw[(size_t)p * nn + o] = dep[p];
        }
    if (np) cem_mode2d_set_debye(m, np, s->pole_f, du, dv, dw);
    free(du); free(dv); free(dw);
    cem_mode2d_set_media(m, eu, ev, ew, pu, pv, pw);
    cem_mode2d_set_loss(m, su, sv, sw);
    cem_mode2d_set_pml(m, 0, s->pml_sc[1] + j0, s->pml_kc[1] + j0, s->pml_ac[1] + j0);
    cem_mode2d_set_pml(m, 1, s->pml_sn[1] + j0, s->pml_kn[1] + j0, s->pml_an[1] + j0);
    cem_mode2d_set_pml(m, 2, s->pml_sc[2] + k0, s->pml_kc[2] + k0, s->pml_ac[2] + k0);
    cem_mode2d_set_pml(m, 3, s->pml_sn[2] + k0, s->pml_kn[2] + k0, s->pml_an[2] + k0);
    free(eu); free(ev); free(ew); free(su); free(sv); free(sw); free(pu); free(pv); free(pw);
    return m;
}
