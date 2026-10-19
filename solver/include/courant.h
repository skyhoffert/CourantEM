/*
 * Courant EM - FDTD solver core (double-precision reference backend).
 *
 * Flat C API so the same core builds natively (tests, CI) and to WebAssembly.
 * All quantities are SI. The grid is a non-uniform Cartesian Yee grid given by
 * its node coordinates on each axis.
 *
 * Yee layout (cell (i,j,k) spans [x_i, x_i+1] x [y_j, y_j+1] x [z_k, z_k+1]):
 *   Ex(i+1/2, j, k)   Ey(i, j+1/2, k)   Ez(i, j, k+1/2)
 *   Hx(i, j+1/2, k+1/2)   Hy(i+1/2, j, k+1/2)   Hz(i+1/2, j+1/2, k)
 * Every component is stored in an (nx+1)(ny+1)(nz+1) array indexed by its
 * integer parts, k fastest.
 *
 * Lifecycle:
 *   cem_create -> geometry/material/ports/probes -> cem_finalize -> cem_run/cem_step
 *   -> read spectra -> cem_destroy
 */
#ifndef COURANT_H
#define COURANT_H

#include <stddef.h>

#if defined(__wasm__)
#define CEM_API __attribute__((used, visibility("default")))
#else
#define CEM_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define CEM_C0 299792458.0
#define CEM_MU0 1.25663706212e-6
#define CEM_EPS0 8.8541878128e-12

enum { CEM_X = 0, CEM_Y = 1, CEM_Z = 2 };
/* faces: 2*axis + side (side 0 = low, 1 = high) */
enum { CEM_XLO = 0, CEM_XHI, CEM_YLO, CEM_YHI, CEM_ZLO, CEM_ZHI };
enum { CEM_BC_PEC = 0, CEM_BC_CPML = 1 };

typedef struct cem_sim cem_sim;

/* ---- creation -------------------------------------------------------- */

/* x, y, z: node coordinates, lengths nx+1, ny+1, nz+1, strictly increasing. */
CEM_API cem_sim *cem_create(int nx, int ny, int nz, const double *x, const double *y, const double *z);
CEM_API void cem_destroy(cem_sim *s);
CEM_API const char *cem_last_error(void);

/* Boundary per face. For CPML, ncells cells adjacent to that face form the layer
 * (they are part of the grid). PEC is the default on every face. */
CEM_API int cem_set_boundary(cem_sim *s, int face, int type, int ncells);

/* ---- geometry and materials (cell based) ------------------------------ */

/* Fill every cell whose centre lies inside the box with eps_r / sigma (S/m). */
CEM_API int cem_material_box(cem_sim *s, double x0, double x1, double y0, double y1, double z0,
                     double z1, double eps_r, double sigma);
/* Debye dispersion. Set the global pole frequencies (Hz, tau = 1/(2 pi f_p))
 * once, before any dispersive material. Then cells in the box get
 *   eps(w) = eps_inf + sum_p deps[p] / (1 + j w tau_p) - j sigma / (w eps0). */
CEM_API int cem_set_debye_poles(cem_sim *s, int npole, const double *f_p);
CEM_API int cem_material_box_debye(cem_sim *s, double x0, double x1, double y0, double y1, double z0, double z1,
                                   double eps_inf, double sigma, const double *deps);
/* Copper loss: every metal surface (PEC cells, sheets, PEC outer faces) gets a
 * surface-impedance boundary Zs = (1+j) sqrt(w mu0 / 2 sigma) * K(w), fitted
 * over [fmin, fmax]. rough: 0 none; 1 Groisse with p1 = Rq (m);
 * 2 Huray single sphere size with p1 = radius (m), p2 = spheres per m^2 of flat area.
 * Edges inside CPML layers and lumped-port edges stay lossless. */
CEM_API int cem_set_conductor_loss(cem_sim *s, double sigma, int rough, double p1, double p2, double fmin,
                                   double fmax);
CEM_API double cem_conductor_fit_error(const cem_sim *s);
CEM_API void cem_conductor_zs(const cem_sim *s, double f, double *re, double *im);
/* Fit a Djordjevic-Sarkar dielectric (Dk, Df at fref) with npole Debye poles for
 * the band [fmin, fmax]. Outputs eps_inf, pole frequencies f_p and weights deps
 * (arrays of npole). Returns the max relative loss-tangent error over the band. */
CEM_API double cem_fit_djordjevic_sarkar(double dk, double df, double fref, double fmin, double fmax, int npole,
                                         double *eps_inf, double *f_p, double *deps);
/* Mark every cell whose centre lies inside the box as solid PEC. */
CEM_API int cem_pec_box(cem_sim *s, double x0, double x1, double y0, double y1, double z0, double z1);
/* Zero-thickness PEC sheet on the grid plane normal to `axis` at coordinate `pos`
 * (snapped to the nearest node plane), covering [a0,a1] x [b0,b1] in the two
 * remaining axes in cyclic order (x: y,z; y: z,x; z: x,y). Every tangential edge
 * lying fully inside the rectangle becomes PEC. */
CEM_API int cem_pec_sheet(cem_sim *s, int axis, double pos, double a0, double a1, double b0,
                  double b1);

/* Rasterised geometry. PEC edges on the node plane m normal to `axis`:
 * mask_u covers edges along the first remaining axis (cyclic), n[ua] x (n[va]+1);
 * mask_v covers edges along the second, (n[ua]+1) x n[va]. Either may be NULL. */
CEM_API int cem_pec_edges_plane(cem_sim *s, int axis, int m, const unsigned char *mask_u,
                                const unsigned char *mask_v);
/* Solid PEC cells from a mask over a cell box (k fastest). */
CEM_API int cem_pec_cells(cem_sim *s, int i0, int j0, int k0, int ni, int nj, int nk, const unsigned char *mask);

/* ---- ports and excitation --------------------------------------------- */

/* Lumped port on a column of E edges along `axis` at node (n1, n2) of the other
 * two axes (cyclic order), spanning nodes [m0, m1) along `axis`.
 * dir = +1: port voltage is the line integral of E in +axis; -1: in -axis.
 * Port voltage V and port current I (flowing out of the port into the structure,
 * against the internal source branch) follow the usual circuit convention so that
 * a = (V + Z I) / (2 sqrt Z), b = (V - Z I) / (2 sqrt Z) with Z = R.
 * Returns port index or -1. */
CEM_API int cem_add_lumped_port(cem_sim *s, int axis, int n1, int n2, int m0, int m1, int dir,
                        double R);
/* Gaussian pulse Vs(t) = amp * exp(-((t - t0)/tau)^2) on the given port. Only one
 * port is excited per run; the others are passive terminations. */
CEM_API int cem_set_gaussian(cem_sim *s, int port, double amp, double t0, double tau);

/* ---- frequency-domain probes ------------------------------------------ */

CEM_API int cem_set_frequencies(cem_sim *s, int nf, const double *f);
/* Voltage probe: line integral of E along `axis` at node (n1, n2), nodes [m0, m1),
 * times dir (+1/-1). Sampled at E times. Returns probe index. */
CEM_API int cem_add_vprobe(cem_sim *s, int axis, int n1, int n2, int m0, int m1, int dir);
/* Current probe: circulation of H around a dual-grid rectangle in the plane
 * normal to `axis` through the half node m + 1/2. The rectangle runs over dual
 * nodes (a0+1/2 .. a1+1/2) x (b0+1/2 .. b1+1/2) in the remaining axes (cyclic
 * order). Result is the current in +axis. Sampled at H times. */
CEM_API int cem_add_iprobe(cem_sim *s, int axis, int m, int a0, int a1, int b0, int b1);

/* ---- running ---------------------------------------------------------- */

/* Compute update coefficients; must be called after all geometry/ports/probes. */
CEM_API int cem_finalize(cem_sim *s, double courant_number);
CEM_API double cem_dt(const cem_sim *s);
CEM_API long cem_step_count(const cem_sim *s);
/* Advance n time steps. */
CEM_API void cem_step(cem_sim *s, long n);
/* Run until total field energy has fallen below `decay` (e.g. 1e-5) of its peak,
 * checking every `check_every` steps, after the source has finished, or until
 * max_steps. Returns steps taken. */
CEM_API long cem_run(cem_sim *s, long max_steps, double decay, int check_every);
/* Electric + magnetic energy (J), computed on demand. */
CEM_API double cem_energy(const cem_sim *s);

/* ---- results (interleaved re/im arrays of length 2*nf) ------------------ */

CEM_API int cem_nfreq(const cem_sim *s);
CEM_API const double *cem_port_v(const cem_sim *s, int port);
CEM_API const double *cem_port_i(const cem_sim *s, int port);
CEM_API const double *cem_source_spectrum(const cem_sim *s); /* DFT of the open-circuit source */
CEM_API const double *cem_vprobe(const cem_sim *s, int probe);
CEM_API const double *cem_iprobe(const cem_sim *s, int probe);
CEM_API double cem_port_impedance(const cem_sim *s, int port);
CEM_API int cem_excited_port(const cem_sim *s);

/* Raw field access for tests / visualisation. comp: 0..2 = Ex,Ey,Ez, 3..5 = Hx,Hy,Hz */
CEM_API double *cem_field(cem_sim *s, int comp);
CEM_API size_t cem_index(const cem_sim *s, int i, int j, int k);
CEM_API void cem_dims(const cem_sim *s, int *nx, int *ny, int *nz);


/* ---- 2D eigenmode solver (port cross-sections) -------------------------- */

typedef struct cem_mode2d cem_mode2d;
/* Cross-section grid: y nodes (ny+1), z nodes (nz+1). Outer boundary is PEC. */
CEM_API cem_mode2d *cem_mode2d_create(int ny, int nz, const double *y, const double *z);
CEM_API void cem_mode2d_destroy(cem_mode2d *m);
/* Relative permittivity and PEC flags per unknown:
 *   E_u (Ey) slots j<ny, k<=nz at index j*(nz+1)+k
 *   E_v (Ez) slots j<=ny, k<nz at index j*nz+k
 *   E_w (Ex) nodes j<=ny, k<=nz at index j*(nz+1)+k */
CEM_API void cem_mode2d_set_media(cem_mode2d *m, const double *epsr_u, const double *epsr_v, const double *epsr_w,
                                  const unsigned char *pec_u, const unsigned char *pec_v, const unsigned char *pec_w);
/* Find the mode with beta^2 nearest (k0 neff_guess)^2. `start` (optional) is an
 * initial transverse E field [E_u | E_v]. Returns 0 on convergence. */
CEM_API int cem_mode2d_solve(cem_mode2d *m, double f, double neff_guess, const double *start, int maxit, double tol);
/* Optional Debye dispersion per unknown: deps arrays are npole blocks, each in
 * the set_media layout for that unknown type. f_p in Hz. */
CEM_API void cem_mode2d_set_debye(cem_mode2d *m, int npole, const double *f_p, const double *deps_u,
                                  const double *deps_v, const double *deps_w);
/* Optional: conductivity (S/m) per unknown, same layouts as set_media. */
CEM_API void cem_mode2d_set_loss(cem_mode2d *m, const double *sig_u, const double *sig_v, const double *sig_w);
/* Optional stretched-coordinate PML s = kappa + sigma/(alpha + j w eps0).
 * which: 0 = y cells (ny), 1 = y nodes (ny+1), 2 = z cells (nz), 3 = z nodes (nz+1). */
CEM_API void cem_mode2d_set_pml(cem_mode2d *m, int which, const double *sigma, const double *kappa,
                                const double *alpha);
/* Propagation constant gamma = alpha + j beta of the last solution (fields ~ exp(-gamma x)). */
CEM_API double cem_mode2d_beta(const cem_mode2d *m);
CEM_API double cem_mode2d_alpha(const cem_mode2d *m);
/* Field arrays are complex, stored as interleaved re/im pairs. */
CEM_API const double *cem_mode2d_eu(const cem_mode2d *m);
CEM_API const double *cem_mode2d_ev(const cem_mode2d *m);
CEM_API const double *cem_mode2d_hu(const cem_mode2d *m);
CEM_API const double *cem_mode2d_hv(const cem_mode2d *m);
CEM_API double cem_mode2d_power(const cem_mode2d *m);
/* Conductor attenuation (Np/m) of the last solved mode for surface resistance Rs,
 * by the perturbation method over all metal surfaces facing a field region. */
CEM_API double cem_mode2d_conductor_alpha(const cem_mode2d *m, double Rs);
/* V and I are complex: each points to 2 doubles (re, im). */
CEM_API void cem_mode2d_vi(const cem_mode2d *m, int j0, int k0, int k1, int ja, int jb, int ka, int kb, double *V,
                           double *I);
/* ---- mode ports (x-directed lines) --------------------------------------- */

typedef struct {
    int dir;                 /* +1 if the device under test lies on the +x side */
    int i_src;               /* x node plane of the modal source (line continues into CPML behind it) */
    int i_mon;               /* x node plane of the monitor (V at i_mon, i_mon+1; I at cell i_mon) */
    int j0, j1, k0, k1;      /* cross-section window in y/z nodes, PEC on its edge */
    int vj, vk0, vk1, vsign; /* voltage path along z at y node vj: V = vsign * sum Ez dz */
    int ja, jb, ka, kb;      /* current loop: dual nodes (ja+1/2..jb+1/2) x (ka+1/2..kb+1/2) */
    double x_ref;            /* reference plane (m) */
    double f0;               /* frequency of the source mode profile (Hz) */
    double neff0;            /* effective index guess for the mode search */
} cem_mode_port_def;

CEM_API int cem_add_mode_port(cem_sim *s, const cem_mode_port_def *d);
/* Same as cem_add_mode_port with the fields passed as arguments (for JS callers). */
CEM_API int cem_add_mode_port_args(cem_sim *s, int dir, int i_src, int i_mon, int j0, int j1, int k0, int k1, int vj,
                                   int vk0, int vk1, int vsign, int ja, int jb, int ka, int kb, double x_ref,
                                   double f0, double neff0);
/* Excite a mode port with a differentiated Gaussian (zero DC). */
CEM_API int cem_set_mode_excitation(cem_sim *s, int port, double t0, double tau);
/* Incident (a) and outgoing (b) waves at the reference plane, normalised to the
 * port's own mode impedance Z0(f); interleaved re/im, 2*nf each. Also returns
 * Z0(f) and beta(f) (real, lossless) if the pointers are non-NULL. Solves the
 * port mode at every output frequency on first call. Returns 0 on success. */
CEM_API int cem_mode_port_waves(cem_sim *s, int port, double *a, double *b, double *z0, double *beta);
/* Attenuation constant alpha(f) in Np/m of the port mode (after cem_mode_port_waves). */
CEM_API int cem_mode_port_alpha(const cem_sim *s, int port, double *alpha);
/* The conductor part of alpha (perturbation estimate; zero without conductor loss). */
CEM_API int cem_mode_port_alpha_c(const cem_sim *s, int port, double *alpha_c);

/* Renormalise a 2-port S matrix (interleaved re/im, order S11 S12 S21 S22) from
 * port reference impedances z1, z2 (real) to zref, in place. */
CEM_API void cem_s2_renormalize(double *S, double z1, double z2, double zref);

/* ---- update plan for external time steppers (WebGPU backend) ------------- */
CEM_API const double *cem_plan_array(const cem_sim *s, int which, int comp);
CEM_API int cem_plan_npml(const cem_sim *s, int face);
CEM_API int cem_plan_has_pml(const cem_sim *s);
CEM_API int cem_plan_nport(const cem_sim *s);
CEM_API int cem_plan_nvprobe(const cem_sim *s);
CEM_API int cem_plan_niprobe(const cem_sim *s);
CEM_API int cem_plan_ndisp(const cem_sim *s);
CEM_API int cem_plan_nsibc(const cem_sim *s);
/* index arrays are returned as doubles (exact up to 2^53) to keep one element type */
CEM_API int cem_plan_port(const cem_sim *s, int p, double *idx, double *vw, double *cs);
CEM_API int cem_plan_vprobe(const cem_sim *s, int p, double *idx, double *w);
CEM_API int cem_plan_iprobe(const cem_sim *s, int p, double *idx, double *w);
CEM_API int cem_plan_mode_source(const cem_sim *s, double *idx, double *w);
CEM_API void cem_plan_waveforms(const cem_sim *s, long n0, int n, double *lumped, double *mode);
CEM_API int cem_plan_excited_lumped(const cem_sim *s);
CEM_API int cem_plan_npole(const cem_sim *s);
CEM_API int cem_plan_disp(const cem_sim *s, double *idx, double *cb, double *beta, double *k);
CEM_API int cem_plan_sibc(const cem_sim *s, double *idx, double *hidx, double *hw, double *c, double *e, int *nterm);
CEM_API void cem_push_samples(cem_sim *s, int n, const double *ip, const double *vold, const double *vnew,
                              const double *vp);

/* ---- checkpoints: save and resume the time-dependent state --------------- */
/* Number of doubles a checkpoint takes (0 if not finalized). */
CEM_API size_t cem_checkpoint_size(const cem_sim *s);
CEM_API int cem_checkpoint_save(const cem_sim *s, double *buf);
/* Restore into an identically built, finalized sim; refuses mismatched layouts. */
CEM_API int cem_checkpoint_load(cem_sim *s, const double *buf);

CEM_API double *cem_state_array(cem_sim *s, int which, int comp, int w);
CEM_API size_t cem_state_len(const cem_sim *s, int which, int comp, int w);

/* Time step in phases for an external thread pool (see sim.c):
 * pre; H over i-slabs [i0,i1) covering 0..nx+1; mid; E over slabs covering 0..nx; post. */
CEM_API void cem_mt_pre(cem_sim *s);
/* Solver version (CalVer from the repo's VERSION file, "dev" if built without it). */
CEM_API const char *cem_version(void);
CEM_API void cem_mt_H(cem_sim *s, int i0, int i1);
CEM_API void cem_mt_mid(cem_sim *s);
CEM_API void cem_mt_E(cem_sim *s, int i0, int i1);
CEM_API void cem_mt_post(cem_sim *s);
CEM_API void cem_set_step(cem_sim *s, long step);

/* Build the cross-section of a finalized 3D sim on the x node plane i, window of
 * y nodes [j0, j1] and z nodes [k0, k1] (PEC on the window edge). Materials,
 * losses and the 3D CPML profiles inside the window are carried over. */
CEM_API cem_mode2d *cem_cross_section_x(const cem_sim *s, int i, int j0, int j1, int k0, int k1);

#ifdef __cplusplus
}
#endif
#endif
