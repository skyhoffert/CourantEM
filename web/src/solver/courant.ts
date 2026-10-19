/**
 * TypeScript wrapper around the Courant EM solver core compiled to WebAssembly
 * (solver/include/courant.h). The module has no imports: it is plain C + libc
 * built with clang for wasm32, so instantiation needs no JS glue.
 *
 * Memory note: wasm memory can grow during any call that allocates, which
 * detaches existing typed-array views. Every accessor below creates a fresh view.
 */

export const Axis = { X: 0, Y: 1, Z: 2 } as const;
export type Axis = (typeof Axis)[keyof typeof Axis];
export const Face = { XLo: 0, XHi: 1, YLo: 2, YHi: 3, ZLo: 4, ZHi: 5 } as const;
export type Face = (typeof Face)[keyof typeof Face];
export const Boundary = { PEC: 0, CPML: 1 } as const;
export type Boundary = (typeof Boundary)[keyof typeof Boundary];

export interface Complex { re: number; im: number }

interface Exports {
  memory: WebAssembly.Memory;
  _initialize?: () => void;
  malloc(n: number): number;
  free(p: number): void;
  cem_last_error(): number;
  cem_create(nx: number, ny: number, nz: number, x: number, y: number, z: number): number;
  cem_destroy(s: number): void;
  cem_set_boundary(s: number, face: number, type: number, n: number): number;
  cem_material_box(s: number, x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, er: number, sigma: number): number;
  cem_pec_box(s: number, x0: number, x1: number, y0: number, y1: number, z0: number, z1: number): number;
  cem_pec_sheet(s: number, axis: number, pos: number, a0: number, a1: number, b0: number, b1: number): number;
  cem_add_lumped_port(s: number, axis: number, n1: number, n2: number, m0: number, m1: number, dir: number, R: number): number;
  cem_set_gaussian(s: number, port: number, amp: number, t0: number, tau: number): number;
  cem_set_frequencies(s: number, nf: number, f: number): number;
  cem_add_vprobe(s: number, axis: number, n1: number, n2: number, m0: number, m1: number, dir: number): number;
  cem_add_iprobe(s: number, axis: number, m: number, a0: number, a1: number, b0: number, b1: number): number;
  cem_finalize(s: number, cfl: number): number;
  cem_dt(s: number): number;
  cem_step_count(s: number): number;   /* C long is 32-bit on wasm32 */
  cem_step(s: number, n: number): void;
  cem_energy(s: number): number;
  cem_nfreq(s: number): number;
  cem_port_v(s: number, p: number): number;
  cem_port_i(s: number, p: number): number;
  cem_source_spectrum(s: number): number;
  cem_vprobe(s: number, p: number): number;
  cem_iprobe(s: number, p: number): number;
  cem_add_mode_port_args(s: number, dir: number, iSrc: number, iMon: number, j0: number, j1: number, k0: number,
    k1: number, vj: number, vk0: number, vk1: number, vsign: number, ja: number, jb: number, ka: number, kb: number,
    xRef: number, f0: number, neff0: number): number;
  cem_set_mode_excitation(s: number, port: number, t0: number, tau: number): number;
  cem_mode_port_waves(s: number, port: number, a: number, b: number, z0: number, beta: number): number;
  cem_set_debye_poles(s: number, n: number, fp: number): number;
  cem_material_box_debye(s: number, x0: number, x1: number, y0: number, y1: number, z0: number, z1: number,
    epsInf: number, sigma: number, deps: number): number;
  cem_set_conductor_loss(s: number, sigma: number, rough: number, p1: number, p2: number, fmin: number,
    fmax: number): number;
  cem_conductor_fit_error(s: number): number;
  cem_fit_djordjevic_sarkar(dk: number, df: number, fref: number, fmin: number, fmax: number, n: number,
    epsInf: number, fp: number, deps: number): number;
  cem_s2_renormalize(S: number, z1: number, z2: number, zref: number): void;
  cem_plan_array(s: number, which: number, comp: number): number;
  cem_plan_npml(s: number, face: number): number;
  cem_plan_has_pml(s: number): number;
  cem_plan_nport(s: number): number;
  cem_plan_nvprobe(s: number): number;
  cem_plan_niprobe(s: number): number;
  cem_plan_ndisp(s: number): number;
  cem_state_len(s: number, which: number, comp: number, w: number): number;
  cem_plan_nsibc(s: number): number;
  cem_plan_port(s: number, p: number, idx: number, vw: number, cs: number): number;
  cem_plan_vprobe(s: number, p: number, idx: number, w: number): number;
  cem_plan_iprobe(s: number, p: number, idx: number, w: number): number;
  cem_plan_mode_source(s: number, idx: number, w: number): number;
  cem_plan_waveforms(s: number, n0: number, n: number, lumped: number, mode: number): void;
  cem_plan_excited_lumped(s: number): number;
  cem_push_samples(s: number, n: number, ip: number, vold: number, vnew: number, vp: number): void;
  cem_field(s: number, comp: number): number;
  cem_plan_npole(s: number): number;
  cem_plan_disp(s: number, idx: number, cb: number, beta: number, k: number): number;
  cem_plan_sibc(s: number, idx: number, hidx: number, hw: number, c: number, e: number, nterm: number): number;
  cem_checkpoint_size(s: number): number;
  cem_checkpoint_save(s: number, buf: number): number;
  cem_checkpoint_load(s: number, buf: number): number;
  cem_state_array(s: number, which: number, comp: number, w: number): number;
}

/** Finalized update plan, as exported for an external time stepper. */
export interface Plan {
  nx: number; ny: number; nz: number; N: number; dt: number;
  ca: Float64Array[]; cb: Float64Array[];
  ide: Float64Array[]; idh: Float64Array[];
  be: Float64Array[]; ce: Float64Array[]; bh: Float64Array[]; ch: Float64Array[];
  d: Float64Array[]; dd: Float64Array[];
  npml: number[]; hasPml: boolean;
  ports: { idx: Float64Array; vw: Float64Array; cs: Float64Array }[];
  excitedLumped: number;
  vprobes: { idx: Float64Array; w: Float64Array }[];
  iprobes: { idx: Float64Array; w: Float64Array }[];
  modeSource: { idx: Float64Array; w: Float64Array };
  ndisp: number; nsibc: number;
  disp?: { npole: number; idx: Float64Array; cb: Float64Array; beta: Float64Array; k: Float64Array };
  sibc?: { idx: Float64Array; hidx: Float64Array; hw: Float64Array; c: Float64Array; e: Float64Array };
}

/** Debye-pole model of a dielectric: eps(w) = epsInf + sum deps_p / (1 + j f / f_p). */
export interface DebyeModel { epsInf: number; fp: number[]; deps: number[]; tanDeltaError: number }

/** Mode port on an x-directed line (see cem_mode_port_def in courant.h). */
export interface ModePortDef {
  dir: 1 | -1; iSrc: number; iMon: number;
  j0: number; j1: number; k0: number; k1: number;
  vj: number; vk0: number; vk1: number; vsign: 1 | -1;
  ja: number; jb: number; ka: number; kb: number;
  xRef: number; f0: number; neff0: number;
}

import { ThreadPool, type MtExports } from './pool.ts';

export class CourantError extends Error {}

export class Courant {
  readonly ex: Exports;
  /** Set when running the multithreaded build. */
  pool?: ThreadPool;
  constructor(instance: WebAssembly.Instance, memory?: WebAssembly.Memory) {
    // the threads build imports its (shared) memory instead of exporting it
    this.ex = (memory ? { ...instance.exports, memory } : instance.exports) as unknown as Exports;
    this.ex._initialize?.();
  }

  /**
   * Load the multithreaded build (courant-mt.wasm) with `threads` threads.
   * Needs crossOriginIsolated and a worker context; results are identical to load().
   */
  static async loadThreaded(source: BufferSource | Response | Promise<Response>, threads: number,
                            workerUrl: string | URL): Promise<Courant> {
    const bytes = source instanceof ArrayBuffer || ArrayBuffer.isView(source) ? source as BufferSource : await (await source).arrayBuffer();
    const module = await WebAssembly.compile(bytes);
    const memory = new WebAssembly.Memory({ initial: 256, maximum: 65536, shared: true });
    const instance = await WebAssembly.instantiate(module, { env: { memory }, wasi_snapshot_preview1: { sched_yield: () => 0 } });
    const c = new Courant(instance, memory);
    if (threads > 1) c.pool = await ThreadPool.create(module, memory, c.ex as unknown as MtExports, threads, workerUrl);
    return c;
  }

  static async load(source: BufferSource | Response | Promise<Response>): Promise<Courant> {
    let instance: WebAssembly.Instance;
    if (source instanceof ArrayBuffer || ArrayBuffer.isView(source)) {
      instance = (await WebAssembly.instantiate(source as BufferSource, {})).instance;
    } else {
      const resp = await source;
      const bytes = await resp.arrayBuffer();
      instance = (await WebAssembly.instantiate(bytes, {})).instance;
    }
    return new Courant(instance);
  }

  lastError(): string {
    const p = this.ex.cem_last_error();
    const mem = new Uint8Array(this.ex.memory.buffer);
    let end = p;
    while (mem[end] !== 0) end++;
    return new TextDecoder().decode(mem.slice(p, end)); // copy: decode() rejects shared memory views
  }

  /** Copy a float64 array into wasm memory; caller frees. */
  allocF64(values: ArrayLike<number>): number {
    const p = this.ex.malloc(values.length * 8);
    if (!p) throw new CourantError('wasm out of memory');
    new Float64Array(this.ex.memory.buffer, p, values.length).set(values);
    return p;
  }

  readF64(ptr: number, n: number): Float64Array {
    return new Float64Array(this.ex.memory.buffer, ptr, n).slice();
  }

  /** Fit a Djordjevic–Sarkar dielectric (Dk, Df at fref) with Debye poles for [fmin, fmax]. */
  fitDjordjevicSarkar(dk: number, df: number, fref: number, fmin: number, fmax: number, npole = 8): DebyeModel {
    const ex = this.ex;
    const pe = ex.malloc(8), pf = ex.malloc(8 * npole), pd = ex.malloc(8 * npole);
    try {
      const err = ex.cem_fit_djordjevic_sarkar(dk, df, fref, fmin, fmax, npole, pe, pf, pd);
      return { epsInf: this.readF64(pe, 1)[0], fp: Array.from(this.readF64(pf, npole)),
        deps: Array.from(this.readF64(pd, npole)), tanDeltaError: err };
    } finally {
      ex.free(pe); ex.free(pf); ex.free(pd);
    }
  }

  /** Renormalise a 2-port S matrix [S11, S12, S21, S22] from z1, z2 to zref. */
  renormalize2(S: [Complex, Complex, Complex, Complex], z1: number, z2: number, zref: number):
    [Complex, Complex, Complex, Complex] {
    const p = this.allocF64(S.flatMap((c) => [c.re, c.im]));
    try {
      this.ex.cem_s2_renormalize(p, z1, z2, zref);
      const v = this.readF64(p, 8);
      return [{ re: v[0], im: v[1] }, { re: v[2], im: v[3] }, { re: v[4], im: v[5] }, { re: v[6], im: v[7] }];
    } finally {
      this.ex.free(p);
    }
  }

  createSim(x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number>): Sim {
    const px = this.allocF64(x), py = this.allocF64(y), pz = this.allocF64(z);
    try {
      const h = this.ex.cem_create(x.length - 1, y.length - 1, z.length - 1, px, py, pz);
      if (!h) throw new CourantError(`cem_create: ${this.lastError()}`);
      return new Sim(this, h, x.length - 1, y.length - 1, z.length - 1);
    } finally {
      this.ex.free(px); this.ex.free(py); this.ex.free(pz);
    }
  }
}

export class Sim {
  private nf = 0;
  private readonly c: Courant;
  private handle: number;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  constructor(c: Courant, handle: number, nx: number, ny: number, nz: number) {
    this.c = c; this.handle = handle; this.nx = nx; this.ny = ny; this.nz = nz;
  }

  private check(rc: number, what: string): number {
    if (rc < 0) throw new CourantError(`${what}: ${this.c.lastError()}`);
    return rc;
  }

  setBoundary(face: Face, type: Boundary, ncells = 0): void {
    this.check(this.c.ex.cem_set_boundary(this.handle, face, type, ncells), 'set_boundary');
  }
  materialBox(x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, epsR: number, sigma = 0): void {
    this.check(this.c.ex.cem_material_box(this.handle, x0, x1, y0, y1, z0, z1, epsR, sigma), 'material_box');
  }
  pecBox(x0: number, x1: number, y0: number, y1: number, z0: number, z1: number): void {
    this.check(this.c.ex.cem_pec_box(this.handle, x0, x1, y0, y1, z0, z1), 'pec_box');
  }
  pecSheet(axis: Axis, pos: number, a0: number, a1: number, b0: number, b1: number): void {
    this.check(this.c.ex.cem_pec_sheet(this.handle, axis, pos, a0, a1, b0, b1), 'pec_sheet');
  }
  addLumpedPort(axis: Axis, n1: number, n2: number, m0: number, m1: number, dir: 1 | -1, R: number): number {
    return this.check(this.c.ex.cem_add_lumped_port(this.handle, axis, n1, n2, m0, m1, dir, R), 'add_lumped_port');
  }
  setGaussian(port: number, amp: number, t0: number, tau: number): void {
    this.check(this.c.ex.cem_set_gaussian(this.handle, port, amp, t0, tau), 'set_gaussian');
  }
  setFrequencies(f: ArrayLike<number>): void {
    const p = this.c.allocF64(f);
    try {
      this.check(this.c.ex.cem_set_frequencies(this.handle, f.length, p), 'set_frequencies');
      this.nf = f.length;
    } finally {
      this.c.ex.free(p);
    }
  }
  addVProbe(axis: Axis, n1: number, n2: number, m0: number, m1: number, dir: 1 | -1): number {
    return this.check(this.c.ex.cem_add_vprobe(this.handle, axis, n1, n2, m0, m1, dir), 'add_vprobe');
  }
  addIProbe(axis: Axis, m: number, a0: number, a1: number, b0: number, b1: number): number {
    return this.check(this.c.ex.cem_add_iprobe(this.handle, axis, m, a0, a1, b0, b1), 'add_iprobe');
  }
  finalize(cfl = 0.99): void {
    if (this.c.ex.cem_finalize(this.handle, cfl) !== 0) throw new CourantError(`finalize: ${this.c.lastError()}`);
    this.nf = this.c.ex.cem_nfreq(this.handle);
  }
  get dt(): number { return this.c.ex.cem_dt(this.handle); }
  get steps(): number { return Number(this.c.ex.cem_step_count(this.handle)); }
  step(n: number): void {
    if (this.c.pool) this.c.pool.step(this.handle, this.nx, n);
    else this.c.ex.cem_step(this.handle, n);
  }
  energy(): number { return this.c.ex.cem_energy(this.handle); }

  private spectrum(ptr: number): Complex[] {
    const a = this.c.readF64(ptr, 2 * this.nf);
    const out: Complex[] = [];
    for (let i = 0; i < this.nf; i++) out.push({ re: a[2 * i], im: a[2 * i + 1] });
    return out;
  }
  portV(p: number): Complex[] { return this.spectrum(this.c.ex.cem_port_v(this.handle, p)); }
  portI(p: number): Complex[] { return this.spectrum(this.c.ex.cem_port_i(this.handle, p)); }
  vprobe(p: number): Complex[] { return this.spectrum(this.c.ex.cem_vprobe(this.handle, p)); }
  iprobe(p: number): Complex[] { return this.spectrum(this.c.ex.cem_iprobe(this.handle, p)); }

  setDebyePoles(fp: number[]): void {
    const p = this.c.allocF64(fp);
    try { this.check(this.c.ex.cem_set_debye_poles(this.handle, fp.length, p), 'set_debye_poles'); }
    finally { this.c.ex.free(p); }
  }
  materialBoxDebye(x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, m: DebyeModel,
                   sigma = 0): void {
    const p = this.c.allocF64(m.deps);
    try {
      this.check(this.c.ex.cem_material_box_debye(this.handle, x0, x1, y0, y1, z0, z1, m.epsInf, sigma, p),
        'material_box_debye');
    } finally { this.c.ex.free(p); }
  }
  /** Copper loss on every metal surface. rough: 0 none, 1 Groisse (p1 = Rq m), 2 Huray (p1 = radius, p2 = density). */
  setConductorLoss(sigma: number, fmin: number, fmax: number, rough = 0, p1 = 0, p2 = 0): void {
    this.check(this.c.ex.cem_set_conductor_loss(this.handle, sigma, rough, p1, p2, fmin, fmax), 'set_conductor_loss');
  }

  addModePort(d: ModePortDef): number {
    return this.check(this.c.ex.cem_add_mode_port_args(this.handle, d.dir, d.iSrc, d.iMon, d.j0, d.j1, d.k0, d.k1,
      d.vj, d.vk0, d.vk1, d.vsign, d.ja, d.jb, d.ka, d.kb, d.xRef, d.f0, d.neff0), 'add_mode_port');
  }
  setModeExcitation(port: number, t0: number, tau: number): void {
    this.check(this.c.ex.cem_set_mode_excitation(this.handle, port, t0, tau), 'set_mode_excitation');
  }
  /** Incident/outgoing waves at the port reference plane (normalised to the port's Z0). */
  modePortWaves(port: number): { a: Complex[]; b: Complex[]; z0: Float64Array; beta: Float64Array } {
    const ex = this.c.ex, nf = this.nf;
    const pa = ex.malloc(16 * nf), pb = ex.malloc(16 * nf), pz = ex.malloc(8 * nf), pbt = ex.malloc(8 * nf);
    try {
      if (ex.cem_mode_port_waves(this.handle, port, pa, pb, pz, pbt) !== 0) {
        throw new CourantError(`mode_port_waves: ${this.c.lastError()}`);
      }
      return { a: this.spectrum(pa), b: this.spectrum(pb), z0: this.c.readF64(pz, nf), beta: this.c.readF64(pbt, nf) };
    } finally {
      ex.free(pa); ex.free(pb); ex.free(pz); ex.free(pbt);
    }
  }

  /** Raw C handle, for importer helpers that upload masks directly. */
  handleForImport(): number { return this.handle; }

  /** Export the finalized plan (coefficients, CPML, sources, probes) for another engine. */
  plan(): Plan {
    const ex = this.c.ex, h = this.handle, n = [this.nx, this.ny, this.nz];
    const N = (this.nx + 1) * (this.ny + 1) * (this.nz + 1);
    const arr = (which: number, comp: number, len: number) => this.c.readF64(ex.cem_plan_array(h, which, comp), len);
    const per = (which: number, extra: number) => [0, 1, 2].map((a) => arr(which, a, n[a] + extra));
    const list = (fn: (p: number, a: number, b: number) => number, count: number) => {
      const out: { idx: Float64Array; w: Float64Array }[] = [];
      for (let p = 0; p < count; p++) {
        const len = fn(p, 0, 0);
        const pa = ex.malloc(8 * Math.max(len, 1)), pb = ex.malloc(8 * Math.max(len, 1));
        fn(p, pa, pb);
        out.push({ idx: this.c.readF64(pa, len), w: this.c.readF64(pb, len) });
        ex.free(pa); ex.free(pb);
      }
      return out;
    };
    const ports: Plan['ports'] = [];
    for (let p = 0; p < ex.cem_plan_nport(h); p++) {
      const len = ex.cem_plan_port(h, p, 0, 0, 0);
      const a = ex.malloc(8 * len), b = ex.malloc(8 * len), c = ex.malloc(8 * len);
      ex.cem_plan_port(h, p, a, b, c);
      ports.push({ idx: this.c.readF64(a, len), vw: this.c.readF64(b, len), cs: this.c.readF64(c, len) });
      ex.free(a); ex.free(b); ex.free(c);
    }
    const ms = list((_p, a, b) => ex.cem_plan_mode_source(h, a, b), 1)[0];
    return {
      nx: this.nx, ny: this.ny, nz: this.nz, N, dt: this.dt,
      ca: [0, 1, 2].map((c) => arr(0, c, N)), cb: [0, 1, 2].map((c) => arr(1, c, N)),
      ide: per(2, 1), idh: per(3, 0), be: per(4, 1), ce: per(5, 1), bh: per(6, 0), ch: per(7, 0),
      d: per(8, 0), dd: per(9, 1),
      npml: [0, 1, 2, 3, 4, 5].map((f) => ex.cem_plan_npml(h, f)), hasPml: ex.cem_plan_has_pml(h) !== 0,
      ports, excitedLumped: ex.cem_plan_excited_lumped(h),
      vprobes: list((p, a, b) => ex.cem_plan_vprobe(h, p, a, b), ex.cem_plan_nvprobe(h)),
      iprobes: list((p, a, b) => ex.cem_plan_iprobe(h, p, a, b), ex.cem_plan_niprobe(h)),
      modeSource: ms,
      ndisp: ex.cem_plan_ndisp(h), nsibc: ex.cem_plan_nsibc(h),
      disp: this.planDisp(), sibc: this.planSibc(),
    };
  }

  private planDisp(): Plan['disp'] {
    const ex = this.c.ex, h = this.handle;
    const nd = ex.cem_plan_ndisp(h), np = ex.cem_plan_npole(h);
    if (!nd) return undefined;
    const a = ex.malloc(8 * nd), b = ex.malloc(8 * nd), c = ex.malloc(8 * nd * np), d = ex.malloc(8 * np);
    ex.cem_plan_disp(h, a, b, c, d);
    const out = { npole: np, idx: this.c.readF64(a, nd), cb: this.c.readF64(b, nd), beta: this.c.readF64(c, nd * np),
      k: this.c.readF64(d, np) };
    ex.free(a); ex.free(b); ex.free(c); ex.free(d);
    return out;
  }

  private planSibc(): Plan['sibc'] {
    const ex = this.c.ex, h = this.handle;
    const ns = ex.cem_plan_nsibc(h);
    if (!ns) return undefined;
    const a = ex.malloc(8 * ns), b = ex.malloc(32 * ns), c = ex.malloc(32 * ns), d = ex.malloc(8 * 16), e = ex.malloc(8 * 16),
      nt = ex.malloc(4);
    ex.cem_plan_sibc(h, a, b, c, d, e, nt);
    const K = new Int32Array(ex.memory.buffer, nt, 1)[0];
    const out = { idx: this.c.readF64(a, ns), hidx: this.c.readF64(b, 4 * ns), hw: this.c.readF64(c, 4 * ns),
      c: this.c.readF64(d, K), e: this.c.readF64(e, K) };
    ex.free(a); ex.free(b); ex.free(c); ex.free(d); ex.free(e); ex.free(nt);
    return out;
  }

  /** Source waveforms (lumped, mode) at H times for steps n0 .. n0+n-1. */
  waveforms(n0: number, n: number): { lumped: Float64Array; mode: Float64Array } {
    const ex = this.c.ex;
    const a = ex.malloc(8 * n), b = ex.malloc(8 * n);
    ex.cem_plan_waveforms(this.handle, n0, n, a, b);
    const out = { lumped: this.c.readF64(a, n), mode: this.c.readF64(b, n) };
    ex.free(a); ex.free(b);
    return out;
  }

  /** Hand samples recorded by an external engine back to the C DFTs (advances the step count). */
  pushSamples(n: number, ip: Float64Array, vold: Float64Array, vnew: Float64Array, vp: Float64Array): void {
    const c = this.c;
    const one = (a: Float64Array) => c.allocF64(a.length ? a : [0]);
    const p1 = one(ip), p2 = one(vold), p3 = one(vnew), p4 = one(vp);
    try { c.ex.cem_push_samples(this.handle, n, p1, p2, p3, p4); }
    finally { c.ex.free(p1); c.ex.free(p2); c.ex.free(p3); c.ex.free(p4); }
  }

  /** Snapshot of all time-dependent state (fields, CPML, Debye, SIBC, DFTs, step). */
  checkpoint(): Float64Array {
    const ex = this.c.ex, n = ex.cem_checkpoint_size(this.handle);
    const p = ex.malloc(8 * n);
    try {
      if (ex.cem_checkpoint_save(this.handle, p)) throw new CourantError(`checkpoint: ${this.c.lastError()}`);
      return this.c.readF64(p, n);
    } finally { ex.free(p); }
  }
  /** Restore a checkpoint taken from an identically built simulation. */
  restore(buf: Float64Array): void {
    const p = this.c.allocF64(buf);
    try {
      if (this.c.ex.cem_checkpoint_load(this.handle, p)) throw new CourantError(`restore: ${this.c.lastError()}`);
    } finally { this.c.ex.free(p); }
  }
  /** Write a field component (0..2 E, 3..5 H) from an external engine. */
  setField(comp: number, data: ArrayLike<number>): void {
    const N = (this.nx + 1) * (this.ny + 1) * (this.nz + 1);
    new Float64Array(this.c.ex.memory.buffer, this.c.ex.cem_field(this.handle, comp), N).set(data);
  }
  /** Length of a state array (see cem_state_len). */
  stateLen(which: number, comp = 0, w = 0): number { return this.c.ex.cem_state_len(this.handle, which, comp, w); }
  /** View of a writable state array (see cem_state_array); null if absent. */
  stateArray(which: number, comp: number, w: number, len: number): Float64Array | null {
    const p = this.c.ex.cem_state_array(this.handle, which, comp, w);
    return p ? new Float64Array(this.c.ex.memory.buffer, p, len) : null;
  }

  /** Copy of a field component (0..2 E, 3..5 H). */
  field(comp: number): Float64Array {
    const N = (this.nx + 1) * (this.ny + 1) * (this.nz + 1);
    return this.c.readF64(this.c.ex.cem_field(this.handle, comp), N);
  }

  destroy(): void {
    if (this.handle) this.c.ex.cem_destroy(this.handle);
    this.handle = 0;
  }
}

/** S-parameters from port spectra with port 1 excited (reference impedance R). */
export function sParams(v1: Complex[], i1: Complex[], v2: Complex[], i2: Complex[], R: number) {
  const sq = Math.sqrt(R);
  const s11: Complex[] = [], s21: Complex[] = [];
  for (let k = 0; k < v1.length; k++) {
    const a1 = { re: (v1[k].re + R * i1[k].re) / (2 * sq), im: (v1[k].im + R * i1[k].im) / (2 * sq) };
    const b1 = { re: (v1[k].re - R * i1[k].re) / (2 * sq), im: (v1[k].im - R * i1[k].im) / (2 * sq) };
    const b2 = { re: (v2[k].re - R * i2[k].re) / (2 * sq), im: (v2[k].im - R * i2[k].im) / (2 * sq) };
    const d = a1.re * a1.re + a1.im * a1.im;
    s11.push({ re: (b1.re * a1.re + b1.im * a1.im) / d, im: (b1.im * a1.re - b1.re * a1.im) / d });
    s21.push({ re: (b2.re * a1.re + b2.im * a1.im) / d, im: (b2.im * a1.re - b2.re * a1.im) / d });
  }
  return { s11, s21 };
}
