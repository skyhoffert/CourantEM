/**
 * WebGPU time stepper (f32). The C core (wasm, f64) builds and finalizes the
 * simulation; this engine uploads the plan, advances the fields on the GPU, and
 * streams probe samples back into the C DFTs with cem_push_samples, so all
 * post-processing is shared with the reference path.
 *
 * Every kernel of the C core runs here too: Yee updates, CPML, lumped and modal
 * sources, probes, Debye dispersion and surface-impedance copper. Buffers are
 * single per quantity for now; the plan's multi-buffer split for >128 MB
 * bindings comes with large problems.
 */
import type { Plan, Sim } from '../solver/courant.ts';
import { dispWGSL, fdtdWGSL, ioWGSL, sibcWGSL } from './fdtd.wgsl.ts';

const WG = 256;
// WebGPU flag values (spec constants), so the code does not depend on lib typings
const BU = { MAP_READ: 1, COPY_SRC: 4, COPY_DST: 8, UNIFORM: 64, STORAGE: 128 } as const;
const MAP_READ_MODE = 1;

function f32(a: ArrayLike<number>): Float32Array { return Float32Array.from(a); }

export class GpuEngine {
  readonly device: GPUDevice;
  private sim: Sim;
  private plan: Plan;
  private chunk: number;
  private pipes: Record<string, GPUComputePipeline> = {};
  private bgIO: GPUBindGroup[] = [];
  private bufs: Record<string, GPUBuffer> = {};
  private nprobe = 0;
  private nsrc = 0;
  private groups = { ip: [0, 0], vold: [0, 0], vnew: [0, 0], vp: [0, 0] };
  private slabCount = { E: [0, 0, 0], H: [0, 0, 0] };
  private psiSlots: { off: number; len: number; which: number; c: number; w: number }[] = [];
  private stepsDone = 0;
  /** Largest |probe sample| in the last chunk (run control without a field read-back). */
  lastChunkPeak = 0;

  private constructor(device: GPUDevice, sim: Sim, plan: Plan, chunk: number) {
    this.device = device; this.sim = sim; this.plan = plan; this.chunk = chunk;
  }

  static async create(sim: Sim, opts: { device?: GPUDevice; chunk?: number } = {}): Promise<GpuEngine> {
    let device = opts.device;
    if (!device) {
      if (!navigator.gpu) throw new Error('WebGPU is not available in this browser');
      const adapter = await navigator.gpu.requestAdapter();
      if (!adapter) throw new Error('no WebGPU adapter');
      device = await adapter.requestDevice({
        requiredLimits: {
          maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize,
          maxBufferSize: adapter.limits.maxBufferSize,
        },
      });
    }
    // the largest bindings are E, H, CA, CB (3 N floats each)
    const N = (sim.nx + 1) * (sim.ny + 1) * (sim.nz + 1);
    const need = 4 * 3 * N, lim = Math.min(device.limits.maxStorageBufferBindingSize, device.limits.maxBufferSize);
    if (need > lim) {
      throw new Error(`this mesh needs ${(need / 2 ** 20).toFixed(0)} MB GPU buffers but this GPU allows ${(lim / 2 ** 20).toFixed(0)} MB; ` +
        `use the CPU backend or a coarser mesh (max ~${(lim / 12 / 1e6).toFixed(1)} M points)`);
    }
    const plan = sim.plan();
    const eng = new GpuEngine(device, sim, plan, opts.chunk ?? 256);
    eng.stepsDone = sim.steps; // resumes continue the source waveforms where they left off
    eng.build();
    eng.buildExtras();
    return eng;
  }

  private buffer(name: string, data: ArrayBufferView | number, usage: number): GPUBuffer {
    const size = typeof data === 'number' ? data : data.byteLength;
    const b = this.device.createBuffer({ size: Math.max(16, Math.ceil(size / 4) * 4), usage, mappedAtCreation: typeof data !== 'number' });
    if (typeof data !== 'number') {
      const dst = new Uint8Array(b.getMappedRange());
      dst.set(new Uint8Array(data.buffer, data.byteOffset, data.byteLength));
      b.unmap();
    }
    this.bufs[name] = b;
    return b;
  }

  private build(): void {
    const p = this.plan, dev = this.device, N = p.N;
    const S = BU.STORAGE, CD = BU.COPY_DST, CS = BU.COPY_SRC;
    // fields start from the C state (normally zero)
    const E = new Float32Array(3 * N), H = new Float32Array(3 * N);
    for (let c = 0; c < 3; c++) { E.set(f32(this.sim.field(c)), c * N); H.set(f32(this.sim.field(3 + c)), c * N); }
    this.buffer('E', E, S | CD | CS);
    this.buffer('H', H, S | CD | CS);
    const CA = new Float32Array(3 * N), CB = new Float32Array(3 * N);
    for (let c = 0; c < 3; c++) { CA.set(f32(p.ca[c]), c * N); CB.set(f32(p.cb[c]), c * N); }
    this.buffer('CA', CA, S);
    this.buffer('CB', CB, S);
    // packed 1D arrays
    const sp: number[] = [];
    const off: Record<string, number[]> = {};
    const push = (name: string, arrs: Float64Array[] | number[][]) => {
      off[name] = [];
      for (const a of arrs) { off[name].push(sp.length); for (const v of a) sp.push(v); }
      off[name].push(0);
    };
    push('ide', p.ide); push('idh', p.idh); push('be', p.be); push('ce', p.ce); push('bh', p.bh); push('ch', p.ch);
    push('iddr', p.dd.map((a) => Array.from(a, (v) => (v > 0 ? 1 / v : 0))));
    push('idr', p.d.map((a) => Array.from(a, (v) => 1 / v)));
    this.buffer('SP', f32(sp), S);
    if (p.hasPml) {
      // CPML auxiliaries start from the C state, so a run can resume from a checkpoint
      // compact layout (as C): slot 2c+w psiE[c][w], 6+2c+w psiH[c][w]
      let tot = 0;
      for (let slot = 0; slot < 12; slot++) {
        const which = slot < 6 ? 0 : 1, c = Math.floor((slot % 6) / 2), w = slot % 2;
        const len = this.sim.stateLen(which, c, w);
        this.psiSlots.push({ off: tot, len, which, c, w });
        tot += len;
      }
      const psi = new Float32Array(Math.max(tot, 4));
      for (const sl of this.psiSlots) {
        const a = this.sim.stateArray(sl.which, sl.c, sl.w, sl.len);
        if (a) psi.set(f32(a), sl.off);
      }
      this.buffer('PSI', psi, S | CS);
    } else {
      this.buffer('PSI', 16, S);
    }
    // uniform block
    const u = new Uint32Array(76);
    const uf = new Float32Array(u.buffer);
    const n = [p.nx, p.ny, p.nz];
    u[0] = p.nx; u[1] = p.ny; u[2] = p.nz; u[3] = N;
    u[4] = (p.ny + 1) * (p.nz + 1); u[5] = p.nz + 1;
    uf[6] = p.dt / 1.25663706212e-6;
    const names = ['ide', 'idh', 'be', 'ce', 'bh', 'ch', 'iddr', 'idr'];
    names.forEach((nm, i) => { for (let a = 0; a < 4; a++) u[8 + 4 * i + a] = off[nm][a]; });
    for (let a = 0; a < 3; a++) {
      const lo = p.npml[2 * a], hi = p.npml[2 * a + 1];
      // E nodes: 1..lo and n-hi..n-1 ; H cells: 0..lo-1 and n-hi..n-1
      const e = [1, lo, n[a] - hi, hi], hh = [0, lo, n[a] - hi, hi];
      for (let q = 0; q < 4; q++) { u[40 + 4 * a + q] = e[q]; u[52 + 4 * a + q] = hh[q]; }
      this.slabCount.E[a] = lo + hi; this.slabCount.H[a] = lo + hi;
    }
    this.psiSlots.forEach((sl, n) => { u[64 + n] = sl.off; });
    this.buffer('U', u, BU.UNIFORM | CD);

    const mod = dev.createShaderModule({ code: fdtdWGSL });
    for (const e of ['updateH', 'updateE', 'pmlHx', 'pmlHy', 'pmlHz', 'pmlEx', 'pmlEy', 'pmlEz']) {
      this.pipes[e] = dev.createComputePipeline({ layout: 'auto', compute: { module: mod, entryPoint: e } });
    }
    // 'auto' layouts keep only the bindings an entry point uses, so each pipeline
    // gets its own bind group
    const entries = ['U', 'E', 'H', 'CA', 'CB', 'SP', 'PSI'].map((nm, i) => ({ binding: i, resource: { buffer: this.bufs[nm] } }));
    this.mainGroups = {};
    for (const [nm, pipe] of Object.entries(this.pipes)) {
      const lay = pipe.getBindGroupLayout(0);
      const used = entries.filter((e) => this.usesBinding(nm, e.binding));
      this.mainGroups[nm] = dev.createBindGroup({ layout: lay, entries: used });
    }

    // ---- sources and probes ----
    const src: number[] = [];
    const pushSrc = (idx: number, coef: number, type: number) => {
      src.push(idx, new Uint32Array(new Float32Array([coef]).buffer)[0], type, 0);
    };
    if (p.excitedLumped >= 0) {
      const pt = p.ports[p.excitedLumped];
      pt.idx.forEach((ix, e) => pushSrc(ix, pt.cs[e], 0));
    }
    p.modeSource.idx.forEach((ix, e) => pushSrc(ix, p.modeSource.w[e], 1));
    this.nsrc = src.length / 4;
    this.buffer('SRC', new Uint32Array(src.length ? src : [0, 0, 0, 0]), S);
    const probes: number[] = [], pent: number[] = [];
    const addProbe = (idx: ArrayLike<number>, w: ArrayLike<number>, field: number) => {
      probes.push(pent.length / 2, idx.length, field, 0);
      for (let e = 0; e < idx.length; e++) pent.push(idx[e], new Uint32Array(new Float32Array([w[e]]).buffer)[0]);
    };
    let k = 0;
    this.groups.ip = [k, k + p.iprobes.length]; for (const q of p.iprobes) addProbe(q.idx, q.w, 1); k += p.iprobes.length;
    this.groups.vold = [k, k + p.ports.length]; for (const q of p.ports) addProbe(q.idx, q.vw, 0); k += p.ports.length;
    this.groups.vnew = [k, k + p.ports.length]; for (const q of p.ports) addProbe(q.idx, q.vw, 0); k += p.ports.length;
    this.groups.vp = [k, k + p.vprobes.length]; for (const q of p.vprobes) addProbe(q.idx, q.w, 0); k += p.vprobes.length;
    this.nprobe = k;
    this.buffer('PROBE', new Uint32Array(probes.length ? probes : [0, 0, 0, 0]), S);
    this.buffer('PENT', new Uint32Array(pent.length ? pent : [0, 0]), S);
    this.buffer('WAVE', 8 * this.chunk, S | CD);
    this.buffer('STATE', 4 * (1 + this.chunk * Math.max(1, this.nprobe)), S | CD | CS);
    this.buffer('READ', 4 * (1 + this.chunk * Math.max(1, this.nprobe)), BU.MAP_READ | CD);
    const io = dev.createShaderModule({ code: ioWGSL });
    for (const e of ['inject', 'probes', 'tick']) {
      this.pipes[e] = dev.createComputePipeline({ layout: 'auto', compute: { module: io, entryPoint: e } });
    }
    const ioU = (p0: number, p1: number) => this.buffer(`IOU${p0}_${p1}`, new Uint32Array([this.nsrc, this.nprobe, p0, p1]),
      BU.UNIFORM);
    const ioEntries = (u: GPUBuffer) => [u, this.bufs.E, this.bufs.H, this.bufs.SRC, this.bufs.WAVE, this.bufs.PROBE,
      this.bufs.PENT, this.bufs.STATE].map((b, i) => ({ binding: i, resource: { buffer: b } }));
    const mk = (pipe: string, u: GPUBuffer) => {
      const lay = this.pipes[pipe].getBindGroupLayout(0);
      const ent = ioEntries(u).filter((e) => this.usesIO(pipe, e.binding));
      return dev.createBindGroup({ layout: lay, entries: ent });
    };
    const uAll = ioU(0, 0);
    this.ioGroups = {
      inject: mk('inject', uAll),
      tick: mk('tick', uAll),
      ip: mk('probes', ioU(this.groups.ip[0], this.groups.ip[1])),
      vold: mk('probes', ioU(this.groups.vold[0], this.groups.vold[1])),
      vnew: mk('probes', ioU(this.groups.vnew[0], this.groups.vnew[1])),
      vp: mk('probes', ioU(this.groups.vp[0], this.groups.vp[1])),
    };
  }

  private extGroups: Record<string, GPUBindGroup> = {};
  private nd = 0;
  private ns = 0;

  private buildExtras(): void {
    const p = this.plan, dev = this.device, S = BU.STORAGE;
    const bits = (v: number) => new Uint32Array(new Float32Array([v]).buffer)[0];
    if (p.disp) {
      const d = p.disp, nd = d.idx.length;
      this.nd = nd;
      const dl = new Uint32Array(2 * nd);
      for (let e = 0; e < nd; e++) { dl[2 * e] = d.idx[e]; dl[2 * e + 1] = bits(d.cb[e]); }
      this.buffer('DL', dl, S);
      this.buffer('DB', f32(d.beta), S);
      const dj = this.sim.stateArray(2, 0, 0, d.beta.length), dold = this.sim.stateArray(3, 0, 0, nd);
      this.buffer('DJ', dj ? f32(dj) : new Float32Array(d.beta.length), S | BU.COPY_SRC);
      this.buffer('DO', dold ? f32(dold) : new Float32Array(nd), S | BU.COPY_SRC);
      const u = new Float32Array(4 + 16);
      const ui = new Uint32Array(u.buffer);
      ui[0] = nd; ui[1] = d.npole; u[2] = 1 / p.dt;
      for (let q = 0; q < d.npole; q++) u[4 + q] = d.k[q];
      this.buffer('DU', u, BU.UNIFORM);
      const mod = dev.createShaderModule({ code: dispWGSL });
      const ents = ['DU', 'E', 'DL', 'DB', 'DJ', 'DO'].map((nm, i) => ({ binding: i, resource: { buffer: this.bufs[nm] } }));
      const use: Record<string, number[]> = { dispSave: [0, 1, 2, 5], dispCorrect: [0, 1, 2, 4], dispUpdate: [0, 1, 2, 3, 4, 5] };
      for (const e of Object.keys(use)) {
        this.pipes[e] = dev.createComputePipeline({ layout: 'auto', compute: { module: mod, entryPoint: e } });
        this.extGroups[e] = dev.createBindGroup({ layout: this.pipes[e].getBindGroupLayout(0),
          entries: ents.filter((x) => use[e].includes(x.binding)) });
      }
    }
    if (p.sibc) {
      const sb = p.sibc, ns = sb.idx.length, K = sb.c.length;
      this.ns = ns;
      const sl = new Uint32Array(12 * ns);
      for (let q = 0; q < ns; q++) {
        sl[12 * q] = sb.idx[q];
        for (let f = 0; f < 4; f++) { sl[12 * q + 1 + f] = sb.hidx[4 * q + f]; sl[12 * q + 5 + f] = bits(sb.hw[4 * q + f]); }
      }
      this.buffer('SL', sl, S);
      const sy = this.sim.stateArray(4, 0, 0, ns * K);
      this.buffer('SY', sy ? f32(sy) : new Float32Array(ns * K), S | BU.COPY_SRC);
      const u = new Float32Array(4 + 16 + 16);
      const ui = new Uint32Array(u.buffer);
      ui[0] = ns; ui[1] = K;
      for (let k = 0; k < K; k++) { u[4 + k] = sb.c[k]; u[20 + k] = sb.e[k]; }
      this.buffer('SU', u, BU.UNIFORM);
      const mod = dev.createShaderModule({ code: sibcWGSL });
      this.pipes.sibc = dev.createComputePipeline({ layout: 'auto', compute: { module: mod, entryPoint: 'sibc' } });
      this.extGroups.sibc = dev.createBindGroup({ layout: this.pipes.sibc.getBindGroupLayout(0),
        entries: ['SU', 'E', 'H', 'SL', 'SY'].map((nm, i) => ({ binding: i, resource: { buffer: this.bufs[nm] } })) });
    }
  }

  private mainGroups: Record<string, GPUBindGroup> = {};
  private ioGroups: Record<string, GPUBindGroup> = {};

  // bindings each entry point actually references ('auto' layouts drop unused ones)
  private usesBinding(entry: string, b: number): boolean {
    const use: Record<string, number[]> = {
      updateH: [0, 1, 2, 5], updateE: [0, 1, 2, 3, 4, 5],
      pmlHx: [0, 1, 2, 5, 6], pmlHy: [0, 1, 2, 5, 6], pmlHz: [0, 1, 2, 5, 6],
      pmlEx: [0, 1, 2, 4, 5, 6], pmlEy: [0, 1, 2, 4, 5, 6], pmlEz: [0, 1, 2, 4, 5, 6],
    };
    return use[entry].includes(b);
  }
  private usesIO(entry: string, b: number): boolean {
    const use: Record<string, number[]> = { inject: [0, 1, 3, 4, 7], probes: [0, 1, 2, 5, 6, 7], tick: [7] };
    return use[entry].includes(b);
  }

  private dispatch(pass: GPUComputePassEncoder, pipe: string, group: GPUBindGroup, threads: number, wg = WG): void {
    if (threads <= 0) return;
    const groups = Math.ceil(threads / wg);
    const x = Math.min(groups, 65535), y = Math.ceil(groups / 65535);
    pass.setPipeline(this.pipes[pipe]);
    pass.setBindGroup(0, group);
    pass.dispatchWorkgroups(x, y);
  }

  private encodeStep(pass: GPUComputePassEncoder): void {
    const p = this.plan, N = p.N, g = this.mainGroups;
    const per = [(p.ny + 1) * (p.nz + 1), (p.nx + 1) * (p.nz + 1), (p.nx + 1) * (p.ny + 1)];
    this.dispatch(pass, 'updateH', g.updateH, N);
    if (p.hasPml) {
      this.dispatch(pass, 'pmlHx', g.pmlHx, this.slabCount.H[0] * per[0]);
      this.dispatch(pass, 'pmlHy', g.pmlHy, this.slabCount.H[1] * per[1]);
      this.dispatch(pass, 'pmlHz', g.pmlHz, this.slabCount.H[2] * per[2]);
    }
    this.dispatch(pass, 'probes', this.ioGroups.ip, this.groups.ip[1] - this.groups.ip[0], 64);
    this.dispatch(pass, 'probes', this.ioGroups.vold, this.groups.vold[1] - this.groups.vold[0], 64);
    if (this.nd) this.dispatch(pass, 'dispSave', this.extGroups.dispSave, this.nd);
    this.dispatch(pass, 'updateE', g.updateE, N);
    if (p.hasPml) {
      this.dispatch(pass, 'pmlEx', g.pmlEx, this.slabCount.E[0] * per[0]);
      this.dispatch(pass, 'pmlEy', g.pmlEy, this.slabCount.E[1] * per[1]);
      this.dispatch(pass, 'pmlEz', g.pmlEz, this.slabCount.E[2] * per[2]);
    }
    if (this.nd) this.dispatch(pass, 'dispCorrect', this.extGroups.dispCorrect, this.nd);
    if (this.ns) this.dispatch(pass, 'sibc', this.extGroups.sibc, this.ns);
    this.dispatch(pass, 'inject', this.ioGroups.inject, this.nsrc, 64);
    this.dispatch(pass, 'probes', this.ioGroups.vnew, this.groups.vnew[1] - this.groups.vnew[0], 64);
    if (this.nd) this.dispatch(pass, 'dispUpdate', this.extGroups.dispUpdate, this.nd);
    this.dispatch(pass, 'probes', this.ioGroups.vp, this.groups.vp[1] - this.groups.vp[0], 64);
    this.dispatch(pass, 'tick', this.ioGroups.tick, 1, 1);
  }

  /** Advance n steps; probe samples are pushed into the C sim after every chunk. */
  async step(n: number, onChunk?: (stepsDone: number) => void): Promise<void> {
    const dev = this.device, p = this.plan;
    let left = n;
    while (left > 0) {
      const m = Math.min(this.chunk, left);
      const w = this.sim.waveforms(this.stepsDone, m);
      const wave = new Float32Array(2 * this.chunk);
      for (let q = 0; q < m; q++) { wave[2 * q] = w.lumped[q]; wave[2 * q + 1] = w.mode[q]; }
      dev.queue.writeBuffer(this.bufs.WAVE, 0, wave);
      dev.queue.writeBuffer(this.bufs.STATE, 0, new Uint32Array([0]));
      const enc = dev.createCommandEncoder();
      const pass = enc.beginComputePass();
      for (let q = 0; q < m; q++) this.encodeStep(pass);
      pass.end();
      const bytes = 4 * (1 + m * Math.max(1, this.nprobe));
      enc.copyBufferToBuffer(this.bufs.STATE, 0, this.bufs.READ, 0, bytes);
      dev.queue.submit([enc.finish()]);
      await this.bufs.READ.mapAsync(MAP_READ_MODE, 0, bytes);
      const raw = new Float32Array(this.bufs.READ.getMappedRange(0, bytes).slice(0));
      this.bufs.READ.unmap();
      const np = this.nprobe;
      const take = (g: number[]) => {
        const len = g[1] - g[0], out = new Float64Array(m * len);
        for (let q = 0; q < m; q++) for (let k = 0; k < len; k++) out[q * len + k] = raw[1 + q * np + g[0] + k];
        return out;
      };
      const ip = take(this.groups.ip), vold = take(this.groups.vold), vnew = take(this.groups.vnew), vp = take(this.groups.vp);
      let pk = 0;
      for (const a of [ip, vnew, vp]) for (let i = 0; i < a.length; i++) pk = Math.max(pk, Math.abs(a[i]));
      this.lastChunkPeak = pk;
      this.sim.pushSamples(m, ip, vold, vnew, vp);
      this.stepsDone += m;
      left -= m;
      onChunk?.(this.stepsDone);
    }
    void p;
  }

  private async readBuffer(name: string, offset: number, bytes: number): Promise<Float32Array> {
    const dev = this.device;
    const tmp = dev.createBuffer({ size: bytes, usage: BU.MAP_READ | BU.COPY_DST });
    const enc = dev.createCommandEncoder();
    enc.copyBufferToBuffer(this.bufs[name], offset, tmp, 0, bytes);
    dev.queue.submit([enc.finish()]);
    await tmp.mapAsync(MAP_READ_MODE);
    const out = new Float32Array(tmp.getMappedRange().slice(0));
    tmp.unmap(); tmp.destroy();
    return out;
  }

  /** Copy the GPU state (fields, CPML, Debye, SIBC) back into the C simulation,
   * e.g. before cem_checkpoint_save. The DFTs are already in C. */
  async syncToSim(): Promise<void> {
    const p = this.plan, N = p.N;
    for (let c = 0; c < 6; c++) this.sim.setField(c, await this.readField(c));
    if (p.hasPml) {
      for (const sl of this.psiSlots) {
        if (sl.len) this.sim.stateArray(sl.which, sl.c, sl.w, sl.len)?.set(await this.readBuffer('PSI', 4 * sl.off, 4 * sl.len));
      }
    }
    if (this.nd && p.disp) {
      this.sim.stateArray(2, 0, 0, p.disp.beta.length)?.set(await this.readBuffer('DJ', 0, 4 * p.disp.beta.length));
      this.sim.stateArray(3, 0, 0, this.nd)?.set(await this.readBuffer('DO', 0, 4 * this.nd));
    }
    if (this.ns && p.sibc) {
      const len = this.ns * p.sibc.c.length;
      this.sim.stateArray(4, 0, 0, len)?.set(await this.readBuffer('SY', 0, 4 * len));
    }
  }

  /** Read back one field component (0..2 E, 3..5 H) as f32. */
  async readField(comp: number): Promise<Float32Array> {
    const N = this.plan.N, dev = this.device;
    const tmp = dev.createBuffer({ size: 4 * N, usage: BU.MAP_READ | BU.COPY_DST });
    const enc = dev.createCommandEncoder();
    enc.copyBufferToBuffer(comp < 3 ? this.bufs.E : this.bufs.H, 4 * N * (comp % 3), tmp, 0, 4 * N);
    dev.queue.submit([enc.finish()]);
    await tmp.mapAsync(MAP_READ_MODE);
    const out = new Float32Array(tmp.getMappedRange().slice(0));
    tmp.unmap(); tmp.destroy();
    return out;
  }

  /** Total field energy proxy (sum of squares) for run control; reads both fields back. */
  async energyProxy(): Promise<number> {
    let s = 0;
    for (let c = 0; c < 6; c++) { const f = await this.readField(c); for (let i = 0; i < f.length; i++) s += f[i] * f[i]; }
    return s;
  }

  destroy(): void { for (const b of Object.values(this.bufs)) b.destroy(); }
}
