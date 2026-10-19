/// <reference lib="webworker" />
/**
 * Solver worker: builds and runs a simulation off the main thread.
 * In:  { type: 'run-gcpw', params, backend, wasmUrl }
 *      { type: 'run-board', input (KiCad text or Gerber files), stackup, net, fmax, df, hmin, lossy, backend, wasmUrl }
 * Out: 'setup' | 'progress' | 'done' (S-parameters, Z0, field snapshot) | 'error'
 */
import { Courant, type Complex, type Sim } from './courant.ts';
import { buildGcpw, gcpwResults, type GcpwParams } from '../scenarios/gcpw.ts';
import { GpuEngine } from '../gpu/engine.ts';
import { loadBoard } from '../import/load.ts';
import { buildBoardSim } from '../import/boardsim.ts';
import { dropCheckpoint, loadCheckpoint, runKey, saveCheckpoint } from './checkpoints.ts';

declare const self: DedicatedWorkerGlobalScope;

interface Snap { x: number[]; y: number[]; k: number; z: number }

/** |E| transverse magnitude on the node plane k (just under the signal layer). */
function slice(sim: Sim, E: (c: number) => ArrayLike<number>, snap: Snap) {
  const nx = sim.nx + 1, ny = sim.ny + 1, nz1 = sim.nz + 1;
  const ey = E(1), ez = E(2);
  const data = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    const q = (i * ny + j) * nz1 + snap.k;
    data[i * ny + j] = Math.hypot(ey[q], ez[q]);
  }
  return { nx, ny, x: snap.x, y: snap.y, z: snap.z, data };
}

/** Advance to energy decay on the chosen backend, grabbing one field snapshot at tSnap.
 * Saves a checkpoint every ~20 s of wall time and resumes from one if present. */
async function runSim(sim: Sim, tStop: number, decay: number, backend: 'cpu' | 'gpu', tSnap: number, snap: Snap,
                      key: string, tMax = 5e-9) {
  const dt = sim.dt, maxSteps = Math.floor(tMax / dt), chunk = 500;
  let lastRatio = 1;
  let peak = 0;
  let shot: ReturnType<typeof slice> | undefined;
  const saved = await loadCheckpoint(key);
  if (saved) {
    try {
      sim.restore(saved.subarray(0, saved.length - 1));
      peak = saved[saved.length - 1];
      self.postMessage({ type: 'resumed', steps: sim.steps });
    } catch (e) { console.warn('checkpoint ignored:', e instanceof Error ? e.message : e); /* stale layout: start over */ }
  }
  let lastSave = performance.now();
  const maybeSave = async (eng?: GpuEngine) => {
    if (performance.now() - lastSave < 20000) return;
    if (eng) await eng.syncToSim();
    const cp = sim.checkpoint();
    const buf = new Float64Array(cp.length + 1);
    buf.set(cp); buf[cp.length] = peak; // run-control peak travels with the state
    const ok = await saveCheckpoint(key, buf);
    self.postMessage({ type: 'checkpoint', steps: sim.steps, ok, mb: buf.byteLength / 1e6 });
    lastSave = performance.now();
  };
  if (backend === 'gpu') {
    const eng = await GpuEngine.create(sim, { chunk });
    while (sim.steps < maxSteps) {
      // the snapshot splits a chunk but never moves the run-control checkpoints
      const toSnap = Math.round(tSnap / dt) - sim.steps;
      if (!shot && toSnap > 0 && toSnap < chunk) {
        await eng.step(toSnap);
        const f1 = await eng.readField(1), f2 = await eng.readField(2);
        shot = slice(sim, (c) => (c === 1 ? f1 : f2), snap);
        await eng.step(chunk - toSnap);
      } else {
        await eng.step(chunk);
      }
      peak = Math.max(peak, eng.lastChunkPeak);
      const ratio = peak > 0 ? (eng.lastChunkPeak / peak) ** 2 : 1;
      lastRatio = ratio;
      self.postMessage({ type: 'progress', steps: sim.steps, time: sim.steps * dt, energyRatio: ratio });
      if (sim.steps * dt > tStop && ratio < decay) break;
      await maybeSave(eng);
    }
    eng.destroy();
  } else {
    while (sim.steps < maxSteps) {
      const toSnap = Math.round(tSnap / dt) - sim.steps;
      if (!shot && toSnap > 0 && toSnap < chunk) {
        sim.step(toSnap);
        shot = slice(sim, (c) => sim.field(c), snap);
        sim.step(chunk - toSnap);
      } else {
        sim.step(chunk);
      }
      const W = sim.energy();
      peak = Math.max(peak, W);
      lastRatio = peak > 0 ? W / peak : 1;
      self.postMessage({ type: 'progress', steps: sim.steps, time: sim.steps * dt, energyRatio: peak > 0 ? W / peak : 1 });
      if (sim.steps * dt > tStop && W < decay * peak) break;
      await maybeSave();
    }
  }
  await dropCheckpoint(key);
  // hit the time limit before the energy decayed: say so, with the level reached
  const capped = sim.steps >= maxSteps ? { energyDb: 10 * Math.log10(Math.max(lastRatio, 1e-30)), tMax } : undefined;
  return { shot, capped };
}

/** Threads: 0 = auto (all cores but one, up to 8), 1 = single-threaded build. The
 * multithreaded build needs cross-origin isolation; without it we fall back. */
async function loadCourant(wasmUrl: string, threads: number): Promise<Courant> {
  const iso = (self as unknown as { crossOriginIsolated?: boolean }).crossOriginIsolated === true;
  const n = threads > 0 ? threads : Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 2) - 1));
  if (!iso || n <= 1) {
    self.postMessage({ type: 'threads', threads: 1, isolated: iso });
    return Courant.load(fetch(wasmUrl));
  }
  const mtUrl = new URL('./courant-mt.wasm', wasmUrl).href;
  const c = await Courant.loadThreaded(fetch(mtUrl), n, new URL('./mtworker.js', import.meta.url));
  self.postMessage({ type: 'threads', threads: n, isolated: iso });
  return c;
}

const div = (x: Complex, y: Complex): Complex => {
  const d = y.re * y.re + y.im * y.im;
  return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d };
};

let cemRef: Courant | undefined;

self.onmessage = async (ev: MessageEvent) => {
  const m = ev.data;
  const key = await runKey(m);
  try {
    const cem = await loadCourant(m.wasmUrl, m.threads ?? 0);
    cemRef = cem;
    const t0 = performance.now();
    if (m.type === 'run-board') {
      const board = loadBoard(m.input, m.stackup).board;
      // one run per excited port: port 1 always; port 2 too for a full (non-symmetric) 2-port
      const nRuns = m.full2 ? 2 : 1;
      const waves: { w1: ReturnType<Sim['modePortWaves']>; w2: ReturnType<Sim['modePortWaves']> }[] = [];
      let freqs: number[] = [], steps = 0, shot: ReturnType<typeof slice> | undefined;
      let capped: { energyDb: number; tMax: number } | undefined;
      for (let run = 0; run < nRuns; run++) {
        const bs = buildBoardSim(cem, board, { fmax: m.fmax, df: m.df, hmin: m.hmin, lossy: m.lossy, thickCopper: m.thick ?? true, yWindow: m.ywin || undefined,
          excite: run, ports: [{ net: m.net, edge: 'xlo', offset: m.offsets?.[0] ?? 0 }, { net: m.net, edge: 'xhi', offset: m.offsets?.[1] ?? 0 }] });
        const [nx, ny, nz] = bs.cells;
        self.postMessage({ type: 'setup', cells: bs.cells, dt: bs.sim.dt, mcells: (nx * ny * nz) / 1e6, run, runs: nRuns });
        const top = board.layers[board.layers.length - 1];
        const kTop = bs.mesh.z.findIndex((v) => Math.abs(v - top.z0) < 1e-12);
        const L = board.outline.max[0] - board.outline.min[0] - (m.offsets?.[0] ?? 0) - (m.offsets?.[1] ?? 0);
        const tSnap = bs.tStop * 0.45 + (0.5 * L + 4e-3) * 1.7 / 299792458.0;
        const rs = await runSim(bs.sim, bs.tStop, m.decay ?? 1e-5, m.backend, tSnap,
          { x: bs.mesh.x, y: bs.mesh.y, k: Math.max(0, kTop - 1), z: top.z0 }, run ? `${key}:p${run + 1}` : key, m.tmax ?? 5e-9);
        if (run === 0) shot = rs.shot;
        capped = capped ?? rs.capped;
        waves.push({ w1: bs.sim.modePortWaves(bs.ports[0]), w2: bs.sim.modePortWaves(bs.ports[1]) });
        freqs = bs.freqs; steps += bs.sim.steps;
        bs.sim.destroy();
      }
      const { w1, w2 } = waves[0];
      const s11 = w1.b.map((b, i) => div(b, w1.a[i])), s21 = w2.b.map((b, i) => div(b, w1.a[i]));
      // with port 2 excited: S12 = b1/a2, S22 = b2/a2; otherwise assume a symmetric, reciprocal line
      const s12 = waves[1] ? waves[1].w1.b.map((b, i) => div(b, waves[1].w2.a[i])) : undefined;
      const s22 = waves[1] ? waves[1].w2.b.map((b, i) => div(b, waves[1].w2.a[i])) : undefined;
      const s11_50: Complex[] = [], s21_50: Complex[] = [], s12_50: Complex[] = [], s22_50: Complex[] = [];
      s11.forEach((a, i) => {
        const r = cem.renormalize2([a, s12?.[i] ?? s21[i], s21[i], s22?.[i] ?? a], w1.z0[i], w2.z0[i], 50);
        s11_50.push(r[0]); s12_50.push(r[1]); s21_50.push(r[2]); s22_50.push(r[3]);
      });
      self.postMessage({ type: 'done', f: freqs, s11, s21, s11_50, s21_50, z0: Array.from(w1.z0), steps,
        ...(s12 ? { s12, s22, s12_50, s22_50 } : {}), wall: (performance.now() - t0) / 1000, snapshot: shot, capped });
      return;
    }
    if (m.type === 'run-gcpw') {
      const p = m.params as GcpwParams;
      const setup = buildGcpw(cem, p);
      const [nx, ny, nz] = setup.cells;
      self.postMessage({ type: 'setup', cells: setup.cells, dt: setup.sim.dt, mcells: (nx * ny * nz) / 1e6 });
      const tSnap = setup.tStop * 0.45 + 0.5 * p.len * 1.65 / 299792458.0 + (p.ports === 'mode' ? 4e-3 * 1.65 / 299792458.0 : 0);
      const { shot, capped } = await runSim(setup.sim, setup.tStop, p.decay, m.backend, tSnap,
        { x: setup.mesh.x, y: setup.mesh.y, k: setup.mesh.kh - 1, z: 254e-6 }, key, (p as GcpwParams & { tmax?: number }).tmax ?? 5e-9);
      const res = gcpwResults(setup, cem);
      self.postMessage({ type: 'done', ...res, steps: setup.sim.steps, wall: (performance.now() - t0) / 1000, snapshot: shot, capped });
      setup.sim.destroy();
    }
  } catch (e) {
    self.postMessage({ type: 'error', message: e instanceof Error ? e.message : String(e) });
  } finally {
    cemRef?.pool?.terminate();
  }
};
