/**
 * In-browser GPU vs CPU test (run by test/gpu.test.mjs in headless Chromium).
 * 1. Parity scenario (every kernel feature except dispersion/SIBC): fields after
 *    N steps, f32 GPU vs f64 C.
 * 2. Small GCPW with mode ports: S-parameters, GPU vs C, against the plan's
 *    tolerance (|S21| +-0.01 dB, phase +-0.1 deg).
 */
import { Axis, Boundary, Courant, type Sim } from '../src/solver/courant.ts';
import { GpuEngine } from '../src/gpu/engine.ts';
import { buildGcpw, defaultGcpw } from '../src/scenarios/gcpw.ts';

const log = (m: string) => { (document.getElementById('log') as HTMLElement).textContent += m + '\n'; console.log(m); };

function parityScene(cem: Courant): Sim {
  const NX = 24, NY = 18, NZ = 14;
  const x: number[] = [], y: number[] = [], z: number[] = [];
  for (let i = 0; i <= NX; i++) x.push(1e-4 * i + 2e-6 * i * i);
  for (let j = 0; j <= NY; j++) y.push(1e-4 * (j - NY / 2) * (1 + 0.01 * Math.abs(j - NY / 2.0)));
  for (let k = 0; k <= NZ; k++) z.push(5e-5 * k + 1e-6 * k * k);
  const s = cem.createSim(x, y, z);
  for (let f = 0; f < 6; f++) if (f !== 4) s.setBoundary(f as 0, Boundary.CPML, 4);
  s.materialBox(-1, 1, -1, 1, 0, z[5], 3.4, 0.02);
  s.pecBox(-1, 1, y[NY - 6], 1, 0, z[5]);
  s.pecSheet(Axis.Z, z[5], x[5], x[19], y[8], y[10]);
  const p1 = s.addLumpedPort(Axis.Z, 5, 9, 0, 5, -1, 50);
  s.addLumpedPort(Axis.Z, 19, 9, 0, 5, -1, 50);
  s.setGaussian(p1, 1.0, 60e-12, 15e-12);
  s.setFrequencies([5e9, 20e9, 40e9]);
  s.addVProbe(Axis.Z, 12, 9, 0, 5, -1);
  s.addIProbe(Axis.X, 12, 7, 11, 4, 5);
  s.finalize(0.95);
  return s;
}

async function main() {
  const cem = await Courant.load(fetch('./courant.wasm'));
  const results: Record<string, number | string> = {};
  // ---- 1: field parity ----
  {
    const a = parityScene(cem), b = parityScene(cem);
    const steps = 600;
    a.step(steps);
    const eng = await GpuEngine.create(b, { chunk: 200 });
    await eng.step(steps);
    let worst = 0, scale = 0;
    for (let c = 0; c < 6; c++) {
      const ref = a.field(c), got = await eng.readField(c);
      for (let i = 0; i < ref.length; i++) { scale = Math.max(scale, Math.abs(ref[i])); }
      for (let i = 0; i < ref.length; i++) worst = Math.max(worst, Math.abs(got[i] - ref[i]));
    }
    // compare per component family relative to the largest field value
    const v1 = a.portV(0), v2 = b.portV(0);
    let dv = 0, sv = 0;
    v1.forEach((v, i) => { dv = Math.max(dv, Math.hypot(v.re - v2[i].re, v.im - v2[i].im)); sv = Math.max(sv, Math.hypot(v.re, v.im)); });
    results.field_rel = worst / scale;
    results.portV_rel = dv / sv;
    log(`parity scene: max |E,H| diff / max = ${(worst / scale).toExponential(2)}, port-1 V spectrum rel diff ${(dv / sv).toExponential(2)}`);
    eng.destroy(); a.destroy(); b.destroy();
  }
  // ---- 2: GCPW with mode ports ----
  {
    const p = { ...defaultGcpw, hmin: 40e-6, dx: 200e-6, len: 4e-3, ports: 'mode' as const };
    const A = buildGcpw(cem, p), B = buildGcpw(cem, p);
    const steps = 2500;
    const t0 = performance.now();
    A.sim.step(steps);
    const tc = performance.now() - t0;
    const eng = await GpuEngine.create(B.sim, { chunk: 250 });
    const t1 = performance.now();
    await eng.step(steps);
    const tg = performance.now() - t1;
    const wa1 = A.sim.modePortWaves(A.ports[0]), wa2 = A.sim.modePortWaves(A.ports[1]);
    const wb1 = B.sim.modePortWaves(B.ports[0]), wb2 = B.sim.modePortWaves(B.ports[1]);
    const div = (x: { re: number; im: number }, y: { re: number; im: number }) => {
      const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d };
    };
    let dmag = 0, dph = 0, ds11 = 0;
    A.freqs.forEach((f, i) => {
      if (f < 1e9) return;
      const sa = div(wa2.b[i], wa1.a[i]), sb = div(wb2.b[i], wb1.a[i]);
      const ma = 20 * Math.log10(Math.hypot(sa.re, sa.im)), mb = 20 * Math.log10(Math.hypot(sb.re, sb.im));
      const r = div(sb, sa);
      dmag = Math.max(dmag, Math.abs(ma - mb));
      dph = Math.max(dph, Math.abs(Math.atan2(r.im, r.re) * 180 / Math.PI));
      const s11a = div(wa1.b[i], wa1.a[i]), s11b = div(wb1.b[i], wb1.a[i]);
      ds11 = Math.max(ds11, Math.hypot(s11a.re - s11b.re, s11a.im - s11b.im));
    });
    results.s21_db = dmag; results.s21_deg = dph; results.s11_abs = ds11;
    results.cpu_ms = tc; results.gpu_ms = tg;
    log(`GCPW mode ports (${A.cells.join('x')} cells, ${steps} steps): max |S21| diff ${dmag.toFixed(5)} dB, phase ${dph.toFixed(4)} deg, |dS11| ${ds11.toExponential(2)}; CPU ${(tc / 1000).toFixed(1)} s, GPU (SwiftShader) ${(tg / 1000).toFixed(1)} s`);
    eng.destroy(); A.sim.destroy(); B.sim.destroy();
  }
  // ---- 3: lossy GCPW (Debye substrate + copper SIBC) ----
  {
    const p = { ...defaultGcpw, hmin: 40e-6, dx: 200e-6, len: 4e-3, ports: 'mode' as const, lossy: true };
    const A = buildGcpw(cem, p), B = buildGcpw(cem, p);
    const steps = 2500;
    A.sim.step(steps);
    const eng = await GpuEngine.create(B.sim, { chunk: 250 });
    await eng.step(steps);
    const wa1 = A.sim.modePortWaves(A.ports[0]), wa2 = A.sim.modePortWaves(A.ports[1]);
    const wb1 = B.sim.modePortWaves(B.ports[0]), wb2 = B.sim.modePortWaves(B.ports[1]);
    const div = (x: { re: number; im: number }, y: { re: number; im: number }) => {
      const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d };
    };
    let dmag = 0, dph = 0;
    A.freqs.forEach((f, i) => {
      if (f < 1e9) return;
      const sa = div(wa2.b[i], wa1.a[i]), sb = div(wb2.b[i], wb1.a[i]), r = div(sb, sa);
      dmag = Math.max(dmag, Math.abs(20 * Math.log10(Math.hypot(r.re, r.im))));
      dph = Math.max(dph, Math.abs(Math.atan2(r.im, r.re) * 180 / Math.PI));
    });
    results.lossy_s21_db = dmag; results.lossy_s21_deg = dph;
    log(`lossy GCPW: max |S21| diff ${dmag.toFixed(5)} dB, phase ${dph.toFixed(4)} deg`);
    // ---- 4: GPU half-run -> syncToSim -> checkpoint -> fresh C sim resumes ----
    const C = buildGcpw(cem, p), D = buildGcpw(cem, p);
    const eng2 = await GpuEngine.create(C.sim, { chunk: 250 });
    await eng2.step(1250);
    await eng2.syncToSim();
    D.sim.restore(C.sim.checkpoint());
    D.sim.step(steps - D.sim.steps);
    const wd1 = D.sim.modePortWaves(D.ports[0]), wd2 = D.sim.modePortWaves(D.ports[1]);
    let rmag = 0, rph = 0;
    A.freqs.forEach((f, i) => {
      if (f < 1e9) return;
      const sa = div(wa2.b[i], wa1.a[i]), sd = div(wd2.b[i], wd1.a[i]), r = div(sd, sa);
      rmag = Math.max(rmag, Math.abs(20 * Math.log10(Math.hypot(r.re, r.im))));
      rph = Math.max(rph, Math.abs(Math.atan2(r.im, r.re) * 180 / Math.PI));
    });
    results.resume_s21_db = rmag; results.resume_s21_deg = rph;
    log(`GPU->checkpoint->CPU resume (lossy): max |S21| diff vs CPU ${rmag.toFixed(5)} dB, phase ${rph.toFixed(4)} deg`);
    eng2.destroy(); C.sim.destroy(); D.sim.destroy();
    eng.destroy(); A.sim.destroy(); B.sim.destroy();
  }
  (window as unknown as { __gpu: unknown }).__gpu = results;
}

main().catch((e) => { log('ERROR ' + (e instanceof Error ? e.stack : String(e))); (window as unknown as { __gpu: unknown }).__gpu = { error: String(e) }; });
