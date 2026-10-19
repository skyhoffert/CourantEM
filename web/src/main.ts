/**
 * Courant EM app shell: one typed store, one-way data flow, panels with
 * mount/update. Solving happens in a Web Worker (CPU f64 wasm or GPU f32 WebGPU).
 */
import { Store, type Panel, type State } from './app/store.ts';
import { PlotPanel } from './app/plots.ts';
import { ModelPanel } from './app/view3d.ts';
import { defaultGcpw, gcpwMesh } from './scenarios/gcpw.ts';
import { estimate } from './mesh/mesher.ts';
import { planBoard } from './import/boardsim.ts';
import { GCPW_NET, gcpwBoard } from './scenarios/gcpwboard.ts';
import { inputFromFiles, loadBoard, type BoardInput } from './import/load.ts';
import type { StackupRow } from './import/gerber.ts';
import { toTouchstone } from './app/touchstone.ts';

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

const store = new Store({
  source: 'gcpw',
  params: { hmin: 20e-6, dx: 100e-6, len: 10e-3, er: 3.66, fmax: 40e9, df: 0.5e9, decay: 1e-5, ports: 'mode', lossy: false, thick: true, full2: false, ywin: 0, tmax: 5e-9 },
  backend: 'cpu',
  threads: 0,
  gpuAvailable: 'gpu' in navigator,
  ref50: true,
  run: { phase: 'idle', message: 'Substrate 254 µm, strip 500 µm, gaps 200 µm, 1 oz copper, via fences (0.3 mm drills, 0.8 mm pitch) 300 µm from the gaps.', progress: 0 },
});

// ---- controls: DOM -> actions; state -> DOM ----
class ControlsPanel implements Panel {
  constructor() {
    const num = (id: string, scale: number, key: keyof State['params']) =>
      $<HTMLInputElement>(id).addEventListener('change', (e) =>
        store.dispatch({ type: 'params', patch: { [key]: Number((e.target as HTMLInputElement).value) * scale } }));
    num('tmax', 1e-9, 'tmax');
    $<HTMLInputElement>('decaydb').addEventListener('change', (e) =>
      store.dispatch({ type: 'params', patch: { decay: 10 ** (Math.min(-10, Number((e.target as HTMLInputElement).value)) / 10) } }));
    num('hmin', 1e-6, 'hmin'); num('dx', 1e-6, 'dx'); num('len', 1e-3, 'len'); num('er', 1, 'er'); num('ywin', 1e-3, 'ywin');
    $<HTMLSelectElement>('ports').addEventListener('change', (e) =>
      store.dispatch({ type: 'params', patch: { ports: (e.target as HTMLSelectElement).value as 'mode' | 'lumped' } }));
    $<HTMLInputElement>('lossy').addEventListener('change', (e) =>
      store.dispatch({ type: 'params', patch: { lossy: (e.target as HTMLInputElement).checked } }));
    $<HTMLInputElement>('full2').addEventListener('change', (e) =>
      store.dispatch({ type: 'params', patch: { full2: (e.target as HTMLInputElement).checked } }));
    $<HTMLInputElement>('thick').addEventListener('change', (e) =>
      store.dispatch({ type: 'params', patch: { thick: (e.target as HTMLInputElement).checked } }));
    $<HTMLInputElement>('ref50').addEventListener('change', (e) =>
      store.dispatch({ type: 'ref50', on: (e.target as HTMLInputElement).checked }));
    $<HTMLSelectElement>('threads').addEventListener('change', (e) =>
      store.dispatch({ type: 'threads', threads: Number((e.target as HTMLSelectElement).value) }));
    $<HTMLSelectElement>('backend').addEventListener('change', (e) =>
      store.dispatch({ type: 'backend', backend: (e.target as HTMLSelectElement).value as 'cpu' | 'gpu' }));
    $<HTMLSelectElement>('source').addEventListener('change', (e) =>
      store.dispatch({ type: 'source', source: (e.target as HTMLSelectElement).value as State['source'] }));
    $<HTMLSelectElement>('net').addEventListener('change', (e) =>
      store.dispatch({ type: 'net', net: (e.target as HTMLSelectElement).value }));
    $<HTMLInputElement>('file').addEventListener('change', async (e) => {
      const files = [...((e.target as HTMLInputElement).files ?? [])];
      if (!files.length) return;
      try {
        loadInput(await inputFromFiles(files));
      } catch (err) {
        store.dispatch({ type: 'status', message: `Import failed: ${err instanceof Error ? err.message : err}` });
      }
    });
    for (const [id, side] of [['offlo', 0], ['offhi', 1]] as const) {
      $<HTMLInputElement>(id).addEventListener('change', (e) => {
        const b = store.get().board;
        if (!b) return;
        const off: [number, number] = [...b.offsets];
        off[side] = Math.max(0, Number((e.target as HTMLInputElement).value) * 1e-3);
        try { loadInput(b.input, b.stackup, off); } catch (err) {
          store.dispatch({ type: 'status', message: `Port plane: ${err instanceof Error ? err.message : err}` });
        }
      });
    }
    $<HTMLButtonElement>('run').addEventListener('click', () => run());
    $<HTMLButtonElement>('save').addEventListener('click', () => {
      const s = store.get();
      if (!s.result) return;
      const name = s.source === 'board' && s.board ? s.board.name.replace(/\.kicad_pcb$/, '').replace(/ \(Gerber\)$/, '') : s.source === 'ref' ? 'gcpw-ref' : 'gcpw';
      const blob = new Blob([toTouchstone(s.result, name)], { type: 'text/plain' });
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `${name}.s2p`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    });
  }
  update(s: State): void {
    $<HTMLSelectElement>('source').value = s.source;
    const board = s.source === 'board';
    const ref = s.source === 'ref';
    // the built-in and board lines are meshed by the rule-based mesher (no dx) with mode ports
    for (const id of ['dx', 'ports']) ($(id) as HTMLInputElement).disabled = !ref;
    for (const id of ['len', 'er']) ($(id) as HTMLInputElement).disabled = board;
    ($('thick') as HTMLInputElement).disabled = ref || (board && !s.board);
    for (const id of ['offlo', 'offhi', 'full2', 'ywin']) ($(id) as HTMLInputElement).disabled = !board || !s.board;
    if (s.board && document.activeElement !== $('offlo')) $<HTMLInputElement>('offlo').value = String(+(s.board.offsets[0] * 1e3).toFixed(4));
    if (s.board && document.activeElement !== $('offhi')) $<HTMLInputElement>('offhi').value = String(+(s.board.offsets[1] * 1e3).toFixed(4));
    const sel = $<HTMLSelectElement>('net');
    const names = s.board?.nets ?? [];
    if (sel.options.length !== names.length || [...sel.options].some((o, i) => o.value !== names[i])) {
      sel.innerHTML = '';
      for (const n of names) sel.add(new Option(n, n));
    }
    if (s.board) sel.value = s.board.net;
    const gpuOpt = $<HTMLSelectElement>('backend').querySelector('option[value="gpu"]') as HTMLOptionElement;
    gpuOpt.disabled = !s.gpuAvailable;
    $<HTMLButtonElement>('run').disabled = s.run.phase === 'running' || (board && !s.board?.net);
    $<HTMLButtonElement>('save').disabled = !s.result;
  }
}

/** (Re)import a board, optionally with an edited stackup, and publish it to the store. */
function loadInput(input: BoardInput, stackup?: StackupRow[], offsets?: [number, number]): void {
  const prev = store.get().board;
  const prevNet = prev?.net;
  const off = offsets ?? (prev && prev.input === input ? prev.offsets : [0, 0]);
  const lb = loadBoard(input, stackup, off);
  const b = lb.board;
  const [x0] = b.outline.min, [x1] = b.outline.max;
  const net = lb.portNets.includes(prevNet ?? '') ? prevNet! : lb.portNets[0] ?? '';
  store.dispatch({ type: 'board', board: { name: input.name, input, stackup: lb.stackup, offsets: off, model: b, nets: lb.portNets, net } });
  const kind = input.kind === 'kicad' ? '' : ` from ${input.files.length} files`;
  store.dispatch({ type: 'status', message: `Loaded ${input.name}${kind}: ${b.copperLayers.length} copper layers, ` +
    `${b.tracks.length} tracks, ${b.vias.length} vias, ${((b.thickness) * 1e3).toFixed(3)} mm thick, board ${((x1 - x0) * 1e3).toFixed(1)} mm long; ` +
    (lb.portNets.length ? 'pick the net to excite and press Run.' : 'no x-directed track crosses both port planes; if the line starts at connector pads, move the port planes past them.') +
    (lb.warnings.length ? ` Warnings: ${lb.warnings.slice(0, 3).join('; ')}${lb.warnings.length > 3 ? ' …' : ''}` : '') });
}

// ---- stackup editor ----
class StackupPanel implements Panel {
  private shape = '';
  constructor() {
    $<HTMLInputElement>('mask').addEventListener('change', (e) => {
      const b = store.get().board;
      if (!b) return;
      const on = (e.target as HTMLInputElement).checked;
      const rows = b.stackup.filter((r) => !r.mask);
      if (on) rows.unshift({ name: 'F.Mask', type: 'dielectric', mask: true, thickness: 20e-6, epsR: 3.5, lossTangent: 0.02 });
      setTimeout(() => {
        try { loadInput(b.input, rows, b.offsets); } catch (err) {
          store.dispatch({ type: 'status', message: `Stackup: ${err instanceof Error ? err.message : err}` });
        }
      });
    });
  }
  private inputs: { inp: HTMLInputElement; get: (r: StackupRow) => number | undefined; scale: number }[] = [];
  update(s: State, prev: State | undefined): void {
    $<HTMLDivElement>('stackcard').hidden = !(s.source === 'board' && s.board);
    if (!s.board || (prev && prev.board?.rev === s.board.rev)) return;
    const rows = s.board.stackup;
    $<HTMLInputElement>('mask').checked = rows.some((r) => r.mask);
    const shape = rows.map((r) => `${r.name}:${r.type}`).join('|');
    if (shape !== this.shape) this.build(rows, shape);
    // refresh values in place (never rebuild under a focused input)
    this.inputs.forEach((f, k) => {
      const r = rows[Math.floor(k / 3)];
      const v = f.get(r);
      if (v !== undefined && document.activeElement !== f.inp) f.inp.value = String(+(v / f.scale).toPrecision(6));
    });
  }
  private build(rows: StackupRow[], shape: string): void {
    this.shape = shape;
    this.inputs = [];
    const tbl = $<HTMLTableElement>('stackup');
    tbl.replaceChildren();
    const head = tbl.insertRow();
    for (const h of ['Layer', 't (µm)', 'εr', 'tan δ', 'Rq (µm)']) { const th = document.createElement('th'); th.textContent = h; head.append(th); }
    rows.forEach((r, i) => {
      const tr = tbl.insertRow();
      tr.className = r.type === 'copper' ? 'cu' : '';
      const name = tr.insertCell(); name.className = 'name'; name.textContent = r.name;
      // three editable slots per row (thickness, then εr/tanδ or roughness), kept in order for refreshes
      const slot = (key: keyof StackupRow | null, scale: number, step: string) => {
        const inp = document.createElement('input');
        inp.type = 'number'; inp.step = step;
        const get = (q: StackupRow) => (key ? ((q[key] as number | undefined) ?? 0) : undefined);
        this.inputs.push({ inp, get, scale });
        if (!key) return; // placeholder keeps three slots per row
        const td = tr.insertCell();
        inp.addEventListener('change', () => {
          const b = store.get().board!;
          const next = b.stackup.map((q, k) => (k === i ? { ...q, [key]: Number(inp.value) * scale } : q));
          setTimeout(() => {
            try { loadInput(b.input, next, b.offsets); } catch (err) {
              store.dispatch({ type: 'status', message: `Stackup: ${err instanceof Error ? err.message : err}` });
            }
          });
        });
        td.append(inp);
      };
      slot('thickness', 1e-6, 'any');
      if (r.type === 'dielectric') { slot('epsR', 1, '0.01'); slot('lossTangent', 1, '0.0001'); tr.insertCell(); }
      else { tr.insertCell(); tr.insertCell(); slot('roughness', 1e-6, '0.1'); slot(null, 1, ''); }
    });
  }
}

// ---- status and progress ----
class StatusPanel implements Panel {
  update(s: State): void {
    $<HTMLParagraphElement>('status').textContent = s.run.message;
    $<HTMLProgressElement>('bar').value = s.run.progress;
  }
}

// ---- pre-run estimate ----
class EstimatePanel implements Panel {
  update(s: State, prev: State | undefined): void {
    if (prev && prev.params === s.params && prev.source === s.source && prev.backend === s.backend && prev.board === s.board &&
      prev.threads === s.threads) return;
    const el = $<HTMLParagraphElement>('estimate');
    let m: { x: number[]; y: number[]; z: number[] }, tSim = 0.35e-9;
    if (s.source !== 'ref') {
      if (s.source === 'board' && !s.board?.net) { el.textContent = ''; return; }
      try {
        const sb = s.source === 'board' ? s.board : undefined;
        const model = sb ? sb.model : gcpwBoard(s.params), net = sb ? sb.net : GCPW_NET;
        const off = sb ? sb.offsets : [0, 0], ywin = sb ? s.params.ywin : 0;
        const p = planBoard(model, { fmax: s.params.fmax, df: s.params.df, hmin: s.params.hmin, thickCopper: s.params.thick,
          yWindow: ywin || undefined,
          ports: [{ net, edge: 'xlo', offset: off[0] }, { net, edge: 'xhi', offset: off[1] }] });
        m = p;
        // the source pulse (10.5 tau, as the run sets it) plus the transit through board and feeds with margin
        const tau = Math.sqrt(6.9) / (Math.PI * 1.25 * s.params.fmax);
        tSim = 10.5 * tau + 1.2 * (p.cx1 - p.cx0 + 8e-3) * Math.sqrt(p.erMax) / 299792458.0;
      } catch (err) { el.textContent = `Estimate: ${err instanceof Error ? err.message : err}`; return; }
    } else {
      m = gcpwMesh({ ...defaultGcpw, ...s.params });
    }
    const e = estimate(m.x, m.y, m.z, { tSim, backend: s.backend, poles: s.params.lossy ? 8 : 0 });
    // cell-steps per second: measured wasm (~50 M per thread, ~75% parallel efficiency), assumed discrete GPU
    const thr = !crossOriginIsolated ? 1 : s.threads || Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 2) - 1));
    const rate = s.backend === 'cpu' ? 50e6 * (thr > 1 ? 0.75 * thr : 1) : 1.5e9;
    el.textContent = `Estimate: ${(e.cells / 1e6).toFixed(2)} M cells, Δt ${(e.dt * 1e15).toFixed(1)} fs, ` +
      `~${e.steps} steps, ${(e.bytes / 1e6).toFixed(0)} MB, ~${Math.max(1, e.cellSteps / rate).toFixed(0)} s on ${s.backend.toUpperCase()}`;
  }
}

let worker: Worker | null = null;

function run(): void {
  const s = store.get();
  worker?.terminate();
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  store.dispatch({ type: 'run-start' });
  worker.onmessage = (ev) => {
    const m = ev.data;
    if (m.type === 'setup') {
      store.dispatch({ type: 'run-setup', cells: m.cells,
        message: `${m.runs > 1 ? `Run ${m.run + 1}/${m.runs}: ` : ''}Mesh ${m.cells.join(' × ')} = ${m.mcells.toFixed(2)} M cells, Δt ${(m.dt * 1e15).toFixed(1)} fs` });
    } else if (m.type === 'resumed') {
      store.dispatch({ type: 'status', message: `Resumed from a checkpoint at step ${m.steps}` });
      (window as unknown as { __resumed: number }).__resumed = m.steps;
    } else if (m.type === 'threads') {
      store.dispatch({ type: 'threads-used', threads: m.threads, isolated: m.isolated });
      (window as unknown as { __threads: number }).__threads = m.threads;
    } else if (m.type === 'checkpoint') {
      console.log(`checkpoint at step ${m.steps}: ${m.ok ? `saved (${m.mb.toFixed(1)} MB)` : 'not saved'}`);
    } else if (m.type === 'progress') {
      const r = Math.max(m.energyRatio, 1e-12);
      store.dispatch({ type: 'run-progress', progress: Math.min(1, Math.log10(1 / r) / 5),
        message: `Step ${m.steps}, t = ${(m.time * 1e12).toFixed(0)} ps, field energy ${(10 * Math.log10(r)).toFixed(1)} dB below peak` });
    } else if (m.type === 'done') {
      const z0 = m.z0 ? ` · line Z0 ${m.z0[0].toFixed(2)}–${m.z0[m.z0.length - 1].toFixed(2)} Ω` : '';
      const st = store.get();
      const how = st.backend === 'gpu' ? 'GPU' : st.threadsUsed && st.threadsUsed > 1 ? `${st.threadsUsed} threads` :
        st.isolated === false ? '1 thread, page not cross-origin isolated' : '1 thread';
      const cap = m.capped ? ` — stopped at the ${(m.capped.tMax * 1e9).toFixed(1)} ns limit with the field energy still ${m.capped.energyDb.toFixed(0)} dB below peak: ` +
        'slowly decaying resonances (often parallel-plate modes between grounds without stitching vias) truncate the spectra; raise the limit, add vias or losses, or widen the band' : '';
      store.dispatch({ type: 'run-done', result: m, message: `Done: ${m.steps} steps in ${m.wall.toFixed(1)} s (${how})${z0}${cap}` });
      (window as unknown as { __result: unknown }).__result = m;
    } else if (m.type === 'error') {
      store.dispatch({ type: 'run-error', message: `Error: ${m.message}` });
    }
  };
  const wasmUrl = new URL('./courant.wasm', document.baseURI).href;
  if (s.source === 'board' && s.board) {
    worker.postMessage({ type: 'run-board', input: s.board.input, stackup: s.board.stackup, offsets: s.board.offsets, net: s.board.net, fmax: s.params.fmax, df: s.params.df,
      hmin: s.params.hmin, lossy: s.params.lossy, thick: s.params.thick, full2: s.params.full2, ywin: s.params.ywin, decay: s.params.decay, tmax: s.params.tmax, backend: s.backend, threads: s.threads, wasmUrl });
  } else if (s.source === 'gcpw') {
    // symmetric line: one run, S22 = S11 and S12 = S21
    worker.postMessage({ type: 'run-board', builtin: { len: s.params.len, er: s.params.er }, offsets: [0, 0], net: GCPW_NET, fmax: s.params.fmax, df: s.params.df,
      hmin: s.params.hmin, lossy: s.params.lossy, thick: s.params.thick, full2: false, ywin: 0, decay: s.params.decay, tmax: s.params.tmax, backend: s.backend, threads: s.threads, wasmUrl });
  } else {
    worker.postMessage({ type: 'run-gcpw', params: { ...defaultGcpw, ...s.params }, backend: s.backend, threads: s.threads, wasmUrl });
  }
}

$('ver').textContent = __COURANT_VERSION__;
store.mount(new ControlsPanel());
store.mount(new StackupPanel());
store.mount(new StatusPanel());
store.mount(new EstimatePanel());
store.mount(new ModelPanel($('model')));
store.mount(new PlotPanel($('plots')));
