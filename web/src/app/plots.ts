/** S-parameter and impedance plots (uPlot). */
import uPlot from 'uplot';
import css from 'uplot/dist/uPlot.min.css';
import type { Complex } from '../solver/courant.ts';
import type { Panel, State } from './store.ts';

const db = (c: Complex) => 10 * Math.log10(c.re * c.re + c.im * c.im + 1e-300);
const deg = (c: Complex) => (Math.atan2(c.im, c.re) * 180) / Math.PI;

function unwrap(p: number[]): number[] {
  const out = [...p];
  for (let i = 1; i < out.length; i++) {
    let d = out[i] - out[i - 1];
    while (d > 180) { out[i] -= 360; d -= 360; }
    while (d < -180) { out[i] += 360; d += 360; }
  }
  return out;
}

let cssDone = false;
function injectCss() {
  if (cssDone) return;
  const st = document.createElement('style');
  st.textContent = css as unknown as string;
  document.head.appendChild(st);
  cssDone = true;
}

export class PlotPanel implements Panel {
  private el: HTMLElement;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private mag?: any;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private second?: any;
  constructor(el: HTMLElement) { this.el = el; injectCss(); }

  update(s: State, prev: State | undefined): void {
    if (prev && s.result === prev.result && s.ref50 === prev.ref50) return;
    this.el.innerHTML = '';
    this.mag?.destroy(); this.second?.destroy();
    const r = s.result;
    if (!r) { this.el.innerHTML = '<p class="muted">Run a simulation to see S-parameters.</p>'; return; }
    const use50 = s.ref50 && r.s11_50;
    const s11 = use50 ? r.s11_50! : r.s11, s21 = use50 ? r.s21_50! : r.s21;
    const s22 = use50 ? r.s22_50 : r.s22;
    const f = r.f.map((v) => v / 1e9);
    const css = getComputedStyle(document.documentElement);
    const fg = css.getPropertyValue('--fg').trim() || '#222', grid = css.getPropertyValue('--grid').trim() || '#ddd';
    const axis = { stroke: fg, grid: { stroke: grid, width: 1 }, ticks: { stroke: grid } };
    const width = Math.max(300, this.el.clientWidth - 8);
    const title = `S-parameters (${use50 ? '50 Ω reference' : 'normalised to the line Z0'})`;
    this.mag = new uPlot({
      title, width, height: 300,
      scales: { x: { time: false }, y: { range: [-80, 5] } },
      axes: [{ ...axis, label: 'GHz' }, { ...axis, label: 'dB' }],
      series: [{ label: 'f (GHz)' }, { label: '|S21|', stroke: '#2a7de1', width: 2 }, { label: '|S11|', stroke: '#e8590c', width: 2 },
        ...(s22 ? [{ label: '|S22|', stroke: '#c2255c', width: 1.5, dash: [6, 4] }] : [])],
    }, [f, s21.map(db), s11.map(db), ...(s22 ? [s22.map(db)] : [])], this.el);
    if (r.z0) {
      this.second = new uPlot({
        title: 'Port Z0 and ∠S21', width, height: 240,
        scales: { x: { time: false }, z: {}, p: {} },
        axes: [{ ...axis, label: 'GHz' }, { ...axis, scale: 'z', label: 'Z0 (Ω)' },
          { ...axis, scale: 'p', side: 1, label: '∠S21 (deg)', grid: { show: false } }],
        series: [{ label: 'f (GHz)' }, { label: 'Z0', stroke: '#2f9e44', width: 2, scale: 'z' },
          { label: '∠S21', stroke: '#7048e8', width: 1.5, scale: 'p' }],
      }, [f, r.z0, unwrap(s21.map(deg))], this.el);
    }
  }
}
