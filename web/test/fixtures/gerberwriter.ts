/** Test helper: write a Board model as KiCad-style Gerbers + Excellon (FS 4.6 mm, integer coordinates). */
import type { Board } from '../../src/import/kicad.ts';

/** Write a Board as KiCad-style Gerbers (FS 4.6 mm, integer coordinates). */
export function toGerbers(b: Board, x2nets: boolean): { name: string; text: string }[] {
  const f = (v: number) => String(Math.round(v * 1e9)); // mm with 6 decimals = nm
  const out: { name: string; text: string }[] = [];
  const cu = [...b.copperLayers].reverse(); // top first
  cu.forEach((layer, k) => {
    const lines = ['%FSLAX46Y46*%', '%MOMM*%', `%TF.FileFunction,Copper,L${k + 1},${k === 0 ? 'Top' : k === cu.length - 1 ? 'Bot' : 'Inr'}*%`];
    const aps = new Map<string, number>();
    const ap = (d: number) => { const key = (d * 1e3).toFixed(6); if (!aps.has(key)) { aps.set(key, 10 + aps.size); lines.push(`%ADD${aps.get(key)}C,${key}*%`); } return aps.get(key)!; };
    const net = (n: number) => { if (x2nets) lines.push(`%TO.N,${b.nets.get(n) ?? ''}*%`); };
    for (const t of b.tracks.filter((t) => t.layer === layer)) {
      const d = ap(t.width); net(t.net);
      lines.push(`D${d}*`, `X${f(t.a[0])}Y${f(t.a[1])}D02*`, `X${f(t.b[0])}Y${f(t.b[1])}D01*`);
      if (x2nets) lines.push('%TD*%');
    }
    for (const v of b.vias) {
      const d = ap(v.size); net(v.net);
      lines.push(`D${d}*`, `X${f(v.at[0])}Y${f(v.at[1])}D03*`);
      if (x2nets) lines.push('%TD*%');
    }
    for (const p of b.polys.filter((p) => p.layer === layer)) {
      net(p.net);
      lines.push('G36*', `X${f(p.pts[0][0])}Y${f(p.pts[0][1])}D02*`);
      for (const q of [...p.pts.slice(1), p.pts[0]]) lines.push(`G01X${f(q[0])}Y${f(q[1])}D01*`);
      lines.push('G37*');
      if (x2nets) lines.push('%TD*%');
    }
    lines.push('M02*');
    out.push({ name: `fx-${layer.replace('.', '_')}.gbr`, text: lines.join('\n') });
  });
  const [x0, y0] = b.outline.min, [x1, y1] = b.outline.max;
  out.push({ name: 'fx-Edge_Cuts.gm1', text: ['%FSLAX46Y46*%', '%MOMM*%', '%TF.FileFunction,Profile,NP*%', '%ADD10C,0.100000*%', 'D10*',
    `X${f(x0)}Y${f(y0)}D02*`, `X${f(x1)}Y${f(y0)}D01*`, `X${f(x1)}Y${f(y1)}D01*`, `X${f(x0)}Y${f(y1)}D01*`, `X${f(x0)}Y${f(y0)}D01*`, 'M02*'].join('\n') });
  const tools = [...new Set(b.vias.map((v) => (v.drill * 1e3).toFixed(3)))];
  out.push({ name: 'fx-PTH.drl', text: ['M48', '; DRILL file {KiCad} ; TYPE=PLATED', 'METRIC', ...tools.map((t, i) => `T${i + 1}C${t}`), '%', 'G90', 'G05',
    ...tools.flatMap((t, i) => [`T${i + 1}`, ...b.vias.filter((v) => (v.drill * 1e3).toFixed(3) === t).map((v) => `X${(v.at[0] * 1e3).toFixed(6)}Y${(v.at[1] * 1e3).toFixed(6)}`)]),
    'M30'].join('\n') });
  return out;
}

