/** Touchstone (.s2p) export of a 2-port result, RI format, 50-ohm reference.
 * With both ports excited (full 2-port) all four S-parameters are written;
 * otherwise the line is assumed symmetric and reciprocal (S22 = S11, S12 = S21),
 * and the header says so. */
import type { Result } from './store.ts';

export function toTouchstone(r: Result, title: string): string {
  const s11 = r.s11_50 ?? r.s11, s21 = r.s21_50 ?? r.s21;
  const full = !!(r.s12 && r.s22);
  const s12 = (r.s11_50 ? r.s12_50 : r.s12) ?? s21, s22 = (r.s11_50 ? r.s22_50 : r.s22) ?? s11;
  const ref = r.s11_50 ? '50' : 'line Z0 (not renormalised)';
  const lines = [
    `! Courant EM ${typeof __COURANT_VERSION__ === 'string' ? __COURANT_VERSION__ : 'dev'} — ${title}`,
    `! ${new Date().toISOString()}; reference ${ref}; ${full ? 'both ports excited' : 'port 1 excited, S22 = S11 and S12 = S21 assumed'}`,
    '# HZ S RI R 50',
  ];
  r.f.forEach((f, i) => {
    const a = s11[i], b = s21[i], c = s12[i], d = s22[i];
    // Touchstone 2-port order: S11 S21 S12 S22
    lines.push([f, a.re, a.im, b.re, b.im, c.re, c.im, d.re, d.im].map((v, k) => (k ? v.toExponential(9) : v.toExponential(6))).join(' '));
  });
  return lines.join('\n') + '\n';
}
