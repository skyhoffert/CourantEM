/** Mesher rules on the milestone-1 GCPW cross-section. Run: node web/test/mesher.test.ts */
import { estimate, lambdaOver20, meshAxis } from '../src/mesh/mesher.ts';

let fail = 0;
const check = (ok: boolean, msg: string) => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${msg}`); if (!ok) fail = 1; };

const h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6, er = 3.66, fmax = 40e9;
const a = w / 2, b = a + s, yf = b + fence, yext = yf + 1e-3;
const hmin = 20e-6;
const sub = lambdaOver20(fmax, er);

// y: strip edges (metal inside |y| < a) and ground edges (metal outside |y| > b)
const Y = meshAxis({
  lo: -yext, hi: yext,
  sheetEdges: [{ x: -b, metal: 'lo' }, { x: -a, metal: 'hi' }, { x: a, metal: 'lo' }, { x: b, metal: 'hi' }],
  planes: [-yf, 0, yf],
  spans: [{ a: -b, b: -a, n: 8 }, { a, b, n: 8 }],
  refine: [-b, -a, a, b],
}, { hmin, hmax: 100e-6 });
const Z = meshAxis({
  lo: 0, hi: h + 1.5e-3, planes: [h], spans: [{ a: 0, b: h, n: 6 }], maxCell: [{ a: 0, b: h, h: sub }], refine: [h],
}, { hmin, hmax: 100e-6 });

for (const [nm, m] of [['y', Y], ['z', Z]] as const) {
  console.log(`${nm}: ${m.report.cells} cells, min ${(m.report.minCell * 1e6).toFixed(1)} µm, max ${(m.report.maxCell * 1e6).toFixed(1)} µm, ratio ${m.report.maxRatio.toFixed(3)}`);
  check(m.report.maxRatio <= 1.3 + 1e-6, `${nm}: neighbour ratio <= 1.3`);
  check(m.report.warnings.length === 0, `${nm}: no warnings ${m.report.warnings.join('; ')}`);
  for (let i = 0; i + 1 < m.nodes.length; i++) if (!(m.nodes[i + 1] > m.nodes[i])) { check(false, `${nm}: strictly increasing`); break; }
}
// thirds rule lines present at every sheet edge
const has = (arr: number[], x: number) => arr.some((v) => Math.abs(v - x) < 1e-12);
for (const [x, side] of [[a, -1], [-a, 1], [b, 1], [-b, -1]] as const) {
  check(has(Y.nodes, x + side * hmin / 3) && has(Y.nodes, x - side * 2 * hmin / 3), `thirds-rule lines at edge ${x * 1e6} µm`);
  check(!has(Y.nodes, x), `no line exactly on the sheet edge ${x * 1e6} µm`);
}
check(has(Z.nodes, h), 'substrate top on a grid line');
check(has(Y.nodes, 0) && has(Y.nodes, yf) && has(Y.nodes, -yf), 'centre and fence planes on grid lines');
const zsub = Z.nodes.filter((v) => v > 0 && v < h).length + 1;
check(zsub >= 6, `>= 6 cells through the substrate (${zsub})`);
check(Z.nodes.every((v, i) => i === 0 || v > h || v - Z.nodes[i - 1] <= sub * (1 + 1e-9)), 'substrate cells <= lambda/20');

const X = meshAxis({ lo: 0, hi: 10e-3, maxCell: [{ a: 0, b: 10e-3, h: 50e-6 }] }, { hmin: 50e-6, hmax: 50e-6 });
const est = estimate(X.nodes, Y.nodes, Z.nodes, { tSim: 0.5e-9, backend: 'gpu' });
console.log(`estimate: ${est.cells} cells, dt ${(est.dt * 1e15).toFixed(1)} fs, ${est.steps} steps, ${(est.bytes / 1e6).toFixed(0)} MB on GPU`);
check(est.cells > 0 && est.dt > 0, 'estimator runs');
console.log(fail ? 'mesher test FAILED' : 'mesher test passed');
process.exit(fail);
