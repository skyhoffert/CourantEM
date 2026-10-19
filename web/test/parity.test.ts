/**
 * Rebuilds solver/apps/parity.c through the TypeScript wrapper and the wasm
 * module, then compares against the native build's output
 * (solver/build/parity_native.txt). Run: node web/test/parity.test.ts
 */
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Axis, Boundary, Courant } from '../src/solver/courant.ts';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..', '..');
const wasm = readFileSync(join(root, 'solver', 'build', 'courant.wasm'));
const native = readFileSync(join(root, 'solver', 'build', 'parity_native.txt'), 'utf8')
  .trim().split('\n').map((l: string) => l.split(/\s+/));

const cem = await Courant.load(wasm);
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
const p2 = s.addLumpedPort(Axis.Z, 19, 9, 0, 5, -1, 50);
s.setGaussian(p1, 1.0, 60e-12, 15e-12);
s.setFrequencies([5e9, 20e9, 40e9]);
const v = s.addVProbe(Axis.Z, 12, 9, 0, 5, -1);
const ip = s.addIProbe(Axis.X, 12, 7, 11, 4, 5);
s.finalize(0.95);
s.step(1500);

const got: [string, number[]][] = [['dt', [s.dt]], ['energy', [s.energy()]]];
const arrs = { V1: s.portV(p1), V2: s.portV(p2), Vp: s.vprobe(v), Ip: s.iprobe(ip) };
for (const [name, a] of Object.entries(arrs)) a.forEach((c, q) => got.push([`${name}[${q}]`, [c.re, c.im]]));
s.destroy();

let worst = 0;
if (got.length !== native.length) throw new Error(`line count ${got.length} vs ${native.length}`);
got.forEach(([name, vals], i) => {
  if (native[i][0] !== name) throw new Error(`label ${name} vs ${native[i][0]}`);
  vals.forEach((val, j) => {
    const ref = Number(native[i][j + 1]);
    worst = Math.max(worst, Math.abs(val - ref) / Math.max(Math.abs(val), Math.abs(ref), 1e-300));
  });
});
console.log(`TS/wasm vs native: max relative difference ${worst.toExponential(2)}`);
if (worst > 1e-12) { console.log('parity FAILED'); process.exit(1); }
console.log('parity passed');
