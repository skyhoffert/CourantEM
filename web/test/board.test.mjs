// Browser test of the board path: upload test/fixtures/msl.kicad_pcb in the demo
// page (or, with SOURCE=gerber, its Gerber + drill set and then enter the stackup
// in the stackup editor), run it, and check the result is a matched ~50 ohm line.
import { chromium } from 'playwright';
import { createServer } from 'node:http';
import { readFileSync, existsSync } from 'node:fs';
import { join, dirname, extname } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const dist = join(here, '..', 'dist');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm', '.map': 'application/json' };
const server = createServer((req, res) => {
  const p = join(dist, new URL(req.url, 'http://x').pathname);
  if (!existsSync(p)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'Content-Type': types[extname(p)] || 'application/octet-stream' });
  res.end(readFileSync(p));
});
await new Promise((r) => server.listen(0, r));
const browser = await chromium.launch({ args: ['--enable-unsafe-webgpu'] });
const page = await browser.newPage({ viewport: { width: 1240, height: 980 }, locale: 'en-US', serviceWorkers: 'block' });
page.on('pageerror', (e) => console.log('[pageerror]', e.message));
page.on('console', (m) => { if (m.type() === 'error') console.log('[console]', m.text()); });
await page.goto(`http://localhost:${server.address().port}/index.html`);
const gerber = process.env.SOURCE === 'gerber';
if (gerber) {
  const g = join(here, 'fixtures', 'gerber');
  await page.setInputFiles('#file', ['msl-F_Cu.gbr', 'msl-B_Cu.gbr', 'msl-Edge_Cuts.gm1', 'msl-PTH.drl'].map((f) => join(g, f)));
  await page.waitForFunction(() => document.getElementById('net').options.length > 0);
  console.log('status:', await page.textContent('#status'));
  // Gerbers carry no stackup: the default is 1.6 mm FR-4; enter 254 um RO4350B in the editor
  const diel = (col) => `#stackup tr:nth-child(3) td:nth-child(${col}) input`;
  // fill() fires only `input`; the app listens for `change` (commit), as a user's Enter or blur does
  await page.fill(diel(2), '254'); await page.dispatchEvent(diel(2), 'change');
  await page.fill(diel(3), '3.66'); await page.dispatchEvent(diel(3), 'change');
  await page.waitForFunction(() => /0\.254 mm thick/.test(document.getElementById('status').textContent) &&
    document.querySelector('#stackup tr:nth-child(3) td:nth-child(3) input').value === '3.66');
}
else {
  await page.setInputFiles('#file', join(here, 'fixtures', 'msl.kicad_pcb'));
  await page.waitForFunction(() => document.getElementById('net').options.length > 0);
}
if (process.env.MASK) { // add a top solder mask in the stackup editor
  await page.check('#mask');
  await page.waitForFunction(() => document.querySelector('#stackup tr:nth-child(2) td.name')?.textContent === 'F.Mask');
}
if (process.env.PLANES) { // move both port planes into the board (mm)
  for (const id of ['#offlo', '#offhi']) { await page.fill(id, process.env.PLANES); await page.dispatchEvent(id, 'change'); }
}
console.log('status:', await page.textContent('#status'));
await page.fill('#hmin', '40'); await page.dispatchEvent('#hmin', 'change');
if (process.env.BACKEND) await page.selectOption('#backend', process.env.BACKEND);
if (process.env.FULL2) await page.check('#full2');
await page.click('#run');
await page.waitForFunction(() => window.__result || /Error|failed/.test(document.getElementById('status').textContent), null, { timeout: 1800000 });
console.log('status:', await page.textContent('#status'));
const r = await page.evaluate(() => window.__result);
if (process.argv[2]) await page.screenshot({ path: process.argv[2] });
await browser.close(); server.close();
if (!r) { console.log('board test FAILED'); process.exit(1); }
const db = (c) => 20 * Math.log10(Math.hypot(c.re, c.im));
const s11 = Math.max(...r.s11.slice(1).map(db)), s21 = Math.max(...r.s21.slice(1).map((c) => Math.abs(db(c))));
console.log(`max |S11| (line Z0) ${s11.toFixed(1)} dB, max ||S21|| dev ${s21.toFixed(3)} dB, Z0 ${r.z0[1].toFixed(2)}..${r.z0[r.z0.length - 1].toFixed(2)} ohm`);
let ok = s11 < -35 && s21 < 0.1 && r.z0[1] > 48 && r.z0[1] < 52;
if (process.env.FULL2) {
  // reciprocity (S12 = S21) and the symmetric line (S22 = S11), from two independent runs
  const d = (a, b) => Math.max(...a.map((c, i) => Math.hypot(c.re - b[i].re, c.im - b[i].im)));
  const recip = r.s12 ? d(r.s12, r.s21) : NaN, sym = r.s22 ? d(r.s22, r.s11) : NaN;
  console.log(`full 2-port: max |S12 - S21| ${recip.toExponential(2)}, max |S22 - S11| ${sym.toExponential(2)}`);
  ok = ok && recip < 1e-3 && sym < 3e-3;
}
console.log(ok ? 'board test passed' : 'board test FAILED');
process.exit(ok ? 0 : 1);
