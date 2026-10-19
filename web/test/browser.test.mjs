// End-to-end browser test: serves web/dist, runs the GCPW demo in headless
// Chromium (Web Worker + wasm), and compares S-parameters with the native app.
// Usage: node web/test/browser.test.mjs <native run dir with sparams.csv>
import { chromium } from 'playwright';
import { createServer } from 'node:http';
import { readFileSync, existsSync } from 'node:fs';
import { join, dirname, extname } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const dist = join(here, '..', 'dist');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm', '.map': 'application/json' };
const server = createServer((req, res) => {
  const p = join(dist, decodeURIComponent(new URL(req.url, 'http://x').pathname));
  if (!existsSync(p)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'Content-Type': types[extname(p)] || 'application/octet-stream' });
  res.end(readFileSync(p));
});
await new Promise((r) => server.listen(0, r));
const url = `http://localhost:${server.address().port}/index.html`;
const browser = await chromium.launch({ executablePath: process.env.CHROME_PATH || undefined, args: ['--enable-unsafe-webgpu'] });
const page = await browser.newPage({ locale: 'en-US', serviceWorkers: 'block' });
page.on('console', (m) => console.log('[page]', m.text()));
await page.goto(url);
const params = JSON.parse(process.env.GCPW_PARAMS || '{}');
for (const [k, v] of Object.entries(params)) {
  if (k === 'ports' || k === 'backend') await page.selectOption('#' + k, String(v));
  else await page.fill('#' + k, String(v));
}
await page.click('#run');
await page.waitForFunction(() => window.__result || /Error/.test(document.getElementById('status').textContent), null, { timeout: 600000 });
const status = await page.textContent('#status');
console.log('status:', status);
const result = await page.evaluate(() => window.__result);
await page.screenshot({ path: join(here, '..', 'dist', 'screenshot.png') });
await browser.close();
server.close();
if (!result) { console.log('browser run FAILED'); process.exit(1); }
if (!process.argv[2]) { console.log('browser run passed (no native run dir given, parity check skipped)'); process.exit(0); }

// native reference: sparams.csv (lumped) or waves.csv (mode ports: S11 = b1/a1, S21 = b2/a1)
let rows;
// Mode-port results go through the 2D eigen solve, whose conditioning at low
// frequency (k0^2 << 1/dx^2) limits platform-to-platform agreement to ~1e-6;
// lumped-port results agree to the CSV's print precision.
const modePorts = existsSync(join(process.argv[2], 'waves.csv'));
const tol = Number(process.env.TOL || (modePorts ? 1e-5 : 1e-9));
if (modePorts) {
  const lines = readFileSync(join(process.argv[2], 'waves.csv'), 'utf8').trim().split('\n');
  const hdr = lines[0].split(',');
  const col = (n) => hdr.indexOf(n);
  rows = lines.slice(1).map((l) => {
    const v = l.split(',').map(Number);
    const a1 = { re: v[col('a1_re')], im: v[col('a1_im')] }, b1 = { re: v[col('b1_re')], im: v[col('b1_im')] };
    const b2 = { re: v[col('b2_re')], im: v[col('b2_im')] };
    const div = (x, y) => { const d = y.re * y.re + y.im * y.im; return { re: (x.re * y.re + x.im * y.im) / d, im: (x.im * y.re - x.re * y.im) / d }; };
    const s11 = div(b1, a1), s21 = div(b2, a1);
    return [v[0], s11.re, s11.im, s21.re, s21.im];
  });
} else {
  rows = readFileSync(join(process.argv[2], 'sparams.csv'), 'utf8').trim().split('\n').slice(1).map((l) => l.split(',').map(Number));
}
const csv = rows;
let worst = 0;
csv.forEach((row, i) => {
  const s21 = result.s21[i], s11 = result.s11[i];
  const d = Math.hypot(s21.re - row[3], s21.im - row[4]) + Math.hypot(s11.re - row[1], s11.im - row[2]);
  worst = Math.max(worst, d);
});
console.log(`browser vs native: max |dS| = ${worst.toExponential(2)} over ${csv.length} frequencies`);
if (worst > tol) { console.log('browser parity FAILED'); process.exit(1); }
console.log('browser parity passed');
