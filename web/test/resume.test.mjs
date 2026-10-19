// Browser test of checkpoint/resume: run the built-in GCPW straight through, then
// start it again, reload the page mid-run (killing the worker after at least one
// IndexedDB checkpoint), start the same run, and check it resumes and gives the
// same S-parameters as the straight run.
import { chromium } from 'playwright';
import { createServer } from 'node:http';
import { readFileSync, existsSync } from 'node:fs';
import { join, dirname, extname } from 'node:path';
import { fileURLToPath } from 'node:url';

const dist = join(dirname(fileURLToPath(import.meta.url)), '..', 'dist');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm' };
const server = createServer((req, res) => {
  const p = join(dist, new URL(req.url, 'http://x').pathname);
  if (!existsSync(p)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'Content-Type': types[extname(p)] || 'application/octet-stream' });
  res.end(readFileSync(p));
});
await new Promise((r) => server.listen(0, r));
const url = `http://localhost:${server.address().port}/index.html`;
const browser = await chromium.launch({ args: ['--enable-unsafe-webgpu'] });
const ctx = await browser.newContext({ locale: 'en-US', serviceWorkers: 'block' });
const page = await ctx.newPage();
page.on('pageerror', (e) => console.log('[pageerror]', e.message));
page.on('console', (m) => { if (/checkpoint/.test(m.text())) console.log('[page]', m.text()); });
page.on('worker', (w) => w.on('console', (m) => console.log('[worker]', m.text())));
const backend = process.env.BACKEND || 'cpu';
const setup = async () => {
  await page.fill('#len', '6'); await page.dispatchEvent('#len', 'change');
  await page.selectOption('#backend', backend);
};
const done = () => page.waitForFunction(() => window.__result || /Error/.test(document.getElementById('status').textContent), null, { timeout: 1800000 });

await page.goto(url); await setup();
await page.click('#run'); await done();
const straight = await page.evaluate(() => window.__result);
console.log(`straight run: ${straight.steps} steps in ${straight.wall.toFixed(1)} s`);

await page.goto(url); await setup();
await page.click('#run');
// checkpoints are written every 20 s of wall time: interrupt after the first one
if (straight.wall < 32) { console.log('run too short to checkpoint before it ends; use a longer line'); process.exit(1); }
await page.waitForTimeout(Math.max(23, 0.55 * straight.wall) * 1000);
const before = await page.textContent('#status');
console.log('interrupting at:', before);
await page.reload(); await setup();
await page.click('#run');
await page.waitForFunction(() => window.__resumed !== undefined || window.__result, null, { timeout: 120000 });
const resumedAt = await page.evaluate(() => window.__resumed);
await done();
const r = await page.evaluate(() => window.__result);
await browser.close(); server.close();
let d = 0;
r.s21.forEach((c, i) => { d = Math.max(d, Math.hypot(c.re - straight.s21[i].re, c.im - straight.s21[i].im) + Math.hypot(r.s11[i].re - straight.s11[i].re, r.s11[i].im - straight.s11[i].im)); });
console.log(`resumed at step ${resumedAt}, finished at ${r.steps} (straight ${straight.steps}); max |dS| vs straight ${d.toExponential(2)}`);
const tol = backend === 'cpu' ? 1e-12 : 1e-5;
const ok = resumedAt > 0 && r.steps === straight.steps && d < tol;
console.log(ok ? 'resume test passed' : 'resume test FAILED');
process.exit(ok ? 0 : 1);
