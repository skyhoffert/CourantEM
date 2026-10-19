// Multithreaded CPU build: serve the app cross-origin isolated (COOP/COEP), run the
// built-in GCPW with 1 thread and with N threads, and require bit-identical
// S-parameters. Prints the speed-up. Run: node test/threads.test.mjs [N]
// With COI=sw the server sends no COOP/COEP headers and isolation must come from
// the bundled coi-serviceworker (as on GitHub Pages), after one automatic reload.
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
  const iso = process.env.COI === 'sw' ? {} : { 'Cross-Origin-Opener-Policy': 'same-origin', 'Cross-Origin-Embedder-Policy': 'require-corp' };
  res.writeHead(200, { 'Content-Type': types[extname(p)] || 'application/octet-stream', ...iso });
  res.end(readFileSync(p));
});
await new Promise((r) => server.listen(0, r));
const N = process.argv[2] || '2';
const browser = await chromium.launch();
const page = await browser.newPage({ locale: 'en-US' });
page.on('pageerror', (e) => console.log('[pageerror]', e.message));
page.on('console', (m) => { if (m.type() === 'error') console.log('[console]', m.text()); });
await page.goto(`http://localhost:${server.address().port}/index.html`);
let iso = false;
for (let k = 0; k < 40 && !iso; k++) {
  iso = await page.evaluate(() => crossOriginIsolated).catch(() => false); // may race the service worker's reload
  if (!iso) await page.waitForTimeout(250);
}
await page.waitForLoadState('load');
console.log('crossOriginIsolated:', iso, process.env.COI === 'sw' ? '(via service worker)' : '(via headers)');
await page.fill('#len', '6'); await page.dispatchEvent('#len', 'change');
const run = async (t) => {
  await page.evaluate(() => { window.__result = undefined; window.__threads = undefined; });
  await page.selectOption('#threads', t);
  await page.click('#run');
  await page.waitForFunction(() => window.__result || /Error/.test(document.getElementById('status').textContent), null, { timeout: 1800000 });
  console.log('status:', await page.textContent('#status'));
  return page.evaluate(() => ({ r: window.__result, threads: window.__threads }));
};
const one = await run('1');
const many = await run(N);
await browser.close(); server.close();
let d = 0;
one.r.s21.forEach((c, i) => { d = Math.max(d, Math.abs(c.re - many.r.s21[i].re), Math.abs(c.im - many.r.s21[i].im),
  Math.abs(one.r.s11[i].re - many.r.s11[i].re), Math.abs(one.r.s11[i].im - many.r.s11[i].im)); });
console.log(`1 thread ${one.r.wall.toFixed(1)} s, ${many.threads} threads ${many.r.wall.toFixed(1)} s (x${(one.r.wall / many.r.wall).toFixed(2)}); max |dS| ${d}`);
const ok = many.threads === Number(N) && d === 0;
console.log(ok ? 'threads test passed' : 'threads test FAILED');
process.exit(ok ? 0 : 1);
