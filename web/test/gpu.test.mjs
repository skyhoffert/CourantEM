// Runs test/gpu-test.ts in headless Chromium with WebGPU (SwiftShader on CI).
import { chromium } from 'playwright';
import { createServer } from 'node:http';
import { readFileSync, existsSync } from 'node:fs';
import { join, dirname, extname } from 'node:path';
import { fileURLToPath } from 'node:url';

const dist = join(dirname(fileURLToPath(import.meta.url)), '..', 'dist');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm', '.map': 'application/json' };
const server = createServer((req, res) => {
  const p = join(dist, new URL(req.url, 'http://x').pathname);
  if (!existsSync(p)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'Content-Type': types[extname(p)] || 'application/octet-stream' });
  res.end(readFileSync(p));
});
await new Promise((r) => server.listen(0, r));
const browser = await chromium.launch({ args: ['--enable-unsafe-webgpu'] });
const page = await browser.newPage({ locale: 'en-US', serviceWorkers: 'block' });
page.on('console', (m) => console.log('[page]', m.text()));
await page.goto(`http://localhost:${server.address().port}/gpu.html`);
await page.waitForFunction(() => window.__gpu, null, { timeout: 1800000 });
const r = await page.evaluate(() => window.__gpu);
await browser.close(); server.close();
console.log(JSON.stringify(r));
const ok = !r.error && r.field_rel < 1e-4 && r.s21_db < 0.01 && r.s21_deg < 0.1 && r.lossy_s21_db < 0.01 && r.lossy_s21_deg < 0.1
  && r.resume_s21_db < 0.01 && r.resume_s21_deg < 0.1;
console.log(ok ? 'gpu test passed' : 'gpu test FAILED');
process.exit(ok ? 0 : 1);
