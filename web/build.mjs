// Bundles the app with esbuild into web/dist and copies the solver wasm.
// Usage: node web/build.mjs   (run `make -C solver wasm` first)
import { build } from 'esbuild';
import { copyFileSync, mkdirSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const dist = join(here, 'dist');
mkdirSync(dist, { recursive: true });
await build({
  entryPoints: { main: join(here, 'src/main.ts'), worker: join(here, 'src/solver/worker.ts'), mtworker: join(here, 'src/solver/mtworker.ts'), 'gpu-test': join(here, 'test/gpu-test.ts') },
  outdir: dist, bundle: true, format: 'esm', target: 'es2022', sourcemap: true, minify: true,
  // packages come from node_modules; web/.vendor is a local fallback when npm is unavailable
  nodePaths: [join(here, '.vendor')],
  loader: { '.css': 'text' },
  // CalVer from the repo's VERSION file, stamped into results and the page
  define: { __COURANT_VERSION__: JSON.stringify(readFileSync(join(here, '..', 'VERSION'), 'utf8').trim()) },
  logLevel: 'info',
});
copyFileSync(join(here, 'index.html'), join(dist, 'index.html'));
copyFileSync(join(here, 'public', 'coi-serviceworker.js'), join(dist, 'coi-serviceworker.js'));
copyFileSync(join(here, 'test', 'gpu.html'), join(dist, 'gpu.html'));
copyFileSync(join(here, '..', 'solver', 'build', 'courant.wasm'), join(dist, 'courant.wasm'));
copyFileSync(join(here, '..', 'solver', 'build', 'courant-mt.wasm'), join(dist, 'courant-mt.wasm'));
console.log('built', dist);
