import { WASI } from 'node:wasi';
import { readFileSync } from 'node:fs';
const wasi = new WASI({ version: 'preview1', args: [], env: {} });
const mod = await WebAssembly.compile(readFileSync(process.argv[2]));
const inst = await WebAssembly.instantiate(mod, wasi.getImportObject());
wasi.start(inst);
