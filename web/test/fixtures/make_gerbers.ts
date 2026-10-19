/** Writes test/fixtures/gerber/msl-*.{gbr,gm1,drl} from msl.kicad_pcb (run: node test/fixtures/make_gerbers.ts). */
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { importKicad } from '../../src/import/kicad.ts';
import { toGerbers } from './gerberwriter.ts';

const here = dirname(fileURLToPath(import.meta.url));
const b = importKicad(readFileSync(join(here, 'msl.kicad_pcb'), 'utf8'));
mkdirSync(join(here, 'gerber'), { recursive: true });
for (const f of toGerbers(b, true)) writeFileSync(join(here, 'gerber', f.name.replace(/^fx-/, 'msl-')), f.text);
console.log('wrote', toGerbers(b, true).length, 'files');
