/** Keep each release's verification records separate from historical evidence. */
import { mkdirSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

export const root = fileURLToPath(new URL('../', import.meta.url));
export const version = JSON.parse(readFileSync(path.join(root, 'package.json'), 'utf8')).version;
if (!/^\d+\.\d+\.\d+$/.test(version)) throw new Error('Unexpected release version');
export const evidenceDirectory = path.join(root, 'evidence', `release-${version}`);
export const artifactDirectory = path.join(root, 'artifacts', version);
mkdirSync(evidenceDirectory, { recursive: true });
mkdirSync(artifactDirectory, { recursive: true });
