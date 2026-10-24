# npm 0.3.4 local release validation

This parity release changes package metadata, the exported `VERSION` literal, documentation and release checks. It adds no JavaScript algorithm change or performance claim. Historical 0.3.3 evidence and the independent Python 0.2.2 oracle remain unchanged.

Both Node 20.11.1 and Node 24.22.1 passed all **27,636 source tests**, with no failures, skips, cancellations or todos. Each runtime also passed **14 installed-consumer checks** against the same inspected tarball: fresh ESM/CJS installs, 27,565 published-wheel oracle checks and 11 exact/profile tests through each entry point, all six packaged TypeScript examples, NodeNext import/require declarations, a browser-targeted ES2022 bundle, an isolated realm and real headless Chrome 155.0.8059.85.

The archive is `artifacts/0.3.4/hashrope-0.3.4.tgz`: **53 files**, **68,195 bytes**, SHA-256 `be10002475ccdb22738ff4430380130f30176a60d01a37eae71d4cae5b9ac580`. `tarball-inspection.json` records every packed file hash. `release-summary.json` binds both source suites, both installed-consumer reports and the runtime-preservation audit to that archive. Binaries and raw TAP logs are retained locally in ignored paths; the JSON reports record the relevant hashes and counts.

`runtime-preservation.json` compares all **43 source/build payloads** to the retained published 0.3.3 tarball, SHA-256 `a6e279f45a3665c6f5366a317affadf5ee4a54f41e8b6de2c2f7dc55dd29bc1c`. Thirty-six files are byte-identical; exactly seven index source/output/declaration/source-map files differ only in the `0.3.3` to `0.3.4` version literal. The old ESM files also match the archived 0.3.3 timing record. No timing was rerun or relabeled for this release.

Commands run from `packages/npm`:

```text
npm ci --ignore-scripts --offline --no-audit --no-fund
node scripts/record-suite.mjs node20
<Node 24.22.1> scripts/record-suite.mjs node24
node scripts/audit-oracle-api.mjs
npm run pack:verified
node scripts/consumers.mjs
<Node 24.22.1> scripts/consumers.mjs
python scripts/verify-preservation.py
python scripts/audit-release-runtime.py
python scripts/verify-release.py
```

`node` was `C:/Program Files/nodejs/node.exe` (20.11.1); Node 24.22.1 was the retained official runtime at `E:/data/code/claudecode/hashrope-js/.scratch/runtimes/node-v24.22.1-win-x64/node.exe`. The compiler and bundler remain pinned to TypeScript 5.9.3 and esbuild 0.28.2. Tests ran on Windows; these are finite checks, not certification of every browser or platform.

Publication is a separate coordinated release step. Publish only the inspected tarball after checking its SHA-256, then verify the public registry artifact and a fresh registry consumer. The local gate does not assert publication.
