# npm 0.3.4 release and provenance

This monorepo package integrates the authoritative standalone source of **published npm hashrope 0.2.2**, rather than an older unpublished worktree. The public registry records 0.2.2 publication on 2026-10-23. Before changes, every one of its 46 tarball files matched the local authoritative source byte for byte.

- Registry: https://registry.npmjs.org/hashrope
- Published source archive: https://registry.npmjs.org/hashrope/-/hashrope-0.2.2.tgz
- Archive SHA-256: `f60c78279264b77a5076d03622961f6b41e17f4ce458a106aaf73e1f7a8580c1`
- Authoritative source location: `E:/data/code/claudecode/hashrope-js`; copied into `packages/npm`, without changing the original.
- Copied-file inventory: `evidence/upstream-0.2.2-provenance.json` (checkout evidence, excluded from npm).

Version 0.3.3 synchronizes the core language packages. `PYTHON_VERSION` deliberately remains 0.2.2 because it identifies the immutable arithmetic oracle, not the current npm version. The default prime/base, polynomial arithmetic, coherent-profile fingerprints, ESM/CJS exports, runtime floor and zero-runtime-dependency contract remain unchanged. Added exact whole/range LCP/equality never use fingerprints. Immutable profile metadata rejects incoherent mixed-profile operations, including no-op paths; this is an intentional compatibility change. The existing `contentEquals` name now uses the exact structural engine.

## Version 0.3.4

This release keeps all six coordinated package versions aligned. JavaScript algorithms, API contracts, defaults and fingerprints are unchanged from published 0.3.3. Changes are limited to release metadata (including the exported `VERSION` literal), documentation and release checks. `PYTHON_VERSION` remains 0.2.2 and continues to identify the immutable arithmetic oracle.

The [coordinated 0.3.4 release notes](https://github.com/jemsbhai/hashrope/blob/master/releases/0.3.4/README.md) describe the retained Rust-only exact-LCP and concat improvements. Remaining Rust measurements include an approximately 2% pinned LCP cost and mixed 100 KB edit/comparison results; no universal speedup is claimed. The historical npm 0.3.3 timing evidence is preserved, not relabeled as a new 0.3.4 experiment.

`audit-release-runtime.py` verifies the retained 0.3.3 artifact SHA-256 and compares every TypeScript source and compiled ESM/CJS payload. It permits only the 0.3.3-to-0.3.4 version literal in the index source, emitted index files, declarations and source maps. All other runtime bytes must match exactly.

## Verification and publication

The user has authorized the coordinated package release. The release owner publishes the final reviewed artifact after package and credential checks. Local tests and authentication checks do not assert registry publication. Credentials are never stored in package source or evidence.

```sh
npm ci --ignore-scripts
node scripts/record-suite.mjs node20 # repeat on Node 24 with label node24
node scripts/audit-oracle-api.mjs
npm run pack:verified
npm run test:consumer # run on both Node versions
python scripts/verify-preservation.py
python scripts/audit-release-runtime.py
python scripts/verify-release.py
npm view hashrope version --registry=https://registry.npmjs.org/
npm whoami --registry=https://registry.npmjs.org/
```

After those gates, the authorized release owner publishes **the inspected tarball**, rather than rebuilding implicitly:

```sh
npm publish artifacts/0.3.4/hashrope-0.3.4.tgz --access public --registry=https://registry.npmjs.org/
```

Confirm the archive checksum against `artifacts/0.3.4/SHA256SUMS` and `evidence/release-0.3.4/tarball-inspection.json`. After publication, verify version/integrity through registry metadata and install the registry version in fresh ESM/CJS consumers. Publication status belongs to the coordinated release record; this file does not assert that 0.3.4 is already published. Earlier 0.3.3 reports remain unchanged; current reports are written under `evidence/release-0.3.4/`.

## Contents and omissions

The explicit npm allowlist contains source TypeScript, compiled ESM/CJS, declarations/maps, README, MIT LICENSE and named documentation. Tests, fixtures, upstream wheel, dependencies, logs, preservation manifests, scratch consumers, browser profiles and npm configuration are excluded. `scripts/pack.mjs` inspects actual archive entries and records their hashes. The standalone 0.2.2 historical TDD evidence remains preserved at its source, not presented as newly executed 0.3.4 evidence.
