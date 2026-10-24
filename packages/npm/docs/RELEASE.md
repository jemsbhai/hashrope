# npm 0.3.3 release and provenance

This monorepo package integrates the authoritative standalone source of **published npm hashrope 0.2.2**, rather than an older unpublished worktree. The public registry records 0.2.2 publication on 2026-10-23. Before changes, every one of its 46 tarball files matched the local authoritative source byte for byte.

- Registry: https://registry.npmjs.org/hashrope
- Published source archive: https://registry.npmjs.org/hashrope/-/hashrope-0.2.2.tgz
- Archive SHA-256: `f60c78279264b77a5076d03622961f6b41e17f4ce458a106aaf73e1f7a8580c1`
- Authoritative source location: `E:/data/code/claudecode/hashrope-js`; copied into `packages/npm`, without changing the original.
- Copied-file inventory: `evidence/upstream-0.2.2-provenance.json` (checkout evidence, excluded from npm).

Version 0.3.3 synchronizes the core language packages. `PYTHON_VERSION` deliberately remains 0.2.2 because it identifies the immutable arithmetic oracle, not the current npm version. The default prime/base, polynomial arithmetic, coherent-profile fingerprints, ESM/CJS exports, runtime floor and zero-runtime-dependency contract remain unchanged. Added exact whole/range LCP/equality never use fingerprints. Immutable profile metadata rejects incoherent mixed-profile operations, including no-op paths; this is an intentional compatibility change. The existing `contentEquals` name now uses the exact structural engine.

## Verification and publication

The user has authorized the coordinated package release. This npm integration does not itself publish: the release owner publishes the final reviewed artifact after package and credential checks. At the initial inventory check on 2026-10-24, npm authentication returned E401. The user subsequently completed the security-key login through the coordinated release session, which reported successful login; the release owner then confirmed `npm whoami` as `jemsbhai`. No credential values were printed or placed in source. The release owner verifies current authentication and publishes; this integration has not attempted publication.

```sh
npm ci --ignore-scripts
npm test
npm run pack:verified
npm run test:consumer
npm view hashrope version --registry=https://registry.npmjs.org/
npm whoami --registry=https://registry.npmjs.org/
```

After those gates, the authorized release owner publishes **the inspected tarball**, rather than rebuilding implicitly:

```sh
npm publish artifacts/hashrope-0.3.3.tgz --access public --registry=https://registry.npmjs.org/
```

Confirm the archive checksum against `artifacts/SHA256SUMS` and `evidence/tarball-inspection.json`. After publication, verify version/integrity through registry metadata and install the registry version in fresh ESM/CJS consumers. Publication status belongs to the coordinated release record; this file does not assert that 0.3.3 is already published.

## Contents and omissions

The explicit npm allowlist contains source TypeScript, compiled ESM/CJS, declarations/maps, README, MIT LICENSE and named documentation. Tests, fixtures, upstream wheel, dependencies, logs, preservation manifests, scratch consumers, browser profiles and npm configuration are excluded. `scripts/pack.mjs` inspects actual archive entries and records their hashes. The standalone 0.2.2 historical TDD evidence remains preserved at its source, not presented as newly executed 0.3.3 evidence.
