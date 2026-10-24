# npm 0.3.4 verification

Version 0.3.4 is a metadata/parity release. `scripts/audit-release-runtime.py` compares every source and ESM/CJS payload to the retained, hash-verified 0.3.3 tarball; only the exported `VERSION` literal and its declarations/source maps may differ. The runtime algorithms, default profile, fingerprints and independent 0.2.2 oracle are unchanged. The new evidence directory preserves all earlier release evidence. No new performance timing is claimed; the table below remains the historical 0.3.3 measurement. See the [coordinated release notes](https://github.com/jemsbhai/hashrope/blob/master/releases/0.3.4/README.md) for Rust-only changes and remaining performance limits.

The source suite contains **27,636 tests**: the unchanged published-wheel oracle has 27,564 records plus one integrity test, with existing unit/property checks and eleven exact/profile regression tests added in 0.3.3. All 27,636 tests passed on both Node 20.11.1 and Node 24.22.1 with zero failures, skips, cancellations or todos. Results are bound to source hashes in `evidence/release-0.3.4/suite-node20.json` and `evidence/release-0.3.4/suite-node24.json`. The release gate `scripts/verify-release.py` requires both complete suites to pass with zero failures, skips, cancellations or todos, and binds the installed-consumer results to the final tarball. Publication must follow that successful gate.

These retained regressions cover default collisions under both moduli; equal bytes with unequal legacy custom-profile fingerprints; comparison without hash calls or leaf-copy getters; 110-bit repetitions of different periods and shifted phases; late edits next to repeated prefixes; independently allocated 100-level shared DAGs; 250 random trees with 3,000 range byte oracles; strict bounds; immutable profiles; mismatch rejection before no-ops; and streaming profile replacement; and word-scan boundaries with every mismatch position, unaligned starts and tails. Tests compare bytes and exact lengths, rather than treating hash agreement as an oracle.

```sh
npm ci --ignore-scripts
npm test
node scripts/record-suite.mjs node20 # use Node 20.11.1
# repeat with Node 24.22.1 and label node24
node scripts/audit-oracle-api.mjs
npm run pack:verified
npm run test:consumer
# run packed consumers with both Node versions, then:
python scripts/verify-preservation.py
python scripts/audit-release-runtime.py
python scripts/verify-release.py
```

`npm test` compiles both module formats under strict TypeScript and runs every test file. The packed consumer harness installs the inspected tarball into fresh ESM, CJS, declaration and browser consumer directories. It executes all 27,565 published-wheel oracle checks plus the exact/profile regressions separately through installed ESM and CommonJS; compiles and runs all six packaged README/API TypeScript examples; checks NodeNext import/require declarations; bundles for ES2022 without Node runtime imports; and runs an isolated VM and a real headless Chrome/Edge when available. Reports state whether the real-browser check ran.

`evidence/release-0.3.4/tarball-inspection.json`, `artifacts/0.3.4/SHA256SUMS` and `evidence/release-0.3.4/consumer-verification-<node-version>.json` bind results to the actual tarball SHA-256. Source suite success alone is not packed-package verification. Evidence/fixtures are excluded from the npm allowlist.

## Runtime and provenance

The preserved compatibility floor is Node 20.11.1; the additional maintained runtime used by the original implementation is Node 24.22.1. The target is ES2022 with native BigInt and Uint8Array, TypeScript 5.9.3 declarations, and esbuild 0.28.2 browser bundling. Results identify the actual runtime/browser tested; they do not certify every browser or earlier TypeScript releases. ESM and CJS nodes use separate constructor identities and must not be mixed in one graph.

The arithmetic oracle remains the independently downloaded PyPI hashrope 0.2.2 wheel, SHA-256 `710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2`. The checked-in corpus and its manifest were copied unchanged from the verified published npm source. See [ORACLE.md](ORACLE.md) for deterministic regeneration and exclusions, and [RELEASE.md](RELEASE.md) for the npm source provenance.

## Limits

Finite oracle/property tests are not an exhaustive proof. Accepted legacy malformed-input arithmetic and substring behavior remain documented compatibility quirks. Forged nodes, reflection, arbitrary subclass semantics, resource exhaustion and nonterminating Python helper cases are outside the contract. Exact comparisons retain traversal/memoization and BigInt costs; unshared content can require scanning every byte. No universal logarithmic, allocation-free or measured-memory bound is claimed. No publication is implied by a local test result.

## Historical 0.3.3 retained-API timing diagnostic

The published 0.2.2 tarball and frozen 0.3.3 candidate were compared on Node 24.22.1 / Intel Core i9-14900HX in three fresh processes. Each process used six within-process rounds, rotated/reversed fixture order, three AB and three BA orders per fixture, five batches, equal repetition counts and excluded construction/setup. There was no CPU affinity or frequency control. The reported statistic is the median of three process medians (each a median of its six round medians), with process ranges below; these are not confidence intervals.

| Fixture | Published 0.2.2 us [process range] | Candidate 0.3.3 us [process range] |
|---|---:|---:|
| Independent equal 64 KiB leaves | 36.346 [36.080, 37.248] | 16.875 [16.874, 17.357] |
| Late mismatch, 64 KiB | 88.407 [87.120, 89.157] | 16.480 [16.410, 16.861] |
| Early mismatch, 64 KiB | 41.651 [41.309, 42.236] | 0.150 [0.142, 0.153] |
| Shared 64 KiB prefix, changed suffix | 89.804 [88.290, 91.899] | 0.691 [0.658, 0.732] |
| Warm concat of two 128-byte leaves | 0.171 [0.169, 0.178] | 0.219 [0.215, 0.227] |
| Construct 128-byte leaf | 12.940 [12.682, 13.172] | 13.189 [13.093, 13.190] |
| Split a 64 KiB leaf in half | 5835.550 [5662.075, 6021.200] | 5908.975 [5679.850, 6007.075] |

The initial scalar exact helper regressed on independent equal leaves (~2x); its source archive and raw run are retained. The final portable DataView word verifier was selected after that observation and boundary/oracle review. Profile checking adds about 0.048 us to the tested tiny concat (~28.5%); same-root equality moved from about 0.004 to 0.007 us. Tiny timings are particularly sensitive to call/measurement overhead. These fixtures do not establish broad performance parity, zero memory overhead, or universal complexity bounds. Stored profiles add references and query memoization allocates.

Huge 101-bit logical Repeat/range/DAG cases completed with exact expected results in all three final processes. These were single completion checks per case, not comparative throughput samples: the old expanded traversal was not run on astronomical logical inputs. All 13 pre-existing polynomial functions/methods have identical TypeScript-normalized source; the profile constructor/getter additions are separate. Default arithmetic and the oracle corpus remain unchanged.

Raw results and all round ranges are in checkout `evidence/retained-api-timing-final-{0,1,2}.json`, the aggregate in `retained-api-timing-summary.json`, and the readable full table in `RETAINED_API_TIMING.md`. Every raw file records frozen runtime and harness hashes; archive hashes are in `retained-api-timing-manifest-final.json`. The rejected initial run is `retained-api-timing-initial.json`. To reproduce, install/extract the verified 0.2.2 and selected candidate tarballs in separate directories, then run `node scripts/benchmark-release.mjs <baseline-package-dir> <candidate-package-dir> <output.json>` with Node 24.22.1 three times and `python scripts/summarize-release-timing.py` with the documented output filenames. Final documentation-only repacking does not change the timed runtime sources.
