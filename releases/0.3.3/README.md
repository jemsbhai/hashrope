# Synchronized 0.3.3 release

All six packages were published at version **0.3.3** on 2026-10-24. Public
downloads match the reviewed artifacts, and fresh registry consumer checks pass
for every package. The three source repositories are pushed and tagged `v0.3.3`.

| Registry / package | Prior published version | Target |
| --- | --- | --- |
| PyPI `hashrope` | 0.2.2 | 0.3.3 |
| crates.io `hashrope` | 0.3.2 | 0.3.3 |
| npm `hashrope` | 0.2.2 | 0.3.3 |
| Maven Central `io.github.jemsbhai:hashrope` | 0.2.2 | 0.3.3 |
| PyPI `hashrope-bio` | 0.1.0 | 0.3.3 |
| crates.io `hashrope-bio` | 0.2.0 | 0.3.3 |

## Correctness and compatibility

The default base-131 collision remains intentional fingerprint compatibility:
bytes `[0,131]` and `[1,0]` both hash to 263. Exact equality and LCP compare bytes,
not fingerprints. All core backends receive whole-content and checked range
APIs, with profile checks at hash-dependent composition/query boundaries.
Existing applications that compare hashes directly must adopt these APIs.
Biology consumers now verify regions, mutations, codons, genes/exons, trajectory
frames and compound candidates instead of using matching hashes as proof.

The Rust candidate integrates the fixes onto published 0.3.2, retaining its
optional arbitrary-precision backend, bounded-window API, overflow guards,
32-bit cache fix and corrected small-modulus arithmetic. It does not revert to
the older source used in the historical optimization studies. Default hashes
and published conformance fixtures remain unchanged. Language-specific types,
lazy-mode support and constructor shapes remain documented rather than forced
into identical APIs.

## Efficiency and limits

Exact queries avoid flattening, splitting and rehashing. Shared ranges and
repetition can be skipped after byte/structural proof; verified DAG ranges are
memoized. Unshared equal bytes can still require a full scan. Arbitrary-precision
arithmetic and caches allocate; no universal logarithmic exact-query bound is
claimed. Profile checks deliberately reject previously unsafe parameter mixing.

The [round-four investigation](../../benchmarks/exact_optimization_round4/README.md)
and all earlier evidence are retained unchanged. Its short-protocol equality
regressions remain unresolved, despite near-unity estimates in separate
controlled diagnostics. Those historical measurements do not establish the
timing of this newly integrated release. Release-specific checks are recorded
separately; parity is not a promise that every language or workload has equal
speed.

The [release integration diagnostic](performance/rust_integration/README.md)
compares the integrated Rust source against the completed round-four exact
implementation, using unchanged workloads in three paired processes. Independent
equal 1 MiB LCP is 10.5% slower by paired ratios, and compressed 2 TiB DAG LCP is
30.8% slower (displayed medians 0.252 to 0.340 microseconds). Both slow in all
three pairs. Independent equality and late LCP are about 1.7% faster; exact churn
changes range from 1.0% to 8.0% slower. The unchanged Ropey controls also vary,
but this does not prove that Hashrope's regressions are noise. Small independent
edits remain slower than Ropey. No universal performance parity is claimed.

The [npm retained-API diagnostic](../../packages/npm/evidence/RETAINED_API_TIMING.md)
compares against published 0.2.2. After optimizing the byte-exact leaf scan, equal
64 KiB comparisons take 16.9 versus 36.3 microseconds in the measured workload.
Profile-checked concatenation adds about 0.05 microseconds. Full samples, source
identities, the rejected scalar implementation, and limits are retained.

## Release validation

- Built Python core and biology wheels/sdists. Installed-wheel tests: **244
  passed**. Installed collision/range/profile/symbolic-repeat release gate passed;
  all four distribution metadata checks passed.
- Merged Rust fixed-width checks: **169 release tests**, **20 debug checks**,
  and release `no_std` build passed.
- Rust bigint: **12 new exact-comparison tests**, **144 library tests**,
  **6 published Python conformance tests**, one doctest, and release
  `no_std + bigint` build passed. These overlap the fixed-width suite and are
  not summed as distinct tests.
- Rust: formatting, strict all-feature/all-target Clippy, warnings-as-errors
  docs, both wasm32 configurations, the actual 32-bit numeric runtime probe,
  and default/all-feature package verification passed. The focused final gate
  has **198 tests**, overlapping the counts above. Large E3/E5 experiments are
  excluded from this local gate.
- Rust biology: library, exact-consumer and benchmark-helper regressions plus
  compile-only doctest passed.
- npm: strict build and **27,636 tests** passed on Node 20 and 24, including pinned published
  Python-oracle fixtures and new exact/profile/range/compressed regressions.
- Java: **10,618 library tests**, **3 installed-JAR consumer tests**, and **26
  packaging/validator tests** passed. Published fixture bytes are unchanged.
  Its retained hash/construction/split/sliding timings are close to published
  0.2.2, with observed variation and no zero-overhead claim. See the
  [Java six-pair diagnostic](https://github.com/jemsbhai/hashrope-java/blob/main/evidence/release-0.3.3/EFFICIENCY.md).
- npm's final packed consumers pass on both tested Node versions, including real
  Chrome. Java's signed bundle passes all 24 file checks and four signature
  checks. Fresh public registry checks are recorded below.

All **14 CI jobs passed**: Python installed wheels on 3.10/3.13; npm source and
packed consumers on Node 20/24; Rust on Linux/macOS/Windows, published reference
fixtures, wasm32 and downstream compatibility; Java on 17/24. Unlike the focused
local Rust gate, the CI native jobs ran the complete release suite including
E3/E5. See [CI identities and results](validation/ci.json).

## Publication

| Package | Public verification |
| --- | --- |
| [Python core](https://pypi.org/project/hashrope/0.3.3/) | Passed: both public distributions match reviewed hashes; fresh PyPI install, exact API gate and combined 244-test suite. |
| [Python biology](https://pypi.org/project/hashrope-bio/0.3.3/) | Passed: both public distributions match reviewed hashes; fresh PyPI install and combined 244-test suite. |
| [Rust core](https://crates.io/crates/hashrope/0.3.3) | Passed: uploaded archive checksum, empty-cache registry resolution, exact/profile/range and bigint tests. |
| [Rust biology](https://crates.io/crates/hashrope-bio/0.3.3) | Passed: uploaded archive checksum, empty-cache registry resolution and exact consumer tests. |
| [npm core](https://www.npmjs.com/package/hashrope/v/0.3.3) | Passed: canonical tarball checksum, normal empty-cache registry installation, all 53 installed files, ESM/CJS exact/profile/compressed smoke tests. |
| [Java core](https://central.sonatype.com/artifact/io.github.jemsbhai/hashrope/0.3.3) | Passed: all 24 public files match signed bundle, four signatures verify, empty-cache Central-only Maven resolution and all three consumer tests pass. |

Python evidence: [public checks](validation/python-public-install.json),
[installed tests](validation/python-public-install.log), and
[reviewed artifacts](validation/python-artifacts.json).
Rust evidence: [registry consumer audit](validation/rust-public-registry/README.md).
The uploaded core crate differs from the earlier precommit package only in
`.cargo_vcs_info.json`: it records the clean release commit. Every other archive
member, including production source, tests, docs and manifests, is identical.

npm evidence: [normal public install](validation/npm-normal-public-install.json)
and [public archive verification](validation/npm-public-registry.json).
Java evidence: [public Central audit](https://github.com/jemsbhai/hashrope-java/tree/main/evidence/release-0.3.3/public-central).
Initial CDN propagation failures and verification-parser corrections are retained;
they are not omitted or represented as passing attempts.

## Published source

| Repository | Source tag | Reviewed release commit |
| --- | --- | --- |
| Core Python/Rust/npm | [v0.3.3](https://github.com/jemsbhai/hashrope/tree/v0.3.3) | `3dda4b5940335785976442c16e8743934fef8092` (merged as `4e2836841748fbc902fb3c9975a87a3a37b7f0af`) |
| Python/Rust biology | [v0.3.3](https://github.com/jemsbhai/hashrope-bio/tree/v0.3.3) | `3d81bc84ea1417837102fafc09407e3079e60651` |
| Java | [v0.3.3](https://github.com/jemsbhai/hashrope-java/tree/v0.3.3) | `ddef675f51c30c725bf36a4cfe1877aa19ad4409` (merged as `bcb06158a27f6896248e8cacc136b03c78f6d952`) |

Later documentation-only commits add public verification records; package source
and fingerprints remain the tested, published release contents.
