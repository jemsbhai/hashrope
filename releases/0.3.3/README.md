# Synchronized 0.3.3 release

Release preparation is in progress. This record does not assert registry
publication until the package-specific public verification is recorded below.

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

## Candidate validation

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
  0.2.2, with observed variation and no zero-overhead claim. See the Java
  repository's `evidence/release-0.3.3/EFFICIENCY.md` for the six-pair diagnostic.
- Final packaged consumers, signing and registry verification are recorded
  below as they complete.

## Publication

Pending. Version numbers, built artifacts, credentials or upload attempts are
not evidence of public availability. Each package requires fresh registry
resolution and an installed-artifact smoke check after publication.
