# Synchronized 0.3.4 release

Version **0.3.4** synchronizes the four core language packages and the Python
and Rust `hashrope-bio` packages. This patch includes the accepted Rust LCP and
concat optimizations. Python, JavaScript/TypeScript, Java and biology algorithms
retain their 0.3.3 behavior; their version and dependency metadata advance to
keep the release family aligned. The rejected round-seven join trial is excluded.

The final release checks and all 14 CI jobs passed. Publication and independent
public-download verification are recorded below as each registry completes.

| Registry / package | Version | Publication and public verification |
|---|---|---|
| [PyPI `hashrope`](https://pypi.org/project/hashrope/0.3.4/) | 0.3.4 | Published; both public distribution hashes, fresh PyPI install, exact API gate and all 244 combined tests passed. |
| [crates.io `hashrope`](https://crates.io/crates/hashrope/0.3.4) | 0.3.4 | Published; exact public archive hash and fresh registry-only consumer checks passed. |
| [npm `hashrope`](https://www.npmjs.com/package/hashrope/v/0.3.4) | 0.3.4 | Published; canonical archive, all 53 installed files and fresh ESM/CJS consumer checks passed. |
| Maven Central `io.github.jemsbhai:hashrope` | 0.3.4 | Pending |
| [PyPI `hashrope-bio`](https://pypi.org/project/hashrope-bio/0.3.4/) | 0.3.4 | Published; both public distribution hashes, fresh PyPI install and all 244 combined tests passed. |
| [crates.io `hashrope-bio`](https://crates.io/crates/hashrope-bio/0.3.4) | 0.3.4 | Published; exact public archive hash and fresh registry-only consumer checks passed. |

## Correctness and compatibility

Exact byte equality and LCP remain available independently of fingerprints,
including cross-profile and compressed inputs, and lazy inputs where the
backend supports them. Hash-profile checks,
overflow checks, fingerprint arithmetic/defaults, public APIs and node layout
are unchanged. The compatible default collision `[0,131]` versus `[1,0]`
still hashes to 263; exact comparison distinguishes those byte strings.
Applications must continue to use the exact APIs when correctness depends on
equality. This release does not make fingerprints collision-free.

The fixed-width Rust runtime matches the accepted round-six source inspected
and validated in the [round-seven study](../../benchmarks/rust_regression_round7/README.md).
Its `rope.rs` SHA-256 is
`50d9441b03d5e6f4acba329bcd449fe5dbc4354b6f09c2d9c28fb54cc916d4e7`.
Release metadata changes to 0.3.4; the rejected direction-argument trial is not
part of this source. The optional bigint engine is unchanged.

## Performance improvements and remaining costs

The accepted Rust changes compare aligned child ranges directly and verify
equal-period compressed ranges at their structurally proven period roots.
Balanced concat reuses checked metadata. Every shortcut retains exact byte or
structural proof; fingerprints are not accepted as equality evidence.

The [round-five study](../../benchmarks/rust_lcp_round5/README.md), against
published 0.3.3, measured 17.0%/13.0% lower independent-equal LCP time at
64 KiB/1 MiB and 58.3% lower self-concat 2 TiB DAG LCP time in its nine-pair
cohort. The [round-six study](../../benchmarks/rust_regression_round6/README.md)
measured 14.5% and 12.2% lower balanced-concat time in separate cohorts. That
concat gain repeats at 15.2% in the latest retained-binary replication and
13.6–15.5% in controlled comparisons. These are fixture-specific results from
separate comparisons, not percentages to add or a new universal release estimate.

**The remaining comparison costs are not dismissed as noise.** The earlier
4.4% independent 100 KB edit/comparison and 6.8% independent 1 MiB LCP slowdown
magnitudes do not repeat consistently. Twelve fresh pairs with the original
executables estimate +2.0% and −0.7%, respectively, with mixed directions and
intervals including one. However, controlled pinned 1 MiB LCP still shows
approximately **+2% in both implementation-slot orders**; one primary interval
is wholly above one. A secondary pinned comparison-only timed mean is +5.2%,
while its median estimate and other placements differ. Their cause remains
unresolved. A/A variation and unchanged comparator instructions do not erase
these observations.

A further join-weight-lookup trial passed correctness checks but did not
establish useful improvement in either target, and some early-mismatch cases
slowed. It was rejected, with all measurements retained. The final retained
core still takes about **1.77× Ropey's time at independent 100 KB and 1.31× at
independent 1 MB**, and is farther behind at 10 KB in the measured ASCII edit
fixtures. Shared-history 1 MB comparisons favor Hashrope in that fixture.
These ratios are process-paired estimates, not universal API comparisons.

All 66 benchmark cells, order/seed effects, uncertainty and rejected-trial data
remain in the [complete investigation](../../benchmarks/rust_regression_round7/README.md).
The measurements concern Windows/default-feature Rust; they do not establish
bigint, other-platform or other-language speedups. Exact comparison can scan
independently stored equal bytes and has no general logarithmic time bound.
Correctness and synchronized version numbers do not imply performance parity.

## Final release validation

Local release gates pass:

- Rust: **221 release tests**, strict all-feature/all-target Clippy, formatting,
  documentation and both wasm32 compilation checks; see
  [final Rust validation](validation/rust-final/results.json).
- Python: both wheels built from sdists, four strict distribution metadata
  checks, **244 installed-package tests** and the exact/profile/range/repeat
  gate; see [artifact evidence](validation/python-artifacts.json).
- Rust biology: **15 tests** against the local 0.3.4 core. The platform CI also
  exercises the current committed biology consumer alongside the older pinned
  compatibility fixture.
- npm: **27,636 source tests** and **14 packed-consumer checks** on each of
  Node 20.11.1 and 24.22.1, including ESM/CJS oracle and exact checks,
  declarations, examples and real Chrome. Runtime preservation permits only
  version metadata changes; see
  [npm release evidence](../../packages/npm/evidence/release-0.3.4/README.md).
- Java: source, installed-JAR, packaging and bundle/signature gates are recorded
  in the [Java release evidence](https://github.com/jemsbhai/hashrope-java/tree/codex/release-0.3.4/evidence/release-0.3.4).

The first local Rust capture saw Cargo update the package version in its ignored
lockfile from 0.3.3 to 0.3.4. All commands passed, but its source-unchanged audit
correctly failed; that attempt is retained in [rust](validation/rust/). The
complete final gate was repeated with the settled lockfile and passes both
command and source-identity checks. No failed evidence is discarded.

All **14 CI jobs passed** on the reviewed release commits: Rust native and
biology on Linux/macOS/Windows, wasm32 execution and published reference
fixtures; Python installed wheels on 3.10/3.13; npm source and packed consumers
on Node 20/24; Java on 17/24. The native Rust jobs ran the complete release suite,
including E3/E5. See [CI identities and results](validation/ci.json).

Historical benchmark sources and raw results remain unchanged.

## Public verification and source

Both Python packages' public wheels and sdists match the reviewed artifacts.
A normal version-pinned PyPI installation into a fresh environment passed
dependency checks, module-origin/version assertions, the exact API release gate
and all 244 combined tests. See [public installation evidence](validation/python-public-install.json)
and [test log](validation/python-public-install.log). Initial metadata and index
propagation failures remain recorded; they were retried with unchanged
verification logic and new empty environments, without republishing.

Rust's clean package and publication dry-run checks passed before upload. Both
public crate archives match the reviewed clean artifacts byte for byte. A fresh
Cargo home resolved both exact 0.3.4 versions solely from crates.io; five default
and seven all-feature consumer tests passed. See the
[Rust public-registry audit](validation/rust-public-registry/README.md).

npm's canonical public archive matches the reviewed tarball. A normal registry
install with empty cache and blank user/global configuration passed all 53 file
checks and ESM/CJS exact/profile/compressed-input checks. See
[public verification](../../packages/npm/evidence/release-0.3.4/public-attempt-2/verification.json).
The initial metadata 404 from registry propagation is retained as
[attempt one](../../packages/npm/evidence/release-0.3.4/public-attempt-1/verification.json).

All three source PRs are merged. Tags identify the reviewed source commits;
later documentation-only commits add publication evidence without changing
package runtime contents.

| Repository | Source tag | Reviewed release commit | Merge commit |
|---|---|---|---|
| Core Python/Rust/npm | [v0.3.4](https://github.com/jemsbhai/hashrope/tree/v0.3.4) | `3f730fc13ff882822d91d3277736f7bde12d487c` | `b2a4c83792f54ba9fd12f0979aba8e353b8cc061` |
| Python/Rust biology | [v0.3.4](https://github.com/jemsbhai/hashrope-bio/tree/v0.3.4) | `8e07e641af532849a629bbb595c18afcfd7ab6ab` | `f48e3540fc1e569cbf06e26db756ec3020332f96` |
| Java | [v0.3.4](https://github.com/jemsbhai/hashrope-java/tree/v0.3.4) | `82e771bdae1f6293f7328bc4ee488e6317abe8c0` | `aa0708905ee629654a4d38654146cf3a8b1d7b42` |

The [source release record](validation/source-release.json) records the exact
heads protected by each merge request. The independent
[Rust publication audit](validation/rust-publication-independent-review.json)
confirms artifact, log and tag identities, and lists the few LF/CRLF-only
differences between validated Windows checkout files and tagged Git blobs.

The [0.3.3 release record](../0.3.3/README.md) preserves the earlier publication
and validation history. It is not substituted for the checks of these 0.3.4
artifacts.
