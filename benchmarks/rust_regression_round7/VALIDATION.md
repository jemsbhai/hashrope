# Round-seven validation

**Both the rejected candidate and the restored production implementation pass
all seven checks below, including 221 release tests each.** The complete
commands, exit codes, timestamps, hashes and logs are recorded separately:

- Trial: [validation/results.json](validation/results.json), with tested Rust
  sources, tests, fixtures, lockfile and driver in
  [validation/sources.zip](validation/sources.zip).
- Retained round-six core, freshly checked after restoration:
  [retained_validation/results.json](retained_validation/results.json), with
  [retained_validation/sources.zip](retained_validation/sources.zip).

| Check | Result in both implementations |
|---|---|
| `cargo fmt --all --check` | Pass |
| `cargo clippy --all-features --all-targets -- -D warnings` | Pass |
| `cargo test --all-features --release` | **221 tests pass** |
| `cargo test --doc` | Pass; zero doctests |
| `cargo check --lib --no-default-features --target wasm32-unknown-unknown` | Pass |
| `cargo check --lib --all-features --target wasm32-unknown-unknown` | Pass |
| `cargo doc --no-deps --all-features`, `RUSTDOCFLAGS=-D warnings` | Pass |

The existing suite covers the default collision, cross-profile and lazy exact
comparisons, repeated/shared DAGs, range/visit-budget fallbacks, scalar-reference
fingerprints, persistent edit traces, full node-sequence goldens, and overflow
before allocation. No new behavior or public API is introduced by the private
heavier-side argument, so no redundant implementation-mirroring test was added.
The optional bigint implementation is unchanged and remains included in the
all-feature tests; it is not covered by the default-feature timing workloads.

The trial's focused preflight also passed 27 release tests and two debug tests;
see [codegen/README.md](codegen/README.md). Trial validation ran before the final
benchmark builds and timings. Restored-core validation ran after the rejection.
Neither overlapped with benchmark timings or agent builds. Both full validation
records verify that their source snapshot remained unchanged throughout checks.

The candidate changes only the private join helper and its two callers from the
frozen round-six runtime: pass the already computed heavier-side predicate
instead of loading both node weights again. The predicate, overflow/profile
guard order, recursive operations, hash calculations, allocation sequence,
defaults, node layout and exact comparison code remain unchanged. Source review
and tests support those invariants; they do not establish a timing improvement.

The rejected candidate `rope.rs` SHA-256 is
`5931fcb08942a8ea887bb38ec25f57e651665583a39a42142af090f687f231d6`.
The final measurement [manifest](final_comparison/manifest.json) freezes the
same trial source. Its benefit was not established, and secondary costs remained;
see [decision.json](decision.json). Restored `rope.rs` has SHA-256
`50d9441b03d5e6f4acba329bcd449fe5dbc4354b6f09c2d9c28fb54cc916d4e7`.
All six core files byte-match frozen round six, independently verified in
[final_review.json](final_review.json) and
[retained_source_identity.json](retained_source_identity.json).
The trial snapshots and measurements remain intact. Package versions remain
0.3.3 and no release is performed by this study.
