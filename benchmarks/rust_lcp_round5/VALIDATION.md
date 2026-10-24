# Candidate correctness and compatibility

The final candidate passed all seven checks in [validation/results.json](validation/results.json)
on Windows x86-64 with Rust 1.94.0. The record contains source/test hashes,
commands, exit codes, and hashes of the complete logs. Its source hashes remained
unchanged throughout validation and match the subsequently frozen timing candidate.
The test dependency lock is retained as
[validation/validation.Cargo.lock](validation/validation.Cargo.lock), matching
the recorded `packages/rust/Cargo.lock` hash.

| Check | Result |
|---|---|
| `cargo fmt --all --check` | Passed |
| `cargo clippy --all-features --all-targets -- -D warnings` | Passed |
| `cargo test --all-features --release` | 219 passed; none failed or ignored, including all E1–E7 experiments and the bigint doctest |
| `cargo test --doc` | Passed; no default-feature doctests defined |
| `cargo check --lib --no-default-features --target wasm32-unknown-unknown` | Passed |
| `cargo check --lib --all-features --target wasm32-unknown-unknown` | Passed |
| `cargo doc --no-deps --all-features` with `RUSTDOCFLAGS=-D warnings` | Passed |

The new five integration regressions cover all phases of several equal periods,
short and clipped ranges, a last-witness-byte mismatch after wrapping, three
phase pieces over enormous nonperiodic compressed units, aligned and mismatched
tree shapes, and shared-prefix/repeat fallback. Two private tests cover odd/even
visit-budget exhaustion and resumption after a verified prefix, and verify that
ordinary period-root witnesses leave the general cache empty.

Existing tests retain the base-131 collision, cross-arena numeric-ID collisions,
different hash profiles and hasher replacement, lazy hashes, near-`u64::MAX`
ranges, non-adjacent shared DAGs, published Python/Pollard fingerprint fixtures,
legacy arithmetic/edit traces, and bigint exact comparisons. The separate
[mechanism probe](mechanism_audit.md) passes sixteen additional fixture/oracle
checks and measures allocation requests without interpreting instrumented time.

Only the fixed-width exact-comparison implementation changed. Arithmetic,
defaults, public APIs, constructors, node metadata layout, the bigint engine,
and package versions are unchanged. The Rust README's bigint dependency example
now names the already published 0.3.3 exact-API release.

This is local validation and wasm compilation, not a fresh Linux/macOS CI run or
a wasm runtime benchmark. The E1–E7 tests are included as the existing suite;
their internal timing output is not used as paired performance evidence here.
The optimization remains unreleased. Reproduce the validation in a fresh study
directory with `python benchmarks/rust_lcp_round5/validate.py`; the script refuses
to overwrite its evidence.
