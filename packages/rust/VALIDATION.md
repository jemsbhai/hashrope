# Rust parity validation

Executed locally on Windows x86_64 with Cargo/Rust 1.94.0, including actual
wasm32 execution under Node. Branch `codex/rust-python-parity`, based on
`100749d1cef9603465f4600b814c76d2ad6e9a6e`. The compatible release version is
`0.3.2`, preserving existing `^0.3` consumers. The baseline local checks below
began on the implementation snapshot before the release version was assigned.

## Crate checks

| Check | Result |
| --- | --- |
| `cargo fmt --all --check` and `git diff --check` | Pass |
| `cargo clippy --all-features --all-targets -- -D warnings` | Pass |
| `cargo test --no-default-features --lib --test python_conformance` | 120 library tests and 4 conformance tests pass |
| `cargo test --all-features --lib --test python_conformance` | 137 library tests and 6 conformance tests pass |
| `cargo test --all-features --release` | Complete release suite passed: 136 library tests at that snapshot, 17 integration/experiment tests, 1 doctest |
| `cargo test --all-features --release --lib --test python_conformance` after deep-traversal repair | Final 137 library tests and 6 conformance tests pass; unchanged experiment targets had already passed |
| `cargo test --doc --all-features`, with `RUSTDOCFLAGS=-D warnings` | Pass |
| `cargo doc --no-deps --all-features`, with `RUSTDOCFLAGS=-D warnings` | Pass; default-feature documentation also passes |
| `cargo build --lib --no-default-features --target wasm32-unknown-unknown` | Pass |
| `cargo build --lib --all-features --target wasm32-unknown-unknown` | Pass after final implementation changes |
| `cargo tree --no-default-features --edges normal` | Only hashrope itself: no default runtime dependencies |
| `cargo package --allow-dirty --no-default-features` | Archive generated and compiled successfully |
| `cargo package --allow-dirty --all-features` | Archive generated and compiled successfully; source, scripts, golden fixtures and documentation included |
| Published-wheel fixture regeneration/checks | Both Python conformance JSON and Pollard JSONL reproduce exactly |

The release experiment tests cover substring scaling, hashing throughput,
concat/split, compressed repetition, geometric sums, sliding throughput,
memory and tree height. Criterion's full benchmark sweep was not run; passing
experiment tests is not a new performance guarantee.

An independent review found a stack overflow when validating a split of a
1024-bit repetition count. The bigint API now keeps traversal/continuation
frames on the heap. A permanent regression exercises deep split, validation,
hashing, rejoin and height. Independent post-fix probes passed 6,000 mixed
operations, 36,000 substring queries across twelve modulus widths, and 12,000
bounded-stream operations. Those probes were separate from the checked-in
golden and unit suites, not a claim of exhaustive formal verification.

## Reproduce the 32-bit runtime checks

From this crate directory, with the wasm32 Rust target and Node installed:

```text
rustc --crate-type cdylib --edition 2021 --target wasm32-unknown-unknown scripts/wasm_numeric_probe.rs -o target/numeric_probe.wasm
node scripts/verify_wasm.cjs target/numeric_probe.wasm
```

This executes full u64 power exponents on a 32-bit target, small-modulus
reduction, `u64::MAX` compressed rope hashing, and deterministic traps for
unaddressable materialization and overflowing repetition. Separate Cargo
builds above verify both `no_std + alloc` library configurations.

## Reproduce published-package conformance

Download/extract the exact wheels named in [PARITY.md](PARITY.md), put their
extracted roots on `PYTHONPATH`, and run:

```text
python scripts/generate_python_vectors.py --wheel <hashrope-0.2.2-wheel> --check
python scripts/verify_pollard_fixture.py --hashrope-wheel <hashrope-0.2.2-wheel> --pollard-wheel <pollard-1.6.0-wheel>
cargo test --all-features --test python_conformance
```

Both Python scripts verify archive SHA-256 and imported source against the
archives before checking outputs. The second script uses Pollard entirely in
memory, without provider calls or writing a store. Hashrope and Pollard are
reference test dependencies only; Rust never executes Python at runtime.

## Downstream checks and adoption limits

- **cdh-sort**: local revision `6d40877ba7675e2a3b41dda5ed524b7644ce2ddb`,
  including preexisting local edits, copied to an isolated audit directory.
  Candidate patched through Cargo configuration; original files untouched.
  `cargo test --lib --tests --offline`: 179 unit + 17 integration + 9 property
  tests pass. `cargo check --all-targets --offline` passes, with a preexisting
  unused import warning in the consumer benchmark. Four lazy-construction
  and 22 filtered sorting benchmark smoke cases pass. An unfiltered benchmark
  smoke run was interrupted during its large workloads; no exhaustive
  downstream performance claim is made.
- **hashrope-bio**: source revision
  `206d1af783ca4c34626cfe1f9ce9a7bd7d4275cb` tested unchanged with an isolated
  candidate-dependency harness. Ten library tests pass and six binaries
  compile. Its original/published `^0.2.1` requirement excludes this crate
  line, and its local relative path is stale. Adoption remains blocked until
  that downstream manifest is separately updated. A downstream doctest
  double-borrows one arena (`E0499`); it also fails against published 0.3.1.
- **Pollard**: three copied store tests pass with published Pollard 1.6.0 and
  Hashrope 0.2.2 wheels. The 785-byte old-version JSONL fixture is checked
  directly by Rust using incremental, reconstructed, chunked and streaming
  paths. Its hash remains `2269240639182000324`.
- **wattwarden**: indirect Pollard Hashrope-extra dependency verified at
  `a1845ae5b0d4307610780849cf5dd1e5a3c4c193`; application/provider execution
  was not performed. No unverified third-party framework runtime is claimed.

The existing modified biology submodule and consumer working trees were
preserved. No npm/Python implementation or credentials were modified. The
`Rust` CI workflow runs library/package tests on Linux, macOS and Windows,
alongside published-wheel conformance and public downstream checks, for code
pushes and pull requests. Native 32-bit runtime tests are not claimed; wasm32
supplies concrete 32-bit numeric regression coverage. The final release PR
records exact checked commit IDs and gate outcomes before merge/publication.
