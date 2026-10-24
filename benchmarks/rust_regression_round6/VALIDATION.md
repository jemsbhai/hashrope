# Correctness and compatibility validation

All seven checks passed on Windows x86-64 with Rust 1.94.0. The
[final record](validation_final/results.json) contains source/test/lock hashes,
commands, exit codes and complete log hashes. Sources remained unchanged during
validation and match every source file in the frozen final benchmark candidate.
The [source archive](validation_final/sources.zip) preserves the validation
driver, full core, tests, fixtures and test dependency lock.

| Check | Result |
|---|---|
| `cargo fmt --all --check` | Passed |
| `cargo clippy --all-features --all-targets -- -D warnings` | Passed |
| `cargo test --all-features --release` | 221 passed, none failed or ignored; includes E1–E7 and bigint doctest |
| `cargo test --doc` | Passed; no default-feature doctests defined |
| `cargo check --lib --no-default-features --target wasm32-unknown-unknown` | Passed |
| `cargo check --lib --all-features --target wasm32-unknown-unknown` | Passed |
| `cargo doc --no-deps --all-features`, `RUSTDOCFLAGS=-D warnings` | Passed |

Two new private tests compare the balanced shortcut with the retained recursive
join path across balanced/unbalanced ordinary and repeated trees, three hash
profiles and eager/lazy modes. They check complete node sequences before and
after lazy hash materialization. Overflow coverage verifies length then weight
rejection before any allocation, including error messages. Existing independent
golden edit traces protect published fingerprints and tree layout.

The full suite retains collision, profile-mutation, cross-arena, range, budget,
lazy, shifted-repeat, huge shared-DAG, near-`u64::MAX`, Python/Pollard golden
fingerprint and bigint regressions, including all seven round-five additions.
No hash arithmetic, defaults, leaf policy, metadata fields, public signatures,
exact comparison code or bigint implementation changed in this round.

The first validation attempt stopped at a strict Clippy warning in a new test:
an unnecessary `unwrap()` on a locally constructed `Some`. Its
[record](validation/results.json), source archive, driver and logs remain.
Replacing that test expression changed no runtime code. The complete final
validation above passed afterward; no failure was discarded.

The separate code-generation prototypes and targeted checks are documented in
[codegen/README.md](codegen/README.md). Their final runtime prefix matches the
measured candidate; their earlier test snapshot precedes the Clippy-only fix.
The controlled diagnostic's byte/fingerprint/oracle checks and the final
benchmark's original oracles supplement, rather than replace, the full suite.

This is local validation plus wasm compilation, not a fresh Linux/macOS CI run
or wasm performance benchmark. E1–E7 internal timing output is retained as part
of the existing tests and is not used as paired performance evidence. Versions
remain 0.3.3; the changes are unreleased.
