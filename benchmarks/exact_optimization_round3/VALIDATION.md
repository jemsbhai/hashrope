# Edit and metadata optimization validation

The immutable baseline is the completed second optimization round. The selected
changes remove three unnecessary leaf-vector clones and group the default
polynomial hash into eight-byte Horner steps. Fingerprint values and defaults
remain the compatibility contract; the hashing implementation intentionally
changes.

| Suite | Passed |
|---|---:|
| Python core | 203 |
| Rust core and selected integrations, release mode | 144 |
| Paper helpers and consumers | 255 |
| Biology Python | 41 |
| Biology Rust, including its compile-only doctest | 15 |
| **Total distinct selected tests** | **658** |

The six new compatibility tests also pass separately in debug mode and against
the frozen baseline library; those repeated executions are not added to the
total. The release `no_std + alloc` check passes. The core doctest command
succeeds with zero doctests discovered. Existing dependency deprecation and
unused-code warnings remain unrelated to this change.

The E3 experiment that materializes a 10 GB repeated buffer and E5 sliding
throughput experiment are explicitly filtered as in prior rounds. All other
core library/integration tests run. Symbolic enormous Repeat/DAG cases remain
covered by the exact-comparison suites.

## Compatibility evidence

The new tests freeze scalar arithmetic independently of the optimized hash
kernel, including accepted custom profiles whose legacy reduction is not
ordinary modular arithmetic. They check block/tail boundaries, long inputs,
maximum block carries and full-u128 reduction boundaries. Edit/split/repeat
traces preserve baseline fingerprints, node layouts and lazy materialization
states across eager/lazy and custom-profile cases. Empty and clamped operation
edges retain their behavior. The `(17, 2)` shape-sensitive hash example remains
unchanged, while exact equality remains independent of it.

The source-scope audit verifies six entire untouched modules and 34 Python
definitions. Of 97 existing Rust production functions, only
`ensure_hash_inner`, `split_inner`, `hash_range` and `PolynomialHash::hash`
change; `hash_default_profile` is the sole new production function. Existing
signatures, `NodeInner`, `Arena`, `PolynomialHash`, constant declarations and
type aliases remain unchanged. Differential tests establish semantic hash
compatibility; source checks alone do not.

The optimized archive also freezes Rust integration tests and their checksums.
Final production source hashes are:

- `rope.rs`: `0eae07ab1c8c4563f64e5babfe017bf6731558c7ea397248c327327cd75ed468`
- `polynomial_hash.rs`: `5cf3dbb650ae183523eec6293ac4951f2079129ba04c78e6106c57c555ecae79`

`unchanged_code.json`, `optimized_manifest.json`,
`final_source_verification.json` and `preserved_results.json` provide the
machine-readable source checks. Python is byte-identical across the measured
snapshots and serves as a control.

## Commands and retained logs

`validation/run_checks.py` runs the Python core, paper and biology suites plus
`cargo test --release` in the biology Rust package. It asserts the local
Hashrope import path, uses offline tokenizer settings and records exact
commands, return codes and before/after core hashes in
`validation/python_and_bio.json`. The source hashes remained unchanged during
all four checks. Full outputs are in the adjacent `.log` files.

The core Rust commands were:

```powershell
cargo test --release --offline --manifest-path packages/rust/Cargo.toml --lib --tests -- --skip e3_repeat_node_scaling --skip e5_sliding_window_throughput
cargo test --offline --manifest-path packages/rust/Cargo.toml --test test_edit_optimization
cargo check --offline --manifest-path packages/rust/Cargo.toml --lib --release
cargo test --offline --manifest-path packages/rust/Cargo.toml --doc
```

Their logs are respectively `validation/core_rust_release.log`,
`validation/edit_compatibility_debug.log`,
`validation/core_rust_no_std_release.log` and
`validation/core_rust_doctests.log`. The same six compatibility tests compiled
against the immutable baseline are recorded in
`validation/edit_compatibility_baseline.log`.


## Final matrix and documentation checks

The complete frozen-source matrix finished with all five context audits
passing: Rust, paper CPU, extended real-pair/K, biology and microbenchmarks.
`summary.json` records `complete_matrix: true`. All byte/token and per-edit
boolean oracles passed. Alternate-output aggregation reproduced the same
summary, and `--kind bio --resume` reused only verified matching results.
Every new Python benchmark/report script parses successfully.

Final production and Rust integration-test sources are byte-identical to the
measured optimized archive. The source-scope audit records exactly four changed
existing functions and one addition; 243 prior benchmark source/artifact files
remain byte-identical. `final_audit.json` records the aggregation/resume checks,
and all three repositories pass `git diff --check` with line-ending warnings
only.

The revised paper, claim ledger and findings distinguish the latest cohort,
unchanged Python controls, changed fingerprint-candidate arms and slower exact
microbenchmark cells. The desktop LaTeX compiler successfully compiled the
final `paper/main.tex`, including the last arena-layout/workload-scope precision
edits. No production or measured workload source changed after freeze.
