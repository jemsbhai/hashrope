# Exact-path optimization validation

Validated on 2026-10-23 before the final timing window. The baseline is the
previous exact implementation captured from the working tree, not repository
HEAD and not the original hash-only implementation. Sources are frozen in
`baseline.zip` and `optimized.zip`; manifests contain their SHA-256 checksums.

## Correctness and build checks

| Suite | Passed |
|---|---:|
| Python core | 203 |
| Rust library and selected integrations | 125 |
| Paper helpers and consumers | 255 |
| Biology Python | 41 |
| Biology Rust, including its compile-only doctest | 15 |
| **Total** | **639** |

Rust's 125 tests include 100 library tests, 11 original exact-comparison
regressions, five additional optimization regressions, and nine existing
integration checks. `cargo check --lib --release` passes with `no_std + alloc`.
`cargo test --doc` for the core passes with zero doctests discovered. The
existing unused `rebalance` warning remains.

Two existing Rust experimental binaries were excluded: E3 materializes a
10 GB repeated buffer, and E5 is a large sliding-throughput experiment. Neither
tests the changed exact-comparison paths. The compressed regression tests
instead compare huge logical ropes without materializing them, including
near-`u64::MAX` lengths.

The paper tests ran with the local optimized Python package and offline
tokenizer settings. Its three dependency deprecation warnings were unrelated
to this change. Biology tests used the same local core packages as the library
tests.

## Regression coverage

- Legacy byte and uint32LE-token collisions still produce the same fingerprints
  and are rejected by exact comparison.
- Byte-oracle comparisons cover equal and unequal independent trees, shared
  edits, different leaf boundaries, shifted ranges, both range starts,
  eager/lazy modes and differing hash profiles.
- Empty, invalid and overflowing ranges are checked before identity shortcuts.
- Huge repeat witnesses, independently allocated self-concatenation DAGs and
  non-adjacent shared DAGs remain compressed.
- Private deterministic Rust tests exhaust the fast traversal's visit budget,
  verify that normalization spends that budget, and check fallback after a
  nonzero matched prefix or at a repeat boundary.
- Python tests cover every mismatch position in short leaves and positions
  around the 64-byte scout and 64 KiB chunk boundaries.
- Comparisons preserve node counts and lazy hash sentinels.

## Compatibility checks

The snapshot comparison verifies unchanged checksums for six hash, sliding and
export modules, identical Python ASTs for 33 definitions outside `_leaf_lcp`,
and byte-identical bodies for 31 Rust constructor, hash and edit functions.
Only the exact comparison paths changed. Node layout, fingerprint arithmetic,
default prime/base, and public API signatures remain compatible.

After timing, the Python `rope_lcp` docstring was qualified to mention cache and
integer-arithmetic costs. `verify_final_core.py` compares the final checkout
with the measured optimized archive: all other core files are byte-identical,
and Python executable ASTs are identical after removing docstrings. The final
digests and this check are recorded in `final_source_verification.json`; the
measured source archive remains intact.

Every benchmark retains its result checks against byte/token oracles. These
fixture checks supplement the adversarial tests; they are not a proof that
hash-only decisions are safe for arbitrary input.
