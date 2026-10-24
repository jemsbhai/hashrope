# Second exact-path optimization validation

The baseline is the previous optimized working tree, captured before this
round's edits. It is neither repository HEAD nor the original hash-only core.
Historical reports and raw results remain intact.

## Correctness and build checks

| Suite | Passed |
|---|---:|
| Python core | 203 |
| Rust core and selected integrations, release mode | 138 |
| Paper helpers and consumers | 255 |
| Biology Python | 41 |
| Biology Rust, including its compile-only doctest | 15 |
| **Total** | **652** |

Rust's 138 tests comprise 103 library tests, 11 original exact-comparison
regressions, five first-round optimization regressions, ten new second-round
regressions, and nine existing integration checks. `cargo check --lib
--release` passes with `no_std + alloc`; the core doctest command succeeds
with zero doctests discovered. Existing unused-code/import and Python
dependency deprecation warnings remain unrelated to this change.

The existing E3 experiment that materializes a 10 GB repeated buffer and E5
sliding-throughput experiment were excluded. The changed paths instead have
symbolic huge-rope tests, including ranges near `u64::MAX`.

Final tested Rust source SHA-256:
`eb9401967c78b0e8eb55943fe457ef6adcdbc2a0a27755c4389fc786fcd1212a`.
The measured optimized archive must contain that same source. Python source
is byte-identical to this round's baseline and was rerun as an unchanged
control with the local package and offline tokenizer settings.

## Additional coverage and review

- Cross-arena long scans, uneven leaf/tree boundaries, clipped ranges and
  same-arena shared ranges agree with independent byte oracles.
- Every first-mismatch position in leaves of lengths 1–512, varied offsets
  and differing bit positions exercise safe word localization.
- Private deterministic tests exhaust visit budgets of zero, one and two,
  fill the 32-frame frontier before and after a verified prefix, and encounter
  a repeat after an 8 KiB common prefix. Fallback resumes at the first
  unverified byte.
- Huge self-concatenation DAGs cover shifted phases, cross-profile queries,
  final-byte edits, collision pairs and a difference at the last byte of the
  Fine–Wilf witness. Non-adjacent shared DAG and overflow checks still pass.
- Exact comparisons preserve lazy sentinels and node counts.

Independent review checked the frontier's byte-order and clipping invariant,
portable little-endian word localization, and structural period proof. Only
explicit repeats and identical-child concatenations establish a smaller
period. The witness compares the original range offsets before skipping any
bytes; saturated arithmetic cannot enable an insufficient witness.

The exact API remains independent of fingerprints. In particular, the
accepted legacy profile `(17, 2)` yields fingerprints 16 and 0 for equal
bytes `[0, 13]` built as one leaf versus two concatenated leaves. A new
regression verifies that exact equality/LCP still return the correct result.
This round does not change that legacy fingerprint behavior or introduce a
hash-inequality rejection shortcut.

## Reproducibility and compatibility

The snapshot and unchanged-code checks preserve Python source, Rust hash and
sliding modules, construction/edit/hash function bodies, node layout, default
parameters, and public API signatures. The benchmark matrix retains the
previous workloads, seeds, inner repetitions, controls and oracle checks.
Fresh process counts and biology's within-process adaptive repetitions are
reported separately; diagnostic timings do not replace final cohort results.

`verify_final_core.py` checks the final checkout against the measured archive;
`verify_preserved.py` verifies historical benchmark artifacts. Their JSON
outputs accompany the benchmark report.

The complete matrix finished with all oracle and provenance checks passing:
18 Rust workers, 18 paper workers, 18 extended workers, two biology workers,
and three alternating microbenchmark rounds for each implementation. This
includes 7,200 real-pair checks and 90 candidate-grid checks across the cores.
Final verification found every measured core file byte-identical to the
checkout and all 144 prior benchmark source/artifact files unchanged.

The updated paper compiled successfully with the desktop LaTeX compiler.
Diff checks pass in the library, paper and biology repositories, with only
line-ending normalization warnings.
