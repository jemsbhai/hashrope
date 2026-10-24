# Predefined controlled comparison diagnostic

This experiment investigates the round-three comparison-only microbenchmark
changes. It does not replace or modify `comparison_safety/bench_rust.rs` or the
final benchmark matrix. The immutable round-two implementation is compared with
the immutable completed round-three implementation. Production is untouched.

The original benchmark uses three operation warmups, fixed case/API ordering,
and 100 iterations per long-comparison sample. It interleaves hash-only work
and construction whose runtime changed greatly in round three. Comparator
source equality alone does not exclude compiler layout, heap placement, cache,
CPU frequency, scheduling, or runtime differences. Variation or overlapping
ranges alone will not be called proof of noise.

## Design fixed before timing

- Build three executables with fixed crate slots: round2/round3, round3/round2,
  and round3/round3 as an identical-core A/A control. Crate copies only change
  package identity; copied Rust sources are checked and recorded.
- Retain the original 64 KiB/1 MiB same-root, independent-equal, early-mismatch,
  late-mismatch, and shared-late inputs, plus the 2 TiB Repeat and self-concat
  DAG cases. Exact byte oracles check ordinary fixtures; known constructed
  periodic contents check compressed fixtures. A collision check runs untimed.
- Six fresh processes per executable and affinity arm, with counterbalanced
  executable order, fixture construction order, shuffled fixture order, and
  alternating LCP/equality order. Each matched fixture/API uses the same
  calibrated iteration count in every A/B block. Calibration and warmup do not
  become measured samples.
- Sustained paired warmup of at least 10 ms, followed by ABBA or BAAB blocks
  targeting 10 ms each. Inputs and outputs are black-boxed. The fixed wrapper
  and loop overhead remain part of tiny-case measurements.
- Preserve the default unpinned process arm. Separately run a controlled arm
  pinned to the lowest allowed logical processor using the Windows process
  affinity API; record requested/actual affinity and processor IDs at both
  ends of each timed block. This does not measure every intervening migration
  or CPU frequency and does not exclude sibling-thread contention.
- Record root-node and first-leaf addresses outside timing to expose placement
  differences, without treating address correlation as a causal explanation.
- Retain every raw block, source/binary hash, compiler identity, and copied
  source archive. The original full-sequence benchmark remains an independent
  control in the final rerun.

## Interpretation fixed before timing

The analysis unit is a fresh process: for each fixture/API take each slot's
median of its two blocks, then the paired round3/round2 ratio. Report all six
ratios, their spread, and geometric mean separately for slot assignment and
affinity mode. Exploratory per-cell 95% log-t intervals describe sampling
uncertainty; they are not adjusted for the many cells and are not universal
performance bounds. Identical-core A/A results quantify role/placement noise
in this diagnostic rather than a presumed zero-variance baseline.

A stable 12–17% penalty should reproduce in direction and approximate magnitude
across both slot assignments and controlled process pairs, beyond corresponding
A/A variability. Failure to reproduce that size limits the original inference;
it does not establish an absolute absence of regression. Smaller effects,
layout-dependent penalties, the original binary's layout, and machine-specific
behavior remain possible. Dual-linking itself changes code layout, so its
results must be read with the unchanged original full-sequence rerun.

No production optimization will be justified solely by unchanged source,
overlapping timing ranges, or a faster selected inner sample.
