# Rust LCP round five: premeasurement protocol

## Question and immutable reference

Does the proposed Rust LCP implementation improve the observed 0.3.3 exact-LCP
costs while preserving exact answers, fingerprint values and the surrounding
edit/comparison workload? The baseline is the published 0.3.3 implementation at
commit `ffba5ee87a4c9401dd8a2e529a6c045fd1c73eff`, captured from Git objects before
implementation edits. `baseline_manifest.json` records all source hashes. The
candidate will be frozen separately only after correctness review and tests.
Neither arm changes its package version or feature configuration for this study.
The committed manifest uses LF and the Windows working manifest uses CRLF;
their normalized bytes must be identical. Both raw snapshots are retained, and
newline-only differences are explicitly recorded rather than editing either arm.

This study measures both arms again. The release-integration and round-four
measurements remain historical context; their timings will not be divided into
new timings as if they were concurrent controls. In particular, round four found
protocol-sensitive equality results that its separate controlled diagnostics did
not resolve. This experiment does not establish a causal explanation for those
earlier results.

## Fixed workloads and process order

The Rust workload source files are copied byte for byte from
`releases/0.3.3/performance/rust_integration/workloads/`. Their source hashes,
original wrappers and lockfile are retained. Actual build lockfiles are the
byte-identical candidate lockfiles retained by that integration study: core
0.3.3, Ropey 1.6.1, smallvec 1.16.4 and str_indices 0.4.4. Builds use
`--locked --offline`; no dependency update is allowed.

1. The unchanged microbenchmark executes all **42 cells**, each with **7 timing
   samples**, plus its three layout/node-count rows. It runs **three fresh
   process pairs**, ordered baseline/candidate, candidate/baseline,
   baseline/candidate. Both arms enable the harness's `exact` feature and use the
   core's default features. Fixture order, three-operation warmups, iteration
   counts, byte/length oracles and the full sequence are unchanged.
2. The unchanged churn benchmark executes all **24 cells**, each with **3 timing
   samples**, in three fresh process pairs with the same alternating order.
   `(seed, invocation)` is `(1,0)`, `(2,0)`, `(3,0)`. Each cell retains its 50
   edits, half equal and half final-byte-different targets, 10K/100K/1M inputs,
   independent/shared histories, four arms and checks outside the timer.

The complete fixed cohort contains **12 fresh processes and 2,196 timing
samples**. No sample is removed, no early stopping is allowed, and interrupted
runs retain their raw output rather than being silently replaced. Hash-only
comparison arms are preserved as historical candidate filters, not exact
alternatives. The Ropey rehash control also calls Hashrope's hash implementation.

## Controls and provenance

All four binaries are built before the timing window using the same Rust
toolchain. The baseline, candidate, workload/analysis driver sources, compiler
identity, exact lockfiles, executable hashes and archived executables are
recorded before measurements. Every subprocess's arguments, return status,
stdout, stderr, sample counts, fixture identities and byte-oracle success are
retained. Source, binary and lock hashes are checked again after measurement.
Existing manifests, archives, build records and raw output cannot be overwritten.

Root coordinates one exclusive CPU timing slot. No competing tests, builds or
other agent timings run during that slot. Processes inherit ordinary unpinned
affinity; CPU frequency, thermal state and background OS activity are not
controlled. The outer Python-microbenchmark interleaving used in the full
round-four matrix is absent, as it was in the release-integration check.

## Predeclared interpretation

The primary LCP outcomes are independent-equal 64KiB and 1MiB inputs and the
independent 2TiB self-concatenation DAG. The long late-mismatch LCP cells are also
reported prominently. All equality cells, shared/early/identity/repeat paths,
core operations, churn arms and layout rows remain visible as controls and
potential regressions. Tiny calls include loop/dispatch overhead; layout sizes
are not measurements of peak or retained process memory.

For every cell the analysis reports all three process medians, their range,
paired candidate/baseline ratios, geometric mean of those paired ratios, ratio
range and number of slower pairs. The median of process medians is a descriptive
time and is not substituted for the paired-ratio estimator. A slowdown exceeding
10% in geometric mean with all three pairs slower is explicitly flagged; smaller
or inconsistent effects remain visible. Unchanged controls help describe run
variability but cannot prove that a change is noise.

Three process pairs provide a bounded diagnostic, not tight statistical
acceptance limits or a universal efficiency guarantee. Improvement requires
consistent measured evidence and a plausible code-path explanation; these data
alone do not prove absence of smaller regressions. If stronger confirmation is
needed, its protocol and fixed count will be declared in a separate artifact
before any extra timing. Correctness and fingerprint compatibility are validated
independently and are required regardless of performance.
