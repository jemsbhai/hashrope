# Round-four pre-execution protocol

Date: 2026-10-24. The user authorized investigation of the roughly 12–17%
comparison-only increases reported in round three, optimization only when
supported by evidence, unchanged correctness/fingerprints, and the same full
benchmark matrix and documentation update.

Before source edits, the completed round-three working state, tests and every
earlier benchmark study were preserved in `baseline.zip`: 497 files,
12,405,734 bytes, SHA-256
`3047745ce88cb5b8a69c37009711d7a7662443ed4afd0eafec12f2f841af2a5b`.
This preserved working baseline is distinct from the comparator selected for
measurement below. Historical files must never be overwritten.

## Hypotheses, diagnostics and decision

The slower round-three same-arena equality/late-LCP cells could reflect code
generation/layout, cache/state effects, hybrid-CPU scheduling, runtime noise or
an implementation issue. Unchanged exact function bodies alone do not identify
the cause. Counterbalanced, longer warmed diagnostic blocks compare immutable
pre-round-three and current sources. Source/binary hashes, arm/fixture/API order,
and any affinity choice must be retained. Identical-source or identical-binary
A/A and swapped-slot controls are separately labeled diagnostics. Instrumented
counters and assembly inspection, if used, remain separate from uninstrumented
final timings. Diagnostics may investigate additional controls, but do not
replace the original final workloads or erase the earlier regressions.

A candidate must preserve exact answers, fingerprints, accepted profiles,
public API semantics and persistence. Acceptance requires passing the existing
correctness/compatibility gates and repeatable focused evidence supporting the
change. No universal speedup or significance claim is pre-registered. Preserve
negative cells and rejected variants. Do not manufacture an optimization if the
observed regression does not justify one.

The final comparator choice and rationale are frozen in `measurement_plan.json`
before final measurement, using one of two explicitly authorized modes:

- **Candidate selected:** stage the completed-round-three working archive as
  `reference.zip`; compare it against the reviewed round-four candidate frozen
  as `optimized.zip`.
- **No production change justified:** stage the exact pre-round-three archive
  from `../exact_optimization_round3/baseline.zip` as `reference.zip`; compare it
  against unchanged current round-three source frozen as `optimized.zip`. This
  repeats the actual old/new question, rather than replacing it with A/A.

The established result keys `baseline` and `optimized` identify these measured
reference/selected cohorts only. `baseline.zip` always remains the untouched
round-four working-state preservation. Context fingerprints bind the plan,
measured archive identities, preservation archive, workloads and runtime.
There is no implicit comparator selection and no plan change after timing starts.

## Final matrix and reporting

After reviewed core/test freeze, required validation and a team CPU pause,
reuse the exact prior consumer runners, corpus inputs, builders, sizes,
match/mismatch mixes, seeds, warmups, repetition counts and correctness oracles:

- Rust edit plus target comparison: seeds 1/2/3 by three fresh processes per
  core, three samples of 50 ASCII edits, half equal and half last-byte mismatch;
  independent and shared targets at 10K/100K/1M bytes, Ropey exact and both
  fingerprint-candidate arms retained.
- Paper CPU: seeds 42/43/44 by three fresh processes per core, full EXP005
  size/fraction grid, controlled EXP017 through 2M tokens and CPU R2 through 1M
  tokens. Query construction remains separately recorded and excluded.
- Extended: nine fresh processes per core, 200 ShareGPT and 200 LMSYS pairs per
  dataset/seed and all five original bounded K cells. Preserve fixed inner
  radix/exact/NumPy timing order and untimed preparation.
- Biology: all 2,000 HIV records plus 16 equal controls, nine within-process
  paired repetitions per core with the existing adaptive batching and recorded
  inner counts. These are not nine fresh-process estimates.
- Structural microbenchmarks: the exact same Python/Rust fixtures, three
  alternating rounds of seven samples, including all operation, memory,
  same/independent/shared roots, early/late mismatch and 2 TiB Repeat/DAG cells.

Rust micro independent roots occupy one Arena and use the recursive path on
1 MiB repeating 0-through-255 bytes. Independent-target churn uses separate
arenas and 1M edited ASCII bytes, with long cross-arena frontier comparisons.
Those timings cannot be substituted or added. Full arena construction and
destruction are timed in the construction microbench; churn setup is excluded.
Both fingerprint-candidate arms use hashing from their measured core. The
pre-round-three/current pairing includes the historical scalar/grouped hash
change, whereas Ropey exact and the Python implementation remain unchanged.

Report means and sample SD across fresh-process medians, ratios of means and
paired win counts. Biology retains within-process paired batch medians;
microbenchmarks retain median of three round medians and the round range.
Dispersion bands/ranges are not confidence intervals. Re-evaluate all prior
claim limits without attributing unchanged Python movements to Rust.

## Integrity, environment and outputs

Final measurements use the same Windows 11 / Intel Core i9-14900HX workstation,
without affinity/frequency locking, to preserve the final protocol. Controlled
or pinned diagnostic runs remain labeled separately. No GPU serving, TTFT,
NVML energy or chromosome-scale experiment is rerun. No heavy team work may
overlap final timing.

Retain all new artifacts under `benchmarks/exact_optimization_round4/`.
Resume/aggregation reject mismatched plan/core/harness/data/workload identities.
Archive tested Rust sources and integration tests; verify final working files
against that archive and previous benchmark artifacts against the preservation
manifest. The production-scope audit includes both function bodies and text
outside them, so inline/code-generation attributes cannot be missed. Explicitly
review every permitted change before freeze. Update the paper and claim ledger
only after the complete final matrix, preserve earlier cohorts, and require a
successful final paper compilation.
