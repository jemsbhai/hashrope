# Round-three pre-execution protocol

Date: 2026-10-23. The user authorized investigation and optimization of edit and
metadata costs behind the remaining Rust-versus-Ropey gap, especially on small
inputs, followed by the same full benchmark matrix. Before any source edits,
the current working implementation, tests and all prior benchmark artifacts
were captured in `baseline.zip`: 363 files, SHA-256
`4c394e6adb99d6cd5f755e6fca139148a04777c2c468e75b715539710428f80d`.
This is the completed second optimization round, not repository HEAD.

## Investigation and acceptance

Use isolated diagnostic copies to separate edit traversal, allocation, leaf
rehashing and internal fingerprint/metadata work from exact comparison costs.
Retain counter evidence and rejected variants. Diagnostic timings explain
candidate changes; they do not replace final paired measurements.

Production changes must preserve exact answers, existing fingerprint values,
defaults, accepted custom-profile behavior and public API semantics. Compare
edit/split/repeat/hash traces against the immutable baseline, including eager,
lazy, cross-profile and adversarial inputs. Record intentionally changed Rust
functions rather than claiming their source bytes are unchanged. Unmodified
modules, Python controls and prior result files remain checksum-verified.
Run the correctness gates before source freeze. Accept no optimization that
fails those gates. Report gains, remaining gaps and regressions from every
final workload without dropping negative cells or requiring a universal win.

## Final matrix

After reviewed source/test freeze and a team CPU pause, capture an immutable
optimized archive. Alternate baseline/optimized process order on the same
workstation and reuse the existing consumer runners, data, builders, sizes,
match/mismatch proportions, seeds, warmups, repetitions and oracle checks:

- Rust edit plus target comparison: three seeds by three fresh processes per
  core, three samples of 50 ASCII edits, 50% equal and 50% last-byte mismatch;
  independent and shared targets at 10K, 100K and 1M bytes. Retain Ropey exact,
  fingerprint-candidate and Ropey rehash controls.
- Paper CPU: nine fresh processes per core; full EXP-005 size/fraction grid,
  controlled EXP-017 through 2M tokens, and R2 CPU through 1M tokens with eager
  query construction separately measured and excluded from query latency.
- Extended comparisons: nine fresh processes per core; original 200-pair
  ShareGPT/LMSYS cells and all five original bounded candidate-set cells.
  Retain the fixed radix/exact/NumPy inner order and untimed preparation.
- Biology: all 2,000 HIV records plus 16 independent equal controls; nine
  within-process paired repetitions per core with the existing adaptive inner
  repetition method and recorded counts. These are not nine fresh processes.
- Python/Rust comparison-safety microbenchmarks: the same three alternating
  rounds of seven samples, including construction/hash/edit, memory, early and
  late mismatches, shared paths, 2 TiB Repeat and self-concatenation DAGs.

The unchanged Python implementation is a control; its timing movements cannot
be attributed to Rust optimization. Fresh-process tables report means and
sample SD across process medians, ratios of means and paired improvement
counts. Microbenchmarks report median of three round medians and their range.
No new statistical-significance or universal-speedup claim is pre-registered.

## Outputs and limits

All new snapshots, diagnostic and final raw results, compiled-binary/locked
dependency provenance, context fingerprints and the report go in
`benchmarks/exact_optimization_round3/`. Resume and aggregation must reject
core, harness, data or workload identity changes. Preserve every previous
benchmark source/result; use fresh reproduction output directories.

The same Windows 11 / Intel Core i9-14900HX workstation is used without CPU
affinity or frequency locking. No GPU serving, TTFT, NVML energy or
chromosome-scale experiment is rerun. Paper claim ledgers and the working
draft will be updated only after the final measurements, with a successful
compile required for the revised paper.
