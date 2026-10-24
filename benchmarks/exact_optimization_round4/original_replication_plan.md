# Retained original-binary replication protocol

This follow-up was fixed after inspecting the controlled dual-core experiment
and before running any follow-up samples. It addresses the dual-link layout
limitation directly, rather than selecting more samples from the first study.

Run exactly 12 fresh process pairs using the two retained original round-three
microbenchmark executables. Their SHA-256 hashes must match the original
`round3/results/microbench.json` values, and their copied Rust harnesses must
match the unchanged `comparison_safety/bench_rust.rs` hash. No rebuilding,
fixture changes, extra warmup, or editing of the original sequence is allowed.

Odd pairs run baseline then optimized; even pairs run optimized then baseline,
giving six pairs in each order. Processes retain the ordinary unpinned affinity
mask. Run all pairs without early stopping, retain every output row and all
seven samples per cell, and check the original executable's exit/oracle status.
No Python microbenchmark is interleaved in this follow-up; this differs from
the historical overall driver while preserving the complete Rust executable.

Primary cells are 1 MiB independent exact equality and 1 MiB late-mismatch exact
LCP, the reported round-three regressions. Report the other comparison cells
too. For each cell compute the optimized/baseline ratio of within-process
medians for each matched process pair, then report all 12 ratios, their
geometric mean, spread, and exploratory 95% log-t interval (11 degrees of
freedom). Inner samples are not independent experimental units. Per-cell
intervals are not adjusted for multiple comparisons or universal bounds.

Agreement across this replication and both controlled slot assignments would
support a repeatable code/runtime penalty. A failure to reproduce the original
12–17% effect limits that inference, while smaller, layout-conditioned,
workload-history, or machine-specific effects can remain. This diagnostic does
not change the final benchmark protocol or justify a speculative core change.
