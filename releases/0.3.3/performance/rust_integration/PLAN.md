# Rust 0.3.3 integration check

This is a targeted comparison of the release integration against the completed
round-four optimized snapshot. The reference already has exact comparisons. It
is not published Rust 0.3.2, which lacks these APIs. The integration adds the
published 0.3.2 arithmetic/overflow/parity work to the prior correctness and
optimization work. No performance-driven change is selected in this check.

The reference archive, candidate sources, workload files, and this driver are
hashed before building or timing. Both arms use fresh release binaries from the
same compiler and default features. The unchanged `bench_rust.rs` workload runs
three fresh process pairs, ordered reference/candidate, candidate/reference,
reference/candidate. Every one of its 42 timing cells and seven samples is
retained, including the hash-only historical comparison arms and memory rows.
The Rust sequence, fixture order, warmup, iteration counts, and byte oracles are
unchanged. The surrounding Python microbenchmark interleaving used in the full
round-four matrix is omitted. These results therefore do not resolve that
matrix's previously disclosed context-sensitive comparison results.

The unchanged `edit_compare.rs` workload also runs three fresh process pairs,
with `(seed, invocation)` equal to `(1,0)`, `(2,0)`, `(3,0)` and the same alternating
arm order. Its 50 edits, 50% equal/50% final-byte-different targets, 10K/100K/1M
sizes, independent/shared histories, four comparison arms, three samples per
process, and full outside-timing result checks remain unchanged. This is three
pairs, not the prior nine-pair cohort. Hash-only candidate arms are explicitly
not exact alternatives. Ropey and its resolved transitive dependencies remain
locked to the original workload versions. The Cargo lock may change only for
the path package's release version and its optional dependency resolution.

All binaries are built before measurements. No concurrent tests, builds, or
other agent timings run during the exclusive CPU window. Processes inherit the
ordinary unpinned affinity; hardware frequency/thermal state is not controlled.
No sample is removed, no early stopping is allowed, and existing result files
cannot be overwritten. Interpretation uses per-process medians, paired ratios,
their spread, and all cells. With only three process pairs, effects are
diagnostic, not tight statistical acceptance bounds or a universal efficiency
guarantee. Recurring material slowdowns will be reported rather than dismissed
as noise.
