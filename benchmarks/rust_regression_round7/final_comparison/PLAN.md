# Round-seven final optimization comparison

Declare this protocol before freezing or measuring an optimization candidate.
Compare frozen round six against the reviewed round-seven candidate, retaining
fingerprints, correctness, defaults, layouts and package versions. This is a
separate optimization comparison from the round-five/round-six diagnostics.

Use the byte-identical original 42-cell microbenchmark (seven samples per cell)
and 24-cell edit/comparison workload (three samples), with unchanged iterations,
setup, warmups, oracles, complete case ordering and locked dependencies. Build
and archive all four binaries before timing. Preserve source/driver/lock/compiler
and executable identities. No Python outer interleaving or process pinning.

Run exactly twelve fresh process pairs per cohort, 48 processes and 8,784 timing
samples. Alternate baseline/candidate then candidate/baseline, six of each.
Churn seeds 1,2,3 repeat four times, invocation zero; each seed receives both
execution orders twice. All micro processes precede all churn processes.

The primary outcomes are independent 100,000-byte exact edit/comparison and
independent-equal 1,048,576-byte exact LCP. Preserve the previous concat gain,
compressed-DAG behavior, every other cell and Ropey/hash-only controls as well.
Hash-only controls are not substitutes for exact APIs. Record all memory-layout
rows. Root holds an exclusive CPU slot with no concurrent agent builds, tests,
disassembly or timings; frequency, temperature and external OS work remain
uncontrolled.

Analyze paired ratios of per-process medians: geometric mean, full range,
direction count, execution-order and seed strata, and exploratory 95% log-t
interval with eleven degrees of freedom. Also retain medians of process medians;
their ratio may differ from the paired estimate. Inner samples and repeated
workload seeds are not independent workload draws. Intervals are unadjusted,
not equivalence bounds or universal performance guarantees. Do not pool earlier
diagnostic cohorts with these measurements.

No source change, early stopping, sample removal, failed-process replacement or
selection of a favorable subset is permitted during this cohort. Retain a
failure's partial output and stop. Assess the entire fixed run and explicit
remaining costs before retaining an optimization; instruction-count changes
alone do not establish speedups. Any further experiment requires a separate
declared protocol and preserved results.

Run `run.py freeze`, `run.py build`, `run.py run`, then `run.py summarize`.
Phases refuse overwrites. Analysis reparses every raw record and validates
orders, arguments, seeds, identities and sample counts. No release is included.
