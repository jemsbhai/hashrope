# Final round-six comparison

Before freezing or measuring, compare the frozen round-five optimization with
the reviewed round-six candidate. This is a separate optimization experiment,
not a continuation or pooling of the published-0.3.3 diagnostic cohorts.

Use byte-identical round-five workloads and locked dependencies: all 42 micro
cells with seven samples and all 24 edit/comparison cells with three samples.
Keep their original iterations, warmups, setup, oracles and complete case order.
Build all four executables before measurements. Preserve source, driver, lock,
compiler and executable hashes and archive the built binaries.

Run exactly nine fresh process pairs per cohort (36 processes, 6,588 timing
samples). Alternate baseline/candidate then candidate/baseline, beginning with
baseline first. Churn seeds are 1,2,3 repeated three times, invocation zero.
Run all micro processes before churn. Record the slight five/four order
imbalance and per-order/per-seed summaries. Seeds repeat workloads, not nine
independent workload draws. No outer Python interleaving or process pinning.

Root holds an exclusive CPU slot during timing: no other agent benchmarks,
builds or tests. CPU frequency, temperature and unrelated OS activity remain
uncontrolled. No early stopping, exclusions, replacement processes, or source
changes during the cohort. Any failure preserves partial output and stops it.

For every cell report per-process medians, the geometric mean of the nine
paired candidate/baseline ratios, range, direction count, and exploratory 95%
log-t interval (eight degrees of freedom). Analyze at the process-pair level,
not as independent inner timing samples. Intervals are not multiplicity
adjusted or equivalence bounds. Keep all 66 cells, including Ropey and
historical hash-only controls, and layout rows. Primary outcomes are balanced
concat and independent 1M exact edit/comparison; preserve the independent LCP
and compressed-DAG gains as well as every remaining regression.

Reparse raw outputs during analysis and verify process order, arguments, seed,
binary identity, sample counts and oracle success. Freeze with `run.py freeze`,
build with `run.py build`, measure with `run.py run`, then audit with
`run.py summarize`. Phases refuse to overwrite retained evidence. No version
bump or publication is part of this experiment.
