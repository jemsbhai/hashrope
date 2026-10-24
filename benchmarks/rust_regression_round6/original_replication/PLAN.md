# Retained round-five executable replication

Declared before any new measurements. Reexecute the exact four archived Windows
executables from round five: published 0.3.3 and the round-five LCP optimization,
each with the original microbenchmark and edit/comparison workload. No rebuild,
source change, dependency change, sample exclusion, or overwritten evidence is
allowed. This checks recurrence in the original executable layout separately
from the new controlled harness.

Run six fresh pairs per cohort, ordered baseline/candidate then
candidate/baseline, alternating throughout. Churn seeds are 1, 2, 3, 1, 2, 3,
always invocation 0. Each seed therefore has both execution orders. Run all
micro processes, then all churn processes, preserving each executable's full
sequence, setup, three-call micro warmup, iterations, samples and byte oracles.

Retain all 24 processes, 42 micro cells with seven samples, 24 churn cells with
three samples, and all layout rows: 4,392 timing samples in total. Primary
outcomes are balanced concat of two 64 KiB inputs and independent 1M-byte exact
edit/comparison. All other cells remain visible, including Ropey and unsafe
historical hash-only controls. A failed process stops the run after retaining
its partial stdout/stderr and status; it is not replaced or silently retried.

Use an exclusive CPU slot without concurrent agent builds/tests/timings. Leave
normal process affinity unpinned. CPU frequency, thermal state and OS activity
are uncontrolled. Outer historical Python interleaving remains absent. Paths
may differ after extraction; executable bytes must match their original hashes.

Analyze paired ratios of per-process medians, retaining their range, slower-pair
count, execution-order and seed breakdowns, and exploratory 95% log-t intervals
over the six process pairs. Intervals are not adjusted for multiple cells, do
not establish equivalence, and do not treat inner samples or repeated seeds as
independent workload draws. These new six pairs stand alone; the earlier nine
remain historical evidence, not data to pool adaptively into a desired result.

Freeze source/parser, plan, archived executables and original provenance hashes
before running. Verify every raw record, expected argument/seed/order, cell and
sample count, and executable hash again during analysis. Run `run.py freeze`,
then `run.py run`, then `run.py summarize`; every phase refuses overwrites.
