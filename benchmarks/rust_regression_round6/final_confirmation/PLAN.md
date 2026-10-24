# Fixed confirmation of the final round-six comparison

Declared after reviewing the initial nine-pair final comparison and before any
confirmation measurement. The initial cohort estimated independent 100K exact
churn at 1.0438 (eight of nine pairs slower; exploratory interval 1.0015–1.0879),
balanced concat at 0.8550, and independent 1M exact churn at 1.0138. This bounded
follow-up checks those three primary outcomes and retains every other cell.
Its selection follows the initial results; it is not an independently selected
universal performance guarantee.

Use exactly the four retained executables in `../final_comparison/binaries.zip`,
with their recorded hashes and unchanged source, dependencies, workloads,
warmups, iterations, setup, byte oracles and full sequence. Baseline is the
frozen round-five optimization; candidate is the reviewed round-six
optimization. No rebuild, source change, sample exclusion, or overwritten
evidence is allowed. This phase does not use the older published-0.3.3 binaries.

Run exactly six new fresh process pairs per cohort, alternating baseline/candidate
then candidate/baseline, three pairs in each order. Churn seeds are 1, 2, 3, 1,
2, 3, always invocation zero; each seed therefore receives both execution
orders. Run all micro processes first, then all churn processes. Keep the normal
unpinned affinity and use an exclusive agent CPU slot with no concurrent builds,
tests or timing. Frequency, thermal state and unrelated OS activity remain
uncontrolled. Historical outer Python interleaving is absent.

Retain all **24 processes, 42 micro cells with seven samples, 24 churn cells
with three samples, and all layout rows: 4,392 timing samples**. Primary cells
are independent 100K exact edit/comparison, balanced concat of two 64KiB inputs,
and independent 1M exact edit/comparison. All remaining cells remain visible,
including Ropey and historical hash-only controls. No stopping early, selective
retiming, failed-process replacement or filtering is permitted. A failed process
stops the fixed cohort after retaining its partial stdout/stderr and status.

Report all 66 cells using paired ratios of per-process medians, their geometric
mean, range, slower-pair count, execution-order and seed breakdowns, and an
exploratory 95% log-t interval over the six process pairs (five degrees of
freedom). The intervals are not adjusted for multiple cells, are not equivalence
bounds, and do not treat inner samples or repeated seeds as independent workload
draws. **Keep this six-pair confirmation separate from the initial nine pairs;
do not pool or replace the initial result.** Report agreement and disagreement
between phases directly.

Before measurement, freeze this plan and driver plus the initial cohort's source,
workload, lock, executable, result and summary identities. Verify every saved raw
record, status, argument, seed, execution order, cell/sample count and executable
hash during analysis. `run.py freeze`, `run.py run`, then `run.py summarize` each
refuse overwrites. There is no build command, version bump, commit or publication
in this confirmation phase.
