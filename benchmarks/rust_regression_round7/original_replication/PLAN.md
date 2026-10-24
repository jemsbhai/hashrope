# Round-seven retained-executable recurrence

Before new measurements, test the two remaining round-six signals with the exact
four archived final-comparison executables. Baseline is the round-five core;
candidate is the round-six concat optimization. There is no rebuild, runtime
source change, dependency change, discarded sample or replacement process.

Run exactly twelve fresh pairs per cohort, alternating baseline/candidate then
candidate/baseline. Churn seeds are 1,2,3 repeated four times, always invocation
zero. Each seed has both execution orders twice. Run all micro processes before
all churn processes. Preserve the complete original sequence, setup, warmups,
iterations, sample counts and byte/fingerprint/result oracles: 42 micro cells
with seven samples and 24 churn cells with three samples, 48 processes and
8,784 timing samples in total.

The two primary outcomes are independent 100,000-byte exact edit/comparison and
independent-equal 1,048,576-byte exact LCP in one arena. Balanced concat and all
other cells remain visible, including Ropey and historical unsafe hash-only
controls. The new twelve-pair cohort is analyzed separately from round six's
nine pairs and six-pair confirmation; there is no adaptive pooling.

Root holds an exclusive CPU slot. No concurrent builds, tests, disassembly or
agent benchmarks. Processes remain unpinned, matching the original protocol;
CPU frequency, temperature and unrelated OS work remain uncontrolled. A failed
process preserves partial stdout/stderr and status, then stops this fixed run.

Analyze paired ratios of per-process medians, their geometric mean, range,
direction counts, execution-order and seed strata. Report exploratory 95% log-t
intervals over twelve process pairs (eleven degrees of freedom), not over the
inner timing samples. These unadjusted intervals are not equivalence tests or
universal bounds; repeated seeds are repeated workloads. Report both the paired
ratio and medians of process medians when they disagree. No early stopping or
sample exclusion is allowed.

Freeze plan/driver/reference hashes first; verify archived executable bytes
before execution and during analysis. Reparse all raw records and verify every
argument, seed, order, cell and sample count. Use `run.py freeze`, `run.py run`,
then `run.py summarize`; phases refuse to overwrite their evidence.
