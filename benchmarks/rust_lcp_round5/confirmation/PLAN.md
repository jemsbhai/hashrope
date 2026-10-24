# Fixed supplemental confirmation: six additional pairs per cohort

This protocol is declared after the original fixed three-pair study and before
any supplemental timing. It investigates unresolved initial estimates:
independent 1M exact churn had an 11.5% slower paired geometric mean (two of
three pairs slower), independent 10K exact churn was 2.55% slower in all three
pairs, and 64KiB independent-equal LCP was faster overall but slower in one pair.
The 1MiB independent-equal and late LCP gains, and the DAG LCP gain, will also
remain prominent. All cells are retained so confirmation cannot select only
favorable operations.

No production source, workload, dependency, feature, executable or compiler
configuration changes. This phase uses the **same four archived executables**
from the original round-five build, extracted without rebuilding and verified
against their recorded SHA-256 hashes. The original source archives, lockfiles,
driver hashes, result file and summary are bound into this phase's manifest.

Run exactly **six additional fresh process pairs for each cohort** (24 new
processes):

| Supplemental pair | Global pair | Order | Churn seed/invocation |
|---:|---:|---|---|
| 0 | 3 | candidate, baseline | 1, 0 |
| 1 | 4 | baseline, candidate | 2, 0 |
| 2 | 5 | candidate, baseline | 3, 0 |
| 3 | 6 | baseline, candidate | 1, 0 |
| 4 | 7 | candidate, baseline | 2, 0 |
| 5 | 8 | baseline, candidate | 3, 0 |

The original three pairs remain intact, and the combined nine-pair sequence
therefore alternates baseline-first/candidate-first throughout. Run all six
micro pairs, then all six churn pairs, as in the original cohort ordering.
Each micro process retains all 42 cells and seven samples; each churn process
retains all 24 cells and three samples. Fixture order, warmup, iterations,
50-edit histories, unsafe historical hash-only arms, and outside-timing oracles
are unchanged. Expected supplemental totals are 3,528 micro timing samples and
864 churn timing samples; pooled totals are 5,292 and 1,296, respectively.

Root grants one exclusive CPU window. No other agent builds, tests or timings
will run concurrently. Affinity remains inherited and unpinned. Thermal state,
CPU frequency and OS background work remain uncontrolled. The outer Python
interleaving is still absent. Run every predefined process without early
stopping, excluding samples or replacing failed output. All stdout, stderr,
arguments, status and executable hashes are retained. Existing artifacts cannot
be overwritten; failure or timeout stops the cohort with evidence retained.

The report must show the supplemental six pairs **separately** from the original
three and also pool **all nine pairs**. It will include every cell, all process
medians, paired ratios, their geometric mean and range, and number of slower
pairs. The original predeclared all-three flag is preserved in the original
summary; this phase reports directions/counts without redefining that outcome.
Ranges are descriptive, not confidence intervals. Nine process pairs with
repeated seeds are not nine independent workload draws, and pooled results
must disclose the two measurement phases.

Consistency across phases can strengthen an estimate; disagreement is retained
as uncertainty. Neither unchanged-control variability nor favorable pooled
averages alone establishes a causal explanation or rules out all smaller
regressions. No additional optimization is selected during confirmation.
