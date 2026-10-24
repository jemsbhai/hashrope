# Round six: regression diagnosis before another optimization

This study investigates two mixed round-five estimates: independent 1M exact
churn was 5.1% slower pooled (five of nine pairs slower; initial/supplemental
ratios 1.115/1.020), and concat was 8.4% slower (six of nine pairs slower).
The source of those differences is not established. Unchanged source, overlap,
or variability in controls cannot by themselves prove noise or rule out effects
from generated code, allocation history, placement or machine state.

`baseline.zip` preserves the reviewed round-five candidate, the optimization
baseline for any final change. `released_033_reference.zip` preserves the
published 0.3.3 source used as the earlier reference. `round5_binaries.zip`
preserves the exact original four executables. Their hashes, source maps and
historical result identities are in `baseline_manifest.json`. Old evidence is
never modified. This protocol precedes new diagnostic builds and measurements.

## A. Retained original-executable recurrence

Run exactly six fresh pairs for each original cohort, using the four original
executables without rebuilding. Preserve all 42 micro cells with seven samples,
all 24 churn cells with three samples, all original oracles and the complete
fixed sequence. Process order alternates released/current then current/released,
three pairs each way. Churn seeds are 1/2/3 repeated twice, invocation 0. Run all
micro pairs then all churn pairs. This is 24 fresh processes and 4,392 timing
samples. It tests recurrence in the original binary/sequence context; it does not
replace the nine retained round-five pairs or pool them as a new controlled run.

## B. Longer, within-process controlled diagnostic

Build three dual-slot executables from immutable copies: released/current,
current/released, and current/current as an identical-code A/A control. Only
diagnostic crate identities/manifests change; copied runtime source bytes must
match their frozen cores. Both slot assignments are needed because dual-linking
changes code layout. A/A is a measured control, not an assumed zero-variance null.

Run six fresh processes for each executable in each of two affinity arms:
ordinary unpinned and Windows-pinned to the lowest allowed logical processor.
The fixed matrix contains **36 processes**. Executable and affinity ordering are
counterbalanced across process seeds 1–6. Record requested/actual affinity and
processor IDs at each complete block's endpoints (including its untimed setup);
endpoints cannot rule out intervening
migration or sibling-thread contention. Record fixture root/leaf addresses and
node counts outside timing as descriptive placement evidence, not causal proof.

Each process contains the following five cases in a seeded shuffled order:

1. **Original concat growth trajectory.** Construct the original two independent
   64KiB roots, perform three untimed concatenations, and retain all seven samples
   of 2,000 concatenations per fresh arena. Do not reset between those samples.
   Repeat complete bounded sessions to obtain a longer block; drop each arena
   outside timing. This preserves allocation-growth positions and the original
   median-of-seven estimator.
2. **Bounded prewarmed concat.** The same roots, but 50,003 untimed concatenations
   before the seven samples. An independent isolated allocator/capacity replay
   must verify the expected 50,513 initial nodes, capacity 65,536, and 64,513 final
   nodes with no timed growth in both cores before this case is launched. No
   counters, capacity hooks or source changes are present in timed cores. This
   control also changes retained memory and operation history; it is not a
   single-variable causal allocation experiment. If the assumption fails, stop
   and revise the protocol before any timing rather than silently adapting it.
3. **Combined independent 1M churn.** Preserve the original 50 seeded edits,
   independent targets, untimed representative edit warmup, exact equality,
   25-match mask and final-byte oracle. Each session creates fresh bounded state,
   with setup/drop outside timing, then times its 50 edit-and-compare operations.
4. **Edit-only 1M churn.** The same edits and representative warmup, with target
   construction and comparisons omitted. Each session times its 50 edits and
   verifies the final bytes afterward. This measures a distinct setup/history;
   subtracting it from combined time is not an isolated comparison cost.
5. **Comparison-only matched edited states.** Prebuild the original independent
   targets and all 50 exact edited source roots, including the same untimed warm
   edits. Time only the 50 exact comparisons, preserving the alternating answer
   mask and verifying all source/target content outside timing. Its access history
   differs from the combined mode, so phase times are not additive components.

Churn fixture seeds are 1/2/3 repeated twice across the six process seeds. Case
order and ABBA/BAAB order are controlled independently of that seed and of each
other; the same fixture/case receives three of each slot order. Construction and
calibration/warmup order are balanced as well.

For each matched case use one common calibrated number of complete sessions for
both slots and every measured block, targeting **25 ms accumulated timed work**
in the slower arm. Calibration and paired warmup of at least **20 ms** precede
measurement and are not samples. Execute one ABBA or BAAB quartet per case.
Black-box operands and results. Setup, oracles, telemetry, output and arena drop
remain outside timing, but their effects on cache/thermal history are a stated
limitation rather than assumed absent.

For concat retain **every seven-sample session trajectory**. Its session statistic
is the median of seven per-call times; the block statistic is the median of its
session medians. Report accumulated timed mean and wall duration separately:
replacing the original estimator with total 14,000-call time would overweight
growth samples. For churn use the median session per-operation time as the block
statistic and also retain aggregate timed means. Every measured session retains
its answer oracle and all durations; no selected fast sample substitutes for it.

## Analysis and decision criteria

The unit of analysis is a fresh process. For each case take the median of each
slot's two block statistics, then form the current/released ratio (B/A for A/A).
Report all six ratios, geometric mean, range and slower-direction count
**separately for slot assignment and affinity mode**. Exploratory 95% log-t
intervals describe process sampling uncertainty; they are not multiplicity
adjusted, equivalence tests, or universal bounds. Inner blocks/sessions are not
independent process replicates. A/A variability remains explicit.

A stable penalty should recur across retained original binaries, both dual-slot
assignments and controlled process samples, with magnitude considered against
A/A variation. Protocol disagreement is retained, not resolved by choosing a
favorable arm. The relevant estimates are approximately 5% churn and 8% concat;
failure to reproduce larger penalties does not establish zero smaller effect.
The source/disassembly and untimed allocation evidence may support a mechanism
but cannot alone establish timing causality. No production change is justified
solely by unchanged source or noisy controls.

If evidence warrants an optimization, freeze its separate candidate only after
correctness/fingerprint checks and review, then rerun the unchanged 42-cell micro
and 24-cell churn workloads against the round-five baseline. That final protocol
and fixed process count will be recorded before its timing. No further core
change is selected during the diagnostic phase.

## Execution and preservation

All builds complete before the exclusive timing slot. Root coordinates a single
timing owner; no concurrent tests, compilation or other agent timing is allowed.
CPU frequency, thermal state and OS background work remain uncontrolled. No
early stopping, sample removal, failed-process replacement or artifact overwrite
is allowed. Record raw stdout/stderr, arguments, status, compiler identity,
source/lock/executable hashes and archived binaries. Timeouts retain captured
output and stop the fixed cohort. Original full-sequence and controlled results
are stored and interpreted separately, with no outer Python-micro interleaving.
