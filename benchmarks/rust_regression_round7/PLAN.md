# Round seven: focused diagnosis of 100K churn and 1MiB LCP

The round-six final comparison and its separate confirmation disagreed on
independent 100K exact edit/comparison (paired ratios 1.0438 and 0.9834) and
same-arena independent-equal 1MiB exact LCP (0.9801 and 1.0680). This study asks
whether either penalty recurs under longer paired measurements, both code-slot
assignments, affinity controls and an identical-code control. Neither unchanged
source nor variable controls can establish that a difference is noise.

`baseline.zip` freezes the current round-six source; `prior_round5.zip` freezes
the previous round-five source. `round6_final_binaries.zip` retains the exact four
executables from the previous final comparison. Source/archive/executable and
prior result identities are recorded in `baseline_manifest.json`. This source
freeze precedes any new implementation change. No production source is edited
by the controlled diagnostic, and no version bump or publication is involved.
The retained full-sequence replication and any final optimization comparison are
separate protocols owned by the parent task; their data are not pooled here.

## Five controlled cases

1. **100K combined independent exact churn.** Use the original 100,000-byte ASCII
   generator, seeds 1/2/3, 50 insert/delete edits, maximum insertion length 32,
   independent target arena, alternating exact-match/final-byte mismatch target
   mask and representative untimed edit warmup. Time each 50-edit sequence with
   exact equality after every edit; verify the 25-match mask and final bytes.
2. **100K edit only.** Use the same initial bytes, edits and representative edit
   warmup, omitting independent targets and comparisons. Time all 50 edits and
   validate final bytes afterward. Its setup/history differs from combined mode.
3. **100K comparison only.** Prebuild all 50 actual edited roots and independent
   targets, preserve the representative warm edits, and verify every state and
   final-byte mismatch outside timing. Time only the 50 cross-arena exact
   comparisons. This follows the cross-arena frontier path; phase timings cannot
   be subtracted or added to infer an isolated comparison cost.
4. **1MiB same-arena independent-equal LCP.** Recreate the original microbenchmark
   byte pattern and its complete local arena construction: first, independent,
   early-mismatch, late-mismatch and shared-late roots, including the split/join.
   Verify byte contents and independent identities outside timing. Preserve the
   target's untimed LCP/equality/hash-prefix checks, then three untimed LCP calls
   and seven samples of 100 LCP calls in one session. Each sample's answer must
   be 1,048,576, and the arena must not grow. The preceding full-sequence
   same-root LCP/equality/hash timed cells are deliberately absent, replaced by
   focused calibration/warmup. Additional byte validation also changes access
   history. This case uses the same-arena recursive comparator; conclusions from
   cross-arena churn cannot simply be transferred to it.
5. **Original-growing balanced concat control.** Construct the original two
   independent 64KiB roots, perform three untimed concatenations, then retain
   seven samples of 2,000 concatenations in the same arena. Check the initial
   513 nodes, final 14,513 nodes, length, weight, fingerprint and byte contents
   outside timing. Do not reset the arena between those seven samples or change
   this estimator into an allocation-growth-weighted mean.

Each complete session creates fresh bounded state and destroys it outside the
local operation timers. No source counters or private capacity hooks are added
to timed cores. Input handles/arenas and results are black-boxed. This diagnostic
has its own wrappers and additional oracles; it is not the unchanged original
full-sequence executable.

## Fixed matrix and measurement

Build three dual-slot executables from byte-identical frozen runtime sources:
prior/current, current/prior, and current/current as A/A. Only the diagnostic
crate identities/manifests differ. Swapped slots address some code-placement
sensitivity, while A/A measures identical-code variation rather than assuming
it is zero. Preserve all source, lock, compiler and executable identities and
archives before timing. Build everything before the exclusive measurement slot.

Run exactly six fresh processes for each build in each affinity arm: ordinary
unpinned and Windows-pinned to the lowest allowed logical processor. This is
**36 processes**. Rotate/reverse build order and alternate affinity order across
process seeds 1–6. Churn seeds are 1/2/3 repeated twice, with opposite ABBA/BAAB
orders for each repeated seed. Shuffle the five case positions using the process
seed; choose ABBA/BAAB independently of shuffled position. Every case receives
three of each slot order in each build/affinity arm. Construction, calibration
and paired warmup initial orders are balanced as well.

For each case calibrate one common integer count of complete sessions for both
slots and every block, targeting approximately 25 ms accumulated timed work in
the slower arm, with at least one session. Follow it with at least 20 ms paired
wall-clock warmup. The target is not a guaranteed minimum; session granularity
and calibration variability are retained. In particular a seven-by-100 1MiB LCP
session can itself exceed 25 ms. Measure one ABBA or BAAB quartet per case,
retaining every sample and session statistic. Do not time selected favorable
cases or stop after an interesting estimate.

For concat and LCP, the session statistic is the median of seven per-call sample
times; the block statistic is the median of session medians. Churn blocks use
median per-edit session times. Also retain accumulated timed means and complete
wall durations as separate estimators, not substitutes for the primary median.
Record node counts, root/leaf addresses, requested/actual affinity and processor
IDs at complete-block endpoints, including setup. Endpoints cannot rule out
intervening migration, and addresses do not establish a placement mechanism.

## Analysis and interpretation

The unit of analysis is a fresh process. For each case use the median of each
slot's two block statistics and form current/prior, or B/A for current/current.
Report all six ratios, geometric mean, range, slower-direction count and an
exploratory 95% log-t interval separately for every slot assignment and affinity.
Keep primary medians, secondary timed means and wall estimates distinct. Inner
sessions are not independent process replicates; repeated seeds are repeated
workloads, not new workload draws. Intervals are not multiplicity adjusted,
equivalence tests or universal bounds.

A stable penalty should reproduce direction and plausible magnitude across the
original-binary recurrence, both slot assignments and paired process samples,
with A/A variation shown explicitly. Protocol disagreement and unfavorable
secondary cells remain visible. Failure to reproduce a 4–7% point estimate does
not establish zero smaller effect; an interval containing one does not prove
noise. No production change is selected solely from source similarity or one
timing protocol. Any subsequent candidate must pass correctness/fingerprint
checks and a separately declared unchanged full-workload comparison against the
frozen round-six baseline.

The parent coordinates an exclusive CPU slot: no concurrent agent builds, tests
or timings. Machine frequency, thermal state and unrelated OS work remain
uncontrolled. No early stopping, sample removal, failed-process replacement,
source change during timing or artifact overwrite is permitted. Retain raw
stdout/stderr, arguments, status and all hashes. A timeout/failure preserves its
partial evidence and stops the fixed cohort. Reparse raw samples and verify
contexts, sidecars, oracles, counts and identities during analysis. No controlled
data are pooled with previous rounds or the separate retained-binary run.
