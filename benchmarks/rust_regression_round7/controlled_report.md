# Round-seven controlled diagnosis

The controlled data retain a modest pinned 1MiB LCP penalty: both slot
assignments estimate about 2%, and the normal-slot primary interval excludes
one. The 100K combined estimates vary by slot assignment and do not establish
a consistent 4–7% penalty, but a secondary pinned comparison-only timed mean
shows a positive signal. These observations cannot be dismissed as noise or
replaced by a blanket no-regression conclusion. The concat control consistently
improves, showing that the protocol can detect the previously introduced concat
optimization in these inputs.

**Current means the frozen round-six source; prior means frozen round five.**
This diagnostic excludes any proposed round-seven implementation change. Both
source archives and the previous exact executables were preserved before new
edits. [baseline_manifest.json](baseline_manifest.json),
[controlled_freeze.json](controlled_freeze.json) and
[controlled_build.json](controlled_build.json) bind the source, protocol, lock,
compiler and binary identities; complete build inputs and executables are in
[controlled_sources.zip](controlled_sources.zip) and
[controlled_binaries.zip](controlled_binaries.zip).

## All 30 primary estimates

Each entry is the geometric mean of six fresh process ratios, an exploratory
95% log-t interval, and the number of processes with a slower numerator.
Prior/current and current/prior headings describe slot A/B placement, while
both report **current / prior**. A/A reports **slot B / slot A**, both containing
current code. Ratios above one are slower. Intervals are not adjusted for
multiple comparisons, are not equivalence bounds, and do not treat inner
sessions or repeated workload seeds as independent process replicates.

### Unpinned

| Case | Prior/current slots | Current/prior slots | Current/current A/A |
|---|---:|---:|---:|
| 100K edit + exact compare | 0.9951 [0.9666, 1.0244]; 3/6 | 1.0011 [0.9674, 1.0359]; 2/6 | 0.9988 [0.9729, 1.0253]; 3/6 |
| 100K edit only | 0.9503 [0.8856, 1.0197]; 0/6 | 0.9950 [0.9805, 1.0096]; 2/6 | 1.0113 [0.9795, 1.0440]; 4/6 |
| 100K cross-arena compare only | 0.9619 [0.8859, 1.0443]; 2/6 | 1.0267 [0.9645, 1.0930]; 4/6 | 1.0031 [0.9666, 1.0409]; 4/6 |
| 1MiB same-arena independent LCP | 0.9778 [0.9453, 1.0115]; 1/6 | 1.0266 [0.9346, 1.1278]; 3/6 | 0.9777 [0.9294, 1.0285]; 2/6 |
| Concat, original growth control | 0.8642 [0.8520, 0.8766]; 0/6 | 0.8508 [0.8008, 0.9040]; 0/6 | 0.9724 [0.9333, 1.0131]; 1/6 |

### Pinned

| Case | Prior/current slots | Current/prior slots | Current/current A/A |
|---|---:|---:|---:|
| 100K edit + exact compare | 0.9787 [0.9349, 1.0246]; 2/6 | 1.0145 [0.9799, 1.0504]; 5/6 | 1.0144 [0.9897, 1.0397]; 4/6 |
| 100K edit only | 0.9834 [0.9676, 0.9994]; 1/6 | 1.0239 [0.9613, 1.0907]; 4/6 | 0.9940 [0.9543, 1.0354]; 4/6 |
| 100K cross-arena compare only | 0.9917 [0.9613, 1.0230]; 1/6 | 1.0198 [0.9936, 1.0466]; 5/6 | 0.9782 [0.9356, 1.0227]; 2/6 |
| 1MiB same-arena independent LCP | 1.0203 [1.0032, 1.0376]; 5/6 | 1.0200 [0.9957, 1.0448]; 4/6 | 1.0114 [0.9991, 1.0238]; 5/6 |
| Concat, original growth control | 0.8520 [0.8316, 0.8730]; 0/6 | 0.8450 [0.8273, 0.8630]; 0/6 | 1.0008 [0.9937, 1.0079]; 1/6 |

## Signals and limits

For same-arena independent 1MiB LCP, pinned current/prior ratios are 1.0203
[1.0032, 1.0376] in normal slots (five of six slower) and 1.0200
[0.9957, 1.0448] after swapping slots (four of six slower). The normal-slot
secondary timed mean, 1.0194 [1.0019, 1.0371], and wall estimate,
1.0221 [1.0017, 1.0430], have the same positive direction. Unpinned primary
ratios are 0.9778 and 1.0266 with wider intervals, so the penalty is not uniform
across all settings. The pinned A/A point estimate is 1.0114 [0.9991, 1.0238],
with five of six slower. This contextual variation cannot be subtracted from
the measured penalty, nor does it identify a cause. A smaller pinned LCP cost
remains unresolved even though the earlier 6.8% point estimate is not uniformly
reproduced.

The 100K combined primary ratios range from 0.9787 to 1.0145 across assignments
and affinity arms, and all four intervals include one. These results do not
establish a stable combined penalty, but do not prove no smaller effect either.
The three phase-specific setup histories differ, so their times are not additive
components and cannot be subtracted to assign a causal comparison cost.
In particular, pinned swapped-slot comparison-only **timed mean** is 1.0516
[1.0182, 1.0862], with all six processes slower; its primary median is 1.0198
[0.9936, 1.0466]. Both estimators are retained. The corresponding normal-slot
primary median is 0.9917, so this positive secondary signal does not reproduce
uniformly across placement. Edit-only has a pinned normal-slot primary ratio
of 0.9834 [0.9676, 0.9994] and an unpinned normal-slot timed-mean ratio of 0.9692
[0.9401, 0.9993]; other assignments do not establish that effect consistently.

Concat's primary current/prior ratios are 0.8450–0.8642, with every comparison
pair faster and all four intervals below one. The timed-mean and wall estimates
also improve in all four comparisons. A/A concat intervals include one. This
is evidence for the existing round-six concat improvement in this diagnostic;
it does not resolve the LCP or churn mechanism or guarantee all concat shapes.

The separate [12-pair original-executable replication](original_replication/summary.json)
retains all original cell sequences and yields 1.0200 [0.9661, 1.0769] for 100K
independent exact churn (eight of twelve slower), 0.9933 [0.9387, 1.0511] for
1MiB independent LCP (seven of twelve slower), and 0.8476 [0.7657, 0.9383] for
concat (one of twelve slower). The earlier round-six initial and confirmation
cohorts remain separate historical evidence. No diagnostic cohorts are pooled.
Differences between original and controlled contexts, and between order,
affinity, estimator and slot assignment, remain observations rather than an
identified thermal, allocator or code-placement cause.

## Measurement and verification

All 36 fixed processes completed successfully: three slot builds × six process
seeds × two affinity modes, **720 measured blocks and 153,552 retained raw
samples**. Child-process wall time totals 121.7 seconds. Every matched case used
one common calibrated session count in both arms and all four ABBA/BAAB blocks.
Each repeated churn seed appears in both slot orders; case position is shuffled
independently of that order. No failed process, sample or unfavorable case was
removed or replaced.

Concat retains three warmup calls and seven samples of 2,000 calls per fresh
arena. LCP retains three warmup calls and seven samples of 100 calls per session.
Their primary session statistic is the median of seven per-call times; the block
statistic is the median of session medians. Churn uses median per-edit session
time. Each process takes the median of its two block statistics per slot before
forming a paired ratio. Accumulated timed means and complete wall estimates are
separate secondary metrics: all **90 estimate vectors**, individual process
ratios, ranges and counts are available in
[summary.json](controlled_results/summary.json). Full session durations, answers,
node/address metadata, statuses and raw hashes are retained in
[results.json](controlled_results/results.json) and adjacent raw TSV/sidecar files.

Calibration targets approximately 25 ms of accumulated timed work, with at least
one complete session. This is not a guaranteed minimum: the shortest block was
11.8 ms, while median block durations were 23.0 ms combined churn, 23.9 ms edit
only, 23.5 ms comparison only, 22.9 ms concat and 40.0 ms LCP. Every LCP block
contained one complete seven-by-100 session; its shortest timed duration was
36.1 ms. Setup, byte oracles, telemetry, output and destruction occur outside
local timers, but may affect cache, allocator, frequency and thermal history.

The 100K generator/50 edits, independent target arena, representative edit warmup,
25-match bitmask and final byte oracle match the original workload. Comparison
only prebuilds and validates every edited state. The 1MiB case preserves the full
original local arena construction, including early/late/shared-late roots, the
target's LCP/equality/hash-prefix checks, independent node and leaf identities,
and exact byte checks. The earlier full-sequence same-root timed cells are not
replayed. Additional materialization and mismatch validation change prior access
history. Therefore this focused harness complements, rather than duplicates,
the original-executable recurrence. Same-arena LCP and cross-arena churn use
different exact traversal paths; results cannot simply be transferred between
them. Fingerprint-collision rejection is checked in both cores outside timing.

Pinned complete-block endpoints remain CPU 0 (zero of 360 endpoint changes).
Unpinned blocks have 199/360 endpoint changes. Endpoints include setup; matching
endpoints do not exclude intervening migration or other machine-state effects.
Addresses describe placement but do not prove a placement mechanism. No counters
or capacity instrumentation were added to timed cores. The summary reparses all
raw samples and verifies context, sidecars, commands, exit status, source/lock/
archive/executable identities and expected counts before accepting the results.
An [independent audit](controlled_review.json), using a separate
[recomputation script](controlled_review.py), also verified all 36 processes,
720 blocks, 153,552 samples and 90 estimate vectors, along with the saved
source/archive/binary/raw identities. The original-executable recurrence was
[audited separately](original_replication_review.json); its data are not pooled.

Protocol and implementation: [PLAN.md](PLAN.md), [controlled.rs](controlled.rs),
[controlled.py](controlled.py). Previous evidence and mechanism limits:
[mechanism_audit.md](mechanism_audit.md). Any round-seven candidate and its final
unchanged-workload comparison are separate experiments; these diagnostics stay
frozen and do not establish a universal no-regression or causal-noise claim.
