# Controlled diagnosis of the round-five regressions

The longer controlled diagnostic did not reproduce a consistent 5.1% independent
1M churn penalty or 8.4% concat penalty across both slot assignments and affinity
modes. This is a protocol-specific result, not a proof that the earlier changes
had no performance cost or that every difference was noise. The separately
retained original-executable replication still estimates concat at 1.0804, with
substantial uncertainty, so its positive direction remains part of the evidence.
No historical result is replaced or pooled into this controlled cohort.

Here **current** means the immutable round-five candidate (`baseline.zip`), and
**released** means the immutable released 0.3.3 reference. Neither term denotes a
subsequently proposed round-six optimization. The full source, protocol, lock and
executable identities were frozen before measurement in
[controlled_freeze.json](controlled_freeze.json),
[controlled_build.json](controlled_build.json),
[controlled_sources.zip](controlled_sources.zip) and
[controlled_binaries.zip](controlled_binaries.zip).

## All controlled estimates

Each entry is the geometric mean of six fresh process ratios, its exploratory
95% log-t interval, and the number of processes whose numerator was slower.
Released/current and current/released columns both report **current / released**;
the two headings describe slot A/B placement. The identical-code control reports
**slot B / slot A**, both containing current code. Ratios above one are slower.
Intervals are not adjusted for multiple comparisons or used as equivalence tests.
Every primary median-statistic interval in the table includes one; this alone
is not evidence of zero effect. Some secondary timed-mean and wall estimates
below exclude one. Inner sessions and blocks are not counted as independent
processes.

### Unpinned

| Case | Released/current slots | Current/released slots | Current/current A/A |
|---|---:|---:|---:|
| Concat, original growth | 1.0050 [0.9588, 1.0534]; 4/6 | 1.0195 [0.9630, 1.0792]; 4/6 | 0.9938 [0.9348, 1.0565]; 2/6 |
| Concat, prewarmed | 0.9929 [0.9660, 1.0206]; 4/6 | 0.9984 [0.9821, 1.0151]; 3/6 | 1.0082 [0.9746, 1.0431]; 2/6 |
| 1M edit + exact compare | 1.0177 [0.9837, 1.0529]; 4/6 | 0.9622 [0.8927, 1.0372]; 2/6 | 1.0159 [0.9176, 1.1247]; 3/6 |
| 1M edit only | 0.9971 [0.9617, 1.0338]; 3/6 | 1.0000 [0.9866, 1.0135]; 3/6 | 0.9982 [0.9869, 1.0097]; 2/6 |
| 1M exact compare only | 1.0079 [0.9936, 1.0223]; 4/6 | 0.9899 [0.9648, 1.0157]; 1/6 | 0.9705 [0.9207, 1.0230]; 2/6 |

### Pinned

| Case | Released/current slots | Current/released slots | Current/current A/A |
|---|---:|---:|---:|
| Concat, original growth | 1.0014 [0.9981, 1.0048]; 4/6 | 1.0070 [0.9886, 1.0256]; 1/6 | 1.0000 [0.9760, 1.0246]; 2/6 |
| Concat, prewarmed | 0.9979 [0.9864, 1.0096]; 1/6 | 0.9873 [0.9662, 1.0088]; 1/6 | 1.0047 [0.9877, 1.0220]; 2/6 |
| 1M edit + exact compare | 0.9945 [0.9826, 1.0065]; 1/6 | 1.0006 [0.9791, 1.0226]; 3/6 | 0.9923 [0.9468, 1.0400]; 3/6 |
| 1M edit only | 1.0013 [0.9986, 1.0040]; 4/6 | 1.0003 [0.9969, 1.0037]; 3/6 | 0.9968 [0.9861, 1.0077]; 3/6 |
| 1M exact compare only | 1.0039 [0.9773, 1.0312]; 3/6 | 0.9925 [0.9706, 1.0149]; 2/6 | 0.9837 [0.9474, 1.0214]; 2/6 |

## Interpretation

The combined 1M churn point estimates range from 0.9622 to 1.0177 across the
four released/current comparisons. The unpinned normal-slot interval still
allows a 5.1% penalty (upper endpoint 1.0529); the other three upper endpoints
are lower. The normal and swapped unpinned estimates have opposite directions,
and the pinned estimates are close to one. Edit-only estimates are 0.9971–1.0013
and comparison-only estimates are 0.9899–1.0079. These distinct setup histories
cannot be subtracted to assign a causal number of nanoseconds to comparison.
They do not isolate a stable edit or comparison slowdown in this diagnostic.

Concat with the original arena-growth trajectory has point estimates of
1.0014–1.0195. Its largest controlled upper interval endpoint is 1.0792.
Prewarmed concat estimates are 0.9873–0.9984; their intervals also include one.
The isolated capacity replay verifies reallocations in original samples 0, 1
and 3, and zero allocator events in all seven prewarmed samples for both cores.
The prewarmed control also changes retained memory and operation history; it is
not a single-variable proof that allocation caused the historical difference.
The A/A results display nonzero variation, including wider churn intervals.
They measure sensitivity of identical code in this environment; they do not
establish that a particular released/current difference is noise.

The separate [original-executable replication](original_replication/summary.json)
retains the unchanged 42-cell micro and 24-cell churn sequences. Its six paired
estimates are 0.9836 [0.9491, 1.0195] for independent 1M exact churn (one of six
slower) and 1.0804 [0.9056, 1.2888] for concat (three of six slower). Thus the
original concat point-estimate direction recurs while the longer controlled
estimates are smaller. The protocols differ in source wrapper, binary layout,
work duration, setup, warmup, operation history and process context. This study
does not isolate which difference matters or erase the original observation.
Evidence for a stable 5.1%/8.4% penalty across all controlled assignments is absent;
smaller effects and original-context effects remain unresolved.

## Workload, estimator and checks

All 36 predeclared processes completed successfully: three slot builds × six
process seeds × two affinity modes, 720 measured blocks and **248,764 retained
raw samples**. Child-process wall time summed to 203.8 seconds. Churn fixture
seeds 1/2/3 repeat with opposite ABBA/BAAB orders. Case order is seeded, slot
order is independent of shuffled position, and each matched case uses one common
calibrated session count for both arms throughout its quartet. No process or
sample was removed or replaced.

Concat sessions preserve the complete three-warmup (or 50,003-prewarm), then
seven × 2,000-call trajectory. The primary block statistic is the median of
within-session medians of seven samples; each process then takes the median of
its two block statistics per slot. This preserves the original seven-sample
estimator rather than silently substituting a growth-weighted timed mean.
Churn blocks use median per-edit session time. Accumulated timed means and wall
per-operation estimates are reported separately in the complete
[summary.json](controlled_results/summary.json). All durations, session counts,
answers, node counts and placement telemetry remain in
[results.json](controlled_results/results.json) and the raw TSV files beside it.

Secondary estimates are not uniformly null: unpinned released/current growing
concat has a timed-mean ratio of 1.0219 [1.0088, 1.0351] and a wall ratio of
1.0171 [1.0039, 1.0305]. Two edit-only wall intervals also exclude one, including
one A/A control. These estimators answer different questions; wall time includes
setup and destruction, and the concat timed mean weights growth samples. Their
unadjusted intervals are retained rather than replaced by the primary medians,
but do not establish a causal allocator or runtime mechanism.

Calibration targets approximately 25 ms accumulated timed work in the slower
arm, using a common integer number of complete sessions, followed by paired
warmup. This is a target, not a guaranteed minimum: observed block durations
ranged down to 13.5 ms, with coarse 50-edit churn sessions and calibration
variation; the median combined and compare-only block durations were 21.3 and
19.7 ms. Concat medians were 24.9/24.7 ms and edit-only 25.0 ms. Setup, complete
byte oracles, telemetry, output and arena destruction occur outside local timers.
Their effects on cache, allocator and thermal history remain possible. This is
a diagnostic harness with additional input black-boxing and validation, not a
byte-identical copy of the original timed executable.

The exact 50-edit generator, independent targets, alternating final-byte mismatch
mask (25 matches), final edited bytes and default fingerprint-collision rejection
are checked. Every concat session verifies node increments, final length, weight,
fingerprint and both byte halves. All saved raw data were reparsed and their
estimators recomputed; sidecars, commands, exit statuses, saved context,
source/archive hashes and executable hashes were verified by the summary step.
No instrumentation or capacity API was added to either timed core.
An [independent audit](controlled_review.md) also recomputed all 248,764 raw
samples and all 90 primary/secondary estimate vectors and verified their saved
identities; its [machine-readable results](controlled_review.json) are retained.

Pinned blocks all begin and end on logical CPU 0 (0/360 endpoint changes).
Unpinned blocks change endpoint processor in 187/360 cases. Endpoint readings
span the complete block, including setup; neither matching endpoints nor pinning
rules out intervening machine-state effects or sibling-thread contention.
Root/leaf addresses are descriptive telemetry, not evidence of causation.

Protocol: [PLAN.md](PLAN.md). Harness and verifier:
[controlled.rs](controlled.rs), [controlled.py](controlled.py). Allocation replay:
[concat_capacity_probe.json](concat_capacity_probe.json). Existing-data and
mechanism limits: [mechanism_audit.md](mechanism_audit.md). Any selected round-six
production change and its unchanged full-workload comparison belong in the
separate final comparison; these controlled data remain unchanged.
