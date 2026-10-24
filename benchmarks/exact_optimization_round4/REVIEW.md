# Independent audit of the round-three comparison timings

Date: 2026-10-24. This independent review audits the preserved round-three
sources and raw measurements, the completed round-four diagnostic follow-ups,
and the selected no-production-change outcome. It introduces no production
changes or additional timed runs. No causal explanation for the original
timing increases is established by this audit.

In this document, **baseline** and **optimized** identify the two round-three
cohorts: the completed round-two implementation and the completed round-three
implementation, respectively. The new round-four baseline is the latter.

## What the original measurements establish

The round-three report correctly preserves two slower comparison-only medians:
independent 1 MiB equality, 58.656 to 68.446 microseconds (+16.69%), and
last-byte-mismatch LCP, 60.029 to 67.337 microseconds (+12.17%). These are medians
of three process/round medians. Each process median summarizes seven sequential
batches of 100 comparisons on the same immutable fixture.

Both affected cells are slower in all three paired rounds, but the size of the
observed difference varies considerably:

| Round | Process order | Independent equality, baseline -> optimized (us) | Change | Late LCP, baseline -> optimized (us) | Change |
| --- | --- | ---: | ---: | ---: | ---: |
| 0 | baseline, optimized | 59.864 -> 68.446 | +14.34% | 60.029 -> 68.763 | +14.55% |
| 1 | optimized, baseline | 58.656 -> 68.978 | +17.60% | 66.136 -> 67.337 | +1.82% |
| 2 | baseline, optimized | 58.386 -> 58.779 | +0.67% | 59.097 -> 61.213 | +3.58% |

The original statistic should not be replaced after seeing these results.
However, three fresh processes per core do not establish a stable 12--17%
penalty across executions or identify its cause. The 21 within-process samples
are not 21 independent process replicates. Overlapping round ranges likewise
do not prove equivalence or justify dismissing the increases as noise.

Other cells do not show a uniform comparison penalty. Independent equal LCP
process medians are 64.310/60.304/58.227 us for the baseline and
57.392/57.467/60.818 us for the optimized core. Late-mismatch equality is slower
in two rounds and faster in the third. The independently measured edit-and-
compare workload improves; it uses different data, construction and comparison
paths, so these results cannot be substituted for one another.

## Retained raw sample observations

Each list below preserves the original temporal order of the seven batches,
in microseconds per call. Full data remain in
[`round-three microbench.json`](../exact_optimization_round3/results/microbench.json).

| Cell | Round / core | Seven batch samples (us/call) |
| --- | --- | --- |
| Independent equality | 0 / baseline | 63.040, 58.763, 57.884, 57.794, 59.864, 60.480, 64.483 |
| Independent equality | 0 / optimized | 57.973, 56.769, 71.150, 68.446, 68.878, 70.019, 58.710 |
| Independent equality | 1 / baseline | 60.229, 58.656, 71.142, 63.906, 56.992, 56.519, 57.681 |
| Independent equality | 1 / optimized | 57.069, 71.122, 69.502, 66.879, 69.482, 68.978, 66.159 |
| Independent equality | 2 / baseline | 57.636, 58.386, 58.977, 59.752, 57.707, 56.393, 59.356 |
| Independent equality | 2 / optimized | 56.899, 57.652, 64.572, 57.234, 58.779, 66.517, 72.147 |
| Late LCP | 0 / baseline | 67.239, 59.335, 59.953, 72.298, 65.568, 60.029, 59.808 |
| Late LCP | 0 / optimized | 73.644, 71.595, 71.504, 68.763, 67.743, 66.211, 66.383 |
| Late LCP | 1 / baseline | 73.689, 67.375, 65.331, 65.892, 65.215, 66.136, 75.862 |
| Late LCP | 1 / optimized | 69.366, 67.474, 67.337, 66.735, 66.694, 68.655, 66.634 |
| Late LCP | 2 / baseline | 59.097, 60.451, 58.097, 59.556, 58.034, 61.299, 58.034 |
| Late LCP | 2 / optimized | 61.867, 61.213, 63.074, 59.636, 67.310, 59.125, 57.300 |

For example, optimized round-zero equality moves from about 57 us to about
69 us and back within a single cell. The source data have no timestamps,
processor/core IDs, frequency, thread CPU time, context-switch, or cache-event
telemetry that could explain this behavior. A thermal, scheduler, allocator,
cache, or generated-code explanation remains a hypothesis.

## Harness and exact-path facts

The original harness is
[`bench_rust.rs`](../comparison_safety/bench_rust.rs), built by
[`bench_rust.py`](../comparison_safety/bench_rust.py). Its order is fixed:

1. Construction/hash/edit primitives run first.
2. The 64 KiB comparison fixtures precede the 1 MiB fixtures.
3. Within a size, the fixture order is same-root, independent equal, early
   mismatch, late mismatch, then shared late mismatch.
4. Within each fixture, all seven LCP batches precede all seven equality
   batches, followed by all seven hash-only LCP batches.
5. Each operation receives only three untimed warmup calls. The affected timed
   batches each contain 100 calls, lasting approximately 6 ms.

The matrix starts fresh processes in baseline/optimized, optimized/baseline,
baseline/optimized order. Each Rust pair follows that round's two Python
microbenchmark processes. Order is counterbalanced only across these three
pairs, not within cells, and there is no logged affinity/frequency control.

The source change that makes hashing faster also changes how much work happens
before the comparator cells. Summing the recorded batch durations preceding the
1 MiB independent LCP gives 234.6--236.8 ms for the baseline and 48.7--51.2 ms for
the optimized core. Before independent equality, it gives 275.3--281.4 ms and
89.9--94.0 ms, respectively. These are **sums of recorded timed intervals**, not
wall-clock timestamps: they exclude construction, warmups, oracle checks and
output. Nevertheless, they demonstrate that unchanged harness order does not
mean identical preceding CPU work after the hash optimization.

Earlier benchmark arenas remain alive while comparison fixtures are constructed.
Removed temporary leaf clones also alter allocator history before those fixtures
are allocated. Same logical node layout and IDs do not prove identical physical
heap addresses/cache placement. The two builds may additionally differ in
generated instruction layout despite unchanged exact Rust function bodies.
Neither possibility is established as the cause by the available artifacts.

Independent microbenchmark roots are separately allocated in **one arena**, with
512-byte leaves and repeated byte values 0 through 255. They use the bounded
recursive same-arena fast path. In contrast, the independent-target churn
workload uses one million edited ASCII bytes and separate target arenas; its
long comparisons use the frontier path. These distinctions preclude attributing
the difference between their measured times to editing alone.

`Arena::equal_across` checks lengths, then calls the same `lcp_across` engine used
for LCP. On these equal-length inputs, no fingerprint is inspected and no rope
is mutated. Thus the source does not contain an equality-specific byte scan
that would by itself explain the larger equality result. Call-site inlining,
code generation and measurement state still require checking rather than an
assumption that the compiled executions must be identical.

## Controlled diagnostic design review

The proposed diagnostic uses both immutable cores in one executable, sustained
warmup, paired ABBA/BAAB blocks, counterbalanced fixture and LCP/equality order,
and six fresh processes. A second executable swaps the core-to-crate-slot mapping
to test sensitivity to linked-code roles/layout. Inputs and plain-byte oracles
remain outside timing; measured aggregate answers are checked.

This addresses the main original-order confounders. Review recommendations were:

- Use a common calibrated iteration count for both arms of a matched
  fixture/operation, recording it and using it for all paired blocks.
- Keep input handles and outputs behind `black_box`; record any function-pointer
  or wrapper dispatch that could dominate tiny comparisons.
- Record processor IDs and thread CPU time outside timed intervals when
  available. Endpoint IDs do not prove there was no migration within a block.
- Retain an original full-sequence, separate-binary control. Dual linking changes
  the original compilation context: disappearance of the difference there would
  weaken a general-regression claim but would not disprove a layout/startup/order
  effect in the original binaries.
- Inspect generated comparator code and retain the source/binary provenance.
  Do not infer a compiler effect solely from different executable hashes.

There is no material design blocker in that proposal.

### Implementation and oracle review before launch

The reviewed implementation uses three fixed crate-slot assignments: old/new,
new/old, and new/new as an identical-source control. Six fresh processes per
assignment run in both unpinned and pinned modes. The build order covers all six
permutations, and pinned/unpinned order alternates. Actual process affinity masks
and block-boundary processor IDs are recorded. These do not measure frequency
or prove the absence of within-block migration in the unpinned mode.

Both arms share each calibrated iteration count. The operation branch is outside
the measured inner loop; arenas, node IDs, and results pass through `black_box`.
Each block's accumulated answer is checked against the expected answer times
its iteration count. Input construction, materialization, and plain-byte oracle
checks occur outside the timer. The ordinary fixtures now explicitly verify
the first input against its source bytes and the second input's full length.
The compressed fixtures and concrete collision have separate exact-answer
assertions. The driver rejects incomplete/malformed block groups and requires a
fresh results directory.

Review found an ordering bug before any timings ran: the original combination
of API-order and arm-order parity always assigned LCP to ABBA and equality to
BAAB. It is corrected: API order now depends on seed plus stable fixture index,
while the initial arm depends on seed plus API identity. Consequently every
fixture/API gets exactly three ABBA and three BAAB sequences over seeds 1--6,
independently of its position in the schedule. The implementation's static
distribution check also confirms three LCP-first and three equality-first
positions for each fixture. Unused pre-review build artifacts are preserved;
all three binaries were rebuilt under `.cache/controlled_v2` before launch.

No remaining prelaunch correctness, balancing, or provenance blocker was found
in the source/driver review. Equivalent-code folding in the identical-source
control and any effect of dual linking must be considered when interpreting its
results; neither an equal result nor a difference alone establishes a cause for
the original separate-binary observation.

## Controlled result review

The completed diagnostic retains 36 fresh-process runs and 3,456 timed blocks.
Independent lightweight checks verified all 36 raw TSV hashes and every block's
aggregate answer against its fixture oracle. Recomputing all 864 process ratios
and 144 summary groups from the block records reproduced the stored geometric
means and exploratory log-t intervals. Every fixture/API/assignment/affinity
group has three ABBA and three BAAB process sequences.

The two primary cells are **independent equality** and **last-byte LCP**, each
on the original 1 MiB fixture. The following values are round-three/round-two
geometric mean ratios with exploratory per-cell 95% log-t intervals. Each row
contains six process pairs; intervals are not multiplicity adjusted.

| Affinity / crate assignment | Independent equality | Last-byte LCP |
| --- | --- | --- |
| Unpinned / old,new | 1.0110 [0.9773, 1.0459] | 0.9955 [0.9685, 1.0232] |
| Unpinned / new,old | 1.0100 [0.9566, 1.0664] | 0.9882 [0.9644, 1.0125] |
| Pinned / old,new | 0.9936 [0.9697, 1.0181] | 1.0374 [0.9804, 1.0976] |
| Pinned / new,old | 0.9868 [0.9629, 1.0113] | 1.0089 [0.9974, 1.0206] |

The original 16.69% independent-equality and 12.17% last-byte-LCP penalties do
not recur in this controlled diagnostic. This is narrower than an equivalence
claim: the pinned old,new last-byte-LCP interval extends to a 9.76% penalty.
The identical-source control's slot-B/slot-A ratios are 0.9933 and 0.9974 for
independent equality (unpinned/pinned), and 0.9979 and 0.9833 for last-byte LCP.
These controls also vary; their results cannot certify universal equality of
execution time or remove a smaller source/layout-conditioned effect.

One secondary cell, pinned new,old **last-byte equality**, has ratio 1.0315
[1.0191, 1.0440], with all six process ratios above one. The corresponding
old,new equality ratio is 0.9963 [0.9456, 1.0496]. This association with crate
assignment is worth preserving, but it does not isolate generated code layout
from allocation state or establish a general version effect. It is also a
different API from the primary last-byte-LCP cell.

The endpoint telemetry records 477 differing processor IDs among 1,728 unpinned
blocks and none among 1,728 pinned blocks; all pinned endpoints are CPU 0.
This confirms the affinity intervention at the sampled boundaries, not a causal
explanation of the original regression. The retained original separate-binary
follow-up below addresses dual linking and the changed fixture schedule and
warmup procedure.

## Retained original-binary replication and selection

The fixed follow-up protocol ran all 12 fresh process pairs, alternating their
order to give six baseline-first and six optimized-first pairs. It used the two
retained original round-three executables without rebuilding or changing their
Rust workload, warmup, or fixture sequence. The ordinary unpinned affinity mask
was retained. Unlike the historical full driver, this follow-up did not
interleave Python microbenchmark processes; that workload-history difference is
explicit in the protocol.

Independent review verified all 24 raw TSV hashes, both executable hashes
against the original round-three JSON, all 1,008 within-process seven-sample
medians, all 504 paired ratios, and all 42 summary geometric means and log-t
intervals. The inherited executable checks the expected exact answers before
timing; this follow-up retains those checks, rather than introducing the
controlled diagnostic's per-block aggregate checks.

| Primary 1 MiB cell | Geometric mean round-three/round-two ratio | Exploratory 95% log-t interval | Range of 12 pair ratios |
| --- | --- | --- | --- |
| Independent exact equality | 0.9959 | [0.9458, 1.0486] | [0.8436, 1.1682] |
| Last-byte exact LCP | 0.9914 | [0.9326, 1.0540] | [0.8158, 1.2226] |

The reported 12--17% penalties do not replicate in either this original-binary
follow-up or the controlled primary cells. Individual pair ratios still vary
substantially. These observations support retaining the existing implementation
and recording the original slower results with their uncertainty; they do not
identify a causal noise source or prove universal performance equivalence.
These two diagnostic protocols alone do not identify an actionable cause. The
complete-matrix result below qualifies their interpretation: equality becomes
slower again under the short full protocol.

The selected outcome is **no production change in round four**. All 16 current
Rust source/test files are byte-for-byte identical to the captured round-four
baseline. The full-protocol matrix compares the immutable pre-round-three
implementation with this unchanged current implementation; it is a repeat of
the earlier optimization comparison, not a new-code gain or an A/A experiment.

## Full-matrix qualification and completed audit

The full matrix is complete. Its original short microbenchmark protocol again
reports slower independent equality. These observations must be retained
alongside the near-unity estimates from the two diagnostic protocols.

| Primary 1 MiB cell | Reference median [minimum, maximum], microseconds | Current median [minimum, maximum], microseconds | Change in reported medians |
| --- | --- | --- | --- |
| Independent exact equality | 54.926 [54.543, 58.861] | 65.895 [55.910, 65.990] | +19.97% |
| Last-byte exact LCP | 57.184 [54.632, 57.449] | 60.059 [57.426, 65.809] | +5.03% |

These are medians of three round medians; the ranges are observed round ranges,
not confidence intervals. The three equality process-pair ratios are 1.2081,
1.1211, and 1.0179, so its direction is consistently slower in this short
replication while its magnitude varies. The late-LCP pair ratios are 0.9996,
1.1508, and 1.0993. Their paired geometric means, 1.1130 and 1.0814, respectively,
are different statistics from the table's ratios of medians and must not be
substituted for them without disclosure.

Within-process variation remains visible: in current equality rounds zero and
one, the first three batches take roughly 56--57 microseconds, followed by four
batches at roughly 66--73 microseconds. This describes the samples, not their
cause. Other secondary negative cells also remain in the full report: 64 KiB
last-byte LCP changes from 3.082 [3.069, 3.244] to 3.695 [3.612, 3.961]
microseconds, and equality from 3.083 [3.078, 3.089] to 3.979 [3.963, 4.013].
They have nonoverlapping observed round ranges, without that fact establishing
a general performance bound or an isolated cause.

The final run restores Python/Rust outer interleaving and uses rebuilt
executables. The final binary hashes, recorded in the raw microbenchmark JSON,
are `f4d7546d041bb450540a9136ea9c614ee0cb8e7aac7c60d0a1708df13bc765d2`
and `ef2bf3295560acbbc7331951fc0c57faf0170d294c0d55307a0f40241280fca2`.
They differ from the retained original-binary pair. Neither outer history nor
linked build context was isolated experimentally. The original equality signal
therefore remains unresolved; the diagnostic results do not erase it or justify
calling it noise. No production change was selected without an isolated cause
and a validated remedy. Existing fingerprints and exact semantics remain
preserved by retaining the tested implementation.

Independent final data checks verified:

- All 490 summary statistic vectors from their retained values, including means,
  sample standard deviations, medians, minima, maxima, and counts.
- All 510 within-process microbenchmark medians against their seven raw samples,
  and every aggregated microbenchmark vector against its three process medians.
- Both biology cohorts' nine-sample medians and recorded byte-oracle success.
- All 95 additional extended, biology, and microbenchmark report rows against
  the completed summary. The previously checked six Rust README rows and 68
  paper mean/SD cells remain consistent.
- All five canonical context fingerprints, 57 worker-to-context bindings,
  31 distinct current harness/helper hashes, and all 20 carried-forward tested
  production/Rust-test file hashes.

The final audit records successful oracle/provenance checks for all cohorts,
unchanged preservation of 376 earlier benchmark artifacts, alternate-output
summary identity, and verified biology resume without rewriting results. Its
complete-matrix flag and five context identities match the final summary.
These checks establish data and source consistency; they do not turn the timing
differences into causal estimates. The interpretation preserves unchanged
Python controls, candidate-versus-exact distinctions, excluded query construction,
and inline layout versus total allocated memory.

## Provenance and validation scope

| Artifact | SHA256 |
| --- | --- |
| Round-four baseline archive | `3047745ce88cb5b8a69c37009711d7a7662443ed4afd0eafec12f2f841af2a5b` |
| Original round-three raw microbenchmark JSON | `7b0945458a12e0a4ba3b6915e81789ac8844af1a4cb6d5c476a25615d343881e` |
| `benchmarks/comparison_safety/bench_rust.rs` | `7f0fcf1c0c9da53767eb116973f92bf605e6c860fee80895c845e6f166e63a70` |
| `benchmarks/comparison_safety/bench_rust.py` | `f14e06ecf535bc075eff402b8c83b033215ea5402f8afb76e1b6f7825d510d84` |
| Round-three `run_matrix.py` | `598ee53aede16e40365277037965348383a2a20aad33f496c883af4ceefe4f6e` |
| Reviewed `compare_controlled.rs` | `072b4b5cba70cdb0c3782b38f2ba7ab4ad0b05c64acb3b2db34707affc8071d6` |
| Reviewed `compare_controlled.py` | `f5c945d60b480a89a2173914bbd906e1eca26db5ed2a03b506553d7ffc837e2b` |
| Controlled initial `results.json` | `8fa307aa739246cc5ec62a67470d847df0b40d42df0de7a445d61be4b75b8a0a` |
| Original-binary replication `results.json` | `334119fdf9c7645236415719b4c2f877f3baee450635d0016f00518c92aae8b8` |
| Final full-protocol `results/microbench.json` | `f5b5930bc8fb1faf5585f44e51d3a6e7bca1d36bd89ab88c35d435688bad0c19` |

No production change was selected. The existing **658 selected-test validation
is carried forward, not rerun in round four**. The archive and working-file
identity checks in [carried_validation.json](carried_validation.json) cover all
20 recorded production and Rust integration-test files from the previously
tested current implementation. The independent Rust-specific check above
verified its 16 source/test files directly. Fresh diagnostic oracle, raw-data,
and summary checks are recorded separately in this review; they are not added
to the carried-forward suite count.

The completed Rust and paper-control portions of the full matrix have also
been reviewed against `summary.json`: all six Rust README rows and their
ratios match, as do all 68 displayed paper mean/SD cells. Recomputed paper
log-log slopes are 0.9628229512 and 0.9984076843, both failing the historical
0.3 gate. The nine-point dispersion-band separation starts at the smallest
measured point, so it does not identify an interior crossover. Python remains
unchanged between the measured cores; its timing variation is a control, not
an effect of the Rust optimization. The final whole-matrix completion and
provenance review is recorded above. The root and Rust README updates preserve
the fresh equality recurrence and the distinct historical diagnostic figures.

Final read-only review of the paper's fourth-round section in `paper/main.tex`,
the updated `CLAIMS.md` rows, and `findings.md` found no material numerical or
interpretive mismatch with the completed data. They retain the full-protocol
primary and secondary slower cells, diagnostic uncertainty, rebuilt-binary and
outer-interleaving differences, unchanged Python controls, excluded query
construction, inline-layout limits, and carried-forward validation. Earlier
third-round measurements remain explicitly historical. The parent separately
confirmed successful compilation of the revised paper.
