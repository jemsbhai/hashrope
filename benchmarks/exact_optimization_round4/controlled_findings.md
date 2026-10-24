# Investigation of the comparison-only timing changes

The comparison timing changes remain **sensitive to the measurement protocol**.
The fresh full matrix again shows slower 1 MiB independent equality under the
original short benchmark (+20.0% by its reported median), while the controlled
paired experiment and 12-pair rerun of retained original executables give
near-unity estimates. The full-matrix late-LCP difference is +5.0%. These results
do not support a blanket claim that the regression disappeared or was merely
noise. No cause or corrective production change has been established. Production
source and the original benchmark remain unchanged; the separate full-matrix
qualification below is essential to interpreting the diagnostic findings.

The two primary cells are 1 MiB independent exact equality and 1 MiB
late-mismatch exact LCP. The historical report used the median of three process
medians: 58.656→68.446 µs (+16.7%) for independent equality and
60.029→67.337 µs (+12.2%) for late LCP. Those original observations remain in
the [round-three raw results](../exact_optimization_round3/results/microbench.json).
The two diagnostic analyses below use matched process ratios, rather than
treating inner timing samples as independent replicates. The final full matrix
retains its original median-of-three reporting convention.

## Controlled paired experiment

The [predefined plan](controlled_plan.md) links immutable round-two and
round-three cores into one executable, repeats with their crate slots swapped,
and includes a round-three/round-three A/A control. Each of these three builds
runs in six fresh processes, separately unpinned and pinned to logical CPU 0:
36 processes and 3,456 timed blocks in total. All exact-result checks passed.

Fixture order, construction order, API order, and ABBA/BAAB order are
counterbalanced. A prelaunch review caught a coupling between API order and
block order; the reviewed v2 harness fixes it and gives each fixture/API exactly
three ABBA and three BAAB processes. The earlier build/source archives are
explicitly retained as `pre_review`; they produced no timed results.

Each matched pair uses a common calibrated iteration count, sustained paired
warmup, and blocks targeting 10 ms. The target is not a hard minimum: actual
blocks range from 1.60 to 21.25 ms because calibration and runtime state vary.
Inputs and outputs are black-boxed. Ordinary fixtures use byte oracles; the
compressed fixtures verify their known constructed periodic contents, lengths,
and exact answers. A concrete default-hash collision is checked untimed.

Values below are geometric means of **six process-paired round3/round2 ratios**,
with exploratory per-cell 95% log-t intervals. A ratio above one means round
three took longer. These intervals are not multiplicity adjusted and are not
universal performance bounds. The two slot builds and affinity arms remain
separate; none are pooled to conceal their differences.

| Affinity / core slots | Independent equality | Late-mismatch LCP |
|---|---:|---:|
| Unpinned, round2 / round3 | 1.0110 [0.9773, 1.0459] | 0.9955 [0.9685, 1.0232] |
| Unpinned, round3 / round2 | 1.0100 [0.9566, 1.0664] | 0.9882 [0.9644, 1.0125] |
| Pinned, round2 / round3 | 0.9936 [0.9697, 1.0181] | 1.0374 [0.9804, 1.0976] |
| Pinned, round3 / round2 | 0.9868 [0.9629, 1.0113] | 1.0089 [0.9974, 1.0206] |

The largest primary-cell upper endpoint is **1.0976 for pinned late LCP**;
the smaller equality bounds must not be generalized to LCP. The point estimates
do not reproduce the earlier 12–17% effect across slot assignments and affinity
arms. They still permit smaller effects under these conditions.

The A/A control compares slot B with slot A, both running the same round-three
source:

| Affinity / identical cores | Independent equality | Late-mismatch LCP |
|---|---:|---:|
| Unpinned | 0.9933 [0.9544, 1.0338] | 0.9979 [0.9753, 1.0211] |
| Pinned | 0.9974 [0.9673, 1.0284] | 0.9833 [0.9415, 1.0270] |

Smaller secondary effects remain. For example, pinned swapped-slot
late-mismatch **equality** is +3.15%, with all six process ratios above one;
the normal-slot pinned counterpart is −0.37%. This is a separate API cell from
the primary late-LCP result. It does not support a slot-independent 12–17%
penalty, but it cautions against claiming absolute runtime equivalence.

Unpinned block endpoints report different logical processors in 477/1,728
blocks. Pinned endpoints are always CPU 0, with zero endpoint changes. Endpoint
sampling cannot detect every intervening migration; pinning does not eliminate
frequency changes, interrupts, or sibling-thread contention. Root-node and
first-leaf addresses are recorded as placement telemetry, not causal proof.

All process ratios, raw blocks, addresses, processor IDs, affinity masks,
source/binary hashes, and fixture metadata are in the
[controlled results](results/controlled_initial/results.json).

## Retained original-binary replication

Dual-linking changes code layout, so the second experiment uses the **exact two
original executables**, recovered from their temporary build directories and
verified against the original round-three SHA-256 records. Their copied Rust
harness hashes also match the unchanged original harness. The executable bytes
are preserved in
[retained_binaries.zip](results/original_replication/retained_binaries.zip).

The [replication protocol](original_replication_plan.md) was saved before this
follow-up started: 12 fresh unpinned process pairs, six old-first and six
new-first, with the entire original Rust sequence and its oracles unchanged.
There was no early stopping or sample removal. All 24 processes completed,
retaining 42 cells × seven samples each, or 7,056 timed samples. Python
microbenchmark runs were not interleaved in this follow-up, so the historical
outer-driver workload history is not exactly reproduced.

| Primary 1 MiB cell | Baseline mean of process medians, µs | Round-three mean, µs | Geometric mean paired ratio [exploratory 95% interval] |
|---|---:|---:|---:|
| Independent exact equality | 59.532 | 59.305 | 0.9959 [0.9458, 1.0486] |
| Late-mismatch exact LCP | 60.803 | 60.284 | 0.9914 [0.9326, 1.0540] |

The same original binary layout therefore also fails to reproduce the earlier
12–17% effect in this replication. The intervals describe this 12-pair study;
their upper endpoints do not guarantee a global maximum penalty. Individual
late-LCP ratios range from 0.816 to 1.223. Its old-first and new-first geometric
means are 0.9703 and 1.0131, respectively, which further cautions against a
causal interpretation from a few short processes.

Other cells are retained rather than omitted. For example, 1 MiB same-root LCP
is about +6.3% in this replication, an absolute change of approximately
0.11 ns, while the self-concatenation DAG equality cell is about −7.4%.
Tiny-case timings include loop/call overhead and timer granularity effects;
these observations do not establish an actionable comparator optimization.
The [full replication data](results/original_replication/results.json) contains
all 42 cell summaries and every process pair.

## Post-matrix qualification: the short-protocol equality slowdown recurs

After these diagnostics, the complete original outer protocol was rerun using
the immutable pre-round-three core and unchanged current core. Its
[raw microbenchmark results](results/microbench.json) and
[aggregate summary](summary.json) show:

| Primary 1 MiB cell | Baseline median [min, max], µs | Current median [min, max], µs | Change in reported medians |
|---|---:|---:|---:|
| Independent exact equality | 54.926 [54.543, 58.861] | 65.895 [55.910, 65.990] | +20.0% |
| Late-mismatch exact LCP | 57.184 [54.632, 57.449] | 60.059 [57.426, 65.809] | +5.0% |

These are medians of three process medians and their observed ranges, not
confidence intervals. Equality is slower in all three matched pairs, although
the magnitude varies. Thus the direction of the equality slowdown recurs in
the original short protocol; the diagnostic near-unity estimates do not negate
that observation. The late-LCP increase is smaller than the original +12.2%.

This final run uses newly rebuilt binaries whose hashes differ from the retained
original pair, and restores the outer driver's Python microbenchmark
interleaving. Neither difference was isolated experimentally. The experiments
therefore establish unresolved protocol/build-context sensitivity, without
identifying a binary-layout, thermal, scheduling, or cache mechanism as its
cause. All diagnostic results above remain valid for their stated protocols;
they are not a general absence-of-regression result. The
[full round-four report](README.md) presents this qualification with the rest
of the final matrix.

## Source and generated-code evidence

The immutable public exact APIs and the entire comparison engine are
byte-identical after newline normalization; the build manifest records the
extracted source hashes. Round three changes default-profile hashing and three
leaf-copy paths, which can alter preceding work and heap placement without
changing the exact comparison algorithm.

Separately emitted release library assembly was inspected for seven comparator
routines. After normalizing crate-slot names, symbol hashes, local labels, and
anonymous relocation labels, the instruction text matches. The remaining
dispatcher text difference is a section ordinal, 45 versus 46. This is **not**
a final linked-binary identity proof: relocation target contents and final
placement are not equated. The preserved
[assembly audit](assembly_audit.json), [difference](assembly_differences.txt),
and [assembly archive](comparison_assembly.zip) document the comparison.

The original harness warms each operation three times, measures only 100 calls
per long-comparison sample, and orders LCP, equality, then hash-only work. Round
three greatly accelerates preceding construction and hashing. Combined with
the observed process/sample variation, this supplies plausible confounders;
neither CPU temperature nor frequency was measured, so no specific thermal or
cache mechanism is established as the cause.

## Decision and reproducibility

The original short-protocol equality slowdown remains an unresolved finding.
The paired controls and retained-binary replication do not identify a stable
comparator change or an actionable cause, so they do not justify a speculative
production optimization. Retain the exact implementation and its compatibility
guarantees while reporting all three protocols together: the full matrix's
recurring equality slowdown, the near-unity diagnostic estimates, and their
different build and workload contexts. No production change was selected, no
cause was proven, and no universal no-regression claim follows.

The controlled harness and driver are
[compare_controlled.rs](compare_controlled.rs) and
[compare_controlled.py](compare_controlled.py); copied source and compiler
provenance are in [controlled_sources.zip](controlled_sources.zip) and
[controlled_build.json](controlled_build.json). The follow-up uses
[replicate_original.py](replicate_original.py). Both timing drivers refuse to
overwrite existing result directories. The controlled driver accepts a fresh
`--output` directory; another original-binary replication should use a separately
recorded driver copy with a fresh output path. Immutable cores and all earlier
studies remain preserved.

[REVIEW.md](REVIEW.md) records independent verification of the prelaunch
ordering correction, 36 controlled raw-file hashes and all 3,456 aggregate
answers, the 24 replication raw hashes and executable identities, and all
1,008 replication medians, 504 paired ratios, and 42 interval summaries.
