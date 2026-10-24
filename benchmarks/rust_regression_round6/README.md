# Rust concat and edit/comparison regression investigation

Retain the narrow concat optimization: it reduces paired concat time by 14.5%
in the primary cohort and 12.2% in confirmation, with all 221 release tests
passing. The original 5.1% edit/comparison penalty is not stably reproduced.
The original concat signal and secondary edit/LCP reversals remain
context-sensitive; this is not a universal performance improvement. Fingerprints,
defaults and exact-comparison semantics are unchanged.

This study follows the [round-five LCP optimization](../rust_lcp_round5/README.md).
It investigates that study's pooled +5.1% independent 1M-byte edit/comparison
estimate and +8.4% balanced-concat estimate. Round-five results remain historical
evidence; new cohorts are kept separate. Package versions remain 0.3.3 and these
optimizations are unreleased.

## Evidence and interpretation

The [original-binary replication](original_replication/summary.json) executes
the exact four archived round-five binaries, without rebuilding, for six fresh
pairs per workload. All 24 processes and 4,392 samples passed their oracles and
were retained. Independent 1M exact edit/comparison has a paired time ratio of
0.9836, faster in five of six pairs, with exploratory 95% log-t interval
[0.9491, 1.0195]. The original 5.1% penalty did not recur in this cohort. This is
evidence against treating that estimate as a stable cost, not proof of zero
effect in every context.

Concat's new ratio is 1.0804, slower in three of six pairs, interval
[0.9056, 1.2888]. Its average penalty recurs, but with substantial variation.
The [existing-data audit](mechanism_audit.md) also finds similar concat estimates
in the original and supplemental cohorts and both process orders. Simple
execution-order effects therefore do not explain the observations.

The unchanged short concat benchmark mixes growth and nongrowth operations:
its seven blocks cause 2, 1, 0, 1, 0, 0, 0 arena reallocations respectively.
An [untimed allocation replay](concat_capacity_probe.json) verifies this for
both cores. It also verifies a bounded prewarmed control with no allocations,
reallocations or deallocations during timing. Prewarming changes memory/history
as well as growth, so it is not a single-variable causal experiment. Equal
allocation trajectories do not explain the historical between-core difference.

The [retained-object inspection](codegen/README.md) finds identical normalized
instructions for the old concat/join/construction routines in the published
and round-five binaries. Final placement and runtime are not thereby identical.
Both versions redundantly check/load metadata along the balanced concat call
chain, providing a concrete optimization independent of the unresolved cause
of the historical gap.

Long cross-arena comparison source was unchanged in round five, but its emitted
dispatcher changed register allocation and gained a spill around leaf comparison.
That is a plausible mechanism, not measured attribution. A scratch-only frontier
outlining experiment changed stack/code organization and added a call; it has
no timing evidence and is not included in production.

The longer controlled run completed all 36 processes, 720 measured blocks and
248,764 raw samples. Neither original penalty persists across the two slot
assignments and affinity settings. Paired round-five/published ratios are:

| Controlled case | Unpinned, original slots | Unpinned, swapped | Pinned, original slots | Pinned, swapped |
|---|---:|---:|---:|---:|
| Combined independent 1M edit/comparison | 1.0177 | 0.9622 | 0.9945 | 1.0006 |
| Concat, original growth trajectory | 1.0050 | 1.0195 | 1.0014 | 1.0070 |

All 30 controlled primary-median exploratory intervals include 1, including the
edit-only, comparison-only, prewarmed concat and A/A controls. This is not an
equivalence result or proof that every observed difference is noise. See the
[complete controlled findings](controlled_findings.md) for intervals, all
separate arms and the distinct accumulated-time/wall estimators. Some secondary
intervals exclude 1: for example, unpinned concat in the original slots has
2.2% higher accumulated timed mean. These results remain visible and are not
substituted for the predeclared primary estimator.

The supported conclusion is narrower: **the 5.1% edit/comparison estimate is
not a stable reproduced penalty; concat's roughly 8% short-protocol signal
recurs but is context-sensitive and its cause remains unisolated.** The data
do not justify changing the exact comparison algorithm. Removing demonstrably
redundant concat work is justified independently, subject to the final rerun.

## Selected optimization and final nine-pair rerun

Balanced public concat now reuses its already checked length/weight values
when constructing the parent. The overflow order, profile guard, wide balance
arithmetic, eager/lazy hashes and node sequence are preserved. The shared
constructor tail is inlined so recursive joins do not acquire another call.
The unbalanced helper reloads weights, allowing the existing join edge to be a
tail jump on this Windows ABI. Exact comparison source is unchanged from round
five. No new cache, metadata field, fingerprint or default is introduced.

The [full validation](VALIDATION.md) passes 221 release tests, strict Clippy,
formatting, wasm compilation in both feature modes and warning-free docs.
Both new regression tests and existing published fingerprint/edit goldens pass.
Layout remains 48 bytes per node, 2,136 inline bytes per arena and 255 nodes
for the 64 KiB fixture; this is not a total-memory measurement.

The final primary cohort retains all 36 processes and 6,588 samples. Times below
are medians of process medians. Changes and intervals use the nine paired time
ratios, so they are not ratios of the displayed medians. The baseline here is
**round five**, not published 0.3.3. See [all 66 cells](final_comparison/all_cells.csv)
and the [verified summary](final_comparison/summary.json).

| Case | Round five | Round six | Paired change | Ratio interval | Faster pairs |
|---|---:|---:|---:|---:|---:|
| Balanced concat, two 64 KiB roots | 15.450 ns | 13.050 ns | −14.5% | [0.7674, 0.9526] | 8/9 |
| Independent 1M edit + exact compare | 136.658 µs | 137.900 µs | +1.4% | [0.9834, 1.0452] | 3/9 |
| Independent 100K edit + exact compare | 17.068 µs | 18.006 µs | +4.4% | [1.0015, 1.0879] | 1/9 |
| Independent 1 MiB exact LCP | 50.454 µs | 49.991 µs | −2.0% | [0.9354, 1.0270] | 6/9 |
| Independent 64 KiB exact LCP | 3.274 µs | 2.766 µs | −1.9% | [0.8578, 1.1220] | 4/9 |
| Self-concat 2 TiB DAG exact LCP | 141.000 ns | 114.000 ns | −7.2% | [0.8204, 1.0507] | 8/9 |
| Build 64 KiB | 27.220 µs | 39.910 µs | +7.2% | [0.8372, 1.3734] | 4/9 |

Concat improves in this cohort, while the original long edit/comparison penalty
does not recur at its previous magnitude. The previous LCP/DAG optimizations
remain in place; the follow-up LCP estimates above do not establish equivalence
or a further universal speedup. The 100K secondary result is a concern: it is
slower in eight pairs and its unadjusted interval barely excludes one. Build
timing is also mixed and highly variable. All outcomes are retained.

### Separate retained-binary confirmation

After observing the 100K result, a [fixed six-pair protocol](final_confirmation/PLAN.md)
was declared before further timing. It uses the exact same four final binaries,
without rebuilding, and reruns every cell: 24 processes and 4,392 samples. Its
[66-cell summary](final_confirmation/summary.json) and
[table](final_confirmation/all_cells.csv) remain separate from the primary nine;
they are not pooled or substituted for the earlier observations.
The [confirmation report](final_confirmation/README.md) includes all 66 cells,
direction counts and execution-order breakdowns.

| Case | Primary nine-pair ratio | Confirmation six-pair ratio | Confirmation interval | Confirmation slower pairs |
|---|---:|---:|---:|---:|
| Balanced concat | 0.8550 | 0.8777 | [0.7252, 1.0624] | 1/6 |
| Independent 100K edit/comparison | 1.0438 | 0.9834 | [0.9228, 1.0479] | 4/6 |
| Independent 1M edit/comparison | 1.0138 | 0.9604 | [0.9070, 1.0169] | 2/6 |
| Independent 1 MiB exact LCP | 0.9801 | 1.0680 | [0.9921, 1.1496] | 5/6 |
| Independent 64 KiB exact LCP | 0.9810 | 1.0776 | [0.9364, 1.2400] | 3/6 |
| Late-mismatch 1 MiB exact LCP | 0.9825 | 1.0410 | [0.9597, 1.1292] | 4/6 |
| Self-concat 2 TiB DAG exact LCP | 0.9284 | 1.0131 | [0.8814, 1.1645] | 2/6 |
| Build 64 KiB | 1.0723 | 0.8958 | [0.6911, 1.1612] | 3/6 |

Concat's point estimate improves by 14.5% and 12.2% in the two cohorts,
supporting retention of the narrow optimization. The confirmation's interval
still includes one. The 100K result does not reproduce as a paired geometric
mean, but four of six confirmation pairs remain slower and the median of
process medians increases from 18.113 to 18.406 µs. Its initial eight-of-nine
direction is not erased: the cost is uncertain/context-sensitive, not proven
noise. These outcomes do not support a universal no-regression claim.

Secondary LCP estimates also reverse direction: independent 1 MiB LCP is
6.8% slower in confirmation (five of six pairs), despite a 2.0% improvement
in the primary cohort. Independent 64 KiB and late-mismatch 1 MiB LCP estimates
are +7.8% and +4.1%. No confirmation slowdown has an unadjusted interval wholly
above one, but that does not make these costs zero or prove all prior LCP gains
survive every context. They remain part of the tradeoff documented here.

An inspection of the actual final churn objects finds identical normalized
instructions for the dispatcher, recursive comparator and frontier expansion
between these final arms; see [the disassembly record](codegen/final_churn_comparison/manifest.json).
That rules out changed instruction bodies in those routines, not linked
placement, data addresses or the effects of changed edit paths. No causal
explanation for the mixed 100K result is established, and no further comparator
or default change was selected.

### Remaining Ropey gaps in the primary cohort

Ropey still wins every independent-target size and both smaller shared-history
cases. These are candidate-cohort µs per cycle; paired Hashrope/Ropey ratios use
the nine corresponding process medians rather than the displayed median ratio:

| Bytes / target history | Round-six Hashrope | Ropey | Paired Hashrope/Ropey ratio |
|---|---:|---:|---:|
| 10,000 / independent | 4.608 | 1.068 | 4.231 |
| 10,000 / shared | 3.586 | 0.950 | 3.890 |
| 100,000 / independent | 18.006 | 9.250 | 1.957 |
| 100,000 / shared | 7.100 | 4.798 | 1.515 |
| 1,000,000 / independent | 137.900 | 98.382 | 1.385 |
| 1,000,000 / shared | 8.214 | 32.498 | 0.241 |

## Protocol and reproducibility

The [diagnostic plan](PLAN.md) predates diagnostic builds and measurements.
The longer diagnostic uses swapped implementation slots, an identical-source
A/A control, separate pinned/unpinned arms, five cases and 36 fresh processes.
It separates combined, edit-only and comparison-only histories; their timings
are not additive. Every concat session retains its full seven-sample trajectory.
The process-level estimator is distinct from accumulated timed mean and wall
time, which are also retained. Inner sessions are not independent replicates.

The [final comparison protocol](final_comparison/PLAN.md) compares the frozen
round-five core with the selected concat optimization using the unchanged
42-cell microbenchmark and 24-cell edit/comparison workload, including Ropey
1.6.1 and historical hash-only controls. Nine pairs per workload give 36 fresh
processes and 6,588 samples. Hash-only controls are not correctness substitutes.

Measurements have an exclusive agent CPU slot, but CPU frequency, temperature
and unrelated OS activity remain uncontrolled. No outer Python interleaving,
sample exclusions, replacement processes or adaptive early stopping are used.
Exploratory intervals are not adjusted for multiple cells and do not establish
equivalence or universal bounds. Repeated seeds represent repeated workloads.
All source snapshots, locks, binaries, hashes, raw output, orders and analysis
scripts are retained. This Windows/default-feature study does not benchmark
other platforms, languages or the bigint backend.
