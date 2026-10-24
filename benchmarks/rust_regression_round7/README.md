# Rust 100 KB edit/comparison and 1 MiB LCP investigation

**Retain the round-six implementation.** The earlier 4.4% 100 KB combined
slowdown and 6.8% 1 MiB LCP slowdown do not repeat at stable magnitudes. However,
controlled pinned LCP still shows about a 2% cost in both implementation-slot
orders, and a secondary comparison-only estimate is slower. These costs are
context-sensitive and unresolved, not established harmless noise.

A narrow join optimization was implemented, validated and measured against
round six, then **rejected**. Neither primary target showed an established
benefit, and some secondary early-mismatch comparisons slowed. All candidate
sources, binaries, results and the [decision](decision.json) remain available.
Production is byte-identical to the six frozen round-six core files; see
[retained_source_identity.json](retained_source_identity.json). Existing exact
APIs, collision handling, fingerprints, defaults and previous optimizations
remain intact. All 221 release tests and tooling/platform compilation checks
pass again on the restored core. No package version or publication changes are
made.

## What repeats with the exact original executables

The [previous-data audit](mechanism_audit.md) rechecks both round-six cohorts,
including the estimator differences and full sample histories. The new
[original-executable protocol](original_replication/PLAN.md) uses the same four
retained binaries, without rebuilding. Twelve balanced process pairs per
workload retain all 42 microbenchmark and 24 edit/comparison cells: **48
processes and 8,784 samples**. The twelve pairs are a separate cohort; they do
not replace or pool with the earlier nine or six.

Ratios below are round six / round five. Intervals are exploratory 95% log-t
intervals over paired process medians, not equivalence bounds. Displayed times
are medians of process medians, so their ratio need not equal the paired ratio.

| Case | Earlier nine-pair ratio | Earlier six-pair ratio | New twelve-pair ratio [interval] | New slower pairs |
|---|---:|---:|---:|---:|
| 100,000-byte independent edit + exact compare | 1.0438 | 0.9834 | 1.0200 [0.9661, 1.0769] | 8/12 |
| 1,048,576-byte independent-equal exact LCP | 0.9801 | 1.0680 | 0.9933 [0.9387, 1.0511] | 7/12 |
| Balanced concat, two 64 KiB roots | 0.8550 | 0.8777 | 0.8476 [0.7657, 0.9383] | 1/12 |

The new 100 KB time medians are 20.368→20.830 µs; the paired estimate is +2.0%.
Its eight slower pairs and earlier eight-of-nine direction still matter even
though the interval includes one. New LCP medians are 55.8565→54.4075 µs; the
paired estimate is −0.7%, with mixed individual directions. Concat's earlier
gain repeats at −15.2%, with eleven of twelve pairs faster.

All cells, individual ratios, order/seed strata and raw samples are retained in
[summary.json](original_replication/summary.json),
[all_cells.csv](original_replication/all_cells.csv) and
[results.json](original_replication/results.json). An independent
[replication audit](original_replication_review.json) verifies all 48 processes,
8,784 samples and 66 summaries. Neither source equality nor a confidence
interval crossing one proves absence of a smaller performance cost.

## Separate phase, placement and affinity diagnostics

The [controlled report](controlled_report.md) preserves all 30 primary cells
and 90 estimator vectors from **36 processes, 720 blocks and 153,552 samples**.
It compares frozen round five and round six, before the proposed new change,
with both implementation-slot assignments and an identical-source A/A control.
Pinned and unpinned processes remain separate. Five cases cover 100 KB combined,
edit-only and comparison-only histories, original-context 1 MiB same-arena LCP,
and the original growing-arena concat control. The three edit/comparison phase
times have different setup histories and are not additive.

The combined 100 KB primary ratios range from 0.9787 to 1.0145 across the four
comparison conditions; all four intervals include one. There is no uniform
4.4% combined penalty. Nevertheless, pinned swapped-slot comparison-only
**timed mean** is 1.0516 [1.0182, 1.0862], with six of six slower. Its primary
median estimate is 1.0198 [0.9936, 1.0466]; normal-slot estimates differ.
This secondary signal remains evidence, not an interchangeable estimate of
the combined cost.

Pinned independent 1 MiB LCP ratios are **1.0203 [1.0032, 1.0376]** in normal
slots and **1.0200 [0.9957, 1.0448]** after swapping, with five and four of six
processes slower. Normal-slot timed mean and wall estimates also have positive
intervals excluding one. Unpinned primary estimates, 0.9778 and 1.0266, vary
more. The pinned A/A estimate, 1.0114 [0.9991, 1.0238], cannot be subtracted
from the measured cost or used to declare it noise. The smaller pinned penalty
is unresolved; no thermal, allocation or placement cause has been isolated.

Concat improves in every comparison process in all four conditions, with
ratios 0.8450–0.8642 and all four primary intervals below one. This supports
the existing round-six concat optimization for this fixture. It does not
resolve the comparison costs or establish universal performance gains.

The [independent diagnostic audit](controlled_review.md) verifies provenance,
all samples, the full estimator vectors and rendered report. Complete protocol
differences and limitations, including extra byte validation and omission of
earlier full-sequence timed cells, are stated in the controlled report.

## Proposed optimization and why it was dropped

Generated-code inspection found real redundant work in the unbalanced join:
both callers already compute node weights and balance, but the helper reloads
both weights to decide which side is heavier. A trial passes that known boolean
in the fourth Windows x64 register argument. It preserves the tail jump and
removes those two lookups without changing the predicate or recursive body.

The tradeoff is visible in [the generated-code record](codegen/README.md).
The helper shrinks from 247 to 241 normalized instruction lines, but `join`
grows from 82 to 91 and saves/restores two extra registers; public concat grows
from 226 to 235. Fewer lookups alone do not establish lower runtime. Exact
comparator instruction bodies are unchanged. In the original round-six
binaries, the dispatcher moves by `0x150`, while the recursive comparator and
frontier expansion remain at the same linked locations. These observations do
not prove a placement explanation for the measured regression.

The trial passed all **221 release tests**, lint/format/docs and both wasm
compilation checks. A [fixed final protocol](final_comparison/PLAN.md) then
compared frozen round six with the trial using the unchanged full workloads:
**12 pairs per workload, 48 fresh processes, 8,784 samples**. All four binaries
were built before timing. No sample, failed process, or unfavorable cell was
removed or replaced. This is a new optimization comparison, separate from the
round-five/round-six diagnostics above.

### Final trial results

Baseline is the **retained round-six core**; candidate is the **rejected
direction-argument trial**. Time changes use paired geometric means. All
intervals are exploratory and unadjusted across the 66 cells.

| Case | Retained baseline | Rejected candidate | Paired time change [ratio interval] | Candidate slower |
|---|---:|---:|---:|---:|
| **100 KB independent edit + exact compare** | 20.197 µs | 20.163 µs | **−1.0% [0.9479, 1.0336]** | 5/12 |
| **1 MiB independent-equal exact LCP** | 54.725 µs | 53.871 µs | **−0.4% [0.9602, 1.0340]** | 6/12 |
| 64 KiB independent-equal exact LCP | 3.027 µs | 2.987 µs | −5.1% [0.8858, 1.0170] | 6/12 |
| 1 MiB late-mismatch exact LCP | 54.679 µs | 53.596 µs | −2.9% [0.9446, 0.9976] | 1/12 |
| 2 TiB self-concat DAG exact LCP | 124 ns | 124 ns | −0.4% [0.9046, 1.0970] | 5/12 |
| Balanced concat, two 64 KiB roots | 13.875 ns | 14.250 ns | −2.9% [0.8037, 1.1734] | 7/12 |
| 64 KiB early-mismatch exact LCP | 38.025 ns | 40.750 ns | +4.9% [1.0121, 1.0874] | 10/12 |
| 1 MiB early-mismatch exact equality | 58.225 ns | 59.550 ns | +3.9% [1.0054, 1.0745] | 9/12 |
| 100 KB shared-history edit + exact compare | 7.549 µs | 7.947 µs | +4.2% [0.9873, 1.0995] | 7/12 |
| 1 MB independent edit + exact compare | 153.603 µs | 152.171 µs | +0.1% [0.9556, 1.0495] | 6/12 |
| 1 MB shared-history edit + exact compare | 9.125 µs | 8.635 µs | −8.5% [0.8091, 1.0347] | 3/12 |
| Split chunked 64 KiB | 444.5 ns | 454.0 ns | +4.9% [0.9813, 1.1218] | 9/12 |

Neither primary target establishes a useful improvement. The 100 KB estimate
changes direction by execution order: baseline-first 1.0082, candidate-first
0.9717; its seed strata are 0.9645, 0.9759 and 1.0302. LCP order estimates are
0.9846 and 1.0084. The concat paired estimate is negative despite seven slower
pairs and a higher displayed candidate median; its ratio range is 0.4927–1.3460
and its order estimates are 0.8137/1.1590. These are mixed data, not a dependable
new concat improvement.

The favorable late-mismatch result and 1 MB fingerprint-candidate result
(−7.2%, interval [0.8638, 0.9959]) are retained too. Six secondary slowdown
intervals exclude one: the two exact early-mismatch cases above, three hash-only
LCP controls (1 MiB independent and same-root, 64 KiB late mismatch), and the
10 KB shared-history fingerprint-candidate workload (+5.7%). Hash-only
controls are not substitutes for correctness-preserving exact APIs. Those
multiple exploratory signals do not prove a universal regression, but there
is insufficient benefit to accept the register/code-size tradeoff. No further
timing was used to seek a favorable subset.

The [complete final table](final_comparison/README.md),
[summary](final_comparison/summary.json), [CSV](final_comparison/all_cells.csv)
and [raw records](final_comparison/results.json) preserve all 66 cells and
orders. Memory layout remains 48-byte `NodeInner`, 2,136-byte inline `Arena`,
and 255 nodes for 64 KiB in both arms. The summary's working-source-match flags
describe the candidate at measurement time; [decision.json](decision.json)
records the subsequent rejection and exact baseline restoration.

### Remaining Ropey gaps in the retained implementation

These measurements use the **baseline arm** of the final comparison, because
that is the implementation retained in production. Times are µs per edit and
comparison cycle. Ratios pair Hashrope/Ropey medians within each process; they
need not equal the displayed time ratios. Ropey is the unchanged locked 1.6.1
control, and these ASCII fixtures do not compare Unicode editing APIs.

| Bytes / target history | Retained Hashrope | Ropey | Paired Hashrope/Ropey time |
|---|---:|---:|---:|
| 10,000 / independent | 4.725 | 1.233 | 3.549× |
| 10,000 / shared | 3.862 | 0.981 | 3.822× |
| 100,000 / independent | 20.197 | 11.506 | 1.770× |
| 100,000 / shared | 7.549 | 5.367 | 1.357× |
| 1,000,000 / independent | 153.603 | 118.153 | 1.312× |
| 1,000,000 / shared | 9.125 | 37.526 | 0.246× |

Behavioral correctness and fingerprint compatibility do not imply efficiency
parity with Ropey. Independently equal content can still require scanning its
bytes; exact comparison has no general logarithmic bound.

## Reproducibility and validation

Frozen original sources/executables, candidate sources, locked dependencies,
compiler details, build logs, raw timing records and analyses are retained.
The final [build record](final_comparison/build.json) identifies Rust 1.94.0,
LLVM 21.1.8, `x86_64-pc-windows-msvc`, Windows 11 and an Intel family 6/model
183 processor with 32 logical processors. No compiler override flags were set.
The original and final twelve-pair cohorts each alternate process order six
times per arm; seeds 1/2/3 each receive both execution orders twice. The
controlled experiment keeps affinity, slot and A/A cells separate. Process
pairs, not thousands of inner samples, are the statistical units. Repeated
seeds are repeated workloads, not independent draws of new content.

Agent builds, tests and timings did not overlap. Frequency, temperature and
external OS activity remain uncontrolled. All intervals are exploratory,
unadjusted, and not equivalence tests. This Windows/default-feature study does
not benchmark bigint, other languages or other platforms.

[Validation](VALIDATION.md) distinguishes the tested rejected trial from the
fresh checks of the restored implementation. The trial remains reproducible
from its archived sources even though it is not selected. The previous
[round-six](../rust_regression_round6/README.md) and
[round-five](../rust_lcp_round5/README.md) studies remain historical evidence;
no earlier claim is silently replaced or promoted to a universal guarantee.
