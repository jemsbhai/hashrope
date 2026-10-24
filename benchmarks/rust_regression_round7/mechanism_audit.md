# Existing-data audit for the two round-seven targets

The round-six results leave two different unresolved patterns. The 100K
independent exact edit/comparison slowdown is strong in the original nine pairs
but changes direction under the confirmation's paired estimator. The 1 MiB
independent-equal LCP slowdown appears in the confirmation, in both execution
orders, after a modest improvement in the first cohort. Neither pattern can be
discarded from source equality, nearby control variation, or intervals crossing
one. Neither identifies a runtime cause by itself.

This audit performs no builds, new timings or production edits.
[existing_data_audit.py](existing_data_audit.py) and its
[JSON output](existing_data_audit.json) independently reparse all 60 processes,
10,980 raw samples and 132 paired statistic vectors from the two round-six
cohorts. Input/source/archive hashes and the four retained executable identities
are bound. Both cohorts used the same four executable bytes; there was no rebuild
between them. All observations remain; no sensitivity subset replaces a primary
estimate. The newer round-seven recurrence and future focused diagnostic are
separate studies, not pooled into this existing-data analysis.

## Target estimates and estimator differences

Times are medians of process medians in microseconds. Ratios and exploratory
95% log-t intervals use paired process medians, not the displayed time ratio.

| Target | Cohort | Baseline | Candidate | Paired ratio [interval] | Slower pairs |
|---|---|---:|---:|---:|---:|
| 100K independent exact edit/comparison | Initial nine | 17.068 | 18.006 | 1.0438 [1.0015, 1.0879] | 8/9 |
| 100K independent exact edit/comparison | Separate six | 18.113 | 18.406 | 0.9834 [0.9228, 1.0479] | 4/6 |
| 1 MiB independent-equal exact LCP | Initial nine | 50.454 | 49.991 | 0.9801 [0.9354, 1.0270] | 3/9 |
| 1 MiB independent-equal exact LCP | Separate six | 53.279 | 55.6685 | 1.0680 [0.9921, 1.1496] | 5/6 |

The 100K confirmation ratios, in process-pair order, are **1.0465, 1.0109,
0.9258, 1.0105, 0.8983, 1.0171**. The two faster pairs have larger proportional
changes than the four slower ones. Their contribution makes the paired geometric
mean fall below one while the ratio of displayed arm medians is 1.0162. This is
an estimator difference, not an error or a reason to choose whichever result
looks preferable. The initial cohort's ratio of displayed medians is 1.0550,
also distinct from its paired 1.0438 estimate.

For 1 MiB LCP, confirmation ratios are **1.0061, 0.9823, 1.1004, 1.1613,
1.1434, 1.0275**. Its displayed-median ratio is 1.0448, versus paired 1.0680.
The point-estimate increase is not confined to one process. Conversely, the
initial cohort includes both increases and decreases, with a paired range of
0.8466–1.0407. These small cohorts and repeated workloads do not establish a
universal distribution or an equivalence bound.

## Execution order, seed and sample position

For 100K, baseline-first/candidate-first geometric means are **1.0477/1.0390**
in the initial cohort and **0.9548/1.0128** in confirmation. Initial seed strata
1/2/3 are **1.0448/1.0229/1.0640**; confirmation strata are
**1.0283/0.9530/0.9704**. Thus a simple uniform first/second-process effect does
not account for the original slowdown. The confirmation is order-sensitive,
but that association does not establish the cause.

Seed is also coupled to the original within-process arm order. With invocation
zero, the exact Hashrope arm is fourth, third or second for seeds 1, 2 or 3.
Its raw 100K independent samples are lines **34–36, 31–33 or 28–30** respectively;
all 10K cases precede them. A seed stratum therefore changes both the edited
bytes and preceding allocations/queries. It is not an isolated input-effect
estimate. Each of its three timed samples builds fresh source/target arenas,
prebuilds 50 independent targets, warms 50 edits, then times another 50 edits
and exact comparisons. Those 50 outcomes alternate equality and a mismatch at
the final byte. Setup, byte validation and destruction are outside its timer,
but still contribute to execution history.

The micro LCP target is always the **25th timing cell**, raw lines **172–178**.
It uses 100 queries per sample after three untimed queries, seven samples total.
The same-root LCP/equality/hash measurements precede it; independent equality
follows it immediately. The arena already contains first, independently equal,
early-mismatch, late-mismatch and shared-late-mismatch roots. Its independent
trees are in **one arena**, using the same-arena recursive exact traversal. The
100K edit workload compares across arenas and uses the long frontier path;
the two target costs cannot be interchanged or added together.

LCP baseline-first/candidate-first means are **0.9726/0.9896** initially and
**1.0818/1.0544** in confirmation. Both orders retain the confirmation slowdown.
The raw batches contain changes in both directions within a process. For example,
confirmation pair 2's candidate starts around 52–54 microseconds and ends around
61–67; pair 3 starts around 67 and falls toward 52–63. Pair 4's baseline has a
77-microsecond final batch after six near 50. There is no uniform
first-sample-only pattern supporting a retrospective warmup exclusion.

The JSON retains every sample vector, raw line position and process timestamp.
Its sample-position ratios and later/earlier summaries are descriptive checks
of correlated samples, not extra independent replicates or replacement
estimators. For example, candidate/baseline LCP ratios by sample position are
1.0590, 1.0387, 0.9757, 0.9711, 0.9783, 0.9828, 0.9768 initially, and
1.0875, 1.0249, 1.0589, 1.0963, 1.0536, 1.0804, 0.9651 in confirmation.

## Nearby controls and limits of attribution

| Nearby cell | Initial ratio | Confirmation ratio |
|---|---:|---:|
| 100K independent Hashrope fingerprint candidate | 0.9747 | 1.0347 |
| 100K independent Ropey exact | 0.9791 | 1.0616 |
| 100K independent Ropey rehash candidate | 1.0047 | 1.0255 |
| 100K shared-history Hashrope exact | 0.9939 | 0.9403 |
| 1M independent Hashrope exact edit/comparison | 1.0138 | 0.9604 |
| 1 MiB independent exact equality | 0.9486 | 1.0307 |
| 1 MiB late-mismatch exact LCP | 0.9825 | 1.0410 |
| 1 MiB late-mismatch exact equality | 1.0196 | 1.0408 |
| 1 MiB independent hash-only LCP | 0.9937 | 0.9818 |
| 64 KiB independent exact LCP | 0.9810 | 1.0776 |

Related exact comparisons share the confirmation's positive direction, while
hash-only LCP does not. The initial 100K slowdown is absent from some nearby
controls, whereas several confirmation controls become slower. This is useful
localization, not a causal adjustment: controls do different work, at different
positions, with different allocations. Dividing a target by a conveniently
selected control would introduce another estimator and would not identify the
cause. Shared-history and cross-arena traversal have different shortcuts.

The retained final churn objects previously had identical normalized comparator
instruction bodies across arms. That does not cover the actual micro call site,
linked addresses, caller/inlining context, allocator addresses, preceding edit
paths or machine state. The four binaries are unchanged across these two
cohorts, so a rebuild is not an explanation for the phase reversal. Static
placement differences between the two arms can still interact with execution
context. Frequency, temperature and interference were not measured; they are
possible confounders, not established explanations.

## Bounded follow-up design and mechanism checks

The root-owned [retained-executable recurrence plan](original_replication/PLAN.md)
fixes 12 new pairs per workload, using the exact same four binaries, full
sequences and all cells. It balances six pairs per execution order, gives each
churn seed both orders twice, and keeps the cohort separate. This preflight
passed review. It can establish recurrence under another full-protocol sample;
it cannot by itself isolate an allocation or instruction mechanism.

A separately frozen focused diagnostic can compare the two archived cores in
swapped crate slots and current/current A/A, with pinned and unpinned modes
reported separately. The proposed five cases are original 100K combined
edit/comparison, edit-only, comparison-only over prebuilt 50 edited states,
original 1 MiB same-arena LCP, and original-growing concat as a control. Retain
the full original LCP arena, describe the replacement of its preceding timed
history by focused warmup, and keep source/target construction outside timers.
Case order and ABBA/BAAB order must be independent; both slots must use one
common calibrated count of bounded sessions. Retain all raw samples, fixed
process counts, original estimators and slower outcomes. Separate setup
histories make the edit-only and comparison-only times nonadditive.

Before any runtime change, useful untimed probes are counts of balanced versus
unbalanced concat/join entries, metadata checks, internal allocations/rotations,
and, for exact queries, recursive/frontier visits, leaf comparisons and spans,
and fallback events. Compare identical seeded edit/node traces and the full
LCP fixture across arms. Equal counts exclude added instances of that measured
work, not different cost per instance. A proposed cursor or constructor change
must preserve range bounds, earliest mismatch, visit-budget fallback, profile
checks, overflow order, node allocation sequence and fingerprints.

Actual retained micro and churn disassembly should distinguish changed
instruction bodies from placement, alignment, call targets and unwind/stack
layout. A controlled placement or factoring ablation must be isolated, preserve
the production source and original evidence, and have a declared measurement
protocol before its timings are observed. No optimization is justified solely
by a favorable subset or by asserting that unchanged source proves noise.

## Review of the bounded direction-flag candidate

The proposed runtime diff against the frozen round-six source only passes
`wl > wr` from the two existing callers into `join_unbalanced` as a boolean,
instead of looking up both weights again inside that helper. Both callers first
apply the unchanged balance test. No mutation intervenes between reading the
weights and calculating the flag, so the flag is the same predicate previously
evaluated inside the helper, including its equality case. The recursive
Internal/Repeat/Leaf branch bodies are unchanged.

The public profile guard, empty handling, checked length-then-weight preflight,
wide balance arithmetic, hash operands, lazy sentinel and allocation sequence
are untouched. This source review found no correctness blocker. It establishes
removal of two source-level weight lookups, not a measured speedup. The inspected
Windows build retains the tail jump with four register arguments, but its
wrapper grows from 82 to 91 static instructions and gains register-save work;
the final complete cohort must assess that tradeoff. Full validation and the
predeclared 12-pair final comparison remain separate from this read-only review.
