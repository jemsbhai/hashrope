# Round six: existing evidence and causal questions

The remaining slowdowns merit investigation. The round-five nine-pair estimate
is **+8.4% for balanced concat** and **+5.1% for independent 1M edit plus exact
equality**. Neither unchanged source nor variable controls proves that these
effects are noise. This document separates observations, source-derived
hypotheses and tests needed before selecting another production change.

## Independent audit of the retained data

[existing_data_audit.py](existing_data_audit.py) independently reads all 36 raw
process records and 6,588 timing samples, checks their recorded hashes and four
archived executables, and reconstructs all 132 supplemental-six/pooled-nine
statistic vectors. All checks passed. Its output retains every selected cell's
raw samples, process medians, paired ratios and phase/order/seed groups in
[existing_data_audit.json](existing_data_audit.json). No observation was excluded
and no new measurement was made by this analysis.

All ratios below are candidate/released-0.3.3 geometric means of paired process
medians. These post hoc strata are descriptive and have only three to six
pairs; they are not significance tests or independent causal experiments.

| Cell / grouping | Original three | Supplemental six | All nine | Slower pairs |
|---|---:|---:|---:|---:|
| Balanced concat | 1.0785 | 1.0865 | 1.0839 | 6/9 |
| 1M independent Hashrope exact churn | 1.1151 | 1.0204 | 1.0511 | 5/9 |
| 1M independent Hashrope fingerprint-only churn | 1.0330 | 0.9729 | 0.9925 | 3/9 |
| 1M independent Ropey exact churn | 1.2986 | 1.0115 | 1.0993 | 5/9 |
| 1M independent Ropey rehash churn | 1.0186 | 0.9873 | 0.9976 | 5/9 |

Concat is slower in both phases by a similar amount. Baseline-first pairs give
1.0859 and candidate-first pairs 1.0813; a simple first-versus-second process
explanation is unsupported. For 1M exact churn those ratios are 1.0473 and
1.0559. Its effect is smaller in the supplemental phase, but it remains positive
in both. The nine exact-churn ratios span 0.933 to 1.307, so the pooled mean must
not be presented as a uniform cost per run.

The 1M exact-churn seed groups give ratios 0.9868, 1.0585 and 1.1117 for seeds
1, 2 and 3. Each group has three repeated processes of the same edit fixture.
Crucially, the harness also rotates comparison-arm order by seed: Hashrope
exact is fourth, third and second respectively. Thus **seed and arm position
are coupled**. These data cannot distinguish input/tree-shape effects from
preceding-work effects. The microbenchmark has no random seed; pair number
modulo three is not treated as a seed here.

The fingerprint-only and Ropey arms expose context variation, but they do not
run the same operation in the same timing block. Their values cannot be
subtracted to isolate exact-comparison cost, and they do not cancel the
Hashrope result. Phase, process order and seed are retained rather than selected
according to which grouping gives the preferred result.

## What balanced concat actually measures

The original microbenchmark constructs two independent 64KiB ropes, each with
255 nodes. It then repeatedly calls `concat(left, right)` with those **same
fixed inputs**: three untimed warmups followed by seven samples of 2,000 joins.
The returned node is black-boxed, but each result remains in the arena. This is
not an ever-growing logical string: each output has 128KiB, 256 leaf weight and
the same content fingerprint. It is an ever-growing **node arena**, receiving
14,003 Internal nodes. Its node count goes from 510 to 14,513.

Each balanced join follows the public profile/overflow checks, the weight
balance check, a new Internal node's metadata/hash calculation, and `Vec::push`.
The power for the fixed right length is warmed before timing. No comparison
engine runs in this cell. The candidate leaves this source algorithm unchanged,
but that does not establish identical linked instructions or execution context.

The raw data have a strong sample-position pattern: positions 0, 1 and 3 are
usually much more expensive than 2, 4, 5 and 6. Under the expected `Vec` doubling
trajectory, the arena reaches length 513/capacity 1,024 after warmup, then grows
at capacities 1,024 and 2,048 in sample 0, at 4,096 in sample 1 and at 8,192 in
sample 3. This predicted reallocation counts `[2,1,0,1,0,0,0]` before the probe.
The subsequent untimed capacity/allocator replay confirmed this exact trajectory
in both immutable cores. Allocator relocation and page behavior still cannot be
inferred from the elapsed samples alone.

When those three growth-heavy samples exceed the other four, the predeclared
median of seven is the **largest of the four growth-free samples**. That helps
explain why a single high growth-free observation can move the reported process
median, but is not evidence that the candidate gap is spurious. The original
estimator and every raw sample must remain intact.

Useful diagnostic controls are the full original seven-sample trajectory and a
separately labeled no-growth setup. A proposed public-API setup performs 50,003
untimed joins before the same seven-by-2,000 loop: expected lengths 50,513 to
64,513 fit within capacity 65,536. The counter-only copies verified those
capacities and zero allocator events in every measured-shaped block. This setup
also changes retained memory, cache and operation history; it is a diagnostic control, not a clean
single-variable proof of allocation causality. Never run one arena until a
long timer expires: that changes the growth workload and retained memory.

The separate [capacity probe](concat_capacity_probe.json) recorded all 32 block
replays (two setups, warmup plus seven blocks, two core versions), with identical
node/capacity trajectories and successful content/metadata oracles. Every join
added one node; every resulting root had length 131,072, weight 256 and default
fingerprint `2138634010386205116`, checked against direct byte hashing outside
the counted interval. The growing case ended at 14,513 nodes/capacity 16,384.
The prewarmed case started at 50,513/capacity 65,536 and ended at 64,513 with the
same capacity. Its seven 2,000-join blocks requested no allocations,
reallocations or deallocations. The counted warmup still made seven node-vector
reallocations, which were outside those blocks. Capacity 65,536 represents
3,145,728 bytes of node-vector storage, not total arena/process memory.

[Probe source](concat_capacity_probe.rs), [driver](concat_capacity_probe.py),
[manifest](concat_capacity_probe_manifest.json), [build record](concat_capacity_probe_build.json),
raw TSVs and [archived copied sources](concat_capacity_probe_sources.zip) retain
the proof. Only a capacity accessor was added to the copied cores; a separate
counting allocator recorded requested traffic. No production or timed-harness
source was instrumented. These counts confirm fixture mechanics, not the cause
or size of a latency difference.

## What independent 1M churn actually measures

Each sample constructs a fresh source arena and a distinct target arena, with
50 independent 1M target ropes. Target bytes and target construction are outside
the timer. Fifty untimed source edits warm representative powers and leave
their nodes in the source arena. The timed section repeats the 50 seeded
delete/insert edits, each followed by equality against its corresponding
prebuilt target. Exactly half are equal and half differ at the final byte;
answer masks and final bytes are checked outside timing.

Every compared range is long and crosses arenas. The dispatcher therefore
uses the fixed-stack **cross-arena frontier**, rather than the round-five
aligned same-arena recursive shortcut. Inputs are ordinary edited trees with
bounded leaves, so the source-level expected path is frontier descent and leaf
comparison, with no periodic witness. That expectation should be verified by
event counters: a hidden fallback must not be assumed absent merely from the
fixture description.

The current frontier reexamines both current nodes each iteration. When one
side has reached a leaf while the other side descends or finishes a shorter
leaf, the stationary leaf can be inspected again. Caching a leaf cursor could
remove source work, but its effect on branching, stack usage and code layout
needs measurement. Counters should distinguish total frontier iterations,
Internal expansions, stationary-leaf observations, leaf comparison calls and
input spans, plus budget/Repeat/depth fallback. Separate query-only and edit-only
sessions are needed to localize cost; neither can simply be added or subtracted
from combined timing because cache and allocation histories differ.

## Bounded diagnostic strategy

The benchmark owner will freeze the actual protocol before timing. The useful
controls, in priority order, are:

1. Retain fixed-count full-sequence process pairs of the original archived
   binaries. This asks whether the same measured executable gap recurs without
   changing code layout or replacing the original workload.
2. Compare released/current cores in a common focused harness, reverse their
   crate slots, and include current/current A/A. Keep the original bounded concat
   trajectory and 50-edit sessions, with setup/drop outside the timer. Longer
   accumulated timing must repeat bounded sessions, not expand one arena.
3. Measure combined edit-plus-equality, edit-only and prebuilt equality-only
   sessions separately. Counterbalance case/arm order independently of seed;
   ABBA/BAAB block order must also be independent of that order. Keep every
   block, using a predeclared process-level estimator.
4. Keep the ordinary unpinned mode and any pinned mode separate. Processor IDs
   at block endpoints and affinity settings describe the run but do not measure
   every migration, frequency change or sibling-thread interference.

The earlier round-four controlled study supplies useful methods, not a verdict
for this new question. It exposed slot/affinity sensitivity and used A/A plus
retained-original-binary checks; the original full protocol later retained a
regression. A focused harness changes linking and preceding work even when its
library source is identical. Failure to reproduce a gap there cannot erase an
effect in the retained original executable cohort.

## Code-generation and correctness review before an optimization

The implementation agent's [retained code-generation evidence](codegen/README.md)
examines construction and frontier/call-site assembly. Normalized COFF
instructions for `concat`, `join`, `make_internal`,
`split_inner` and `ExactFrontier::expand` match across the round-five executables'
retained library objects. The balanced path nevertheless repeats checked
length/weight work through the three construction routines. The inlined frontier
inside `exact_lcp` changes from 490 to 507 disassembled instruction lines, with
different register/stack handling around the leaf-comparison call. That is
generated-code evidence supporting investigation of an outlining ablation;
it is **not** a measured attribution of the 5.1% churn gap. Uniformly normalized
function text is weaker than final instruction/address identity, and matching
function bodies alone do not bind their callees, allocator state or placement.
No cause is assigned solely from source equality.

A prechecked balanced-concat constructor could avoid repeated length/weight
loads, if assembly confirms them. It must preserve the public profile guard,
checked length then weight overflow behavior, failure before allocating output
nodes, exact polynomial arithmetic, and the balanced/unbalanced structure. The
unbalanced path, custom profiles and boundary-sized lengths need independent
regression coverage; unchecked arithmetic is not an acceptable shortcut.

A cached frontier must retain byte order, exact first-mismatch location,
cross-arena identity rules, a work bound in stored-node terms, and Repeat/depth
fallback at the first unverified byte. A smaller event count is not sufficient
if it silently bypasses the compact-DAG bound. Its stack and allocation costs
must also remain bounded. These are review criteria, not a recommendation to
ship either design before evidence.

No production source is changed by this audit, and it performs no new timing.
