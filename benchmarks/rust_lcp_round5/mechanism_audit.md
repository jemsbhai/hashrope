# Rust LCP mechanism audit

Two guarded changes remove measurable structural work while preserving the
byte proofs. They are plausible optimizations of the current implementation;
they do **not** establish why release integration slowed some earlier cells.
The uninstrumented paired benchmark must establish their actual timing effect.

## What the historical comparison can establish

The exact engine from completed round four and committed 0.3.3 is identical
after formatting both sections with the same `rustfmt`. This includes the
dispatcher, cache-free recursive/frontier traversal, mismatch localization,
structural-period discovery and memoized comparison. The formatted sections
both have SHA256 `9fc58bbce89ae77dd8ecb2163010d013be2dd4f1452005c75eca7a586def0441`.
[source_comparison_audit.json](source_comparison_audit.json) records both archive
hashes and the formatter; [audit_historical_source.py](audit_historical_source.py)
reproduces the check without changing either archive.

This is source evidence, not identical machine-code or final-link evidence.
The surrounding 0.3.3 integration changed arithmetic, overflow handling and
optional bigint support. Those changes can affect construction history,
inlining or linked code placement even when the comparison algorithm is
unchanged. This audit does not isolate any such cause. In particular, the
earlier three-pair release diagnostic had independent-equal LCP and DAG LCP
slower in all three pairs, while the corresponding equality cells behaved
differently. Neither shared traversal nor variable controls invalidate those
observations. The prior full-matrix protocol-sensitive recurrence also remains
unresolved; this new study uses a fresh concurrent 0.3.3 baseline.

## Instrumented source-event evidence

[source_work_probe.py](source_work_probe.py) instruments **copies** of both cores;
[source_work_probe.rs](source_work_probe.rs) constructs eight fixtures and makes
one oracle-checked exact range-LCP query per fixture and arm. Production source
is never instrumented. The baseline comes from the immutable Git-object
snapshot. The captured candidate `rope.rs` hash is
`8439fe38d04bde7a79ac5981ce3a4c6f5a5b4f727d9eeeaa278e2d286608b222`.
Both use the narrow backend and the same Rust compiler. The candidate manifest
differs from the Git snapshot only in CRLF/LF bytes; its contents are unchanged.

All sixteen oracle checks passed. Counts below come from actual execution of
the instrumented source, with construction outside the counted interval.
Atomic counters perturb generated code, so **no elapsed-time interpretation is
made**. These are source events rather than CPU instructions, cache misses or
physical memory loads. `leaf_compare_bytes` is the sum of lengths submitted to
the byte-comparison expression, not an assertion about how many bytes the
compiler's comparison routine reads. Allocator values count requested heap
traffic during the query, including its temporary cache's destruction; they do
not measure retained or peak process memory.

| Fixture | Source event | 0.3.3 baseline | Candidate |
|---|---|---:|---:|
| 64KiB aligned independent roots | Normalizer node visits | 510 | 256 |
| 64KiB aligned independent roots | Boundary calls | 254 | 0 |
| 1MiB aligned independent roots | Normalizer node visits | 8,190 | 4,096 |
| 1MiB aligned independent roots | Boundary calls | 4,094 | 0 |
| 1MiB aligned independent roots | Leaf comparisons / input bytes | 2,048 / 1,048,576 | 2,048 / 1,048,576 |
| 2TiB self-concat DAG | Period-chain node inspections | 82 | 82 |
| 2TiB self-concat DAG | Normalizer node visits | 84 | 2 |
| 2TiB self-concat DAG | Memo lookups / insert calls | 1 / 2 | 0 / 0 |
| 2TiB self-concat DAG | Allocations / requested bytes | 1 / 368 | 0 / 0 |
| Shifted 2TiB self-concat DAG | Normalizer node visits | 94 | 4 |
| Shifted 2TiB self-concat DAG | Allocations | 1 | 0 |
| Three-piece compressed unit | Normalizer node visits | 1,065 | 453 |
| Three-piece compressed unit | Memo lookups / allocations | 122 / 6 | 89 / 5 |
| Last-witness-byte mismatch in compressed unit | Normalizer node visits | 1,526 | 914 |
| Last-witness-byte mismatch in compressed unit | Allocations | 16 | 16 |

The two ordinary late-mismatch fixtures have the same traversal counts as
their equal counterparts, with answers 65,535 and 1,048,575. Both ordinary
versions allocate nothing during comparison. The two DAG fixtures still submit
only two bytes to leaf comparisons. The nonperiodic compressed unit is
`a^(2^32) b^(2^32)` (8GiB), paired with an independently built rotation by
17 bytes. Nonzero phases 18 and 1 require three valid pieces of one period.
The 16GiB query is answered without expansion; changing the other unit's byte
immediately before its start phase gives the expected mismatch at `2^33 - 1`.
That negative test still allocates sixteen cache nodes, so the optimization
does not make general compressed comparison allocation-free.

Full counts, answers and allocator requests are in
[source_work_probe.json](source_work_probe.json). The
[manifest](source_work_probe_manifest.json), [build record](source_work_probe_build.json),
raw per-arm TSVs and instrumented-source archive retain source, harness,
compiler and binary identities. Build logs are separate. A reproduction uses
the same instrumented source archive or the recorded candidate source with
`source_work_probe.py prepare`, `build`, then `run` in a fresh study directory;
the script refuses to overwrite its artifacts.

## Correctness and cost review of the aligned-child shortcut

The shortcut requires zero offsets, two Internal nodes and identical positive
left-child lengths strictly smaller than the remaining comparison. The shared
identity check still runs first. It proves the complete left prefix through
exact comparison, then moves to both right children. A mismatch or a bounded
fallback propagates its already verified prefix, including any earlier bytes
from the caller. Both parent visits debit the same visit budget used by the
general path, so aligned compact DAGs cannot silently bypass the bound.

For two perfect independent trees with `L` leaves, the baseline normalizes
`4L - 2` nodes and calls `boundary` `2(L - 1)` times. The candidate's successful
aligned path removes Internal normalization and boundary rediscovery, leaving
`2L` leaf-normalizer visits. It does not change the number of leaf comparisons
or their byte spans. The added branch can cost something for unaligned,
early-mismatch or shared cases; only the complete benchmark can assess that
tradeoff. This is a constant-factor traversal change, not a new asymptotic
bound or a claim that every independent tree has aligned boundaries.

## Correctness and cost review of period-root comparison

Period discovery now returns the underlying node, following only Repeat edges
and Internal nodes whose two child IDs are identical. The final unit is a Leaf
or an ordinary Internal; no hash or inferred byte periodicity is involved.
For equal periods `p`, one full period proves equality of the two phase-shifted
streams. Their offsets are reduced modulo `p`, and each proof piece ends before
either unit boundary. Across `p` bytes each phase wraps at most once, producing
at most three pieces. A differing byte returns its exact position; skipping the
rest is permitted only after all `p` witness bytes are proved equal.

Root reduction only runs when the witness is shorter than the query. Zero
length, nonperiodic roots, unequal periods, or queries no longer than the
witness retain their existing paths. Unequal periods keep the general
`p + q - gcd(p,q)` witness and saturating arithmetic. Public range validation
and cross-arena identity rules remain unchanged. Periods are nonzero because
allocated leaves and structural children are nonempty; every piece is bounded
by `period - phase`, keeping its range and additions valid even near `u64::MAX`.

For ordinary units, each bounded delegated comparison uses the existing fast
path with general fallback. The unit itself cannot immediately re-enter the
same root-period reduction, and at most three calls are made. For obviously
compressed units (logical weight exceeding combined stored-node count), the
candidate deliberately keeps **one shared memo table across all phase pieces**.
This avoids discarding reusable proofs on large nonperiodic units. Less obvious
sharing still receives a bounded fast traversal and memoized fallback; the
weight check is a cost heuristic rather than a correctness premise.

The eight diagnostic fixtures support the intended work reductions and cover
the large-unit cache concern. They supplement the implementation agent's
adversarial regression tests; they do not replace the full correctness suite
or prove a universal performance gain. No additional generated-code build or
timed ablation was used to select these conclusions.
