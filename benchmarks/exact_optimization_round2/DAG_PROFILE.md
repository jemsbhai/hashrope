# Compressed-DAG investigation

These are diagnostic measurements, separate from the final alternating
benchmark cohorts in [README.md](README.md). The fixture represents 2 TiB
using two independently built chains of 40 self-concatenations of `ab`.
No expanded buffer is constructed.

The previous round's memo engine was unchanged, yet its LCP fixture changed
from 1.431 to 1.646 microseconds. Removing the ordinary-tree dispatcher in a
controlled copy gives 1.586 microseconds versus 1.599 with it. That ablation
does **not** explain the earlier 15% change. Code layout, compilation and
run-to-run effects remain possible; the evidence does not isolate a cause.

The profile does identify removable work. A concatenation of the identical
child proves a power of that child, just as an explicit `Repeat` does.
Following only these edges finds a smaller proven period. Two periodic
query ranges need their actual first `p + q - gcd(p, q)` bytes checked;
after that exact witness succeeds, their remaining bytes agree. The query
retains both original offsets, so the witness includes any phase difference.
No fingerprint is consulted, and an ordinary unequal-child concatenation
does not establish a smaller period.

| Diagnostic, microseconds | Before | Query-level reduction |
|---|---:|---:|
| Equal 2 TiB self-concatenation DAG | 1.599 | 0.339 |
| Shifted equal window in those DAGs | 6.404 | 0.332 |
| Last-byte mismatch after editing one DAG | 2.947 | 2.962 |
| Explicit 2 TiB repeats | 0.098 | 0.077 |
| Nonperiodic shared DAG | 4.108 | 4.172 |

These are medians of nine batches in one process per variant, with 2,000
queries per batch. They are not confidence intervals or independent-process
speedup estimates. `dag_before.json`, `dag_bypass.json` and
`dag_periodic_query.json` retain the samples. The same queries are checked
against analytically known results before and during measurement.

Instrumented work for the equal 2 TiB DAG:

| Work | Before | Query-level reduction |
|---|---:|---:|
| General memo-engine calls | 41 | 1 |
| Normalization node visits | 404 | 84 |
| Structural-period node visits | 0 | 82 |
| Memo lookups | 81 | 1 |
| Query heap allocations | 7 | 1 |
| Leaf bytes verified | 2 | 2 |

The new period traversal is included: the combined visit count falls from
404 to 166. Counter-build timings include instrumentation overhead and are
not used for speed estimates. Full counters are in `dag_final_counters.json`.

An initial experiment applied period reduction at every recursive boundary.
Although equal-DAG latency improved, the late-mismatch diagnostic grew to
8.405 microseconds because it repeatedly walked the same chains. That
experiment was rejected. The shipped approach reduces periods once at query
entry, then delegates to the original memoized engine. Its late-mismatch
counters remain unchanged. The rejected diagnostic source and results are
preserved in `dag_rejected_each_level.rs.txt` and
`dag_structural_period.json` for review.

`dag_profile.py` builds isolated copies of a supplied core; its optional
`--bypass` flag is the dispatcher ablation. Instrumentation is never added to
production source. Reproduce a diagnostic with:

```text
python benchmarks/exact_optimization_round2/dag_profile.py --core PATH_TO_ARCHIVED_CORE --label before --output diagnostic.json
```
