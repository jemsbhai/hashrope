# Edit and metadata decomposition

This diagnostic identifies work behind the small-input Rust gap. It is separate
from the unchanged nine-process final benchmark matrix. The selected core keeps
the default fingerprint values and exact-comparison implementation, groups the
default polynomial evaluation into eight-byte blocks, and removes redundant
full-leaf copies before splitting, partial hashing, and lazy materialization.

The principal evidence is [baseline_final](edit_profile.baseline_final.json)
and [selected](edit_profile.selected.json). The earlier baseline-only diagnostic
is preserved in [baseline](edit_profile.baseline.json). Every result records
source and binary SHA-256 values, the exact locked dependencies, compiler
version, and all raw samples. Each final core also has a `.sources.zip` archive
containing the clean core, diagnostic Rust/Python source, and build manifest.

## Protocol and attribution limits

The generator and edit operations reproduce the original 50-step seeded ASCII
churn workload: seeds 1/2/3, equal and final-byte-different targets in equal
proportions, independent targets and shared ancestors. The original 10K/100K/1M
sizes remain; **1K is an additional exploratory cell**. The baseline and selected
core alternate order by seed. Values below are the mean and sample standard
deviation of three fresh-process medians, each based on three timed samples.
They are diagnostic dispersion, not confidence intervals or a nine-process
confirmatory result.

Edit-only, comparison-only, and combined are separately prepared experiments.
Comparison-only retains all 50 edited roots, then repeats their 50 queries three
times. This changes locality and, for Ropey, copy-on-write behavior during
untimed preparation. It is not valid to add the phase times or subtract them
to obtain an exact attribution of combined latency. The diagnostic also omits
the original Ropey arm's unused untimed target-fingerprint preprocessing.
Setup, byte-oracle verification, and final arena/rope teardown are outside the
timer; allocations and deallocations performed by an edit remain inside it.
Expected equality masks and final/intermediate byte states are checked.

Timing binaries contain neither the counting allocator nor the hash counters.
Separate instrumented binaries count one sample per seed. Counts cover the
actual hash entry calls and bytes, concatenation hashes, inline/ordered-cache
power hits and misses, leaf splits and their redundant source-leaf clone bytes,
new arena nodes, and allocator requests. Allocation sizes are requested traffic,
**not peak live memory, retained memory, or RSS**. Reallocation counts and full
requested replacement sizes are recorded separately.

Lazy Hashrope suppresses materialization of fingerprint metadata and is included
only as a diagnostic. Eager-minus-lazy is not an isolated hash-cost measurement:
the modes also change cache and allocation behavior. Query-only controls and
Ropey show variation even when their relevant source is unchanged.

## Small-input costs

Times are microseconds per operation, mean ± sample SD across three process
medians. Both phases use the selected run's Ropey control for comparison.

| Bytes / target history | Eager edit, baseline → selected | Ropey edit | Combined exact, baseline → selected | Ropey combined |
|---|---:|---:|---:|---:|
| 1K / independent, exploratory | 2.284 ± 0.040 → 1.373 ± 0.215 | 0.401 ± 0.019 | 3.272 ± 0.237 → 2.556 ± 0.144 | 0.619 ± 0.011 |
| 1K / shared, exploratory | 2.333 ± 0.348 → 1.586 ± 0.076 | 0.351 ± 0.019 | 3.665 ± 0.252 → 2.737 ± 0.225 | 0.573 ± 0.021 |
| 10K / independent | 5.065 ± 0.236 → 1.983 ± 0.068 | 0.437 ± 0.009 | 7.881 ± 0.328 → 4.661 ± 0.402 | 1.416 ± 0.167 |
| 10K / shared | 4.759 ± 0.123 → 1.901 ± 0.014 | 0.431 ± 0.049 | 6.645 ± 0.311 → 4.019 ± 0.086 | 1.114 ± 0.093 |
| 100K / independent | 8.090 ± 0.085 → 3.733 ± 0.085 | 0.781 ± 0.021 | 24.517 ± 1.643 → 19.019 ± 0.223 | 9.863 ± 0.845 |
| 100K / shared | 8.753 ± 0.612 → 4.057 ± 0.118 | 1.190 ± 0.256 | 13.243 ± 0.182 → 8.139 ± 0.324 | 5.188 ± 0.073 |
| 1M / independent | 10.573 ± 0.322 → 5.411 ± 0.090 | 1.212 ± 0.065 | 152.771 ± 13.256 → 144.747 ± 5.835 | 106.997 ± 4.700 |
| 1M / shared | 8.620 ± 0.511 → 3.394 ± 0.135 | 0.927 ± 0.051 | 16.049 ± 1.636 → 9.854 ± 0.107 | 34.711 ± 0.238 |

The selected edit path is about 2.55× faster at 10K independent targets and
2.50× faster at 10K shared targets in this diagnostic. The small-input gap to
Ropey remains. The 1K exploratory combined cases improve, but remain about
4.13× and 4.77× slower than their Ropey controls. These ratios are ratios of
reported means, not confidence bounds.

Comparison-only at 10K independent targets changes from 2.238 ± 0.061 to
2.116 ± 0.079 microseconds, while the shared case changes from 1.717 ± 0.106
to 1.778 ± 0.025. At 1M independent it changes from 126.822 ± 3.440 to
127.153 ± 4.544. The exact comparator did not change in this round; these
controls do not establish a new comparison speedup. Lazy edit-only at 10K
independent changes from 1.105 ± 0.097 to 1.002 ± 0.118 microseconds; other lazy
cells have small gains or losses. The major consistent change is in eager edit
work, while long independent comparisons still dominate the larger workload.

## Deterministic work counters

These are mean counts **per edit** across seeds, for independent eager inputs.
The shared-history inputs have identical hash/split/node counts. All hash,
power, split, and new-node counts below are unchanged by the selected core;
the temporary clone allocations disappear.

| Bytes | New nodes | Leaf hash calls / bytes | Concat hashes | Power inline / tree hits / misses | Full-leaf clone bytes, before → after |
|---|---:|---:|---:|---:|---:|
| 1K, exploratory | 36.07 | 6.63 / 199.95 | 29.45 | 23.18 / 6.27 / 0 | 182.99 → 0 |
| 10K | 40.41 | 6.92 / 708.95 | 33.49 | 9.89 / 23.60 / 0 | 692.40 → 0 |
| 100K | 49.98 | 6.99 / 1,213.83 | 42.99 | 3.93 / 39.06 / 0 | 1,197.74 → 0 |
| 1M | 65.47 | 7.00 / 1,316.19 | 58.47 | 2.47 / 56.01 / 0 | 1,299.49 → 0 |

| Bytes | Fresh allocation requests, before → after | Fresh allocation bytes, before → after | Reallocation requests / replacement bytes, unchanged |
|---|---:|---:|---:|
| 1K, exploratory | 9.44 → 6.63 | 382.95 → 199.95 | 0.02 / 3,932.16 |
| 10K | 9.88 → 6.92 | 1,401.35 → 708.95 | 0 / 0 |
| 100K | 9.98 → 6.99 | 2,411.57 → 1,213.83 | 0.02 / 7,864.32 |
| 1M | 10.00 → 7.00 | 2,615.69 → 1,316.19 | 0.02 / 15,728.64 |

At 10K, eliminating the source-leaf clone removes 2.96 allocations and 692.40
requested bytes per edit, approximately 30% of fresh allocation requests.
Arena-node creation and node-array growth do not decrease. Larger cells retain
the same amortized reallocation traffic, so the fresh-allocation byte reduction
must not be presented as an equivalent reduction in total memory traffic or
memory footprint.

The original untimed edit warmup produces **zero power-cache misses** in every
measured edit cell. This is evidence against optimizing cold exponentiation for
this workload. Ordered-cache hits remain, but their counts alone do not isolate
their latency. A repeatedly accessed tiny-cache concat primitive is not a proxy
for the editing workload's larger, mixed-exponent cache.

## Fingerprint primitive

Both cores evaluate the same default polynomial. The selected eight-byte
evaluation leaves the generic/custom-profile path unchanged. The independent
scalar and immutable-baseline regression tests supply the compatibility oracle;
the diagnostic's same-hasher result check is not by itself that proof.

| Bytes hashed | Baseline ns/call | Selected ns/call | Baseline / selected |
|---:|---:|---:|---:|
| 1 | 1.830 | 2.253 | 0.81× |
| 16 | 61.644 | 5.077 | 12.14× |
| 256 | 1,071.032 | 86.221 | 12.42× |
| 512 | 2,168.233 | 168.893 | 12.84× |

There is a small counterexample to a universal speedup: the one-byte primitive
is about **0.423 ns slower** in this diagnostic. The 16–512-byte gains and the
unchanged 709–1,316 hashed bytes per edit support optimizing leaf hashing; the
paired edit measurements demonstrate the resulting workload improvement.
They do not imply that all fingerprint operations improve by the primitive's
factor. Existing exact-comparison behavior and custom-profile fingerprints
remain separate compatibility requirements.

## Reproduction and follow-up

Build each diagnostic core separately, then alternate the binaries. Use fresh
labels and an empty output directory for a new run; do not overwrite completed
or interrupted diagnostic artifacts.

```text
python edit_profile.py --core-root BASELINE_CORE --label baseline_reproduction --output-dir NEW_DIRECTORY --build-only
python edit_profile.py --core-root SELECTED_CORE --label selected_reproduction --output-dir NEW_DIRECTORY --build-only
python edit_profile_pair.py --labels baseline_reproduction selected_reproduction --output-dir NEW_DIRECTORY
```

`edit_profile.py` copies and instruments the supplied crate in a temporary
project. Production files are never changed. Instrumentation is enabled only in
the counter binary; the timing binary uses the ordinary allocator. The pair
driver verifies input and compiled-binary hashes before executing. The original
final benchmark protocol remains unchanged and supplies the confirmatory
performance results in the main report.

The counters support the selected default-profile grouping and borrowed-leaf
copy removal. Further metadata optimization should first measure the warm
ordered-cache and node-construction work that remains; the present evidence
does not justify changing fingerprint parameters, tree shape, persistence, or
the exact API's hash-independent decisions.
