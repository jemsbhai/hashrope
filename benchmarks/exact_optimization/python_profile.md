Python profiling found a narrow, measurable bottleneck: `_leaf_lcp` scanned the
first unequal chunk byte by byte in Python. For a 64 KiB leaf differing only at
the last byte, this consumed over 99% of profiled query time. Construction and
hashing were outside the profile.

The change uses native byte comparisons to narrow a known unequal interval.
The first byte returns directly; an exact 64-byte prefix probe keeps small
offsets out of the large-interval search. Equal halves are skipped, so total
probe work remains linear in the bounded chunk size. The tree traversal,
periodic witness, verified-range cache, profiles, and fingerprints are unchanged.

The isolated leaf-offset diagnostic measured the following medians in
microseconds. These are single-process diagnostic samples, not substitutes for
the paired consumer benchmark reruns.

| Chunk / first mismatch | Baseline | Optimized |
|---|---:|---:|
| 4 KiB / byte 0 | 0.364 | 0.216 |
| 4 KiB / byte 1 | 0.388 | 0.396 |
| 4 KiB / byte 4095 | 121.721 | 1.859 |
| 64 KiB / byte 0 | 0.363 | 0.217 |
| 64 KiB / byte 1 | 0.398 | 0.440 |
| 64 KiB / byte 64 | 2.161 | 4.091 |
| 64 KiB / byte 65535 | 1979.246 | 7.180 |

There is a small local regression at the scout boundary, and no claim that
every offset becomes faster. On the diagnostic stored-corpus workload, the
2 MiB late-mismatch query moved from 1.024 ms to 0.883 ms. Independent equal
trees still spend most of their time in traversal and normalization; this
change does not optimize that path or change its complexity. Small differences
on unchanged paths in these profiles are measurement variation.

`python_profile.baseline.{json,txt}` and `python_profile.optimized.{json,txt}`
record source/corpus hashes, query samples, and cProfile counts. The isolated
offset sweep is in `python_profile.localization.json`; it records both source
hashes and all samples. The `initial_bisection` and `localization_initial` files
retain the diagnostic evidence that motivated the short-prefix probe.

Reproduce after unpacking the recorded baseline archive as described by the
optimization benchmark harness:

```powershell
python benchmarks/exact_optimization/python_profile.py --label baseline --package-root benchmarks/exact_optimization/.cache/baseline/hashrope/packages/python
python benchmarks/exact_optimization/python_profile.py --label optimized
python benchmarks/exact_optimization/python_profile.localize.py
```

Validation: the complete Python core suite passed, **203 tests**, including
all mismatch positions for several short leaf lengths, offsets around 64 KiB
chunk boundaries, multiple differences, checked ranges, cross-profile content,
the default collisions, trillion-copy repeats, and independent self-concat DAGs.
The bounded interval argument was also independently reviewed. Temporary byte
payload is bounded by 192 KiB, including both chunk buffers and probe slices.
