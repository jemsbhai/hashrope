# Exact comparison and profile-check benchmarks

Measured on 2026-10-23. These are local microbenchmarks, not universal performance guarantees.

Baseline: `100749d1cef9603465f4600b814c76d2ad6e9a6e`. Current source SHA-256 values and all raw samples are in `results.json`.
The measured sources remained unchanged during the run: **True**.

Environment: Intel(R) Core(TM) i9-14900HX; Windows-11-10.0.26300-SP0; Python 3.12.2; rustc 1.94.0 (4a4ef493e 2026-03-02). Rust uses Cargo release defaults. No CPU pinning, frequency locking, or process isolation was applied.

3 rounds alternate baseline/current and current/baseline. Each round has seven samples after warm-up. Tables give the median of the round medians; brackets show the minimum and maximum round median. Times are microseconds per operation. Fixture creation is excluded unless the operation is construction. GC is disabled only during Python timing. All comparison fixtures are checked against known exact answers before timing.

## Existing operations

Python default construction remains one leaf; its other nonrepeat fixtures are manually built from 512-byte leaves. Rust construction already chunks to 512-byte leaves. Python construction reuses a hasher; Rust construction includes a fresh arena and its disposal. Do not compare absolute timings across languages. Concat/split retain their newly allocated nodes in Rust's arena; Python discards each returned result. Substring and exponent caches are warmed. These results do not cover lazy first-query materialization.

### Python

| Operation | Baseline µs | Current µs | Current / baseline |
|---|---:|---:|---:|
| from_bytes_64KiB | 13217.050 [12370.500–13251.850] | 12987.250 [12907.750–13238.450] | 0.983× |
| chunked_build_64KiB | 13319.050 [13024.450–13378.000] | 13296.050 [12929.200–13445.400] | 0.998× |
| concat_balanced_64KiB | 0.944 [0.916–0.970] | 1.167 [1.151–1.195] | 1.236× |
| split_chunked_64KiB | 115.855 [109.070–117.390] | 115.160 [114.895–117.885] | 0.994× |
| substr_warm_64KiB | 217.670 [215.470–221.215] | 208.495 [206.510–214.850] | 0.958× |
| root_hash | 0.048 [0.045–0.049] | 0.046 [0.046–0.046] | 0.959× |
| repeat_join_2pow32 | 597.650 [556.625–601.900] | 579.550 [577.000–609.725] | 0.970× |

### Rust

| Operation | Baseline µs | Current µs | Current / baseline |
|---|---:|---:|---:|
| from_bytes_64KiB | 310.550 [310.365–311.210] | 318.405 [304.210–323.995] | 1.025× |
| concat_balanced_64KiB | 0.013 [0.013–0.015] | 0.015 [0.014–0.021] | 1.093× |
| split_chunked_64KiB | 2.394 [2.384–2.586] | 2.587 [2.580–2.616] | 1.081× |
| substr_warm_64KiB | 4.489 [4.471–4.818] | 4.782 [4.648–4.817] | 1.065× |
| root_hash | 0.001 [0.001–0.001] | 0.001 [0.001–0.001] | 0.960× |
| repeat_join_2pow32 | 6.287 [6.244–6.690] | 6.376 [6.305–6.473] | 1.014× |

## Exact comparisons

Hash-only LCP is the existing prefix-hash binary-search strategy using current code; its speed does not make it an exact correctness substitute. The equality API is also timed in the raw results. Byte inputs deterministically cycle through 0–255; results need validation on application traces before deployment decisions. Shared late mismatch reuses the original prefix after replacing its last byte. The Repeat case compares independent `ab × 2^40` with `abab × 2^39`; the DAG case independently self-concatenates two seeds 40 times. Both represent 2 TiB without expansion.

### Python

| Input | Exact LCP µs | Hash-only LCP µs |
|---|---:|---:|
| 1048576/same_root | 0.245 [0.240–0.315] | 1764.033 [1705.400–1820.467] |
| 1048576/independent_equal | 3230.005 [3103.875–3238.815] | 1768.000 [1727.100–1836.133] |
| 1048576/early_mismatch | 9.630 [8.845–9.780] | 233.167 [233.133–258.233] |
| 1048576/late_mismatch | 3209.520 [3195.085–3260.995] | 1801.133 [1670.000–1868.533] |
| 1048576/shared_late_mismatch | 14.675 [13.610–14.770] | 1690.233 [1578.700–1708.700] |
| repeat_2TiB | 2.697 [2.695–2.720] | 3020.450 [2988.600–3051.750] |
| self_concat_DAG_2TiB | 52.729 [52.011–54.883] | 928.850 [907.900–931.650] |

### Rust

| Input | Exact LCP µs | Hash-only LCP µs |
|---|---:|---:|
| 1048576/same_root | 0.003 [0.002–0.003] | 41.695 [38.763–41.724] |
| 1048576/independent_equal | 392.279 [392.114–404.656] | 41.825 [39.634–41.950] |
| 1048576/early_mismatch | 0.150 [0.135–0.151] | 5.429 [5.177–5.720] |
| 1048576/late_mismatch | 396.050 [392.863–406.338] | 39.129 [39.031–40.078] |
| 1048576/shared_late_mismatch | 0.399 [0.395–0.446] | 36.608 [36.559–37.248] |
| repeat_2TiB | 0.103 [0.103–0.104] | 36.327 [35.989–37.643] |
| self_concat_DAG_2TiB | 1.527 [1.379–1.528] | 23.126 [22.351–23.265] |

## Memory

Python `tracemalloc` starts after creating the input bytes. Construction figures include new nodes, chunk slices, hasher/cache, and profiles; the single-leaf input buffer is reused and therefore excluded. These are live traced Python allocations, not process RSS.

| Python construction | Baseline retained bytes | Current retained bytes | Change |
|---|---:|---:|---:|
| single_leaf_64KiB | 636 | 700 | +64 (+10.06%) |
| chunked_64KiB | 105660 | 107756 | +2096 (+1.98%) |

The profile adds one pointer per Python node; the shared profile/hasher has a small fixed cost. Scratch peaks below include temporary memo/stack entries and slice buffers. They exclude existing ropes and warm caches.

| Python query | Current peak scratch bytes |
|---|---:|
| 1048576/independent_equal/exact_lcp | 492808 |
| 1048576/independent_equal/hash_only_lcp | 1612 |
| 1048576/shared_late_mismatch/exact_lcp | 4152 |
| 1048576/shared_late_mismatch/hash_only_lcp | 1607 |
| repeat_2TiB/exact_lcp | 1096 |
| self_concat_DAG_2TiB/exact_lcp | 17632 |

Rust reports layout/node counts only, not allocator totals or exact-query scratch peaks:

| Rust metric | Baseline | Current |
|---|---:|---:|
| node_inner_bytes | 48 | 48 |
| arena_inline_bytes | 2120 | 2136 |
| nodes_for_64KiB | 255 | 255 |

## Reproduce

From the repository root:

```powershell
python benchmarks/comparison_safety/run_benchmarks.py --baseline-ref 100749d1cef9603465f4600b814c76d2ad6e9a6e
python benchmarks/comparison_safety/summarize.py
```

The runner extracts the baseline into a fresh temporary directory and builds separate Rust binaries. It does not change either checkout or inspect the computational-biology submodule. The individual `bench_python.py` and `bench_rust.py` launchers also accept `--package PATH` and `--exact`, allowing focused reruns. Remove temporary benchmark directories manually when no longer needed.
