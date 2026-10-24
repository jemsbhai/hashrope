# Exact comparison optimization: measured correction costs

The frozen baseline is the already-correct exact implementation, including the prior collision fixes. Fingerprints and defaults are identical. These results measure optimization, separately from the earlier hash-only-to-exact migration.

## Evidence and method

The archives and manifests capture both complete measured core versions. `unchanged_code.json` verifies six hash/sliding/export modules, 33 Python definitions outside `_leaf_lcp`, and 31 Rust construction/hash/edit functions stayed unchanged. [VALIDATION.md](VALIDATION.md) records 639 passing tests, the no_std release check, and explicit experimental exclusions.

Rust profiling of independent equal 1M-byte trees found 483 query allocations requesting 184,176 bytes, 192,304 ordered-memo key comparisons and zero memo hits. The optimized ordinary-tree path performs zero heap allocations on this fixture, uses no memo entries, and reduces normalization visits from 19,534 to 15,626 while still checking 1M bytes. The compressed-repeat/DAG engine remains available through the bounded fallback. `rust_profile*.json` records diagnostic counters and ablations; their single-process timings are explanatory, not the final performance estimates.

Python profiling isolated the unequal leaf's byte-by-byte scan. Native comparisons now narrow that unequal interval. The tree/cache/periodic logic is unchanged. [python_profile.md](python_profile.md) discloses the local tradeoff: a 64 KiB leaf mismatching at byte 64 increased from 2.161 to 4.091 µs, while its last-byte mismatch fell from 1,979.246 to 7.180 µs. Those are diagnostic offset samples; not every mismatch position improves.

Final measurements alternate frozen baseline and optimized cores by seed/invocation, in fresh processes. The same data, builders, sizes, match/mismatch mix, inner arm order, warmups, oracles and repetitions are retained. Tables report mean ± sample SD over nine process medians, in milliseconds unless labeled otherwise; speedups are ratios of means. Real-pair and R2 queries are first averaged within each process. Construction is excluded from query latency. Biology retains nine paired repetitions within one process per core and its original adaptive batching. Microbenchmarks retain three alternating rounds and seven samples per round. Team load-generating tests were paused; this Windows 11 / Intel Core i9-14900HX workstation had no CPU affinity or frequency locking.

## Rust: edit plus exact target comparison

Each batch performs 50 seeded churn edits and compares each edited value with a prebuilt target, 50% equal and 50% different at the last byte. All 50 boolean outcomes are checked. Fixtures are ASCII, so Hashrope byte offsets and Ropey character offsets coincide. Target preparation and destruction remain outside timing. Ropey is pinned to 1.6.1.

| Bytes / target history | Baseline exact | Optimized exact | Speedup | Ropey exact, optimized cohort |
|---|---:|---:|---:|---:|
| 10,000 / independent | 0.0245 ± 0.0009 | 0.0085 ± 0.0005 | 2.87× | 0.0013 ± 0.0001 |
| 10,000 / shared ancestry | 0.0186 ± 0.0011 | 0.0068 ± 0.0004 | 2.73× | 0.0010 ± 0.0001 |
| 100,000 / independent | 0.1078 ± 0.0031 | 0.0308 ± 0.0015 | 3.50× | 0.0117 ± 0.0007 |
| 100,000 / shared ancestry | 0.0359 ± 0.0020 | 0.0134 ± 0.0006 | 2.68× | 0.0057 ± 0.0005 |
| 1,000,000 / independent | 0.8938 ± 0.0244 | 0.1888 ± 0.0082 | 4.73× | 0.1114 ± 0.0059 |
| 1,000,000 / shared ancestry | 0.0589 ± 0.0016 | 0.0169 ± 0.0007 | 3.49× | 0.0374 ± 0.0046 |

At 1M bytes, optimized Hashrope remains **1.70× slower than Ropey** on independent targets, and is **2.22× faster** with shared ancestry. This is an improvement with explicit limits, not a universal equality win. Independent candidate-only control time is 0.01068→0.01050 ms; shared candidate control is 0.00874→0.00885 ms.

## Paper CPU consumers

EXP-005, mismatching halfway through independently built corpus-derived ropes:

| Bytes | Baseline exact | Optimized exact |
|---:|---:|---:|
| 64,000 | 0.1255 ± 0.0452 | 0.0157 ± 0.0010 |
| 256,000 | 0.0800 ± 0.0047 | 0.0526 ± 0.0033 |
| 1,000,000 | 0.2090 ± 0.0135 | 0.2097 ± 0.0287 |
| 2,000,000 | 0.4492 ± 0.0501 | 0.4391 ± 0.0639 |
| 4,000,000 | 0.8487 ± 0.0366 | 0.8367 ± 0.0261 |
| 8,000,000 | 1.7137 ± 0.0772 | 1.6514 ± 0.0904 |

The exact f=0.5 log-log slope fitted to the six cohort means is 0.591→0.975. The former slope ≤ 0.3 gate fails for the optimized core. Optimization does not establish a logarithmic exact-comparison bound. Fractions 0 and 0.99 are retained in the raw/summary data.

Controlled EXP-017 prefix queries:

| Shared tokens | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |
|---:|---:|---:|---:|---:|
| 1,000 | 0.1219 ± 0.0059 | 0.0038 ± 0.0002 | 0.0347 ± 0.0020 | 0.0023 ± 0.0002 |
| 4,000 | 0.1142 ± 0.0058 | 0.0086 ± 0.0007 | 0.1227 ± 0.0042 | 0.0029 ± 0.0002 |
| 16,000 | 0.1012 ± 0.0046 | 0.0281 ± 0.0010 | 0.4711 ± 0.0170 | 0.0051 ± 0.0003 |
| 64,000 | 0.1735 ± 0.0448 | 0.1037 ± 0.0062 | 2.1016 ± 0.1067 | 0.0180 ± 0.0010 |
| 128,000 | 0.2027 ± 0.0123 | 0.2016 ± 0.0157 | 4.3559 ± 0.7728 | 0.0397 ± 0.0061 |
| 256,000 | 0.4110 ± 0.0252 | 0.4197 ± 0.0271 | 7.8484 ± 0.2705 | 0.0768 ± 0.0103 |
| 512,000 | 0.8323 ± 0.0384 | 0.8333 ± 0.0525 | 15.8612 ± 0.5615 | 0.1847 ± 0.0514 |
| 1,000,000 | 1.8315 ± 0.1606 | 1.6856 ± 0.0790 | 30.3253 ± 0.6414 | 0.3372 ± 0.0960 |
| 2,000,000 | 3.6350 ± 0.2514 | 3.6690 ± 0.2048 | 60.7689 ± 0.8450 | 1.7129 ± 0.1259 |

Optimized exact wins 9/9 paired comparisons at 2M tokens. The original strict crossover criterion passes on this grid; separated mean ± 1 SD bands start at 1000 tokens. This empirical grid result is separate from the retired logarithmic hash-call guard.
The separation threshold equals the smallest measured point. The criterion's below-threshold condition is therefore vacuous; this shows superiority over the sampled grid, not a measured interior crossover.

EXP-015 R2 CPU identification:

| Cached tokens | Baseline exact | Optimized exact | Query build, optimized cohort |
|---:|---:|---:|---:|
| 1,000 | 0.0014 ± 0.0003 | 0.0013 ± 0.0001 | 1.6235 ± 0.0879 |
| 16,000 | 0.0265 ± 0.0054 | 0.0276 ± 0.0088 | 13.4671 ± 0.3200 |
| 64,000 | 0.0975 ± 0.0034 | 0.0967 ± 0.0036 | 51.6162 ± 1.5187 |
| 256,000 | 0.4224 ± 0.0275 | 0.4131 ± 0.0115 | 212.4969 ± 2.8874 |
| 1,000,000 | 1.6755 ± 0.0938 | 1.6876 ± 0.0811 | 803.7662 ± 8.0603 |

Eager query construction is part of an integration that starts with flat tokens; the small prebuilt-query latency does not include it. No 4-A100 serving, TTFT or NVML energy measurement was rerun, and no new GPU claim follows.

## Real conversational pairs and candidate sets

Original 200 pairs per dataset per seed, with independently built ropes:

| Dataset | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |
|---|---:|---:|---:|---:|
| sharegpt_sample.jsonl | 0.0702 ± 0.0109 | 0.0038 ± 0.0004 | 0.0329 ± 0.0029 | 0.0022 ± 0.0004 |
| lmsys_sample.jsonl | 0.0469 ± 0.0082 | 0.0029 ± 0.0004 | 0.0209 ± 0.0040 | 0.0019 ± 0.0003 |

The optimized aggregate means reverse the previous real-pair radix advantage on both datasets; NumPy remains faster. These are query-only results for the selected conversational pairs, not a universal indexing advantage.

K candidates retain one radix lookup versus K pairwise Hashrope/NumPy queries; prefix values are shared but the Hashrope roots are independently constructed.

| Shared tokens / K | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |
|---|---:|---:|---:|---:|
| 64000/1 | 0.2193 ± 0.0719 | 0.1480 ± 0.0510 | 3.1494 ± 1.7192 | 0.0241 ± 0.0149 |
| 64000/10 | 1.9940 ± 0.4565 | 1.5126 ± 0.4380 | 2.6332 ± 0.6128 | 0.1992 ± 0.0614 |
| 64000/100 | 20.9065 ± 4.4562 | 11.9176 ± 0.5368 | 2.8955 ± 0.8027 | 3.6866 ± 0.2873 |
| 512000/1 | 1.0172 ± 0.2221 | 0.8775 ± 0.0921 | 17.9653 ± 3.4465 | 0.1941 ± 0.0909 |
| 512000/10 | 10.2742 ± 1.3644 | 11.3337 ± 2.2663 | 17.6282 ± 2.8181 | 3.9220 ± 1.3471 |

At 64K/K=100, exact remains 4.12× slower than the indexed radix lookup. At 512K/K=10 the optimized exact mean increases, with overlapping sample variability; this cell does not establish an optimization gain.

## Biology consumers

All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. These are whole-corpus median batch times over nine within-process paired repetitions per core, not nine fresh-process estimates.

| Workload | Baseline exact ms | Optimized exact ms | Speedup | Inner repetitions baseline / optimized |
|---|---:|---:|---:|---:|
| regions | 9.1388 | 8.8234 | 1.04× | 3 / 4 |
| gene_reports | 10.4843 | 11.8006 | 0.89× | 2 / 1 |
| panels | 222.6273 | 207.5836 | 1.07× | 1 / 1 |

Adaptive batching is unchanged from the earlier protocol; values normalize the recorded inner repetitions. The gene-report median increases, while the old-hash control also varies substantially. These samples do not support a uniform biology-consumer speedup. Construction is excluded. These timings do not validate clinical efficacy or reproduce chromosome-scale claims.

## Existing structural and operation microbenchmarks

Same fixtures, three alternating rounds × seven samples. Values below are the median [minimum, maximum] of the three within-process medians, matching the original microbenchmark aggregation. These round ranges are not confidence intervals. Full memory/scratch samples remain in `results/microbench.json`.

Python (µs per operation):

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 13331.6500 [12993.2500, 13344.0500] | 13067.5000 [12690.3500, 13621.3000] |
| core/chunked_build_64KiB | 13385.9000 [13173.2500, 13532.6000] | 13232.3000 [12954.0000, 13383.5500] |
| core/concat_balanced_64KiB | 1.2158 [1.2145, 1.2210] | 1.2060 [1.1378, 1.2158] |
| core/split_chunked_64KiB | 119.4750 [119.0800, 119.6750] | 115.9250 [112.7200, 116.4900] |
| core/substr_warm_64KiB | 214.8450 [214.1550, 218.0900] | 207.5750 [207.0150, 209.7350] |
| core/root_hash | 0.0473 [0.0458, 0.0483] | 0.0464 [0.0449, 0.0467] |
| core/repeat_join_2pow32 | 594.5750 [574.3750, 600.5000] | 578.0750 [575.3000, 591.2750] |
| compare/65536/same_root/exact_lcp | 0.2050 [0.1950, 0.2600] | 0.2350 [0.2150, 0.2950] |
| compare/65536/same_root/exact_equal | 0.0850 [0.0850, 0.1050] | 0.0950 [0.0850, 0.1450] |
| compare/65536/same_root/hash_only_lcp | 1694.5667 [1611.6667, 1709.2000] | 1608.6000 [1603.0667, 1716.4667] |
| compare/65536/independent_equal/exact_lcp | 192.7400 [187.3200, 193.8550] | 186.0100 [184.5400, 190.3800] |
| compare/65536/independent_equal/exact_equal | 194.2550 [192.4100, 196.3600] | 186.4050 [184.1450, 193.8200] |
| compare/65536/independent_equal/hash_only_lcp | 1669.6000 [1609.9333, 1773.0000] | 1614.2667 [1606.4000, 1639.1333] |
| compare/65536/early_mismatch/exact_lcp | 6.3050 [6.2750, 6.7350] | 6.0000 [5.9100, 6.3800] |
| compare/65536/early_mismatch/exact_equal | 6.7250 [6.5850, 6.9850] | 6.0550 [5.8900, 6.1950] |
| compare/65536/early_mismatch/hash_only_lcp | 221.0000 [220.0000, 230.6333] | 217.2667 [215.6667, 231.9667] |
| compare/65536/late_mismatch/exact_lcp | 202.2950 [198.8050, 203.0750] | 187.1250 [186.9050, 188.9950] |
| compare/65536/late_mismatch/exact_equal | 201.4800 [200.6200, 210.2450] | 187.3100 [187.2250, 187.3600] |
| compare/65536/late_mismatch/hash_only_lcp | 1723.5667 [1654.0333, 1755.5667] | 1632.6333 [1623.4000, 1662.3333] |
| compare/65536/shared_late_mismatch/exact_lcp | 10.3400 [10.2250, 10.3500] | 10.2050 [10.1000, 10.3950] |
| compare/65536/shared_late_mismatch/exact_equal | 10.3850 [10.0300, 10.8100] | 10.3050 [10.0800, 10.3300] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 1561.8333 [1544.5667, 1604.9000] | 1526.9000 [1522.5333, 1534.0000] |
| compare/1048576/same_root/exact_lcp | 0.2700 [0.2250, 0.4900] | 0.2700 [0.2400, 0.2950] |
| compare/1048576/same_root/exact_equal | 0.1450 [0.1200, 0.2850] | 0.1000 [0.0850, 0.1050] |
| compare/1048576/same_root/hash_only_lcp | 2081.9667 [1779.4667, 3211.7333] | 1795.2667 [1761.8333, 2515.6333] |
| compare/1048576/independent_equal/exact_lcp | 3792.7950 [3342.9750, 3916.5750] | 3590.2950 [3561.3150, 3816.9250] |
| compare/1048576/independent_equal/exact_equal | 3627.6500 [3400.4600, 4026.6650] | 3412.2150 [3233.9050, 3432.3800] |
| compare/1048576/independent_equal/hash_only_lcp | 1775.0333 [1656.8333, 1841.2333] | 1898.8667 [1770.1000, 2054.0667] |
| compare/1048576/early_mismatch/exact_lcp | 8.9350 [8.8050, 11.0100] | 9.6150 [8.4350, 9.7350] |
| compare/1048576/early_mismatch/exact_equal | 9.0200 [8.9500, 9.8950] | 8.9050 [8.5600, 9.3650] |
| compare/1048576/early_mismatch/hash_only_lcp | 232.9333 [230.1667, 257.3667] | 245.1667 [233.8000, 245.4667] |
| compare/1048576/late_mismatch/exact_lcp | 3378.5300 [3268.7500, 3448.2100] | 3492.3150 [3424.2500, 3658.0550] |
| compare/1048576/late_mismatch/exact_equal | 3839.9250 [3525.9400, 3991.6900] | 3364.7150 [3337.1000, 3519.9450] |
| compare/1048576/late_mismatch/hash_only_lcp | 1709.3667 [1678.8667, 2070.4667] | 1858.2000 [1729.3667, 3184.9333] |
| compare/1048576/shared_late_mismatch/exact_lcp | 20.2750 [15.2900, 30.1200] | 14.8100 [13.3000, 25.0350] |
| compare/1048576/shared_late_mismatch/exact_equal | 14.3900 [14.3600, 15.0950] | 15.3650 [14.6900, 24.9200] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 1664.9667 [1638.1333, 2031.2667] | 1718.0667 [1678.4000, 1779.6667] |
| compare/repeat_2TiB/exact_lcp | 2.8720 [2.7495, 4.2460] | 2.9305 [2.8610, 3.0435] |
| compare/repeat_2TiB/exact_equal | 3.0350 [2.7780, 3.6025] | 2.9025 [2.8795, 3.0700] |
| compare/repeat_2TiB/hash_only_lcp | 3501.1500 [2889.2500, 3888.4500] | 3089.8000 [3061.1000, 3401.6000] |
| compare/self_concat_DAG_2TiB/exact_lcp | 53.3310 [50.7520, 90.9330] | 53.1200 [51.5380, 53.3410] |
| compare/self_concat_DAG_2TiB/exact_equal | 52.4050 [50.5210, 58.0910] | 55.5330 [54.2810, 77.3460] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 903.2000 [872.5000, 1244.1000] | 927.3500 [862.1000, 978.2000] |
Rust (µs per operation):

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 300.6000 [294.2750, 301.0300] | 303.4700 [298.9500, 317.8700] |
| core/concat_balanced_64KiB | 0.0194 [0.0151, 0.0211] | 0.0151 [0.0147, 0.0221] |
| core/split_chunked_64KiB | 2.3490 [2.3250, 2.6220] | 2.6020 [2.5250, 2.6030] |
| core/substr_warm_64KiB | 4.4840 [4.3861, 4.5757] | 4.7860 [4.5271, 4.8614] |
| core/root_hash | 0.0006 [0.0005, 0.0006] | 0.0006 [0.0006, 0.0006] |
| core/repeat_join_2pow32 | 6.5590 [6.4330, 6.5980] | 6.4040 [6.2250, 6.9280] |
| compare/65536/same_root/exact_lcp | 0.0025 [0.0025, 0.0025] | 0.0019 [0.0019, 0.0019] |
| compare/65536/same_root/exact_equal | 0.0035 [0.0035, 0.0035] | 0.0039 [0.0038, 0.0039] |
| compare/65536/same_root/hash_only_lcp | 39.6610 [37.5030, 39.8040] | 38.4680 [37.5610, 39.7380] |
| compare/65536/independent_equal/exact_lcp | 11.0140 [10.9800, 11.9850] | 3.6420 [3.4650, 3.9730] |
| compare/65536/independent_equal/exact_equal | 10.9350 [10.7700, 11.2340] | 3.6320 [3.3930, 4.0170] |
| compare/65536/independent_equal/hash_only_lcp | 36.7260 [36.5600, 37.9310] | 37.7210 [37.7050, 38.2270] |
| compare/65536/early_mismatch/exact_lcp | 0.0912 [0.0883, 0.0937] | 0.0615 [0.0612, 0.0616] |
| compare/65536/early_mismatch/exact_equal | 0.0942 [0.0847, 0.0945] | 0.0635 [0.0633, 0.0638] |
| compare/65536/early_mismatch/hash_only_lcp | 5.0830 [4.8540, 5.3360] | 5.2210 [5.1990, 5.3100] |
| compare/65536/late_mismatch/exact_lcp | 10.7270 [10.5940, 11.5580] | 3.9330 [3.4760, 4.1930] |
| compare/65536/late_mismatch/exact_equal | 10.3270 [10.2830, 12.4960] | 3.7650 [3.7110, 4.2430] |
| compare/65536/late_mismatch/hash_only_lcp | 37.0630 [36.8430, 37.1720] | 37.4160 [37.1350, 38.0960] |
| compare/65536/shared_late_mismatch/exact_lcp | 0.2380 [0.2370, 0.2400] | 0.1220 [0.1100, 0.1240] |
| compare/65536/shared_late_mismatch/exact_equal | 0.2440 [0.2430, 0.2440] | 0.1260 [0.1120, 0.1270] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 34.3760 [33.4670, 35.0230] | 35.8830 [35.3030, 35.8940] |
| compare/1048576/same_root/exact_lcp | 0.0025 [0.0022, 0.0025] | 0.0019 [0.0018, 0.0019] |
| compare/1048576/same_root/exact_equal | 0.0035 [0.0032, 0.0035] | 0.0037 [0.0035, 0.0039] |
| compare/1048576/same_root/hash_only_lcp | 41.2190 [39.1530, 41.7220] | 39.6720 [39.4560, 40.9710] |
| compare/1048576/independent_equal/exact_lcp | 390.8570 [389.8690, 394.1980] | 63.5930 [62.2130, 69.0780] |
| compare/1048576/independent_equal/exact_equal | 391.3790 [388.9130, 400.7000] | 62.7580 [62.1300, 63.3770] |
| compare/1048576/independent_equal/hash_only_lcp | 41.1300 [38.6170, 41.9450] | 38.9060 [38.4210, 42.1170] |
| compare/1048576/early_mismatch/exact_lcp | 0.1424 [0.1368, 0.1503] | 0.0822 [0.0820, 0.0912] |
| compare/1048576/early_mismatch/exact_equal | 0.1520 [0.1482, 0.1523] | 0.0853 [0.0852, 0.0935] |
| compare/1048576/early_mismatch/hash_only_lcp | 5.2640 [5.2320, 5.6700] | 5.2150 [5.2140, 5.6830] |
| compare/1048576/late_mismatch/exact_lcp | 396.8400 [390.5180, 409.2660] | 62.5330 [62.4080, 63.5860] |
| compare/1048576/late_mismatch/exact_equal | 395.9010 [385.0920, 397.4190] | 62.7620 [62.0230, 63.3500] |
| compare/1048576/late_mismatch/hash_only_lcp | 39.3850 [38.2530, 39.3970] | 39.8280 [38.8450, 40.1570] |
| compare/1048576/shared_late_mismatch/exact_lcp | 0.4150 [0.3940, 0.4520] | 0.1730 [0.1590, 0.1750] |
| compare/1048576/shared_late_mismatch/exact_equal | 0.4130 [0.4010, 0.4670] | 0.1780 [0.1620, 0.1810] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 36.6580 [36.1190, 37.0930] | 37.8990 [35.9090, 39.0790] |
| compare/repeat_2TiB/exact_lcp | 0.1017 [0.0939, 0.1057] | 0.1070 [0.0977, 0.1083] |
| compare/repeat_2TiB/exact_equal | 0.1056 [0.1050, 0.1075] | 0.1047 [0.1010, 0.1093] |
| compare/repeat_2TiB/hash_only_lcp | 38.4950 [36.3030, 38.5030] | 37.1150 [36.7060, 38.5800] |
| compare/self_concat_DAG_2TiB/exact_lcp | 1.4310 [1.4200, 1.5360] | 1.6460 [1.5630, 1.6480] |
| compare/self_concat_DAG_2TiB/exact_equal | 1.5380 [1.4080, 1.5760] | 1.5070 [1.4730, 1.6550] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 22.8860 [22.8140, 24.0860] | 23.9650 [23.0130, 24.0340] |

Rust layout controls are identical: NodeInner 48 bytes, Arena inline 2136 bytes, 255 nodes for 64 KiB. Python retained/root-object measurements and comparison scratch samples are retained in the summary and raw data. Hash/edit source is unchanged; control timings still vary with measurement conditions and compiled layout. Source equivalence does not establish identical timing. The 2 TiB paths retain their compact behavior, but Rust DAG LCP, split and substring controls show higher optimized medians in these rounds. Python early/late/shared and same-root measurements also move in both directions. This is not a universal no-regression result.

## Provenance and reproduction

[final_source_verification.json](final_source_verification.json) records the post-timing source check against the measured snapshot, including any documentation-only edits.

The original in-memory driver is preserved in `executed_matrix_source.txt`; `run_matrix.py` adds hardened guards for future runs. `results/contexts/*.json` binds each result's executed driver, timed harness/helpers, data hashes, workloads, runtimes and core archives. The initial `execution_context_*.json` files are retained as captured; the final paper context additionally records the verified startup-only source variant. Each final result references its context fingerprint; Rust additionally records compiled-binary hashes and Cargo.lock. The first baseline paper worker ran before three untimed startup import-path-check statements were added. Both exact harness versions are preserved in `paper_cpu_pre_import_guards.txt` and `paper_cpu_import_guards.txt`, matching the recorded source hashes. The paper context explicitly binds both versions, verifies the two timed nested-function ASTs are identical and every remaining source byte matches, and verifies that worker's recorded import path selects the baseline snapshot. All later workers assert the package/module paths at startup. Future resume attempts reject driver, harness, data, workload or core mismatches. Use a new output directory for a fresh reproduction with the hardened driver:

```powershell
python benchmarks/exact_optimization/run_matrix.py --output-root benchmarks/exact_optimization/reproduction_results
python benchmarks/exact_optimization/summarize.py --results-root benchmarks/exact_optimization/reproduction_results --output-dir benchmarks/exact_optimization/reproduction_report
```

The runner restores saved source archives and verifies/reuses the prior token caches, regenerating missing caches with the offline GPT-2 tokenizer. Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies. `preserved_results.json` verifies 25 original consumer_exactness and comparison_safety report/result files remain byte-for-byte preserved. Post-timing documentation-only source edits, if any, are recorded separately with AST/body equivalence to the measured snapshot.
