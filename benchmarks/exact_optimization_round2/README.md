# Exact comparison optimization: round two

The baseline is the completed first optimization round. This investigation targets the remaining Rust-versus-Ropey gap and the compressed-DAG LCP regression. The [round-one report](../exact_optimization/README.md), its raw results and both original source archives remain preserved.

In the repeated 1M-byte edit-and-compare workload, independent-target exact time falls 19.5%, but remains 1.50× Ropey's time. Shared-history time falls 14.2% and is 2.96× faster than Ropey in this fixture. The 2 TiB self-concatenation DAG's LCP changes from 1.622 to 0.282 µs (median of three round medians).

## Method and scope

Both measured core versions are immutable archives with file manifests. The same consumer runners, corpus inputs, builders, sizes, seed/invocation order, match/mismatch mix, warmups and fixed repetition counts are retained. Nine-process tables report mean ± sample SD across process medians, in milliseconds. Real-pair and R2 queries are first averaged within each process. Speedup ratios use the reported means. Biology retains nine paired repetitions within one process per core, with recorded adaptive batching; its inner counts may differ across cores. Microbenchmarks report the median of three round medians and their minimum/maximum; these ranges are not confidence intervals.

Python is unchanged between the two core snapshots and is rerun as a control. Differences in Python consumer timings therefore do not measure the Rust optimization. Preparation remains outside query timing. CPU-heavy team work is paused for final timing, on the same Windows 11 / Intel Core i9-14900HX workstation, without CPU affinity or frequency locking. No GPU serving, TTFT, NVML energy or chromosome-scale experiment is rerun.

[VALIDATION.md](VALIDATION.md) records the correctness checks and explicit experimental exclusions.

[unchanged_code.json](unchanged_code.json) verifies 7 complete modules, 34 Python definitions and 31 Rust construction/hash/edit functions against the baseline.

## Profiling evidence

[DAG_PROFILE.md](DAG_PROFILE.md) records the compressed-DAG investigation and rejected per-recursion shortcut. The dispatcher-only ablation does not explain the previous 15% regression, so its cause remains unresolved. The selected query-level structural-period proof reduces the equal-DAG diagnostic from 404 normalization visits to 84 normalization plus 82 period-proof visits, 81 memo lookups to one and seven query allocations to one, while verifying the same two-byte witness. These work counters explain a new optimization; single-process diagnostic timings are separate from the final cohorts below.

[fast_path_findings.md](fast_path_findings.md) records the independent-range investigation: direct-child traversal, portable eight-byte mismatch localization, and a fixed-stack frontier for cross-arena spans of at least 8 KiB. Same-arena comparisons retain recursive shared-subtree skips. The frontier uses 1,600 bytes of frontier storage on the stack on this 64-bit target, excluding other locals, and no heap allocation; its diagnostic runs are separate from the final matrix. A 512-byte leaf with a first-byte mismatch measured about 0.7 ns slower in the combined diagnostic variant, so the leaf-localization improvement is not universal.

Exact comparisons remain independent of fingerprints. Existing accepted custom profiles can produce different stored fingerprints for identical bytes in different tree shapes: with modulus 17 and base 2, bytes [0, 13] hash to 16 as one leaf and 0 as a concatenation under the legacy reduction. Rejecting equality on a fingerprint mismatch would therefore break existing exact semantics. This round preserves those profiles and fingerprint values.

## Rust edit plus exact target comparison

Each timed batch performs 50 seeded churn edits and compares the edited rope with a prebuilt target: 50% equal and 50% different at the last byte. Every boolean result is checked. ASCII fixtures make byte and Ropey character offsets coincide. Target preparation/destruction is excluded. Ropey is pinned to 1.6.1. Independent targets and shared ancestry remain separate regimes.

| Bytes / history | Baseline exact | Optimized exact | Baseline / optimized | Paired improvements | Ropey, optimized cohort |
|---|---:|---:|---:|---:|---:|
| 10,000 / independent | 0.0085 ± 0.0005 | 0.0077 ± 0.0004 | 1.11× | 8/9 | 0.0013 ± 0.0002 |
| 10,000 / shared | 0.0070 ± 0.0004 | 0.0068 ± 0.0004 | 1.03× | 6/9 | 0.0010 ± 0.0001 |
| 100,000 / independent | 0.0301 ± 0.0011 | 0.0247 ± 0.0022 | 1.22× | 8/9 | 0.0108 ± 0.0007 |
| 100,000 / shared | 0.0136 ± 0.0006 | 0.0124 ± 0.0007 | 1.09× | 8/9 | 0.0055 ± 0.0004 |
| 1,000,000 / independent | 0.2029 ± 0.0307 | 0.1633 ± 0.0243 | 1.24× | 9/9 | 0.1090 ± 0.0090 |
| 1,000,000 / shared | 0.0175 ± 0.0009 | 0.0150 ± 0.0006 | 1.17× | 9/9 | 0.0443 ± 0.0098 |

At 1M bytes (independent targets), optimized exact Hashrope/Ropey is 1.499. The separate fingerprint-candidate control changes from 0.01162 to 0.01097 ms/cycle.

At 1M bytes (shared ancestry), optimized exact Hashrope/Ropey is 0.338. The separate fingerprint-candidate control changes from 0.00857 to 0.00933 ms/cycle.

For 10K independent targets, the edit-plus-fingerprint-candidate control is 0.0050 ± 0.0003→0.0050 ± 0.0003 ms, versus 0.0013 ± 0.0002 ms for Ropey edit plus exact comparison. The small-document gap therefore includes edit and metadata costs, not only comparison traversal.

These fixture-specific ratios and paired win counts do not establish statistical significance, a universal equality advantage or identical timings on every path.

## Unchanged Python paper controls

EXP-005, mismatch halfway through independently built corpus-derived ropes; all other fractions remain in the raw results:

| Bytes | Baseline core | Optimized core |
|---:|---:|---:|
| 64,000 | 0.0155 ± 0.0004 | 0.0154 ± 0.0007 |
| 256,000 | 0.0567 ± 0.0149 | 0.0556 ± 0.0099 |
| 1,000,000 | 0.2079 ± 0.0128 | 0.2133 ± 0.0296 |
| 2,000,000 | 0.4974 ± 0.1350 | 0.4725 ± 0.1309 |
| 4,000,000 | 0.9450 ± 0.1895 | 0.8549 ± 0.0359 |
| 8,000,000 | 1.8473 ± 0.4158 | 1.7181 ± 0.1032 |

Log-log slopes fitted to the six cohort means are 1.000 and 0.983. The legacy slope gate remains descriptive; no logarithmic exact-LCP bound follows.

Controlled EXP-017:

| Shared tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0081 ± 0.0091 | 0.0045 ± 0.0014 | 0.0363 ± 0.0034 | 0.0024 ± 0.0001 |
| 4,000 | 0.0089 ± 0.0014 | 0.0129 ± 0.0129 | 0.1561 ± 0.0759 | 0.0039 ± 0.0019 |
| 16,000 | 0.0309 ± 0.0074 | 0.0331 ± 0.0098 | 0.6022 ± 0.1999 | 0.0068 ± 0.0027 |
| 64,000 | 0.1766 ± 0.2253 | 0.1261 ± 0.0728 | 3.6645 ± 1.1248 | 0.0269 ± 0.0023 |
| 128,000 | 0.1930 ± 0.0072 | 0.2417 ± 0.0670 | 6.4730 ± 2.9433 | 0.1490 ± 0.1465 |
| 256,000 | 0.4840 ± 0.1544 | 0.4920 ± 0.1209 | 10.5923 ± 2.3529 | 0.1014 ± 0.0318 |
| 512,000 | 1.1253 ± 0.2912 | 1.2539 ± 0.8187 | 17.6651 ± 1.9882 | 0.2162 ± 0.1098 |
| 1,000,000 | 1.6865 ± 0.2245 | 2.2832 ± 0.9061 | 36.6723 ± 5.1902 | 0.4369 ± 0.1671 |
| 2,000,000 | 3.7581 ± 0.5552 | 3.8726 ± 0.4638 | 71.0269 ± 7.7577 | 2.0314 ± 0.5115 |

Optimized-cohort exact wins 9/9 paired comparisons against Radix at 2M tokens. The literal original crossover-grid condition passes; stable band separation starts at 1000 tokens.
That is the smallest measured point, so the below-threshold condition is vacuous; an interior crossover is not observed.

EXP-015 R2 CPU identification:

| Cached tokens | Baseline exact | Optimized-core exact | Excluded query build, optimized cohort |
|---:|---:|---:|---:|
| 1,000 | 0.0015 ± 0.0004 | 0.0015 ± 0.0004 | 1.9220 ± 0.4020 |
| 16,000 | 0.0294 ± 0.0103 | 0.0270 ± 0.0037 | 17.0146 ± 2.8288 |
| 64,000 | 0.1404 ± 0.0592 | 0.1159 ± 0.0349 | 59.9238 ± 4.7303 |
| 256,000 | 0.4857 ± 0.0930 | 0.5044 ± 0.1077 | 252.8374 ± 18.2924 |
| 1,000,000 | 1.9027 ± 0.3677 | 1.8795 ± 0.1730 | 952.9690 ± 43.1254 |

An integration starting with flat tokens must account for query construction; the prebuilt-query figures exclude it.

## Unchanged Python real-pair and candidate-set controls

Original 200 ShareGPT/LMSYS pairs per dataset/seed. Inner timing order remains fixed: radix, exact Hashrope, NumPy. K candidates retain one indexed radix lookup versus K pairwise comparisons with independently built Hashrope roots.

Conversational pairs:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy |
|---|---:|---:|---:|---:|
| sharegpt_sample.jsonl | 0.0039 ± 0.0004 | 0.0041 ± 0.0004 | 0.0361 ± 0.0033 | 0.0024 ± 0.0004 |
| lmsys_sample.jsonl | 0.0030 ± 0.0004 | 0.0033 ± 0.0005 | 0.0227 ± 0.0050 | 0.0022 ± 0.0004 |

Shared tokens / K:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy |
|---|---:|---:|---:|---:|
| 64000/1 | 0.2453 ± 0.0873 | 0.3401 ± 0.1994 | 3.2155 ± 0.8665 | 0.0275 ± 0.0179 |
| 64000/10 | 1.8001 ± 0.6079 | 1.5459 ± 0.4195 | 3.1008 ± 0.7380 | 0.2159 ± 0.0685 |
| 64000/100 | 13.8306 ± 2.6112 | 13.3491 ± 2.2234 | 3.2121 ± 0.9941 | 3.6887 ± 0.4640 |
| 512000/1 | 1.2663 ± 0.4530 | 1.0929 ± 0.3005 | 16.9679 ± 1.2354 | 0.1849 ± 0.0465 |
| 512000/10 | 9.8868 ± 0.7644 | 11.4587 ± 2.5154 | 18.4316 ± 2.8434 | 3.4748 ± 0.6020 |

## Unchanged Python biology controls

All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. Whole-corpus median batch times below are milliseconds, from nine within-process paired repetitions; they are not nine fresh-process estimates. Adaptive inner counts are normalized in the timings.

| Workload | Baseline exact | Optimized-core exact | Old-hash control, baseline / optimized | Inner repetitions, baseline / optimized |
|---|---:|---:|---:|---:|
| regions | 8.5498 | 9.0015 | 0.7834 / 0.9429 | 3 / 3 |
| gene_reports | 9.6823 | 11.9568 | 1.8061 / 2.9430 | 1 / 1 |
| panels | 204.9282 | 199.1552 | 138.3379 / 137.8355 | 1 / 1 |

Construction is excluded. These are timing controls, not clinical-efficacy or chromosome-scale measurements.

## Structural and operation microbenchmarks

Same fixtures, three alternating rounds × seven samples. Values are median [minimum, maximum] of three round medians, in microseconds. These round ranges are not confidence intervals; all individual samples and scratch-memory measurements remain in results/microbench.json.

Python:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 12771.7500 [12755.4000, 12786.0000] | 12646.2000 [12629.2000, 12995.7500] |
| core/chunked_build_64KiB | 12808.9500 [12478.8000, 13008.2500] | 12846.3000 [12608.0000, 13055.6500] |
| core/concat_balanced_64KiB | 1.1227 [1.1090, 1.1558] | 1.1362 [1.1350, 1.2379] |
| core/split_chunked_64KiB | 115.6350 [115.0450, 119.9800] | 112.7250 [110.8850, 116.2650] |
| core/substr_warm_64KiB | 196.1550 [195.8200, 197.4600] | 200.1100 [197.1700, 203.3550] |
| core/root_hash | 0.0448 [0.0440, 0.0449] | 0.0449 [0.0442, 0.0451] |
| core/repeat_join_2pow32 | 572.4750 [556.7250, 577.3000] | 568.8750 [562.4250, 574.7500] |
| compare/65536/same_root/exact_lcp | 0.2000 [0.1900, 0.2500] | 0.2150 [0.1900, 0.2150] |
| compare/65536/same_root/exact_equal | 0.0950 [0.0900, 0.1000] | 0.0900 [0.0850, 0.0950] |
| compare/65536/same_root/hash_only_lcp | 1627.6333 [1588.9333, 1709.8667] | 1607.7333 [1589.6333, 1637.2000] |
| compare/65536/independent_equal/exact_lcp | 182.0650 [181.3650, 191.4250] | 187.1950 [182.8100, 187.4300] |
| compare/65536/independent_equal/exact_equal | 181.4700 [179.9700, 185.7500] | 184.7300 [183.1300, 187.9400] |
| compare/65536/independent_equal/hash_only_lcp | 1614.5667 [1604.9333, 1625.6667] | 1658.8000 [1598.9333, 1670.7667] |
| compare/65536/early_mismatch/exact_lcp | 5.8000 [5.7650, 5.8050] | 5.8600 [5.7400, 5.8850] |
| compare/65536/early_mismatch/exact_equal | 5.9650 [5.8200, 6.0400] | 5.8950 [5.8500, 6.0200] |
| compare/65536/early_mismatch/hash_only_lcp | 211.2333 [209.5333, 213.5333] | 216.4000 [214.1667, 227.8000] |
| compare/65536/late_mismatch/exact_lcp | 183.5150 [180.8300, 191.8800] | 186.2200 [186.1750, 195.1350] |
| compare/65536/late_mismatch/exact_equal | 185.5700 [183.7150, 189.0700] | 186.8500 [182.8850, 189.9850] |
| compare/65536/late_mismatch/hash_only_lcp | 1596.0667 [1585.5000, 1655.4667] | 1592.5333 [1581.5000, 1644.5667] |
| compare/65536/shared_late_mismatch/exact_lcp | 9.3650 [9.2300, 9.9750] | 9.4300 [9.2600, 9.6600] |
| compare/65536/shared_late_mismatch/exact_equal | 9.2800 [9.2150, 9.3550] | 9.4500 [9.4250, 9.4700] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 1488.8333 [1482.0333, 1495.1000] | 1563.4667 [1504.2000, 1574.8667] |
| compare/1048576/same_root/exact_lcp | 0.2100 [0.1950, 0.2850] | 0.2600 [0.2000, 0.2900] |
| compare/1048576/same_root/exact_equal | 0.1000 [0.0950, 0.1150] | 0.0950 [0.0850, 0.1150] |
| compare/1048576/same_root/hash_only_lcp | 1697.4000 [1689.1333, 1722.3000] | 1870.9000 [1685.5333, 3141.6667] |
| compare/1048576/independent_equal/exact_lcp | 3502.8350 [3365.2800, 4040.6350] | 3660.8850 [3243.6500, 4088.6000] |
| compare/1048576/independent_equal/exact_equal | 3380.2050 [3319.0300, 3612.5250] | 3423.4850 [3181.8050, 3489.7800] |
| compare/1048576/independent_equal/hash_only_lcp | 1716.5000 [1691.7667, 3438.5667] | 1731.2667 [1696.6667, 1765.4333] |
| compare/1048576/early_mismatch/exact_lcp | 8.5250 [8.4750, 8.9900] | 8.6400 [8.6150, 8.7950] |
| compare/1048576/early_mismatch/exact_equal | 8.9450 [8.8250, 8.9850] | 8.6100 [8.5750, 8.7900] |
| compare/1048576/early_mismatch/hash_only_lcp | 249.6000 [234.5667, 519.8333] | 237.5333 [232.3333, 246.0667] |
| compare/1048576/late_mismatch/exact_lcp | 3319.0350 [3189.1600, 3325.1500] | 3490.0450 [3395.1500, 3607.6400] |
| compare/1048576/late_mismatch/exact_equal | 3376.6950 [3247.3550, 3379.9800] | 3300.1800 [3290.6100, 3920.9800] |
| compare/1048576/late_mismatch/hash_only_lcp | 1748.1333 [1738.5000, 3837.1667] | 1706.6333 [1677.2333, 2398.2000] |
| compare/1048576/shared_late_mismatch/exact_lcp | 14.0650 [13.5850, 14.3350] | 14.0200 [13.3400, 15.2750] |
| compare/1048576/shared_late_mismatch/exact_equal | 13.8450 [13.4400, 32.7100] | 14.4550 [13.1650, 15.6600] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 1618.0667 [1562.4333, 1657.8333] | 1727.6667 [1588.4000, 3037.6000] |
| compare/repeat_2TiB/exact_lcp | 2.7820 [2.7180, 3.0790] | 2.6535 [2.6395, 5.4795] |
| compare/repeat_2TiB/exact_equal | 3.0060 [2.8035, 3.0520] | 2.8390 [2.7850, 3.0625] |
| compare/repeat_2TiB/hash_only_lcp | 3572.9000 [2935.9000, 3605.7500] | 2932.3000 [2852.0000, 3098.0500] |
| compare/self_concat_DAG_2TiB/exact_lcp | 51.9300 [49.5090, 57.3350] | 50.4080 [49.3020, 71.0200] |
| compare/self_concat_DAG_2TiB/exact_equal | 51.3560 [51.0770, 52.4920] | 50.7330 [49.9600, 60.9900] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 881.6000 [856.6500, 1556.6000] | 890.5500 [854.0500, 891.9500] |

Rust:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 297.3100 [292.2450, 318.9850] | 298.3100 [297.9000, 301.1450] |
| core/concat_balanced_64KiB | 0.0150 [0.0133, 0.0152] | 0.0162 [0.0133, 0.0206] |
| core/split_chunked_64KiB | 2.3660 [2.3270, 2.6000] | 2.5260 [2.3950, 2.7730] |
| core/substr_warm_64KiB | 4.5108 [4.4263, 4.8640] | 4.5253 [4.4296, 4.5446] |
| core/root_hash | 0.0005 [0.0005, 0.0006] | 0.0006 [0.0005, 0.0006] |
| core/repeat_join_2pow32 | 6.2990 [6.0130, 6.5780] | 6.4240 [6.3330, 6.4370] |
| compare/65536/same_root/exact_lcp | 0.0019 [0.0019, 0.0019] | 0.0019 [0.0018, 0.0019] |
| compare/65536/same_root/exact_equal | 0.0039 [0.0037, 0.0039] | 0.0037 [0.0033, 0.0039] |
| compare/65536/same_root/hash_only_lcp | 40.0810 [39.6610, 40.1210] | 38.7140 [36.8710, 39.8160] |
| compare/65536/independent_equal/exact_lcp | 4.1940 [3.6610, 4.2050] | 3.3640 [3.2220, 3.7820] |
| compare/65536/independent_equal/exact_equal | 4.1920 [3.6380, 4.1970] | 3.3460 [3.0800, 3.9130] |
| compare/65536/independent_equal/hash_only_lcp | 38.1820 [38.0610, 40.1800] | 37.2300 [37.0560, 37.9010] |
| compare/65536/early_mismatch/exact_lcp | 0.0614 [0.0555, 0.0617] | 0.0569 [0.0513, 0.0570] |
| compare/65536/early_mismatch/exact_equal | 0.0638 [0.0601, 0.0643] | 0.0532 [0.0531, 0.0592] |
| compare/65536/early_mismatch/hash_only_lcp | 5.3650 [5.2630, 5.3710] | 4.8090 [4.8080, 5.0360] |
| compare/65536/late_mismatch/exact_lcp | 3.7880 [3.4060, 4.2360] | 3.1860 [3.0370, 3.3060] |
| compare/65536/late_mismatch/exact_equal | 3.7680 [3.5190, 3.9110] | 3.2480 [3.0720, 3.3610] |
| compare/65536/late_mismatch/hash_only_lcp | 37.8830 [37.7670, 40.0820] | 37.4200 [37.3510, 37.8410] |
| compare/65536/shared_late_mismatch/exact_lcp | 0.1220 [0.1100, 0.1230] | 0.0860 [0.0840, 0.0940] |
| compare/65536/shared_late_mismatch/exact_equal | 0.1250 [0.1240, 0.1260] | 0.0870 [0.0860, 0.0960] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 35.6040 [35.3340, 37.7630] | 34.5440 [34.4670, 34.9380] |
| compare/1048576/same_root/exact_lcp | 0.0018 [0.0018, 0.0019] | 0.0017 [0.0016, 0.0024] |
| compare/1048576/same_root/exact_equal | 0.0039 [0.0035, 0.0039] | 0.0035 [0.0035, 0.0037] |
| compare/1048576/same_root/hash_only_lcp | 39.4270 [39.3430, 40.7400] | 39.6850 [37.8560, 39.9390] |
| compare/1048576/independent_equal/exact_lcp | 64.2420 [61.8110, 71.2910] | 56.5560 [56.0300, 57.0090] |
| compare/1048576/independent_equal/exact_equal | 65.2740 [63.5590, 67.6650] | 57.1200 [56.1460, 57.5090] |
| compare/1048576/independent_equal/hash_only_lcp | 39.4330 [38.8990, 42.1060] | 38.9740 [38.8310, 40.0690] |
| compare/1048576/early_mismatch/exact_lcp | 0.0915 [0.0897, 0.0917] | 0.0840 [0.0833, 0.0841] |
| compare/1048576/early_mismatch/exact_equal | 0.0933 [0.0850, 0.0938] | 0.0861 [0.0775, 0.0866] |
| compare/1048576/early_mismatch/hash_only_lcp | 5.5150 [5.4170, 5.6830] | 5.3020 [5.1140, 5.5080] |
| compare/1048576/late_mismatch/exact_lcp | 66.7510 [63.3530, 67.1430] | 57.1900 [56.5680, 66.8350] |
| compare/1048576/late_mismatch/exact_equal | 63.6230 [61.5260, 63.9210] | 56.9310 [55.6500, 62.4140] |
| compare/1048576/late_mismatch/hash_only_lcp | 39.0960 [38.2320, 40.9380] | 38.6060 [38.4470, 38.6290] |
| compare/1048576/shared_late_mismatch/exact_lcp | 0.1580 [0.1570, 0.1750] | 0.1130 [0.1120, 0.1130] |
| compare/1048576/shared_late_mismatch/exact_equal | 0.1590 [0.1580, 0.1770] | 0.1140 [0.1140, 0.1150] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 36.9440 [35.7730, 37.7540] | 36.3750 [36.2490, 39.7200] |
| compare/repeat_2TiB/exact_lcp | 0.1057 [0.0958, 0.1157] | 0.0796 [0.0740, 0.0800] |
| compare/repeat_2TiB/exact_equal | 0.1064 [0.0983, 0.1135] | 0.0825 [0.0812, 0.0827] |
| compare/repeat_2TiB/hash_only_lcp | 35.9990 [35.4560, 38.5160] | 35.5010 [35.3070, 35.9320] |
| compare/self_concat_DAG_2TiB/exact_lcp | 1.6220 [1.6090, 1.7150] | 0.2820 [0.2510, 0.2830] |
| compare/self_concat_DAG_2TiB/exact_equal | 1.6270 [1.4570, 1.7170] | 0.3420 [0.3050, 0.3420] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 22.9900 [22.4160, 23.8920] | 22.5410 [22.0120, 22.7570] |

Rust layout remains NodeInner 48 bytes / Arena inline 2136 bytes, with 255 nodes for 64 KiB. Source compatibility does not imply identical operation timing. Compressed-path results and unchanged construction/hash/edit controls must be assessed with their recorded round variability.

## Provenance and reproduction

Every result binds its immutable core, driver, consumer helpers, input hashes, workload and runtime context. Resume and aggregation reject mismatched identities. Rust results also reference compiled-binary hashes and locked dependency manifests. Python worker imports select and check the intended frozen package. Previous reports/results are retained unchanged.

[preserved_results.json](preserved_results.json) verifies 144 prior benchmark source/artifact files remain byte-identical.

[final_source_verification.json](final_source_verification.json) verifies the final working source against the measured optimized archive, including any documentation-only edits.

Use a fresh output directory for a full reproduction:

```powershell
python benchmarks/exact_optimization_round2/run_matrix.py --output-root benchmarks/exact_optimization_round2/reproduction_results
python benchmarks/exact_optimization_round2/provenance.py --results-root benchmarks/exact_optimization_round2/reproduction_results
python benchmarks/exact_optimization_round2/summarize.py --results-root benchmarks/exact_optimization_round2/reproduction_results --output-dir benchmarks/exact_optimization_round2/reproduction_report
```

The runner restores saved core archives, verifies corpus-token caches and can regenerate missing caches with the offline GPT-2 tokenizer. The original local paper/biology datasets remain required. Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies.
