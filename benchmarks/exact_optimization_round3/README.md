# Edit and metadata optimization: round three

The baseline is the completed second optimization round. This investigation isolates edit and metadata costs behind the remaining Rust-versus-Ropey gap, especially on small inputs. The [round-one report](../exact_optimization/README.md) and [round-two report](../exact_optimization_round2/README.md), their raw results and all original source archives remain preserved.

At 10K bytes, independent-target edit-and-compare changes from 7.803 to 4.709 µs (1.66× speedup), and shared history from 7.108 to 3.795 µs (1.87×), with 9/9 and 9/9 paired improvements respectively. Exact comparisons plus edits still take 3.92× and 3.82× Ropey's time at 10K bytes. In the repeated 1M-byte workload, independent-target exact time falls 3.6% and is 1.35× Ropey's time. Shared-history time falls 34.0% and takes 25.9% of Ropey's time in this fixture. The 2 TiB self-concatenation DAG's LCP changes from 0.281 to 0.251 µs (median of three round medians).

## Method and scope

Both measured core versions are immutable archives with file manifests. The same consumer runners, corpus inputs, builders, sizes, seed/invocation order, match/mismatch mix, warmups and fixed repetition counts are retained. Nine-process tables report mean ± sample SD across process medians, in milliseconds. Real-pair and R2 queries are first averaged within each process. Speedup ratios use the reported means. Biology retains nine paired repetitions within one process per core, with recorded adaptive batching; its inner counts may differ across cores. Microbenchmarks report the median of three round medians and their minimum/maximum; these ranges are not confidence intervals.

Python is unchanged between the two core snapshots and is rerun as a control. Differences in Python consumer timings therefore do not measure the Rust optimization. Preparation remains outside query timing. CPU-heavy team work is paused for final timing, on the same Windows 11 / Intel Core i9-14900HX workstation, without CPU affinity or frequency locking. No GPU serving, TTFT, NVML energy or chromosome-scale experiment is rerun.

[VALIDATION.md](VALIDATION.md) records the correctness checks and explicit experimental exclusions.

[unchanged_code.json](unchanged_code.json) verifies 6 complete modules, 34 Python definitions and storage/default declarations against the baseline. It audits 97 existing Rust production functions, recording 4 selected changes and 1 addition explicitly; all other production function bodies and existing signatures remain unchanged. Source checks establish change scope; differential baseline tests establish the preserved fingerprint semantics.

## Profiling evidence

[edit_profile.md](edit_profile.md) separates eager/lazy edits, exact comparisons, combined operations and fingerprint primitives. Its three-seed process diagnostics retain the original sizes and add a clearly labeled exploratory 1K cell. Separately prepared phases have different locality and retained-state behavior; their timings cannot be added or subtracted for exact attribution.

At 10K independent targets, diagnostic eager editing falls from 5.065 to 1.983 µs while comparison-only changes from 2.238 to 2.116 µs. The exact comparator is unchanged. Leaf hashing still processes about 709 bytes in 6.92 calls per edit; eliminating redundant source-leaf copies removes 2.96 fresh allocation requests and 692.40 requested bytes per edit. Node creation and arena growth do not decrease. Those allocation figures describe requested traffic, not peak or retained memory. Warmed power-cache misses are zero, so cold exponentiation is not the measured bottleneck.

Default-profile eight-byte polynomial grouping reduces the diagnostic 512-byte hash from 2,168 to 169 ns, but the one-byte hash increases from 1.830 to 2.253 ns (0.423 ns). The generic/custom-profile path remains scalar. [REVIEW.md](REVIEW.md) gives the arithmetic bounds and independent compatibility review. [hash_kernel_diagnostic.json](hash_kernel_diagnostic.json) and [edit_variants.json](edit_variants.json) preserve primitive/combined ablations, including the rejected extra modular-arithmetic specialization. Their diagnostic repetitions and source provenance remain distinct from the final matrix below.

Exact comparisons remain independent of fingerprints. Existing accepted custom profiles can produce different stored fingerprints for identical bytes in different tree shapes: with modulus 17 and base 2, bytes [0, 13] hash to 16 as one leaf and 0 as a concatenation under the legacy reduction. Rejecting equality on a fingerprint mismatch would therefore break existing exact semantics. This round preserves those profiles and fingerprint values.

## Rust edit plus exact target comparison

Each timed batch performs 50 seeded churn edits and compares the edited rope with a prebuilt target: 50% equal and 50% different at the last byte. Every boolean result is checked. ASCII fixtures make byte and Ropey character offsets coincide. Target preparation/destruction is excluded. Ropey is pinned to 1.6.1. Independent targets and shared ancestry remain separate regimes. Both fingerprint-candidate arms use Hashrope from their measured core version, including Ropey's rehash-candidate arm, so they can benefit from this round's hash optimization. Ropey exact uses its unchanged equality implementation.

| Bytes / history | Baseline exact | Optimized exact | Baseline / optimized | Paired improvements | Ropey, optimized cohort |
|---|---:|---:|---:|---:|---:|
| 10,000 / independent | 0.0078 ± 0.0003 | 0.0047 ± 0.0002 | 1.66× | 9/9 | 0.0012 ± 0.0001 |
| 10,000 / shared | 0.0071 ± 0.0002 | 0.0038 ± 0.0002 | 1.87× | 9/9 | 0.0010 ± 0.0001 |
| 100,000 / independent | 0.0246 ± 0.0022 | 0.0194 ± 0.0006 | 1.27× | 9/9 | 0.0109 ± 0.0004 |
| 100,000 / shared | 0.0124 ± 0.0008 | 0.0078 ± 0.0005 | 1.60× | 9/9 | 0.0051 ± 0.0002 |
| 1,000,000 / independent | 0.1491 ± 0.0027 | 0.1437 ± 0.0042 | 1.04× | 8/9 | 0.1069 ± 0.0048 |
| 1,000,000 / shared | 0.0144 ± 0.0009 | 0.0095 ± 0.0004 | 1.51× | 9/9 | 0.0367 ± 0.0047 |

At 1M bytes (independent targets), optimized exact Hashrope/Ropey is 1.345. The separate fingerprint-candidate control changes from 0.01018 to 0.00518 ms/cycle.

At 1M bytes (shared ancestry), optimized exact Hashrope/Ropey is 0.259. The separate fingerprint-candidate control changes from 0.00855 to 0.00322 ms/cycle.

For 10K independent targets, the edit-plus-fingerprint-candidate arm is 0.0051 ± 0.0003→0.0020 ± 0.0001 ms, versus 0.0012 ± 0.0001 ms for Ropey edit plus exact comparison. The decomposition profiles separately measure edit and metadata costs.

These fixture-specific ratios and paired win counts do not establish statistical significance, a universal equality advantage or identical timings on every path.

Fingerprint-candidate arms from the same runs (milliseconds per cycle); these are distinct from exact equality:

| Bytes / history | Hashrope candidate, baseline | Hashrope candidate, optimized | Ropey rehash candidate, baseline | Ropey rehash candidate, optimized |
|---|---:|---:|---:|---:|
| 10,000 / independent | 0.0051 ± 0.0003 | 0.0020 ± 0.0001 | 0.0469 ± 0.0015 | 0.0054 ± 0.0003 |
| 10,000 / shared | 0.0049 ± 0.0002 | 0.0019 ± 0.0002 | 0.0466 ± 0.0006 | 0.0050 ± 0.0005 |
| 100,000 / independent | 0.0086 ± 0.0005 | 0.0037 ± 0.0002 | 0.4385 ± 0.0170 | 0.0429 ± 0.0010 |
| 100,000 / shared | 0.0088 ± 0.0007 | 0.0037 ± 0.0002 | 0.4392 ± 0.0164 | 0.0422 ± 0.0012 |
| 1,000,000 / independent | 0.0102 ± 0.0004 | 0.0052 ± 0.0002 | 4.3588 ± 0.1600 | 0.4291 ± 0.0247 |
| 1,000,000 / shared | 0.0085 ± 0.0004 | 0.0032 ± 0.0002 | 4.3025 ± 0.0431 | 0.4103 ± 0.0149 |

## Unchanged Python paper controls

EXP-005, mismatch halfway through independently built corpus-derived ropes; all other fractions remain in the raw results:

| Bytes | Baseline core | Optimized core |
|---:|---:|---:|
| 64,000 | 0.0151 ± 0.0007 | 0.0150 ± 0.0008 |
| 256,000 | 0.0525 ± 0.0030 | 0.0522 ± 0.0028 |
| 1,000,000 | 0.2138 ± 0.0549 | 0.2173 ± 0.0302 |
| 2,000,000 | 0.4081 ± 0.0180 | 0.4481 ± 0.0935 |
| 4,000,000 | 0.8633 ± 0.2004 | 0.8272 ± 0.0430 |
| 8,000,000 | 1.6273 ± 0.0672 | 1.6441 ± 0.1035 |

Log-log slopes fitted to the six cohort means are 0.979 and 0.982. The historical ≤0.3 gate fails for the baseline and fails for the optimized cohort. This gate remains descriptive; no logarithmic exact-LCP bound follows.

Controlled EXP-017:

| Shared tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0041 ± 0.0004 | 0.0042 ± 0.0010 | 0.0405 ± 0.0125 | 0.0029 ± 0.0014 |
| 4,000 | 0.0082 ± 0.0003 | 0.0101 ± 0.0044 | 0.1205 ± 0.0033 | 0.0029 ± 0.0001 |
| 16,000 | 0.0277 ± 0.0010 | 0.0327 ± 0.0096 | 0.6637 ± 0.3129 | 0.0066 ± 0.0028 |
| 64,000 | 0.1005 ± 0.0045 | 0.1085 ± 0.0296 | 2.5025 ± 0.6814 | 0.0198 ± 0.0059 |
| 128,000 | 0.1927 ± 0.0062 | 0.1992 ± 0.0068 | 5.0478 ± 2.0782 | 0.0402 ± 0.0057 |
| 256,000 | 0.4448 ± 0.1000 | 0.4068 ± 0.0230 | 7.9662 ± 0.3222 | 0.1048 ± 0.0663 |
| 512,000 | 0.8311 ± 0.0709 | 0.8159 ± 0.0458 | 16.0834 ± 1.1422 | 0.1737 ± 0.0488 |
| 1,000,000 | 1.6315 ± 0.1154 | 1.7658 ± 0.4428 | 30.9313 ± 2.0240 | 0.4753 ± 0.4322 |
| 2,000,000 | 3.4247 ± 0.2175 | 3.5261 ± 0.3160 | 64.4321 ± 6.7456 | 1.7897 ± 0.3307 |

Optimized-cohort exact wins 9/9 paired comparisons against Radix at 2M tokens. The literal original crossover-grid condition passes; stable mean ± one sample SD band separation starts at 1000 tokens. These are dispersion bands, not confidence intervals or a significance test.
That is the smallest measured point, so the below-threshold condition is vacuous; an interior crossover is not observed.

EXP-015 R2 CPU identification:

| Cached tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | Excluded query build, optimized cohort |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0012 ± 0.0001 | 0.0012 ± 0.0000 | 0.0349 ± 0.0025 | 1.7231 ± 0.3002 |
| 16,000 | 0.0260 ± 0.0025 | 0.0270 ± 0.0035 | 0.5146 ± 0.0695 | 13.8848 ± 0.8696 |
| 64,000 | 0.1193 ± 0.0374 | 0.1074 ± 0.0255 | 2.7463 ± 0.9148 | 53.7571 ± 3.9221 |
| 256,000 | 0.4249 ± 0.0465 | 0.4099 ± 0.0132 | 8.7402 ± 1.5370 | 224.4793 ± 21.0266 |
| 1,000,000 | 1.7198 ± 0.2106 | 1.8923 ± 0.6131 | 32.4642 ± 3.7963 | 850.9699 ± 70.2959 |

An integration starting with flat tokens must account for query construction; the prebuilt-query figures exclude it.

## Unchanged Python real-pair and candidate-set controls

Original 200 ShareGPT/LMSYS pairs per dataset/seed. Inner timing order remains fixed: radix, exact Hashrope, NumPy. K candidates retain one indexed radix lookup versus K pairwise comparisons with independently built Hashrope roots.

Conversational pairs:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---|---:|---:|---:|---:|
| sharegpt_sample.jsonl | 0.0039 ± 0.0003 | 0.0042 ± 0.0005 | 0.0371 ± 0.0045 | 0.0026 ± 0.0005 |
| lmsys_sample.jsonl | 0.0031 ± 0.0006 | 0.0033 ± 0.0006 | 0.0231 ± 0.0040 | 0.0022 ± 0.0007 |

Shared tokens / K:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---|---:|---:|---:|---:|
| 64000/1 | 0.2989 ± 0.1837 | 0.2102 ± 0.1236 | 2.9610 ± 0.7499 | 0.0370 ± 0.0301 |
| 64000/10 | 2.3046 ± 1.7140 | 1.6801 ± 0.6414 | 3.6836 ± 1.9491 | 0.2255 ± 0.0802 |
| 64000/100 | 14.2944 ± 2.9450 | 14.5547 ± 4.1186 | 3.3229 ± 0.7224 | 4.3851 ± 1.1917 |
| 512000/1 | 1.0889 ± 0.3636 | 1.2704 ± 0.5288 | 19.2216 ± 4.4009 | 0.2403 ± 0.1390 |
| 512000/10 | 12.2294 ± 4.0391 | 10.3849 ± 1.6810 | 17.9846 ± 2.7567 | 3.7742 ± 1.0887 |

## Unchanged Python biology controls

All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. Whole-corpus median batch times below are milliseconds, from nine within-process paired repetitions; they are not nine fresh-process estimates. Adaptive inner counts are normalized in the timings.

| Workload | Baseline exact | Optimized-core exact | Old-hash control, baseline / optimized | Inner repetitions, baseline / optimized |
|---|---:|---:|---:|---:|
| regions | 8.8537 | 8.8427 | 0.9487 / 0.9649 | 3 / 3 |
| gene_reports | 10.7233 | 10.9026 | 2.2239 / 2.2757 | 1 / 1 |
| panels | 222.3583 | 223.9911 | 145.2110 / 151.3969 | 1 / 1 |

Construction is excluded. These are timing controls, not clinical-efficacy or chromosome-scale measurements.

## Structural and operation microbenchmarks

Same fixtures, three alternating rounds × seven samples. Values are median [minimum, maximum] of three round medians, in microseconds. These round ranges are not confidence intervals; all individual samples and scratch-memory measurements remain in results/microbench.json.

Rust microbench independent roots are built separately within one Arena and use the recursive same-arena exact path. Their 1 MiB input repeats the byte sequence 0 through 255; the independent-target churn arm uses 1M edited ASCII bytes and separately allocated target arenas, with the frontier path for long cross-arena ranges. These comparison timings are not interchangeable or additive. The construction microbench includes Arena::new and arena destruction inside timing, whereas churn excludes initial setup.

The changed Rust construction/hash/edit paths improve in these fixtures: 64 KiB full arena construction/destruction 311.035→45.925 µs, split 2.334→0.466 µs, warm substring fingerprint 4.503→0.512 µs. The exact engines are unchanged, and their timing movements are mixed: independent 1 MiB LCP 60.304→57.467 µs, equality 58.656→68.446 µs, late-mismatch LCP 60.029→67.337 µs. Overlapping round ranges and unchanged engine source do not establish either a new exact-comparison optimization or timing equivalence. Compressed DAG/Repeat cases remain symbolic. Hash-only LCP rows are retained performance controls and remain collision-unsafe.

Python:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 13033.4500 [12906.5000, 13123.0000] | 12961.3500 [12867.2000, 12988.1500] |
| core/chunked_build_64KiB | 12932.1500 [12891.6500, 13034.9000] | 13133.9000 [12754.2500, 13265.9500] |
| core/concat_balanced_64KiB | 1.1550 [1.1481, 1.2072] | 1.1718 [1.1507, 1.2450] |
| core/split_chunked_64KiB | 113.5250 [112.8200, 120.5900] | 114.6850 [113.5050, 118.2450] |
| core/substr_warm_64KiB | 207.2350 [206.2450, 209.2900] | 211.2250 [203.4800, 214.6800] |
| core/root_hash | 0.0458 [0.0445, 0.0480] | 0.0456 [0.0452, 0.0484] |
| core/repeat_join_2pow32 | 603.2250 [583.7000, 610.5000] | 589.4250 [577.1250, 589.9000] |
| compare/65536/same_root/exact_lcp | 0.2350 [0.2100, 0.2950] | 0.2250 [0.2250, 0.2500] |
| compare/65536/same_root/exact_equal | 0.1100 [0.0850, 0.1250] | 0.0850 [0.0800, 0.1000] |
| compare/65536/same_root/hash_only_lcp | 1673.3333 [1598.4333, 1704.0333] | 1652.7333 [1620.6333, 1668.1667] |
| compare/65536/independent_equal/exact_lcp | 198.6150 [193.3150, 205.9800] | 195.2750 [190.2800, 199.9450] |
| compare/65536/independent_equal/exact_equal | 196.4650 [190.2200, 204.6250] | 192.4150 [190.4400, 194.8200] |
| compare/65536/independent_equal/hash_only_lcp | 1649.3667 [1597.6333, 1669.4000] | 1694.4333 [1679.1000, 1700.0667] |
| compare/65536/early_mismatch/exact_lcp | 6.2600 [6.1650, 6.3150] | 6.0850 [5.9300, 6.6000] |
| compare/65536/early_mismatch/exact_equal | 6.0600 [5.9350, 6.1600] | 6.1950 [5.9400, 6.2100] |
| compare/65536/early_mismatch/hash_only_lcp | 214.5667 [212.7000, 223.0000] | 218.5333 [217.9333, 229.3333] |
| compare/65536/late_mismatch/exact_lcp | 187.6250 [186.8850, 192.9950] | 197.4950 [192.0900, 202.3000] |
| compare/65536/late_mismatch/exact_equal | 190.9050 [187.9800, 202.1450] | 198.3300 [195.9150, 205.4600] |
| compare/65536/late_mismatch/hash_only_lcp | 1615.4333 [1606.9000, 1652.6000] | 1678.5667 [1653.1000, 1686.9667] |
| compare/65536/shared_late_mismatch/exact_lcp | 9.3500 [9.3150, 10.3850] | 9.7800 [9.7700, 9.9750] |
| compare/65536/shared_late_mismatch/exact_equal | 9.8950 [9.8900, 10.1300] | 9.6200 [9.4200, 9.7700] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 1570.7667 [1566.6333, 1644.2667] | 1563.0333 [1542.7000, 1602.4000] |
| compare/1048576/same_root/exact_lcp | 0.3500 [0.3350, 0.4050] | 0.3550 [0.3300, 0.4400] |
| compare/1048576/same_root/exact_equal | 0.1700 [0.1150, 0.2400] | 0.1700 [0.1600, 0.1900] |
| compare/1048576/same_root/hash_only_lcp | 1874.2000 [1795.3000, 2382.9000] | 1926.7667 [1877.0667, 1930.3333] |
| compare/1048576/independent_equal/exact_lcp | 4018.0950 [3963.1050, 4060.2200] | 3916.5200 [3841.4200, 3948.1700] |
| compare/1048576/independent_equal/exact_equal | 3934.5850 [3563.0550, 4426.2700] | 3986.5400 [3939.6100, 4217.2850] |
| compare/1048576/independent_equal/hash_only_lcp | 1969.6333 [1818.6333, 3169.6333] | 1922.8000 [1783.3000, 1958.0667] |
| compare/1048576/early_mismatch/exact_lcp | 9.6900 [9.2800, 16.3100] | 9.5100 [8.9600, 9.7400] |
| compare/1048576/early_mismatch/exact_equal | 9.5200 [9.5050, 9.7400] | 9.0850 [8.9350, 9.4850] |
| compare/1048576/early_mismatch/hash_only_lcp | 254.6333 [235.5000, 254.9667] | 250.2667 [246.2333, 254.5667] |
| compare/1048576/late_mismatch/exact_lcp | 3717.1050 [3716.7400, 4550.2650] | 3758.4200 [3665.2100, 4082.4600] |
| compare/1048576/late_mismatch/exact_equal | 4188.1550 [3953.9250, 4493.0350] | 3838.6850 [3626.7750, 4187.3100] |
| compare/1048576/late_mismatch/hash_only_lcp | 1776.3667 [1772.9000, 3317.9667] | 1895.1000 [1872.1667, 1922.6667] |
| compare/1048576/shared_late_mismatch/exact_lcp | 16.8250 [15.1300, 26.1000] | 14.6450 [14.5150, 15.0250] |
| compare/1048576/shared_late_mismatch/exact_equal | 14.8850 [14.0150, 25.7100] | 14.4250 [14.0250, 14.6500] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 1659.7333 [1647.1667, 1727.5333] | 1746.5667 [1728.1333, 1804.9000] |
| compare/repeat_2TiB/exact_lcp | 3.0065 [2.9205, 5.3940] | 2.8745 [2.7945, 3.0905] |
| compare/repeat_2TiB/exact_equal | 3.2000 [2.8975, 5.5205] | 3.1685 [2.9950, 3.2070] |
| compare/repeat_2TiB/hash_only_lcp | 3466.1500 [3141.8000, 3482.0500] | 3642.1000 [3248.9000, 3958.5000] |
| compare/self_concat_DAG_2TiB/exact_lcp | 59.5950 [57.4490, 61.9740] | 55.8620 [55.8290, 67.6000] |
| compare/self_concat_DAG_2TiB/exact_equal | 65.1400 [58.5280, 90.3820] | 54.8700 [52.1800, 63.7290] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 935.6500 [901.4000, 1577.1000] | 1003.9000 [959.1000, 1020.4000] |

Rust:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 311.0350 [295.4350, 312.3000] | 45.9250 [28.3750, 46.4750] |
| core/concat_balanced_64KiB | 0.0143 [0.0137, 0.0152] | 0.0152 [0.0136, 0.0156] |
| core/split_chunked_64KiB | 2.3340 [2.3260, 2.5340] | 0.4660 [0.4580, 0.4700] |
| core/substr_warm_64KiB | 4.5035 [4.4940, 4.5190] | 0.5124 [0.4803, 0.5169] |
| core/root_hash | 0.0005 [0.0005, 0.0006] | 0.0005 [0.0005, 0.0006] |
| core/repeat_join_2pow32 | 6.3580 [6.1880, 6.4130] | 6.2840 [6.2350, 6.9290] |
| compare/65536/same_root/exact_lcp | 0.0019 [0.0016, 0.0019] | 0.0019 [0.0019, 0.0019] |
| compare/65536/same_root/exact_equal | 0.0038 [0.0034, 0.0038] | 0.0037 [0.0037, 0.0039] |
| compare/65536/same_root/hash_only_lcp | 38.4070 [37.2300, 38.7530] | 4.9950 [4.9880, 5.0130] |
| compare/65536/independent_equal/exact_lcp | 3.5390 [3.3770, 3.5490] | 3.3510 [3.0500, 3.9320] |
| compare/65536/independent_equal/exact_equal | 3.7760 [3.2260, 3.9180] | 3.5090 [3.3430, 3.9350] |
| compare/65536/independent_equal/hash_only_lcp | 37.9320 [37.1890, 38.1390] | 4.9840 [4.9800, 4.9850] |
| compare/65536/early_mismatch/exact_lcp | 0.0537 [0.0512, 0.0571] | 0.0571 [0.0570, 0.0573] |
| compare/65536/early_mismatch/exact_equal | 0.0531 [0.0531, 0.0593] | 0.0535 [0.0532, 0.0595] |
| compare/65536/early_mismatch/hash_only_lcp | 4.9190 [4.9170, 5.0960] | 0.6520 [0.6490, 0.7130] |
| compare/65536/late_mismatch/exact_lcp | 3.4430 [3.1030, 3.4740] | 3.8020 [3.3830, 3.9650] |
| compare/65536/late_mismatch/exact_equal | 3.2210 [3.1470, 3.6110] | 3.1020 [3.0540, 3.7320] |
| compare/65536/late_mismatch/hash_only_lcp | 37.9810 [37.5780, 38.1860] | 4.5650 [4.5520, 4.8750] |
| compare/65536/shared_late_mismatch/exact_lcp | 0.0880 [0.0840, 0.0950] | 0.0930 [0.0850, 0.0950] |
| compare/65536/shared_late_mismatch/exact_equal | 0.0960 [0.0870, 0.0970] | 0.0960 [0.0880, 0.5010] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 35.1560 [35.0570, 40.6650] | 4.7100 [4.3880, 4.7520] |
| compare/1048576/same_root/exact_lcp | 0.0018 [0.0017, 0.0019] | 0.0018 [0.0017, 0.0019] |
| compare/1048576/same_root/exact_equal | 0.0037 [0.0035, 0.0039] | 0.0035 [0.0035, 0.0039] |
| compare/1048576/same_root/hash_only_lcp | 39.3930 [39.3830, 40.9070] | 6.3870 [6.2710, 6.5620] |
| compare/1048576/independent_equal/exact_lcp | 60.3040 [58.2270, 64.3100] | 57.4670 [57.3920, 60.8180] |
| compare/1048576/independent_equal/exact_equal | 58.6560 [58.3860, 59.8640] | 68.4460 [58.7790, 68.9780] |
| compare/1048576/independent_equal/hash_only_lcp | 40.2690 [39.5280, 40.7220] | 6.2370 [6.1770, 6.7190] |
| compare/1048576/early_mismatch/exact_lcp | 0.0837 [0.0785, 0.0843] | 0.0806 [0.0757, 0.0838] |
| compare/1048576/early_mismatch/exact_equal | 0.0775 [0.0771, 0.0804] | 0.0788 [0.0779, 0.0794] |
| compare/1048576/early_mismatch/hash_only_lcp | 5.3760 [5.1670, 5.4490] | 0.8850 [0.8700, 0.9670] |
| compare/1048576/late_mismatch/exact_lcp | 60.0290 [59.0970, 66.1360] | 67.3370 [61.2130, 68.7630] |
| compare/1048576/late_mismatch/exact_equal | 58.7460 [57.7360, 65.9320] | 66.8740 [57.1990, 67.2500] |
| compare/1048576/late_mismatch/hash_only_lcp | 39.3860 [38.8490, 41.0230] | 6.5390 [6.4900, 6.8810] |
| compare/1048576/shared_late_mismatch/exact_lcp | 0.1140 [0.1130, 0.1250] | 0.1260 [0.1250, 0.1270] |
| compare/1048576/shared_late_mismatch/exact_equal | 0.1180 [0.1150, 0.1280] | 0.1280 [0.1270, 0.1280] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 36.9030 [36.7160, 38.7160] | 6.1850 [6.1590, 6.4340] |
| compare/repeat_2TiB/exact_lcp | 0.0813 [0.0801, 0.0888] | 0.0810 [0.0803, 0.0817] |
| compare/repeat_2TiB/exact_equal | 0.0820 [0.0819, 0.0880] | 0.0837 [0.0814, 0.0841] |
| compare/repeat_2TiB/hash_only_lcp | 38.5370 [38.3960, 38.9050] | 36.8150 [35.9760, 36.9810] |
| compare/self_concat_DAG_2TiB/exact_lcp | 0.2810 [0.2520, 0.2870] | 0.2510 [0.2510, 0.2800] |
| compare/self_concat_DAG_2TiB/exact_equal | 0.3400 [0.3050, 0.3410] | 0.2810 [0.2730, 0.2970] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 23.1460 [22.9630, 24.0750] | 22.6670 [22.4590, 22.8660] |

Rust layout remains NodeInner 48 bytes / Arena inline 2136 bytes, excluding heap allocations, with 255 nodes for 64 KiB. Source compatibility does not imply identical operation timing. Compressed-path, construction/hash/edit and unchanged Python results must be assessed with their recorded round variability.

## Provenance and reproduction

Every result binds its immutable core, driver, consumer helpers, input hashes, workload and runtime context. Resume and aggregation reject mismatched identities. Rust results also reference compiled-binary hashes and locked dependency manifests. Python worker imports select and check the intended frozen package. Previous reports/results are retained unchanged.

[preserved_results.json](preserved_results.json) verifies 243 prior benchmark source/artifact files remain byte-identical.

[final_source_verification.json](final_source_verification.json) verifies the final working source against the measured optimized archive, including any documentation-only edits.

[final_audit.json](final_audit.json) records the complete matrix, context identities, verified resume and alternate-output aggregation checks. The final revised paper compiled successfully; the completion record is in VALIDATION.md and the paper LOGBOOK.

Use a fresh output directory for a full reproduction:

```powershell
python benchmarks/exact_optimization_round3/run_matrix.py --output-root benchmarks/exact_optimization_round3/reproduction_results
python benchmarks/exact_optimization_round3/provenance.py --results-root benchmarks/exact_optimization_round3/reproduction_results
python benchmarks/exact_optimization_round3/summarize.py --results-root benchmarks/exact_optimization_round3/reproduction_results --output-dir benchmarks/exact_optimization_round3/reproduction_report
```

The runner restores saved core archives, verifies corpus-token caches and can regenerate missing caches with the offline GPT-2 tokenizer. The original local paper/biology datasets remain required. Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies.
