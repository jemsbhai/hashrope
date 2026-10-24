# Comparison regression investigation: round four

The independent-equality slowdown recurs in the final short, full-protocol rerun, although the longer controlled and retained-original-binary follow-ups did not establish a stable primary old/new penalty. The evidence is protocol-sensitive and its mechanism remains unresolved. No round-four production change was selected without an isolated cause and validated remedy. [controlled_findings.md](controlled_findings.md) records the estimates, exploratory intervals and limits. The [round-one report](../exact_optimization/README.md), [round-two report](../exact_optimization_round2/README.md) and [round-three report](../exact_optimization_round3/README.md), their raw results and original source archives remain preserved.

The frozen [measurement plan](measurement_plan.json) compares **pre-round-three (completed round two)** with **unchanged completed round three**. The result keys baseline/optimized identify these two measured cohorts; baseline.zip separately preserves the completed-round-three working state. No round-four production optimization was selected. This is a fresh pre-round-three versus unchanged-current retest, not an A/A replacement for that comparison. Any identical-source or identical-binary controls remain separately labeled diagnostics.

Same-arena independent 1 MiB equality measures 54.926→65.895 µs; late-mismatch LCP measures 57.184→60.059 µs. These are medians of three round medians; full ranges and all controls are reported below. The 2 TiB self-concatenation DAG's LCP is 0.280→0.280 µs.

At 10K bytes, independent-target edit-and-compare changes from 7.428 to 4.468 µs (1.66× speedup), and shared history from 6.643 to 3.703 µs (1.79×), with 9/9 and 9/9 paired improvements respectively. Exact comparisons plus edits still take 3.82× and 4.06× Ropey's time at 10K bytes. In the repeated 1M-byte workload, independent-target exact time falls 6.0% and is 1.35× Ropey's time. Shared-history time falls 37.6% and takes 26.8% of Ropey's time in this fixture.

## Method and scope

Both measured core versions are immutable archives with file manifests. The same consumer runners, corpus inputs, builders, sizes, seed/invocation order, match/mismatch mix, warmups and fixed repetition counts are retained. Nine-process tables report mean ± sample SD across process medians, in milliseconds. Real-pair and R2 queries are first averaged within each process. Speedup ratios use the reported means. Biology retains nine paired repetitions within one process per core, with recorded adaptive batching; its inner counts may differ across cores. Microbenchmarks report the median of three round medians and their minimum/maximum; these ranges are not confidence intervals.

Python is unchanged between the two core snapshots and is rerun as a control. Differences in Python consumer timings therefore do not measure a Rust implementation effect. Preparation remains outside query timing. CPU-heavy team work is paused for final timing, on the same Windows 11 / Intel Core i9-14900HX workstation, without CPU affinity or frequency locking. No GPU serving, TTFT, NVML energy or chromosome-scale experiment is rerun.

[VALIDATION.md](VALIDATION.md) records the correctness checks and explicit experimental exclusions.

[unchanged_code.json](unchanged_code.json) verifies 6 complete modules, 34 Python definitions and storage/default declarations against the preserved completed-round-three working state. It audits 98 existing Rust production functions against the preserved completed-round-three working state, recording 0 selected changes and 0 additions explicitly; all other production function bodies and existing signatures remain unchanged. The audit also compares production text outside function bodies, including inline/code-generation attributes. It is distinct from the measured reference choice. Source checks establish change scope; compatibility tests and the retained arithmetic review support fingerprint semantics.

Exact comparisons remain independent of fingerprints. Existing accepted custom profiles can produce different stored fingerprints for identical bytes in different tree shapes: with modulus 17 and base 2, bytes [0, 13] hash to 16 as one leaf and 0 as a concatenation under the legacy reduction. Rejecting equality on a fingerprint mismatch would therefore break existing exact semantics. This round preserves those profiles and fingerprint values.

## Why no further production change was selected

The original round-three observations remain in their original report: independent 1 MiB equality rose from 58.656 to 68.446 µs, and late-mismatch LCP from 60.029 to 67.337 µs, using medians of three round medians. [controlled_findings.md](controlled_findings.md) records two preregistered follow-ups, their raw evidence and limitations; [REVIEW.md](REVIEW.md) records the independent correctness and measurement audit.

The controlled diagnostic used 36 fresh processes and 3,456 timed blocks, with counterbalanced core, fixture, API and block order, swapped crate slots, identical-current A/A controls, and separate unpinned and CPU-0-pinned arms. Across old/new slot and affinity combinations, independent-equality current/reference geometric ratios were 0.9868–1.0110, and late-LCP ratios were 0.9882–1.0374. The largest primary-cell exploratory 95% interval upper endpoint was 1.0976 for pinned late LCP. A separate pinned/swapped late-equality cell, distinct from late LCP, was 3.15% slower; its normal-slot counterpart did not show the same effect. These checks do not prove equivalence or exclude smaller context-dependent costs.

The retained-original-binary replication verified the exact two earlier executable hashes, then ran a fixed 12 fresh unpinned pairs with six old-first and six current-first pairs, retaining all 7,056 raw samples. Independent-equality current/reference paired geometric ratio was 0.9959 (exploratory 95% interval 0.9458–1.0486); late LCP was 0.9914 (0.9326–1.0540). These paired geometric estimates differ from the original three-round median statistic. The original outer Python microbenchmarks were not interleaved in this diagnostic, which is one reason for the separate full-protocol matrix below. See the [controlled data](results/controlled_initial/results.json) and [original-binary replication](results/original_replication/results.json).

Comparator source was unchanged, and separately emitted normalized instruction text matched for seven comparator routines; the assembly audit retains the dispatcher section-ordinal difference. This does not establish final linked-binary or relocation-target identity. Those two follow-ups did not reproduce the original primary 12–17% penalties as a stable old/new effect, but the full-protocol equality recurrence leaves the cause unresolved, so no validated remedy was selected. Heap placement, code layout, preceding work, migration and frequency effects remain possible influences, without an established causal explanation. Intervals are exploratory, per cell and not multiplicity adjusted.

## Rust edit plus exact target comparison

Each timed batch performs 50 seeded churn edits and compares the edited rope with a prebuilt target: 50% equal and 50% different at the last byte. Every boolean result is checked. ASCII fixtures make byte and Ropey character offsets coincide. Target preparation/destruction is excluded. Ropey is pinned to 1.6.1. Independent targets and shared ancestry remain separate regimes. Both fingerprint-candidate arms use Hashrope from their measured core version, including Ropey's rehash-candidate arm. Ropey exact uses its unchanged equality implementation. The selected pre-round-three reference includes the earlier scalar hash kernel; completed-round-three and later cores include the reviewed grouped kernel, so the candidate arms are not unchanged controls across that historical boundary.

| Bytes / history | Baseline exact | Optimized exact | Baseline / optimized | Paired improvements | Ropey, optimized cohort |
|---|---:|---:|---:|---:|---:|
| 10,000 / independent | 0.0074 ± 0.0006 | 0.0045 ± 0.0002 | 1.66× | 9/9 | 0.0012 ± 0.0002 |
| 10,000 / shared | 0.0066 ± 0.0005 | 0.0037 ± 0.0001 | 1.79× | 9/9 | 0.0009 ± 0.0001 |
| 100,000 / independent | 0.0236 ± 0.0011 | 0.0192 ± 0.0019 | 1.23× | 9/9 | 0.0112 ± 0.0016 |
| 100,000 / shared | 0.0120 ± 0.0007 | 0.0074 ± 0.0006 | 1.63× | 9/9 | 0.0054 ± 0.0006 |
| 1,000,000 / independent | 0.1534 ± 0.0147 | 0.1442 ± 0.0073 | 1.06× | 7/9 | 0.1067 ± 0.0080 |
| 1,000,000 / shared | 0.0147 ± 0.0008 | 0.0092 ± 0.0005 | 1.60× | 9/9 | 0.0342 ± 0.0027 |

At 1M bytes (independent targets), optimized exact Hashrope/Ropey is 1.351. The separate fingerprint-candidate control changes from 0.01145 to 0.00512 ms/cycle.

At 1M bytes (shared ancestry), optimized exact Hashrope/Ropey is 0.268. The separate fingerprint-candidate control changes from 0.00994 to 0.00309 ms/cycle.

The unchanged Ropey exact control also varies between cohorts: at 1M bytes with shared history it measures 0.0440 ± 0.0107→0.0342 ± 0.0027 ms. This control variation limits causal attribution of timing movements and makes the cohort-specific comparison important.


For 10K independent targets, the edit-plus-fingerprint-candidate arm is 0.0048 ± 0.0003→0.0019 ± 0.0002 ms, versus 0.0012 ± 0.0002 ms for Ropey edit plus exact comparison. The retained round-three decomposition separately measured edit and metadata costs.

These fixture-specific ratios and paired win counts do not establish statistical significance, a universal equality advantage or identical timings on every path.

Fingerprint-candidate arms from the same runs (milliseconds per cycle); these are distinct from exact equality:

| Bytes / history | Hashrope candidate, baseline | Hashrope candidate, optimized | Ropey rehash candidate, baseline | Ropey rehash candidate, optimized |
|---|---:|---:|---:|---:|
| 10,000 / independent | 0.0048 ± 0.0003 | 0.0019 ± 0.0002 | 0.0437 ± 0.0018 | 0.0049 ± 0.0003 |
| 10,000 / shared | 0.0045 ± 0.0002 | 0.0018 ± 0.0001 | 0.0435 ± 0.0017 | 0.0049 ± 0.0006 |
| 100,000 / independent | 0.0084 ± 0.0005 | 0.0034 ± 0.0002 | 0.4276 ± 0.0103 | 0.0402 ± 0.0020 |
| 100,000 / shared | 0.0083 ± 0.0009 | 0.0036 ± 0.0004 | 0.4249 ± 0.0074 | 0.0407 ± 0.0036 |
| 1,000,000 / independent | 0.0115 ± 0.0036 | 0.0051 ± 0.0003 | 4.5387 ± 0.3766 | 0.4055 ± 0.0125 |
| 1,000,000 / shared | 0.0099 ± 0.0037 | 0.0031 ± 0.0003 | 5.5574 ± 0.5083 | 0.4162 ± 0.0218 |

## Unchanged Python paper controls

EXP-005, mismatch halfway through independently built corpus-derived ropes; all other fractions remain in the raw results:

| Bytes | Baseline core | Optimized core |
|---:|---:|---:|
| 64,000 | 0.0170 ± 0.0040 | 0.0157 ± 0.0008 |
| 256,000 | 0.0678 ± 0.0205 | 0.0589 ± 0.0136 |
| 1,000,000 | 0.2127 ± 0.0131 | 0.2314 ± 0.0632 |
| 2,000,000 | 0.4694 ± 0.1377 | 0.4929 ± 0.0986 |
| 4,000,000 | 0.8458 ± 0.0564 | 0.8404 ± 0.0406 |
| 8,000,000 | 1.9239 ± 0.4077 | 2.0507 ± 0.5292 |

Log-log slopes fitted to the six cohort means are 0.963 and 0.998. The historical ≤0.3 gate fails for the baseline and fails for the optimized cohort. This gate remains descriptive; no logarithmic exact-LCP bound follows.

Controlled EXP-017:

| Shared tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0045 ± 0.0010 | 0.0050 ± 0.0019 | 0.0393 ± 0.0100 | 0.0029 ± 0.0014 |
| 4,000 | 0.0110 ± 0.0042 | 0.0095 ± 0.0030 | 0.1475 ± 0.0455 | 0.0039 ± 0.0018 |
| 16,000 | 0.0479 ± 0.0372 | 0.0432 ± 0.0209 | 0.6136 ± 0.2251 | 0.0064 ± 0.0025 |
| 64,000 | 0.1320 ± 0.0367 | 0.1100 ± 0.0259 | 3.3188 ± 1.0617 | 0.0273 ± 0.0154 |
| 128,000 | 0.2070 ± 0.0137 | 0.2545 ± 0.0631 | 6.5265 ± 2.6186 | 0.0867 ± 0.0888 |
| 256,000 | 0.5258 ± 0.1258 | 0.4990 ± 0.2148 | 10.0934 ± 3.0933 | 0.1459 ± 0.1266 |
| 512,000 | 0.9807 ± 0.2487 | 0.9417 ± 0.2090 | 18.6884 ± 1.5118 | 0.2570 ± 0.1378 |
| 1,000,000 | 2.0488 ± 0.3958 | 2.1717 ± 0.5821 | 35.0915 ± 3.4419 | 0.5915 ± 0.2826 |
| 2,000,000 | 4.3893 ± 0.6572 | 4.1296 ± 0.5370 | 75.3032 ± 8.9531 | 1.9083 ± 0.2610 |

Optimized-cohort exact wins 9/9 paired comparisons against Radix at 2M tokens. The literal original crossover-grid condition passes; stable mean ± one sample SD band separation starts at 1000 tokens. These are dispersion bands, not confidence intervals or a significance test.
That is the smallest measured point, so the below-threshold condition is vacuous; an interior crossover is not observed.

EXP-015 R2 CPU identification:

| Cached tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | Excluded query build, optimized cohort |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0015 ± 0.0003 | 0.0016 ± 0.0008 | 0.0421 ± 0.0094 | 1.8599 ± 0.2641 |
| 16,000 | 0.0288 ± 0.0043 | 0.0317 ± 0.0084 | 0.6481 ± 0.0920 | 17.3144 ± 1.9619 |
| 64,000 | 0.1206 ± 0.0121 | 0.1247 ± 0.0269 | 3.6316 ± 0.9295 | 60.5874 ± 5.7008 |
| 256,000 | 0.5373 ± 0.0611 | 0.5141 ± 0.1387 | 10.6803 ± 1.4550 | 263.1952 ± 13.7981 |
| 1,000,000 | 2.1766 ± 0.6161 | 1.7993 ± 0.1694 | 38.1418 ± 4.1017 | 999.2188 ± 35.3970 |

An integration starting with flat tokens must account for query construction; the prebuilt-query figures exclude it.

## Unchanged Python real-pair and candidate-set controls

Original 200 ShareGPT/LMSYS pairs per dataset/seed. Inner timing order remains fixed: radix, exact Hashrope, NumPy. K candidates retain one indexed radix lookup versus K pairwise comparisons with independently built Hashrope roots.

Conversational pairs:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---|---:|---:|---:|---:|
| sharegpt_sample.jsonl | 0.0039 ± 0.0004 | 0.0039 ± 0.0005 | 0.0340 ± 0.0039 | 0.0022 ± 0.0004 |
| lmsys_sample.jsonl | 0.0032 ± 0.0005 | 0.0031 ± 0.0005 | 0.0214 ± 0.0041 | 0.0020 ± 0.0004 |

Shared tokens / K:

| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |
|---|---:|---:|---:|---:|
| 64000/1 | 0.2809 ± 0.3044 | 0.2354 ± 0.2299 | 2.8210 ± 0.7074 | 0.0246 ± 0.0118 |
| 64000/10 | 1.5400 ± 0.5786 | 1.4360 ± 0.3892 | 3.1041 ± 1.6970 | 0.1807 ± 0.0522 |
| 64000/100 | 12.7485 ± 1.3204 | 13.8696 ± 2.4057 | 2.9453 ± 0.6905 | 3.9263 ± 0.9100 |
| 512000/1 | 1.1428 ± 0.4853 | 0.9910 ± 0.2353 | 18.5511 ± 3.1610 | 0.2071 ± 0.0533 |
| 512000/10 | 10.3699 ± 1.3409 | 11.6710 ± 3.6026 | 18.2668 ± 3.8485 | 3.5338 ± 1.0850 |

## Unchanged Python biology controls

All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. Whole-corpus median batch times below are milliseconds, from nine within-process paired repetitions; they are not nine fresh-process estimates. Adaptive inner counts are normalized in the timings.

| Workload | Baseline exact | Optimized-core exact | Old-hash control, baseline / optimized | Inner repetitions, baseline / optimized |
|---|---:|---:|---:|---:|
| regions | 7.6802 | 7.8064 | 0.7583 / 0.7505 | 4 / 4 |
| gene_reports | 9.4610 | 9.7551 | 1.7666 / 1.8173 | 2 / 2 |
| panels | 213.9133 | 194.1372 | 143.5535 / 141.9041 | 1 / 1 |

Construction is excluded. These are timing controls, not clinical-efficacy or chromosome-scale measurements.

## Structural and operation microbenchmarks

Same fixtures, three alternating rounds × seven samples. Values are median [minimum, maximum] of three round medians, in microseconds. These round ranges are not confidence intervals; all individual samples and scratch-memory measurements remain in results/microbench.json.

Rust microbench independent roots are built separately within one Arena and use the recursive same-arena exact path. Their 1 MiB input repeats the byte sequence 0 through 255; the independent-target churn arm uses 1M edited ASCII bytes and separately allocated target arenas, with the frontier path for long cross-arena ranges. These comparison timings are not interchangeable or additive. The construction microbench includes Arena::new and arena destruction inside timing, whereas churn excludes initial setup.

Rust operation medians in these fixtures: 64 KiB full arena construction/destruction 288.730→47.020 µs, split 2.553→0.460 µs, warm substring fingerprint 4.491→0.492 µs. The comparison-only medians are: independent 1 MiB LCP 61.702→58.778 µs, equality 54.926→65.895 µs, late-mismatch LCP 57.184→60.059 µs. Round ranges alone do not establish statistical significance or timing equivalence, and source equality alone does not isolate the cause of runtime movements. Compressed DAG/Repeat cases remain symbolic. Hash-only LCP rows are retained performance controls and remain collision-unsafe.

The final short-protocol equality median is 20.0% higher, and late-LCP median 5.0% higher. The equality penalty's direction therefore recurs in the complete outer protocol, rather than disappearing from all replications. This differs from the paired geometric diagnostic estimates above; those estimates do not invalidate these observations. The full matrix interleaves the original Python/Rust rounds and retains the short fixed sequence. It also rebuilds the executables; their hashes differ from the retained-original-binary pair, so both outer execution history and linked build context differ. The unchanged source/harness therefore does not isolate the outer sequence as the cause. No controlled experiment yet isolates the responsible mechanism, so a universal slowdown, timing equivalence, or a dismissal as noise would each overstate the evidence.

Other slower cells remain visible: 64 KiB late-mismatch LCP is 3.082→3.695 µs (+19.9%), and late-mismatch equality 3.083→3.979 µs (+29.1%), with nonoverlapping observed three-round ranges. These secondary short-protocol observations also require qualification; the 1 MiB equality cell is not the only recorded increase.

Final Rust micro executable SHA-256 values: reference `f4d7546d041bb450540a9136ea9c614ee0cb8e7aac7c60d0a1708df13bc765d2`, current `ef2bf3295560acbbc7331951fc0c57faf0170d294c0d55307a0f40241280fca2`. The earlier retained executable identities are recorded in the [replication data](results/original_replication/results.json).

Python:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 12658.2500 [12605.5500, 13639.8500] | 12607.2500 [12375.1500, 12712.6000] |
| core/chunked_build_64KiB | 13098.3000 [12608.2500, 14123.4000] | 12899.4500 [12580.0000, 13045.8000] |
| core/concat_balanced_64KiB | 1.1565 [1.1204, 1.2602] | 1.1307 [1.1160, 1.2450] |
| core/split_chunked_64KiB | 113.8550 [111.6250, 122.8150] | 114.5450 [110.1000, 115.0250] |
| core/substr_warm_64KiB | 214.0400 [198.1950, 220.8100] | 204.7000 [197.6950, 206.3750] |
| core/root_hash | 0.0473 [0.0462, 0.0478] | 0.0465 [0.0447, 0.0468] |
| core/repeat_join_2pow32 | 578.4500 [575.3500, 586.1250] | 595.2500 [582.3000, 649.1750] |
| compare/65536/same_root/exact_lcp | 0.2200 [0.1900, 0.2250] | 0.2150 [0.1800, 0.2550] |
| compare/65536/same_root/exact_equal | 0.0850 [0.0850, 0.0850] | 0.0850 [0.0750, 0.0900] |
| compare/65536/same_root/hash_only_lcp | 1663.3333 [1595.7333, 1727.4667] | 1598.9333 [1584.4333, 1681.2000] |
| compare/65536/independent_equal/exact_lcp | 190.7750 [182.9300, 199.4300] | 184.3350 [178.4100, 184.6850] |
| compare/65536/independent_equal/exact_equal | 192.9850 [191.5000, 193.5250] | 183.8900 [179.2200, 185.6550] |
| compare/65536/independent_equal/hash_only_lcp | 1654.9000 [1614.8000, 1676.8333] | 1633.5000 [1564.8333, 1667.0000] |
| compare/65536/early_mismatch/exact_lcp | 5.9700 [5.9200, 6.5400] | 5.9600 [5.9200, 5.9950] |
| compare/65536/early_mismatch/exact_equal | 5.9200 [5.8850, 6.8150] | 6.3600 [5.9900, 6.5500] |
| compare/65536/early_mismatch/hash_only_lcp | 213.0000 [212.5667, 222.9333] | 212.6667 [211.7667, 231.1000] |
| compare/65536/late_mismatch/exact_lcp | 189.3050 [184.6400, 193.7850] | 190.1150 [181.0100, 196.3600] |
| compare/65536/late_mismatch/exact_equal | 192.5650 [189.5400, 195.7100] | 191.8500 [180.2500, 196.0150] |
| compare/65536/late_mismatch/hash_only_lcp | 1627.4667 [1595.3333, 1661.4333] | 1661.4333 [1573.0000, 1716.5333] |
| compare/65536/shared_late_mismatch/exact_lcp | 9.4950 [9.3750, 10.4450] | 10.3600 [9.4100, 10.3800] |
| compare/65536/shared_late_mismatch/exact_equal | 9.5600 [9.4500, 10.5050] | 9.6450 [9.3450, 10.4650] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 1521.0000 [1504.5333, 1578.0667] | 1548.0667 [1478.8667, 1609.0000] |
| compare/1048576/same_root/exact_lcp | 0.2700 [0.2650, 0.3100] | 0.2400 [0.2100, 0.4050] |
| compare/1048576/same_root/exact_equal | 0.1150 [0.0850, 0.1300] | 0.1000 [0.0900, 0.1850] |
| compare/1048576/same_root/hash_only_lcp | 1758.8333 [1738.1667, 1932.3667] | 1717.3667 [1710.5667, 1721.9667] |
| compare/1048576/independent_equal/exact_lcp | 3753.6850 [3525.5900, 4288.4150] | 3601.4250 [3284.1650, 4121.1450] |
| compare/1048576/independent_equal/exact_equal | 3667.0300 [3569.6300, 4128.6700] | 3371.7050 [3265.7900, 3932.7600] |
| compare/1048576/independent_equal/hash_only_lcp | 2015.0667 [1717.1667, 2828.0333] | 1965.2333 [1696.5333, 3188.3000] |
| compare/1048576/early_mismatch/exact_lcp | 9.6700 [9.3200, 9.6750] | 8.8500 [8.6350, 8.8900] |
| compare/1048576/early_mismatch/exact_equal | 9.6200 [9.2100, 16.4650] | 8.9850 [8.3950, 16.7500] |
| compare/1048576/early_mismatch/hash_only_lcp | 250.2000 [243.4667, 291.7333] | 242.0000 [228.2667, 248.8333] |
| compare/1048576/late_mismatch/exact_lcp | 3709.5850 [3606.5550, 4040.5700] | 3508.4850 [3265.2800, 3702.4050] |
| compare/1048576/late_mismatch/exact_equal | 3949.3400 [3607.8200, 4386.8200] | 3417.9950 [3322.4400, 3662.1750] |
| compare/1048576/late_mismatch/hash_only_lcp | 3221.3000 [1786.4000, 3238.7333] | 1721.0667 [1657.6667, 1844.1667] |
| compare/1048576/shared_late_mismatch/exact_lcp | 25.0200 [13.7400, 26.4150] | 13.7650 [13.6450, 13.8700] |
| compare/1048576/shared_late_mismatch/exact_equal | 24.8850 [13.4050, 27.7450] | 14.2000 [13.6850, 14.4850] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 1783.4000 [1650.8667, 1795.8000] | 1700.0333 [1640.2667, 2793.1333] |
| compare/repeat_2TiB/exact_lcp | 2.9635 [2.6695, 3.2800] | 2.7220 [2.6875, 2.9935] |
| compare/repeat_2TiB/exact_equal | 3.1170 [2.8500, 5.5120] | 2.7595 [2.7260, 3.1920] |
| compare/repeat_2TiB/hash_only_lcp | 3082.9000 [2864.9500, 4550.1000] | 3132.8000 [2916.0500, 3368.1500] |
| compare/self_concat_DAG_2TiB/exact_lcp | 51.2480 [49.8170, 62.7260] | 53.5610 [52.5340, 61.6070] |
| compare/self_concat_DAG_2TiB/exact_equal | 51.6170 [49.7580, 52.6840] | 51.4480 [49.8320, 55.8090] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 879.5000 [877.1500, 995.3500] | 900.1500 [889.1500, 902.5500] |

Rust:

| Fixture | Baseline | Optimized |
|---|---:|---:|
| core/from_bytes_64KiB | 288.7300 [288.2750, 309.4500] | 47.0200 [42.6400, 51.9050] |
| core/concat_balanced_64KiB | 0.0150 [0.0133, 0.0161] | 0.0150 [0.0135, 0.0173] |
| core/split_chunked_64KiB | 2.5530 [2.3860, 2.5760] | 0.4600 [0.4220, 0.5070] |
| core/substr_warm_64KiB | 4.4907 [4.3618, 4.5901] | 0.4923 [0.4803, 0.5147] |
| core/root_hash | 0.0005 [0.0005, 0.0005] | 0.0005 [0.0005, 0.0006] |
| core/repeat_join_2pow32 | 6.3370 [5.6230, 6.4850] | 6.0240 [5.7250, 6.5040] |
| compare/65536/same_root/exact_lcp | 0.0016 [0.0016, 0.0018] | 0.0017 [0.0016, 0.0019] |
| compare/65536/same_root/exact_equal | 0.0035 [0.0033, 0.0035] | 0.0033 [0.0033, 0.0037] |
| compare/65536/same_root/hash_only_lcp | 37.7270 [36.9830, 38.0030] | 4.7420 [4.5560, 4.9370] |
| compare/65536/independent_equal/exact_lcp | 3.5360 [3.0370, 3.5730] | 3.5440 [3.2840, 3.6180] |
| compare/65536/independent_equal/exact_equal | 3.5500 [2.9960, 3.5520] | 3.5240 [3.4180, 3.6120] |
| compare/65536/independent_equal/hash_only_lcp | 36.8590 [36.1230, 39.8100] | 4.5640 [4.5480, 4.7490] |
| compare/65536/early_mismatch/exact_lcp | 0.0512 [0.0512, 0.0512] | 0.0521 [0.0521, 0.0572] |
| compare/65536/early_mismatch/exact_equal | 0.0531 [0.0530, 0.0590] | 0.0544 [0.0543, 0.0594] |
| compare/65536/early_mismatch/hash_only_lcp | 4.9650 [4.8800, 5.0290] | 0.7110 [0.6510, 0.7120] |
| compare/65536/late_mismatch/exact_lcp | 3.0820 [3.0690, 3.2440] | 3.6950 [3.6120, 3.9610] |
| compare/65536/late_mismatch/exact_equal | 3.0830 [3.0780, 3.0890] | 3.9790 [3.9630, 4.0130] |
| compare/65536/late_mismatch/hash_only_lcp | 37.4040 [36.0300, 38.1610] | 4.7120 [4.5560, 4.7670] |
| compare/65536/shared_late_mismatch/exact_lcp | 0.0950 [0.0840, 0.0950] | 0.0940 [0.0870, 0.0970] |
| compare/65536/shared_late_mismatch/exact_equal | 0.0960 [0.0860, 0.0960] | 0.0900 [0.0880, 0.0960] |
| compare/65536/shared_late_mismatch/hash_only_lcp | 35.2000 [33.9650, 35.7130] | 4.3650 [4.3630, 4.6640] |
| compare/1048576/same_root/exact_lcp | 0.0017 [0.0016, 0.0018] | 0.0018 [0.0016, 0.0018] |
| compare/1048576/same_root/exact_equal | 0.0035 [0.0035, 0.0035] | 0.0035 [0.0033, 0.0039] |
| compare/1048576/same_root/hash_only_lcp | 38.7880 [37.7270, 40.3990] | 6.2900 [6.2680, 6.8640] |
| compare/1048576/independent_equal/exact_lcp | 61.7020 [54.7610, 69.5240] | 58.7780 [57.9440, 66.0950] |
| compare/1048576/independent_equal/exact_equal | 54.9260 [54.5430, 58.8610] | 65.8950 [55.9100, 65.9900] |
| compare/1048576/independent_equal/hash_only_lcp | 39.7160 [37.8530, 40.2930] | 6.4460 [6.3700, 6.5580] |
| compare/1048576/early_mismatch/exact_lcp | 0.0771 [0.0752, 0.0838] | 0.0833 [0.0784, 0.0843] |
| compare/1048576/early_mismatch/exact_equal | 0.0775 [0.0772, 0.0862] | 0.0805 [0.0782, 0.0866] |
| compare/1048576/early_mismatch/hash_only_lcp | 5.4570 [5.1680, 5.5700] | 0.9700 [0.9670, 1.0490] |
| compare/1048576/late_mismatch/exact_lcp | 57.1840 [54.6320, 57.4490] | 60.0590 [57.4260, 65.8090] |
| compare/1048576/late_mismatch/exact_equal | 56.6840 [55.6120, 58.0690] | 59.7710 [57.3960, 63.8940] |
| compare/1048576/late_mismatch/hash_only_lcp | 38.5550 [38.5170, 39.8580] | 6.2590 [6.1850, 6.8430] |
| compare/1048576/shared_late_mismatch/exact_lcp | 0.1130 [0.1130, 0.1130] | 0.1260 [0.1150, 0.1260] |
| compare/1048576/shared_late_mismatch/exact_equal | 0.1150 [0.1150, 0.1150] | 0.1280 [0.1170, 0.1280] |
| compare/1048576/shared_late_mismatch/hash_only_lcp | 37.3790 [36.0630, 37.8690] | 6.2870 [6.0590, 6.4080] |
| compare/repeat_2TiB/exact_lcp | 0.0790 [0.0716, 0.0839] | 0.0820 [0.0729, 0.0849] |
| compare/repeat_2TiB/exact_equal | 0.0832 [0.0739, 0.0873] | 0.0777 [0.0767, 0.0884] |
| compare/repeat_2TiB/hash_only_lcp | 35.7440 [35.2620, 35.7490] | 37.2000 [35.1230, 38.0440] |
| compare/self_concat_DAG_2TiB/exact_lcp | 0.2800 [0.2510, 0.2800] | 0.2800 [0.2560, 0.2830] |
| compare/self_concat_DAG_2TiB/exact_equal | 0.3070 [0.3060, 0.3330] | 0.3130 [0.2720, 0.3140] |
| compare/self_concat_DAG_2TiB/hash_only_lcp | 23.4690 [22.0470, 23.7500] | 23.4620 [22.7700, 23.6440] |

Rust inline layout remains NodeInner 48 bytes / Arena 2136 bytes, excluding heap allocations, with 255 nodes for 64 KiB. Source compatibility does not imply identical operation timing. Compressed-path, construction/hash/edit and unchanged Python results must be assessed with their recorded round variability.

## Provenance and reproduction

Every result binds its immutable core, driver, consumer helpers, input hashes, workload and runtime context. Resume and aggregation reject mismatched identities. Rust results also reference compiled-binary hashes and locked dependency manifests. Python worker imports select and check the intended frozen package. Previous reports/results are retained unchanged.

[preserved_results.json](preserved_results.json) verifies 376 prior benchmark source/artifact files remain byte-identical.

[final_source_verification.json](final_source_verification.json) verifies measured production and integration-test source against the optimized archive; subsequent report and paper edits are outside the timed source scope.

[final_audit.json](final_audit.json) records the complete matrix, context identities, verified resume and alternate-output aggregation checks. The final validation and paper status are recorded in VALIDATION.md and the paper LOGBOOK.

Use a fresh output directory for a full reproduction:

```powershell
python benchmarks/exact_optimization_round4/run_matrix.py --output-root benchmarks/exact_optimization_round4/reproduction_results
python benchmarks/exact_optimization_round4/provenance.py --results-root benchmarks/exact_optimization_round4/reproduction_results
python benchmarks/exact_optimization_round4/summarize.py --results-root benchmarks/exact_optimization_round4/reproduction_results --output-dir benchmarks/exact_optimization_round4/reproduction_report
```

The runner restores saved core archives, verifies corpus-token caches and can regenerate missing caches with the offline GPT-2 tokenizer. The original local paper/biology datasets remain required. Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies.
