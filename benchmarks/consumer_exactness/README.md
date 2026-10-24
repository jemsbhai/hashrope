# Exact consumer reruns: 2026-10-23

This correction preserves fingerprints and evaluates exact consumer decisions. All historical paper result files remain untouched.

## Protocol and scope

Paper CPU measurements use three corpus seeds × three fresh interpreter processes, five warmed query repetitions per cell, and mean ± sample standard deviation across the nine process medians. Old and exact arms run on the same independently built ropes; the old arm is loaded from paper commit `c9d8902e0f443f1a5ba12da78f087d5dbd2cad68`. Both use current Hashrope, isolating consumer semantics. No shared-root shortcut is supplied. The original 4096-byte leaf builder, controlled divergent-tail construction, vendored radix implementation, and NumPy int64 comparator are retained. Old/exact order reverses on alternating invocations. This differs from the earlier 512-byte-leaf microbenchmarks; leaf size materially changes costs.

All corpora were re-tokenized offline using the cached GPT-2 tokenizer. Their stored canonical token prefixes match exactly. Full corpus hashes, source hashes, runtime versions, per-query samples and oracle expectations are in `paper_results/`. Each observed result was asserted against its oracle before timing data were accepted. The pre-tokenization manifest predates a final core range-wrapper refactor; each timed worker independently records its actual sources. Timed worker source consistency: **True**. The workstation is Windows 11, Intel Core i9-14900HX, Python 3.12.2; no CPU affinity/frequency locking was used. Team CPU-heavy tests were held during timed workers after the startup barrier; a sub-second preliminary biology probe may overlap the earliest work. These measurements are local CPU results, not a reproduction of the historical GPU server.

Construction is excluded from query latency, as in the original methodology. EXP015 query-rope construction is reported separately. Its three queries per prefix are averaged within each process before computing uncertainty. The full original ShareGPT/LMSYS real-pair and bounded K-candidate sweeps were rerun separately in nine fresh processes. They retain the original fixed timed arm order: radix, exact, NumPy.

## EXP005: exact LCP correctness survives; logarithmic exact-query claim does not

Times below are milliseconds, at a mismatch halfway through independent corpus-derived buffers.

| Bytes | Old hash-only | Exact | Python byte oracle |
|---:|---:|---:|---:|
| 64,000 | 17.7320 ± 0.4560 | 0.1152 ± 0.0053 | 0.8122 ± 0.0739 |
| 256,000 | 9.8738 ± 1.9517 | 0.0902 ± 0.0218 | 3.8978 ± 1.2418 |
| 1,000,000 | 10.2946 ± 2.7004 | 0.2587 ± 0.0788 | 13.3176 ± 1.7511 |
| 2,000,000 | 11.4460 ± 3.5314 | 0.4763 ± 0.0804 | 25.8967 ± 1.6714 |
| 4,000,000 | 14.1869 ± 1.9494 | 1.0075 ± 0.2342 | 62.6101 ± 9.3023 |
| 8,000,000 | 22.2575 ± 4.8108 | 2.2238 ± 0.5442 | 114.1726 ± 21.5212 |

Log-log slopes over 64K–8M bytes: old hash-only **0.042**, exact **0.646**, byte oracle **1.018**. The historical ≤0.3 exact-scaling gate fails; there is no general logarithmic exactness guarantee. Zero hash calls now measures fingerprint independence, not work. Fractions 0 and 0.99 are also present in the raw and summary data.

## EXP017: controlled competitive win remains, with different crossover

| Shared tokens | Old hash-only ms | Exact ms | Vendored radix ms | NumPy ms |
|---:|---:|---:|---:|---:|
| 1,000 | 21.5846 ± 4.3032 | 0.1310 ± 0.0155 | 0.0393 ± 0.0105 | 0.0028 ± 0.0013 |
| 4,000 | 20.1641 ± 1.5563 | 0.1319 ± 0.0362 | 0.1435 ± 0.0575 | 0.0033 ± 0.0010 |
| 16,000 | 18.4033 ± 3.8634 | 0.1261 ± 0.0514 | 0.5091 ± 0.0928 | 0.0066 ± 0.0035 |
| 64,000 | 16.8613 ± 3.5010 | 0.1848 ± 0.0800 | 3.5089 ± 1.2541 | 0.0547 ± 0.0792 |
| 128,000 | 19.7211 ± 3.5788 | 0.2158 ± 0.0334 | 5.6653 ± 1.2785 | 0.0570 ± 0.0302 |
| 256,000 | 19.9361 ± 3.7527 | 0.4627 ± 0.0978 | 10.5528 ± 3.0603 | 0.1380 ± 0.0705 |
| 512,000 | 19.0475 ± 3.3587 | 1.0725 ± 0.3790 | 18.9455 ± 1.9858 | 0.2125 ± 0.0543 |
| 1,000,000 | 20.3631 ± 2.9101 | 2.1540 ± 0.4602 | 35.6583 ± 3.2891 | 0.5094 ± 0.1982 |
| 2,000,000 | 15.5921 ± 2.4319 | 4.0392 ± 0.4225 | 75.2928 ± 8.7130 | 2.2386 ± 0.3887 |

At 2M tokens exact beats radix in **9/9** paired runs. The ratio of mean radix time to mean exact time is **18.64×**. Mean ± 1 SD bands remain separated from **16,000 tokens** upward on this grid. This is a measured CPU implementation result, not evidence for sublinear exact comparison. NumPy remains the flat-array reference. The original strict crossover gate **fails**: it also demands that radix's mean win at every smaller point; at 4K the means favor exact but the bands overlap. Thus this is a qualified long-context result, not a renewed full historical B1 promotion.

### Real conversational pairs and K candidates

The original 200 pairs per dataset per seed were rerun with exact, radix, and NumPy helpers. The table averages the per-pair median across each process's 200 pairs, then reports mean ± sample SD across nine fresh processes. Cached/query ropes are independently constructed. The same pairs recur across invocations of a seed; the 3,600 comparisons are not 3,600 independent samples.

| Dataset | Mean shared tokens | Exact ms/pair | Radix ms/pair | NumPy ms/pair |
|---|---:|---:|---:|---:|
| sharegpt_sample.jsonl | 956.4 | 0.0683 ± 0.0056 | 0.0336 ± 0.0028 | 0.0021 ± 0.0003 |
| lmsys_sample.jsonl | 554.6 | 0.0474 ± 0.0047 | 0.0206 ± 0.0030 | 0.0019 ± 0.0002 |

For K candidates, radix performs one tree lookup; Hashrope and NumPy perform K pairwise comparisons. The original workload shares token values but builds separate Hashrope roots; it does not add a Hashrope candidate index or share prefix nodes. Times below cover the whole lookup across K candidates, with construction excluded.

| Shared tokens / K | Exact ms | Radix ms | NumPy ms |
|---|---:|---:|---:|
| 64,000 / 1 | 0.2710 ± 0.0849 | 3.0999 ± 0.7014 | 0.0305 ± 0.0129 |
| 64,000 / 10 | 2.6117 ± 0.7488 | 2.9849 ± 0.5554 | 0.2180 ± 0.0583 |
| 64,000 / 100 | 20.5336 ± 4.8691 | 3.0623 ± 1.0679 | 4.0590 ± 0.5703 |
| 512,000 / 1 | 1.2100 ± 0.3584 | 18.6280 ± 3.8707 | 0.2408 ± 0.0491 |
| 512,000 / 10 | 10.4803 ± 2.0759 | 17.2971 ± 2.9295 | 3.9497 ± 1.2061 |

All real-pair and K-cell oracle gates passed. These cells test the corrected implementation directly; the historical hash-only penalties should not be carried forward as exact performance estimates.

## EXP015: CPU identification only

| Cached tokens | Old hash-only ms | Exact ms | Radix ms | NumPy ms | Query-rope build ms |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 18.5452 ± 1.5027 | 0.0015 ± 0.0003 | 0.0360 ± 0.0041 | 0.0020 ± 0.0002 | 1.9059 ± 0.1863 |
| 16,000 | 16.7420 ± 1.5239 | 0.0277 ± 0.0046 | 0.6260 ± 0.1502 | 0.0059 ± 0.0010 | 16.1006 ± 2.3400 |
| 64,000 | 15.2869 ± 1.9974 | 0.1360 ± 0.0539 | 3.6032 ± 0.5923 | 0.0447 ± 0.0287 | 61.0543 ± 5.0629 |
| 256,000 | 22.8357 ± 1.7649 | 0.5513 ± 0.1456 | 10.2416 ± 1.3453 | 0.1168 ± 0.0433 | 252.2974 ± 10.6192 |
| 1,000,000 | 20.4899 ± 2.2471 | 2.1341 ± 0.3438 | 37.7532 ± 4.2591 | 0.5918 ± 0.2091 | 968.4544 ± 43.9292 |

At 1M cached tokens, query-rope construction alone costs **968.4544 ± 43.9292 ms**, excluded from the **2.1341 ± 0.3438 ms** exact query. This eager construction cost is part of any end-to-end integration that builds each incoming query from flat tokens. Prebuilt exact identification can remain competitive while eager query-rope construction is much more expensive. The 4-A100 model-serving, TTFT and NVML energy measurements were not rerun; their historical engine-cache observations must not be presented as measurements of the corrected exact identification path. No new GPU/energy claim is supported here.

## New edit-plus-target comparison workload

This is separate from EXP020's maintained-fingerprint read. Each batch performs 50 seeded ASCII churn edits and compares each edited value with a prebuilt target: 50% equal, 50% different at the last byte. Independent targets and targets sharing edit ancestry are separate regimes. Target construction/destruction and warmup are untimed. Seeds 1/2/3 × three fresh processes, three batches per process; per-process medians are aggregated below. Every 50-bit equality-result sequence is checked exactly. Arm order rotates across seed/invocation. Ropey 1.6.1 uses its own exact equality; rehash uses Hashrope's polynomial hash over materialized Ropey bytes.

| Bytes / target history | Hashrope exact ms/cycle | Ropey exact | Hashrope candidate | Ropey rehashed candidate |
|---|---:|---:|---:|---:|
| 10,000 / independent | 0.0250 ± 0.0009 | 0.0012 ± 0.0001 | 0.0050 ± 0.0003 | 0.0461 ± 0.0015 |
| 10,000 / shared ancestry | 0.0190 ± 0.0008 | 0.0010 ± 0.0002 | 0.0047 ± 0.0002 | 0.0459 ± 0.0018 |
| 100,000 / independent | 0.1046 ± 0.0053 | 0.0116 ± 0.0008 | 0.0087 ± 0.0005 | 0.4357 ± 0.0122 |
| 100,000 / shared ancestry | 0.0356 ± 0.0021 | 0.0052 ± 0.0004 | 0.0089 ± 0.0004 | 0.4286 ± 0.0038 |
| 1,000,000 / independent | 0.9415 ± 0.1087 | 0.1124 ± 0.0124 | 0.0117 ± 0.0027 | 4.6212 ± 0.2905 |
| 1,000,000 / shared ancestry | 0.0610 ± 0.0082 | 0.0462 ± 0.0085 | 0.0090 ± 0.0011 | 5.1037 ± 0.1713 |

At 1M bytes, the ratio of mean Hashrope exact time to mean candidate time is **80.75×** for independent targets and **6.74×** for targets with shared ancestry. These ratios include the same edits in both arms; they measure this explicit workload, not a general equality bound. Ropey exact equality is faster in both 1M-byte regimes: Hashrope exact takes **8.37×** its time for independent targets and **1.32×** with shared ancestry.

A maintained fingerprint read is still cheap, but it is not an exact identity decision. The historical EXP020 speedup applies to its original maintained-versus-recomputed fingerprint workload; exact equality needs the separate target-dependent measurements above.

## Biology consumers

`bio_results.json` and `bio_paired.py` compare old/exact region, gene-report and resistance-panel decisions on 2,000 stored HIV patient records plus 16 independent equal controls against byte oracles. These are nine paired repetitions in one process, not nine independent process estimates. Construction is excluded; fixed 64-byte leaves and the data selection are recorded. The results do not reproduce the old chr22 gene-scale speedup claims or establish clinical efficacy.

| Whole corpus batch | Old hash-only median ms | Exact median ms | Exact / old |
|---|---:|---:|---:|
| regions | 0.8074 | 8.5084 | 10.54× |
| gene reports | 1.8406 | 10.2893 | 5.59× |
| panels | 135.5568 | 207.7366 | 1.53× |

Constructing the 2,016 cases took 946.75 ms and is excluded from these batch timings. The longer exact decisions are a real correction cost on these independently built inputs; tiny-codon panel queries amortize that cost differently.

## Unchanged or unavailable claims

EXP019 branch sharing/COW, raw edits, and compression do not depend on treating matching fingerprints as content equality. Their structural claims are not invalidated by this collision; their full historical benchmarks were not rerun here. Changes to default hash parameters, all large-scale/real-serving effects, and an end-to-end production reuse integration remain outside these measurements.

## Reproduce

The Python environment needs NumPy, Torch, Transformers and a locally cached GPT-2 tokenizer; tokenization deliberately runs offline. Rust requires Cargo and rustc. From the Hashrope checkout, capture the pinned paper baseline into a new temporary folder:

```powershell
$paperRepo = 'E:/data/code/claudecode/hashropeictaipaper'
$paperBaseline = Join-Path $env:TEMP ('hashrope-paper-baseline-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $paperBaseline | Out-Null
$paperArchive = Join-Path $paperBaseline 'baseline.zip'
git -C $paperRepo archive --format=zip --output=$paperArchive c9d8902e0f443f1a5ba12da78f087d5dbd2cad68
Expand-Archive -LiteralPath $paperArchive -DestinationPath $paperBaseline
python benchmarks/consumer_exactness/paper_cpu.py --paper $paperRepo --baseline $paperBaseline
python benchmarks/consumer_exactness/competitive_extended.py --paper $paperRepo
python benchmarks/consumer_exactness/edit_compare.py
python benchmarks/consumer_exactness/bio_paired.py --samples 2000 --output benchmarks/consumer_exactness/bio_results.json
python benchmarks/consumer_exactness/summarize.py
```

The Rust runner uses the checked-in dependency lock, fetches pinned Ropey 1.6.1 when necessary, and records Cargo.lock, binary hash and compile-time source hashes. A reused --binary must match its build manifest and current sources. Token arrays are ignored, regenerable preparation artifacts; their SHA-256 values remain in provenance. No historical latest JSON or figure is overwritten by these commands.
