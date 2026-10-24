# Complete final trial results

**The candidate was rejected; production retains the baseline.** See the
[study report](../README.md) and [independent audit](../final_review.md).

All 48 processes, 8,784 samples and 66 cells are retained. Each ratio is the
geometric mean of twelve paired process-median ratios; intervals are
exploratory 95% log-t intervals, unadjusted across cells. Time columns are
medians of process medians and need not yield the same ratio. Slower counts
refer to the rejected candidate. Order columns identify which implementation
ran first. Full individual ratios, seed strata and sample vectors remain in
[summary.json](summary.json) and [results.json](results.json).

Baseline is frozen round six; candidate is the direction-argument trial.
The source-match flags in the summary precede the documented restoration.
No cells or unfavorable samples are removed. Hash-only rows are historical
candidate controls, not substitutes for exact byte comparison.

Summary SHA-256: `7cecad18844d0fd483e74bd640fc7d732eb55018b3bef3669dab025a29a90d35`.

## Microbenchmark

| Cell | Baseline | Rejected candidate | Ratio [interval] | Slower | Ratio range | Baseline-first / candidate-first |
|---|---:|---:|---:|---:|---:|---:|
| `compare/1048576/early_mismatch/exact_equal` | 58.2250 ns | 59.5500 ns | 1.0394 [1.0054, 1.0745] | 9/12 | 0.9614–1.1570 | 1.0400 / 1.0388 |
| `compare/1048576/early_mismatch/exact_lcp` | 57.2000 ns | 57.4000 ns | 1.0345 [0.9986, 1.0717] | 7/12 | 0.9520–1.1662 | 1.0541 / 1.0152 |
| `compare/1048576/early_mismatch/hash_only_lcp` | 955.0000 ns | 956.0000 ns | 1.0032 [0.9603, 1.0480] | 6/12 | 0.8694–1.1168 | 0.9842 / 1.0225 |
| `compare/1048576/independent_equal/exact_equal` | 53361.5000 ns | 54192.0000 ns | 1.0217 [0.9673, 1.0791] | 8/12 | 0.8801–1.2720 | 1.0720 / 0.9737 |
| `compare/1048576/independent_equal/exact_lcp` | 54724.5000 ns | 53870.5000 ns | 0.9964 [0.9602, 1.0340] | 6/12 | 0.8653–1.0753 | 0.9846 / 1.0084 |
| `compare/1048576/independent_equal/hash_only_lcp` | 6289.5000 ns | 6661.0000 ns | 1.0387 [1.0106, 1.0676] | 11/12 | 0.9583–1.0925 | 1.0328 / 1.0447 |
| `compare/1048576/late_mismatch/exact_equal` | 53157.5000 ns | 53896.0000 ns | 1.0196 [0.9860, 1.0543] | 7/12 | 0.9397–1.1534 | 1.0094 / 1.0298 |
| `compare/1048576/late_mismatch/exact_lcp` | 54679.0000 ns | 53595.5000 ns | 0.9707 [0.9446, 0.9976] | 1/12 | 0.8722–1.0377 | 0.9695 / 0.9720 |
| `compare/1048576/late_mismatch/hash_only_lcp` | 6650.0000 ns | 6575.5000 ns | 0.9910 [0.9665, 1.0162] | 6/12 | 0.9435–1.0614 | 1.0101 / 0.9723 |
| `compare/1048576/same_root/exact_equal` | 3.6500 ns | 3.7000 ns | 0.9626 [0.8599, 1.0776] | 5/12 | 0.5620–1.1000 | 0.9121 / 1.0159 |
| `compare/1048576/same_root/exact_lcp` | 1.7500 ns | 1.8500 ns | 0.9693 [0.8725, 1.0767] | 4/12 | 0.5846–1.1176 | 0.9446 / 0.9945 |
| `compare/1048576/same_root/hash_only_lcp` | 6353.0000 ns | 6695.0000 ns | 1.0324 [1.0024, 1.0632] | 8/12 | 0.9598–1.0912 | 1.0105 / 1.0547 |
| `compare/1048576/shared_late_mismatch/exact_equal` | 99.0000 ns | 97.5000 ns | 0.9432 [0.8530, 1.0430] | 6/12 | 0.5904–1.1031 | 1.0105 / 0.8804 |
| `compare/1048576/shared_late_mismatch/exact_lcp` | 104.0000 ns | 97.5000 ns | 0.9498 [0.8761, 1.0298] | 3/12 | 0.7007–1.1556 | 0.9942 / 0.9075 |
| `compare/1048576/shared_late_mismatch/hash_only_lcp` | 6595.5000 ns | 6473.5000 ns | 0.9803 [0.9443, 1.0176] | 6/12 | 0.9108–1.0710 | 1.0157 / 0.9461 |
| `compare/65536/early_mismatch/exact_equal` | 39.8250 ns | 42.6500 ns | 1.0364 [0.9975, 1.0769] | 7/12 | 0.9239–1.1260 | 1.0169 / 1.0563 |
| `compare/65536/early_mismatch/exact_lcp` | 38.0250 ns | 40.7500 ns | 1.0491 [1.0121, 1.0874] | 10/12 | 0.9186–1.1075 | 1.0164 / 1.0829 |
| `compare/65536/early_mismatch/hash_only_lcp` | 691.5000 ns | 700.0000 ns | 1.0099 [0.9708, 1.0505] | 8/12 | 0.9193–1.1106 | 1.0051 / 1.0146 |
| `compare/65536/independent_equal/exact_equal` | 3113.5000 ns | 3005.5000 ns | 0.9504 [0.8666, 1.0424] | 7/12 | 0.7440–1.1005 | 0.8948 / 1.0095 |
| `compare/65536/independent_equal/exact_lcp` | 3027.0000 ns | 2987.0000 ns | 0.9491 [0.8858, 1.0170] | 6/12 | 0.7905–1.1061 | 0.9064 / 0.9939 |
| `compare/65536/independent_equal/hash_only_lcp` | 4581.5000 ns | 4702.0000 ns | 1.0179 [0.9891, 1.0475] | 8/12 | 0.9393–1.0957 | 0.9991 / 1.0370 |
| `compare/65536/late_mismatch/exact_equal` | 2883.0000 ns | 2989.5000 ns | 0.9901 [0.9403, 1.0426] | 7/12 | 0.8360–1.0954 | 0.9612 / 1.0198 |
| `compare/65536/late_mismatch/exact_lcp` | 3060.5000 ns | 3042.5000 ns | 0.9867 [0.9295, 1.0473] | 6/12 | 0.8277–1.1685 | 0.9506 / 1.0241 |
| `compare/65536/late_mismatch/hash_only_lcp` | 4430.5000 ns | 4720.0000 ns | 1.0347 [1.0069, 1.0633] | 10/12 | 0.9566–1.1046 | 1.0231 / 1.0464 |
| `compare/65536/same_root/exact_equal` | 3.7000 ns | 3.7000 ns | 0.9508 [0.8118, 1.1136] | 4/12 | 0.4431–1.1449 | 1.0069 / 0.8979 |
| `compare/65536/same_root/exact_lcp` | 1.8250 ns | 1.8500 ns | 1.0136 [0.9769, 1.0517] | 4/12 | 0.9189–1.1176 | 1.0147 / 1.0126 |
| `compare/65536/same_root/hash_only_lcp` | 4677.0000 ns | 4789.5000 ns | 1.0329 [0.9956, 1.0717] | 8/12 | 0.9364–1.1053 | 1.0073 / 1.0592 |
| `compare/65536/shared_late_mismatch/exact_equal` | 71.0000 ns | 77.0000 ns | 1.0261 [0.9918, 1.0615] | 7/12 | 0.9091–1.1143 | 1.0134 / 1.0389 |
| `compare/65536/shared_late_mismatch/exact_lcp` | 70.0000 ns | 73.5000 ns | 1.0166 [0.9792, 1.0554] | 7/12 | 0.9000–1.0946 | 1.0206 / 1.0126 |
| `compare/65536/shared_late_mismatch/hash_only_lcp` | 4428.0000 ns | 4566.5000 ns | 1.0248 [0.9888, 1.0621] | 9/12 | 0.9258–1.0964 | 1.0048 / 1.0453 |
| `compare/repeat_2TiB/exact_equal` | 82.8750 ns | 82.4750 ns | 0.9933 [0.9197, 1.0728] | 6/12 | 0.8273–1.3163 | 1.0421 / 0.9468 |
| `compare/repeat_2TiB/exact_lcp` | 81.5750 ns | 80.6250 ns | 0.9592 [0.8964, 1.0263] | 5/12 | 0.7433–1.1128 | 0.9785 / 0.9402 |
| `compare/repeat_2TiB/hash_only_lcp` | 27982.5000 ns | 27786.0000 ns | 0.9928 [0.9629, 1.0236] | 6/12 | 0.9267–1.0705 | 1.0067 / 0.9790 |
| `compare/self_concat_DAG_2TiB/exact_equal` | 145.0000 ns | 142.0000 ns | 0.9797 [0.8897, 1.0788] | 3/12 | 0.7278–1.3621 | 0.9982 / 0.9615 |
| `compare/self_concat_DAG_2TiB/exact_lcp` | 124.0000 ns | 124.0000 ns | 0.9962 [0.9046, 1.0970] | 5/12 | 0.7914–1.3964 | 0.9758 / 1.0170 |
| `compare/self_concat_DAG_2TiB/hash_only_lcp` | 24566.0000 ns | 24472.5000 ns | 0.9854 [0.9548, 1.0171] | 7/12 | 0.8968–1.0321 | 0.9862 / 0.9846 |
| `core/concat_balanced_64KiB` | 13.8750 ns | 14.2500 ns | 0.9711 [0.8037, 1.1734] | 7/12 | 0.4927–1.3460 | 0.8137 / 1.1590 |
| `core/from_bytes_64KiB` | 30312.5000 ns | 30082.5000 ns | 0.9725 [0.8142, 1.1616] | 8/12 | 0.6237–1.4407 | 0.9184 / 1.0298 |
| `core/repeat_join_2pow32` | 5386.0000 ns | 5423.5000 ns | 1.0042 [0.9707, 1.0389] | 6/12 | 0.9212–1.1032 | 0.9900 / 1.0186 |
| `core/root_hash` | 0.5530 ns | 0.5637 ns | 1.0050 [0.9660, 1.0455] | 9/12 | 0.8893–1.0772 | 1.0055 / 1.0044 |
| `core/split_chunked_64KiB` | 444.5000 ns | 454.0000 ns | 1.0492 [0.9813, 1.1218] | 9/12 | 0.8364–1.2793 | 0.9861 / 1.1162 |
| `core/substr_warm_64KiB` | 476.2000 ns | 474.5000 ns | 1.0077 [0.9783, 1.0381] | 7/12 | 0.9335–1.0866 | 0.9774 / 1.0390 |

## Edit/comparison workload

| Cell | Baseline | Rejected candidate | Ratio [interval] | Slower | Ratio range | Baseline-first / candidate-first |
|---|---:|---:|---:|---:|---:|---:|
| `10000/false/hashrope_candidate` | 2.2600 µs | 2.2490 µs | 0.9952 [0.9404, 1.0532] | 5/12 | 0.8775–1.1499 | 1.0189 / 0.9721 |
| `10000/false/hashrope_exact` | 4.7250 µs | 4.8090 µs | 1.0178 [0.9652, 1.0734] | 7/12 | 0.8890–1.1894 | 1.0186 / 1.0171 |
| `10000/false/ropey_exact` | 1.2330 µs | 1.1660 µs | 0.9425 [0.8369, 1.0613] | 4/12 | 0.7139–1.5033 | 0.9825 / 0.9041 |
| `10000/false/ropey_rehash_candidate` | 4.8360 µs | 5.1380 µs | 1.0463 [0.9855, 1.1108] | 7/12 | 0.9100–1.1840 | 1.0429 / 1.0497 |
| `10000/true/hashrope_candidate` | 1.8970 µs | 2.0400 µs | 1.0567 [1.0112, 1.1043] | 11/12 | 0.9389–1.1912 | 1.0598 / 1.0535 |
| `10000/true/hashrope_exact` | 3.8620 µs | 4.0130 µs | 1.0209 [0.9531, 1.0936] | 7/12 | 0.8312–1.1808 | 1.0346 / 1.0074 |
| `10000/true/ropey_exact` | 0.9810 µs | 0.9680 µs | 1.0006 [0.9515, 1.0522] | 7/12 | 0.8933–1.1677 | 0.9897 / 1.0116 |
| `10000/true/ropey_rehash_candidate` | 4.7340 µs | 4.8100 µs | 1.0218 [0.9706, 1.0758] | 9/12 | 0.8437–1.1780 | 1.0129 / 1.0308 |
| `100000/false/hashrope_candidate` | 3.8790 µs | 3.8260 µs | 1.0024 [0.9589, 1.0480] | 6/12 | 0.8992–1.1362 | 1.0160 / 0.9891 |
| `100000/false/hashrope_exact` | 20.1970 µs | 20.1630 µs | 0.9898 [0.9479, 1.0336] | 5/12 | 0.8628–1.0896 | 1.0082 / 0.9717 |
| `100000/false/ropey_exact` | 11.5060 µs | 11.5770 µs | 1.0160 [0.9464, 1.0908] | 5/12 | 0.8959–1.3147 | 1.0454 / 0.9874 |
| `100000/false/ropey_rehash_candidate` | 39.9170 µs | 40.4890 µs | 1.0208 [0.9933, 1.0490] | 6/12 | 0.9594–1.1172 | 1.0285 / 1.0131 |
| `100000/true/hashrope_candidate` | 3.9130 µs | 3.9380 µs | 1.0214 [0.9690, 1.0766] | 6/12 | 0.9459–1.2777 | 1.0537 / 0.9901 |
| `100000/true/hashrope_exact` | 7.5490 µs | 7.9470 µs | 1.0419 [0.9873, 1.0995] | 7/12 | 0.9252–1.2005 | 1.0612 / 1.0230 |
| `100000/true/ropey_exact` | 5.3670 µs | 5.5050 µs | 1.0167 [0.9270, 1.1151] | 7/12 | 0.8119–1.3193 | 1.0006 / 1.0330 |
| `100000/true/ropey_rehash_candidate` | 39.6890 µs | 40.5560 µs | 1.0217 [0.9866, 1.0580] | 7/12 | 0.9533–1.1266 | 1.0119 / 1.0315 |
| `1000000/false/hashrope_candidate` | 6.1460 µs | 5.6760 µs | 0.9275 [0.8638, 0.9959] | 2/12 | 0.7405–1.1115 | 0.9381 / 0.9171 |
| `1000000/false/hashrope_exact` | 153.6030 µs | 152.1710 µs | 1.0014 [0.9556, 1.0495] | 6/12 | 0.8677–1.1155 | 0.9830 / 1.0202 |
| `1000000/false/ropey_exact` | 118.1530 µs | 115.7770 µs | 0.9886 [0.9262, 1.0551] | 5/12 | 0.8763–1.2318 | 1.0095 / 0.9680 |
| `1000000/false/ropey_rehash_candidate` | 415.0570 µs | 414.5770 µs | 0.9605 [0.8816, 1.0465] | 7/12 | 0.6337–1.0452 | 0.9121 / 1.0115 |
| `1000000/true/hashrope_candidate` | 3.6080 µs | 3.6500 µs | 1.0212 [0.9619, 1.0842] | 7/12 | 0.8874–1.2713 | 1.0470 / 0.9961 |
| `1000000/true/hashrope_exact` | 9.1250 µs | 8.6350 µs | 0.9150 [0.8091, 1.0347] | 3/12 | 0.5183–1.1001 | 0.8510 / 0.9837 |
| `1000000/true/ropey_exact` | 37.5260 µs | 35.8890 µs | 0.9666 [0.8617, 1.0843] | 5/12 | 0.7105–1.3519 | 0.9727 / 0.9605 |
| `1000000/true/ropey_rehash_candidate` | 421.9440 µs | 415.9740 µs | 0.9847 [0.9613, 1.0087] | 4/12 | 0.9189–1.0348 | 0.9762 / 0.9934 |

## Layout

Both arms report a 48-byte `NodeInner`, 2,136-byte inline `Arena`,
and 255 nodes for the 64 KiB construction. No layout change occurred.
