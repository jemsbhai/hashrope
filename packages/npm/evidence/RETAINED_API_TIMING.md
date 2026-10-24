# npm retained-API timing diagnostic

median of three process medians; each process is median of six within-process round medians, each from five batches. Min/max process and all-round ranges retained. No confidence/significance claim.

six within-process rounds; fixture order rotated/reversed; each fixture has three AB and three BA arm orders; five batches per arm/round; equal repetitions; setup excluded; bounded warmup; no affinity/frequency control

Runtime: v24.22.1; CPU: Intel(R) Core(TM) i9-14900HX. Baseline is the verified published 0.2.2 tarball; candidate is frozen 0.3.3 runtime code. Units below are microseconds per operation.

| Fixture | 0.2.2 median [process range] | 0.3.3 median [process range] | Candidate/baseline time |
|---|---:|---:|---:|
| hash_4096B | 361.360 [343.975, 369.002] | 359.887 [339.238, 366.370] | 0.996x |
| construct_leaf_128B | 12.940 [12.682, 13.172] | 13.189 [13.093, 13.190] | 1.019x |
| concat_two_128B_warm | 0.171 [0.169, 0.178] | 0.219 [0.215, 0.227] | 1.285x |
| split_64KiB_leaf_middle | 5835.550 [5662.075, 6021.200] | 5908.975 [5679.850, 6007.075] | 1.013x |
| substring_64B_in_64KiB_leaf | 25.524 [25.415, 25.680] | 25.488 [25.299, 27.343] | 0.999x |
| repeat_metadata_100bit | 14.467 [14.343, 14.796] | 14.430 [14.325, 14.466] | 0.997x |
| equal_independent_64KiB | 36.346 [36.080, 37.248] | 16.875 [16.874, 17.357] | 0.464x |
| equal_late_mismatch_64KiB | 88.407 [87.120, 89.157] | 16.480 [16.410, 16.861] | 0.186x |
| equal_early_mismatch_64KiB | 41.651 [41.309, 42.236] | 0.150 [0.142, 0.153] | 0.004x |
| equal_shared_prefix_64KiB | 89.804 [88.290, 91.899] | 0.691 [0.658, 0.732] | 0.008x |
| equal_same_root | 0.004 [0.004, 0.004] | 0.007 [0.006, 0.007] | 1.632x |

These bounded, warmed fixtures are diagnostic evidence, not broad package performance parity. Setup is excluded and there is no affinity/frequency control. Tiny operations are particularly sensitive to call/measurement overhead. The independent-equality candidate initially regressed (~2x) with a scalar leaf helper; the initial raw run and archive are retained separately. A portable DataView word scan was selected only after that observation and byte-oracle review. Unchanged arithmetic timings also vary, so no causal interpretation is assigned to their small movements. Profile safety adds checks and stored profile references; zero overhead is not claimed.

Huge Repeat and independent DAG cases each completed with the expected exact result in all three processes. They are single completion observations per case, not comparable throughput samples; the old implementation was not asked to traverse astronomical expanded content. Exact comparisons still have general traversal, memoization, BigInt and unshared-byte costs.

Raw samples: `evidence/retained-api-timing-final-{0,1,2}.json`; frozen runtime and harness digests are embedded in every file, with archive provenance in `evidence/retained-api-timing-manifest-final.json`. `evidence/retained-api-timing-initial.json` preserves the rejected scalar-helper diagnostic. Final documentation-only packaging changes are outside the measured runtime source scope.
