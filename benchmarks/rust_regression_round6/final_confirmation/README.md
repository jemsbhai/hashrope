# Fixed final confirmation: six new pairs, separate results

All 24 predeclared retained-binary processes completed successfully, preserving
4,392 raw timing samples and all 66 cells. No executable was rebuilt, no source
changed, and no observation was replaced or removed. Baseline is the frozen
round-five optimization; candidate is the reviewed round-six optimization.
This is a separate six-pair follow-up to the initial nine-pair final comparison.
It was declared after observing that cohort's 100K result, and is not pooled
with or substituted for it.

The concat point estimate improves in both phases. The initial 100K independent
exact-churn slowdown does not repeat as a paired geometric mean, but the
confirmation still has four of six slower pairs and a higher candidate median
across processes. The results do not support a blanket no-regression claim.

## Primary outcomes

Each entry is candidate/baseline paired geometric mean, exploratory 95% log-t
interval, and slower-direction count. These process-level intervals are not
multiplicity adjusted or equivalence bounds.

| Outcome | Initial nine pairs | Separate confirmation six pairs |
|---|---:|---:|
| Balanced concat | 0.8550 [0.7674, 0.9526]; 1/9 slower | 0.8777 [0.7252, 1.0624]; 1/6 slower |
| Independent 100K exact churn | 1.0438 [1.0015, 1.0879]; 8/9 slower | 0.9834 [0.9228, 1.0479]; 4/6 slower |
| Independent 1M exact churn | 1.0138 [0.9834, 1.0452]; 6/9 slower | 0.9604 [0.9070, 1.0169]; 2/6 slower |

Confirmation medians of process medians are 17.100→14.725 ns for concat,
18.113→18.406 µs for 100K exact churn, and 146.035→140.095 µs for 1M exact churn.
The 100K median ratio and paired geometric mean have different directions:
they are different estimators, and the paired ratios range from 0.8983 to 1.0465.
Neither estimator is discarded. Baseline-first/candidate-first 100K ratio
geomeans are 0.9548/1.0128; concat ratios are 0.7843/0.9822. The order breakdown
is retained, not used to select a preferred subset.

Secondary comparison results also reverse between phases. Independent 1MiB
exact LCP changes from an initial ratio of 0.9801 to confirmation 1.0680
[0.9921, 1.1496], with five of six slower pairs. Independent 64KiB exact LCP
changes from 0.9810 to 1.0776 [0.9364, 1.2400], with three of six slower pairs;
late-mismatch 1MiB LCP changes from 0.9825 to 1.0410 [0.9597, 1.1292]. These
observations remain unresolved performance limitations. No confirmation
slowdown interval lies wholly above one, but that does not establish zero cost,
prove noise, or make the initial observations disappear. Unchanged-control
variation and source/disassembly similarity do not identify a timing cause.

## Every retained confirmation cell

Times below are medians of the six process medians in each arm. Micro times
are **ns**; churn times are **ms per edit/comparison**. Ratios use the six matched
process medians, so they need not equal the ratio of the two displayed times.
All per-process values, ratio ranges, order/seed breakdowns and layout rows are
in [summary.json](summary.json); [all_cells.csv](all_cells.csv) contains all 66
rows. The first nine-pair results remain in
[the separate initial summary](../final_comparison/summary.json).

### Microbenchmark (ns)

| Cell | Baseline | Candidate | Paired ratio [95% interval]; slower |
|---|---:|---:|---:|
| `compare/1048576/early_mismatch/exact_equal` | 59.35 | 59.475 | 0.9964 [0.9776, 1.0155]; 2/6 slower |
| `compare/1048576/early_mismatch/exact_lcp` | 57.225 | 57.225 | 0.9927 [0.9613, 1.0251]; 3/6 slower |
| `compare/1048576/early_mismatch/hash_only_lcp` | 964 | 955.5 | 1.0303 [0.9478, 1.1200]; 3/6 slower |
| `compare/1048576/independent_equal/exact_equal` | 52698 | 53346 | 1.0307 [0.9567, 1.1105]; 4/6 slower |
| `compare/1048576/independent_equal/exact_lcp` | 53279 | 55668.5 | 1.0680 [0.9921, 1.1496]; 5/6 slower |
| `compare/1048576/independent_equal/hash_only_lcp` | 6631.5 | 6363.5 | 0.9818 [0.9315, 1.0347]; 2/6 slower |
| `compare/1048576/late_mismatch/exact_equal` | 51778.5 | 53982.5 | 1.0408 [0.9115, 1.1883]; 5/6 slower |
| `compare/1048576/late_mismatch/exact_lcp` | 52757 | 54775 | 1.0410 [0.9597, 1.1292]; 4/6 slower |
| `compare/1048576/late_mismatch/hash_only_lcp` | 6681.5 | 6637.5 | 0.9964 [0.9557, 1.0389]; 3/6 slower |
| `compare/1048576/same_root/exact_equal` | 3.8 | 3.65 | 0.9752 [0.9396, 1.0122]; 1/6 slower |
| `compare/1048576/same_root/exact_lcp` | 1.85 | 1.8 | 0.9779 [0.8851, 1.0804]; 2/6 slower |
| `compare/1048576/same_root/hash_only_lcp` | 6645.5 | 6607 | 1.0005 [0.9441, 1.0603]; 3/6 slower |
| `compare/1048576/shared_late_mismatch/exact_equal` | 105 | 105.5 | 0.9856 [0.9238, 1.0515]; 2/6 slower |
| `compare/1048576/shared_late_mismatch/exact_lcp` | 100.5 | 102.5 | 1.0099 [0.9337, 1.0924]; 3/6 slower |
| `compare/1048576/shared_late_mismatch/hash_only_lcp` | 6506.5 | 6443 | 1.0092 [0.9437, 1.0793]; 2/6 slower |
| `compare/65536/early_mismatch/exact_equal` | 42.7 | 42.775 | 0.9924 [0.9693, 1.0161]; 3/6 slower |
| `compare/65536/early_mismatch/exact_lcp` | 40.8 | 40.625 | 0.9840 [0.9618, 1.0068]; 1/6 slower |
| `compare/65536/early_mismatch/hash_only_lcp` | 690 | 693 | 0.9947 [0.9718, 1.0181]; 5/6 slower |
| `compare/65536/independent_equal/exact_equal` | 2950.5 | 2983.5 | 1.0253 [0.8563, 1.2275]; 3/6 slower |
| `compare/65536/independent_equal/exact_lcp` | 2946 | 3030.5 | 1.0776 [0.9364, 1.2400]; 3/6 slower |
| `compare/65536/independent_equal/hash_only_lcp` | 4758.5 | 4758 | 1.0015 [0.9375, 1.0699]; 3/6 slower |
| `compare/65536/late_mismatch/exact_equal` | 3011.5 | 3008 | 1.0735 [0.9583, 1.2027]; 3/6 slower |
| `compare/65536/late_mismatch/exact_lcp` | 3002 | 3017 | 1.0747 [0.9619, 1.2006]; 3/6 slower |
| `compare/65536/late_mismatch/hash_only_lcp` | 4913 | 4713.5 | 0.9709 [0.9033, 1.0435]; 1/6 slower |
| `compare/65536/same_root/exact_equal` | 3.725 | 3.85 | 1.0135 [0.9531, 1.0778]; 2/6 slower |
| `compare/65536/same_root/exact_lcp` | 1.9 | 1.85 | 0.9604 [0.9384, 0.9829]; 0/6 slower |
| `compare/65536/same_root/hash_only_lcp` | 4809 | 4737.5 | 0.9853 [0.9373, 1.0357]; 1/6 slower |
| `compare/65536/shared_late_mismatch/exact_equal` | 78 | 77 | 0.9891 [0.9281, 1.0541]; 1/6 slower |
| `compare/65536/shared_late_mismatch/exact_lcp` | 75 | 75 | 1.0072 [0.9592, 1.0576]; 2/6 slower |
| `compare/65536/shared_late_mismatch/hash_only_lcp` | 4575.5 | 4614.5 | 1.0079 [0.9562, 1.0624]; 2/6 slower |
| `compare/repeat_2TiB/exact_equal` | 82.825 | 78.75 | 0.9488 [0.8803, 1.0225]; 1/6 slower |
| `compare/repeat_2TiB/exact_lcp` | 80.075 | 80.075 | 1.0383 [0.8824, 1.2219]; 2/6 slower |
| `compare/repeat_2TiB/hash_only_lcp` | 28627.5 | 27828 | 0.9645 [0.8797, 1.0573]; 2/6 slower |
| `compare/self_concat_DAG_2TiB/exact_equal` | 158 | 158 | 0.9933 [0.8723, 1.1310]; 2/6 slower |
| `compare/self_concat_DAG_2TiB/exact_lcp` | 125 | 124.5 | 1.0131 [0.8814, 1.1645]; 2/6 slower |
| `compare/self_concat_DAG_2TiB/hash_only_lcp` | 25405 | 24971 | 0.9708 [0.9149, 1.0302]; 2/6 slower |
| `core/concat_balanced_64KiB` | 17.1 | 14.725 | 0.8777 [0.7252, 1.0624]; 1/6 slower |
| `core/from_bytes_64KiB` | 42527.5 | 41930 | 0.8958 [0.6911, 1.1612]; 3/6 slower |
| `core/repeat_join_2pow32` | 5321.5 | 5388 | 1.0160 [0.9915, 1.0412]; 5/6 slower |
| `core/root_hash` | 0.55935 | 0.55895 | 1.0029 [0.9769, 1.0295]; 3/6 slower |
| `core/split_chunked_64KiB` | 451.5 | 452 | 0.9930 [0.9261, 1.0648]; 2/6 slower |
| `core/substr_warm_64KiB` | 488.95 | 462.8 | 0.9729 [0.9366, 1.0106]; 1/6 slower |

### Edit/comparison workload (ms per operation)

| Cell | Baseline | Candidate | Paired ratio [95% interval]; slower |
|---|---:|---:|---:|
| `10000/false/hashrope_candidate` | 0.002088 | 0.002202 | 1.0375 [0.9493, 1.1339]; 4/6 slower |
| `10000/false/hashrope_exact` | 0.004669 | 0.004688 | 1.0151 [0.9504, 1.0841]; 3/6 slower |
| `10000/false/ropey_exact` | 0.001115 | 0.001134 | 1.0333 [0.8298, 1.2866]; 3/6 slower |
| `10000/false/ropey_rehash_candidate` | 0.004818 | 0.004956 | 1.0314 [0.9365, 1.1358]; 3/6 slower |
| `10000/true/hashrope_candidate` | 0.001892 | 0.001905 | 1.0118 [0.8888, 1.1517]; 3/6 slower |
| `10000/true/hashrope_exact` | 0.00381 | 0.003781 | 1.0319 [0.9595, 1.1099]; 5/6 slower |
| `10000/true/ropey_exact` | 0.000912 | 0.000952 | 1.0541 [0.9487, 1.1712]; 4/6 slower |
| `10000/true/ropey_rehash_candidate` | 0.004562 | 0.004699 | 1.1304 [0.8657, 1.4760]; 4/6 slower |
| `100000/false/hashrope_candidate` | 0.0037 | 0.003772 | 1.0347 [0.9579, 1.1177]; 4/6 slower |
| `100000/false/hashrope_exact` | 0.018113 | 0.018406 | 0.9834 [0.9228, 1.0479]; 4/6 slower |
| `100000/false/ropey_exact` | 0.010654 | 0.011108 | 1.0616 [0.9273, 1.2155]; 3/6 slower |
| `100000/false/ropey_rehash_candidate` | 0.038665 | 0.039235 | 1.0255 [0.9803, 1.0727]; 4/6 slower |
| `100000/true/hashrope_candidate` | 0.004162 | 0.004088 | 1.0026 [0.8840, 1.1372]; 2/6 slower |
| `100000/true/hashrope_exact` | 0.007619 | 0.007284 | 0.9403 [0.8694, 1.0170]; 1/6 slower |
| `100000/true/ropey_exact` | 0.005165 | 0.00492 | 0.9491 [0.8793, 1.0244]; 2/6 slower |
| `100000/true/ropey_rehash_candidate` | 0.039513 | 0.038916 | 0.9811 [0.9229, 1.0429]; 1/6 slower |
| `1000000/false/hashrope_candidate` | 0.005677 | 0.005425 | 0.9766 [0.8965, 1.0640]; 2/6 slower |
| `1000000/false/hashrope_exact` | 0.146035 | 0.140095 | 0.9604 [0.9070, 1.0169]; 2/6 slower |
| `1000000/false/ropey_exact` | 0.110423 | 0.108114 | 0.9705 [0.8678, 1.0853]; 2/6 slower |
| `1000000/false/ropey_rehash_candidate` | 0.400253 | 0.415836 | 1.0082 [0.9249, 1.0990]; 4/6 slower |
| `1000000/true/hashrope_candidate` | 0.003609 | 0.003341 | 0.9135 [0.7988, 1.0446]; 2/6 slower |
| `1000000/true/hashrope_exact` | 0.008999 | 0.008689 | 0.9601 [0.8906, 1.0350]; 1/6 slower |
| `1000000/true/ropey_exact` | 0.037181 | 0.035102 | 0.9607 [0.8473, 1.0893]; 2/6 slower |
| `1000000/true/ropey_rehash_candidate` | 0.419994 | 0.401037 | 0.9641 [0.9211, 1.0092]; 1/6 slower |

## Protocol and provenance

The fixed [PLAN.md](PLAN.md) names the three primary outcomes before timing,
uses six balanced alternating process pairs for each cohort, and repeats churn
seeds 1/2/3 twice with invocation zero. The exact four final-comparison archived
executables preserve the original 42×7 micro and 24×3 churn cell/sample counts,
complete sequences, warmups, iterations, byte oracles and layout rows. All
micro processes precede churn; affinity is unpinned. CPU frequency, temperature
and unrelated OS activity remain uncontrolled, while other agent builds, tests
and timings were paused. No outer Python timing interleaving was introduced.

[manifest.json](manifest.json) binds the source, workload, lock, executable,
initial-result, initial-summary and confirmation-driver identities. Every raw
stdout/stderr and process sidecar remains under [raw](raw). The analysis reparses
all samples and verifies their status, argument/seed/order, cell counts, expected
binary hashes and sidecars before producing [results.json](results.json) and the
separate confirmation summary. The retained driver has no build command and
refuses output overwrites. No version bump, commit or publication occurs here.

An [independent review](../final_review.md) verified all 24 raw processes, 4,392
samples, 66 summary vectors and rendered table rows, all four retained binaries
and source/parent bindings. Its [audit results](../confirmation_review.json) and
[recomputation script](../confirmation_review.py) are retained separately.
