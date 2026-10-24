# Independent diagnostic and recurrence review

The controlled results and report pass independent review. The sole report
correction requested was rounding concat's median block duration, 22.92225 ms,
to 22.9 rather than 23.0 ms. No material numerical, provenance or interpretation
issue was found. No builds or timings were run by this audit.

[controlled_review.py](controlled_review.py) and
[controlled_review.json](controlled_review.json) independently verify all
36 processes, 720 blocks, 153,552 raw samples and 90 estimate vectors, including
the source/archive/build/input/executable identities, saved context, commands,
sidecars, raw hashes, process medians, ratios and exploratory intervals. All
360 pinned endpoint pairs match the required processor. Additional read-only
checks reproduced the seeded case permutation, ABBA/BAAB parity, common session
counts, exact operation/sample counts, answer masks and node-growth constraints.
All 30 rendered primary estimates agree with the data. The 121.7-second wall
sum, 11.8-ms minimum block and 199 unpinned endpoint changes also agree.

[original_replication_review.py](original_replication_review.py) and its
[JSON result](original_replication_review.json) separately verify all 48 retained
executable processes, 8,784 samples, 66 paired vectors and CSV rows. The exact
four old executable hashes and frozen reference identities match. Its 100K
ratio is 1.0200 [0.9661, 1.0769], with 8/12 slower; 1 MiB LCP is 0.9933
[0.9387, 1.0511], with 7/12 slower; concat is 0.8476 [0.7657, 0.9383], with
11/12 faster. These observations remain separate from the prior nine and six
pairs and from the controlled data.

Here controlled **current** is frozen round six and **prior** is round five;
neither arm includes the proposed round-seven direction flag. Pinned LCP retains
a roughly 2% increase in both slot assignments: 1.0203 [1.0032, 1.0376] and
1.0200 [0.9957, 1.0448]. The normal-slot accumulated-time and wall intervals
also lie above one. The A/A estimate, 1.0114 [0.9991, 1.0238], is distinct
evidence of contextual variation. It cannot be subtracted as a measured causal
correction or used to label the approximately 2% increase noise.

Combined 100K estimates are mixed, but the pinned swapped-slot comparison-only
accumulated mean is 1.0516 [1.0182, 1.0862], with all six ratios slower. Its
primary median is 1.0198 with an interval crossing one. Both are correctly kept;
the different setup histories do not permit subtracting the phase timings.
The report also retains the consistent existing concat improvement without
claiming it explains the LCP/churn differences. No equivalence, universal
no-regression or identified thermal/allocator/placement cause is established.
