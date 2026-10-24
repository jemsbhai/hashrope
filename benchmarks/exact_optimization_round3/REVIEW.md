# Independent review of the selected round-three changes

Reviewed 2026-10-23 against the captured round-three baseline. This review covers
the selected production diff and the edit diagnostic's workload, correctness
checks and allocation accounting. It does not substitute for the final paired
benchmark results. No production or benchmark source was changed during this
review, and no additional timed work was run.

**Conclusion:** no correctness blocker found in the selected production changes,
and no allocation-accounting defect found in the diagnostic. The reporting
limits below are material when interpreting its results.

## Production behavior and arithmetic

The production diff is limited to the default-profile byte-hash implementation
and borrowing leaf data in `ensure_hash_inner`, `split_inner`, and `hash_range`.
The generic reducer, power cache, accepted profiles, node representation,
balancing decisions and exact comparison engines are unchanged.

For the default profile, grouping eight byte updates preserves the same
polynomial modulo `p = 2^61 - 1`: the next value is
`h * 131^8 + block`, where `block` evaluates the same eight `byte + 1`
coefficients without intermediate reduction. The relevant exact bounds are:

- `131^8 = 86,730,203,469,006,241`.
- The largest block is `256 * (131^8 - 1) / 130 = 170,791,785,292,812,288`,
  below `2^58`, so its `u64` arithmetic cannot overflow.
- With `h <= p - 1`, the largest combined accumulator is
  `199,986,233,356,689,308,921,606,600,506,754,238`, below `2^118`, so the
  `u128` multiplication and addition cannot overflow.
- At that bound, the existing two Mersenne folds followed by one conditional
  subtraction produce the usual residue below `p`. Thus the bound on `h`
  continues to hold after every full block. Remaining bytes use the original
  scalar recurrence, and empty input still returns zero.

The specialization requires both the original default modulus and base 131.
That restriction is essential: accepted legacy non-Mersenne profiles do not
necessarily implement ordinary modular arithmetic. Regrouping their operations
could change fingerprints even if a conventional polynomial derivation appeared
valid. In particular, profile `(17, 2)` retains the existing hash 16 for `[0, 13]`
in one leaf and hash 0 for the same bytes after split/rejoin.

The borrowed leaf changes remove temporary copies without exposing mutable
aliases. Splitting still creates the left leaf, then the right leaf, using the
same slices and `make_leaf` calls. Internal child IDs are copied before recursive
mutation. Hash materialization and partial-leaf queries read existing bytes;
their cached values and node allocation behavior remain the same. Public profile
checks, clamped split behavior and empty-range behavior remain intact.

## Compatibility evidence

The six new tests in `packages/rust/tests/test_edit_optimization.rs` passed
against both the captured baseline and the selected implementation. They include
an independent frozen scalar recurrence, full-width reducer boundary cases,
eight-byte block boundaries, every byte value, unaligned slices and inputs over
1 MiB. The selected implementation also passed these tests with debug overflow
checks enabled.

Five baseline edit traces pin 19 root fingerprints each plus a digest of every
allocated node, its child IDs, bytes, length, weight, repetition count and cached
hash before/after substring queries. They cover default eager/lazy modes, another
Mersenne base, profiles `(17, 2)` and `(31, 3)`, persistent edits, repeats beyond
`2^32`, split/rejoin, and lazy materialization. These are regression oracles,
alongside byte comparisons and structural validation; the digest is not a claim
of cryptographic identity.

The selected release suite passed 144 tests. Two experimental tests were
explicitly filtered: E3's 10 GB materialization and E5's sliding throughput test.
The `no_std` release library check passed. The documentation test command passed
with zero doctests. Logs are retained under `validation/`.

## Diagnostic workload and oracle review

`edit_profile.rs` uses the same RNG transitions, ASCII byte generation, 50 edit
operations, deletion/insertion choices and alternating equal/final-byte-different
targets as `benchmarks/consumer_exactness/edit_compare.rs`. ASCII makes Ropey's
character indices equivalent to the byte indices used by Hashrope. Independent
and shared-history targets are kept separate. The original untimed Hashrope edit
warmup is retained; it warms power-cache entries and leaves historical arena
nodes allocated, so this is the existing warmed workload rather than a cold
construction benchmark.

For combined comparison, the full 50-bit expected answer pattern is checked after
timing, rather than checking only a count of matches. The final edited content is
also checked against the separately maintained byte-vector model. Comparison-only
preparation retains all 50 intermediate states and verifies their bytes after
measurement. Its immutable comparison batch is repeated three times and the
three result masks are XORed. Edit-only checks the final byte state. Shared
targets are constructed using each library's edit routine, so the combined arm
does not independently flatten-check every intermediate target; avoid describing
it as doing so. The separate scalar and baseline regression tests provide the
stronger implementation-compatibility evidence.

Primitive hash and concat probes derive their expected values from the same
hasher before timing. Those checks detect unstable probe results; they are not
independent proofs of fingerprint correctness. Similarly, `fingerprint_only` is
a noncryptographic diagnostic control, not an exact equality implementation.

Timing and counter executables are separate. Instrumentation is confined to an
isolated copied core and compiled out of the timing build. The recorded timing
build includes the selected production operations, not an alternate algorithm.
Build metadata pins source hashes, binaries, locked dependencies and compiler
version; paired runs check those identities before executing the binaries.

## Allocation accounting and interpretation

The global allocator counts successful `alloc`/`alloc_zeroed` calls and their
requested byte sizes, successful `realloc` calls and their full new requested
sizes, and explicit deallocations while measurement is active. It delegates to
`System` and uses allocation-free atomic counters. The single-threaded workload
makes the relaxed atomic accounting adequate here. Temporary allocations and
deallocations inside an edit are included; final arena/Rope destruction, setup,
reporting and byte-oracle flattening occur outside the counter interval.

`new_nodes` is the increase in the arena's stored node count. Byte verification
does not add rope nodes. The separately instrumented `split_leaf_clone_bytes`
counter charges the baseline's explicit full-leaf clone and correctly charges
zero when the borrowed-leaf implementation is selected. It does not count the
two output slices as eliminated copies: those allocations remain visible in the
allocator totals. Hash byte counters measure submitted bytes, not reduction
instructions, so grouping eight bytes does not artificially reduce that count.

These counters are cumulative requested traffic, **not** peak live memory, RSS,
allocator overhead, retained capacity, or a claim that Ropey and Hashrope have
the same ownership model. Reallocation requested bytes are not net growth and
must not be added to allocation bytes and called retained memory.

Separate edit-only, comparison-only and combined runs are not additive timing
components. Retaining all intermediate comparison states changes locality and
Ropey's copy-on-write conditions. Lazy mode suppresses hash work but also changes
cache and allocation behavior. Therefore eager-minus-lazy is not a precise
percentage of runtime spent hashing, and combined-minus-comparison-only is not
a clean edit-cost estimate. The diagnostic's three seed/process summaries are
exploratory evidence; final nine-process measurements support the end-to-end
performance claims.

The standalone `edit_profile.py` runner refuses an existing final JSON but does
not refuse every partial TSV/build-metadata pathname after an interrupted run;
use a fresh output directory when reproducing that standalone diagnostic. The
paired `edit_profile_pair.py` runner explicitly refuses existing raw TSV paths
as well as final JSON paths. The current paired run uses fresh output labels and
is unaffected. This is a limitation of standalone reproduction, not evidence
that existing results were overwritten in the reviewed run.

## Reviewed source identities

Baseline archive SHA256:
`4c394e6adb99d6cd5f755e6fca139148a04777c2c468e75b715539710428f80d`.

| File | SHA256 |
| --- | --- |
| `packages/rust/src/rope.rs` | `0eae07ab1c8c4563f64e5babfe017bf6731558c7ea397248c327327cd75ed468` |
| `packages/rust/src/polynomial_hash.rs` | `5cf3dbb650ae183523eec6293ac4951f2079129ba04c78e6106c57c555ecae79` |
| `packages/rust/tests/test_edit_optimization.rs` | `9d9e228797983b4516b7246db5453006d96f7ef0900ded5eb6ab288dcf27dd48` |
| `edit_profile.rs` | `4651b3d5491dde8d98523ad228b21d0e3502f42966d83d681d28b3430658e27b` |
| `edit_profile.py` | `97c5c6c31a5f7e62e367950a847ce8a8990bf90fc13cc244d81608178a297c90` |
| `edit_profile_pair.py` | `b61dfe34ab9909f594f9dc27487ae55ebb191c66b6888fb20dc904aa1fb5db66` |
