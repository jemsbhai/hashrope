# Changelog

## [0.1.0] - 2026-04-07

### Added
- Initial release extracted from UHC (Unified Hash-Compression Engine).
- PolynomialHash: Mersenne-prime polynomial hashing with O(log n) power caching.
- phi: Geometric accumulator via inverse-free repeated doubling in O(log q).
- Leaf, Internal, RepeatNode: Immutable hash rope nodes (Definition 6).
- BB[2/7] weight-balanced tree with single/double rotations (Definition 7).
- rope_concat: O(k * |h_L - h_R|) balanced concatenation (Theorem 10).
- rope_split: O(k * log w) split at arbitrary byte position (Theorem 7).
- rope_repeat: O(k * log q) repetition via RepeatNode (Theorem 8).
- rope_substr_hash: O(k * log w) allocation-free substring hashing (Theorem 9).
- SlidingWindow: Bounded-memory streaming hash with eviction (Theorem 11).
- validate_rope: Invariant checker for all nine rope invariants.
- rope_from_bytes, rope_to_bytes: Construction and reconstruction helpers.
