# Python 0.2.2 to JavaScript parity matrix

This audit covers every name in the pinned wheel's `hashrope.__all__`, public class constructors, fields/properties and methods, module version, and rope balance constants. The JavaScript package root is the supported import surface. Type-only exports are erased at runtime. Deliberate language and resource differences are part of the contract in [COMPATIBILITY.md](COMPATIBILITY.md).

The primary expected-value source is the independently downloaded, SHA-256-pinned PyPI wheel, not the original checkout or Java implementation. The Java [matrix](https://github.com/jemsbhai/hashrope-java/blob/main/docs/PARITY-MATRIX.md) is supporting research. Fixture counts, actual pass/fail totals, compiler/runtime versions and package-consumer results belong to [VERIFICATION.md](VERIFICATION.md); this table identifies coverage and does not substitute for execution evidence.

## Evidence groups

| Group | Source / purpose |
|---|---|
| Arithmetic | `test/polynomial.test.mjs`: constants, exact integer arithmetic, profiles, signed inputs, errors, collisions |
| Nodes | `test/nodes.test.mjs`: metadata, ownership, structural equality, validation, materialization limits |
| Rope | `test/rope.test.mjs`: persistent balanced joins, splits and repeated-child decomposition |
| Query | `test/rope-query.test.mjs`: repetition edge cases, all valid small ranges, malformed ranges, huge intervals and metadata fast paths |
| Edits/equality | `test/edits-equality.test.mjs`: clamped insert/delete, exact byte comparison despite collisions, profile/shape differences and deterministic edit models |
| Review | `test/review.test.mjs`: actual typed-array storage, detached copies, retained callback isolation, profile/shape equality and 150-bit ranges |
| Streaming | `test/sliding.test.mjs`: complete state, copies, eviction, errors, inherited quirks and huge logical sizes |
| Oracle | `test/oracle.test.mjs`: verifies `test/fixtures/manifest.json` and executes every published-wheel record in `test/fixtures/oracle.json`; hashes, shape digests, metadata and bytes are compared |
| Package | `test/package-contract.test.mjs`, `scripts/pack.mjs` and `scripts/consumers.mjs`: allowlist, clean tarball installations, ESM/CommonJS imports and TypeScript declarations |
| Browser | `scripts/consumers.mjs`: browser-target build and execution environment described in VERIFICATION; no Node-only runtime imports |
| TDD | [TDD.md](TDD.md) and `evidence/`: recorded red, green and refactor runs |

## Every top-level Python export

| Python export | JavaScript mapping | Contract and coverage |
|---|---|---|
| `PolynomialHash` | `PolynomialHash` | Exact BigInt arithmetic; methods below. Arithmetic, Oracle |
| `phi` | `phi` | Inverse-free geometric accumulator; signed-count behavior and termination boundary documented. Arithmetic, Oracle |
| `mersenne_mod` | `mersenneMod` | Positive-modulus Python folding, including signed operands and accepted non-Mersenne masks. Arithmetic, Oracle |
| `mersenne_mul` | `mersenneMul` | Exact product followed by reference folding. Arithmetic, Oracle |
| `MERSENNE_61` | `MERSENNE_61` | BigInt `2^61-1`. Arithmetic |
| `MERSENNE_127` | `MERSENNE_127` | BigInt `2^127-1`; does not imply cryptographic security. Arithmetic |
| `Leaf` | `Leaf` | Nonempty owned bytes; length/hash metadata; weight one. Nodes, Oracle |
| `Internal` | `Internal` | Shared immutable children, direct metadata join; constructor does not rebalance. Nodes, Oracle |
| `RepeatNode` | `RepeatNode` | Shared child, reps >= 2, exact expanded metadata, no expansion. Nodes, Oracle |
| `Node` | `Node` (TypeScript type) | Node union including null; no runtime Python typing object. Package |
| `rope_concat` | `ropeConcat` | Persistent BB[2/7] balanced join with null identity. Rope, Oracle |
| `rope_split` | `ropeSplit` | Byte offsets, endpoint clamping, repeated-period splitting, `[left,right]` tuple. Rope, Oracle |
| `rope_repeat` | `ropeRepeat` | Count zero -> null; one -> identity; otherwise RepeatNode. Query, Oracle |
| `rope_substr_hash` | `ropeSubstrHash` | Exact valid byte intervals, metadata fast paths, inherited malformed-range behavior. Query, Oracle |
| `rope_len` | `ropeLen` | Stored BigInt logical byte length; null -> 0n. Nodes, Oracle |
| `rope_hash` | `ropeHash` | Stored BigInt fingerprint; null -> 0n; no rehash. Nodes, Oracle |
| `rope_height` | `ropeHeight` | Actual structural height as Number; null/leaf zero, repeat child+1, internal max+1. Nodes, Oracle |
| `rope_from_bytes` | `ropeFromBytes` | Empty Uint8Array -> null; nonempty -> one owned leaf. Nodes, Rope |
| `rope_to_bytes` | `ropeToBytes` | Fresh expanded Uint8Array; deliberate 64 MiB default budget, checked override and runtime limits. Nodes, Oracle |
| `validate_rope` | `validateRope` | Structural and profile checks; optional hasher must match; no rehashing. Deliberate 0.3.3 safety change. Nodes, Rope, Oracle |
| `SlidingWindow` | `SlidingWindow` | Mutable stream state; methods and fields below. Streaming, Oracle |

## Version and balance metadata

| Python name | JavaScript mapping | Contract |
|---|---|---|
| `hashrope.__version__` | `PYTHON_VERSION` | Pinned reference version `"0.2.2"` |
| Package release version | `VERSION` and package.json `version` | JavaScript release currently `"0.3.4"`; separate ecosystem release identity |
| `hashrope.rope.ALPHA_NUM` | `ALPHA_NUM` | BigInt `2n` |
| `hashrope.rope.ALPHA_DEN` | `ALPHA_DEN` | BigInt `7n` |

## PolynomialHash constructor, properties and methods

| Python member | JavaScript mapping | Contract / cases |
|---|---|---|
| `PolynomialHash(prime=MERSENNE_61,base=131)` | `new PolynomialHash(prime?,base?)` | Exact interval check; no primality/Mersenne proof; safe Number arguments additive |
| `.prime` | `.prime` | Read-only BigInt modulus |
| `.base` | `.base` | Read-only BigInt base |
| `.power(n)` | `.power(n)` | Cached modular exponentiation; huge/negative exponents and noninvertible-base error |
| `.hash(data)` | `.hash(data)` | Uint8Array bytes, byte+1 coefficients, empty zero, all unsigned bytes |
| `.hash_concat(h_a,len_b,h_b)` | `.hashConcat(hA,lenB,hB)` | Composable hash, exact large/signed lengths |
| `.hash_repeat(h_s,d,q)` | `.hashRepeat(hS,d,q)` | Whole-period hash; zero/one/huge counts, alpha=1 degeneracy |
| `.hash_overlap(h_p,d,l,h_prefix)` | `.hashOverlap(hP,d,length,hPrefix)` | Signed floor division, whole periods plus supplied prefix; zero-period error |
| Method overrides | Ordinary JavaScript class overrides | Dispatch retained; custom implementations own their semantics |
| `_power_cache`, `_p`, `_x` | Private implementation state | Not public APIs; internal cache identity/order not matched |

Arithmetic tests and oracle records cover these methods. The two default-base collision inputs `[0,131]` and `[1,0]` each produce `263n` under both exported moduli. The test is a regression requirement, not an unexpected failure.

## Node fields, ownership and equality

| Python field/behavior | JavaScript mapping | Contract / coverage |
|---|---|---|
| `Leaf.data` | `Leaf.data` getter | Fresh owned Uint8Array on each read; Buffer input and callback mutation isolation. Nodes |
| Every node `.len` | `.len` | Exact BigInt expanded byte length. Nodes, Oracle |
| Every node `.weight` | `.weight` | Exact BigInt expanded leaf count, not byte count. Nodes, Rope, Oracle |
| Every node `.hash_val` | `.hashVal` | Stored BigInt fingerprint; collisions are permitted. Nodes, Oracle |
| `Internal.left` / `.right` | `.left` / `.right` | Immutable shared child references. Nodes, Rope |
| `RepeatNode.child` | `.child` | Immutable shared child, including nested repeats. Nodes, Rope |
| `RepeatNode.reps` | `.reps` | BigInt repetition count, >=2; no safe-integer truncation. Nodes, Rope |
| Frozen dataclass attributes | Frozen node instances and readonly TypeScript state | Native strict-mode write failure; private byte ownership. Nodes |
| Dataclass equality | `structuralEquals(a,b)` | Explicit additive helper compares structure, bytes and stored metadata; not hash-only. Nodes |
| Python object hash/repr | No exact equivalent | JavaScript identity and declared accessors; no process-specific hash/repr parity claim |
| Split tuple | Readonly TypeScript tuple `[Node,Node]` | Destructurable JavaScript array; no Python tuple runtime identity |
| Optional node union | `Node`, `NonEmptyNode` types | Null empty value; no arbitrary Python duck-typed node contract |

`contentEquals` is a separate JavaScript helper for exact represented-byte comparison, including collision pairs and same bytes with differing structure or hash profiles. `test/edits-equality.test.mjs` covers these cases. It aliases ropeEqual: identical shared ranges, already verified ranges, and finite periodic witnesses can skip traversal. It never substitutes fingerprint agreement for byte equality and may require expanded-length work.

## SlidingWindow constructor, properties and methods

| Python member | JavaScript mapping | Contract / coverage |
|---|---|---|
| `SlidingWindow(d_max=32768,m_max=258,prime=MERSENNE_61,base=131)` | `new SlidingWindow(dMax?,mMax?,prime?,base?)` | BigInt defaults; positive capacities checked in Python order; custom profiles. Streaming, Oracle |
| `.h` | `.h` | Mutable hasher field |
| `.d_max` | `.dMax` | Mutable BigInt retention target |
| `.m_max` | `.mMax` | Mutable BigInt capacity component |
| `.W` | `.W` | Mutable BigInt eviction threshold; initially dMax+mMax |
| `.h_prefix` | `.hPrefix` | Mutable BigInt fingerprint of evicted prefix |
| `.l_prefix` | `.lPrefix` | Mutable BigInt prefix length |
| `.r_window` | `.rWindow` | Mutable root field; root value itself is persistent/immutable |
| `.pos` | `.pos` | Mutable BigInt decoded counter |
| `.window_len` | `.windowLen` getter | Derived BigInt retained logical length |
| `.current_hash()` | `.currentHash()` | Full-stream fingerprint from prefix and window |
| `.final_hash()` | `.finalHash()` | Calls currentHash; no close/freeze barrier, override dispatch retained |
| `.append_bytes(data)` | `.appendBytes(data)` | Owned literal leaf, empty Uint8Array no-op; null is deliberately TypeError |
| `.append_copy(offset,length)` | `.appendCopy(offset,length)` | Nonoverlap and periodic overlap, exact counters, bounds/errors and signed quirks |
| `_evict()` | Internal implementation | Observable threshold strictly above W, retention to dMax, prefix folding |
| Assignable slots | Corresponding public mutable fields | No automatic consistency repair; integer fields require BigInt state |

Streaming tests and oracle snapshots cover the complete observable stream state, not only final hashes. dMax and mMax do not impose per-copy bounds. Huge compressed copies must be checked through metadata and limited substring queries without flattening the logical output.

## Explicit adaptations and limits

| Area | Decision and verification |
|---|---|
| JS numeric inputs | Safe integer Number convenience; all unsafe/nonfinite/fractional numbers rejected, no silent rounding. Arithmetic, Nodes, Rope, Streaming |
| Byte types | Uint8Array only; Buffer copied; null/string/array inputs rejected. Arithmetic, Nodes, Streaming |
| Hash-profile mixing | Deliberate 0.3.3 change: immutable profile metadata and mismatch rejection for fingerprint operations; exact comparisons allow different profiles. Nodes, Rope |
| Structural validation | Does not prove hash correctness; independent hash assertions supplement it. Nodes, Oracle |
| Malformed substring ranges | Python shape-dependent behavior retained; documented valid range for application use. Query, Oracle |
| Negative copy quirks | Retained integer behavior can make pos inconsistent; documented as outside valid stream semantics. Streaming, Oracle |
| Nonpositive helper modulus | Explicit RangeError avoids possible Python nontermination; short circuits documented. Arithmetic |
| Materialization | Default 64 MiB budget; exact length check before Number conversion, configurable limit and engine allocation bounds. Nodes |
| > Number.MAX_SAFE_INTEGER | Logical lengths, weights, offsets, counts and counters stay BigInt; compressed tests avoid materialization. Nodes, Rope, Streaming |
| Equality / LCP | Exact whole/range LCP and equality; structural helper remains separate. Collision, custom-profile and compressed-DAG regressions. Nodes, Edits/equality |
| Insertion/deletion | `ropeInsert` clamps its position; `ropeDelete` clamps start, rejects negative length and removes up to the remaining length. Both compose split/join and are not Python exports. Edits/equality |
| Exception types | RangeError/TypeError/RopeValidationError mapping; Python runtime diagnostics not promised. All unit groups |
| Runtime limits | BigInt/array/heap/stack bounds differ from Python; exhaustion is not numerical truncation. Documented, not deliberate resource exhaustion tests |
| Complexity claims | Byte-sized leaf costs, BigInt arithmetic, cache growth, retained snapshots and expansion costs explicitly qualified. Source review and targeted behavior tests |
| ESM/CommonJS | Separate entry points and declarations from the tarball, tested in clean consumers. Constructor/private-state identities differ; keep one module flavor and package instance per graph. Package |
| Browser | ES2022 build, zero runtime dependencies and no Node-only runtime imports; executed environment listed separately. Browser |
| Serialization / monkeypatching | Python repr/hash/pickle, forged nodes and arbitrary dynamic object behavior are not port contracts |

Private Python helpers such as `_join`, `_balance`, `_decompose` and `_make_repeat` are implementation details. Their observable balancing and repetition behavior is exercised through public operations rather than exposed as additional package APIs. Internal module exports such as integer-adaptation helpers are not package-root public exports.

The 0.3.3 additions `HashProfile`, `ropeEqual`, `ropeLcp`, `ropeRangeEqual` and `ropeRangeLcp` are covered by `test/exact.test.mjs`, including installed ESM/CJS consumers. The pinned published Python 0.2.2 oracle remains unchanged.
