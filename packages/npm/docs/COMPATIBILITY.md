# Compatibility and exact-comparison contracts

Version 0.3.4 retains the runtime behavior of npm 0.3.3, which extends published npm 0.2.2 while retaining its Python **hashrope 0.2.2** arithmetic oracle. The JavaScript API uses camelCase names, `bigint` integer results and `Uint8Array` bytes. Compatibility means the documented numerical, persistent-rope and streaming contracts below; it does not mean identical dynamic-language behavior, exceptions or resource limits. The [parity matrix](PARITY-MATRIX.md) audits every Python export and public member.

## Reference and provenance

The primary oracle is the independently downloaded PyPI wheel `hashrope-0.2.2-py3-none-any.whl`, SHA-256:

```text
710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2
```

The fixture generator verifies that digest and imports the wheel itself. It does not import the original Python checkout or obtain expected results from this TypeScript implementation. [PyPI metadata](https://pypi.org/pypi/hashrope/json) identified 0.2.2 as the latest stable release during this task's registry check on 2026-10-23. `PYTHON_VERSION` identifies that reference and `VERSION` identifies this package release; `PYTHON_VERSION` remains `"0.2.2"`; `VERSION` is `"0.3.4"`.

The upstream [Java parity matrix](https://github.com/jemsbhai/hashrope-java/blob/main/docs/PARITY-MATRIX.md) and [Java compatibility notes](https://github.com/jemsbhai/hashrope-java/blob/main/docs/COMPATIBILITY.md) are supporting research, not the independent oracle. MIT licensing and the upstream copyright notice are preserved in [LICENSE](../LICENSE). This port adds no claim that the paper's theoretical statements cover every host-language allocation or malformed input.

## Runtime and integer representation

The declared Node minimum is **20.11.1**. The runtime is compiled to **ES2022**, with separate ESM and CommonJS entry points and corresponding TypeScript declarations. It has zero runtime dependencies and no Node-only runtime imports. Browser use requires native ES2022 semantics, `bigint` and `Uint8Array`; BigInt is not polyfilled. Browser-target bundling and the executed verification environments are reported in [VERIFICATION.md](VERIFICATION.md). A successful bundle is not certification of every browser version. The development compiler is pinned in `package-lock.json`; declaration-consumer checks use that compiler rather than guessing an earlier TypeScript compatibility floor.

Choose **one module flavor and one installed package instance per node graph**. The ESM and CommonJS builds have separate class constructors and private ownership maps. Their nodes are not interchangeable: an ESM `Leaf` is not an instance of the CommonJS `Leaf`, and mixed-build operations may throw `TypeError`. Sharing bytes and rebuilding a node through the receiving API is supported; sharing node objects across these package instances is not. This does not affect using either build independently.

Hashes, primes, bases, powers, lengths, weights, repetition counts, byte positions and streaming counters are `bigint`. Integer arguments use the additive `Integer = bigint | number` type. A `number` is accepted only when `Number.isSafeInteger` is true, then converted exactly. Fractions, infinities, NaN and numbers outside the safe integer range throw `RangeError`; strings, booleans, null and other types throw `TypeError`. Construct a BigInt literal or explicitly parse a string before calling the API. A number already rounded by JavaScript cannot be recovered by a library conversion.

Compressed logical lengths beyond `Number.MAX_SAFE_INTEGER` remain exact. They are not converted to Number for rope arithmetic, splitting, repetition or substring addressing. `ropeHeight` is a Number describing actual structural depth, which is separately limited by traversal stack and available memory. BigInt and typed arrays still have finite engine and memory limits; arbitrary precision is not unlimited resources or a promise of Python's exact limits. JSON does not serialize BigInt by default: consumers should deliberately encode decimal strings when needed.

## Byte semantics and ownership

Bytes are unsigned: `0x00` contributes coefficient 1 and `0xff` contributes coefficient 256. All core indices count **bytes**, not Unicode characters, code points, grapheme clusters or tokens. Encode strings explicitly, for example with `TextEncoder`; a valid split can cut through a UTF-8 multibyte character. There is no normalization or implicit string encoding.

Byte entry points require a `Uint8Array` recognized by the local realm, including Node Buffer instances. A typed array from another JavaScript realm can fail the `instanceof Uint8Array` check; copy it with the receiving realm's `new Uint8Array(foreignBytes)` before passing it in. Unlike Python's incidental acceptance of integer sequences or some falsy values, plain arrays, strings, `Int8Array`, null and undefined are rejected. In particular, `appendBytes(null)` is a `TypeError`, not Python's falsy-input no-op. An empty `Uint8Array` is accepted by `hash`, `ropeFromBytes` and `appendBytes`; it produces zero, null and a no-op respectively. A direct `Leaf` cannot be empty.

A leaf copies its incoming bytes and its `data` getter returns a fresh `Uint8Array`. Buffer inputs are copied using typed-array construction, not Buffer's view-producing `slice`. The bytes passed to an overridable hashing callback do not expose the retained private leaf buffer. Public node state is frozen, with immutable shared child references; TypeScript readonly annotations are accompanied by runtime freezing. A custom hasher may intentionally return a nonstandard fingerprint, but it cannot obtain the retained leaf storage through that callback.

`null` is the empty rope. Concatenation with an empty side, repetition by one and endpoint splits preserve the existing root reference. Splitting at an internal boundary returns the children. Prior roots remain valid after later persistent operations. Mutations through monkeypatching, reflection, forged nodes, or arbitrary Python-style node subclasses are outside this API contract.

## Fingerprints, profiles and equality

The default profile is prime `MERSENNE_61 = 2^61-1`, base 131. `MERSENNE_127 = 2^127-1` is also exported; changing the modulus does not change the default base. Use the same prime/base pair for every node and operation in a rope. Since npm 0.3.3, nodes and hashers expose immutable `profile` metadata. Direct constructors, fingerprint queries and persistent edits reject mismatched profiles, including no-op and cached-hash boundaries. `validateRope(node,h?)` checks child profiles and a supplied hasher. Separate hashers with the same prime/base remain compatible. Stored hashes are preserved rather than silently recomputed with the most recently supplied hasher. Algebraic fingerprint comparisons are meaningful only under compatible profiles and the documented arithmetic domain.

**Fingerprints do not prove content equality and are not cryptographic digests.** The default byte coefficients have a deterministic collision before any modular reduction:

```text
H([0, 131]) = 1 * 131 + 132 = 263
H([1,   0]) = 2 * 131 +   1 = 263
```

Both exported moduli therefore produce the same value for these different equal-length inputs. A larger modulus does not remove this collision, provide authentication, or make adversarial comparisons safe. Explicit tests retain this behavior for parity. The exact APIs never inspect hashes, even to reject unequal hashes: accepted legacy non-Mersenne profiles can give different fingerprints to equal bytes in different tree shapes.

JavaScript `===` compares node identity. The additive `structuralEquals(a,b)` compares node kinds, retained bytes, child structure, repetitions and stored metadata, corresponding to Python dataclass value equality. Different structures representing the same bytes can compare structurally unequal. Stored metadata participates, so same-content nodes made with different profiles can also be structurally unequal. Python object hash values and repr text are not reproduced.

`ropeEqual(a,b)` and the compatibility alias `contentEquals(a,b)` compare actual represented bytes, independent of fingerprints and profiles. `ropeLcp(a,b)` returns their exact byte LCP as a BigInt. `ropeRangeEqual(a,aStart,b,bStart,length)` and `ropeRangeLcp` check equally sized ranges without splitting or rehashing. Negative or out-of-bounds ranges throw `RangeError`, including empty ranges whose start exceeds the rope length; bounds precede identity shortcuts.

The engine verifies owned leaf bytes without exposing/copying buffers, skips identical shared ranges, memoizes already verified subproblems, and uses the Fine-Wilf witness of length `p+q-gcd(p,q)` for two periodic streams. This handles large repeated strings and independent shared DAGs without blindly expanding them. It is not a universal logarithmic or compressed-size bound: unshared bytes may require a full scan, traversal and memoization consume storage, and BigInt costs depend on operand size. It does not mutate ropes or hash caches.

## Polynomial arithmetic

`PolynomialHash(prime,base)` checks the Python interval `2 <= base < prime`. It does not prove primality or Mersenne shape. For positive Mersenne moduli and nonnegative operands, `mersenneMod` implements modular reduction. For accepted positive non-Mersenne masks and negative operands, it deliberately reproduces Python's bit-folding algorithm; this need not equal the mathematical remainder. Algebraic guarantees outside their mathematical assumptions should not be inferred from acceptance of an input.

`power(n)` caches modular powers and supports negative exponents through modular inversion. An unavailable inverse throws `RangeError`. Cache growth follows distinct requested exponents and is not bounded by window capacity. `phi` uses the inverse-free geometric algorithm and preserves the release's bit-scan behavior for signed counts; negative counts do not describe a mathematical negative-length sequence. `hashOverlap` uses Python floor-division/remainder semantics for signed dimensions. A zero period raises `RangeError`, even when the requested length is zero. The caller-provided prefix hash is not independently verified.

**Termination boundary:** nonpositive helper moduli are rejected with `RangeError` rather than reproducing potentially nonterminating Python folding. The `phi` zero/one short circuits are retained after JavaScript integer-argument checks. Python nonterminating cases are not run as oracle fixtures. Runtime-specific division and type-error wording is not a cross-language compatibility promise.

`PolynomialHash` and `SlidingWindow` remain extensible. Ordinary method overrides are observable. In particular, RepeatNode metadata follows the reference's direct `power` / `phi` / multiply sequence rather than delegating to `hashRepeat`, and `finalHash` calls the overridable `currentHash`. Subclass authors remain responsible for the semantics of their returned hashes and for mutable state.

## Rope operations, balancing and ranges

`ropeFromBytes` constructs one leaf for nonempty input. A leaf's weight is one regardless of byte count. An internal node's weight is the sum of child weights; a repeat's weight is child weight times repetition count. Weight counts **expanded leaves**, not bytes or allocated node objects. `ALPHA_NUM = 2n` and `ALPHA_DEN = 7n` describe the BB[2/7] balance condition. For an internal node with total weight greater than two:

```text
2 * totalWeight <= 7 * leftWeight <= 5 * totalWeight
```

`ropeConcat` uses the reference's persistent balanced join, including repetition-count decomposition and checked rotations. Direct construction of `Internal` combines metadata without rebalancing, so a manually constructed tree may fail validation. Repeats are not expanded when their hash or logical size is calculated; nested repeat structure is preserved.

`ropeSplit(node,pos,h)` clamps endpoint positions: `pos <= 0` gives `[null,node]`; `pos >= ropeLen(node)` gives `[node,null]`. Interior repeat splits handle both whole-period boundaries and offsets inside a period. It returns a JavaScript tuple rather than a Python tuple object.

`ropeRepeat(node,0,h)` is null and count one returns the input. Other counts use RepeatNode, whose count must be at least two. Repeating null zero or one times is null; repeating null two or more times raises a mapped type error. `ropeHeight` returns zero for null and leaves, child height plus one for a repeat, and maximum child height plus one for an internal.

For valid ranges `0 <= start`, `0 <= length`, `start + length <= ropeLen(node)`, `ropeSubstrHash` returns the exact interval fingerprint without constructing new rope nodes. Complete coverage returns stored metadata. Partial repeated intervals combine a tail, whole periods and a head.

The pinned Python release does not consistently reject malformed substring bounds. The port retains its shape-dependent behavior rather than presenting it as a general clipping API. Leaf queries use Python-style slicing; repeat queries may hash periodic bytes beyond the represented length; internal queries can compose hashes using requested rather than available lengths. For example, querying start 6, length 1 on `Leaf("ababab")` yields zero, while the same query on `RepeatNode(Leaf("ab"),3)` yields the hash of `"a"` (98). Negative lengths can also produce nonzero hashes. Applications should use valid ranges. Empty input or zero requested length follows the reference short circuit.

`ropeInsert(node,pos,value,h)` and `ropeDelete(node,start,length,h)` are additive JavaScript editing helpers, built from split and concatenation. Insertion clamps its byte position using `ropeSplit`. Deletion clamps start, rejects a negative length, and removes up to the available remaining bytes; deleting beyond the end does not add data. Both preserve prior roots. Python 0.2.2 exports no standalone insert/delete operation. See [API.md](API.md) for examples.

## Validation

`validateRope(node,h?)` checks the Python structural invariants: leaf length/weight; internal length, weight and balance; repeat count, length and weight; and descendants. The optional hasher must match the root profile; child profiles are also checked. It **does not recompute hashes** or prove content equality. Independent unit/property/oracle tests check numerical hashes in addition to the compatibility validator.

A violation throws the exported `RopeValidationError`, whose `name` is `"AssertionError"`. This is an ordinary runtime check and is not disabled by a Python-style optimized execution flag. The error class is a JavaScript adaptation, not Python exception inheritance.

## Sliding-window semantics

Defaults are `dMax=32768n` and `mMax=258n`; `W=dMax+mMax`. Construction requires positive capacities, checking dMax first. `appendBytes` appends one owned literal leaf. `appendCopy(offset,length)` can copy a nonoverlapping source interval or repeat a source period for an overlapping back-reference.

Copy validation checks offset against total decoded `pos` before retained `windowLen`, including for zero-length requests. `(0,0)` is a no-op, while zero offset with positive length triggers a division `RangeError`. dMax and mMax are not per-operation caps: a distance above dMax is accepted while retained, and a copy can exceed mMax. A literal append can exceed W.

Eviction happens strictly above W. At exactly W, no bytes are evicted. Crossing the threshold folds the removed prefix into `hPrefix` and retains dMax represented bytes. `currentHash` covers the full stream, including evicted bytes. `finalHash` is an alias, not a finalization barrier; later appends are allowed.

The Python public mutable slots map to `h`, `dMax`, `mMax`, `W`, `hPrefix`, `lPrefix`, `rWindow` and `pos`. `windowLen` is read-only derived state. Direct assignments do not repair consistency or recompute other fields; assign BigInts to integer state. The constructor captures the stream profile: replacing `h` with a different profile or attaching an incompatible `rWindow` is rejected before hash/append operations, including empty appends. A new stream is required to change profiles. The returned root is an immutable snapshot, while SlidingWindow itself is mutable.

The reference's negative-copy quirks are retained for integer inputs. They can make `pos` inconsistent with represented bytes: after appending `abc`, offset 1 and length -1 reduce pos to two while retaining three bytes. Such inputs are compatibility behavior, not valid decoded-stream operations. Under valid operations and unmodified state, `pos-lPrefix == windowLen`, the retained represented length is at most W, and currentHash matches the decoded byte stream.

## Resource boundaries and complexity

`ropeToBytes(node,maxBytes=64*1024*1024)` materializes a fresh array with a default **64 MiB** budget. An explicit integer budget between zero and `Number.MAX_SAFE_INTEGER` can permit more, but the output length must also fit that budget and the engine's typed-array/allocation limits. A budget above the safe integer limit is rejected even when the rope itself is small. Oversize logical output is rejected before Number conversion and allocation. An allowed size can still exhaust memory. This explicit budget is a deliberate JavaScript difference from Python's attempt to allocate its full bytes result. There is no truncation. A supplied budget does not preallocate that amount or restrict compressed logical lengths.

Tree path operation counts do not imply fixed-time arithmetic or allocation:

- Balanced join and split traverse/rebuild logarithmic paths in expanded leaf weight under their preconditions. Splitting a large leaf additionally copies and hashes its byte slices. Python and this port impose no small fixed leaf-size limit.
- Repetition adds one node and performs a number of modular operations proportional to the repetition count's bit length, plus power computation. BigInt multiplication and exponentiation costs depend on operand bit lengths.
- A substring query constructs no new rope nodes, but a partial-leaf query obtains an owned copy of that leaf before hashing the slice, and arithmetic/cache objects can allocate. "No rope nodes" does not mean allocation-free JavaScript execution.
- Flattening necessarily requires expanded-output time and memory. Exact content comparison may require expanded-length traversal. Structural height and validation traverse child references without a visited-node cache; a shared subtree may be visited more than once. Manually deep trees remain limited by the JavaScript stack.
- Window eviction bounds retained represented bytes under valid operations, not total process memory. Power caches, allocated arithmetic temporaries, persistent roots retained by callers and transient copies can exceed that bound.

## Error and language mapping

| Python behavior | JavaScript behavior |
|---|---|
| Explicit library `ValueError`, invalid numeric ranges, noninvertible negative power | `RangeError` with stable library messages where mapped |
| `ZeroDivisionError` | `RangeError` |
| Wrong dynamic types or invalid/null child access | `TypeError`; engine-generated wording is not promised |
| Validation `AssertionError` | `RopeValidationError`, named `AssertionError` |
| Arbitrary-precision integer | `bigint`; lossless safe-Number argument convenience |
| `bytes` and incidental duck-typed inputs | Explicit `Uint8Array` with defensive ownership |
| Frozen dataclasses and Optional node union | Frozen node classes, TypeScript `Node`, null empty rope |
| Dataclass value equality | Explicit `structuralEquals`; `===` remains identity |
| Python repr / object hash / serialization | Not reproduced; use declared data and hash accessors |
| Python allocation/recursion limits | Documented materialization budget and host resource limits |

See [VERIFICATION.md](VERIFICATION.md) for executed checks and remaining limits. Neither a matrix entry nor the existence of an API is by itself a claim that a test has passed.
