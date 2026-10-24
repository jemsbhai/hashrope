# hashrope

Hashrope **0.3.3** for Node.js and browsers, preserving the published npm 0.2.2 fingerprint contracts and adding profile-safe metadata operations and exact whole/range comparisons. Persistent byte ropes carry polynomial fingerprints through concatenation, splitting, repetition, substring queries and streaming copies. Zero runtime dependencies. MIT licensed.

```sh
npm install hashrope
```

The starting implementation is the published npm `hashrope@0.2.2` artifact. All 46 published files matched the authoritative standalone source before integration here; see [release provenance](docs/RELEASE.md). This checkout prepares version 0.3.3; publication status is recorded separately by the release process.

```ts
import { PolynomialHash, Leaf, ropeConcat, ropeSplit, ropeHash,
  ropeSubstrHash, ropeRepeat, ropeLen, ropeToBytes } from 'hashrope';

const encode = (s: string) => new TextEncoder().encode(s);
const h = new PolynomialHash();
const rope = ropeConcat(new Leaf(encode('hello '), h), new Leaf(encode('world'), h), h);
console.assert(ropeHash(rope) === h.hash(encode('hello world')));
const [left, right] = ropeSplit(rope, 6n, h);
console.assert(ropeSubstrHash(rope, 0n, 5n, h) === h.hash(encode('hello')));
console.log(new TextDecoder().decode(ropeToBytes(right))); // world

const huge = ropeRepeat(new Leaf(encode('ab'), h), 1n << 100n, h);
console.assert(ropeLen(huge) === 2n << 100n); // no expansion
```

CommonJS consumers can use `const { PolynomialHash } = require('hashrope')`. Both module forms include TypeScript declarations. Runtime support is pinned to **Node >=20.11.1** and browsers with native **ES2022**, BigInt and Uint8Array. No Node-only modules are imported by the runtime. Browser compatibility is checked by bundling the packed package for a browser and executing its examples in the documented verification environment; this does not certify every browser version. The verified browser baseline is Chrome 155.0.8059.85. See [verification](docs/VERIFICATION.md).

## Streaming

```ts
import { SlidingWindow, PolynomialHash } from 'hashrope';
const sw = new SlidingWindow(1024n, 256n);
sw.appendBytes(new TextEncoder().encode('hello world'));
sw.appendCopy(5n, 5n); // copies "world"
console.assert(sw.finalHash() === new PolynomialHash().hash(
  new TextEncoder().encode('hello worldworld')));
```

`dMax + mMax` is the eviction threshold; after crossing it, `dMax` represented bytes remain. These parameters are not per-operation limits. `currentHash()` covers every byte ingested, including evicted bytes. `finalHash()` is an alias and does not close the stream.

## Contracts that matter

- All core positions are **byte offsets**, not Unicode character or token indices. Encode text explicitly. A split may divide a UTF-8 character.
- Hashes, lengths, weights, repetitions and counters are **bigint**. Integer arguments also accept safe integer numbers; unsafe, fractional or nonfinite numbers are rejected. Use `12345678901234567890n` for large values.
- Bytes use **Uint8Array**. Leaves copy input bytes, including Node Buffer inputs; the `data` getter returns a fresh copy. Persistent nodes are frozen and share immutable children. `null` is the empty rope.
- Nodes retain an immutable `HashProfile` (prime/base). Fingerprint operations require compatible profiles, including no-op boundaries. Exact comparisons may compare roots made with different profiles. This deliberate safety change rejects previously accepted incoherent mixed-profile construction.
- A fingerprint is **not proof of content equality** and is not a cryptographic digest. At the default base 131, `[0,131]` and `[1,0]` both hash to **263**, before reduction, under both exported moduli. Increasing the modulus does not fix that collision. `ropeEqual` and its compatibility alias `contentEquals` compare represented bytes. `ropeLcp`, `ropeRangeEqual` and `ropeRangeLcp` are exact and never use fingerprints. `structuralEquals` retains the existing structure/bytes/stored-metadata semantics.
- Exact comparisons skip identical shared ranges, cache verified subproblems, and prove periodic equality using finite byte witnesses. Independent content can still require a full scan; traversal, memoization and BigInt work allocate. There is no universal logarithmic bound.
- Repetition can represent enormous expanded lengths with one additional node. Flattening necessarily takes expanded-output time and memory. `ropeToBytes` defaults to a 64 MiB materialization budget; an explicit larger budget is possible, subject to safe indexing and the runtime's allocation limits.
- Leaves may be arbitrarily large within a typed array's limits. Splitting or partially hashing a leaf costs time proportional to its bytes. Tree path bounds alone do not establish blanket logarithmic runtime. Substring hashing creates no rope nodes, but byte slices, BigInt operations and caches may allocate. Power caches and retained earlier roots can grow independently of window capacity.

## API and parity

[API examples](docs/API.md) cover every operation. [Compatibility](docs/COMPATIBILITY.md) explains valid ranges, inherited malformed-input quirks, errors and memory limits. The [parity matrix](docs/PARITY-MATRIX.md) maps every Python public export and member. `ropeInsert`, `ropeDelete`, `contentEquals` and `structuralEquals` are explicit JavaScript conveniences; Python's published export list has no standalone insert/delete functions.

The reference is the independently downloaded PyPI wheel, SHA-256 `710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2`. Oracle generation verifies this hash and imports the wheel rather than a source checkout. The Java port's [parity matrix](https://github.com/jemsbhai/hashrope-java/blob/main/docs/PARITY-MATRIX.md) and [compatibility notes](https://github.com/jemsbhai/hashrope-java/blob/main/docs/COMPATIBILITY.md) are supporting references.

## Development and release

```sh
npm ci
npm test
npm run build
```

See [TDD evidence](docs/TDD.md), [verification](docs/VERIFICATION.md), and [release steps](docs/RELEASE.md) for reproducible oracle, consumer and tarball checks. Tests, downloaded upstream code, scratch directories, logs and credentials are excluded from the npm file allowlist. Publishing uses the separately verified tarball and existing release authorization; see the release checklist.
