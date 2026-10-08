# Rust parity with published Python

Audit baseline: repository commit `100749d1cef9603465f4600b814c76d2ad6e9a6e`
(verified current `master` on 2026-10-08), PyPI `hashrope 0.2.2`, and crates.io
`hashrope 0.3.1`. Version numbers alone were not used as evidence of parity.
This document records the compatible Rust 0.3.2 release against those baselines.

The default Rust API keeps its existing types and default-profile byte hashes. Arbitrary
precision is additive through `features = ["bigint"]` and the `bigint` module.
No Python, npm, or downstream implementation is modified by this change.

## Feature and behavior matrix

| Feature or contract | Published Python 0.2.2 | Rust candidate | Status and compatibility choice |
| --- | --- | --- | --- |
| Default profile | `2^61-1`, base 131, Horner hash with byte+1 | Identical for both APIs | Already supported; frozen by published golden vectors |
| Empty bytes | Hash 0, empty rope `None` | Hash 0, empty `Node=None` | Already supported |
| Binary and Unicode representation | Raw bytes; no text normalization | `&[u8]`, explicit `str::as_bytes()` for UTF-8 | Already supported; fixtures include all 256 bytes, invalid UTF-8, NUL, emoji, combining marks and splits inside code points |
| Custom hash parameters | Mersenne modulus, base in `[2,p)` | Existing `u64` API; arbitrary-size `bigint::PolynomialHash` | Wide profiles implemented. `MERSENNE_127: u128` added at crate root and re-exported in bigint module |
| Moduli beyond 127 bits | Python arbitrary integers | Optional `BigUint` arithmetic | Implemented; 521-bit profile and 1200-bit reduction inputs in fixture |
| `mersenne_mod`, `mersenne_mul`, `phi` | Arbitrary integer arithmetic | Existing fixed-width functions plus bigint equivalents | Fixed small-modulus/full-u128 reduction; tested zero/one alpha and huge q |
| `power`, concat/repeat/overlap hash algebra | Arbitrary exponents/counts | Existing u64 methods plus bigint equivalents | Fixed 32-bit cache truncation; wide exponents/counts implemented |
| Byte construction and reconstruction | `rope_from_bytes`, `rope_to_bytes` | Arena methods and equivalent wrappers in both modules | Already supported for u64, implemented for bigint. Reconstruction remains limited by real memory |
| Leaf chunking | `rope_from_bytes` creates one Leaf; consumers can chunk inputs before concat | Existing default 512-byte arena construction; bigint also offers `from_bytes_chunked` | Rust-specific capability preserved. Tree shape is not a byte/hash contract |
| Persistent concat, split, repeat | Immutable nodes with sharing | Immutable arena nodes with stable handles | Already supported; wide equivalents implemented. No mutation of earlier roots |
| Insert, delete, replace | Compose split/concat primitives | Same composition | Already supported; edit trace checks bytes/hash and old snapshots |
| Slicing and indexing | Byte `split(pos)`, `substr_hash(start,length)` | Same valid byte ranges in both modules | Already supported; wide equivalent implemented. No character/grapheme indexing implied |
| Search | No public find/contains API | No new search API | Absent in Python; no artificial parity requirement |
| Equality | Frozen dataclass structural equality; no dedicated content-equality API | Arena-local NodeId identity; compare bytes for exact content | Intentionally language-specific. Equal fingerprints alone do not prove equality |
| Serialization | No core wire format; bytes can be reconstructed | No new core wire format or serde contract | Already applicable via bytes. Applications own their persistence formats; big integers should not be truncated to floating-point JSON numbers |
| Public nodes | `Leaf`, `Internal`, `RepeatNode` dataclasses | Existing `NodeInner::{Leaf,Internal,Repeat}` unchanged; separate bigint enum | Language-specific ownership. No variants/fields/types added to existing enum |
| Lazy hashing | No lazy API | Existing `new_lazy`, `ensure_hash`, `LAZY_SENTINEL` unchanged | Rust-specific feature preserved. `hash(&self)` deliberately returns sentinel before materialization; bigint is eager |
| Length/weight limits | Arbitrary nonnegative integers | Checked u64 old API; `BigUint` new API | Implemented wide metadata. Fixed overflow and u128 balance comparisons; arena node IDs stay u32 |
| Streaming append/copy and prefix hash | `SlidingWindow` with eviction, overlap and non-overlap copies | Same hash behavior in both modules | Already supported narrow; wide implementation added and golden-tested |
| Streaming storage reclamation | Unreachable nodes collectible | New opt-in `new_bounded` constructors copy retained graph after eviction | Implemented without changing old stable-handle contract. `arena_mut()` disables reclamation; legacy constructors retain evicted allocations |
| `rope_height`, `validate_rope` | Structural inspection | Methods and free wrappers | Already supported, wide equivalents implemented. Heights/weights need not match across different chunk sizes |
| Invalid input behavior | Python exceptions; inconsistent malformed substring ranges across node types | Existing Rust behavior retained except deterministic representational-overflow checks; new bigint methods validate documented ranges | Intentionally not standardized. Fixtures use valid inputs only; non-Mersenne constructor quirks are outside the documented contract |
| Dependency/runtime surface | Pure Python, no dependencies | Default Rust zero runtime dependencies; bigint uses optional num-bigint/num-traits | `no_std + alloc` preserved; no framework integration dependency added |

## Published reference provenance

- Python wheel SHA-256:
  `710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2`.
  Registry: <https://pypi.org/pypi/hashrope/0.2.2/json>.
  `__init__.py` and `polynomial_hash.py` match repository bytes; sliding source
  matches after CRLF normalization. Repository rope changes since the wheel
  add byte-index/UTF-8 docstrings, not runtime changes.
- Rust crate checksum:
  `c4dcb7477fcea9253c7142b950db236ca20f6d2bb6d43c259b030ca7861ea884`.
  Registry: <https://crates.io/api/v1/crates/hashrope/0.3.1>.
  Published VCS commit: `98a9209475c10ea4cd00b610a6456f16950ff6b0`.
- `tests/fixtures/python-0.2.2.json` was generated by the exact downloaded
  wheel, not by candidate Rust calculations. Integers are decimal strings and
  bytes are hex, making the fixture reusable across language implementations.
  It contains eight profiles: moduli 3/7/31, 61-bit with bases 131/257, 127-bit
  with bases 131/p-2, and 521-bit with base 131. Existing Rust tests exercise
  all representable values; bigint tests exercise the complete fixture.
- `scripts/generate_python_vectors.py --wheel <downloaded-wheel> --check`
  reproduces and verifies the fixture with the extracted wheel on `PYTHONPATH`.
  The script verifies the wheel checksum, version, and imported source bytes.
  Its output covers hashes, powers, geometric sums, overlap hashes, byte splits,
  substring ranges, repeats, persistent edits and sliding streams. The fixture
  schema is a test format, not a newly imposed serialization contract.

## Verified downstream contracts

| Consumer | Verified contract | Result / adoption constraint |
| --- | --- | --- |
| cdh-sort | `hashrope="0.3"`; public Arena/Node/NodeId; exhaustive NodeInner fields and variants; eager and lazy construction | Candidate passes unchanged consumer unit, integration and property suites. Exhaustive matching compiles. Existing enum layout and lazy sentinel semantics preserved |
| hashrope-bio Rust | `hashrope="0.2.1"` plus stale path `../../../../packages/rust`; published dependency also excludes 0.3.x | Existing downstream declaration cannot adopt this candidate. Unchanged library source passes under an audit-only manifest pointing to candidate. No downstream dependency was broadened. A preexisting doctest borrows the same arena mutably twice and fails on both published 0.3.1 and candidate |
| hashrope-bio Python | `PolynomialHash`, byte split/substr start+length, custom chunking, repeated sequence operations | Covered by core byte/range/repeat fixtures; Python implementation was not edited or replaced with Rust |
| Pollard Python | Six Hashrope imports: PolynomialHash, rope_concat/from_bytes/hash/to_bytes, validate_rope. Hashes exact persisted UTF-8 JSONL | Published Pollard 1.6.0 + Hashrope 0.2.2 golden log retained byte-for-byte; candidate Rust incremental, rebuilt, chunked, and streaming hashes match old wheel output |
| Pollard Rust / npm | No direct Hashrope dependency in inspected manifests | No package adoption implied; Python persistence contract is the confirmed integration |
| wattwarden | `pollard[openai,mcp,hashrope,tokenmaster]>=1.5.1` | Indirect dependency verified; no end-to-end app/provider integration claimed |
| Document-diff / RAG examples | Repository design documents describe planned use cases | No implemented framework integration verified; no speculative adapter or compatibility claim added |

The Pollard fixture is `tests/fixtures/pollard-1.6.0-hashrope-0.2.2.jsonl`:
785 exact bytes, three deterministic put/put/meta operations, SHA-256
`45174ce9b17341c5bd9efb800a403308d148ce799b60862896341c248a37e623`,
default-profile hash `2269240639182000324`. Its Unicode includes composed and
decomposed accents, CJK and emoji; newlines inside payloads remain escaped.
Published Pollard wheel SHA-256:
`569fb5f130a82c9be327b8dcbd285e3be063200bd9773ca15c5d6bb62edd627f`.

For external applications, preserve exact encoded bytes, modulus, base, and
numeric representation together. Hashrope does not silently migrate old stored
hashes, choose a new default, normalize JSON, or turn fingerprints into
cryptographic identity proofs.

The reduction bug fix deliberately corrects previously wrong results for small
Mersenne moduli and unrestricted full-width reduction inputs: for example,
`mersenne_mod(256, 3)` changes from 13 to the Python/mathematical result 1.
Applications persisting results from those affected custom configurations
must identify their generating version before comparing or rebuilding hashes.
No stored values are rewritten. The default 61-bit profile, including Pollard's
persisted log hash, is unchanged.

Bounded streaming refers to retained arena storage after eviction. A large
append temporarily allocates its input graph, and compaction temporarily holds
old and new graphs. Arbitrary-precision counters/hashes also scale with their
bit widths. No constant-total-memory or per-append O(log w) claim is made for
the optional graph-copy reclamation step.

## Validation and limits

Validation commands and final outcomes are recorded in `VALIDATION.md`.
Cross-platform release gates and the publication process are separate from
the recorded baseline audit; no downstream source modification is required.
Adoption by hashrope-bio requires a separate authorized dependency correction.
Rust allocator/platform capacity and u32 arena node counts remain finite even
when logical lengths and hashes use arbitrary-precision integers.
