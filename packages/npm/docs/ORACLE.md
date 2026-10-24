# Independent Python oracle

The differential corpus was generated with Python 3.12.2 from the published
`hashrope-0.2.2-py3-none-any.whl`, SHA-256
`710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2`.
The generator verifies the digest before importing, adds the wheel itself to
`sys.path`, and asserts that every imported Hashrope module originated inside
that wheel. It does not import the sibling Python checkout, Java source, or JS
implementation. Run the script with Python's isolated-interpreter flag:

```sh
python -I scripts/test-oracle-generator.py
python -I scripts/generate-oracle.py
npm run build
node scripts/audit-oracle-api.mjs
node --test test/oracle.test.mjs
```

The wheel is a local development input, excluded from the npm tarball. To
reproduce in a fresh checkout, download that exact wheel from PyPI into
`oracle/` first; the script rejects a missing or mismatched artifact. The pinned wheel digest and original import/source provenance are retained in `test/fixtures/manifest.json`; current npm source provenance is in `evidence/upstream-0.2.2-provenance.json`.

`test/fixtures/manifest.json` records the original 21-name `__all__` inventory,
public class members, balance constants, source-file digests, seed 20261023,
wheel/import provenance, corpus digest, per-kind counts, and upstream license
digest. All integers in numerical inputs and results are decimal strings; byte
strings are hexadecimal. Tests decode integers directly to BigInt.

The corpus has **27,564 rows**, including **42 expected exception rows**:

| Kind | Records |
| --- | ---: |
| byte hashes | 155 |
| modular powers | 60 |
| geometric accumulator | 595 |
| modular folding | 110 |
| modular multiplication | 440 |
| concatenated hashes | 180 |
| repeated hashes | 180 |
| overlap hashes | 180 |
| rope metadata/bytes/tree shape | 176 |
| split pairs | 4,672 |
| substring hashes | 20,259 |
| complete ordinary stream snapshots | 520 |
| stream edge/error snapshots | 16 |
| invalid hash constructors | 5 |
| invalid window constructors | 5 |
| invalid node constructors | 7 |
| unbalanced-tree validation | 4 |

Every row is registered as an individual Node test. The manifest test verifies
corpus bytes, count, unique IDs, kind inventory and expected-error inventory.
The adapter fails on an unknown record kind. Tree state checks include exact
byte length, expanded leaf weight, stored fingerprint, structural height,
materialized bytes when small, structural validity, and a SHA-256 of the exact
canonical tree shape. That SHA-256 is a test-only structural digest; it is not
part of the library API or an alternative content-equality guarantee.

Profiles include both documented moduli with base 131, the 127-bit modulus with
base 257, modulus 31/base 3, and arithmetic-only positive non-Mersenne modulus
100/base 3. Cases include signed arithmetic, negative modular exponents, the
base-131 collision, UTF-8 split by bytes, append/prepend balancing histories,
nested repetition, lengths above 2^100 without materialization, malformed
substring behavior, eviction thresholds, overlapping copies and negative-copy
quirks. Invalid positive-domain arguments and Python error messages are tested;
VM-dependent null/division diagnostics use documented JS error-category mapping.

Nonpositive-modulus Python folding may not terminate and is deliberately not
invoked. Resource exhaustion, Python duck types, object repr/hash, arbitrary
reflection, and every possible subclass are outside the fixture contract.
Dedicated JS unit/property tests cover ownership, safe Number rejection,
materialization limits, exact equality, dispatch and persistence beyond this
published-wheel corpus. Corpus size alone is not a claim of proof or complete
all-language identity.

The following historical files remain in the original standalone 0.2.2 checkout: `evidence/oracle-generator-red.log`, `oracle-generator-green.log`, and
`oracle-generator-refactor.log` record the provenance TDD cycle. The initial
implementation syntax failure is retained separately. `oracle-js-red.log`
records the differential acceptance suite before the public JS entry point was
built. Final pass/fail outcomes belong to the verification report, rather than
the fixed corpus counts above.
