# Historical development evidence and current checks

The chronology below describes the original standalone npm 0.2.2 implementation. Historical logs/snapshots remain in the preserved authoritative standalone checkout; they were not regenerated or copied wholesale into this monorepo. It is not a claim that the 0.3.3 correction followed those historical red/green stages. Current 0.3.4 checks are recorded in [VERIFICATION.md](VERIFICATION.md), and source provenance is in [RELEASE.md](RELEASE.md).

Development proceeded through test-first stages in parallel modules. Red-green-refactor was applied at the feature level; the overall work was not one globally serialized cycle. For example, the query red run occurred before the first complete rope green while node dependencies were still arriving. The preserved timestamps make this overlap explicit; a claim that every stage finished green before any other stage began would be inaccurate. Test sources were written before feature implementations; recorded runs show missing APIs or explicit throwing scaffolds, followed by implementations and green/refactor runs. The initial package/arithmetic/node/oracle baselines include missing-module failures; later behavioral reds fail on the actual unimplemented method or assertion. These distinctions are retained rather than presenting every red as a deep semantic failure.

The shared `scripts/tdd.mjs` runner records UTC time, exact command, exit code, output, and SHA-256 source/test snapshots under `evidence/tdd/<label>/`. Initial baselines before the runner existed have module-specific log files. Evidence files and upstream downloads are source-checkout artifacts, deliberately excluded from the npm tarball. No historical commit chronology or branch history is fabricated.

| Area | Recorded stages | Evidence prefixes |
|---|---|---|
| Package | metadata/export/allowlist contract red, minimal package green; documentation allowlist extension red/green | `tooling-*` |
| Arithmetic | integer/folding primitives; profile/hash/power; phi/concat/repeat/overlap; refactor | `polynomial-01-*`, `polynomial-02-*`, `polynomial-03-*`, `polynomial-refactor*` |
| Nodes | owned leaves; internal/repeat metadata; materialization budget; validator; structural equality | `nodes-*` |
| Rope | throwing concat/split scaffold, repeated-node paths; implementation green | `rope-01-*` |
| Queries | repetition and substring red; numerical/fast-path/malformed range green | `rope-02-queries-*` |
| JS additions | persistent edits, exact equality and seeded model red; implementation green | `rope-03-edits-equality-*` |
| Streaming | constructor/accessors; literals/eviction; back-reference copies and signed quirks | `sliding-01-*`, `sliding-02-*`, `sliding-03-*` |
| Oracle | failing missing-generator provenance tests, verified-wheel generator green and refactor; public package differential run | `oracle-*` |
| Consumers | absent-artifact red; packed package and clean installation verification | `tooling-*`, `consumer-verification-*` |
| Independent review | overridden byte-iterator discrepancy red; intrinsic byte iteration fix green plus ownership/equality/huge-range checks | `review-01-byte-iterator-*` |

The first full oracle run passed all nonstreaming records while 525 streaming records failed against a build without `appendCopy`. That failed run is retained as `evidence/oracle-js-full-attempt1.log`. An initial rope green attempt also failed on node validation/materialization stubs being developed concurrently. Neither was called a passing result. See [verification](VERIFICATION.md) for the final integrated outcome.

To repeat tests from the final source checkout:

```sh
npm ci
npm test
python -I scripts/test-oracle-generator.py
python -I scripts/generate-oracle.py
npm run pack:verified
npm run test:consumer
```

The committed fixtures are independent expected values from the pinned published Python wheel; JS output was never used to manufacture expected oracle values. Deterministic byte/modular models additionally check edits and streaming. Full TDD snapshots can be inspected with their command/output/manifest; they are provenance evidence, not a claim that every intermediate dependency combination was a releasable package.
