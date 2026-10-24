# Python and biology release preparation

Version **0.3.4** is published and verified from public PyPI for both Python
packages. The preparation source commit for biology is
`8e07e641af532849a629bbb595c18afcfd7ab6ab`, pushed on `codex/release-0.3.4`.

Both Python wheels were built from their source distributions. Strict Twine
checks pass for all four distributions. A fresh environment imports both
packages from site-packages, passes `pip check`, all **244 installed tests**,
and the collision/range/profile/symbolic-repeat release gate. Biology's Rust
library, consumer and benchmark-helper regressions and doctest total **15
passing tests** against local core 0.3.4. No benchmark workload was run.

Only versions, dependency requirements, release documentation and benchmark
version provenance changed in these packages. The build driver verifies Python
AST equivalence after allowing the `__version__` metadata update. Biology Rust
consumer algorithms are unchanged; its benchmark report names the new required
core version. API/version alignment does not establish equal performance across
languages or remove the unresolved core comparison costs documented in the
[release notes](../README.md).

- [Artifact hashes, source hashes and commands](python-artifacts.json)
- [Installed Python tests](python-final-installed-tests.log)
- [Exact release gate](python-final-gate.log)
- [Strict distribution checks](python-twine-check.log)
- [Biology Rust validation](bio-rust-tests.json) and [log](bio-rust-tests.log)
- [Final preparation identity](python-bio-ready.json)
- [Reproducible build driver](build-python-bio.py)
- [Post-publication verification driver](verify-python-registry.py)

Artifacts are retained under the ignored `.release/0.3.4/python-dist` directory;
their exact paths, sizes and SHA-256 values are recorded above. Standard PyPI
credentials are available through the existing user configuration; only
availability booleans are recorded, never secret values. The reviewed artifacts
were uploaded after all coordinated release gates and source merges passed.

All four canonical public downloads match the reviewed hashes. A new empty
environment installed `hashrope==0.3.4` and `hashrope-bio==0.3.4` through the
ordinary PyPI index with isolated pip configuration and no cache. `pip check`,
installed origin/version assertions, the exact gate and **244 tests** passed
(34.99 seconds for the test suite). See the [public result](python-public-install.json),
[log](python-public-install.log) and [publication record](python-publication.json).

Initial canonical metadata requests returned 404 during propagation, followed
by two ordinary pip resolver attempts that could not yet see one or both new
versions. All failed logs and environments remain preserved. Cache-bypass
requests were diagnostic only; they did not supply the final installation.
The verifier's logic stayed byte-identical, with wrappers choosing new empty
environment paths. The [attempt lineage](python-public-lineage.json) binds the
failures, probes, wrapper sources and final successful result by SHA-256.
