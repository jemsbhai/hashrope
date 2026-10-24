# Python and biology release preparation

Version **0.3.4** is built and verified locally. These records do not claim
publication. The biology source commit is
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
availability booleans are recorded, never secret values. Publication must use
the reviewed artifacts after the coordinated release gates pass. The public
verification driver then checks registry byte hashes, creates a fresh PyPI-only
installation, and reruns the exact gate and combined installed test suite.
