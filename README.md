# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

Zero dependencies. Formally proven. Multiple language implementations.

## Packages

| Language | Path | Status |
|----------|------|--------|
| Python | [`packages/python`](packages/python) | 88 tests passing |
| Rust | `packages/rust` | Planned |
| C | `packages/c` | Planned |

## What it does

hashrope is a rope data structure where every node carries a polynomial hash, enabling O(log n) concat, split, substring hashing, and O(log q) repetition encoding -- all with strict BB[2/7] balance guarantees and Mersenne prime modular arithmetic.

See [`packages/python/README.md`](packages/python/README.md) for full API documentation.

## License

MIT
