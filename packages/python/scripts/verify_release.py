"""Release gate for the hashrope Python package.

Run this against the FRESHLY BUILT WHEEL installed in a scratch venv
BEFORE `twine upload`, and again AFTER publishing against the real
PyPI install. Exits non-zero on any failure.

    & "$env:TEMP\hr-verify\Scripts\python.exe" scripts\verify_release.py

The first output line shows exactly which installation is being
exercised -- confirm it points at the venv's site-packages, not the
local source tree.
"""
import os
import random
import sys
import time

import hashrope
from hashrope import PolynomialHash, rope_concat, rope_from_bytes, rope_substr_hash

CHUNK = 4096
N = 500_000
HALF_RANGE_BUDGET_MS = 50.0  # 0.2.2 fast path: ~1 ms typical; pre-fix: ~120 ms


def build(data: bytes, h: PolynomialHash):
    leaves = [rope_from_bytes(data[i:i + CHUNK], h)
              for i in range(0, len(data), CHUNK)]
    while len(leaves) > 1:
        leaves = [
            rope_concat(leaves[i], leaves[i + 1], h)
            if i + 1 < len(leaves) else leaves[i]
            for i in range(0, len(leaves), 2)
        ]
    return leaves[0]


def main() -> int:
    print(f"hashrope {hashrope.__version__} "
          f"loaded from: {os.path.dirname(hashrope.__file__)}")

    rng = random.Random(20260610)
    h = PolynomialHash()
    data = os.urandom(N)
    rope = build(data, h)

    # Gate 1 -- correctness: substring hashes match the byte-level oracle.
    for _ in range(50):
        length = rng.randint(1, 4000)
        start = rng.randint(0, N - length)
        if rope_substr_hash(rope, start, length, h) != h.hash(data[start:start + length]):
            print("FAIL: substr hash mismatch vs byte oracle", file=sys.stderr)
            return 1
    print("gate 1 OK: 50/50 random ranges match byte oracle")

    # Gate 2 -- Theorem 9: half-range query must be metadata-bound, not O(n).
    times = []
    for _ in range(100):
        t0 = time.perf_counter()
        rope_substr_hash(rope, 0, N // 2, h)
        times.append((time.perf_counter() - t0) * 1e3)
    times.sort()
    p50 = times[50]
    print(f"gate 2: half-range substr_hash p50 = {p50:.3f} ms "
          f"(budget {HALF_RANGE_BUDGET_MS:.0f} ms)")
    if p50 > HALF_RANGE_BUDGET_MS:
        print("FAIL: Theorem 9 fast path not active in this build", file=sys.stderr)
        return 1

    print("OK: release gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
