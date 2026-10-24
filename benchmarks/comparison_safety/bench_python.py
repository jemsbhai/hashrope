"""Reproducible baseline/current benchmark; no third-party dependencies.

Run from any directory with --package pointing at packages/python. Each process
imports exactly that source tree. Timings include return-value disposal, exclude
fixture preparation, and use warmed hash-power caches. JSON is written to stdout.
"""
from __future__ import annotations

import argparse
import gc
import json
import platform
import statistics
import sys
import time
import tracemalloc
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--exact", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.package.resolve()))
    import hashrope.rope as r
    from hashrope.polynomial_hash import PolynomialHash

    measurements = []
    scratch_memory = []

    def scratch(name, function):
        gc.collect()
        tracemalloc.start()
        function()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        scratch_memory.append({"name": name, "retained_bytes": current,
                               "peak_bytes": peak})

    def measure(name, function, iterations):
        for _ in range(min(iterations, 3)):
            function()
        timings = []
        for _ in range(args.samples):
            gc.collect()
            gc.disable()
            try:
                start = time.perf_counter_ns()
                for _ in range(iterations):
                    function()
                timings.append((time.perf_counter_ns() - start) / iterations)
            finally:
                gc.enable()
        measurements.append({"name": name, "iterations_per_sample": iterations,
                             "ns_per_op": timings,
                             "median_ns": statistics.median(timings)})

    def chunked(data, h):
        nodes = [r.rope_from_bytes(data[i:i + 512], h)
                 for i in range(0, len(data), 512)]
        while len(nodes) > 1:
            nodes = [r.rope_concat(nodes[i], nodes[i + 1], h)
                     if i + 1 < len(nodes) else nodes[i]
                     for i in range(0, len(nodes), 2)]
        return nodes[0] if nodes else None

    def hash_lcp(a, b, h):
        lo, hi = 0, min(r.rope_len(a), r.rope_len(b))
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if r.rope_substr_hash(a, 0, mid, h) == r.rope_substr_hash(b, 0, mid, h):
                lo = mid
            else:
                hi = mid - 1
        return lo

    data = bytes(range(256)) * 256
    h = PolynomialHash()
    left = chunked(data, h)
    right = chunked(data, h)
    measure("core/from_bytes_64KiB", lambda: r.rope_from_bytes(data, h), 2)
    measure("core/chunked_build_64KiB", lambda: chunked(data, h), 2)
    measure("core/concat_balanced_64KiB", lambda: r.rope_concat(left, right, h), 2000)
    measure("core/split_chunked_64KiB", lambda: r.rope_split(left, 32769, h), 20)
    measure("core/substr_warm_64KiB", lambda: r.rope_substr_hash(left, 1, len(data)-2, h), 20)
    measure("core/root_hash", lambda: r.rope_hash(left), 100000)
    repeated = r.rope_repeat(r.rope_from_bytes(b"ab", h), 2**32, h)
    suffix = r.rope_from_bytes(b"!", h)
    measure("core/repeat_join_2pow32", lambda: r.rope_concat(repeated, suffix, h), 4)

    memory = []
    for name, factory in [("single_leaf_64KiB", lambda ph: r.rope_from_bytes(data, ph)),
                          ("chunked_64KiB", lambda ph: chunked(data, ph))]:
        runs = []
        for _ in range(3):
            gc.collect()
            tracemalloc.start()
            ph = PolynomialHash()
            node = factory(ph)
            current, peak = tracemalloc.get_traced_memory()
            runs.append({"retained_bytes": current, "peak_bytes": peak,
                         "root_object_bytes": sys.getsizeof(node)})
            tracemalloc.stop()
            del node, ph
        memory.append({"name": name, "runs": runs})

    if args.exact:
        # Small and medium buffers separate early-exit costs from full scans.
        for size in (65536, 1048576):
            raw = bytes(range(256)) * (size // 256)
            first = chunked(raw, h)
            independent = chunked(raw, h)
            early = chunked(bytes([raw[0] ^ 1]) + raw[1:], h)
            late = chunked(raw[:-1] + bytes([raw[-1] ^ 1]), h)
            prefix, _ = r.rope_split(first, size - 1, h)
            shared_late = r.rope_concat(prefix, r.rope_from_bytes(bytes([raw[-1] ^ 1]), h), h)
            scenarios = [("same_root", first, first, size),
                         ("independent_equal", first, independent, size),
                         ("early_mismatch", first, early, 0),
                         ("late_mismatch", first, late, size - 1),
                         ("shared_late_mismatch", first, shared_late, size - 1)]
            for label, a, b, expected in scenarios:
                assert r.rope_lcp(a, b) == expected, label
                assert r.rope_equal(a, b) == (expected == size), label
                assert hash_lcp(a, b, h) == expected, label
                prefix_name = f"compare/{size}/{label}"
                measure(prefix_name + "/exact_lcp", lambda a=a, b=b: r.rope_lcp(a, b), 20)
                measure(prefix_name + "/exact_equal", lambda a=a, b=b: r.rope_equal(a, b), 20)
                measure(prefix_name + "/hash_only_lcp", lambda a=a, b=b: hash_lcp(a, b, h), 3)
                if size == 1048576 and label in ("independent_equal", "shared_late_mismatch"):
                    scratch(prefix_name + "/exact_lcp", lambda a=a, b=b: r.rope_lcp(a, b))
                    scratch(prefix_name + "/hash_only_lcp", lambda a=a, b=b: hash_lcp(a, b, h))
        # 2 TiB represented in a handful of nodes, with different pattern sizes.
        a = r.rope_repeat(r.rope_from_bytes(b"ab", h), 2**40, h)
        b = r.rope_repeat(r.rope_from_bytes(b"abab", h), 2**39, h)
        assert r.rope_lcp(a, b) == 2**41
        assert r.rope_equal(a, b)
        measure("compare/repeat_2TiB/exact_lcp", lambda: r.rope_lcp(a, b), 200)
        measure("compare/repeat_2TiB/exact_equal", lambda: r.rope_equal(a, b), 200)
        measure("compare/repeat_2TiB/hash_only_lcp", lambda: hash_lcp(a, b, h), 2)
        scratch("compare/repeat_2TiB/exact_lcp", lambda: r.rope_lcp(a, b))
        a = r.rope_from_bytes(b"ab", h)
        b = r.rope_from_bytes(b"ab", h)
        for _ in range(40):
            a = r.rope_concat(a, a, h)
            b = r.rope_concat(b, b, h)
        assert r.rope_lcp(a, b) == 2**41
        assert r.rope_equal(a, b)
        measure("compare/self_concat_DAG_2TiB/exact_lcp", lambda: r.rope_lcp(a, b), 100)
        measure("compare/self_concat_DAG_2TiB/exact_equal", lambda: r.rope_equal(a, b), 100)
        measure("compare/self_concat_DAG_2TiB/hash_only_lcp", lambda: hash_lcp(a, b, h), 2)
        scratch("compare/self_concat_DAG_2TiB/exact_lcp", lambda: r.rope_lcp(a, b))

    print(json.dumps({"runtime": platform.python_version(), "platform": platform.platform(),
                      "package": str(args.package.resolve()), "samples": args.samples,
                      "measurements": measurements, "memory": memory,
                      "scratch_memory": scratch_memory}, indent=2))


if __name__ == "__main__":
    main()
