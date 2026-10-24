"""Profile exact Python comparisons with construction outside measurement.

This diagnostic is separate from the paired benchmark reruns. It records
cProfile call counts and costs, source identity, and unprofiled query medians.
"""
from __future__ import annotations

import argparse
import cProfile
import hashlib
import io
import json
from pathlib import Path
import platform
import pstats
import statistics
import sys
import time


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=root / "packages/python")
    parser.add_argument("--label", required=True)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--corpus", type=Path,
                        default=root.parent / "hashropeictaipaper/data/raw/corpus_s42.txt")
    args = parser.parse_args()
    sys.path.insert(0, str(args.package_root.resolve()))
    import hashrope.rope as module
    from hashrope import (Internal, Leaf, PolynomialHash, rope_concat, rope_lcp,
                          rope_range_lcp, rope_repeat, rope_split)

    corpus = args.corpus.read_bytes()
    payload = corpus[:2 * 1024 * 1024]
    assert len(payload) == 2 * 1024 * 1024
    h = PolynomialHash()

    def build(data, width=4096):
        nodes = [Leaf(data[i:i + width], h) for i in range(0, len(data), width)]
        while len(nodes) > 1:
            nodes = [Internal(nodes[i], nodes[i + 1], h) if i + 1 < len(nodes)
                     else nodes[i] for i in range(0, len(nodes), 2)]
        return nodes[0]

    a = build(payload)
    b = build(bytes(bytearray(payload)))
    changed = payload[:-1] + bytes([payload[-1] ^ 1])
    b_changed = build(changed)
    prefix, _ = rope_split(a, len(payload) - 1, h)
    shared_changed = rope_concat(prefix, Leaf(changed[-1:], h), h)
    shifted_a, shifted_b = build(b"abc" + payload), build(b"12345" + payload)
    leaf_a = Leaf(payload[:65536], h)
    leaf_late = Leaf(payload[:65535] + bytes([payload[65535] ^ 1]), h)
    leaf_early = Leaf(bytes([payload[0] ^ 1]) + payload[1:65536], h)
    repeat_a = rope_repeat(Leaf(b"abc", h), 10**12, h)
    repeat_b = rope_repeat(Leaf(b"abcabc", h), 5 * 10**11, h)
    dag_a, dag_b = Leaf(b"abc", h), Leaf(bytes(bytearray(b"abc")), h)
    for _ in range(40):
        dag_a, dag_b = Internal(dag_a, dag_a, h), Internal(dag_b, dag_b, h)

    cases = [
        ("independent_2MiB_equal_4KiB_leaves", lambda: rope_lcp(a, b), len(payload)),
        ("independent_2MiB_late_mismatch", lambda: rope_lcp(a, b_changed), len(payload) - 1),
        ("shared_2MiB_late_edit", lambda: rope_lcp(a, shared_changed), len(payload) - 1),
        ("independent_shifted_equal_ranges", lambda: rope_range_lcp(
            shifted_a, 3, shifted_b, 5, len(payload)), len(payload)),
        ("64KiB_leaf_last_byte_mismatch", lambda: rope_lcp(leaf_a, leaf_late), 65535),
        ("64KiB_leaf_first_byte_mismatch", lambda: rope_lcp(leaf_a, leaf_early), 0),
        ("independent_trillion_repeats", lambda: rope_lcp(repeat_a, repeat_b), 3 * 10**12),
        ("independent_self_concat_DAG", lambda: rope_lcp(dag_a, dag_b), 3 * 2**40),
    ]
    report = {
        "label": args.label, "python": platform.python_version(),
        "source": str(Path(module.__file__).resolve()),
        "source_sha256": hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
        "corpus": str(args.corpus.resolve()), "corpus_sha256": hashlib.sha256(corpus).hexdigest(),
        "repetitions": args.reps, "construction_timed": False,
        "note": "Diagnostic single-process timings; use paired rerun results for performance claims.",
        "cases": [],
    }
    text = io.StringIO()
    for name, query, expected in cases:
        assert query() == expected, name
        samples = []
        for _ in range(9):
            start = time.perf_counter_ns()
            for _ in range(args.reps):
                actual = query()
            samples.append((time.perf_counter_ns() - start) / args.reps / 1e6)
            assert actual == expected, name
        profile = cProfile.Profile()
        profile.enable()
        for _ in range(args.reps):
            query()
        profile.disable()
        stats = pstats.Stats(profile, stream=text).strip_dirs().sort_stats("tottime")
        text.write(f"\n{name}: {args.reps} queries, expected={expected}\n")
        stats.print_stats(16)
        report["cases"].append({"name": name, "expected": expected,
                                "median_ms": statistics.median(samples),
                                "samples_ms": samples,
                                "profile_total_calls": stats.total_calls,
                                "profile_total_seconds": stats.total_tt})
        print(f"{name}: {statistics.median(samples):.6f} ms", flush=True)
    prefix = Path(__file__).with_name(f"python_profile.{args.label}")
    prefix.with_suffix(prefix.suffix + ".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    prefix.with_suffix(prefix.suffix + ".txt").write_text(text.getvalue(), encoding="utf-8")


if __name__ == "__main__":
    main()
