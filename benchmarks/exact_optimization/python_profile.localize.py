"""Isolated before/after leaf-offset sweep; diagnostic, not a final benchmark."""
import ast
import hashlib
import json
from pathlib import Path
import statistics
import time

root = Path(__file__).resolve().parents[2]
base = root / "benchmarks/exact_optimization/.cache/baseline/hashrope/packages/python/hashrope/rope.py"
current = root / "packages/python/hashrope/rope.py"


def load(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_leaf_lcp")
    namespace = {}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    return namespace["_leaf_lcp"]


old, new = load(base), load(current)
rows = []
for length in (4096, 65536):
    left = b"x" * length
    for position in (0, 1, 2, 8, 15, 16, 63, 64, length // 2, length - 1):
        right = left[:position] + b"y" + left[position + 1:]
        row = {"length": length, "position": position}
        for label, function in (("baseline", old), ("optimized", new)):
            assert function(left, 0, right, 0, length) == position
            values = []
            for _ in range(7):
                start = time.perf_counter_ns()
                for _ in range(100):
                    function(left, 0, right, 0, length)
                values.append((time.perf_counter_ns() - start) / 100 / 1000)
            row[label + "_us"] = statistics.median(values)
            row[label + "_samples_us"] = values
        rows.append(row)
        print(length, position, row["baseline_us"], row["optimized_us"])
report = {"baseline_sha256": hashlib.sha256(base.read_bytes()).hexdigest(),
          "optimized_sha256": hashlib.sha256(current.read_bytes()).hexdigest(),
          "repetitions": 100, "samples": 7,
          "note": "Isolated diagnostic timings, not final consumer benchmark evidence.",
          "cases": rows}
(root / "benchmarks/exact_optimization/python_profile.localization.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8")
