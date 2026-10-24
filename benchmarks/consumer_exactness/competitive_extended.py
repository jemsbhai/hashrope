"""Rerun EXP017 real-pair and K-candidate cells with exact helpers.

Calls the original paper benchmark functions directly. Preparation is excluded,
as before; fixed timed arm order is radix, Hashrope exact, then NumPy. Nine fresh
processes retain the original 200 pairs/dataset/seed and full bounded K grid.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("competitive_extended_results"))
    parser.add_argument("--tokens", type=Path, default=Path(__file__).with_name("paper_results"))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--invocation", type=int)
    args = parser.parse_args()
    root = Path(os.environ.get("HASHROPE_BENCH_CORE", Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(args.paper.resolve()))
    sys.path.insert(0, str(root / "packages/python"))
    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.worker:
        for seed in [42, 43, 44]:
            for invocation in range(3):
                print(f"Extended competitive worker seed={seed} invocation={invocation}", flush=True)
                subprocess.run([sys.executable, __file__, "--paper", str(args.paper),
                    "--tokens", str(args.tokens), "--output", str(args.output), "--worker",
                    "--seed", str(seed), "--invocation", str(invocation)], check=True)
        return

    import numpy as np
    import torch
    import hashrope
    assert Path(hashrope.__file__).resolve() == (root / "packages/python/hashrope/__init__.py").resolve()
    assert Path(sys.modules["hashrope.rope"].__file__).resolve() == (root / "packages/python/hashrope/rope.py").resolve()
    from scripts.exp017_bench import (DEFAULT_K_SWEEP, run_k_sweep_cell,
                                     run_real_pair_cells, verify_radix_sha)
    verify_radix_sha(args.paper)
    sources = [Path(__file__), root / "packages/python/hashrope/rope.py",
               args.paper / "scripts/exp017_bench.py", args.paper / "src/lcp.py",
               args.paper / "src/competitive.py", args.paper / "src/realpairs.py",
               args.paper / "src/flatten.py",
               args.paper / "third_party/sglang_radix_cache/radix_cache.py"]
    datasets = [args.paper / "data/canonical" / name
                for name in ["sharegpt_sample.jsonl", "lmsys_sample.jsonl"]]
    tokens_file = args.tokens / f"tokens_s{args.seed}.npy"
    result = {"comparison_mode": "exact_structural", "seed": args.seed,
              "invocation": args.invocation, "reps": 5, "real_pairs_per_dataset": 200,
              "method": "original EXP017 real-pair and K cells; construction excluded; fixed radix/exact/NumPy timing order",
              "source_sha256": {str(p.resolve()): sha(p) for p in sources},
              "data_sha256": {str(p.resolve()): sha(p) for p in [*datasets, tokens_file]},
              "python": platform.python_version(), "platform": platform.platform(), "hashrope_file": hashrope.__file__,
              "numpy": np.__version__, "torch": torch.__version__, "torch_threads": torch.get_num_threads(),
              "real_pairs": {}, "k_sweep": {}}
    output = args.output / f"run_s{args.seed}_i{args.invocation}.json"
    for path in datasets:
        print(f"  Real pairs: {path.name}", flush=True)
        cell = run_real_pair_cells(str(path), 200, args.seed, 5, args.paper)
        assert cell["n_pairs"] == 200 and cell["all_correct"] and cell["mismatches"] == 0
        result["real_pairs"][path.name] = cell
    tokens = np.load(tokens_file).tolist()
    for length, ks in DEFAULT_K_SWEEP.items():
        for k in ks:
            print(f"  K sweep: L={length} K={k}", flush=True)
            cell = run_k_sweep_cell(tokens, length, k, args.seed, 5, args.paper)
            assert cell["correctness_ok"]
            result["k_sweep"][f"{length}/{k}"] = cell
    assert result["source_sha256"] == {str(p.resolve()): sha(p) for p in sources}
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
