"""Run three alternating baseline/current rounds and save all raw measurements.

Example (from repository root):
  python benchmarks/comparison_safety/run_benchmarks.py --baseline-ref 100749d

The baseline is extracted from git into a fresh temporary directory. No checkout,
tracked source, or submodule is modified. The Rust harness builds offline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

from bench_rust import parse_output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-ref", required=True)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    scripts = Path(__file__).resolve().parent
    root = scripts.parent.parent
    baseline = Path(tempfile.mkdtemp(prefix="hashrope-benchmark-baseline-"))
    archive = baseline / "source.tar"
    subprocess.run(["git", "archive", "--format=tar", "--output", str(archive), args.baseline_ref,
                    "packages/python", "packages/rust"], cwd=root, check=True)
    with tarfile.open(archive) as source:
        source.extractall(baseline, filter="data")
    sources = {"baseline": baseline, "current": root}
    results = {"baseline_ref": subprocess.check_output(["git", "rev-parse", args.baseline_ref], cwd=root, text=True).strip(),
               "platform": platform.platform(),
               "python": platform.python_version(),
               "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
               "rounds": args.rounds, "order": [], "python_runs": [], "rust_runs": []}
    if sys.platform == "win32":
        results["cpu"] = json.loads(subprocess.check_output([
            "powershell", "-NoProfile", "-Command",
            "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json -Compress"
        ], text=True))
    paths = ["packages/python/hashrope/rope.py", "packages/python/hashrope/polynomial_hash.py",
             "packages/python/hashrope/sliding.py", "packages/rust/src/rope.rs",
             "packages/rust/src/polynomial_hash.rs"]
    results["current_source_sha256"] = {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in paths}
    binaries = {}
    for variant, directory in sources.items():
        command = [sys.executable, str(scripts / "bench_rust.py"), "--package",
                   str(directory / "packages/rust"), "--build-only"]
        if variant == "current":
            command.append("--exact")
        binaries[variant] = subprocess.check_output(command, text=True).strip()
    for round_id in range(args.rounds):
        order = ["baseline", "current"] if round_id % 2 == 0 else ["current", "baseline"]
        results["order"].append(order)
        for variant in order:
            print(f"Round {round_id + 1}: Python {variant}", file=sys.stderr)
            command = [sys.executable, str(scripts / "bench_python.py"), "--package",
                       str(sources[variant] / "packages/python")]
            if variant == "current":
                command.append("--exact")
            run = json.loads(subprocess.check_output(command, text=True, timeout=180))
            run.update(round=round_id + 1, variant=variant)
            results["python_runs"].append(run)
        for variant in order:
            print(f"Round {round_id + 1}: Rust {variant}", file=sys.stderr)
            run = parse_output(subprocess.check_output([binaries[variant]], text=True, timeout=180))
            run.update(round=round_id + 1, variant=variant)
            results["rust_runs"].append(run)
        args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    results["source_unchanged_during_run"] = all(
        hashlib.sha256((root / p).read_bytes()).hexdigest() == digest
        for p, digest in results["current_source_sha256"].items())
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
