"""Build a separate Rust/Ropey exact-comparison workload and run 3x3 processes.

Does not change EXP020's scratch project, Criterion output, or historical JSON.
Hashrope and Ropey both perform each seeded edit before querying a concrete
prebuilt target. Targets are 50% equal, 50% different in the final byte.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--core-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = args.core_root or here.parents[1]
    sources = [here / "edit_compare.rs", root / "packages/rust/Cargo.toml",
               *sorted((root / "packages/rust/src").glob("*.rs"))]
    source_hashes = {str(p): sha256(p) for p in sources}
    if args.binary:
        binary = args.binary
        build = json.loads(binary.with_suffix(".build.json").read_text(encoding="utf-8"))
        assert build["source_sha256"] == source_hashes, "Sources changed after compilation; rebuild"
        assert build["binary_sha256"] == sha256(binary), "Compiled binary changed"
    else:
        project = Path(tempfile.mkdtemp(prefix="hashrope-exact-edit-"))
        (project / "src").mkdir()
        shutil.copyfile(here / "edit_compare.rs", project / "src/main.rs")
        dependency = json.dumps((root / "packages/rust").as_posix())
        (project / "Cargo.toml").write_text(
            '[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n'
            '[dependencies]\nhashrope={path=' + dependency + '}\nropey="=1.6.1"\n', encoding="utf-8")
        shutil.copyfile(here / "edit_compare.Cargo.lock", project / "Cargo.lock")
        subprocess.run(["cargo", "build", "--release", "--locked", "--manifest-path", str(project / "Cargo.toml")], check=True)
        binary = project / "target/release/exact-edit-comparison.exe"
        assert source_hashes == {str(p): sha256(p) for p in sources}, "Sources changed during compilation"
        build = {"source_sha256": source_hashes, "binary_sha256": sha256(binary),
                 "cargo_manifest": (project / "Cargo.toml").read_text(encoding="utf-8"),
                 "cargo_lock": (project / "Cargo.lock").read_text(encoding="utf-8"),
                 "cargo_lock_sha256": sha256(project / "Cargo.lock"),
                 "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip()}
        binary.with_suffix(".build.json").write_text(json.dumps(build, indent=2), encoding="utf-8")
    if args.build_only:
        print(binary)
        return
    results = {"method": "new edit-plus-target-comparison workload, not historical EXP020; 50 cycles, 50% equal / 50% late mismatch",
               "seeds": [1, 2, 3], "invocations": 3, "samples_per_process": 3,
               "arm_order": "rotated by (seed+invocation) modulo 4",
               "sizes": [10000,100000,1000000], "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
               "build": build, "source_sha256": source_hashes, "runs": []}
    for seed in results["seeds"]:
        for invocation in range(3):
            print(f"Edit-plus-exact comparison seed={seed} invocation={invocation}", flush=True)
            output = subprocess.check_output([str(binary), str(seed), str(invocation)], text=True)
            cells = {}
            for line in output.splitlines():
                _, n, shared, arm, _, ms = line.split("\t")
                key = f"{n}/{shared}/{arm}"
                cell = cells.setdefault(key, {"bytes": int(n), "shared_history": shared == "true",
                                              "arm": arm, "samples_ms_per_cycle": []})
                cell["samples_ms_per_cycle"].append(float(ms))
            for cell in cells.values():
                cell["median_ms_per_cycle"] = statistics.median(cell["samples_ms_per_cycle"])
            results["runs"].append({"seed": seed, "invocation": invocation, "cells": list(cells.values())})
            destination = args.output or here / "edit_results.json"
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
