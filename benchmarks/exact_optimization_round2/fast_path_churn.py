"""Run the unchanged edit-plus-comparison harness against fast-path ablations.

Preliminary diagnostics only: seed 1, one fresh process per variant, three
samples per cell. The final benchmark matrix supplies independent repetitions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess

HERE = Path(__file__).resolve().parent
CONSUMERS = HERE.parent / "consumer_exactness"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(variant):
    project = HERE / ".cache" / ("fast_churn_" + variant)
    library = HERE / ".cache" / ("fast_path_" + variant) / "hashrope"
    (project / "src").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(CONSUMERS / "edit_compare.rs", project / "src/main.rs")
    manifest = ('[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n'
                '[dependencies]\nhashrope={path=' + json.dumps(library.as_posix()) + '}\nropey="=1.6.1"\n')
    (project / "Cargo.toml").write_text(manifest)
    shutil.copyfile(CONSUMERS / "edit_compare.Cargo.lock", project / "Cargo.lock")
    subprocess.run(["cargo", "build", "--release", "--offline", "--locked", "--manifest-path", str(project / "Cargo.toml")], capture_output=True, check=True)
    binary = project / "target/release/exact-edit-comparison.exe"
    metadata = {"variant": variant, "source_sha256": {str(path): sha(path) for path in sorted((library / "src").glob("*.rs"))},
                "harness_sha256": sha(project / "src/main.rs"), "lock_sha256": sha(project / "Cargo.lock"),
                "binary_sha256": sha(binary), "binary": str(binary)}
    binary.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
    return metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--variants", default="baseline,children_words,children_words_frontier,children_words_frontier_cross")
    parser.add_argument("--output", type=Path, default=HERE / "fast_path_churn.json")
    args = parser.parse_args()
    variants = args.variants.split(",")
    builds = [json.loads((HERE / ".cache" / ("fast_churn_" + variant) / "target/release/exact-edit-comparison.json").read_text())
              if args.run else build(variant) for variant in variants]
    if args.build_only:
        print(json.dumps(builds, indent=2))
        return
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite {args.output}")
    result = {"purpose": __doc__, "seed": 1, "invocation": 0, "samples_per_cell": 3, "builds": builds, "runs": []}
    for metadata in builds:
        binary = Path(metadata["binary"])
        assert sha(binary) == metadata["binary_sha256"]
        cells = {}
        for line in subprocess.check_output([str(binary), "1", "0"], text=True).splitlines():
            _, size, shared, arm, _, ms = line.split("\t")
            cell = cells.setdefault((size, shared, arm), {"bytes": int(size), "shared": shared == "true", "arm": arm, "samples_ms": []})
            cell["samples_ms"].append(float(ms))
        for cell in cells.values():
            cell["median_ms"] = statistics.median(cell["samples_ms"])
        result["runs"].append({"variant": metadata["variant"], "cells": list(cells.values())})
        print(metadata["variant"], flush=True)
    args.output.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
