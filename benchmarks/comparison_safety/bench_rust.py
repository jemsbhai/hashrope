"""Compile/run bench_rust.rs against a specified packages/rust directory.

Only a temporary Cargo project is written. No dependencies beyond hashrope/std.
Use --build-only to return a binary path for alternating baseline/current runs.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile


def parse_output(output):
    records, memory = {}, {}
    for line in output.splitlines():
        parts = line.split("\t")
        if parts[0] == "time":
            _, name, _, iterations, ns = parts
            record = records.setdefault(name, {"name": name, "iterations_per_sample": int(iterations),
                                               "ns_per_op": []})
            record["ns_per_op"].append(float(ns))
        elif parts[0] == "memory":
            memory[parts[1]] = int(parts[2])
    for record in records.values():
        record["median_ns"] = statistics.median(record["ns_per_op"])
    return {"measurements": list(records.values()), "memory": memory}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--build-only", action="store_true")
    args = parser.parse_args()
    directory = Path(tempfile.mkdtemp(prefix="hashrope-benchmark-"))
    (directory / "src").mkdir()
    shutil.copyfile(Path(__file__).with_suffix(".rs"), directory / "src/main.rs")
    # Forward slashes and JSON string syntax are valid TOML basic strings here.
    package_path = json.dumps(args.package.resolve().as_posix())
    (directory / "Cargo.toml").write_text(
        '[package]\nname="hashrope-comparison-benchmark"\nversion="0.0.0"\nedition="2021"\n'
        '[features]\nexact=[]\n[dependencies]\nhashrope={path=' + package_path + '}\n', encoding="utf-8")
    command = ["cargo", "build", "--release", "--offline", "--manifest-path", str(directory / "Cargo.toml")]
    if args.exact:
        command += ["--features", "exact"]
    subprocess.run(command, check=True)
    binary = directory / "target/release" / ("hashrope-comparison-benchmark.exe" if os.name == "nt" else "hashrope-comparison-benchmark")
    if args.build_only:
        print(binary)
        return
    result = parse_output(subprocess.check_output([str(binary)], text=True))
    result["package"] = str(args.package.resolve())
    result["runtime"] = subprocess.check_output(["rustc", "--version"], text=True).strip()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
