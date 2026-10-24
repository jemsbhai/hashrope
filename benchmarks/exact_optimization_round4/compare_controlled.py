"""Build and run a separate paired comparison diagnostic; never edits core sources."""
import argparse
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CORES = {
    "round2": ROOT / "benchmarks/exact_optimization_round3/.cache/baseline/hashrope/packages/rust",
    "round3": HERE / ".cache/baseline/hashrope/packages/rust",
}
BUILDS = {"old_new": {"a": "round2", "b": "round3"},
          "new_old": {"a": "round3", "b": "round2"},
          "same_same": {"a": "round3", "b": "round3"}}
CACHE = HERE / ".cache/controlled_v2"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_map(root):
    return {str(path.relative_to(root)).replace("\\", "/"): sha(path)
            for path in sorted(root.rglob("*.rs")) if "target" not in path.parts}


def build():
    if CACHE.exists():
        raise SystemExit(f"Preserve existing build: {CACHE}")
    CACHE.mkdir(parents=True)
    metadata = {"runtime": subprocess.check_output(["rustc", "-Vv"], text=True),
                "harness_sha256": {p.name: sha(p) for p in [Path(__file__), HERE / "compare_controlled.rs"]},
                "cores": {label: {"path": str(path), "source_sha256": source_map(path),
                                   "cargo_sha256": sha(path / "Cargo.toml")} for label, path in CORES.items()},
                "builds": {}}
    old = (CORES["round2"] / "src/rope.rs").read_text()
    new = (CORES["round3"] / "src/rope.rs").read_text()
    portions = {"public_exact_apis": ("    pub fn equal(&self", "    /// Weight (number of leaves)"),
                "comparison_engine": ("struct RangeStart", "#[cfg(test)]")}
    audit = {}
    for name, (start, end) in portions.items():
        left = old[old.index(start):old.index(end, old.index(start))]
        right = new[new.index(start):new.index(end, new.index(start))]
        audit[name] = {"identical": left == right, "round2_sha256": hashlib.sha256(left.encode()).hexdigest(),
                       "round3_sha256": hashlib.sha256(right.encode()).hexdigest(), "characters": len(left)}
        assert left == right, f"Expected unchanged comparator source: {name}"
    metadata["source_audit"] = audit
    for label, mapping in BUILDS.items():
        project = CACHE / label
        (project / "src").mkdir(parents=True)
        shutil.copyfile(HERE / "compare_controlled.rs", project / "src/main.rs")
        for slot, core in mapping.items():
            destination = project / f"slot_{slot}"
            shutil.copytree(CORES[core] / "src", destination / "src")
            # Fixed crate slot names permit a swapped assignment and an A/A null.
            (destination / "Cargo.toml").write_text(
                f'[package]\nname="hashrope-diagnostic-slot-{slot}"\nversion="0.0.0"\nedition="2021"\n', encoding="utf-8")
        (project / "Cargo.toml").write_text(
            '[package]\nname="hashrope-controlled-comparison"\nversion="0.0.0"\nedition="2021"\n'
            '[dependencies]\nslot_a={package="hashrope-diagnostic-slot-a",path="slot_a"}\n'
            'slot_b={package="hashrope-diagnostic-slot-b",path="slot_b"}\n', encoding="utf-8")
        subprocess.run(["cargo", "build", "--release", "--offline", "--manifest-path", str(project / "Cargo.toml")], check=True)
        executable = project / "target/release" / ("hashrope-controlled-comparison.exe" if os.name == "nt" else "hashrope-controlled-comparison")
        metadata["builds"][label] = {"mapping": mapping, "binary": str(executable), "binary_sha256": sha(executable),
            "copied_core_sha256": {slot: source_map(project / f"slot_{slot}") for slot in mapping},
            "cargo_lock_sha256": sha(project / "Cargo.lock")}
    (HERE / "controlled_build.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with zipfile.ZipFile(HERE / "controlled_sources.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for path in [Path(__file__), HERE / "compare_controlled.rs", HERE / "controlled_plan.md", HERE / "controlled_build.json"]:
            archive.write(path, path.name)
        for label, core in CORES.items():
            archive.write(core / "Cargo.toml", f"{label}/Cargo.toml")
            for path in sorted((core / "src").rglob("*.rs")):
                archive.write(path, f"{label}/{path.relative_to(core).as_posix()}")
    print("All diagnostic builds ready; no measurements have run.")


def parse(raw):
    result = {"metadata": [], "fixtures": [], "blocks": []}
    for line in raw.splitlines():
        parts = line.split("\t")
        if parts[0] == "meta":
            result["metadata"].append(parts[1:])
        elif parts[0] == "fixture":
            _, slot, case, nodes, length, expected, root, leaf = parts
            result["fixtures"].append(dict(slot=slot, case=case, nodes=int(nodes), length=int(length),
                expected=int(expected), root_address=int(root), first_leaf_address=int(leaf)))
        elif parts[0] == "time":
            _, case, api, case_order, api_order, position, slot, iterations, ns, answer, cpu0, cpu1 = parts
            result["blocks"].append(dict(case=case, api=api, case_order=int(case_order), api_order=int(api_order),
                position=int(position), slot=slot, iterations=int(iterations), elapsed_ns=int(ns),
                ns_per_op=int(ns)/int(iterations), answer=int(answer), processor_before=int(cpu0), processor_after=int(cpu1)))
        else:
            raise ValueError(line)
    assert len(result["fixtures"]) == 24
    assert len(result["blocks"]) == 96
    for case in {row["case"] for row in result["blocks"]}:
        for api in ["lcp", "equal"]:
            rows = [row for row in result["blocks"] if row["case"] == case and row["api"] == api]
            assert sorted(row["slot"] for row in rows) == ["a", "a", "b", "b"]
            assert len({row["iterations"] for row in rows}) == 1
    return result


def summarize(runs):
    output = []
    for mode in ["unpinned", "pinned"]:
        for build_label, mapping in BUILDS.items():
            matching = [run for run in runs if run["build"] == build_label and run["affinity_mode"] == mode]
            if not matching:
                continue
            assert sorted(run["seed"] for run in matching) == list(range(1, 7))
            keys = sorted({(row["case"], row["api"]) for row in matching[0]["blocks"]})
            for case, api in keys:
                values = []
                for run in sorted(matching, key=lambda row: row["seed"]):
                    rows = [row for row in run["blocks"] if row["case"] == case and row["api"] == api]
                    medians = {slot: statistics.median(row["ns_per_op"] for row in rows if row["slot"] == slot) for slot in ["a", "b"]}
                    old_slot, new_slot = ("b", "a") if build_label == "new_old" else ("a", "b")
                    values.append({"seed": run["seed"], "denominator_ns": medians[old_slot],
                                   "numerator_ns": medians[new_slot], "ratio": medians[new_slot] / medians[old_slot]})
                logs = [math.log(value["ratio"]) for value in values]
                mean, sd = statistics.mean(logs), statistics.stdev(logs)
                margin = 2.570581835636305 * sd / math.sqrt(len(logs))  # t(5), 97.5th percentile
                output.append({"affinity_mode": mode, "build": build_label, "case": case, "api": api,
                    "ratio_definition": "slot_b / slot_a (identical core)" if build_label == "same_same" else "round3 / round2",
                    "process_pairs": values, "geomean_ratio": math.exp(mean),
                    "ratio_min": min(value["ratio"] for value in values), "ratio_max": max(value["ratio"] for value in values),
                    "exploratory_95pct_log_t_interval": [math.exp(mean-margin), math.exp(mean+margin)],
                    "numerator_mean_ns": statistics.mean(value["numerator_ns"] for value in values),
                    "denominator_mean_ns": statistics.mean(value["denominator_ns"] for value in values)})
    return output


def available_affinity():
    if os.name != "nt":
        raise SystemExit("This diagnostic's pinned arm currently requires Windows.")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    kernel.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
    process, system = ctypes.c_size_t(), ctypes.c_size_t()
    if not kernel.GetProcessAffinityMask(kernel.GetCurrentProcess(), ctypes.byref(process), ctypes.byref(system)):
        raise ctypes.WinError(ctypes.get_last_error())
    return process.value, system.value


def run(output, block_ms, warm_ms):
    if output.exists():
        raise SystemExit(f"Use a fresh diagnostic output directory: {output}")
    output.mkdir(parents=True)
    build_metadata = json.loads((HERE / "controlled_build.json").read_text())
    assert all(sha(HERE / name) == value for name, value in build_metadata["harness_sha256"].items())
    for label, core in build_metadata["cores"].items():
        assert source_map(Path(core["path"])) == core["source_sha256"]
    for label, item in build_metadata["builds"].items():
        assert sha(item["binary"]) == item["binary_sha256"]
    available, system = available_affinity()
    pin_mask = available & -available
    results = {"protocol": "Separate controlled diagnostic; original benchmark unchanged.",
        "build": build_metadata, "source_archive_sha256": sha(HERE / "controlled_sources.zip"),
        "plan_sha256": sha(HERE / "controlled_plan.md"), "runtime": sys.version,
        "platform": platform.platform(), "initial_affinity_mask": available,
        "system_affinity_mask": system, "pinned_mask": pin_mask,
        "block_target_ms": block_ms, "paired_warmup_target_ms": warm_ms,
        "uncertainty": "Six process-paired ratios per build/affinity; exploratory per-cell log-t intervals, not multiplicity adjusted; inner blocks are not independent replicates.",
        "runs": []}
    for seed in range(1, 7):
        builds = list(BUILDS)
        rotation = (seed - 1) % len(builds)
        builds = builds[rotation:] + builds[:rotation]
        if seed > 3:
            builds.reverse()
        modes = ["unpinned", "pinned"] if seed % 2 else ["pinned", "unpinned"]
        for mode in modes:
            for label in builds:
                print(f"Controlled diagnostic seed={seed} affinity={mode} build={label}", flush=True)
                command = [build_metadata["builds"][label]["binary"], str(seed),
                           str(pin_mask if mode == "pinned" else 0), str(block_ms), str(warm_ms)]
                raw = subprocess.check_output(command, text=True)
                filename = f"{label}_{mode}_s{seed}.tsv"
                (output / filename).write_text(raw, encoding="utf-8")
                result = parse(raw)
                result.update(seed=seed, affinity_mode=mode, build=label, raw_file=filename,
                              raw_sha256=sha(output / filename), command=command)
                results["runs"].append(result)
                (output / "partial.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    results["summary"] = summarize(results["runs"])
    (output / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("Controlled diagnostic complete; all blocks and oracles retained.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["build", "run"])
    parser.add_argument("--output", type=Path, default=HERE / "results/controlled_initial")
    parser.add_argument("--block-ms", type=int, default=10)
    parser.add_argument("--warm-ms", type=int, default=10)
    args = parser.parse_args()
    if args.action == "build":
        build()
    else:
        run(args.output.resolve(), args.block_ms, args.warm_ms)


if __name__ == "__main__":
    main()
