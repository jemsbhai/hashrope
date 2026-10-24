"""Fixed 12-pair replication of verified retained original binaries."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import subprocess
import zipfile
from compare_controlled import available_affinity

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    output = HERE / "results/original_replication"
    if output.exists():
        raise SystemExit(f"Preserve existing result: {output}")
    metadata = json.loads((HERE / "retained_original_binaries.json").read_text())
    prior = json.loads((ROOT / "benchmarks/exact_optimization_round3/results/microbench.json").read_text())
    parser_path = ROOT / "benchmarks/comparison_safety/bench_rust.py"
    original_harness = parser_path.with_suffix(".rs")
    spec = importlib.util.spec_from_file_location("original_parser", parser_path)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    for variant, item in metadata.items():
        assert sha(item["binary"]) == item["binary_sha256"] == prior["binary_sha256"][variant]
        assert item["rust_harness_sha256"] == sha(original_harness) == prior["harness_sha256"]["bench_rust.rs"]
    output.mkdir(parents=True)
    # Preserve the exact machine code used, not just a future rebuild recipe.
    with zipfile.ZipFile(output / "retained_binaries.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for variant, item in metadata.items():
            archive.write(item["binary"], f"{variant}.exe")
        for path in [Path(__file__), HERE / "compare_controlled.py", HERE / "original_replication_plan.md", HERE / "retained_original_binaries.json", original_harness, parser_path]:
            archive.write(path, path.name)
    process_mask, system_mask = available_affinity()
    results = {"protocol_sha256": sha(HERE / "original_replication_plan.md"), "driver_sha256": sha(__file__),
               "parent_affinity_mask": process_mask, "system_affinity_mask": system_mask,
               "affinity_helper_sha256": sha(HERE / "compare_controlled.py"),
               "original_micro_sha256": sha(ROOT / "benchmarks/exact_optimization_round3/results/microbench.json"),
               "binary_archive_sha256": sha(output / "retained_binaries.zip"), "binaries": metadata,
               "uncertainty": "12 process pairs; exploratory per-cell 95% log-t intervals, not adjusted for multiple comparisons.",
               "runs": []}
    for pair in range(1, 13):
        order = ["baseline", "optimized"] if pair % 2 else ["optimized", "baseline"]
        for variant in order:
            print(f"Original binary replication pair={pair}/12 variant={variant}", flush=True)
            raw = subprocess.check_output([metadata[variant]["binary"]], text=True)
            path = output / f"pair_{pair:02}_{variant}.tsv"
            path.write_text(raw, encoding="utf-8")
            result = parser.parse_output(raw)
            assert len(result["measurements"]) == 42
            assert all(len(row["ns_per_op"]) == 7 for row in result["measurements"])
            result.update(pair=pair, order=order, variant=variant, raw_sha256=sha(path), raw_file=path.name)
            results["runs"].append(result)
            (output / "partial.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    summary = []
    names = [row["name"] for row in results["runs"][0]["measurements"]]
    for name in names:
        values = []
        for pair in range(1, 13):
            medians = {run["variant"]: next(row["median_ns"] for row in run["measurements"] if row["name"] == name)
                       for run in results["runs"] if run["pair"] == pair}
            values.append(dict(pair=pair, baseline_ns=medians["baseline"], optimized_ns=medians["optimized"],
                               ratio=medians["optimized"] / medians["baseline"]))
        logs = [math.log(row["ratio"]) for row in values]
        mean = statistics.mean(logs)
        margin = 2.200985160082949 * statistics.stdev(logs) / math.sqrt(len(logs))
        summary.append({"name": name, "process_pairs": values, "geomean_ratio": math.exp(mean),
                        "ratio_min": min(row["ratio"] for row in values), "ratio_max": max(row["ratio"] for row in values),
                        "exploratory_95pct_log_t_interval": [math.exp(mean-margin), math.exp(mean+margin)],
                        "baseline_mean_ns": statistics.mean(row["baseline_ns"] for row in values),
                        "optimized_mean_ns": statistics.mean(row["optimized_ns"] for row in values)})
    results["summary"] = summary
    (output / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("All 12 original-binary pairs complete; every original sample retained.")


if __name__ == "__main__":
    main()
