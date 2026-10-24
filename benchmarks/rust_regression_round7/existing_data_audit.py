"""Read-only analysis of both retained round-six cohorts; performs no timing."""
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics as stats
import zipfile

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "rust_regression_round6"
TARGETS = {
    "100K_independent_exact_churn": ("churn", "100000/false/hashrope_exact"),
    "1MiB_independent_equal_lcp": ("micro", "compare/1048576/independent_equal/exact_lcp"),
}
CONTROLS = [
    ("churn", "100000/false/hashrope_candidate"),
    ("churn", "100000/false/ropey_exact"),
    ("churn", "100000/false/ropey_rehash_candidate"),
    ("churn", "100000/true/hashrope_exact"),
    ("churn", "10000/false/hashrope_exact"),
    ("churn", "1000000/false/hashrope_exact"),
    ("micro", "compare/1048576/independent_equal/exact_equal"),
    ("micro", "compare/1048576/late_mismatch/exact_lcp"),
    ("micro", "compare/1048576/late_mismatch/exact_equal"),
    ("micro", "compare/1048576/independent_equal/hash_only_lcp"),
    ("micro", "compare/65536/independent_equal/exact_lcp"),
    ("micro", "core/concat_balanced_64KiB"),
]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geomean(values):
    return math.exp(stats.mean(math.log(x) for x in values))


def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12), (a, b)


def main():
    built = read(PREVIOUS / "final_comparison/build.json")
    source_manifest = read(PREVIOUS / "final_comparison/manifest.json")
    inputs = {}
    for directory in ["final_comparison", "final_confirmation"]:
        for name in ["manifest.json", "results.json", "summary.json", "run.py", "PLAN.md"]:
            inputs[f"{directory}/{name}"] = sha(PREVIOUS / directory / name)
    for name in ["baseline.zip", "candidate.zip", "binaries.zip", "build.json", "drivers.zip"]:
        inputs[f"final_comparison/{name}"] = sha(PREVIOUS / "final_comparison" / name)
    for name, digest in source_manifest["archives"].items():
        assert sha(PREVIOUS / "final_comparison" / name) == digest
    with zipfile.ZipFile(PREVIOUS / "final_comparison/binaries.zip") as archive:
        for label in ["baseline", "candidate"]:
            for kind in ["micro", "churn"]:
                assert hashlib.sha256(archive.read(f"{label}/{kind}.exe")).hexdigest() == built[label][kind]["binary_sha256"]
    confirmation_manifest = read(PREVIOUS / "final_confirmation/manifest.json")
    for name, digest in confirmation_manifest["reference_artifacts"].items():
        assert sha(PREVIOUS / "final_comparison" / name) == digest
    output = {"scope": "Separate descriptive audit of nine-pair and six-pair cohorts; no pooling, filtering, new timing or causal adjustment.",
              "script_sha256": sha(Path(__file__)), "input_sha256": inputs,
              "binary_sha256": {label: {kind: built[label][kind]["binary_sha256"] for kind in ["micro", "churn"]}
                                  for label in ["baseline", "candidate"]},
              "phases": {}}
    total_samples = total_processes = total_vectors = 0
    for directory, pairs, t in [("final_comparison", 9, 2.306004135204166),
                                 ("final_confirmation", 6, 2.5705818366147395)]:
        path = PREVIOUS / directory
        results, summary = read(path / "results.json"), read(path / "summary.json")
        assert summary["results_sha256"] == sha(path / "results.json")
        assert results["manifest_sha256"] == sha(path / "manifest.json")
        assert len(results["runs"]) == pairs * 4
        index, positions = {}, {}
        for row in results["runs"]:
            kind, pair, label = row["kind"], row["pair"], row["label"]
            raw = path / row["raw"]
            assert sha(raw) == row["raw_sha256"] and row == read(raw.with_suffix(".json"))
            assert row["returncode"] == 0 and not row["timed_out"]
            assert raw.with_suffix(".stderr").read_bytes() == b""
            assert row["binary_sha256"] == built[label][kind]["binary_sha256"]
            assert row["order"] == (["baseline", "candidate"] if pair % 2 == 0 else ["candidate", "baseline"])
            assert row["seed"] == (pair % 3 + 1 if kind == "churn" else None)
            parsed, lines = {}, defaultdict(list)
            for line_number, line in enumerate(raw.read_text().splitlines(), 1):
                fields = line.split("\t")
                if kind == "micro":
                    if fields[0] == "memory":
                        continue
                    tag, name, sample, iterations, value = fields
                    assert tag == "time" and int(iterations) == row["cells"][name]["iterations"]
                else:
                    seed, size, shared, arm, sample, value = fields
                    assert int(seed) == row["seed"]
                    name = f"{size}/{shared}/{arm}"
                values = parsed.setdefault(name, [])
                assert int(sample) == len(values)
                values.append(float(value))
                lines[name].append(line_number)
                total_samples += 1
            assert list(parsed) == list(row["cells"]) and len(parsed) == (42 if kind == "micro" else 24)
            unit = "ns" if kind == "micro" else "ms"
            for name, values in parsed.items():
                assert len(values) == (7 if kind == "micro" else 3)
                assert values == row["cells"][name]["samples_" + unit]
                assert stats.median(values) == row["cells"][name]["median_" + unit]
            index[kind, pair, label] = row
            positions[kind, pair, label] = lines
            total_processes += 1
        for kind, cells in summary["cells"].items():
            unit = "ns" if kind == "micro" else "ms"
            for name, estimate in cells.items():
                values = {label: [index[kind, pair, label]["cells"][name]["median_" + unit] for pair in range(pairs)]
                          for label in ["baseline", "candidate"]}
                assert values == estimate["process_medians"]
                ratios = [b / a for a, b in zip(values["baseline"], values["candidate"])]
                assert ratios == estimate["paired_ratios"]
                close(geomean(ratios), estimate["paired_ratio_geomean"])
                logs = list(map(math.log, ratios))
                margin = t * stats.stdev(logs) / math.sqrt(pairs)
                for actual, expected in zip(estimate["exploratory_95_log_t_interval"],
                                             [math.exp(stats.mean(logs) - margin), math.exp(stats.mean(logs) + margin)]):
                    close(actual, expected)
                assert estimate["candidate_slower_pairs"] == sum(x > 1 for x in ratios)
                total_vectors += 1
        phase = {"process_pairs_per_cohort": pairs, "targets": {}, "nearby_controls": {}}
        for target, (kind, name) in TARGETS.items():
            unit = "ns" if kind == "micro" else "ms"
            factor_us = 0.001 if kind == "micro" else 1000
            estimate = summary["cells"][kind][name]
            target_result = {"cell": name, "original_unit": unit, "summary": estimate, "pairs": []}
            raw_vectors = {label: [] for label in ["baseline", "candidate"]}
            for pair in range(pairs):
                item = {"pair": pair, "first_arm": index[kind, pair, "baseline"]["order"][0],
                        "seed": index[kind, pair, "baseline"]["seed"], "arms": {}}
                for label in ["baseline", "candidate"]:
                    row = index[kind, pair, label]
                    values = row["cells"][name]["samples_" + unit]
                    raw_vectors[label].append(values)
                    arm = {"samples_us": [x * factor_us for x in values], "median_us": stats.median(values) * factor_us,
                           "timing_cell_position_one_based": list(row["cells"]).index(name) + 1,
                           "raw_sample_lines_one_based": positions[kind, pair, label][name],
                           "started_utc": row["started_utc"], "finished_utc": row["finished_utc"]}
                    if kind == "churn":
                        arm["position_within_100K_independent_arms_one_based"] = [x for x in row["cells"] if x.startswith("100000/false/")].index(name) + 1
                    item["arms"][label] = arm
                item["paired_ratio"] = estimate["paired_ratios"][pair]
                target_result["pairs"].append(item)
            target_result["ratio_of_arm_medians"] = estimate["median_of_process_medians"]["candidate"] / estimate["median_of_process_medians"]["baseline"]
            target_result["raw_sample_position_paired_geomeans"] = [
                geomean(b[j] / a[j] for a, b in zip(raw_vectors["baseline"], raw_vectors["candidate"]))
                for j in range(len(raw_vectors["baseline"][0]))]
            target_result["within_process_later_earlier_geomeans"] = {
                label: geomean((stats.median(v[3:]) / stats.median(v[:3]) if kind == "micro" else v[2] / v[0]) for v in values)
                for label, values in raw_vectors.items()}
            target_result["sample_position_scope"] = "Descriptive correlated inner samples only; no independent-replicate inference, sample removal, replacement primary estimator or causal attribution. Micro later/earlier is median(samples3..6)/median(samples0..2); churn is sample2/sample0."
            phase["targets"][target] = target_result
        for kind, name in CONTROLS:
            phase["nearby_controls"][name] = summary["cells"][kind][name]
        output["phases"][directory] = phase
    assert (total_processes, total_samples, total_vectors) == (60, 10980, 132)
    output["verified"] = dict(processes=total_processes, samples=total_samples, paired_vectors=total_vectors)
    with (HERE / "existing_data_audit.json").open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps(output["verified"]))


if __name__ == "__main__":
    main()
