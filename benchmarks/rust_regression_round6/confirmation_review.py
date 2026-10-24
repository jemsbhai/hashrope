"""Independent audit of the separate six-pair final-binary confirmation."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics as stats
import zipfile
from final_review import read, sha, close

HERE = Path(__file__).resolve().parent
STUDY = HERE / "final_confirmation"
REFERENCE = HERE / "final_comparison"


def audit():
    manifest, results, summary = [read(STUDY / (name + ".json")) for name in ["manifest", "results", "summary"]]
    built = read(REFERENCE / "build.json")
    for name, digest in manifest["reference_artifacts"].items():
        assert sha(REFERENCE / name) == digest
    for name, digest in manifest["drivers"].items():
        assert sha(STUDY / name) == digest
    assert results["manifest_sha256"] == sha(STUDY / "manifest.json")
    assert summary["results_sha256"] == sha(STUDY / "results.json")
    assert summary["analysis_sha256"] == sha(STUDY / "run.py")
    assert built["binaries_archive_sha256"] == sha(REFERENCE / "binaries.zip")
    with zipfile.ZipFile(REFERENCE / "binaries.zip") as archive:
        assert len(archive.namelist()) == 4
        for label in ["baseline", "candidate"]:
            for kind in ["micro", "churn"]:
                digest = built[label][kind]["binary_sha256"]
                assert hashlib.sha256(archive.read(f"{label}/{kind}.exe")).hexdigest() == digest
                assert sha(STUDY / ".cache" / label / (kind + ".exe")) == digest
    assert manifest["protocol"]["churn_seeds"] == [1, 2, 3, 1, 2, 3]
    orders = [["baseline", "candidate"] if pair % 2 == 0 else ["candidate", "baseline"] for pair in range(6)]
    assert manifest["protocol"]["orders"] == orders
    assert [(r["kind"], r["pair"], r["label"]) for r in results["runs"]] == [
        (kind, pair, label) for kind in ["micro", "churn"] for pair, order in enumerate(orders) for label in order]
    index, samples = {}, 0
    for row in results["runs"]:
        kind, pair, label = row["kind"], row["pair"], row["label"]
        path = STUDY / row["raw"]
        assert sha(path) == row["raw_sha256"] and row == read(path.with_suffix(".json"))
        assert sha(path.with_suffix(".stderr")) == row["stderr_sha256"]
        assert path.with_suffix(".stderr").read_bytes() == b""
        assert row["returncode"] == 0 and not row["timed_out"] and row["order"] == orders[pair]
        assert row["seed"] == (pair % 3 + 1 if kind == "churn" else None)
        assert row["binary_sha256"] == built[label][kind]["binary_sha256"]
        assert row["arguments"] == [str(STUDY / ".cache" / label / (kind + ".exe"))] + ([str(row["seed"]), "0"] if kind == "churn" else [])
        cells, memory = {}, {}
        for line in path.read_text().splitlines():
            fields = line.split("\t")
            if kind == "micro":
                if fields[0] == "memory":
                    assert fields[1] not in memory
                    memory[fields[1]] = int(fields[2])
                    continue
                tag, name, sample, iterations, value = fields
                assert tag == "time" and int(iterations) == row["cells"][name]["iterations"]
            else:
                seed, length, shared, arm, sample, value = fields
                assert int(seed) == row["seed"]
                name = f"{length}/{shared}/{arm}"
            values = cells.setdefault(name, [])
            assert int(sample) == len(values)
            value = float(value)
            assert math.isfinite(value) and value > 0
            values.append(value)
            samples += 1
        assert cells.keys() == row["cells"].keys() and len(cells) == (42 if kind == "micro" else 24)
        unit = "ns" if kind == "micro" else "ms"
        for name, values in cells.items():
            assert len(values) == (7 if kind == "micro" else 3)
            assert values == row["cells"][name]["samples_" + unit]
            assert stats.median(values) == row["cells"][name]["median_" + unit]
            index[kind, pair, label, name] = stats.median(values)
        if kind == "micro":
            assert memory == row["memory"] == summary["memory"][label]
    cells_verified = 0
    for kind, cells in summary["cells"].items():
        for name, estimate in cells.items():
            values = {label: [index[kind, pair, label, name] for pair in range(6)] for label in ["baseline", "candidate"]}
            assert values == estimate["process_medians"]
            for label in values:
                assert stats.median(values[label]) == estimate["median_of_process_medians"][label]
            ratios = [b / a for a, b in zip(values["baseline"], values["candidate"])]
            assert ratios == estimate["paired_ratios"]
            logs = [math.log(x) for x in ratios]
            center = stats.mean(logs)
            margin = 2.5705818366147395 * stats.stdev(logs) / math.sqrt(6)
            close(math.exp(center), estimate["paired_ratio_geomean"])
            assert [min(ratios), max(ratios)] == estimate["paired_ratio_range"]
            assert sum(x > 1 for x in ratios) == estimate["candidate_slower_pairs"]
            for a, b in zip([math.exp(center - margin), math.exp(center + margin)], estimate["exploratory_95_log_t_interval"]):
                close(a, b)
            for label, value in estimate["order_geomeans"].items():
                close(math.exp(stats.mean(math.log(ratios[pair]) for pair in range(6) if orders[pair][0] == label)), value)
            if kind == "churn":
                for seed, value in estimate["seed_geomeans"].items():
                    close(math.exp(stats.mean(math.log(ratios[pair]) for pair in range(6) if pair % 3 + 1 == int(seed))), value)
            cells_verified += 1
    rows = list(csv.DictReader((STUDY / "all_cells.csv").open(encoding="utf-8", newline="")))
    assert len(rows) == cells_verified == 66
    for row in rows:
        e = summary["cells"][row["cohort"]][row["cell"]]
        for field, value in [("baseline_median", e["median_of_process_medians"]["baseline"]),
                             ("candidate_median", e["median_of_process_medians"]["candidate"]),
                             ("paired_ratio", e["paired_ratio_geomean"]), ("ci_low", e["exploratory_95_log_t_interval"][0]),
                             ("ci_high", e["exploratory_95_log_t_interval"][1]), ("slower_pairs", e["candidate_slower_pairs"])]:
            close(float(row[field]), value)
    assert len(results["runs"]) == summary["processes"] == 24 and samples == summary["samples"] == 4392
    output = {"audit": "Independent raw/provenance/statistical verification of the separate confirmation; no builds or timings.",
              "script_sha256": sha(Path(__file__)), "helper_sha256": sha(HERE / "final_review.py"),
              "summary_sha256": sha(STUDY / "summary.json"),
              "verified": {"processes": 24, "samples": samples, "paired_vectors": cells_verified, "csv_rows": len(rows), "retained_executables": 4}}
    with (HERE / "confirmation_review.json").open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps(output["verified"]))


if __name__ == "__main__":
    audit()
