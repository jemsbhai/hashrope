"""Independent raw/provenance/statistical audit of the final nine-pair cohort."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics as stats
import zipfile

HERE = Path(__file__).resolve().parent
STUDY = HERE / "final_comparison"
ROOT = HERE.parents[1]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12), (a, b)


def audit():
    manifest, built, results, summary = [read(STUDY / (name + ".json"))
                                         for name in ["manifest", "build", "results", "summary"]]
    assert sha(HERE / "baseline_manifest.json") == manifest["parent_manifest_sha256"]
    assert sha(HERE.parent / "rust_lcp_round5/run.py") == manifest["imported_builder_sha256"]
    for name, digest in manifest["archives"].items():
        assert sha(STUDY / name) == digest
    with zipfile.ZipFile(STUDY / "drivers.zip") as archive:
        for name, digest in manifest["drivers"].items():
            assert sha(STUDY / name) == digest == hashlib.sha256(archive.read(name)).hexdigest()
    for label, sources in manifest["sources"].items():
        with zipfile.ZipFile(STUDY / (label + ".zip")) as archive:
            assert set(archive.namelist()) == {"hashrope/" + name for name in sources}
            for name, digest in sources.items():
                assert hashlib.sha256(archive.read("hashrope/" + name)).hexdigest() == digest
    assert results["manifest_sha256"] == built["environment"]["manifest_sha256"] == sha(STUDY / "manifest.json")
    assert results["build_sha256"] == sha(STUDY / "build.json")
    assert summary["results_sha256"] == sha(STUDY / "results.json")
    assert summary["analysis_sha256"] == sha(STUDY / "run.py")
    assert built["binaries_archive_sha256"] == sha(STUDY / "binaries.zip")
    with zipfile.ZipFile(STUDY / "binaries.zip") as archive:
        assert len(archive.namelist()) == 4
        for label in ["baseline", "candidate"]:
            for kind in ["micro", "churn"]:
                info = built[label][kind]
                assert hashlib.sha256(archive.read(f"{label}/{kind}.exe")).hexdigest() == info["binary_sha256"]
                assert sha(Path(info["binary"])) == info["binary_sha256"]
                assert info["source_sha256"] == manifest["drivers"]["workloads/" + ("bench_rust.rs" if kind == "micro" else "edit_compare.rs")]
                assert info["cargo_lock_sha256"] == manifest["drivers"][f"workloads/{kind}.Cargo.lock"]
                project = STUDY / ".cache/build" / f"{label}_{kind}"
                for name, key in [("src/main.rs", "source_sha256"), ("Cargo.toml", "cargo_manifest_sha256"), ("Cargo.lock", "cargo_lock_sha256")]:
                    assert sha(project / name) == info[key]
    expected_order = [(kind, pair, label) for kind in ["micro", "churn"]
                      for pair, order in enumerate(manifest["protocol"]["orders"]) for label in order]
    assert [(r["kind"], r["pair"], r["label"]) for r in results["runs"]] == expected_order
    index, samples = {}, 0
    for row in results["runs"]:
        kind, pair, label = row["kind"], row["pair"], row["label"]
        path = STUDY / row["raw"]
        assert sha(path) == row["raw_sha256"] and row == read(path.with_suffix(".json"))
        assert sha(path.with_suffix(".stderr")) == row["stderr_sha256"]
        assert path.with_suffix(".stderr").read_bytes() == b""
        assert row["returncode"] == 0 and not row["timed_out"]
        assert row["order"] == manifest["protocol"]["orders"][pair]
        assert row["seed"] == (pair % 3 + 1 if kind == "churn" else None)
        assert row["binary_sha256"] == built[label][kind]["binary_sha256"]
        assert row["arguments"] == [built[label][kind]["binary"]] + ([str(row["seed"]), "0"] if kind == "churn" else [])
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
            values = {label: [index[kind, pair, label, name] for pair in range(9)] for label in ["baseline", "candidate"]}
            assert values == estimate["process_medians"]
            for label in values:
                assert stats.median(values[label]) == estimate["median_of_process_medians"][label]
            ratios = [b / a for a, b in zip(values["baseline"], values["candidate"])]
            assert ratios == estimate["paired_ratios"]
            logs = [math.log(x) for x in ratios]
            center = stats.mean(logs)
            margin = 2.306004135204166 * stats.stdev(logs) / 3
            close(math.exp(center), estimate["paired_ratio_geomean"])
            assert [min(ratios), max(ratios)] == estimate["paired_ratio_range"]
            assert sum(x > 1 for x in ratios) == estimate["candidate_slower_pairs"]
            for a, b in zip([math.exp(center - margin), math.exp(center + margin)], estimate["exploratory_95_log_t_interval"]):
                close(a, b)
            for label, value in estimate["order_geomeans"].items():
                close(math.exp(stats.mean(math.log(ratios[pair]) for pair in range(9)
                                          if manifest["protocol"]["orders"][pair][0] == label)), value)
            if kind == "churn":
                for seed, value in estimate["seed_geomeans"].items():
                    close(math.exp(stats.mean(math.log(ratios[pair]) for pair in range(9) if pair % 3 + 1 == int(seed))), value)
            cells_verified += 1
    csv_rows = list(csv.DictReader((STUDY / "all_cells.csv").open(encoding="utf-8", newline="")))
    assert len(csv_rows) == cells_verified == 66
    for row in csv_rows:
        e = summary["cells"][row["cohort"]][row["cell"]]
        for field, value in [("baseline_median", e["median_of_process_medians"]["baseline"]),
                             ("candidate_median", e["median_of_process_medians"]["candidate"]),
                             ("paired_ratio", e["paired_ratio_geomean"]), ("ci_low", e["exploratory_95_log_t_interval"][0]),
                             ("ci_high", e["exploratory_95_log_t_interval"][1]), ("slower_pairs", e["candidate_slower_pairs"])]:
            close(float(row[field]), value)
    working = {name: sha(ROOT / name) == digest for name, digest in manifest["sources"]["candidate"].items()}
    assert all(working.values()) and working == summary["working_source_matches_candidate"]
    assert len(results["runs"]) == summary["processes"] == 36
    assert samples == summary["samples"] == 6588
    output = {"audit": "Independent final cohort raw/provenance/statistical verification; no builds or timings.",
              "script_sha256": sha(Path(__file__)), "summary_sha256": sha(STUDY / "summary.json"),
              "verified": {"processes": 36, "samples": samples, "paired_vectors": cells_verified, "csv_rows": len(csv_rows),
                           "archived_executables": 4, "source_files_per_core": 6}, "working_source_matches_candidate": working}
    with (HERE / "final_review.json").open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps(output["verified"]))


if __name__ == "__main__":
    audit()
