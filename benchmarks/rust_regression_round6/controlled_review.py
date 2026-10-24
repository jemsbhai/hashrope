"""Independent, untimed audit of the frozen controlled diagnostic's raw data."""
import hashlib
import json
import math
from pathlib import Path
import statistics as stats
import zipfile

HERE = Path(__file__).resolve().parent
OUT = HERE / "controlled_results"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12), (a, b)


def audit():
    frozen = read(HERE / "controlled_freeze.json")
    built = read(HERE / "controlled_build.json")
    data = read(OUT / "results.json")
    summary = read(OUT / "summary.json")
    assert summary["results_sha256"] == sha(OUT / "results.json")
    assert data["freeze_sha256"] == built["freeze_sha256"] == sha(HERE / "controlled_freeze.json")
    assert data["build_sha256"] == sha(HERE / "controlled_build.json")
    assert data["context_sha256"] == sha(OUT / "context.json")
    assert read(OUT / "context.json") == {k: v for k, v in data.items() if k not in {"runs", "context_sha256"}}
    for name, digest in frozen["files"].items():
        assert sha(HERE / name) == digest
    for name, digest in frozen["parent_archives"].items():
        assert sha(HERE / name) == digest
    assert sha(HERE / "controlled_sources.zip") == frozen["source_archive_sha256"]
    assert sha(HERE / "controlled_binaries.zip") == built["binary_archive_sha256"]
    with zipfile.ZipFile(HERE / "controlled_sources.zip") as archive:
        for label, record in frozen["builds"].items():
            for name, digest in record["input_sha256"].items():
                assert hashlib.sha256(archive.read(f"builds/{label}/{name}")).hexdigest() == digest
                assert sha(HERE / ".cache/controlled" / label / name) == digest
            for slot, arm in record["mapping"].items():
                for name, digest in frozen["source_sha256"][arm].items():
                    assert record["input_sha256"][f"slot_{slot}/{name}"] == digest
    with zipfile.ZipFile(HERE / "controlled_binaries.zip") as archive:
        assert len(archive.namelist()) == 3
        for label, record in built["builds"].items():
            assert hashlib.sha256(archive.read(f"{label}/hashrope-controlled-regression.exe")).hexdigest() == record["binary_sha256"]
            assert sha(Path(record["binary"])) == record["binary_sha256"]
    assert data["schedule"] == frozen["schedule"] and len(data["runs"]) == 36
    counts = {"processes": 0, "blocks": 0, "samples": 0, "pinned_blocks": 0, "estimate_vectors": 0}
    calculated = {}
    for index, run in enumerate(data["runs"]):
        assert {k: run[k] for k in frozen["schedule"][index]} == frozen["schedule"][index]
        assert run["index"] == index and run["exit_status"] == 0
        raw = OUT / run["raw_file"]
        assert sha(raw) == run["raw_sha256"]
        sidecar = read(raw.with_suffix(".json"))
        assert sidecar == {k: run[k] for k in sidecar}
        assert sha(OUT / run["stderr_file"]) == run["stderr_sha256"]
        assert (OUT / run["stderr_file"]).read_bytes() == b""
        mask = data["pinned_mask"] if run["affinity_mode"] == "pinned" else 0
        assert run["command"] == [built["builds"][run["build"]]["binary"], str(run["seed"]), str(mask)]
        lines = iter(raw.read_text().splitlines())
        tag, seed, churn_seed, requested, initial, actual = next(lines).split("\t")
        assert tag == "meta" and int(seed) == run["seed"] and int(churn_seed) == run["churn_seed"]
        assert int(requested) == mask and int(actual) == (mask or int(initial))
        assert mask & int(initial) == mask
        assert run["metadata"] == dict(seed=int(seed), churn_seed=int(churn_seed), requested_affinity=mask,
                                        initial_affinity=int(initial), actual_affinity=int(actual))
        for order, case in enumerate(run["cases"]):
            assert next(lines).split("\t") == ["calibration", case["case"], str(order), str(case["sessions_per_block"]), case["pattern"]]
            for position, block in enumerate(case["blocks"]):
                assert block["slot"] == case["pattern"][position].lower()
                durations, session_medians = [], []
                for session_index, session in enumerate(block["sessions"]):
                    sample_durations = []
                    for sample_index, duration in enumerate(session["durations_ns"]):
                        fields = next(lines).split("\t")
                        expected = ["sample", case["case"], str(order), str(position), block["slot"],
                                    str(session_index), str(sample_index), str(session["operations_per_sample"])]
                        assert fields[:8] == expected
                        assert list(map(int, fields[8:])) == [duration, session["answer"], session["nodes_start"],
                                                              session["nodes_end"], session["root_address"], session["leaf_address"]]
                        assert duration > 0
                        sample_durations.append(duration)
                        counts["samples"] += 1
                    value = stats.median(sample_durations) / session["operations_per_sample"]
                    close(value, session["statistic_ns"])
                    session_medians.append(value)
                    durations.extend(sample_durations)
                fields = next(lines).split("\t")
                assert fields[:5] == ["block", case["case"], str(order), str(position), block["slot"]]
                assert int(fields[5]) == len(block["sessions"]) == case["sessions_per_block"]
                assert int(fields[6]) == block["operations"]
                close(float(fields[7]), stats.median(session_medians))
                close(block["statistic_ns"], stats.median(session_medians))
                assert int(fields[8]) == sum(durations) == block["timed_ns"]
                assert int(fields[9]) == block["wall_ns"] >= block["timed_ns"]
                assert list(map(int, fields[10:])) == [block["processor_before"], block["processor_after"]]
                if mask:
                    assert block["processor_before"] == block["processor_after"] == mask.bit_length() - 1
                    counts["pinned_blocks"] += 1
                metrics = {"statistic_ns": stats.median(session_medians), "timed_mean_ns": sum(durations) / block["operations"],
                           "wall_per_operation_ns": block["wall_ns"] / block["operations"]}
                for metric, value in metrics.items():
                    close(value, block[metric])
                    calculated.setdefault((run["build"], run["affinity_mode"], case["case"], run["seed"], metric, block["slot"]), []).append(value)
                counts["blocks"] += 1
        assert next(lines, None) is None
        counts["processes"] += 1
    highlights = []
    for cell in summary["summary"]:
        numerator, denominator = ("a", "b") if cell["build"] == "current_released" else ("b", "a")
        for metric, estimate in cell["estimates"].items():
            ratios, a_values, b_values = [], [], []
            for seed, pair in enumerate(estimate["process_pairs"], 1):
                prefix = (cell["build"], cell["affinity_mode"], cell["case"], seed, metric)
                a = stats.median(calculated[prefix + (denominator,)])
                b = stats.median(calculated[prefix + (numerator,)])
                assert pair["seed"] == seed
                close(a, pair["denominator_ns"])
                close(b, pair["numerator_ns"])
                close(b / a, pair["ratio"])
                a_values.append(a)
                b_values.append(b)
                ratios.append(b / a)
            logs = [math.log(x) for x in ratios]
            center = stats.mean(logs)
            margin = 2.570581835636305 * stats.stdev(logs) / math.sqrt(6)
            close(math.exp(center), estimate["geomean_ratio"])
            close(min(ratios), estimate["ratio_min"])
            close(max(ratios), estimate["ratio_max"])
            close(stats.median(a_values), estimate["denominator_median_ns"])
            close(stats.median(b_values), estimate["numerator_median_ns"])
            assert sum(x > 1 for x in ratios) == estimate["numerator_slower_count"]
            for expected, actual in zip([math.exp(center - margin), math.exp(center + margin)], estimate["exploratory_95pct_log_t_interval"]):
                close(expected, actual)
            counts["estimate_vectors"] += 1
        e = cell["estimates"]["statistic_ns"]
        highlights.append({k: cell[k] for k in ["build", "affinity_mode", "case"]} |
                          {k: e[k] for k in ["geomean_ratio", "exploratory_95pct_log_t_interval", "numerator_slower_count"]})
    assert counts == dict(processes=36, blocks=720, samples=248764, pinned_blocks=360, estimate_vectors=90)
    assert summary["processes"] == counts["processes"] and summary["blocks"] == counts["blocks"] and summary["samples"] == counts["samples"]
    output = {"audit": "Independent scalar recomputation and frozen artifact/source/raw identity checks; no timings or builds.",
              "script_sha256": sha(Path(__file__)), "results_sha256": sha(OUT / "results.json"),
              "summary_sha256": sha(OUT / "summary.json"), "verified": counts, "primary_estimates": highlights}
    with (HERE / "controlled_review.json").open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    audit()
