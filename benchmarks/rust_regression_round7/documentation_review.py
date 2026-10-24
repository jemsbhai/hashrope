"""Read-only final documentation/validation audit; no build or timing work."""
import hashlib
import json
import math
import re
import statistics
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8"))


def numbers(text):
    return [float(n) for n in re.findall(r"[-+]?\d+(?:\.\d+)?", text.replace("−", "-"))]


def close(actual, expected, tolerance=0.000050001):
    assert abs(actual - expected) <= tolerance, (actual, expected, tolerance)


def main():
    summary_path = HERE / "final_comparison/summary.json"
    summary = load("final_comparison/summary.json")
    table = (HERE / "final_comparison/README.md").read_text(encoding="utf-8")
    assert digest(summary_path.read_bytes()) in table
    cells = {name: row for group in summary["cells"].values() for name, row in group.items()}
    rendered = [line for line in table.splitlines() if line.startswith("| `")]
    assert len(rendered) == len(cells) == 66
    seen = set()
    for line in rendered:
        cols = [s.strip() for s in line.strip("|").split("|")]
        name = cols[0].strip("`")
        assert name not in seen
        seen.add(name)
        row = cells[name]
        factor = 1 if row["unit"] == "ns" else 1000
        unit = "ns" if factor == 1 else "µs"
        vectors = row["process_medians"]
        ratios = [c / b for b, c in zip(vectors["baseline"], vectors["candidate"])]
        assert len(ratios) == 12
        ratio = math.exp(statistics.mean(map(math.log, ratios)))
        for col, arm in zip(cols[1:3], ("baseline", "candidate")):
            assert col.endswith(unit)
            close(numbers(col)[0], statistics.median(vectors[arm]) * factor)
        ratio_ci = numbers(cols[3])
        for actual, expected in zip(ratio_ci, [ratio] + row["exploratory_95_log_t_interval"]):
            close(actual, expected)
        assert numbers(cols[4]) == [sum(r > 1 for r in ratios), 12]
        for actual, expected in zip(numbers(cols[5]), [min(ratios), max(ratios)]):
            close(actual, expected)
        for actual, expected in zip(numbers(cols[6]), row["order_geomeans"].values()):
            close(actual, expected)

    report = (HERE / "README.md").read_text(encoding="utf-8")
    retained_table = report.split("### Remaining Ropey gaps", 1)[1].split("## Reproducibility", 1)[0]
    ropey_rows = [line for line in retained_table.splitlines() if re.match(r"\| [\d,]+ /", line)]
    assert len(ropey_rows) == 6
    ropey_checks = []
    for line in ropey_rows:
        cols = [s.strip() for s in line.strip("|").split("|")]
        size, history = cols[0].replace(",", "").split(" / ")
        key = size + "/" + str(history == "shared").lower() + "/"
        h = cells[key + "hashrope_exact"]["process_medians"]["baseline"]
        r = cells[key + "ropey_exact"]["process_medians"]["baseline"]
        expected = [statistics.median(h) * 1000, statistics.median(r) * 1000,
                    math.exp(statistics.mean(math.log(a / b) for a, b in zip(h, r)))]
        for actual, value in zip([numbers(c)[0] for c in cols[1:]], expected):
            close(actual, value, 0.000500001)
        ropey_checks.append({"cell": key, "baseline_hashrope_ropey_ratio": expected[2]})

    main_map = [
        ("100000/false/hashrope_exact", 1000),
        ("compare/1048576/independent_equal/exact_lcp", .001),
        ("compare/65536/independent_equal/exact_lcp", .001),
        ("compare/1048576/late_mismatch/exact_lcp", .001),
        ("compare/self_concat_DAG_2TiB/exact_lcp", 1),
        ("core/concat_balanced_64KiB", 1),
        ("compare/65536/early_mismatch/exact_lcp", 1),
        ("compare/1048576/early_mismatch/exact_equal", 1),
        ("100000/true/hashrope_exact", 1000),
        ("1000000/false/hashrope_exact", 1000),
        ("1000000/true/hashrope_exact", 1000),
        ("core/split_chunked_64KiB", 1),
    ]
    main_table = report.split("### Final trial results", 1)[1].split("Neither primary", 1)[0]
    main_rows = [line for line in main_table.splitlines() if line.startswith("| ")][1:]
    assert len(main_rows) == len(main_map) == 12
    for line, (name, factor) in zip(main_rows, main_map):
        cols = [c.strip().replace("**", "") for c in line.strip("|").split("|")]
        row = cells[name]
        for col, arm in zip(cols[1:3], ("baseline", "candidate")):
            close(numbers(col)[0], row["median_of_process_medians"][arm] * factor, .000500001)
        change = numbers(cols[3])
        close(change[0], (row["paired_ratio_geomean"] - 1) * 100, .050000001)
        for actual, expected in zip(change[1:], row["exploratory_95_log_t_interval"]):
            close(actual, expected)
        assert numbers(cols[4]) == [row["candidate_slower_pairs"], 12]
    fingerprint_percent = (cells["1000000/false/hashrope_candidate"]["paired_ratio_geomean"] - 1) * 100
    assert f"{fingerprint_percent:.1f}%".replace("-", "−") in report

    validation = {}
    for phase in ("validation", "retained_validation"):
        result = load(phase + "/results.json")
        assert result["success"] and result["source_unchanged"]
        assert len(result["commands"]) == 7
        counts = [0, 0, 0]
        for command in result["commands"]:
            assert command["returncode"] == 0
            raw = (HERE / command["log"]).read_bytes()
            assert digest(raw) == command["sha256"]
            parsed = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;", raw.decode("utf-8"))
            count = [sum(int(row[i]) for row in parsed) for i in range(3)]
            assert count[0] == command["passed_tests"]
            counts = [a + b for a, b in zip(counts, count)]
        assert counts == [221, 0, 0], counts
        archive = HERE / phase / "sources.zip"
        assert digest(archive.read_bytes()) == result["source_archive_sha256"]
        with zipfile.ZipFile(archive) as z:
            assert set(z.namelist()) == set(result["source_sha256"]) | {"validation_driver.py"}
            assert digest(z.read("validation_driver.py")) == result["driver_sha256"]
            for name, expected in result["source_sha256"].items():
                assert digest(z.read(name)) == expected
                if phase == "retained_validation":
                    assert digest((ROOT / name).read_bytes()) == expected
        validation[phase] = {"checks": 7, "passed_tests": counts[0],
                             "archived_source_files": len(result["source_sha256"]),
                             "results_sha256": digest((HERE / phase / "results.json").read_bytes())}

    identity = load("retained_source_identity.json")
    baseline = HERE / "final_comparison/baseline.zip"
    assert digest(baseline.read_bytes()) == identity["baseline_archive_sha256"]
    with zipfile.ZipFile(baseline) as z:
        assert len(identity["matches_frozen_round6"]) == 6
        for name, matches in identity["matches_frozen_round6"].items():
            assert matches and (ROOT / name).read_bytes() == z.read("hashrope/" + name)

    docs = [HERE / "README.md", HERE / "VALIDATION.md", HERE / "final_comparison/README.md",
            ROOT / "README.md", ROOT / "packages/rust/README.md",
            ROOT / "packages/rust/PARITY.md", ROOT / "packages/rust/CHANGELOG.md"]
    output = {"status": "passed", "builds_or_timings_run": False,
              "final_table_cells": 66, "main_trial_rows": 12, "retained_ropey_rows": ropey_checks,
              "validation": validation, "current_core_files_equal_baseline": 6,
              "fingerprint_candidate_time_change_percent": fingerprint_percent,
              "documentation_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): digest(p.read_bytes()) for p in docs}}
    with (HERE / "documentation_review.json").open("x", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
        f.write("\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
