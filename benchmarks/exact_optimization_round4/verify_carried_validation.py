"""Verify the source identity supporting carried-forward round-three tests."""
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "exact_optimization_round3"
ROOT = HERE.parents[1]


def main():
    old = json.loads((PREVIOUS / "optimized_manifest.json").read_text())
    current = json.loads((HERE / "optimized_manifest.json").read_text())
    checks = {}
    with zipfile.ZipFile(PREVIOUS / "optimized.zip") as previous_zip, zipfile.ZipFile(HERE / "optimized.zip") as selected_zip:
        for group in ["core_sha256", "test_sha256"]:
            assert old[group] == current[group]
            for name, digest in current[group].items():
                archive_name = "hashrope/" + name
                before = previous_zip.read(archive_name)
                after = selected_zip.read(archive_name)
                working = (ROOT / name).read_bytes()
                assert before == after == working, name
                assert hashlib.sha256(working).hexdigest() == digest, name
                checks[name] = {"sha256": digest, "group": group, "byte_identical": True}
    assert len(checks) == 20
    document = {
        "validation_status": "658 selected tests carried forward; unchanged suites not rerun in round four",
        "previous_report": "../exact_optimization_round3/VALIDATION.md",
        "previous_archive_sha256": hashlib.sha256((PREVIOUS / "optimized.zip").read_bytes()).hexdigest(),
        "selected_archive_sha256": hashlib.sha256((HERE / "optimized.zip").read_bytes()).hexdigest(),
        "all_source_and_test_files_identical": True,
        "files": checks,
    }
    (HERE / "carried_validation.json").write_text(json.dumps(document, indent=2), encoding="utf-8")
    print(f"{len(checks)} tested source/integration files are byte-identical to round three and the working tree")


if __name__ == "__main__":
    main()
