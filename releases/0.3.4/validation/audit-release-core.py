"""Independent read-only release identity/metadata audit; no builds or timings."""
import ast
import hashlib
import json
from pathlib import Path
import re
import tarfile
import tomllib
import zipfile

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "packages/python/pyproject.toml").is_file())
HERE = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    versions = {}
    for rel, field in (("packages/python/pyproject.toml", "project"),
                       ("packages/rust/Cargo.toml", "package"),
                       ("usecases/computational-biology/code/pyproject.toml", "project"),
                       ("usecases/computational-biology/code/rust/Cargo.toml", "package")):
        data = tomllib.loads((ROOT / rel).read_text(encoding="utf-8"))
        versions[rel] = data[field]["version"]
        assert versions[rel] == "0.3.4"
        if "biology" in rel:
            if field == "project":
                assert data["project"]["dependencies"] == ["hashrope>=0.3.4,<0.4"]
            else:
                assert data["dependencies"]["hashrope"]["version"] == "0.3.4"
    for rel in ("packages/python/hashrope/__init__.py", "usecases/computational-biology/code/src/hashrope_bio/__init__.py"):
        tree = ast.parse((ROOT / rel).read_text(encoding="utf-8"))
        value = next(n.value.value for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "__version__" for t in n.targets))
        assert value == "0.3.4"
        versions[rel] = value
    npm = ROOT / "packages/npm"
    for rel in ("package.json", "package-lock.json"):
        data = read(npm / rel)
        assert data["version"] == "0.3.4"
        versions["packages/npm/" + rel] = data["version"]
        if rel == "package-lock.json":
            assert data["packages"][""]["version"] == "0.3.4"
    assert "export const VERSION = '0.3.4';" in (npm / "src/index.ts").read_text(encoding="utf-8")
    versions["packages/npm/src/index.ts"] = "0.3.4"

    record = read(HERE / "rust-final/results.json")
    assert record["success"] and record["source_unchanged"] and len(record["commands"]) == 7
    tests = 0
    for command in record["commands"]:
        assert command["returncode"] == 0
        log = (HERE.parent / command["log"]).read_bytes()
        assert sha(log) == command["sha256"]
        counts = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;", log.decode("utf-8"))
        assert not any(int(f) or int(i) for _, f, i in counts)
        total = sum(int(p) for p, _, _ in counts)
        assert total == command["passed_tests"]
        tests += total
    assert tests == 221
    archive = HERE / "rust-final/sources.zip"
    assert sha(archive.read_bytes()) == record["source_archive_sha256"]
    with zipfile.ZipFile(archive) as z:
        assert sha(z.read("validation_driver.py")) == record["driver_sha256"]
        for rel, expected in record["source_sha256"].items():
            assert sha(z.read(rel)) == sha((ROOT / rel).read_bytes()) == expected
    old_archive = ROOT / "benchmarks/rust_regression_round7/retained_validation/sources.zip"
    changes = []
    with zipfile.ZipFile(old_archive) as z:
        for rel in record["source_sha256"]:
            old, current = z.read(rel), (ROOT / rel).read_bytes()
            if old == current:
                continue
            changes.append(rel)
            old = old.decode("utf-8").replace("\r\n", "\n")
            current = current.decode("utf-8").replace("\r\n", "\n")
            if rel.endswith("Cargo.toml"):
                assert old.replace('version = "0.3.3"', 'version = "0.3.4"', 1) == current
            elif rel.endswith("Cargo.lock"):
                assert old.replace('name = "hashrope"\nversion = "0.3.3"', 'name = "hashrope"\nversion = "0.3.4"') == current
            else:
                raise AssertionError("Unexpected accepted-source change: " + rel)
    assert sorted(changes) == ["packages/rust/Cargo.lock", "packages/rust/Cargo.toml"]
    rope_sha = sha((ROOT / "packages/rust/src/rope.rs").read_bytes())
    assert rope_sha == "50d9441b03d5e6f4acba329bcd449fe5dbc4354b6f09c2d9c28fb54cc916d4e7"
    assert "fn join_unbalanced(&mut self, left: NodeId, right: NodeId)" in (ROOT / "packages/rust/src/rope.rs").read_text(encoding="utf-8")
    first = read(HERE / "rust/results.json")
    assert not first["success"] and not first["source_unchanged"]

    preserve = read(npm / "evidence/release-0.3.4/runtime-preservation.json")
    baseline_tar = npm / "artifacts/hashrope-0.3.3.tgz"
    current_tar = npm / "artifacts/0.3.4/hashrope-0.3.4.tgz"
    inspection = read(npm / "evidence/release-0.3.4/tarball-inspection.json")
    assert sha(baseline_tar.read_bytes()) == preserve["baselineArtifactSha256"]
    assert sha(current_tar.read_bytes()) == inspection["sha256"]
    identical, metadata = 0, 0
    with tarfile.open(baseline_tar) as before, tarfile.open(current_tar) as after:
        for entry in preserve["files"]:
            old = before.extractfile("package/" + entry["path"]).read()
            new = after.extractfile("package/" + entry["path"]).read()
            assert sha(old) == entry["baselineSha256"] and sha(new) == entry["currentSha256"]
            if old == new:
                identical += 1
            else:
                assert old.replace(b"0.3.3", b"0.3.4") == new
                metadata += 1
    assert (identical, metadata) == (36, 7)
    output = {"status": "passed", "versions": versions, "rust_validation_checks": 7, "rust_tests": tests,
              "rust_current_source_and_fixture_files": len(record["source_sha256"]), "changes_from_accepted_source": changes,
              "retained_rope_sha256": rope_sha, "rejected_direction_trial_absent": True,
              "initial_lock_update_guard_failure_preserved": True,
              "npm_runtime_payloads_identical": identical, "npm_version_only_payloads": metadata,
              "release_notes_sha256": sha((HERE.parent / "README.md").read_bytes()),
              "scope": "Core Python/Rust/npm and biology metadata plus Rust accepted-source/gate identity; Java audited separately by its release owner. No publication, builds or timings."}
    with (HERE / "core-release-review.json").open("x", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
        f.write("\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
