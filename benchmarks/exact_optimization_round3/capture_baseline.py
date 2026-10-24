"""Snapshot the completed round-two exact implementation before round-three edits."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAPER = ROOT.parent / "hashropeictaipaper"
BIO = ROOT / "usecases/computational-biology/code"


def main():
    destination = HERE / "baseline.zip"
    if destination.exists():
        raise SystemExit("Baseline already exists; never replace a captured baseline")
    selection = [
        (ROOT, "hashrope", "packages/python/hashrope"),
        (ROOT, "hashrope", "packages/python/tests"),
        (ROOT, "hashrope", "packages/python/pyproject.toml"),
        (ROOT, "hashrope", "packages/rust/src"),
        (ROOT, "hashrope", "packages/rust/tests"),
        (ROOT, "hashrope", "packages/rust/Cargo.toml"),
        (ROOT, "hashrope", "packages/rust/README.md"),
        (ROOT, "hashrope", "benchmarks/consumer_exactness"),
        (ROOT, "hashrope", "benchmarks/comparison_safety"),
        (ROOT, "hashrope", "benchmarks/exact_optimization"),
        (ROOT, "hashrope", "benchmarks/exact_optimization_round2"),
        (BIO, "bio", "src"), (BIO, "bio", "rust/src"),
        (BIO, "bio", "rust/Cargo.toml"), (BIO, "bio", "rust/Cargo.lock"),
        (BIO, "bio", "pyproject.toml"),
        (PAPER, "paper", "src"), (PAPER, "paper", "scripts"),
        (PAPER, "paper", "tests"), (PAPER, "paper", "third_party"),
    ]
    files = {}
    for base, prefix, relative in selection:
        path = base / relative
        candidates = sorted(path.rglob("*")) if path.is_dir() else [path]
        for file in candidates:
            if not file.is_file() or ("__pycache__" in file.parts or ".cache" in file.parts) or file.suffix in {".npy", ".pyc"}:
                continue
            name = prefix + "/" + file.relative_to(base).as_posix()
            files[name] = file.read_bytes()
    manifest = {
        "captured_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "purpose": "Current completed round-two exact implementation; baseline for round three, not git HEAD",
        "git_heads": {name: subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
                      for name, path in [("hashrope", ROOT), ("paper", PAPER), ("bio", BIO)]},
        "files": {name: {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)} for name, data in files.items()},
        "excluded": "Generated token caches and raw datasets remain at original paths; dataset hashes are in captured benchmark provenance",
    }
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    manifest["archive_sha256"] = hashlib.sha256(destination.read_bytes()).hexdigest()
    (HERE / "baseline_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    with zipfile.ZipFile(destination) as archive:
        archive.extractall(HERE / ".cache/baseline")
    print(json.dumps({"files": len(files), "archive": str(destination), "bytes": destination.stat().st_size,
                      "sha256": manifest["archive_sha256"]}))


if __name__ == "__main__":
    main()
