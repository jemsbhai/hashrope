#!/usr/bin/env python3
"""Validate pinned downstream source against a local Hashrope crate.

Python 3.10+, Cargo and (for --*-source) Git are required. Public snapshots
are downloaded without credentials. Original consumer repositories are never
edited, and local snapshots come from git archive at the pinned commit, not
the working tree. See ../DOWNSTREAM_TESTING.md for scope and known limitations.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import zipfile

PINS = {
    "cdh-sort": "6d40877ba7675e2a3b41dda5ed524b7644ce2ddb",
    "hashrope-bio": "206d1af783ca4c34626cfe1f9ce9a7bd7d4275cb",
}
BIO_DOC_LIMIT = (
    "The pinned gene_diff.rs doctest passes &mut arena twice to diff_genes "
    "and fails E0499 with both the released 0.3.1 crate and the candidate. "
    "The unchanged-source gate uses --all-targets, which excludes doctests."
)


class AuditError(RuntimeError):
    """An actionable source or compatibility failure."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def package_version(manifest: Path) -> str:
    text = manifest.read_text(encoding="utf-8-sig")
    package = re.search(r"(?ms)^\[package\]\s*\n(.*?)(?=^\[|\Z)", text)
    version = re.search(r'^version\s*=\s*"([^"]+)"', package.group(1), re.M) if package else None
    if not version:
        raise AuditError(f"Cannot read package version from {manifest}")
    return version.group(1)


def candidate_fingerprint(candidate: Path) -> dict[str, str]:
    files = [candidate / "Cargo.toml", *sorted((candidate / "src").rglob("*.rs"))]
    return {str(path.relative_to(candidate)).replace("\\", "/"): sha256(path.read_bytes()) for path in files}


def run(command: list[str], cwd: Path, log: Path) -> None:
    """Stream a command to both the console and an artifact; fail on nonzero exit."""
    print(f"\n[{cwd.name}] {subprocess.list2cmdline(command)}", flush=True)
    with log.open("w", encoding="utf-8", newline="\n") as output:
        output.write(json.dumps({"command": command, "cwd": str(cwd)}) + "\n")
        output.flush()
        with subprocess.Popen(
            command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        ) as process:
            assert process.stdout is not None
            for line in process.stdout:
                print(line, end="", flush=True)
                output.write(line)
            code = process.wait()
    if code:
        raise AuditError(f"Command failed with exit {code}; see {log}")


def wanted(consumer: str, name: str) -> bool:
    path = PurePosixPath(name)
    if consumer == "hashrope-bio":
        return name in {"rust/Cargo.toml", "rust/Cargo.lock"} or name.startswith("rust/src/")
    return (
        name in {"Cargo.toml", "Cargo.lock", "README.md", "LICENSE"}
        or (path.parts and path.parts[0] in {"src", "tests", "benches", "examples"})
    )


def write_member(destination: Path, name: str, content: bytes) -> None:
    """Copy only relative regular files; never trust an archive's target paths."""
    path = PurePosixPath(name)
    if not path.parts or path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise AuditError(f"Unsafe archive path: {name!r}")
    if "\\" in name or ":" in name:
        raise AuditError(f"Nonportable archive path: {name!r}")
    target = destination.joinpath(*path.parts).resolve()
    if not target.is_relative_to(destination.resolve()):
        raise AuditError(f"Archive path escapes destination: {name!r}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)


def snapshot(consumer: str, local_source: Path | None, destination: Path) -> dict[str, Any]:
    revision = PINS[consumer]
    destination.mkdir(parents=True)
    info: dict[str, Any] = {"repository": f"jemsbhai/{consumer}", "commit": revision}
    if local_source is not None:
        # git archive reads the named commit without checking out or altering any files.
        source = local_source.resolve()
        check = subprocess.run(
            ["git", "-C", str(source), "rev-parse", "--verify", f"{revision}^{{commit}}"],
            check=True, capture_output=True, text=True,
        )
        if check.stdout.strip() != revision:
            raise AuditError(f"{source} does not contain the required commit {revision}")
        archive = subprocess.run(
            ["git", "-C", str(source), "archive", "--format=tar", revision],
            check=True, capture_output=True,
        ).stdout
        info.update(source="local git archive (working tree ignored)", archive_sha256=sha256(archive))
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
            for member in tar.getmembers():
                if wanted(consumer, member.name):
                    if member.isdir():
                        continue
                    if not member.isfile():
                        raise AuditError(f"Unsupported archive member: {member.name}")
                    content = tar.extractfile(member)
                    if content is None:
                        raise AuditError(f"Unreadable archive member: {member.name}")
                    write_member(destination, member.name, content.read())
    else:
        url = f"https://codeload.github.com/jemsbhai/{consumer}/zip/{revision}"
        print(f"Downloading pinned public snapshot: {url}", flush=True)
        try:
            with urlopen(Request(url, headers={"User-Agent": "hashrope-downstream-verifier"}), timeout=90) as response:
                archive = response.read()
        except (HTTPError, URLError, TimeoutError) as error:
            hint = " Use --cdh-source with an authorized local checkout, or --consumer hashrope-bio for the public gate." if consumer == "cdh-sort" else ""
            raise AuditError(f"Cannot retrieve public {consumer} commit {revision}: {error}.{hint}") from error
        info.update(source=url, archive_sha256=sha256(archive))
        expected_prefix = f"{consumer}-{revision}/"
        with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
            for member in bundle.infolist():
                if not member.filename.startswith(expected_prefix):
                    raise AuditError(f"Unexpected archive root: {member.filename!r}")
                relative = member.filename[len(expected_prefix):]
                if member.is_dir() or not wanted(consumer, relative):
                    continue
                if stat.S_ISLNK(member.external_attr >> 16):
                    raise AuditError(f"Archive contains a symlink: {member.filename}")
                write_member(destination, relative, bundle.read(member))
    manifest = destination / ("rust/Cargo.toml" if consumer == "hashrope-bio" else "Cargo.toml")
    if not manifest.is_file():
        raise AuditError(f"Pinned {consumer} snapshot did not contain {manifest.name}")
    info["files_sha256"] = {
        str(path.relative_to(destination)).replace("\\", "/"): sha256(path.read_bytes())
        for path in sorted(destination.rglob("*")) if path.is_file()
    }
    return info


def cargo_command(cargo: str, *args: str, offline: bool = False) -> list[str]:
    result = [cargo, *args]
    if offline:
        result.append("--offline")
    return result


def verify_cdh(source: Path, candidate: Path, version: str, output: Path, cargo: str, offline: bool) -> None:
    # Only the temporary checkout's configuration and lockfile are changed.
    config = source / ".cargo" / "config.toml"
    config.parent.mkdir()
    config.write_text(
        "[patch.crates-io]\nhashrope = { path = " + json.dumps(candidate.as_posix()) + " }\n",
        encoding="utf-8",
    )
    run(cargo_command(cargo, "update", "-p", "hashrope", "--precise", version, offline=offline), source, output / "cdh-update.log")
    run(cargo_command(cargo, "test", "--lib", "--tests", offline=offline), source, output / "cdh-tests.log")
    run(cargo_command(cargo, "check", "--all-targets", offline=offline), source, output / "cdh-all-targets.log")
    command = cargo_command(cargo, "test", "--bench", "diag_benchmark", offline=offline)
    run([*command, "--", "build_rope_lazy", "--test"], source, output / "cdh-lazy-smoke.log")
    # Validate resolution explicitly; a green test of a registry crate is not this gate.
    metadata = subprocess.run(
        cargo_command(cargo, "metadata", "--format-version", "1", offline=offline),
        cwd=source, check=True, capture_output=True, text=True,
    )
    packages = json.loads(metadata.stdout)["packages"]
    resolved = [p for p in packages if p["name"] == "hashrope"]
    if len(resolved) != 1 or Path(resolved[0]["manifest_path"]).resolve() != candidate / "Cargo.toml":
        raise AuditError("cdh-sort did not resolve the local candidate Hashrope crate")
    (output / "cdh-resolved-candidate.json").write_text(json.dumps(resolved[0], indent=2) + "\n", encoding="utf-8")
    shutil.copyfile(source / "Cargo.lock", output / "cdh-Cargo.lock")


def verify_bio(source: Path, candidate: Path, output: Path, cargo: str, offline: bool) -> None:
    original = source / "rust" / "Cargo.toml"
    shutil.copyfile(original, output / "bio-original-Cargo.toml")
    # This is a new test harness, not an edit to the consumer's dependency contract.
    # Its source paths point directly to the unchanged archived library and binaries.
    harness = source / "compatibility-harness"
    harness.mkdir()
    manifest = [
        "[package]", 'name = "hashrope-bio-compatibility-audit"',
        'version = "0.0.0"', 'edition = "2021"', "autobins = false", "",
        "[lib]", 'name = "hashrope_bio"', 'path = "../rust/src/lib.rs"', "",
        "[dependencies]", "hashrope = { path = " + json.dumps(candidate.as_posix()) + " }", "",
    ]
    for binary in sorted((source / "rust" / "src" / "bin").glob("*.rs")):
        manifest.extend(["[[bin]]", f'name = "{binary.stem}"', f'path = "../rust/src/bin/{binary.name}"', ""])
    (harness / "Cargo.toml").write_text("\n".join(manifest), encoding="utf-8")
    shutil.copyfile(harness / "Cargo.toml", output / "bio-harness-Cargo.toml")
    run(cargo_command(cargo, "test", "--all-targets", offline=offline), harness, output / "bio-tests.log")
    shutil.copyfile(harness / "Cargo.lock", output / "bio-Cargo.lock")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consumer", choices=["all", *PINS], default="all")
    parser.add_argument("--candidate", type=Path, default=Path(__file__).resolve().parents[1], help="Hashrope crate directory")
    parser.add_argument("--cdh-source", type=Path, help="Authorized cdh-sort Git checkout containing the pinned commit; ignores dirty files")
    parser.add_argument("--bio-source", type=Path, help="Optional local hashrope-bio Git checkout containing the pinned commit")
    parser.add_argument("--output-dir", type=Path, default=Path("downstream-results"), help="Directory for command logs and report.json")
    parser.add_argument("--cargo", default="cargo", help="Cargo executable path")
    parser.add_argument("--offline", action="store_true", help="Use Cargo's offline mode; public source downloads still require network")
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {"status": "failed", "consumers": {}, "limitations": [BIO_DOC_LIMIT]}
    try:
        if not (candidate / "src" / "lib.rs").is_file():
            raise AuditError(f"Not a Hashrope crate directory: {candidate}")
        version = package_version(candidate / "Cargo.toml")
        fingerprint = candidate_fingerprint(candidate)
        cargo_version = subprocess.run([args.cargo, "--version"], check=True, capture_output=True, text=True).stdout.strip()
        report.update(candidate=str(candidate), candidate_version=version, candidate_files_sha256=fingerprint,
                      cargo_version=cargo_version, python_version=sys.version, platform=sys.platform)
        selected = list(PINS) if args.consumer == "all" else [args.consumer]
        with tempfile.TemporaryDirectory(prefix="hashrope-downstream-") as temporary:
            work = Path(temporary).resolve()
            # Cleanup is limited to the fresh temporary directory created above.
            if work.parent != Path(tempfile.gettempdir()).resolve():
                raise AuditError("Unexpected temporary directory parent")
            for consumer in selected:
                local = args.cdh_source if consumer == "cdh-sort" else args.bio_source
                source = work / consumer
                info = snapshot(consumer, local, source)
                report["consumers"][consumer] = info
                if consumer == "cdh-sort":
                    verify_cdh(source, candidate, version, output, args.cargo, args.offline)
                else:
                    info["adoption_limit"] = "Original ^0.2.1 dependency and stale path exclude 0.3.x; only the isolated harness selects the candidate."
                    info["doctest_limit"] = BIO_DOC_LIMIT
                    verify_bio(source, candidate, output, args.cargo, args.offline)
                info["status"] = "passed"
        if candidate_fingerprint(candidate) != fingerprint:
            raise AuditError("Candidate sources changed during verification; rerun against a stable checkout")
        report["status"] = "passed"
        print("\nAll selected downstream compatibility gates passed.", flush=True)
        return 0
    except (AuditError, OSError, subprocess.CalledProcessError, zipfile.BadZipFile, tarfile.TarError) as error:
        report["error"] = str(error)
        print(f"\nDownstream verification failed: {error}", file=sys.stderr, flush=True)
        return 1
    finally:
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"Evidence: {output / 'report.json'}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
