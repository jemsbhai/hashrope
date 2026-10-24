"""Build and validate synchronized Python artifacts without publishing."""
import ast
import configparser
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import venv

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "packages/python/pyproject.toml").is_file())
HERE = Path(__file__).resolve().parent
STAGE = ROOT / ".release/0.3.4"
DIST = STAGE / "python-dist"
ENV = STAGE / "python-installed"
BIO = ROOT / "usecases/computational-biology/code"
VERSION = "0.3.4"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_check():
    records = {}
    for repository, subpath in ((ROOT, "packages/python/hashrope"), (BIO, "src/hashrope_bio")):
        for path in sorted((repository / subpath).rglob("*.py")):
            relative = path.relative_to(repository).as_posix()
            before = subprocess.check_output(["git", "show", "HEAD:" + relative], cwd=repository).decode("utf-8")
            after = path.read_text(encoding="utf-8")
            if path.name == "__init__.py":
                before = before.replace('__version__ = "0.3.3"', '__version__ = "0.3.4"')
            assert ast.dump(ast.parse(before)) == ast.dump(ast.parse(after)), relative
            records[path.relative_to(ROOT).as_posix()] = sha(path)
    return records


def main():
    assert not DIST.exists() and not ENV.exists(), "Use a fresh stage; do not overwrite reviewed artifacts"
    DIST.mkdir(parents=True)
    env = os.environ.copy()
    env["PIP_INDEX_URL"] = "https://pypi.org/simple"
    env.pop("PIP_EXTRA_INDEX_URL", None)
    commands = []

    def run(name, command, cwd=ROOT):
        path = HERE / name
        assert not path.exists()
        print(name, flush=True)
        with path.open("w", encoding="utf-8") as log:
            result = subprocess.run(list(map(str, command)), cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
        commands.append({"arguments": list(map(str, command)), "cwd": str(cwd), "returncode": result.returncode,
                         "log": path.name, "log_sha256": sha(path)})
        if result.returncode:
            (HERE / "python-build-failure.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")
            raise SystemExit(f"Failed: {name}; retained log")

    source = source_check()
    run("python-core-build.log", [sys.executable, "-m", "build", "--outdir", DIST, ROOT / "packages/python"])
    run("python-bio-build.log", [sys.executable, "-m", "build", "--outdir", DIST, BIO])
    artifacts = sorted(DIST.iterdir())
    assert len(artifacts) == 4
    run("python-twine-check.log", [sys.executable, "-m", "twine", "check", "--strict", *artifacts])
    venv.create(ENV, with_pip=True)
    python = ENV / "Scripts/python.exe"
    run("python-wheel-install.log", [python, "-m", "pip", "install", "--no-cache-dir", "--index-url", "https://pypi.org/simple", "pytest", *sorted(DIST.glob("*.whl"))], ENV)
    run("python-pip-check.log", [python, "-m", "pip", "check"], ENV)
    run("python-installed-origin.log", [python, "-I", "-c", "import hashrope,hashrope_bio,sys,importlib.metadata; from pathlib import Path; assert hashrope.__version__==hashrope_bio.__version__=='0.3.4'; assert all(Path(m.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()) for m in (hashrope,hashrope_bio)); print(hashrope.__file__,hashrope_bio.__file__); print(importlib.metadata.version('hashrope'),importlib.metadata.version('hashrope-bio'))"], ENV)
    run("python-final-gate.log", [python, "-I", ROOT / "packages/python/scripts/verify_release.py"], ENV)
    run("python-final-installed-tests.log", [python, "-I", "-m", "pytest", "--import-mode=importlib", "--rootdir", ENV, ROOT / "packages/python/tests", BIO / "tests", "-q"], ENV)
    assert "244 passed" in (HERE / "python-final-installed-tests.log").read_text(encoding="utf-8")
    assert source == source_check()
    report = {"status": "built and verified locally; not published", "version": VERSION,
              "artifacts": [{"filename": p.name, "path": str(p), "bytes": p.stat().st_size, "sha256": sha(p)} for p in artifacts],
              "installed_tests": 244, "exact_collision_profile_repeat_gate": "pass", "metadata_check": "all four pass",
              "scope": "Wheels built from source distributions; fresh isolated site-packages imports. No algorithm change from 0.3.3; __version__ metadata only.",
              "source_sha256": source, "commands": commands}
    (HERE / "python-artifacts.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    config = configparser.ConfigParser()
    config.read(Path.home() / ".pypirc")
    auth = {"twine_username_env_present": bool(os.environ.get("TWINE_USERNAME")),
            "twine_password_env_present": bool(os.environ.get("TWINE_PASSWORD")),
            "pypirc_present": (Path.home() / ".pypirc").is_file(),
            "pypirc_pypi_password_present": bool(config.get("pypi", "password", fallback=""))}
    (HERE / "python-auth-availability.json").write_text(json.dumps(auth, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "artifacts": report["artifacts"], "installed_tests": 244, "auth": auth}, indent=2))


if __name__ == "__main__":
    main()
