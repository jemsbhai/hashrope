"""After authorized publication, verify public bytes and fresh installed consumers."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import urllib.error
import urllib.request
import venv

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "packages/python/pyproject.toml").is_file())
HERE = Path(__file__).resolve().parent
VERSION = "0.3.4"
EXPECTED = {x["filename"]: x for x in json.loads((HERE / "python-artifacts.json").read_text(encoding="utf-8"))["artifacts"]}
ENV = ROOT / ".release/0.3.4/public-python"


def main():
    assert not ENV.exists(), "Use the original preserved verification environment; do not overwrite it"
    records, attempts = [], []
    for package in ("hashrope", "hashrope-bio"):
        url = f"https://pypi.org/pypi/{package}/{VERSION}/json"
        for delay in (0, 5, 10, 20):
            time.sleep(delay)
            try:
                with urllib.request.urlopen(url, timeout=30) as response:
                    metadata = json.load(response)
                attempts.append({"url": url, "status": 200})
                break
            except urllib.error.HTTPError as error:
                attempts.append({"url": url, "status": error.code})
                (HERE / "python-public-attempts.json").write_text(json.dumps(attempts, indent=2) + "\n", encoding="utf-8")
                if error.code != 404 or delay == 20:
                    raise
        assert metadata["info"]["version"] == VERSION
        for entry in metadata["urls"]:
            reference = EXPECTED[entry["filename"]]
            content = urllib.request.urlopen(entry["url"], timeout=30).read()
            digest = hashlib.sha256(content).hexdigest()
            assert digest == entry["digests"]["sha256"] == reference["sha256"]
            assert len(content) == reference["bytes"]
            records.append({"package": package, "filename": entry["filename"], "url": entry["url"], "sha256": digest, "bytes": len(content)})
    assert len(records) == 4
    venv.create(ENV, with_pip=True)
    python = ENV / "Scripts/python.exe"
    commands = [
        [python, "-m", "pip", "--isolated", "install", "--no-cache-dir", "--index-url", "https://pypi.org/simple", "hashrope==0.3.4", "hashrope-bio==0.3.4", "pytest"],
        [python, "-m", "pip", "check"],
        [python, "-I", "-c", "import hashrope,hashrope_bio,sys; from pathlib import Path; assert hashrope.__version__==hashrope_bio.__version__=='0.3.4'; assert all(Path(m.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()) for m in (hashrope,hashrope_bio)); print(hashrope.__file__,hashrope_bio.__file__)"],
        [python, "-I", ROOT / "packages/python/scripts/verify_release.py"],
        [python, "-I", "-m", "pytest", "--import-mode=importlib", "--rootdir", ENV, ROOT / "packages/python/tests", ROOT / "usecases/computational-biology/code/tests", "-q"],
    ]
    with (HERE / "python-public-install.log").open("x", encoding="utf-8") as log:
        for command in commands:
            result = subprocess.run(list(map(str, command)), cwd=ENV, stdout=log, stderr=subprocess.STDOUT)
            log.flush()
            assert result.returncode == 0, "Public install check failed; preserved log"
    assert "244 passed" in (HERE / "python-public-install.log").read_text(encoding="utf-8")
    report = {"status": "public artifacts match reviewed builds; fresh registry installation and 244 tests passed", "version": VERSION,
              "artifacts": records, "metadata_attempts": attempts,
              "commands": [list(map(str, command)) for command in commands],
              "log_sha256": hashlib.sha256((HERE / "python-public-install.log").read_bytes()).hexdigest()}
    (HERE / "python-public-install.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(report["status"])


if __name__ == "__main__":
    main()
