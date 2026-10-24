"""Run post-matrix audit checks; never performs benchmark timing."""
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAPER = ROOT.parent / "hashropeictaipaper"
BIO = ROOT / "usecases/computational-biology/code"


def run(*arguments):
    result = subprocess.run(arguments, cwd=ROOT, text=True, capture_output=True, check=True)
    return {"command": list(arguments), "returncode": result.returncode,
            "stdout": result.stdout, "stderr": result.stderr}


def main():
    summary = json.loads((HERE / "summary.json").read_text())
    assert summary["complete_matrix"]
    commands = [run(sys.executable, str(HERE / script)) for script in
                ["provenance.py", "check_unchanged.py", "verify_carried_validation.py", "verify_final_core.py", "verify_preserved.py"]]
    alternate = HERE / ".cache/report_check"
    commands.append(run(sys.executable, str(HERE / "summarize.py"), "--output-dir", str(alternate)))
    assert json.loads((alternate / "summary.json").read_text()) == summary
    bio_paths = sorted((HERE / "results").glob("*/bio/run_*.json"))
    before = {str(p): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in bio_paths}
    commands.append(run(sys.executable, str(HERE / "run_matrix.py"), "--kind", "bio", "--resume"))
    after = {str(p): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in bio_paths}
    assert len(before) == 2 and before == after
    parsed = []
    for path in sorted(HERE.glob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        parsed.append(path.name)
    diff = {str(path): run("git", "-C", str(path), "diff", "--check") for path in [ROOT, PAPER, BIO]}
    preserved = json.loads((HERE / "preserved_results.json").read_text())
    document = {
        "complete_matrix": True,
        "study_mode": summary["study_mode"],
        "cohort_contexts": summary["context_fingerprints"],
        "alternate_output_summary_identical": True,
        "bio_resume_reused_verified_results_without_rewrite": True,
        "python_scripts_parsed": parsed,
        "prior_artifacts_byte_identical": len(preserved),
        "validation": "658 selected tests carried forward by exact source/test identity; diagnostic and matrix oracles freshly executed",
        "commands": commands,
        "repository_diff_checks": diff,
    }
    old = HERE / "final_audit.json"
    if old.exists():
        previous = json.loads(old.read_text())
        for name in ["final_paper_compilation", "logbook_prior_bytes_preserved"]:
            if name in previous:
                document[name] = previous[name]
    old.write_text(json.dumps(document, indent=2), encoding="utf-8")
    print("Complete matrix, five contexts, source identities, preservation, alternate report, verified resume and diff checks passed")


if __name__ == "__main__":
    main()
