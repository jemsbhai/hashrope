"""Repeat the unchanged verifier after preserving PyPI propagation failures."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
source = HERE / "verify-python-registry.py"
failed_log = HERE / "python-public-install.log"
preserved_log = HERE / "python-public-install-initial.log"
assert failed_log.is_file() and not preserved_log.exists()
assert "No matching distribution found for hashrope==0.3.4" in failed_log.read_text(encoding="utf-8")
failed_log.rename(preserved_log)
(HERE / "python-public-index-failure.json").write_text(json.dumps({
    "returncode": 1,
    "failure": "Initial public-index installation could not see hashrope 0.3.4; later canonical index probes list both packages.",
    "log": preserved_log.name,
    "log_sha256": hashlib.sha256(preserved_log.read_bytes()).hexdigest(),
    "verifier_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "retry": "Unchanged verification logic, new empty environment; ordinary isolated PyPI install remains mandatory.",
}, indent=2) + "\n", encoding="utf-8")
spec = importlib.util.spec_from_file_location("public_verification", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.ENV = module.ROOT / ".release/0.3.4/public-python-final"
module.main()
