"""Run the unchanged verifier after both normal pip index probes see 0.3.4."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
probe = json.loads((HERE / "python-resolver-propagation.json").read_text(encoding="utf-8"))
assert len(probe) >= 2
assert {item["package"] for item in probe[-2:]} == {"hashrope", "hashrope-bio"}
assert all(item["contains_0_3_4"] and item["returncode"] == 0 for item in probe[-2:])
source = HERE / "verify-python-registry.py"
original = json.loads((HERE / "python-public-initial-failure.json").read_text(encoding="utf-8"))
assert hashlib.sha256(source.read_bytes()).hexdigest() == original["script_sha256"]
assert not (HERE / "python-public-install.log").exists()
spec = importlib.util.spec_from_file_location("public_verification", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.ENV = module.ROOT / ".release/0.3.4/public-python-confirmed"
module.main()
