"""Verify the retained JSONL fixture with published Pollard/Hashrope wheels.

Put the extracted wheels on PYTHONPATH first. This script does not install
packages, change a store on disk, or contact any external model/provider.
"""

import argparse
import hashlib
from pathlib import Path
import zipfile

import hashrope
import pollard
from pollard import HashRopeStore
from pollard.tree import Node, NodeKind


def verify_wheel(path, module, expected):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    package = module.__name__
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.startswith(package + "/") and name.endswith(".py"):
                loaded = Path(module.__file__).parent / name.removeprefix(package + "/")
                assert loaded.read_bytes() == archive.read(name), f"Not published source: {loaded}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hashrope-wheel", required=True, type=Path)
    parser.add_argument("--pollard-wheel", required=True, type=Path)
    args = parser.parse_args()
    verify_wheel(args.hashrope_wheel, hashrope, "710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2")
    verify_wheel(args.pollard_wheel, pollard, "569fb5f130a82c9be327b8dcbd285e3be063200bd9773ca15c5d6bb62edd627f")
    store = HashRopeStore()
    parent = Node.make(kind=NodeKind.ROOT, parent=None,
                       payload={"run": "parity-golden-2026", "label": "caf\u00e9 \u6f22\u5b57 \U0001f600"})
    child = Node.make(kind=NodeKind.MODEL_CALL, parent=parent.id,
                      payload={"model": "mock-1", "prompt": "line\nbreak"},
                      result={"text": "na\u00efve \U0001f980"})
    store.put(parent)
    store.put(child)
    store.update_meta(child.id, {"label": "kept", "bytes": [0, 127, 255], "unicode": "e\u0301"})
    store.validate_log()
    data = store.to_bytes()
    expected = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "pollard-1.6.0-hashrope-0.2.2.jsonl"
    assert data == expected.read_bytes()
    assert len(data) == 785
    assert hashlib.sha256(data).hexdigest() == "45174ce9b17341c5bd9efb800a403308d148ce799b60862896341c248a37e623"
    replayed = HashRopeStore(data)
    assert replayed.content_hash() == store.content_hash() == hashrope.PolynomialHash().hash(data) == 2269240639182000324
    assert replayed.to_bytes() == data
    print("Published Pollard JSONL matches retained bytes and hash")


if __name__ == "__main__":
    main()
