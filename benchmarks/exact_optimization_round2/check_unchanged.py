"""Verify fingerprint modules and construction/edit code against the snapshot."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
BASE = HERE / ".cache/baseline/hashrope"
OPT = HERE / ".cache/optimized/hashrope"


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value,str) else value).hexdigest()


def rust_function(source,name):
    match = re.search(r"^    (?:pub )?fn " + re.escape(name) + r"\(",source,re.M)
    assert match,name
    start = source.index("{",match.start())
    depth = 1
    end = start+1
    while depth:
        depth += (source[end]=="{") - (source[end]=="}")
        end += 1
    return source[match.start():end]


def main():
    modules = ["packages/python/hashrope/polynomial_hash.py","packages/python/hashrope/sliding.py",
               "packages/python/hashrope/__init__.py","packages/python/hashrope/rope.py","packages/rust/src/polynomial_hash.rs",
               "packages/rust/src/sliding.rs","packages/rust/src/lib.rs"]
    result = {"module_checksums":{},"python_ast":{},"rust_functions":{}}
    for name in modules:
        before,after = digest((BASE/name).read_bytes()),digest((OPT/name).read_bytes())
        result["module_checksums"][name] = {"before":before,"after":after,"unchanged":before==after}
        assert before==after,name
    python = [ast.parse((root/"packages/python/hashrope/rope.py").read_text()) for root in [BASE,OPT]]
    for node in python[0].body:
        if isinstance(node,(ast.FunctionDef,ast.ClassDef)) :
            current=next(n for n in python[1].body if isinstance(n,type(node)) and n.name==node.name)
            before,after=digest(ast.dump(node)),digest(ast.dump(current))
            result["python_ast"][node.name]={"before":before,"after":after,"unchanged":before==after}
            assert before==after,node.name
    rust=[(root/"packages/rust/src/rope.rs").read_text() for root in [BASE,OPT]]
    names=["new","new_lazy","with_hash","with_hasher","hasher","hasher_mut","is_lazy","node_count",
           "node","node_len","node_weight","ensure_hash","make_leaf","make_internal","make_repeat_node",
           "make_repeat","len","hash","weight","height","node_height","balance","join","concat","split",
           "repeat","substr_hash","from_bytes","to_bytes","validate","hash_bytes"]
    for name in names:
        before,after=[digest(rust_function(s,name)) for s in rust]
        result["rust_functions"][name]={"before":before,"after":after,"unchanged":before==after}
        assert before==after,name
    (HERE/"unchanged_code.json").write_text(json.dumps(result,indent=2))
    print(f"Unchanged: {len(modules)} modules, {len(result['python_ast'])} Python definitions, {len(names)} Rust construction/hash/edit functions")


if __name__=="__main__":
    main()
