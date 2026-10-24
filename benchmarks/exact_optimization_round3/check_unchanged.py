"""Audit the selected edit/hash changes and verify all other production bodies."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
BASE = HERE / ".cache/baseline/hashrope"
OPT = HERE / ".cache/optimized/hashrope"

# This explicit list is reviewed before the optimized snapshot is frozen.
EXPECTED_CHANGES = {
    "rope.rs": {"ensure_hash_inner", "split_inner", "hash_range"},
    "polynomial_hash.rs": {"hash"},
}
EXPECTED_ADDITIONS = {"rope.rs": set(), "polynomial_hash.rs": {"hash_default_profile"}}


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def block_end(source, start):
    depth = 1
    end = start + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return end


def rust_functions(source):
    """Extract production functions; occurrence IDs distinguish separate impls."""
    source = source.split("#[cfg(test)]", 1)[0]
    functions = {}
    counts = {}
    pattern = r"^[ \t]*(?:pub(?:\([^)]*\))? )?(?:const )?fn ([a-zA-Z_][a-zA-Z_0-9]*)(?:<[^\n]*?>)?\("
    for match in re.finditer(pattern, source, re.M):
        name = match.group(1)
        occurrence = counts.get(name, 0)
        counts[name] = occurrence + 1
        start = source.index("{", match.end())
        functions[f"{name}#{occurrence}"] = {
            "name": name,
            "signature": source[match.start():start].strip(),
            "source": source[match.start():block_end(source, start)],
        }
    return functions


def declaration(source, name):
    match = re.search(r"^(?:pub )?(?:struct|enum) " + re.escape(name) + r"\b", source, re.M)
    assert match, name
    start = source.index("{", match.start())
    return source[match.start():block_end(source, start)]


def main():
    modules = ["packages/python/hashrope/polynomial_hash.py", "packages/python/hashrope/sliding.py",
               "packages/python/hashrope/__init__.py", "packages/python/hashrope/rope.py",
               "packages/rust/src/sliding.rs", "packages/rust/src/lib.rs"]
    result = {"module_checksums": {}, "python_ast": {}, "rust_functions": {},
              "rust_declarations": {}, "intentional_changes": [], "new_functions": [],
              "note": "Source checks establish change scope, not semantic fingerprint equivalence; baseline differential tests supply that evidence."}
    for name in modules:
        before, after = digest((BASE/name).read_bytes()), digest((OPT/name).read_bytes())
        result["module_checksums"][name] = {"before": before, "after": after, "unchanged": before == after}
        assert before == after, name
    python = [ast.parse((root/"packages/python/hashrope/rope.py").read_text()) for root in [BASE, OPT]]
    for node in python[0].body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            current = next(n for n in python[1].body if isinstance(n, type(node)) and n.name == node.name)
            before, after = digest(ast.dump(node)), digest(ast.dump(current))
            result["python_ast"][node.name] = {"before": before, "after": after, "unchanged": before == after}
            assert before == after, node.name
    for filename in ["rope.rs", "polynomial_hash.rs"]:
        texts = [(root/"packages/rust/src"/filename).read_text() for root in [BASE, OPT]]
        before, after = map(rust_functions, texts)
        assert set(before) <= set(after), (filename, "removed production function", set(before)-set(after))
        for key, original in before.items():
            current = after[key]
            assert original["signature"] == current["signature"], (filename, key, "signature changed")
            old_hash, new_hash = digest(original["source"]), digest(current["source"])
            same = old_hash == new_hash
            allowed = original["name"] in EXPECTED_CHANGES[filename]
            result["rust_functions"][filename+"::"+key] = {"before": old_hash, "after": new_hash,
                                                        "unchanged": same, "intentional_change_allowed": allowed}
            assert same or allowed, (filename, key, "unexpected production change")
            if not same:
                result["intentional_changes"].append(filename+"::"+key)
        for key in sorted(set(after)-set(before)):
            assert after[key]["name"] in EXPECTED_ADDITIONS[filename], (filename, key, "unexpected function added")
            result["new_functions"].append({"name": filename+"::"+key, "sha256": digest(after[key]["source"])})
        names = ["NodeInner", "Arena"] if filename == "rope.rs" else ["PolynomialHash"]
        for name in names:
            old_hash, new_hash = [digest(declaration(text, name)) for text in texts]
            assert old_hash == new_hash, (filename, name, "storage declaration changed")
            result["rust_declarations"][filename+"::"+name] = {"before": old_hash, "after": new_hash, "unchanged": True}
        constants = [re.findall(r"^(?:pub )?(?:const|type) [^\n]+", text, re.M) for text in texts]
        assert constants[0] == constants[1], (filename, "constant or type-alias declaration changed")
        result["rust_declarations"][filename+"::constants_and_aliases"] = {"before": digest(str(constants[0])), "after": digest(str(constants[1])), "unchanged": True}
    (HERE/"unchanged_code.json").write_text(json.dumps(result, indent=2))
    print(f"Verified {len(modules)} unchanged modules, {len(result['python_ast'])} Python definitions; "
          f"audited {len(result['rust_functions'])} existing Rust production functions, "
          f"{len(result['intentional_changes'])} intentional changes and {len(result['new_functions'])} "
          f"addition{'s' if len(result['new_functions']) != 1 else ''}")


if __name__ == "__main__":
    main()
