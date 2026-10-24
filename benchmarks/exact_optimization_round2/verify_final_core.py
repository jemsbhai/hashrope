"""Verify that final documentation edits preserve the measured core code.

The optimized archive remains the authoritative benchmark source. Python
docstrings/comments may change after timing; all other AST nodes must match.
Rust and every other core file must remain byte-identical.
"""
import ast
import hashlib
import json
from pathlib import Path
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class WithoutDocstrings(ast.NodeTransformer):
    def generic_visit(self, node):
        super().generic_visit(node)
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr):
                expression = node.body[0].value
                if isinstance(expression, ast.Constant) and isinstance(expression.value, str):
                    node.body = node.body[1:]
        return node


def code_ast(source):
    return ast.dump(WithoutDocstrings().visit(ast.parse(source)), include_attributes=False)


def main():
    manifest = json.loads((HERE / "optimized_manifest.json").read_text())
    rows = []
    with zipfile.ZipFile(HERE / "optimized.zip") as archive:
        for name, expected in manifest["core_sha256"].items():
            frozen = archive.read("hashrope/" + name)
            final = (ROOT / name).read_bytes()
            frozen_sha = hashlib.sha256(frozen).hexdigest()
            final_sha = hashlib.sha256(final).hexdigest()
            assert frozen_sha == expected, f"Snapshot mismatch: {name}"
            same = frozen == final
            if not same:
                assert name.endswith(".py"), f"Non-Python source changed: {name}"
                assert code_ast(frozen) == code_ast(final), f"Executable AST changed: {name}"
            rows.append({"path": name, "frozen_sha256": frozen_sha,
                         "final_sha256": final_sha, "byte_identical": same,
                         "executable_ast_identical": True})
    result = {"method": __doc__, "all_measured_code_preserved": True, "files": rows}
    (HERE / "final_source_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("Measured core preserved; documentation-only changes:",
          [row["path"] for row in rows if not row["byte_identical"]])


if __name__ == "__main__":
    main()
