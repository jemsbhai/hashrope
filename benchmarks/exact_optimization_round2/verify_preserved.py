"""Verify that prior benchmark sources/results remain byte-for-byte preserved."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    manifest=json.loads((HERE/'baseline_manifest.json').read_text())
    checks={}
    for name,record in manifest['files'].items():
        if not name.startswith('hashrope/benchmarks/'):
            continue
        path=ROOT/name.removeprefix('hashrope/')
        current=hashlib.sha256(path.read_bytes()).hexdigest()
        checks[name]={'captured_sha256':record['sha256'],'current_sha256':current,'unchanged':current==record['sha256']}
    assert checks
    assert all(c['unchanged'] for c in checks.values()),[name for name,c in checks.items() if not c['unchanged']]
    (HERE/'preserved_results.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(f'{len(checks)} prior benchmark source/artifact files are byte-identical')


if __name__=='__main__':
    main()
