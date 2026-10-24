"""Freeze committed 0.3.3 sources before candidate edits; never overwrite evidence."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REVISION = 'ffba5ee'


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args])


def main():
    archive = HERE / 'baseline.zip'
    manifest = HERE / 'baseline_manifest.json'
    if archive.exists() or manifest.exists():
        raise RuntimeError('Baseline evidence already exists; refusing overwrite')
    revision = git('rev-parse', REVISION).decode().strip()
    paths = git('ls-tree', '-r', '--name-only', revision, '--',
                'packages/rust/src', 'packages/rust/Cargo.toml').decode().splitlines()
    assert len(paths) >= 6 and 'packages/rust/Cargo.toml' in paths
    files = {path: git('show', f'{revision}:{path}') for path in paths}
    assert b'version = "0.3.3"' in files['packages/rust/Cargo.toml']
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as output:
        for name, data in files.items():
            output.writestr('hashrope/' + name, data)
    with zipfile.ZipFile(archive) as snapshot:
        snapshot.extractall(HERE / '.cache/baseline')
    evidence = {
        'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_revision': revision,
        'source_tree': git('rev-parse', revision + '^{tree}').decode().strip(),
        'label': 'published Rust 0.3.3 committed implementation',
        'source_method': 'git show for every src file and Cargo.toml; ignores working-tree changes',
        'scope': 'Source snapshot; driver lockfiles will come unchanged from release integration evidence.',
        'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
        'source_sha256': {name: hashlib.sha256(data).hexdigest() for name, data in files.items()},
        'freeze_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    with manifest.open('x', encoding='utf-8', newline='\n') as output:
        output.write(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(evidence, indent=2))


if __name__ == '__main__':
    main()
