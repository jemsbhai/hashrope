"""Record the final candidate's correctness/tooling checks before measurement."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CRATE = ROOT / 'packages/rust'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    destination = HERE / 'validation'
    destination.mkdir()
    commands = [
        ['cargo', 'fmt', '--all', '--check'],
        ['cargo', 'clippy', '--all-features', '--all-targets', '--', '-D', 'warnings'],
        ['cargo', 'test', '--all-features', '--release'],
        ['cargo', 'test', '--doc'],
        ['cargo', 'check', '--lib', '--no-default-features', '--target', 'wasm32-unknown-unknown'],
        ['cargo', 'check', '--lib', '--all-features', '--target', 'wasm32-unknown-unknown'],
        ['cargo', 'doc', '--no-deps', '--all-features'],
    ]
    paths = [CRATE / 'Cargo.toml', CRATE / 'Cargo.lock']
    paths += sorted((CRATE / 'src').rglob('*.rs'))
    paths += sorted((CRATE / 'tests').rglob('*'))
    source_hashes = {p.relative_to(ROOT).as_posix(): sha(p) for p in paths if p.is_file()}
    with zipfile.ZipFile(destination / 'sources.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for name in source_hashes:
            archive.write(ROOT / name, name)
        archive.write(Path(__file__), 'validation_driver.py')
    record = {
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'driver_sha256': sha(Path(__file__)),
        'source_sha256': source_hashes,
        'source_archive_sha256': sha(destination / 'sources.zip'),
        'commands': [],
    }
    environment = dict(os.environ, RUSTDOCFLAGS='-D warnings', CARGO_TERM_COLOR='never')
    for index, arguments in enumerate(commands):
        print('Running ' + ' '.join(arguments), flush=True)
        path = destination / f'{index:02d}.log'
        with path.open('xb') as output:
            result = subprocess.run(arguments, cwd=CRATE, env=environment,
                                    stdout=output, stderr=subprocess.STDOUT)
        text = path.read_text(encoding='utf-8', errors='replace')
        passed = sum(map(int, re.findall(r'test result: ok\. (\d+) passed', text)))
        item = {'arguments': arguments, 'working_directory': str(CRATE),
                'returncode': result.returncode, 'log': path.relative_to(HERE).as_posix(),
                'sha256': sha(path), 'passed_tests': passed}
        record['commands'].append(item)
        print(f'Exit {result.returncode}; passed tests {passed}; log {path.name}', flush=True)
        if result.returncode:
            break
    record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    record['source_unchanged'] = all(sha(ROOT / name) == digest for name, digest in source_hashes.items())
    record['success'] = (len(record['commands']) == len(commands)
                         and all(row['returncode'] == 0 for row in record['commands'])
                         and record['source_unchanged'])
    with (destination / 'results.json').open('x', encoding='utf-8', newline='\n') as output:
        json.dump(record, output, indent=2)
        output.write('\n')
    sys.exit(0 if record['success'] else 1)


if __name__ == '__main__':
    main()
