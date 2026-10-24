"""Run only after publication; verify public crate hashes and consumer contracts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--core-sha256', required=True)
    parser.add_argument('--bio-sha256', required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    run = Path(tempfile.mkdtemp(prefix='hashrope-public-registry-033-'))
    project = run / 'consumer'
    (project / 'src').mkdir(parents=True)
    for relative in ['Cargo.toml', 'src/lib.rs']:
        shutil.copyfile(source / relative, project / relative)
    env = dict(os.environ)
    for key in list(env):
        if key.startswith('CARGO_'):
            del env[key]
    env['CARGO_HOME'] = str(run / 'cargo-home')
    env['CARGO_TARGET_DIR'] = str(run / 'target')
    env['CARGO_REGISTRIES_CRATES_IO_PROTOCOL'] = 'sparse'
    env['CARGO_REGISTRIES_CRATES_IO_INDEX'] = 'sparse+https://index.crates.io/'
    cargo = shutil.which('cargo')
    assert cargo, 'cargo not found'
    evidence = {'directory': str(run), 'cargo': cargo, 'commands': [],
                'consumer_sha256': {name: sha(source / name) for name in ['Cargo.toml', 'src/lib.rs', 'verify.py']}}

    def command(name, arguments):
        result = subprocess.run([cargo, *arguments], cwd=project, env=env,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (run / f'{name}.log').write_text(result.stdout + '\n--- stderr ---\n' + result.stderr, encoding='utf-8')
        evidence['commands'].append({'arguments': arguments, 'returncode': result.returncode,
                                     'log': f'{name}.log', 'sha256': sha(run / f'{name}.log')})
        print(f'{name}: exit {result.returncode}', flush=True)
        if result.returncode:
            raise RuntimeError(f'{name} failed; see {run / (name + ".log")}')
        return result.stdout

    try:
        metadata = json.loads(command('metadata', ['metadata', '--format-version', '1', '--all-features']))
        dependencies = [p for p in metadata['packages'] if p['name'] != 'hashrope-registry-release-smoke']
        assert dependencies
        assert all(p['source'] == 'registry+https://github.com/rust-lang/crates.io-index' for p in dependencies), dependencies
        evidence['registry_packages'] = [{key: p[key] for key in ['name', 'version', 'source', 'manifest_path']} for p in dependencies]
        expected = {'hashrope': args.core_sha256.lower(), 'hashrope-bio': args.bio_sha256.lower()}
        evidence['crate_sha256'] = {}
        for name, digest in expected.items():
            matches = list((run / 'cargo-home' / 'registry' / 'cache').glob(f'*/{name}-0.3.3.crate'))
            assert len(matches) == 1, matches
            actual = sha(matches[0])
            assert actual == digest, (name, actual, digest)
            evidence['crate_sha256'][name] = actual
            package = [p for p in dependencies if p['name'] == name]
            assert len(package) == 1 and package[0]['version'] == '0.3.3', package
        command('default-tests', ['test', '--locked', '--release', '--no-default-features'])
        command('bigint-tests', ['test', '--locked', '--release', '--all-features'])
        command('dependency-tree', ['tree', '--locked', '--all-features'])
        evidence['passed'] = True
    finally:
        if (project / 'Cargo.lock').exists():
            evidence['cargo_lock_sha256'] = sha(project / 'Cargo.lock')
        (run / 'evidence.json').write_text(json.dumps(evidence, indent=2) + '\n', encoding='utf-8')
        print(f'Evidence retained: {run}', flush=True)


if __name__ == '__main__':
    main()
