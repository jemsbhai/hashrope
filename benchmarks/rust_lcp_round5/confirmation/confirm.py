"""Fixed supplemental study using the original archived executables, no builds."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
from run import LABELS, KINDS, now, sha, read, save, verify, parse_micro, parse_churn


def freeze():
    verify()
    for name in ['manifest.json', 'drivers.zip']:
        if (HERE / name).exists():
            raise RuntimeError(f'Refusing to overwrite {name}')
    files = [HERE / name for name in ['PLAN.md', 'confirm.py', 'summarize_confirmation.py']]
    with zipfile.ZipFile(HERE / 'drivers.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, path.name)
    dependencies = ['baseline.zip', 'candidate.zip', 'drivers.zip', 'manifest.json',
                    'build.json', 'binaries.zip', 'results.json', 'summary.json', 'run.py', 'summarize.py']
    manifest = {'declared_utc': now(),
                'parent_artifacts': {name: sha(PARENT / name) for name in dependencies},
                'drivers': {path.name: sha(path) for path in files},
                'drivers_archive_sha256': sha(HERE / 'drivers.zip'),
                'protocol': {'pairs_per_cohort': 6, 'global_pairs': list(range(3, 9)),
                             'order': [list(reversed(LABELS)) if i % 2 == 0 else list(LABELS) for i in range(6)],
                             'churn_seed_invocations': [[i % 3 + 1, 0] for i in range(6)]}}
    save(HERE / 'manifest.json', manifest)
    verify_all()
    print('Supplemental six-pair protocol frozen; original binaries verified; no timing.', flush=True)


def verify_all():
    verify()
    manifest = read(HERE / 'manifest.json')
    for name, digest in manifest['parent_artifacts'].items():
        assert sha(PARENT / name) == digest, name
    for name, digest in manifest['drivers'].items():
        assert sha(HERE / name) == digest, name
    assert sha(HERE / 'drivers.zip') == manifest['drivers_archive_sha256']
    built = read(PARENT / 'build.json')
    assert sha(PARENT / 'binaries.zip') == built['binaries_archive_sha256']
    binaries = {}
    with zipfile.ZipFile(PARENT / 'binaries.zip') as archive:
        for label in LABELS:
            binaries[label] = {}
            for kind in KINDS:
                data = archive.read(f'{label}/{kind}.exe')
                assert hashlib.sha256(data).hexdigest() == built[label][kind]['binary_sha256']
                assert hashlib.sha256(built[label][kind]['cargo_lock'].encode()).hexdigest() == sha(PARENT / 'workloads' / f'{kind}.Cargo.lock')
                path = HERE / '.cache' / label / f'{kind}.exe'
                if not path.exists():
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with path.open('xb') as output:
                        output.write(data)
                assert sha(path) == built[label][kind]['binary_sha256']
                binaries[label][kind] = str(path)
    return manifest, built, binaries


def run():
    manifest, built, binaries = verify_all()
    raw = HERE / 'raw'
    raw.mkdir()
    output = {'started_utc': now(), 'manifest_sha256': sha(HERE / 'manifest.json'),
              'parent_build_sha256': sha(PARENT / 'build.json'), 'runs': []}
    save(HERE / 'run_started.json', {key: value for key, value in output.items() if key != 'runs'})
    for kind in KINDS:
        for pair, order in enumerate(manifest['protocol']['order']):
            seed, invocation = manifest['protocol']['churn_seed_invocations'][pair]
            for label in order:
                print(f'{kind} global_pair={pair + 3} {label}', flush=True)
                arguments = [binaries[label][kind]]
                if kind == 'churn':
                    arguments += [str(seed), str(invocation)]
                started = now()
                timed_out = False
                try:
                    result = subprocess.run(arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            encoding='utf-8', errors='strict', timeout=180)
                except subprocess.TimeoutExpired as failure:
                    def decode(value):
                        return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value or ''
                    result = subprocess.CompletedProcess(arguments, None, decode(failure.stdout), decode(failure.stderr))
                    timed_out = True
                target = raw / f'{kind}_{pair}_{label}.tsv'
                for path, text in [(target, result.stdout), (target.with_suffix('.stderr'), result.stderr)]:
                    with path.open('x', encoding='utf-8', newline='\n') as handle:
                        handle.write(text)
                record = {'kind': kind, 'pair': pair, 'global_pair': pair + 3, 'label': label,
                          'order': order, 'seed': seed if kind == 'churn' else None,
                          'invocation': invocation if kind == 'churn' else None,
                          'arguments': arguments, 'started_utc': started, 'finished_utc': now(),
                          'returncode': result.returncode, 'timed_out': timed_out,
                          'raw': target.relative_to(HERE).as_posix(), 'raw_sha256': sha(target),
                          'stderr_sha256': sha(target.with_suffix('.stderr')),
                          'binary_sha256': built[label][kind]['binary_sha256']}
                if result.returncode == 0:
                    record.update(parse_micro(result.stdout) if kind == 'micro' else parse_churn(result.stdout, seed))
                save(target.with_suffix('.json'), record)
                assert result.returncode == 0, (target, result.stderr)
                output['runs'].append(record)
    output['finished_utc'] = now()
    verify_all()
    save(HERE / 'results.json', output)
    print('All24 supplemental processes complete; all samples and oracles retained.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['freeze', 'run', 'verify'])
    args = parser.parse_args()
    {'freeze': freeze, 'run': run, 'verify': verify_all}[args.phase]()
