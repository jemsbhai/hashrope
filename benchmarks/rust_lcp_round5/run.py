"""Freeze, build and run paired 0.3.3 LCP diagnostics without overwriting evidence."""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LABELS = ('baseline', 'candidate')
KINDS = ('micro', 'churn')


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as output:
        output.write(json.dumps(value, indent=2) + '\n')


def hashes(root):
    return {p.relative_to(root).as_posix(): sha(p) for p in sorted(root.rglob('*')) if p.is_file()}


def core(label):
    return HERE / '.cache' / label / 'hashrope'


def verify_inputs():
    baseline = read(HERE / 'baseline_manifest.json')
    assert sha(HERE / 'baseline.zip') == baseline['archive_sha256']
    assert sha(HERE / 'freeze_baseline.py') == baseline['freeze_script_sha256']
    workloads = read(HERE / 'workloads_manifest.json')
    for name, digest in workloads['workload_sha256'].items():
        assert sha(HERE / name) == digest, name
    return baseline, workloads


def freeze():
    baseline, _ = verify_inputs()
    if core('candidate').exists():
        raise RuntimeError('Candidate cache already exists; refusing to reuse unverified files')
    for name in ['candidate.zip', 'drivers.zip', 'manifest.json']:
        if (HERE / name).exists():
            raise RuntimeError(f'Refusing to overwrite {name}')
    with zipfile.ZipFile(HERE / 'candidate.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for relative in baseline['source_sha256']:
            archive.write(ROOT / relative, 'hashrope/' + relative)
    for label in LABELS:
        if not core(label).exists():
            with zipfile.ZipFile(HERE / f'{label}.zip') as archive:
                archive.extractall(core(label).parent)
    files = [HERE / name for name in ['PLAN.md', 'run.py', 'summarize.py', 'freeze_baseline.py',
                                      'baseline_manifest.json', 'workloads_manifest.json']]
    files += sorted((HERE / 'workloads').iterdir())
    with zipfile.ZipFile(HERE / 'drivers.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, path.relative_to(HERE).as_posix())
    sources = {label: hashes(core(label)) for label in LABELS}
    assert sources['baseline'] == baseline['source_sha256']
    assert sources['candidate'].keys() == sources['baseline'].keys()
    normalized = {label: {name: hashlib.sha256((core(label) / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                          for name in sources[label]} for label in LABELS}
    # Git stores this manifest with LF while a Windows working checkout can use
    # CRLF. Preserve both raw snapshots, and require exact normalized content.
    assert normalized['candidate']['packages/rust/Cargo.toml'] == normalized['baseline']['packages/rust/Cargo.toml']
    manifest = {
        'frozen_utc': now(), 'baseline_revision': baseline['source_revision'],
        'baseline_label': baseline['label'], 'candidate_label': 'reviewed Rust LCP candidate; unchanged 0.3.3 version',
        'archives': {name: sha(HERE / name) for name in ['baseline.zip', 'candidate.zip', 'drivers.zip']},
        'sources': sources, 'drivers': {path.relative_to(HERE).as_posix(): sha(path) for path in files},
        'changed_source_files': [name for name in sources['baseline'] if sources['baseline'][name] != sources['candidate'][name]],
        'newline_normalized_source_sha256': normalized,
        'newline_only_differences': [name for name in sources['baseline'] if sources['baseline'][name] != sources['candidate'][name] and normalized['baseline'][name] == normalized['candidate'][name]],
        'newline_normalized_changed_source_files': [name for name in sources['baseline'] if normalized['baseline'][name] != normalized['candidate'][name]],
        'protocol': {'pairs_per_cohort': 3, 'micro_cells': 42, 'micro_samples': 7,
                     'churn_cells': 24, 'churn_samples': 3, 'churn_seed_invocations': [[1, 0], [2, 0], [3, 0]],
                     'order': [list(LABELS), list(reversed(LABELS)), list(LABELS)], 'python_interleaving': False}}
    save(HERE / 'manifest.json', manifest)
    print(json.dumps({'frozen': True, 'changed_source_files': manifest['changed_source_files']}, indent=2))


def verify():
    verify_inputs()
    manifest = read(HERE / 'manifest.json')
    for name, digest in manifest['archives'].items():
        assert sha(HERE / name) == digest, name
    for name, digest in manifest['drivers'].items():
        assert sha(HERE / name) == digest, name
    for label in LABELS:
        if not core(label).exists():
            with zipfile.ZipFile(HERE / f'{label}.zip') as archive:
                archive.extractall(core(label).parent)
        assert hashes(core(label)) == manifest['sources'][label], label
    return manifest


def command(arguments, log):
    if log.exists():
        raise RuntimeError(f'Refusing to overwrite {log}')
    result = subprocess.run(arguments, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            encoding='utf-8', errors='replace')
    with log.open('x', encoding='utf-8', newline='\n') as output:
        output.write(result.stdout)
    if result.returncode:
        raise RuntimeError(f'Command exited {result.returncode}; see {log}')
    return result.stdout


def build():
    verify()
    if os.environ.get('CARGO_TARGET_DIR'):
        raise RuntimeError('Unset CARGO_TARGET_DIR so each frozen arm keeps its own target directory')
    save(HERE / 'build_started.json', {'started_utc': now(), 'manifest_sha256': sha(HERE / 'manifest.json')})
    output = {}
    for label in LABELS:
        output[label] = {}
        for kind in KINDS:
            project = HERE / '.cache/build' / f'{label}_{kind}'
            (project / 'src').mkdir(parents=True)
            source = HERE / 'workloads' / ('bench_rust.rs' if kind == 'micro' else 'edit_compare.rs')
            shutil.copyfile(source, project / 'src/main.rs')
            name = 'hashrope-comparison-benchmark' if kind == 'micro' else 'exact-edit-comparison'
            cargo = '[package]\nname=' + json.dumps(name) + '\nversion="0.0.0"\nedition="2021"\n'
            cargo += '[features]\nexact=[]\n[dependencies]\nhashrope={path=' + json.dumps((core(label) / 'packages/rust').as_posix()) + '}\n'
            if kind == 'churn':
                cargo += 'ropey="=1.6.1"\n'
            (project / 'Cargo.toml').write_text(cargo, encoding='utf-8', newline='\n')
            lock = HERE / 'workloads' / f'{kind}.Cargo.lock'
            shutil.copyfile(lock, project / 'Cargo.lock')
            arguments = ['cargo', 'build', '--release', '--offline', '--locked', '--manifest-path', str(project / 'Cargo.toml')]
            if kind == 'micro':
                arguments += ['--features', 'exact']
            command(arguments, HERE / f'build_{label}_{kind}.log')
            assert sha(project / 'Cargo.lock') == sha(lock), (label, kind, 'lock changed')
            binary = project / 'target/release' / (name + ('.exe' if os.name == 'nt' else ''))
            output[label][kind] = {'binary': str(binary), 'binary_sha256': sha(binary),
                'source_sha256': sha(project / 'src/main.rs'), 'cargo_manifest': cargo,
                'cargo_manifest_sha256': sha(project / 'Cargo.toml'), 'cargo_lock': lock.read_text(encoding='utf-8'),
                'cargo_lock_sha256': sha(lock), 'build_arguments': arguments}
    verify()
    output['environment'] = {
        'rustc': subprocess.check_output(['rustc', '--version', '--verbose'], text=True),
        'cargo': subprocess.check_output(['cargo', '--version'], text=True).strip(),
        'python': sys.version, 'platform': platform.platform(),
        'processor': os.environ.get('PROCESSOR_IDENTIFIER'), 'logical_processors': os.cpu_count(),
        'affinity': 'inherited unpinned; no process pinning',
        'compiler_env': {name: os.environ.get(name) for name in ['RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'CARGO_BUILD_TARGET']},
        'manifest_sha256': sha(HERE / 'manifest.json')}
    with zipfile.ZipFile(HERE / 'binaries.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for label in LABELS:
            for kind in KINDS:
                archive.write(output[label][kind]['binary'], f'{label}/{kind}.exe')
    output['binaries_archive_sha256'] = sha(HERE / 'binaries.zip')
    save(HERE / 'build.json', output)
    print('All four binaries built, locked and archived. No timings have run.', flush=True)


def parse_micro(output):
    cells, memory = {}, {}
    for line in output.splitlines():
        parts = line.split('\t')
        if parts[0] == 'time':
            _, name, sample, iterations, ns = parts
            row = cells.setdefault(name, {'samples_ns': [], 'iterations': int(iterations)})
            assert int(sample) == len(row['samples_ns']) and int(iterations) == row['iterations']
            assert math.isfinite(float(ns)) and float(ns) > 0
            row['samples_ns'].append(float(ns))
        elif parts[0] == 'memory':
            assert parts[1] not in memory
            memory[parts[1]] = int(parts[2])
        else:
            raise ValueError(line)
    expected = {'core/from_bytes_64KiB': 20, 'core/concat_balanced_64KiB': 2000,
                'core/split_chunked_64KiB': 100, 'core/substr_warm_64KiB': 1000,
                'core/root_hash': 1000000, 'core/repeat_join_2pow32': 100}
    for size in [65536, 1048576]:
        for label in ['same_root', 'independent_equal', 'early_mismatch', 'late_mismatch', 'shared_late_mismatch']:
            for operation in ['exact_lcp', 'exact_equal', 'hash_only_lcp']:
                expected[f'compare/{size}/{label}/{operation}'] = 2000 if operation != 'hash_only_lcp' and label in ['same_root', 'early_mismatch'] else 100
    for label in ['repeat_2TiB', 'self_concat_DAG_2TiB']:
        for operation in ['exact_lcp', 'exact_equal', 'hash_only_lcp']:
            expected[f'compare/{label}/{operation}'] = 2000 if label == 'repeat_2TiB' and operation != 'hash_only_lcp' else 100
    assert cells.keys() == expected.keys()
    assert memory.keys() == {'node_inner_bytes', 'arena_inline_bytes', 'nodes_for_64KiB'}
    for name, row in cells.items():
        assert row['iterations'] == expected[name]
        assert len(row['samples_ns']) == 7
        row['median_ns'] = statistics.median(row['samples_ns'])
    return {'cells': cells, 'memory': memory}


def parse_churn(output, seed):
    cells = {}
    for line in output.splitlines():
        tag, n, shared, arm, sample, ms = line.split('\t')
        assert int(tag) == seed
        assert int(n) in [10_000, 100_000, 1_000_000] and shared in ['false', 'true']
        assert arm in ['hashrope_exact', 'hashrope_candidate', 'ropey_exact', 'ropey_rehash_candidate']
        row = cells.setdefault(f'{n}/{shared}/{arm}', {'samples_ms': []})
        assert int(sample) == len(row['samples_ms']) and math.isfinite(float(ms)) and float(ms) > 0
        row['samples_ms'].append(float(ms))
    assert len(cells) == 24
    for row in cells.values():
        assert len(row['samples_ms']) == 3
        row['median_ms'] = statistics.median(row['samples_ms'])
    return {'cells': cells}


def verify_build():
    manifest = verify()
    built = read(HERE / 'build.json')
    assert built['environment']['manifest_sha256'] == sha(HERE / 'manifest.json')
    assert built['binaries_archive_sha256'] == sha(HERE / 'binaries.zip')
    for label in LABELS:
        for kind in KINDS:
            row = built[label][kind]
            assert sha(Path(row['binary'])) == row['binary_sha256']
            project = HERE / '.cache/build' / f'{label}_{kind}'
            for path, digest in [('src/main.rs', row['source_sha256']), ('Cargo.lock', row['cargo_lock_sha256']), ('Cargo.toml', row['cargo_manifest_sha256'])]:
                assert sha(project / path) == digest, (label, kind, path)
    return manifest, built


def run():
    manifest, built = verify_build()
    raw = HERE / 'raw'
    raw.mkdir()
    output = {'started_utc': now(), 'manifest_sha256': sha(HERE / 'manifest.json'),
              'build_sha256': sha(HERE / 'build.json'), 'runs': []}
    save(HERE / 'run_started.json', {key: value for key, value in output.items() if key != 'runs'})
    for kind in KINDS:
        for pair, order in enumerate(manifest['protocol']['order']):
            for label in order:
                print(f'{kind} pair={pair} {label}', flush=True)
                arguments = [built[label][kind]['binary']]
                if kind == 'churn':
                    arguments += [str(pair + 1), '0']
                started = now()
                timed_out = False
                try:
                    result = subprocess.run(arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            encoding='utf-8', errors='strict', timeout=180)
                except subprocess.TimeoutExpired as failure:
                    # Preserve captured output and an explicit unknown exit code;
                    # do not retry or substitute a completed process's samples.
                    def decode(value):
                        return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value or ''
                    result = subprocess.CompletedProcess(arguments, None, decode(failure.stdout), decode(failure.stderr))
                    timed_out = True
                target = raw / f'{kind}_{pair}_{label}.tsv'
                for path, text in [(target, result.stdout), (target.with_suffix('.stderr'), result.stderr)]:
                    with path.open('x', encoding='utf-8', newline='\n') as handle:
                        handle.write(text)
                record = {'kind': kind, 'pair': pair, 'label': label, 'order': order,
                          'arguments': arguments, 'started_utc': started, 'finished_utc': now(),
                          'returncode': result.returncode, 'timed_out': timed_out, 'raw': target.relative_to(HERE).as_posix(),
                          'raw_sha256': sha(target), 'stderr_sha256': sha(target.with_suffix('.stderr')),
                          'binary_sha256': built[label][kind]['binary_sha256']}
                if result.returncode == 0:
                    record.update(parse_micro(result.stdout) if kind == 'micro' else parse_churn(result.stdout, pair + 1))
                save(target.with_suffix('.json'), record)
                assert result.returncode == 0, (target, result.stderr)
                output['runs'].append(record)
    output['finished_utc'] = now()
    verify_build()
    save(HERE / 'results.json', output)
    print('All 12 fresh processes completed and verified; no samples removed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['freeze', 'build', 'run', 'verify'])
    arguments = parser.parse_args()
    {'freeze': freeze, 'build': build, 'run': run, 'verify': verify}[arguments.phase]()
