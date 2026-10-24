"""Build copied-source capacity/allocator probes only; no timed measurements."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
CACHE = HERE / '.cache/concat_capacity_probe'
ARMS = {'released': 'released_033_reference.zip', 'round5': 'baseline.zip'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2)
        handle.write('\n')


def prepare():
    CACHE.mkdir(parents=True, exist_ok=False)
    baseline = read(HERE / 'baseline_manifest.json')
    manifest = {'scope': 'Untimed allocator and private capacity observations in isolated source copies; not timings',
                'script_sha256': sha(Path(__file__)), 'harness_sha256': sha(HERE / 'concat_capacity_probe.rs'),
                'arms': {}}
    for label, archive in ARMS.items():
        assert sha(HERE / archive) == baseline['archives'][archive]
        directory = CACHE / label
        with zipfile.ZipFile(HERE / archive) as zipped:
            zipped.extractall(directory)
        core = directory / 'hashrope'
        sources = {p.relative_to(core).as_posix(): sha(p) for p in sorted(core.rglob('*')) if p.is_file()}
        rope = core / 'packages/rust/src/rope.rs'
        text = rope.read_text(encoding='utf-8')
        marker = '    pub fn node_count(&self) -> usize {\n        self.nodes.len()\n    }'
        assert text.count(marker) == 1
        text = text.replace(marker, marker + '\n\n    /// Diagnostic copy only.\n    pub fn diagnostic_node_capacity(&self) -> usize { self.nodes.capacity() }')
        rope.write_text(text, encoding='utf-8', newline='\n')
        project = directory / 'probe'
        (project / 'src').mkdir(parents=True)
        shutil.copyfile(HERE / 'concat_capacity_probe.rs', project / 'src/main.rs')
        (project / 'Cargo.toml').write_text('[package]\nname="hashrope-concat-capacity-probe"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="../hashrope/packages/rust"}\n', encoding='utf-8')
        manifest['arms'][label] = {'archive': archive, 'archive_sha256': sha(HERE / archive), 'original_sources': sources,
                                   'instrumented_rope_sha256': sha(rope), 'cargo_sha256': sha(project / 'Cargo.toml')}
    save(HERE / 'concat_capacity_probe_manifest.json', manifest)
    print('Prepared isolated capacity accessor copies; no builds or timings.')


def build():
    manifest = read(HERE / 'concat_capacity_probe_manifest.json')
    builds = {'compiler': subprocess.check_output(['rustc', '-Vv'], text=True), 'arms': {}}
    for label in ARMS:
        project = CACHE / label / 'probe'
        args = ['cargo', 'build', '--release', '--offline', '--manifest-path', str(project / 'Cargo.toml')]
        result = subprocess.run(args, capture_output=True, text=True)
        with (HERE / f'concat_capacity_build_{label}.log').open('x', encoding='utf-8') as handle:
            handle.write(result.stdout + result.stderr)
        assert result.returncode == 0
        binary = project / 'target/release/hashrope-concat-capacity-probe.exe'
        assert sha(CACHE / label / 'hashrope/packages/rust/src/rope.rs') == manifest['arms'][label]['instrumented_rope_sha256']
        builds['arms'][label] = {'arguments': args, 'binary': str(binary), 'binary_sha256': sha(binary),
                                'cargo_lock': (project / 'Cargo.lock').read_text()}
    save(HERE / 'concat_capacity_probe_build.json', builds)


def run():
    built = read(HERE / 'concat_capacity_probe_build.json')
    result = {'scope': 'Untimed requested allocation traffic and Vec capacities, not elapsed time or retained/peak process memory',
              'manifest_sha256': sha(HERE / 'concat_capacity_probe_manifest.json'),
              'build_sha256': sha(HERE / 'concat_capacity_probe_build.json'), 'arms': {}}
    for label in ARMS:
        binary = Path(built['arms'][label]['binary'])
        assert sha(binary) == built['arms'][label]['binary_sha256']
        command = subprocess.run([str(binary)], capture_output=True, text=True)
        with (HERE / f'concat_capacity_{label}.tsv').open('x', encoding='utf-8') as output:
            output.write(command.stdout)
        assert command.returncode == 0, command.stderr
        rows = []
        for line in command.stdout.splitlines():
            name, sample, count, before, cap_before, after, cap_after, root_len, root_weight, fingerprint, allocations, events = line.split('\t')
            counts = json.loads(allocations)
            rows.append({'setup': name, 'sample': int(sample), 'concat_calls': int(count),
                         'nodes_before': int(before), 'capacity_before': int(cap_before), 'nodes_after': int(after), 'capacity_after': int(cap_after),
                         'root_length': int(root_len), 'root_weight': int(root_weight), 'fingerprint': int(fingerprint),
                         'allocator': dict(zip(['allocations', 'reallocations', 'deallocations', 'requested_bytes', 'realloc_old_bytes', 'realloc_new_bytes'], counts)),
                         'realloc_events_old_bytes_new_bytes_moved': json.loads(events)})
        assert len(rows) == 16
        for setup, final in [('original_growth', 14513), ('prewarmed_no_growth', 64513)]:
            selected = [r for r in rows if r['setup'] == setup]
            assert [r['sample'] for r in selected] == list(range(-1, 7))
            assert selected[-1]['nodes_after'] == final
            assert all(r['nodes_after'] - r['nodes_before'] == r['concat_calls'] for r in selected)
        assert [r['allocator']['reallocations'] for r in rows if r['setup'] == 'original_growth' and r['sample'] >= 0] == [2, 1, 0, 1, 0, 0, 0]
        assert all(r['allocator']['reallocations'] == 0 for r in rows if r['setup'] == 'prewarmed_no_growth' and r['sample'] >= 0)
        result['arms'][label] = rows
    assert {r['fingerprint'] for rows in result['arms'].values() for r in rows}.__len__() == 1
    save(HERE / 'concat_capacity_probe.json', result)
    archive = HERE / 'concat_capacity_probe_sources.zip'
    files = []
    for label in ARMS:
        files.extend(p for p in (CACHE / label / 'hashrope').rglob('*') if p.is_file())
        files.extend(CACHE / label / 'probe' / name for name in ['src/main.rs', 'Cargo.toml', 'Cargo.lock'])
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as zipped:
        for path in files:
            zipped.write(path, path.relative_to(CACHE).as_posix())
    save(HERE / 'concat_capacity_probe_archive.json', {'archive_sha256': sha(archive), 'sources': {p.relative_to(CACHE).as_posix(): sha(p) for p in files}})
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('phase', choices=['prepare', 'build', 'run'])
    {'prepare': prepare, 'build': build, 'run': run}[p.parse_args().phase]()
