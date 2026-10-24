"""Count exact-comparison source events and allocator calls; never benchmark time.

Instrumentation lives only in .cache/source_work_probe. The baseline comes from
the immutable study archive; the candidate is captured from the working source.
The atomic counters perturb generated code, so their elapsed time is irrelevant.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = HERE / '.cache/source_work_probe'
NAMES = [
    'dispatch_calls', 'fast_normalize_calls', 'fast_normalize_steps',
    'fast_recursive_calls', 'fast_recursive_iterations', 'boundary_calls',
    'period_root_steps', 'memo_normalize_calls', 'memo_normalize_steps',
    'memo_lcp_calls', 'memo_lcp_inner_calls', 'memo_cache_lookups',
    'memo_cache_insert_calls', 'leaf_compare_calls', 'leaf_compare_bytes',
    'aligned_branch_taken', 'equal_period_pieces',
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as output:
        json.dump(value, output, indent=2)
        output.write('\n')


def tick(name, amount='1'):
    return f'source_probe_tick({NAMES.index(name)}, {amount});'


def function_edit(text, name, transform, occurrence=0):
    """These known functions have balanced braces and no brace string literals."""
    needle = 'fn ' + name + '('
    start = -1
    for _ in range(occurrence + 1):
        start = text.index(needle, start + 1)
    begin = text.index('{', start)
    depth = 1
    end = begin + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[:begin + 1] + transform(text[begin + 1:end - 1]) + text[end - 1:]


def once(body, old, new):
    assert body.count(old) == 1, (old, body.count(old))
    return body.replace(old, new)


def instrument(source):
    prefix, comparison = source.split('// Exact comparison only reads', 1)
    comparison, suffix = comparison.split('// Free-function wrappers', 1)
    comparison = function_edit(comparison, 'exact_lcp', lambda b: tick('dispatch_calls') + b)
    comparison = function_edit(comparison, 'normalize', lambda b: tick('fast_normalize_calls') + once(b, 'loop {', 'loop {' + tick('fast_normalize_steps')), 0)
    comparison = function_edit(comparison, 'normalize', lambda b: tick('memo_normalize_calls') + once(b, 'loop {', 'loop {' + tick('memo_normalize_steps')), 1)
    comparison = function_edit(comparison, 'lcp_recursive', lambda b: tick('fast_recursive_calls') + once(b, 'while matched < length {', 'while matched < length {' + tick('fast_recursive_iterations')))
    comparison = function_edit(comparison, 'boundary', lambda b: tick('boundary_calls') + b)
    period_name = 'structural_period_root' if 'fn structural_period_root(' in comparison else 'structural_period'
    comparison = function_edit(comparison, period_name, lambda b: once(b, 'loop {', 'loop {' + tick('period_root_steps')))
    # FastComparison::lcp is occurrence 0; ExactComparison::lcp is occurrence 1.
    comparison = function_edit(comparison, 'lcp', lambda b: tick('memo_lcp_calls') + b, 1)
    comparison = function_edit(comparison, 'lcp_inner', lambda b: tick('memo_lcp_inner_calls') + b)
    comparison = once(comparison, 'if self.equal_ranges.contains(&key) {', tick('memo_cache_lookups') + 'if self.equal_ranges.contains(&key) {')
    assert comparison.count('self.equal_ranges.insert(') == 3
    comparison = comparison.replace('self.equal_ranges.insert(', tick('memo_cache_insert_calls') + 'self.equal_ranges.insert(')
    comparison = once(comparison, 'return Ok(if x == y {', tick('leaf_compare_calls') + tick('leaf_compare_bytes', 'x.len() as u64') + 'return Ok(if x == y {')
    assert comparison.count('if x != y {') == 2
    comparison = comparison.replace('if x != y {', tick('leaf_compare_calls') + tick('leaf_compare_bytes', 'x.len() as u64') + 'if x != y {')
    if 'let (al, ar, bl, br) =' in comparison:
        comparison = comparison.replace('let (al, ar, bl, br) =', tick('aligned_branch_taken') + 'let (al, ar, bl, br) =', 1)
    if 'fn equal_period_lcp(' in comparison:
        comparison = function_edit(comparison, 'equal_period_lcp', lambda b: once(b, 'while matched < period {', 'while matched < period {' + tick('equal_period_pieces')))
    support = f'''
static SOURCE_PROBE: [core::sync::atomic::AtomicU64; {len(NAMES)}] =
    [const {{ core::sync::atomic::AtomicU64::new(0) }}; {len(NAMES)}];
fn source_probe_tick(index: usize, amount: u64) {{
    SOURCE_PROBE[index].fetch_add(amount, core::sync::atomic::Ordering::Relaxed);
}}
pub fn source_probe_reset() {{
    for counter in &SOURCE_PROBE {{ counter.store(0, core::sync::atomic::Ordering::Relaxed); }}
}}
pub fn source_probe_counts() -> [u64; {len(NAMES)}] {{
    core::array::from_fn(|i| SOURCE_PROBE[i].load(core::sync::atomic::Ordering::Relaxed))
}}
'''
    return prefix + support + '// Exact comparison only reads' + comparison + '// Free-function wrappers' + suffix


def prepare():
    CACHE.mkdir(parents=True, exist_ok=False)
    baseline = json.loads((HERE / 'baseline_manifest.json').read_text())
    assert sha((HERE / 'baseline.zip').read_bytes()) == baseline['archive_sha256']
    manifest = {'scope': 'Source-event counts and query-time allocator traffic, not elapsed time or hardware performance counters',
                'counter_names': NAMES, 'harness_sha256': sha((HERE / 'source_work_probe.rs').read_bytes()),
                'instrumenter_sha256': sha(Path(__file__).read_bytes()), 'arms': {}}
    for label in ['baseline', 'candidate']:
        root = CACHE / label / 'hashrope'
        with zipfile.ZipFile(HERE / 'baseline.zip') as archive:
            for name in baseline['source_sha256']:
                data = archive.read('hashrope/' + name) if label == 'baseline' else (ROOT / name).read_bytes()
                dest = root / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
        sources = {p.relative_to(root).as_posix(): sha(p.read_bytes()) for p in sorted(root.rglob('*')) if p.is_file()}
        source_file = root / 'packages/rust/src/rope.rs'
        source_file.write_text(instrument(source_file.read_text()), encoding='utf-8', newline='\n')
        project = CACHE / label / 'probe'
        (project / 'src').mkdir(parents=True)
        shutil.copyfile(HERE / 'source_work_probe.rs', project / 'src/main.rs')
        (project / 'Cargo.toml').write_text('[package]\nname="hashrope-source-work-probe"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="../hashrope/packages/rust"}\n', encoding='utf-8')
        manifest['arms'][label] = {'sources': sources, 'instrumented_rope_sha256': sha(source_file.read_bytes()),
                                   'cargo_sha256': sha((project / 'Cargo.toml').read_bytes())}
    save(HERE / 'source_work_probe_manifest.json', manifest)
    print('Prepared source-event probes; no compiler or timings executed.')


def build():
    manifest = json.loads((HERE / 'source_work_probe_manifest.json').read_text())
    builds = {}
    for label in ['baseline', 'candidate']:
        project = CACHE / label / 'probe'
        arguments = ['cargo', 'build', '--release', '--offline', '--manifest-path', str(project / 'Cargo.toml')]
        result = subprocess.run(arguments, capture_output=True, text=True)
        with (HERE / f'source_work_probe_build_{label}.log').open('x', encoding='utf-8') as out:
            out.write(result.stdout + result.stderr)
        assert result.returncode == 0
        binary = project / 'target/release/hashrope-source-work-probe.exe'
        builds[label] = {'arguments': arguments, 'binary': str(binary), 'binary_sha256': sha(binary.read_bytes()),
                         'cargo_lock': (project / 'Cargo.lock').read_text(),
                         'instrumented_rope_sha256': sha((CACHE / label / 'hashrope/packages/rust/src/rope.rs').read_bytes())}
        assert builds[label]['instrumented_rope_sha256'] == manifest['arms'][label]['instrumented_rope_sha256']
    builds['compiler'] = subprocess.check_output(['rustc', '-Vv'], text=True)
    save(HERE / 'source_work_probe_build.json', builds)


def run():
    builds = json.loads((HERE / 'source_work_probe_build.json').read_text())
    manifest = json.loads((HERE / 'source_work_probe_manifest.json').read_text())
    result = {'scope': manifest['scope'], 'manifest_sha256': sha((HERE / 'source_work_probe_manifest.json').read_bytes()),
              'build_sha256': sha((HERE / 'source_work_probe_build.json').read_bytes()), 'arms': {}}
    for label in ['baseline', 'candidate']:
        binary = Path(builds[label]['binary'])
        assert sha(binary.read_bytes()) == builds[label]['binary_sha256']
        command = subprocess.run([str(binary)], capture_output=True, text=True)
        with (HERE / f'source_work_probe_{label}.tsv').open('x', encoding='utf-8') as output:
            output.write(command.stdout)
        assert command.returncode == 0, command.stderr
        rows = {}
        for line in command.stdout.splitlines():
            name, answer, expected, alloc, realloc, dealloc, alloc_bytes, counters = line.split('\t')
            values = json.loads(counters)
            assert len(values) == len(NAMES) and int(answer) == int(expected)
            rows[name] = {'lcp': int(answer), 'expected': int(expected), 'allocations': int(alloc),
                          'reallocations': int(realloc), 'deallocations': int(dealloc), 'requested_allocation_bytes': int(alloc_bytes),
                          'source_events': dict(zip(NAMES, values))}
        result['arms'][label] = rows
    assert result['arms']['baseline'].keys() == result['arms']['candidate'].keys()
    save(HERE / 'source_work_probe.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['prepare', 'build', 'run'])
    {'prepare': prepare, 'build': build, 'run': run}[parser.parse_args().phase]()
