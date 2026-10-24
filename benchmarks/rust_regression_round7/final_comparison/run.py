"""Freeze, build, measure and audit the final round-seven optimization against round six."""
import argparse
import csv
import hashlib
import importlib.util
import io
import math
from pathlib import Path
import statistics
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parents[1] / 'rust_lcp_round5'
ROOT = HERE.parents[2]
import shutil
spec = importlib.util.spec_from_file_location('round5_driver', REFERENCE / 'run.py')
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
read, save, sha, now = driver.read, driver.save, driver.sha, driver.now
LABELS, KINDS = driver.LABELS, driver.KINDS
ORDERS = [list(LABELS) if i % 2 == 0 else list(reversed(LABELS)) for i in range(12)]


def core(label):
    return HERE / '.cache' / label / 'hashrope'


def freeze():
    for name in ['candidate.zip', 'baseline.zip', 'drivers.zip', 'manifest.json', '.cache']:
        if (HERE / name).exists():
            raise RuntimeError(f'Refusing to overwrite partial or completed freeze: {name}')
    parent = read(HERE.parent / 'baseline_manifest.json')
    assert sha(HERE.parent / 'baseline.zip') == parent['archives']['baseline.zip']
    original_workloads = read(REFERENCE / 'workloads_manifest.json')
    for name, digest in original_workloads['workload_sha256'].items():
        assert sha(HERE / name) == digest, name
    with zipfile.ZipFile(HERE / 'candidate.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for name in parent['source_sha256']['baseline']:
            archive.write(ROOT / name, 'hashrope/' + name)
    shutil.copyfile(HERE.parent / 'baseline.zip', HERE / 'baseline.zip')
    for label in LABELS:
        with zipfile.ZipFile(HERE / f'{label}.zip') as archive:
            archive.extractall(core(label).parent)
    sources = {label: driver.hashes(core(label)) for label in LABELS}
    assert sources['baseline'] == parent['source_sha256']['baseline']
    assert sources['candidate'].keys() == sources['baseline'].keys()
    normalized = {label: {name: hashlib.sha256((core(label) / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                          for name in sources[label]} for label in LABELS}
    changed = [name for name in sources['baseline'] if normalized['baseline'][name] != normalized['candidate'][name]]
    assert changed == ['packages/rust/src/rope.rs'], changed
    driver_paths = [HERE / 'PLAN.md', HERE / 'run.py'] + sorted((HERE / 'workloads').iterdir())
    with zipfile.ZipFile(HERE / 'drivers.zip', 'x', zipfile.ZIP_DEFLATED) as archive:
        for path in driver_paths:
            archive.write(path, path.relative_to(HERE).as_posix())
    save(HERE / 'manifest.json', {
        'frozen_utc': now(),
        'labels': {'baseline': 'frozen round-six optimization', 'candidate': 'round-seven optimization; version unchanged'},
        'parent_manifest_sha256': sha(HERE.parent / 'baseline_manifest.json'),
        'imported_builder_sha256': sha(REFERENCE / 'run.py'),
        'archives': {name: sha(HERE / name) for name in ['baseline.zip', 'candidate.zip', 'drivers.zip']},
        'sources': sources, 'normalized_changed_source_files': changed,
        'drivers': {path.relative_to(HERE).as_posix(): sha(path) for path in driver_paths},
        'protocol': {'orders': ORDERS, 'churn_seeds': [1, 2, 3] * 4,
                     'invocation': 0, 'pairs_per_cohort': 12}})
    verify_sources()
    print('Final source and twelve-pair protocol frozen; no timings.', flush=True)


def verify_sources():
    manifest = read(HERE / 'manifest.json')
    assert sha(HERE.parent / 'baseline_manifest.json') == manifest['parent_manifest_sha256']
    assert sha(REFERENCE / 'run.py') == manifest['imported_builder_sha256']
    for name, digest in manifest['archives'].items():
        assert sha(HERE / name) == digest, name
    for name, digest in manifest['drivers'].items():
        assert sha(HERE / name) == digest, name
    for label in LABELS:
        if not core(label).exists():
            with zipfile.ZipFile(HERE / f'{label}.zip') as archive:
                archive.extractall(core(label).parent)
        assert driver.hashes(core(label)) == manifest['sources'][label], label
    assert manifest['protocol']['orders'] == ORDERS
    assert manifest['protocol']['churn_seeds'] == [1, 2, 3] * 4
    return manifest


def build():
    # Reuse the reviewed locked build procedure; all paths and verification are
    # rebound to this study. Parsers remain byte-identical to round five.
    driver.HERE, driver.ROOT = HERE, ROOT
    driver.verify = verify_sources
    driver.build()


def verify():
    manifest = verify_sources()
    built = read(HERE / 'build.json')
    assert built['environment']['manifest_sha256'] == sha(HERE / 'manifest.json')
    assert built['binaries_archive_sha256'] == sha(HERE / 'binaries.zip')
    binaries = {}
    with zipfile.ZipFile(HERE / 'binaries.zip') as archive:
        for label in LABELS:
            for kind in KINDS:
                key = f'{label}/{kind}.exe'
                data = archive.read(key)
                assert hashlib.sha256(data).hexdigest() == built[label][kind]['binary_sha256']
                path = Path(built[label][kind]['binary'])
                if not path.exists():
                    path = HERE / '.cache/archived' / key
                    if not path.exists():
                        path.parent.mkdir(parents=True, exist_ok=True)
                        with path.open('xb') as output:
                            output.write(data)
                assert sha(path) == built[label][kind]['binary_sha256']
                binaries[label, kind] = str(path)
    return manifest, built, binaries


def parse(kind, text, seed):
    return driver.parse_micro(text) if kind == 'micro' else driver.parse_churn(text, seed)


def run():
    _, built, binaries = verify()
    raw = HERE / 'raw'
    raw.mkdir()
    output = {'started_utc': now(), 'manifest_sha256': sha(HERE / 'manifest.json'), 'build_sha256': sha(HERE / 'build.json'), 'runs': []}
    save(HERE / 'run_started.json', {k: v for k, v in output.items() if k != 'runs'})
    for kind in KINDS:
        for pair, order in enumerate(ORDERS):
            seed = pair % 3 + 1
            for label in order:
                arguments = [binaries[label, kind]] + ([str(seed), '0'] if kind == 'churn' else [])
                print(f'{kind} pair={pair} {label}', flush=True)
                started, timed_out = now(), False
                try:
                    result = subprocess.run(arguments, capture_output=True, encoding='utf-8',
                                            errors='strict', timeout=180)
                except subprocess.TimeoutExpired as failure:
                    def decode(value):
                        return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value or ''
                    result = subprocess.CompletedProcess(arguments, None, decode(failure.stdout), decode(failure.stderr))
                    timed_out = True
                path = raw / f'{kind}_{pair}_{label}.tsv'
                for target, value in [(path, result.stdout), (path.with_suffix('.stderr'), result.stderr)]:
                    with target.open('x', encoding='utf-8', newline='\n') as handle:
                        handle.write(value)
                row = {'kind': kind, 'pair': pair, 'label': label, 'order': order,
                       'seed': seed if kind == 'churn' else None,
                       'arguments': arguments, 'started_utc': started, 'finished_utc': now(),
                       'returncode': result.returncode, 'timed_out': timed_out,
                       'raw': path.relative_to(HERE).as_posix(), 'raw_sha256': sha(path),
                       'stderr_sha256': sha(path.with_suffix('.stderr')),
                       'binary_sha256': built[label][kind]['binary_sha256']}
                if result.returncode == 0:
                    row.update(parse(kind, result.stdout, seed))
                save(path.with_suffix('.json'), row)
                assert result.returncode == 0, (path, result.stderr)
                output['runs'].append(row)
    verify()
    output['finished_utc'] = now()
    save(HERE / 'results.json', output)
    print('All 48 final-comparison processes completed; all samples retained.', flush=True)


def geomean(values):
    return math.exp(statistics.mean(map(math.log, values)))


def summarize():
    _, built, binaries = verify()
    results = read(HERE / 'results.json')
    assert results['manifest_sha256'] == sha(HERE / 'manifest.json')
    assert results['build_sha256'] == sha(HERE / 'build.json')
    expected = [(kind, pair, label) for kind in KINDS for pair, order in enumerate(ORDERS) for label in order]
    assert [(r['kind'], r['pair'], r['label']) for r in results['runs']] == expected
    indexed = {}
    for row in results['runs']:
        kind, pair, label = row['kind'], row['pair'], row['label']
        path = HERE / row['raw']
        assert row == read(path.with_suffix('.json'))
        assert sha(path) == row['raw_sha256'] and sha(path.with_suffix('.stderr')) == row['stderr_sha256']
        assert row['returncode'] == 0 and not row['timed_out']
        assert row['order'] == ORDERS[pair]
        assert row['seed'] == (pair % 3 + 1 if kind == 'churn' else None)
        assert row['arguments'][0] in [built[label][kind]['binary'],
            str(HERE / '.cache/archived' / label / f'{kind}.exe')]
        assert row['arguments'][1:] == ([str(row['seed']), '0'] if kind == 'churn' else [])
        assert row['binary_sha256'] == built[label][kind]['binary_sha256']
        parsed = parse(kind, path.read_text(encoding='utf-8'), row['seed'])
        for name, value in parsed.items():
            assert row[name] == value
        indexed[kind, pair, label] = row
    summary = {'scope': 'Twelve fresh unpinned pairs per cohort: round-six baseline versus round-seven candidate; unchanged workloads',
               'results_sha256': sha(HERE / 'results.json'), 'analysis_sha256': sha(Path(__file__)),
               'processes': 48, 'samples': 8784,
               'working_source_matches_candidate': {name: sha(ROOT / name) == digest for name, digest in read(HERE / 'manifest.json')['sources']['candidate'].items()}, 'cells': {}, 'memory': {}}
    for kind in KINDS:
        summary['cells'][kind] = {}
        unit = 'ns' if kind == 'micro' else 'ms'
        names = indexed[kind, 0, 'baseline']['cells'].keys()
        for name in sorted(names):
            values = {label: [indexed[kind, pair, label]['cells'][name]['median_' + unit]
                              for pair in range(12)] for label in LABELS}
            ratios = [b / a for a, b in zip(values['baseline'], values['candidate'])]
            logs = list(map(math.log, ratios))
            center = statistics.mean(logs)
            half_width = 2.200985160082949 * statistics.stdev(logs) / math.sqrt(12)
            summary['cells'][kind][name] = {
                'unit': unit, 'process_medians': values,
                'median_of_process_medians': {label: statistics.median(v) for label, v in values.items()},
                'paired_ratios': ratios, 'paired_ratio_geomean': math.exp(center),
                'paired_ratio_range': [min(ratios), max(ratios)],
                'exploratory_95_log_t_interval': [math.exp(center - half_width), math.exp(center + half_width)],
                'candidate_slower_pairs': sum(x > 1 for x in ratios),
                'order_geomeans': {order: geomean(ratios[i] for i in range(12) if ORDERS[i][0] == order)
                                   for order in LABELS},
                **({'seed_geomeans': {str(seed): geomean(ratios[i] for i in range(12) if i % 3 + 1 == seed)
                                      for seed in [1, 2, 3]}} if kind == 'churn' else {})}
    for label in LABELS:
        values = [indexed['micro', pair, label]['memory'] for pair in range(12)]
        assert all(v == values[0] for v in values)
        summary['memory'][label] = values[0]
    save(HERE / 'summary.json', summary)
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(['cohort', 'cell', 'unit', 'baseline_median', 'candidate_median', 'paired_ratio', 'ci_low', 'ci_high', 'slower_pairs', 'ratios'])
    for kind, rows in summary['cells'].items():
        for name, row in rows.items():
            writer.writerow([kind, name, row['unit'], *[row['median_of_process_medians'][x] for x in LABELS],
                             row['paired_ratio_geomean'], *row['exploratory_95_log_t_interval'],
                             row['candidate_slower_pairs'], row['paired_ratios']])
    with (HERE / 'all_cells.csv').open('x', encoding='utf-8', newline='\n') as output:
        output.write(stream.getvalue())
    print('All 48 processes and 66 paired summaries independently reaudited.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['freeze', 'build', 'run', 'summarize', 'verify'])
    args = parser.parse_args()
    {'freeze': freeze, 'build': build, 'run': run, 'summarize': summarize, 'verify': verify}[args.phase]()
