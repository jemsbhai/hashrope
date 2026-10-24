"""Reexecute the exact round-five binaries under a fixed six-pair protocol."""
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
spec = importlib.util.spec_from_file_location('round5_driver', REFERENCE / 'run.py')
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
read, save, sha, now = driver.read, driver.save, driver.sha, driver.now
LABELS, KINDS = driver.LABELS, driver.KINDS
ORDERS = [list(LABELS) if i % 2 == 0 else list(reversed(LABELS)) for i in range(6)]


def freeze():
    driver.verify()
    paths = ['baseline.zip', 'candidate.zip', 'drivers.zip', 'manifest.json',
             'build.json', 'binaries.zip', 'run.py', 'results.json', 'summary.json',
             'confirmation/results.json', 'confirmation/summary.json']
    save(HERE / 'manifest.json', {
        'declared_utc': now(),
        'labels': {'baseline': 'published 0.3.3', 'candidate': 'round-five LCP optimization'},
        'reference_artifacts': {name: sha(REFERENCE / name) for name in paths},
        'drivers': {name: sha(HERE / name) for name in ['PLAN.md', 'run.py']},
        'protocol': {'orders': ORDERS, 'churn_seeds': [1, 2, 3, 1, 2, 3],
                     'invocation': 0, 'pairs_per_cohort': 6}})
    verify()
    print('Protocol frozen; original executable bytes verified; no timings.', flush=True)


def verify():
    driver.verify()
    manifest = read(HERE / 'manifest.json')
    for name, digest in manifest['reference_artifacts'].items():
        assert sha(REFERENCE / name) == digest, name
    for name, digest in manifest['drivers'].items():
        assert sha(HERE / name) == digest, name
    assert manifest['protocol']['orders'] == ORDERS
    assert manifest['protocol']['churn_seeds'] == [1, 2, 3, 1, 2, 3]
    built = read(REFERENCE / 'build.json')
    binaries = {}
    with zipfile.ZipFile(REFERENCE / 'binaries.zip') as archive:
        for label in LABELS:
            for kind in KINDS:
                key = f'{label}/{kind}.exe'
                data = archive.read(key)
                assert hashlib.sha256(data).hexdigest() == built[label][kind]['binary_sha256']
                path = HERE / '.cache' / key
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
    output = {'started_utc': now(), 'manifest_sha256': sha(HERE / 'manifest.json'), 'runs': []}
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
    print('All 24 original-binary processes completed; all samples retained.', flush=True)


def geomean(values):
    return math.exp(statistics.mean(map(math.log, values)))


def summarize():
    _, built, binaries = verify()
    results = read(HERE / 'results.json')
    assert results['manifest_sha256'] == sha(HERE / 'manifest.json')
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
        assert row['arguments'] == [binaries[label, kind]] + ([str(row['seed']), '0'] if kind == 'churn' else [])
        assert row['binary_sha256'] == built[label][kind]['binary_sha256']
        parsed = parse(kind, path.read_text(encoding='utf-8'), row['seed'])
        for name, value in parsed.items():
            assert row[name] == value
        indexed[kind, pair, label] = row
    summary = {'scope': 'Six new unpinned original-executable pairs per cohort; no pooling with round five',
               'results_sha256': sha(HERE / 'results.json'), 'analysis_sha256': sha(Path(__file__)),
               'processes': 24, 'samples': 4392, 'cells': {}, 'memory': {}}
    for kind in KINDS:
        summary['cells'][kind] = {}
        unit = 'ns' if kind == 'micro' else 'ms'
        names = indexed[kind, 0, 'baseline']['cells'].keys()
        for name in sorted(names):
            values = {label: [indexed[kind, pair, label]['cells'][name]['median_' + unit]
                              for pair in range(6)] for label in LABELS}
            ratios = [b / a for a, b in zip(values['baseline'], values['candidate'])]
            logs = list(map(math.log, ratios))
            center = statistics.mean(logs)
            half_width = 2.5705818366147395 * statistics.stdev(logs) / math.sqrt(6)
            summary['cells'][kind][name] = {
                'unit': unit, 'process_medians': values,
                'median_of_process_medians': {label: statistics.median(v) for label, v in values.items()},
                'paired_ratios': ratios, 'paired_ratio_geomean': math.exp(center),
                'paired_ratio_range': [min(ratios), max(ratios)],
                'exploratory_95_log_t_interval': [math.exp(center - half_width), math.exp(center + half_width)],
                'candidate_slower_pairs': sum(x > 1 for x in ratios),
                'order_geomeans': {order: geomean(ratios[i] for i in range(6) if ORDERS[i][0] == order)
                                   for order in LABELS},
                **({'seed_geomeans': {str(seed): geomean(ratios[i] for i in range(6) if i % 3 + 1 == seed)
                                      for seed in [1, 2, 3]}} if kind == 'churn' else {})}
    for label in LABELS:
        values = [indexed['micro', pair, label]['memory'] for pair in range(6)]
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
    print('All 24 processes and 66 paired summaries independently reaudited.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['freeze', 'run', 'summarize', 'verify'])
    args = parser.parse_args()
    {'freeze': freeze, 'run': run, 'summarize': summarize, 'verify': verify}[args.phase]()
