"""Verify all frozen/raw records and retain every paired diagnostic outcome."""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import zipfile

from run import HERE, ROOT, LABELS, KINDS, read, save, sha, verify, parse_micro, parse_churn


def main():
    for name in ['summary.json', 'all_cells.csv']:
        if (HERE / name).exists():
            raise RuntimeError(f'Refusing to overwrite {name}')
    manifest = verify()
    built = read(HERE / 'build.json')
    results = read(HERE / 'results.json')
    assert results['manifest_sha256'] == built['environment']['manifest_sha256'] == sha(HERE / 'manifest.json')
    assert results['build_sha256'] == sha(HERE / 'build.json')
    assert built['binaries_archive_sha256'] == sha(HERE / 'binaries.zip')
    with zipfile.ZipFile(HERE / 'binaries.zip') as archive:
        for label in LABELS:
            for kind in KINDS:
                assert hashlib.sha256(archive.read(f'{label}/{kind}.exe')).hexdigest() == built[label][kind]['binary_sha256']
                assert hashlib.sha256(built[label][kind]['cargo_lock'].encode()).hexdigest() == sha(HERE / 'workloads' / f'{kind}.Cargo.lock')
    expected = [(kind, pair, label) for kind in KINDS for pair, order in enumerate(manifest['protocol']['order']) for label in order]
    assert [(r['kind'], r['pair'], r['label']) for r in results['runs']] == expected
    indexed = {}
    for record in results['runs']:
        path = HERE / record['raw']
        assert record == read(path.with_suffix('.json'))
        assert sha(path) == record['raw_sha256'] and sha(path.with_suffix('.stderr')) == record['stderr_sha256']
        assert record['returncode'] == 0 and record['timed_out'] is False
        assert record['binary_sha256'] == built[record['label']][record['kind']]['binary_sha256']
        assert record['order'] == manifest['protocol']['order'][record['pair']]
        expected_arguments = [built[record['label']][record['kind']]['binary']]
        if record['kind'] == 'churn':
            expected_arguments += [str(record['pair'] + 1), '0']
        assert record['arguments'] == expected_arguments
        raw = path.read_text(encoding='utf-8')
        parsed = parse_micro(raw) if record['kind'] == 'micro' else parse_churn(raw, record['pair'] + 1)
        for name, value in parsed.items():
            assert record[name] == value
        indexed[record['kind'], record['pair'], record['label']] = record
    summary = {'scope': 'Paired published Rust 0.3.3 versus frozen LCP candidate; default core features; unchanged micro/churn; no outer Python interleaving',
               'pairs_per_cohort': 3, 'result_sha256': sha(HERE / 'results.json'),
               'analysis_sha256': sha(Path(__file__)), 'cells': {}, 'memory': {},
               'working_source_matches_frozen_candidate': {name: sha(ROOT / name) == digest for name, digest in manifest['sources']['candidate'].items()}}
    for kind in KINDS:
        unit = 'ns' if kind == 'micro' else 'ms'
        median_key = 'median_' + unit
        names = set(indexed[kind, 0, 'baseline']['cells'])
        assert all(set(indexed[kind, pair, label]['cells']) == names for pair in range(3) for label in LABELS)
        rows = {}
        for name in sorted(names):
            values = {label: [indexed[kind, pair, label]['cells'][name][median_key] for pair in range(3)] for label in LABELS}
            ratios = [candidate / baseline for baseline, candidate in zip(values['baseline'], values['candidate'])]
            rows[name] = {
                'unit': unit, 'process_medians': values,
                'median_of_process_medians': {label: statistics.median(v) for label, v in values.items()},
                'process_median_range': {label: [min(v), max(v)] for label, v in values.items()},
                'paired_ratios_candidate_over_baseline': ratios,
                'paired_ratio_geomean': math.exp(statistics.mean(map(math.log, ratios))),
                'paired_ratio_range': [min(ratios), max(ratios)], 'candidate_slower_pairs': sum(r > 1 for r in ratios)}
        summary['cells'][kind] = rows
    for label in LABELS:
        memory = [indexed['micro', pair, label]['memory'] for pair in range(3)]
        assert all(value == memory[0] for value in memory)
        summary['memory'][label] = memory[0]
    summary['verified_counts'] = {'fresh_processes': 12, 'micro_cells': 42, 'churn_cells': 24,
                                  'micro_timing_samples': 1764, 'churn_timing_samples': 432}
    summary['recurring_slowdowns_over_10_percent'] = {
        kind: [name for name, row in rows.items() if row['paired_ratio_geomean'] > 1.1 and row['candidate_slower_pairs'] == 3]
        for kind, rows in summary['cells'].items()}
    save(HERE / 'summary.json', summary)
    output = io.StringIO()
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(['cohort', 'cell', 'unit', 'baseline_median', 'candidate_median', 'paired_ratio_geomean',
                     'ratio_min', 'ratio_max', 'slower_pairs', 'baseline_pair_0', 'baseline_pair_1', 'baseline_pair_2',
                     'candidate_pair_0', 'candidate_pair_1', 'candidate_pair_2'])
    for kind, rows in summary['cells'].items():
        for name, row in rows.items():
            writer.writerow([kind, name, row['unit'], *[row['median_of_process_medians'][label] for label in LABELS],
                             row['paired_ratio_geomean'], *row['paired_ratio_range'], row['candidate_slower_pairs'],
                             *row['process_medians']['baseline'], *row['process_medians']['candidate']])
    with (HERE / 'all_cells.csv').open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(output.getvalue())
    print(json.dumps({'verified_counts': summary['verified_counts'], 'memory': summary['memory'],
                      'recurring_slowdowns_over_10_percent': summary['recurring_slowdowns_over_10_percent']}, indent=2))


if __name__ == '__main__':
    main()
