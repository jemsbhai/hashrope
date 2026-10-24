"""Recompute all paired cells and verify retained source/raw/binary provenance."""
import csv
import argparse
import io
import json
import math
from pathlib import Path
import statistics
import zipfile

from run import HERE, ROOT, LABELS, read, save, sha, verify, parse_micro, parse_churn


def main(csv_only=False):
    manifest = verify()
    build = read(HERE / 'build.json')
    result = read(HERE / 'results.json')
    assert result['manifest_sha256'] == sha(HERE / 'manifest.json')
    assert result['build_sha256'] == sha(HERE / 'build.json')
    assert sha(HERE / 'binaries.zip') == build['binaries_archive_sha256']
    with zipfile.ZipFile(HERE / 'binaries.zip') as binaries:
        import hashlib
        for label in LABELS:
            for kind in ('micro', 'churn'):
                assert hashlib.sha256(binaries.read(f'{label}/{kind}.exe')).hexdigest() == build[label][kind]['binary_sha256']
    expected = {(kind, pair, label) for kind in ('micro', 'churn') for pair in range(3) for label in LABELS}
    assert len(result['runs']) == len(expected)
    assert {(r['kind'], r['pair'], r['label']) for r in result['runs']} == expected
    indexed = {}
    for record in result['runs']:
        path = HERE / record['raw']
        assert sha(path) == record['raw_sha256']
        assert record['binary_sha256'] == build[record['label']][record['kind']]['binary_sha256']
        assert record['order'] == manifest['protocol']['order'][record['pair']]
        text = path.read_text(encoding='utf-8')
        if record['kind'] == 'churn':
            assert all(int(line.split('\t')[0]) == record['pair'] + 1 for line in text.splitlines())
        parsed = (parse_micro if record['kind'] == 'micro' else parse_churn)(text)
        for key, value in parsed.items(): assert record[key] == value
        indexed[record['kind'], record['pair'], record['label']] = record
    summary = {'scope': 'Targeted Rust 0.3.3 integration against completed round-four optimized sources; default features; no Python interleaving',
               'process_pairs_per_cohort': 3, 'source_sha256': sha(HERE / 'results.json'),
               'analysis_sha256': sha(Path(__file__)), 'cells': {}, 'memory': {},
               'rust_working_source_matches_candidate_snapshot': {}}
    for relative, digest in manifest['sources']['candidate'].items():
        summary['rust_working_source_matches_candidate_snapshot'][relative] = sha(ROOT / relative) == digest
    for kind in ('micro', 'churn'):
        unit = 'ns' if kind == 'micro' else 'ms'
        median_key = 'median_' + unit
        names = set(indexed[kind, 0, 'reference']['cells'])
        assert all(set(indexed[kind, pair, label]['cells']) == names for pair in range(3) for label in LABELS)
        rows = {}
        for name in sorted(names):
            values = {label: [indexed[kind, pair, label]['cells'][name][median_key] for pair in range(3)] for label in LABELS}
            ratios = [b / a for a, b in zip(values['reference'], values['candidate'])]
            rows[name] = {'unit': unit, 'process_medians': values,
                         'median_of_process_medians': {label: statistics.median(v) for label, v in values.items()},
                         'process_median_range': {label: [min(v), max(v)] for label, v in values.items()},
                         'paired_ratios_candidate_over_reference': ratios,
                         'paired_ratio_geomean': math.exp(statistics.mean(map(math.log, ratios))),
                         'paired_ratio_range': [min(ratios), max(ratios)],
                         'candidate_slower_pairs': sum(r > 1 for r in ratios)}
        summary['cells'][kind] = rows
    for label in LABELS:
        memory = [indexed['micro', pair, label]['memory'] for pair in range(3)]
        assert all(v == memory[0] for v in memory)
        summary['memory'][label] = memory[0]
    summary['verified_counts'] = {'fresh_processes': 12, 'micro_cells': 42, 'churn_cells': 24,
                                  'micro_timing_samples': 1764, 'churn_timing_samples': 432}
    if csv_only:
        previous = read(HERE / 'summary.json')
        old_analysis_hash = previous.pop('analysis_sha256')
        comparison = dict(summary)
        comparison.pop('analysis_sha256')
        assert previous == comparison, 'Derived numbers changed; CSV-only regeneration refused'
    else:
        save(HERE / 'summary.json', summary)
    output = io.StringIO()
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(['cohort', 'cell', 'unit', 'reference_median', 'candidate_median', 'ratio_geomean',
                     'ratio_min', 'ratio_max', 'slower_pairs', 'reference_pair_0', 'reference_pair_1', 'reference_pair_2',
                     'candidate_pair_0', 'candidate_pair_1', 'candidate_pair_2'])
    for kind, rows in summary['cells'].items():
        for name, row in rows.items():
            writer.writerow([kind, name, row['unit'], *[row['median_of_process_medians'][label] for label in LABELS],
                             row['paired_ratio_geomean'], *row['paired_ratio_range'], row['candidate_slower_pairs'],
                             *row['process_medians']['reference'], *row['process_medians']['candidate']])
    target = HERE / 'all_cells.csv'
    assert csv_only or not target.exists()
    target.write_text(output.getvalue(), encoding='utf-8', newline='\n')
    if csv_only:
        path = HERE / 'summary.json'
        text = path.read_text(encoding='utf-8')
        assert text.count(old_analysis_hash) == 1
        path.write_text(text.replace(old_analysis_hash, summary['analysis_sha256']), encoding='utf-8')
        print('Reverified all raw inputs; regenerated only CSV and refreshed its analysis-script hash.')
        return
    print(json.dumps({'memory': summary['memory'], 'verified_counts': summary['verified_counts'],
                      'recurring_slowdowns_over_10_percent': {kind: {name: row for name, row in rows.items()
                          if row['paired_ratio_geomean'] > 1.1 and row['candidate_slower_pairs'] == 3}
                          for kind, rows in summary['cells'].items()}}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv-only', action='store_true', help='Reverify raw inputs and replace only the derived CSV and analysis-script hash')
    main(csv_only=parser.parse_args().csv_only)
