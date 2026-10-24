"""Report supplemental six and all nine pairs separately, retaining original data."""
import csv
import io
import json
import math
import statistics

from confirm import HERE, PARENT, LABELS, KINDS, verify_all, read, save, sha, parse_micro, parse_churn


def audit(folder, results, orders, supplemental):
    expected = [(kind, pair, label) for kind in KINDS for pair, order in enumerate(orders) for label in order]
    assert [(r['kind'], r['pair'], r['label']) for r in results['runs']] == expected
    built = read(PARENT / 'build.json')
    indexed = {}
    for row in results['runs']:
        path = folder / row['raw']
        assert row == read(path.with_suffix('.json'))
        assert sha(path) == row['raw_sha256'] and sha(path.with_suffix('.stderr')) == row['stderr_sha256']
        assert row['returncode'] == 0 and row['timed_out'] is False
        assert row['binary_sha256'] == built[row['label']][row['kind']]['binary_sha256']
        assert row['order'] == orders[row['pair']]
        global_pair = row['pair'] + (3 if supplemental else 0)
        seed = row['pair'] % 3 + 1
        binary = str(HERE / '.cache' / row['label'] / f'{row["kind"]}.exe') if supplemental else built[row['label']][row['kind']]['binary']
        assert row['arguments'] == [binary] + ([str(seed), '0'] if row['kind'] == 'churn' else [])
        if supplemental:
            assert row['global_pair'] == global_pair
            assert row['seed'] == (seed if row['kind'] == 'churn' else None)
            assert row['invocation'] == (0 if row['kind'] == 'churn' else None)
        text = path.read_text(encoding='utf-8')
        parsed = parse_micro(text) if row['kind'] == 'micro' else parse_churn(text, seed)
        for key, value in parsed.items():
            assert row[key] == value
        indexed[row['kind'], global_pair, row['label']] = row
    return indexed


def aggregate(indexed, pairs):
    cells = {}
    for kind in KINDS:
        unit = 'ns' if kind == 'micro' else 'ms'
        names = indexed[kind, pairs[0], 'baseline']['cells'].keys()
        assert all(indexed[kind, pair, label]['cells'].keys() == names for pair in pairs for label in LABELS)
        cells[kind] = {}
        for name in sorted(names):
            values = {label: [indexed[kind, pair, label]['cells'][name]['median_' + unit] for pair in pairs] for label in LABELS}
            ratios = [b / a for a, b in zip(values['baseline'], values['candidate'])]
            cells[kind][name] = {'unit': unit, 'global_pairs': pairs, 'process_medians': values,
                'median_of_process_medians': {label: statistics.median(v) for label, v in values.items()},
                'process_median_range': {label: [min(v), max(v)] for label, v in values.items()},
                'paired_ratios_candidate_over_baseline': ratios,
                'paired_ratio_geomean': math.exp(statistics.mean(map(math.log, ratios))),
                'paired_ratio_range': [min(ratios), max(ratios)], 'candidate_slower_pairs': sum(r > 1 for r in ratios)}
    return cells


def main():
    if (HERE / 'summary.json').exists() or (HERE / 'all_cells.csv').exists():
        raise RuntimeError('Refusing to overwrite derived evidence')
    manifest, _, _ = verify_all()
    original_manifest = read(PARENT / 'manifest.json')
    original = read(PARENT / 'results.json')
    additional = read(HERE / 'results.json')
    assert original['manifest_sha256'] == sha(PARENT / 'manifest.json')
    assert original['build_sha256'] == additional['parent_build_sha256'] == sha(PARENT / 'build.json')
    assert additional['manifest_sha256'] == sha(HERE / 'manifest.json')
    indexed = audit(PARENT, original, original_manifest['protocol']['order'], False)
    indexed.update(audit(HERE, additional, manifest['protocol']['order'], True))
    original_recomputed = aggregate(indexed, list(range(3)))
    old_summary = read(PARENT / 'summary.json')
    for kind in KINDS:
        for name, row in original_recomputed[kind].items():
            expected = dict(row)
            expected.pop('global_pairs')
            assert expected == old_summary['cells'][kind][name]
    summary = {'scope': 'Same original executables; supplemental six pairs separately and pooled all nine; repeated churn seeds across two phases',
               'source_sha256': {'original_results': sha(PARENT / 'results.json'), 'supplemental_results': sha(HERE / 'results.json')},
               'analysis_sha256': sha(HERE / 'summarize_confirmation.py'),
               'supplemental_six': aggregate(indexed, list(range(3, 9))),
               'pooled_nine': aggregate(indexed, list(range(9))), 'memory': {},
               'verified_counts': {'original_processes': 12, 'supplemental_processes': 24,
                    'supplemental_micro_samples': 3528, 'supplemental_churn_samples': 864,
                    'pooled_micro_samples': 5292, 'pooled_churn_samples': 1296}}
    for label in LABELS:
        values = [indexed['micro', pair, label]['memory'] for pair in range(9)]
        assert all(v == values[0] for v in values)
        summary['memory'][label] = values[0]
    save(HERE / 'summary.json', summary)
    output = io.StringIO()
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(['phase', 'cohort', 'cell', 'unit', 'baseline_median', 'candidate_median', 'paired_ratio_geomean', 'ratio_min', 'ratio_max', 'slower_pairs', 'total_pairs', 'paired_ratios'])
    for phase in ['supplemental_six', 'pooled_nine']:
        for kind, rows in summary[phase].items():
            for name, row in rows.items():
                writer.writerow([phase, kind, name, row['unit'], *[row['median_of_process_medians'][label] for label in LABELS],
                    row['paired_ratio_geomean'], *row['paired_ratio_range'], row['candidate_slower_pairs'], len(row['global_pairs']),
                    json.dumps(row['paired_ratios_candidate_over_baseline'])])
    with (HERE / 'all_cells.csv').open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(output.getvalue())
    print(json.dumps(summary['verified_counts'], indent=2))


if __name__ == '__main__':
    main()
