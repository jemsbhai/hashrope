"""Read-only reanalysis of every retained round-five pair; no new measurements."""
from pathlib import Path
import hashlib
import json
import math
import statistics
import zipfile

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'rust_lcp_round5'
LABELS = ('baseline', 'candidate')
KINDS = ('micro', 'churn')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def geomean(values):
    return math.exp(statistics.mean(map(math.log, values)))


def stats(values):
    return {'count': len(values), 'geomean': geomean(values),
            'range': [min(values), max(values)], 'above_one': sum(x > 1 for x in values)}


def collect():
    built = read(PRIOR / 'build.json')
    rows = {}
    samples = 0
    identities = []
    for phase, folder, count, shift in [('initial_three', PRIOR, 3, 0), ('supplemental_six', PRIOR / 'confirmation', 6, 3)]:
        results = read(folder / 'results.json')
        expected_order = [(kind, pair, label) for kind in KINDS for pair in range(count)
                          for label in (LABELS if (pair + shift) % 2 == 0 else tuple(reversed(LABELS)))]
        assert [(r['kind'], r['pair'], r['label']) for r in results['runs']] == expected_order
        for record in results['runs']:
            raw = folder / record['raw']
            assert sha(raw) == record['raw_sha256']
            assert sha(raw.with_suffix('.stderr')) == record['stderr_sha256']
            assert record == read(raw.with_suffix('.json'))
            assert record['returncode'] == 0 and not record['timed_out']
            assert record['binary_sha256'] == built[record['label']][record['kind']]['binary_sha256']
            cells = {}
            for line in raw.read_text(encoding='utf-8').splitlines():
                fields = line.split('\t')
                if fields[0] == 'memory':
                    continue
                if record['kind'] == 'micro':
                    assert fields[0] == 'time'
                    name, sample, value = fields[1], int(fields[2]), float(fields[4])
                else:
                    assert int(fields[0]) == record['pair'] % 3 + 1
                    name, sample, value = '/'.join(fields[1:4]), int(fields[4]), float(fields[5])
                values = cells.setdefault(name, [])
                assert sample == len(values) and math.isfinite(value) and value > 0
                values.append(value)
                samples += 1
            assert len(cells) == (42 if record['kind'] == 'micro' else 24)
            assert all(len(v) == (7 if record['kind'] == 'micro' else 3) for v in cells.values())
            global_pair = record['pair'] + shift
            rows[record['kind'], global_pair, record['label']] = {
                'phase': phase, 'global_pair': global_pair, 'cells': cells,
                'order': record['order'], 'started_utc': record['started_utc'], 'finished_utc': record['finished_utc'],
                'seed': record['pair'] % 3 + 1 if record['kind'] == 'churn' else None,
            }
            identities.append({'raw': str(raw.relative_to(PRIOR)), 'sha256': sha(raw)})
    with zipfile.ZipFile(PRIOR / 'binaries.zip') as archive:
        for label in LABELS:
            for kind in KINDS:
                assert hashlib.sha256(archive.read(f'{label}/{kind}.exe')).hexdigest() == built[label][kind]['binary_sha256']
    return rows, samples, identities


def cell(rows, kind, name):
    output = []
    for pair in range(9):
        a, b = [rows[kind, pair, label] for label in LABELS]
        medians = {label: statistics.median(rows[kind, pair, label]['cells'][name]) for label in LABELS}
        row = {'global_pair': pair, 'phase': a['phase'], 'first_process': a['order'][0],
               'seed': a['seed'], 'samples': {label: rows[kind, pair, label]['cells'][name] for label in LABELS},
               'process_medians': medians, 'ratio': medians['candidate'] / medians['baseline']}
        if kind == 'churn':
            # The retained fixture rotates [exact, hash-only, Ropey, Ropey-rehash]
            # by seed%4. Invocation is zero in every original and repeated run.
            arms = ['hashrope_exact', 'hashrope_candidate', 'ropey_exact', 'ropey_rehash_candidate']
            at = a['seed'] % 4
            row['within_process_arm_order'] = arms[at:] + arms[:at]
        output.append(row)
    fields = ['phase', 'first_process'] + (['seed'] if kind == 'churn' else [])
    strata = {field: {str(value): stats([r['ratio'] for r in output if r[field] == value])
                      for value in dict.fromkeys(r[field] for r in output)} for field in fields}
    return {'kind': kind, 'cell': name, 'unit': 'ns' if kind == 'micro' else 'ms',
            'pairs': output, 'all_pairs': stats([r['ratio'] for r in output]), 'strata': strata}


def audit_all_summary_vectors(rows):
    summaries = read(PRIOR / 'confirmation/summary.json')
    vectors = 0
    for phase, pairs in [('supplemental_six', list(range(3, 9))), ('pooled_nine', list(range(9)))]:
        for kind, cells in summaries[phase].items():
            for name, expected in cells.items():
                medians = {label: [statistics.median(rows[kind, pair, label]['cells'][name]) for pair in pairs] for label in LABELS}
                ratios = [b / a for a, b in zip(medians['baseline'], medians['candidate'])]
                assert medians == expected['process_medians']
                assert ratios == expected['paired_ratios_candidate_over_baseline']
                assert geomean(ratios) == expected['paired_ratio_geomean']
                assert [min(ratios), max(ratios)] == expected['paired_ratio_range']
                assert sum(r > 1 for r in ratios) == expected['candidate_slower_pairs']
                for label in LABELS:
                    assert statistics.median(medians[label]) == expected['median_of_process_medians'][label]
                    assert [min(medians[label]), max(medians[label])] == expected['process_median_range'][label]
                vectors += 1
    return vectors


def main():
    rows, samples, identities = collect()
    vectors = audit_all_summary_vectors(rows)
    targets = [('micro', 'core/concat_balanced_64KiB')]
    targets += [('churn', f'1000000/false/{arm}') for arm in ['hashrope_exact', 'hashrope_candidate', 'ropey_exact', 'ropey_rehash_candidate']]
    targets += [('churn', '10000/false/hashrope_exact')]
    cells = {name: cell(rows, kind, name) for kind, name in targets}
    concat = cells['core/concat_balanced_64KiB']
    concat['per_sample_position'] = {
        'warning': 'Exploratory position summaries retain every sample. Positions mix arena growth, allocation and elapsed context; they are not replacements for the predeclared estimator.',
        'positions': [{'sample': i,
                       'arm_median_ns': {label: statistics.median([r['samples'][label][i] for r in concat['pairs']]) for label in LABELS},
                       'paired_ratio_geomean': geomean([r['samples']['candidate'][i] / r['samples']['baseline'][i] for r in concat['pairs']])}
                      for i in range(7)],
        'source_derived_capacity_hypothesis': {
            'status': 'Expected doubling trajectory, not observed capacities; confirm with separate allocator/capacity probe.',
            'nodes_after_inputs': 510, 'untimed_warmup_pushes': 3, 'timed_pushes': 14000, 'nodes_after_concat_cell': 14513,
            'assumed_capacity_after_warmup': 1024, 'predicted_reallocation_counts_by_sample': [2, 1, 0, 1, 0, 0, 0],
            'note': 'Three growth-heavy and four growth-free positions can make median-of-seven equal the largest growth-free value, conditional on their ordering. No samples are excluded.'}}
    output = {'scope': 'Descriptive existing-data audit only; no new timings, causal test, significance test or excluded samples',
              'script_sha256': sha(Path(__file__)),
              'input_sha256': {str(p.relative_to(PRIOR)): sha(p) for p in [PRIOR / 'results.json', PRIOR / 'summary.json', PRIOR / 'confirmation/results.json', PRIOR / 'confirmation/summary.json', PRIOR / 'binaries.zip', PRIOR / 'workloads/bench_rust.rs', PRIOR / 'workloads/edit_compare.rs']},
              'verified_processes': len(rows), 'verified_samples': samples, 'verified_summary_vectors': vectors,
              'raw_identities': identities, 'cells': cells,
              'limits': ['Two timing phases and repeated seeds are retained as groups, not independent workload draws.',
                         'Churn seed and within-process arm order are coupled; these strata cannot isolate input shape from arm position.',
                         'Microbenchmark has no random seed; pair index modulo3 is not a seed and is not analyzed as one.',
                         'Control variability is descriptive; it neither cancels a regression nor establishes noise.']}
    with (HERE / 'existing_data_audit.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(output, handle, indent=2)
        handle.write('\n')
    print(json.dumps({'verified_processes': len(rows), 'verified_samples': samples, 'verified_summary_vectors': vectors,
                      'cells': {name: {'all_pairs': row['all_pairs'], 'strata': row['strata']} for name, row in cells.items()}}, indent=2))


if __name__ == '__main__':
    main()
