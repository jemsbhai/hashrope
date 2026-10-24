"""Render every final trial cell; this performs no new timing or inference."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
data_path = HERE / 'final_comparison/summary.json'
data = json.loads(data_path.read_text(encoding='utf-8'))
lines = [
    '# Complete final trial results', '',
    '**The candidate was rejected; production retains the baseline.** See the',
    '[study report](../README.md) and [independent audit](../final_review.md).', '',
    'All 48 processes, 8,784 samples and 66 cells are retained. Each ratio is the',
    'geometric mean of twelve paired process-median ratios; intervals are',
    'exploratory 95% log-t intervals, unadjusted across cells. Time columns are',
    'medians of process medians and need not yield the same ratio. Slower counts',
    'refer to the rejected candidate. Order columns identify which implementation',
    'ran first. Full individual ratios, seed strata and sample vectors remain in',
    '[summary.json](summary.json) and [results.json](results.json).', '',
    'Baseline is frozen round six; candidate is the direction-argument trial.',
    'The source-match flags in the summary precede the documented restoration.',
    'No cells or unfavorable samples are removed. Hash-only rows are historical',
    'candidate controls, not substitutes for exact byte comparison.', '',
    'Summary SHA-256: `' + hashlib.sha256(data_path.read_bytes()).hexdigest() + '`.', '',
]
for kind, rows in data['cells'].items():
    lines += ['## ' + ('Microbenchmark' if kind == 'micro' else 'Edit/comparison workload'), '',
              '| Cell | Baseline | Rejected candidate | Ratio [interval] | Slower | Ratio range | Baseline-first / candidate-first |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for name, row in rows.items():
        med = row['median_of_process_medians']
        factor = 1 if row['unit'] == 'ns' else 1000
        unit = 'ns' if row['unit'] == 'ns' else 'µs'
        ratio = row['paired_ratio_geomean']
        low, high = row['exploratory_95_log_t_interval']
        rlow, rhigh = row['paired_ratio_range']
        orders = row['order_geomeans']
        lines.append(f"| `{name}` | {med['baseline'] * factor:.4f} {unit} | "
                     f"{med['candidate'] * factor:.4f} {unit} | {ratio:.4f} [{low:.4f}, {high:.4f}] | "
                     f"{row['candidate_slower_pairs']}/12 | {rlow:.4f}–{rhigh:.4f} | "
                     f"{orders['baseline']:.4f} / {orders['candidate']:.4f} |")
    lines.append('')
lines += ['## Layout', '',
          'Both arms report a 48-byte `NodeInner`, 2,136-byte inline `Arena`,',
          'and 255 nodes for the 64 KiB construction. No layout change occurred.', '']
(HERE / 'final_comparison/README.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
