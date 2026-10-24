"""Summarize three fresh diagnostic processes without treating batches as replicates."""
from pathlib import Path
import json,statistics,hashlib
ROOT=Path(__file__).resolve().parent.parent
runs=[json.loads((ROOT/'evidence'/f'retained-api-timing-final-{i}.json').read_text()) for i in range(3)]
for key in ['harness_sha256','baseline_runtime_sha256','candidate_runtime_sha256','protocol']:
    assert all(r[key]==runs[0][key] for r in runs)
rows=[]
for first in runs[0]['summary']:
    name=first['fixture'];r=[next(x for x in run['summary'] if x['fixture']==name) for run in runs]
    b=[x['baseline_median_ns'] for x in r];c=[x['candidate_median_ns'] for x in r]
    rows.append({'fixture':name,'baseline_median_ns':statistics.median(b),'candidate_median_ns':statistics.median(c),'baseline_process_range_ns':[min(b),max(b)],'candidate_process_range_ns':[min(c),max(c)],'baseline_all_round_range_ns':[min(x['baseline_round_range_ns'][0] for x in r),max(x['baseline_round_range_ns'][1] for x in r)],'candidate_all_round_range_ns':[min(x['candidate_round_range_ns'][0] for x in r),max(x['candidate_round_range_ns'][1] for x in r)],'candidate_over_baseline':statistics.median(c)/statistics.median(b),'paired_process_ratios':[y/x for x,y in zip(b,c)]})
result={'statistic':'median of three process medians; each process is median of six within-process round medians, each from five batches. Min/max process and all-round ranges retained. No confidence/significance claim.','protocol':runs[0]['protocol'],'node':runs[0]['node'],'cpu':runs[0]['cpu'],'harness_sha256':runs[0]['harness_sha256'],'baseline_runtime_sha256':runs[0]['baseline_runtime_sha256'],'candidate_runtime_sha256':runs[0]['candidate_runtime_sha256'],'rows':rows,'compressed_completion_by_process':[r['compressed_completion'] for r in runs]}
(ROOT/'evidence/retained-api-timing-summary.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# npm retained-API timing diagnostic','',result['statistic'],'',result['protocol'],'',f"Runtime: {result['node']}; CPU: {result['cpu']}. Baseline is the verified published 0.2.2 tarball; candidate is frozen 0.3.3 runtime code. Units below are microseconds per operation.",'','| Fixture | 0.2.2 median [process range] | 0.3.3 median [process range] | Candidate/baseline time |','|---|---:|---:|---:|']
for x in rows:
    def fmt(prefix):
        lo,hi=x[f'{prefix}_process_range_ns'];return f"{x[f'{prefix}_median_ns']/1000:.3f} [{lo/1000:.3f}, {hi/1000:.3f}]"
    lines.append(f"| {x['fixture']} | {fmt('baseline')} | {fmt('candidate')} | {x['candidate_over_baseline']:.3f}x |")
lines+=['','These bounded, warmed fixtures are diagnostic evidence, not broad package performance parity. Setup is excluded and there is no affinity/frequency control. Tiny operations are particularly sensitive to call/measurement overhead. The independent-equality candidate initially regressed (~2x) with a scalar leaf helper; the initial raw run and archive are retained separately. A portable DataView word scan was selected only after that observation and byte-oracle review. Unchanged arithmetic timings also vary, so no causal interpretation is assigned to their small movements. Profile safety adds checks and stored profile references; zero overhead is not claimed.','', 'Huge Repeat and independent DAG cases each completed with the expected exact result in all three processes. They are single completion observations per case, not comparable throughput samples; the old implementation was not asked to traverse astronomical expanded content. Exact comparisons still have general traversal, memoization, BigInt and unshared-byte costs.','', 'Raw samples: `evidence/retained-api-timing-final-{0,1,2}.json`; frozen runtime and harness digests are embedded in every file, with archive provenance in `evidence/retained-api-timing-manifest-final.json`. `evidence/retained-api-timing-initial.json` preserves the rejected scalar-helper diagnostic. Final documentation-only packaging changes are outside the measured runtime source scope.']
(ROOT/'evidence/RETAINED_API_TIMING.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
