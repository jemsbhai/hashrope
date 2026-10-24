import { readFileSync, writeFileSync, readdirSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { pathToFileURL, fileURLToPath } from 'node:url';
import path from 'node:path';
import os from 'node:os';
const [baselineDir, candidateDir, output] = process.argv.slice(2);
if (!output) throw new Error('Usage: node scripts/benchmark-release.mjs <baseline-package-dir> <candidate-package-dir> <output.json>');
const baseline = await import(pathToFileURL(path.resolve(baselineDir, 'dist/esm/index.js')).href);
const candidate = await import(pathToFileURL(path.resolve(candidateDir, 'dist/esm/index.js')).href);
const sha = s => createHash('sha256').update(s).digest('hex');
const median = xs => { const x = [...xs].sort((a,b)=>a-b); const m=x.length>>1; return x.length%2 ? x[m] : (x[m-1]+x[m])/2; };
const bytes = Uint8Array.from({length:65536},(_,i)=>i%251);
const changed = bytes.slice(); changed[changed.length-1]^=1;
const early = bytes.slice(); early[0]^=1;
function fixtures(api) {
  const h=new api.PolynomialHash();
  const big=new api.Leaf(bytes,h), independent=new api.Leaf(bytes,h), late=new api.Leaf(changed,h), first=new api.Leaf(early,h);
  const small=new api.Leaf(bytes.subarray(0,128),h), other=new api.Leaf(bytes.subarray(128,256),h);
  const branchA=new api.Internal(big,new api.Leaf(Uint8Array.of(1),h),h), branchB=new api.Internal(big,new api.Leaf(Uint8Array.of(2),h),h);
  const cases = [
    ['hash_4096B',20,()=>h.hash(bytes.subarray(0,4096))],
    ['construct_leaf_128B',100,()=>new api.Leaf(bytes.subarray(0,128),h).hashVal],
    ['concat_two_128B_warm',1000,()=>api.ropeConcat(small,other,h).hashVal],
    ['split_64KiB_leaf_middle',2,()=>api.ropeSplit(big,32768,h)[0].hashVal],
    ['substring_64B_in_64KiB_leaf',100,()=>api.ropeSubstrHash(big,123,64,h)],
    ['repeat_metadata_100bit',100,()=>api.ropeRepeat(small,1n<<100n,h).hashVal],
    ['equal_independent_64KiB',100,()=>api.contentEquals(big,independent)],
    ['equal_late_mismatch_64KiB',100,()=>api.contentEquals(big,late)],
    ['equal_early_mismatch_64KiB',100,()=>api.contentEquals(big,first)],
    ['equal_shared_prefix_64KiB',100,()=>api.contentEquals(branchA,branchB)],
    ['equal_same_root',100000,()=>api.contentEquals(big,big)],
  ];
  return cases.map(([name,repetitions,run])=>({name,repetitions,run}));
}
const arms={baseline:fixtures(baseline),candidate:fixtures(candidate)};
for(let i=0;i<arms.baseline.length;i++) {
  if(arms.baseline[i].run()!==arms.candidate[i].run()) throw new Error(`Retained operation mismatch: ${arms.baseline[i].name}`);
  for(const arm of Object.values(arms)) for(let j=0;j<Math.min(20,arm[i].repetitions);j++) arm[i].run();
}
const raw=[];
let sink;
for(let round=0;round<6;round++) {
  const indexes=Array.from({length:arms.baseline.length},(_,i)=>(i+round)%arms.baseline.length);
  if(round%2) indexes.reverse();
  for(const index of indexes) {
    const order=(round+index)%2 ? ['candidate','baseline'] : ['baseline','candidate'];
    for(const arm of order) {
      const fixture=arms[arm][index], samples=[];
      for(let batch=0;batch<5;batch++) {
        const start=process.hrtime.bigint();
        for(let rep=0;rep<fixture.repetitions;rep++) sink=fixture.run();
        samples.push(Number(process.hrtime.bigint()-start)/fixture.repetitions);
      }
      raw.push({round,arm,fixture:fixture.name,order,repetitions:fixture.repetitions,ns_per_call:samples,median_ns:median(samples)});
    }
  }
}
const summary=arms.baseline.map(({name})=>{
  const b=raw.filter(x=>x.fixture===name&&x.arm==='baseline').map(x=>x.median_ns);
  const c=raw.filter(x=>x.fixture===name&&x.arm==='candidate').map(x=>x.median_ns);
  return {fixture:name,baseline_median_ns:median(b),candidate_median_ns:median(c),baseline_round_range_ns:[Math.min(...b),Math.max(...b)],candidate_round_range_ns:[Math.min(...c),Math.max(...c)],candidate_over_baseline:median(c)/median(b)};
});
const h=new candidate.PolynomialHash(), q=1n<<100n;
const a=new candidate.RepeatNode(new candidate.Leaf(Uint8Array.of(1,2),h),q,h);
const b=new candidate.RepeatNode(new candidate.Leaf(Uint8Array.of(1,2,1,2),h),q/2n,h);
let da=new candidate.Leaf(Uint8Array.of(1,2),h), db=new candidate.Leaf(Uint8Array.of(1,2),h);
for(let i=0;i<100;i++){da=new candidate.Internal(da,da,h);db=new candidate.Internal(db,db,h);}
const compressed=[];
for(const [name,run,expected] of [
  ['equal_repeat_different_period_101bit',()=>candidate.ropeEqual(a,b),true],
  ['lcp_repeat_shifted_range_101bit',()=>candidate.ropeRangeLcp(a,1n,b,3n,a.len-3n),a.len-3n],
  ['lcp_independent_dag_101bit',()=>candidate.ropeLcp(da,db),da.len],
]) {
  const start=process.hrtime.bigint(), actual=run(), elapsed=Number(process.hrtime.bigint()-start);
  if(actual!==expected) throw new Error(`Compressed oracle failed: ${name}`);
  compressed.push({name,elapsed_ns:elapsed,expected:String(expected),passed:true});
}
function hashes(dir) {
  return Object.fromEntries(readdirSync(path.join(dir,'dist/esm')).filter(f=>f.endsWith('.js')).map(f=>[f,sha(readFileSync(path.join(dir,'dist/esm',f)))]));
}
const result={protocol:'six within-process rounds; fixture order rotated/reversed; each fixture has three AB and three BA arm orders; five batches per arm/round; equal repetitions; setup excluded; bounded warmup; no affinity/frequency control',statistic:'median of six round medians, with min/max round medians; diagnostic only, not independent processes or confidence intervals',node:process.version,platform:process.platform,architecture:process.arch,cpu:os.cpus()[0]?.model,logical_cpus:os.cpus().length,baseline_version:baseline.VERSION,candidate_version:candidate.VERSION,harness_sha256:sha(readFileSync(fileURLToPath(import.meta.url))),baseline_runtime_sha256:hashes(baselineDir),candidate_runtime_sha256:hashes(candidateDir),raw,summary,compressed_completion:compressed,compressed_note:'one untimed-setup completion per case; no old-engine huge run because its traversal can expand logical content; these are finite completion observations, not comparative throughput or universal bounds',sink_type:typeof sink};
writeFileSync(output,JSON.stringify(result,null,2)+'\n');
for(const row of summary) console.log(`${row.fixture}: ${(row.baseline_median_ns/1000).toFixed(3)} -> ${(row.candidate_median_ns/1000).toFixed(3)} us (${row.candidate_over_baseline.toFixed(3)}x time)`);
