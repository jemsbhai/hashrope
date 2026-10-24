"""Diagnostic ablations using the unchanged final Rust churn harness.

Three alternating rounds, seed 1 only. These diagnostics select code paths;
the full nine-process matrix provides final performance estimates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CONSUMERS = HERE.parent / 'consumer_exactness'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build(variant):
    project = HERE / '.cache' / ('edit_variant_' + variant)
    library = project / 'hashrope'
    label = 'baseline' if variant == 'baseline' else 'optimized'
    archive = HERE / (label + '.zip')
    frozen = HERE / '.cache' / label
    if not archive.exists():
        raise RuntimeError('Capture the selected immutable archive before rebuilding diagnostics: '+str(archive))
    with zipfile.ZipFile(archive) as z:
        # Re-extract the preserved archive so a stale diagnostic copy is not used.
        z.extractall(frozen)
    source = frozen / 'hashrope/packages/rust'
    shutil.copytree(source / 'src', library / 'src', dirs_exist_ok=True)
    for name in ['Cargo.toml', 'README.md']:
        shutil.copyfile(source / name, library / name)
    if variant == 'no_clone':
        with zipfile.ZipFile(HERE / 'baseline.zip') as z:
            (library/'src/polynomial_hash.rs').write_bytes(z.read('hashrope/packages/rust/src/polynomial_hash.rs'))
    if variant == 'grouped_fastmod':
        path = library / 'src/polynomial_hash.rs'
        text = path.read_text()
        marker = 'pub fn mersenne_mod(a: u128, p: u64) -> u64 {\n'
        assert text.count(marker) == 1
        text = text.replace(marker, marker + '''    if p == MERSENNE_61 {
        let mut v = (a >> 61) + (a & MERSENNE_61 as u128);
        v = (v >> 61) + (v & MERSENNE_61 as u128);
        let r = v as u64;
        return if r >= MERSENNE_61 { r - MERSENNE_61 } else { r };
    }
''')
        path.write_text(text)
    (project / 'src').mkdir(exist_ok=True)
    shutil.copyfile(CONSUMERS / 'edit_compare.rs', project / 'src/main.rs')
    (project / 'Cargo.toml').write_text('[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="hashrope"}\nropey="=1.6.1"\n')
    shutil.copyfile(CONSUMERS / 'edit_compare.Cargo.lock', project / 'Cargo.lock')
    result = subprocess.run(['cargo','build','--release','--offline','--locked','--manifest-path',str(project / 'Cargo.toml')],capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    binary = project / 'target/release/exact-edit-comparison.exe'
    metadata = {'variant':variant,'source_sha256':{p.name:sha(p) for p in sorted((library/'src').glob('*.rs'))},
                'harness_sha256':sha(project/'src/main.rs'),'lock_sha256':sha(project/'Cargo.lock'),
                'binary':str(binary),'binary_sha256':sha(binary)}
    binary.with_suffix('.json').write_text(json.dumps(metadata,indent=2))
    return metadata

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--variants',default='baseline,no_clone,grouped,grouped_fastmod')
    p.add_argument('--build-only',action='store_true')
    p.add_argument('--run',action='store_true')
    p.add_argument('--output',type=Path,default=HERE/'edit_variants.json')
    args = p.parse_args()
    builds = [json.loads((HERE/'.cache'/('edit_variant_'+v)/'target/release/exact-edit-comparison.json').read_text())
              if args.run else build(v) for v in args.variants.split(',')]
    if args.build_only:
        print(json.dumps(builds,indent=2))
        return
    if args.output.exists():
        raise SystemExit('Refusing to overwrite '+str(args.output))
    result = {'purpose':__doc__,'builds':builds,'runs':[]}
    for round_id in range(3):
        for meta in builds if round_id%2==0 else list(reversed(builds)):
            binary=Path(meta['binary'])
            assert sha(binary)==meta['binary_sha256']
            cells={}
            for line in subprocess.check_output([str(binary),'1','0'],text=True).splitlines():
                _,size,shared,arm,_,ms=line.split('\t')
                c=cells.setdefault((size,shared,arm),{'bytes':int(size),'shared':shared=='true','arm':arm,'samples_ms':[]})
                c['samples_ms'].append(float(ms))
            for c in cells.values():
                c['median_ms']=statistics.median(c['samples_ms'])
            result['runs'].append({'variant':meta['variant'],'round':round_id,'cells':list(cells.values())})
            print(meta['variant'],round_id,flush=True)
    args.output.write_text(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
