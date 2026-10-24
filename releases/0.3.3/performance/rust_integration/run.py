"""Freeze/build/run a bounded release integration study; never publish or overwrite."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import tomllib
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
LABELS = ('reference', 'candidate')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    if path.exists():
        raise RuntimeError(f'Refusing to overwrite {path}')
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def hashes(root):
    return {str(p.relative_to(root)).replace('\\', '/'): sha(p)
            for p in sorted(root.rglob('*')) if p.is_file()}


def core(label):
    return HERE / '.cache' / label / 'hashrope'


def freeze(historical):
    if (HERE / 'manifest.json').exists() or (HERE / 'candidate.zip').exists():
        raise RuntimeError('A frozen snapshot already exists')
    reference = historical / 'benchmarks/exact_optimization_round4/optimized.zip'
    expected = read(reference.with_name('optimized_manifest.json'))['archive_sha256']
    assert sha(reference) == expected
    shutil.copyfile(reference, HERE / 'reference.zip')
    with zipfile.ZipFile(HERE / 'candidate.zip', 'w', zipfile.ZIP_DEFLATED) as out:
        for relative in ['packages/rust/src', 'packages/rust/Cargo.toml', 'packages/rust/Cargo.lock']:
            source = ROOT / relative
            for file in sorted(source.rglob('*')) if source.is_dir() else [source]:
                if file.is_file():
                    out.write(file, 'hashrope/' + file.relative_to(ROOT).as_posix())
    for label in LABELS:
        with zipfile.ZipFile(HERE / f'{label}.zip') as source:
            source.extractall(core(label).parent)
    workloads = HERE / 'workloads'
    workloads.mkdir()
    for name in ['bench_rust.py', 'bench_rust.rs']:
        shutil.copyfile(historical / 'benchmarks/comparison_safety' / name, workloads / name)
    for name in ['edit_compare.py', 'edit_compare.rs', 'edit_compare.Cargo.lock']:
        shutil.copyfile(historical / 'benchmarks/consumer_exactness' / name, workloads / name)
    with zipfile.ZipFile(HERE / 'drivers.zip', 'w', zipfile.ZIP_DEFLATED) as out:
        for file in [HERE/'PLAN.md', HERE/'run.py', *sorted(workloads.iterdir())]:
            out.write(file, file.relative_to(HERE).as_posix())
    manifest = {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'reference_origin': str(reference), 'reference_archive_sha256': expected,
        'candidate_root': str(ROOT), 'reference_label': 'completed round-four optimized core',
        'candidate_label': 'Rust 0.3.3 upstream integration',
        'archives': {name: sha(HERE/name) for name in ['reference.zip','candidate.zip','drivers.zip']},
        'sources': {label: hashes(core(label)) for label in LABELS},
        'drivers': {p.relative_to(HERE).as_posix(): sha(p) for p in [HERE/'run.py',HERE/'PLAN.md',*sorted(workloads.iterdir())]},
        'protocol': {'micro_pairs':3,'micro_samples':7,'churn_pairs':3,'churn_samples':3,
                     'order': [list(LABELS),list(reversed(LABELS)),list(LABELS)],
                     'churn_seed_invocations':[[1,0],[2,0],[3,0]],'python_interleaving':False}}
    save(HERE/'manifest.json',manifest)


def verify():
    manifest = read(HERE/'manifest.json')
    for name,digest in manifest['archives'].items():
        assert sha(HERE/name)==digest, name
    for name,digest in manifest['drivers'].items():
        assert sha(HERE/name)==digest, name
    for label in LABELS:
        if not core(label).exists():
            with zipfile.ZipFile(HERE/f'{label}.zip') as source:
                source.extractall(core(label).parent)
        assert hashes(core(label))==manifest['sources'][label], label
    return manifest


def command(args, log):
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.write_text(result.stdout,encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); see {log}')
    return result.stdout


def build():
    manifest=verify()
    if (HERE/'build.json').exists(): raise RuntimeError('Build already recorded')
    output={}
    original_lock=tomllib.loads((HERE/'workloads/edit_compare.Cargo.lock').read_text())
    retained={p['name']:p for p in original_lock['package'] if p['name'] not in ['hashrope','exact-edit-comparison']}
    for label in LABELS:
        output[label]={}
        for kind in ['micro','churn']:
            project=HERE/'.cache/build'/f'{label}_{kind}'
            project.mkdir(parents=True)
            (project/'src').mkdir()
            source=HERE/'workloads'/('bench_rust.rs' if kind=='micro' else 'edit_compare.rs')
            shutil.copyfile(source,project/'src/main.rs')
            name='hashrope-comparison-benchmark' if kind=='micro' else 'exact-edit-comparison'
            cargo='[package]\nname='+json.dumps(name)+'\nversion="0.0.0"\nedition="2021"\n'
            cargo+='[features]\nexact=[]\n[dependencies]\nhashrope={path='+json.dumps((core(label)/'packages/rust').as_posix())+'}\n'
            if kind=='churn':
                cargo+='ropey="=1.6.1"\n'
                shutil.copyfile(HERE/'workloads/edit_compare.Cargo.lock',project/'Cargo.lock')
            (project/'Cargo.toml').write_text(cargo,encoding='utf-8')
            lock_args=['cargo','update','--offline','--manifest-path',str(project/'Cargo.toml')]
            if kind=='churn': lock_args+=['-p','hashrope']
            command(lock_args,HERE/f'build_{label}_{kind}_lock.log')
            if kind=='churn':
                lock=tomllib.loads((project/'Cargo.lock').read_text())
                resolved={p['name']:p for p in lock['package']}
                for dep,record in retained.items(): assert resolved[dep]==record,dep
            args=['cargo','build','--release','--offline','--locked','--manifest-path',str(project/'Cargo.toml')]
            if kind=='micro': args+=['--features','exact']
            command(args,HERE/f'build_{label}_{kind}.log')
            binary=project/'target/release'/(name+('.exe' if os.name=='nt' else ''))
            output[label][kind]={'binary':str(binary),'binary_sha256':sha(binary),
                'source_sha256':sha(project/'src/main.rs'),'cargo_manifest':cargo,
                'cargo_lock':(project/'Cargo.lock').read_text(), 'cargo_lock_sha256':sha(project/'Cargo.lock')}
    verify()
    output['environment']={'rustc':subprocess.check_output(['rustc','--version','--verbose'],text=True),
        'cargo':subprocess.check_output(['cargo','--version'],text=True).strip(),
        'python':sys.version,'platform':platform.platform(),'processor':os.environ.get('PROCESSOR_IDENTIFIER'),
        'logical_processors':os.cpu_count(),'affinity':'inherited unpinned; no process pinning',
        'manifest_sha256':sha(HERE/'manifest.json')}
    with zipfile.ZipFile(HERE/'binaries.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for label in LABELS:
            for kind in ['micro','churn']:
                path=Path(output[label][kind]['binary'])
                archive.write(path,f'{label}/{kind}.exe')
    output['binaries_archive_sha256']=sha(HERE/'binaries.zip')
    save(HERE/'build.json',output)


def parse_micro(output):
    cells={}; memory={}
    for line in output.splitlines():
        part=line.split('\t')
        if part[0]=='time':
            _,name,sample,iterations,ns=part
            row=cells.setdefault(name,{'samples_ns':[],'iterations':int(iterations)})
            assert int(sample)==len(row['samples_ns']) and int(iterations)==row['iterations']
            row['samples_ns'].append(float(ns))
        elif part[0]=='memory': memory[part[1]]=int(part[2])
        else: raise ValueError(line)
    assert len(cells)==42 and len(memory)==3
    for row in cells.values():
        assert len(row['samples_ns'])==7
        row['median_ns']=statistics.median(row['samples_ns'])
    return {'cells':cells,'memory':memory}


def parse_churn(output):
    cells={}
    for line in output.splitlines():
        tag,n,shared,arm,sample,ms=line.split('\t')
        assert int(tag) in (1,2,3)
        key=f'{n}/{shared}/{arm}'
        row=cells.setdefault(key,{'samples_ms':[]})
        assert int(sample)==len(row['samples_ms'])
        row['samples_ms'].append(float(ms))
    assert len(cells)==24
    for row in cells.values():
        assert len(row['samples_ms'])==3
        row['median_ms']=statistics.median(row['samples_ms'])
    return {'cells':cells}


def run():
    manifest=verify(); built=read(HERE/'build.json')
    assert built['environment']['manifest_sha256']==sha(HERE/'manifest.json')
    assert built['binaries_archive_sha256']==sha(HERE/'binaries.zip')
    for label in LABELS:
        for kind in ['micro','churn']:
            assert sha(Path(built[label][kind]['binary']))==built[label][kind]['binary_sha256']
    raw=HERE/'raw'; raw.mkdir()
    output={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'manifest_sha256':sha(HERE/'manifest.json'),'build_sha256':sha(HERE/'build.json'),'runs':[]}
    for kind in ['micro','churn']:
        for pair,order in enumerate(manifest['protocol']['order']):
            for label in order:
                print(f'{kind} pair={pair} {label}',flush=True)
                args=[built[label][kind]['binary']]
                if kind=='churn': args += [str(pair+1),'0']
                result=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=180)
                target=raw/f'{kind}_{pair}_{label}.tsv'
                target.write_text(result.stdout,encoding='utf-8')
                target.with_suffix('.stderr').write_text(result.stderr,encoding='utf-8')
                assert result.returncode==0,(target,result.stderr)
                data=(parse_micro if kind=='micro' else parse_churn)(result.stdout)
                output['runs'].append({'kind':kind,'pair':pair,'label':label,'order':order,
                    'raw':target.relative_to(HERE).as_posix(),'raw_sha256':sha(target),
                    'binary_sha256':built[label][kind]['binary_sha256'],**data})
    output['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    verify()
    save(HERE/'results.json',output)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=['freeze','build','run'])
    parser.add_argument('--historical',type=Path)
    args=parser.parse_args()
    if args.phase=='freeze':
        assert args.historical is not None
        freeze(args.historical)
    elif args.phase=='build': build()
    else: run()
