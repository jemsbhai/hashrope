"""Reproduce the standalone hash-kernel diagnostic; never replace prior output."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=HERE/'hash_kernel_reproduction.json')
    p.add_argument('--existing-raw',type=Path)
    args=p.parse_args()
    if args.output.exists():
        raise SystemExit('Refusing to overwrite '+str(args.output))
    src=HERE/'hash_kernel_profile.rs'
    binary=HERE/'.cache/hash_kernel_profile.exe'
    command=['rustc','-O',str(src),'-o',str(binary)]
    if args.existing_raw:
        raw=args.existing_raw.read_text()
    else:
        binary.parent.mkdir(exist_ok=True)
        subprocess.run(command,check=True)
        raw=subprocess.check_output([str(binary)],text=True)
    rows=[json.loads(line) for line in raw.splitlines()]
    assert len(rows)==13*6*7
    result={'purpose':__doc__,'build_command':command,
            'rustc':subprocess.check_output(['rustc','--version'],text=True).strip(),
            'source_sha256':sha(src),'binary_sha256':sha(binary),
            'existing_raw_sha256':sha(args.existing_raw) if args.existing_raw else None,
            'method':'One process, fixed kernel order, seven samples/cell. Diagnostic only; no affinity or frequency locking.',
            'rows':rows}
    args.output.write_text(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
