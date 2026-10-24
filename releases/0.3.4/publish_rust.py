"""Publish a reviewed 0.3.4 crate only after the coordinated release gate."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('package', choices=['hashrope', 'hashrope-bio'])
    parser.add_argument('--release-gate-passed', action='store_true', required=True)
    args = parser.parse_args()
    core = args.package == 'hashrope'
    crate = ROOT / ('packages/rust' if core else 'usecases/computational-biology/code/rust')
    directory = HERE / 'validation' / ('rust-core-publication' if core else 'rust-bio-publication')
    directory.mkdir()
    assert 'version = "0.3.4"' in (crate / 'Cargo.toml').read_text()
    status = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no', '--', '.'], cwd=crate, text=True)
    assert not status, status
    record = {'package': args.package, 'version': '0.3.4',
              'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=crate, text=True).strip(),
              'commands': [], 'published': False}

    def command(name, arguments):
        with (directory / f'{name}.log').open('x', encoding='utf-8', newline='\n') as log:
            result = subprocess.run(arguments, cwd=crate, stdout=log, stderr=subprocess.STDOUT, text=True)
        record['commands'].append({'arguments': arguments, 'returncode': result.returncode,
                                   'log': name + '.log', 'sha256': sha(directory / f'{name}.log')})
        print(f'{args.package} {name}: exit {result.returncode}', flush=True)
        if result.returncode:
            raise RuntimeError(f'{name} failed; inspect {directory / (name + ".log")} before retrying')

    try:
        command('package', ['cargo', 'package', '--all-features'])
        command('dry-run', ['cargo', 'publish', '--dry-run'])
        archive = crate / 'target/package' / f'{args.package}-0.3.4.crate'
        reviewed = directory / 'reviewed.crate'
        shutil.copyfile(archive, reviewed)
        record['reviewed_archive_sha256'] = sha(reviewed)
        with tarfile.open(reviewed) as packed:
            vcs = json.load(packed.extractfile(f'{args.package}-0.3.4/.cargo_vcs_info.json'))
            assert vcs['git']['sha1'] == record['source_commit'] and not vcs['git'].get('dirty', False), vcs
        if core:
            expected = json.loads((HERE / 'validation/rust-clean-package.json').read_text())
            assert record['reviewed_archive_sha256'] == expected['sha256'], 'Reviewed core artifact changed'
        command('publish', ['cargo', 'publish'])
        record['published'] = True
        record['upload_archive_sha256'] = sha(archive)
        assert record['upload_archive_sha256'] == record['reviewed_archive_sha256'], 'Upload archive changed'
        url = f'https://static.crates.io/crates/{args.package}/{args.package}-0.3.4.crate'
        request = urllib.request.Request(url, headers={'User-Agent': 'hashrope-release-verification/0.3.4'})
        with urllib.request.urlopen(request, timeout=60) as response:
            public = response.read()
        (directory / 'public.crate').write_bytes(public)
        record['public_url'] = url
        record['public_archive_sha256'] = sha(directory / 'public.crate')
        assert record['public_archive_sha256'] == record['reviewed_archive_sha256'], 'Public artifact differs'
        record['public_archive_matches_reviewed'] = True
        print(f'{args.package} public artifact matches: {record["public_archive_sha256"]}', flush=True)
    finally:
        record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (directory / 'results.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
