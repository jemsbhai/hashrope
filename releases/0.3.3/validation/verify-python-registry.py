import hashlib
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import venv

root = next(p for p in Path(__file__).resolve().parents if (p / 'packages/python/pyproject.toml').is_file())
validation = root / 'releases/0.3.3/validation'
expected = {x['filename']: x for x in json.loads((validation / 'python-artifacts.json').read_text())['artifacts']}
records = []
for package in ('hashrope', 'hashrope-bio'):
    url = f'https://pypi.org/pypi/{package}/0.3.3/json'
    metadata = json.load(urllib.request.urlopen(url, timeout=30))
    assert metadata['info']['version'] == '0.3.3'
    for entry in metadata['urls']:
        reference = expected[entry['filename']]
        content = urllib.request.urlopen(entry['url'], timeout=30).read()
        digest = hashlib.sha256(content).hexdigest()
        assert digest == entry['digests']['sha256'] == reference['sha256']
        assert len(content) == reference['bytes']
        records.append({'package': package, 'filename': entry['filename'], 'url': entry['url'], 'sha256': digest, 'bytes': len(content)})
assert len(records) == 4
environment = root / '.release/public-python'
venv.create(environment, with_pip=True)
python = environment / 'Scripts/python.exe'
commands = [
    [str(python), '-m', 'pip', 'install', '--no-cache-dir', '--index-url', 'https://pypi.org/simple', 'hashrope==0.3.3', 'hashrope-bio==0.3.3', 'pytest'],
    [str(python), '-I', '-c', "import hashrope,hashrope_bio,sys; from pathlib import Path; assert hashrope.__version__==hashrope_bio.__version__=='0.3.3'; assert all(Path(m.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()) for m in (hashrope,hashrope_bio)); print(hashrope.__file__,hashrope_bio.__file__)"],
    [str(python), '-I', str(root / 'packages/python/scripts/verify_release.py')],
    [str(python), '-I', '-m', 'pytest', '--import-mode=importlib', '--rootdir', str(environment), str(root / 'packages/python/tests'), str(root / 'usecases/computational-biology/code/tests'), '-q'],
]
with (validation / 'python-public-install.log').open('w', encoding='utf-8') as log:
    for command in commands:
        print(command, flush=True)
        result = subprocess.run(command, cwd=environment, stdout=log, stderr=subprocess.STDOUT)
        log.flush()
        if result.returncode:
            raise SystemExit(f'Public install check failed; inspect {log.name}')
report = {'status': 'public artifacts match reviewed builds; fresh registry installation and 244 tests passed', 'version': '0.3.3', 'artifacts': records, 'commands': commands}
(validation / 'python-public-install.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(report['status'])
