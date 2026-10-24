"""Require a metadata-only 0.3.4 runtime against the verified 0.3.3 tarball."""
from pathlib import Path
import hashlib
import json
import tarfile

ROOT = Path(__file__).resolve().parent.parent
VERSION = json.loads((ROOT / 'package.json').read_text())['version']
assert VERSION == '0.3.4', 'Review the previous release and permitted changes before reuse'
EVIDENCE = ROOT / 'evidence' / f'release-{VERSION}'
EVIDENCE.mkdir(parents=True, exist_ok=True)
BASELINE = ROOT / 'artifacts/hashrope-0.3.3.tgz'
BASELINE_SHA256 = 'a6e279f45a3665c6f5366a317affadf5ee4a54f41e8b6de2c2f7dc55dd29bc1c'
sha = lambda data: hashlib.sha256(data).hexdigest()
assert sha(BASELINE.read_bytes()) == BASELINE_SHA256
historical = json.loads((ROOT / 'evidence/tarball-inspection.json').read_text())
assert historical['passed'] and historical['sha256'] == BASELINE_SHA256
with tarfile.open(BASELINE, 'r:gz') as archive:
    baseline = {
        member.name.removeprefix('package/'): archive.extractfile(member).read()
        for member in archive.getmembers()
        if member.isfile() and member.name.startswith(('package/src/', 'package/dist/'))
    }
current = {
    p.relative_to(ROOT).as_posix(): p.read_bytes()
    for folder in ['src', 'dist'] for p in (ROOT / folder).rglob('*') if p.is_file()
}
assert baseline.keys() == current.keys(), 'Runtime/source file inventory differs'
permitted = {
    'src/index.ts',
    'dist/esm/index.js', 'dist/esm/index.d.ts', 'dist/esm/index.js.map',
    'dist/cjs/index.js', 'dist/cjs/index.d.ts', 'dist/cjs/index.js.map',
}
files = []
for name, old in sorted(baseline.items()):
    new = current[name]
    identical = old == new
    if not identical:
        assert name in permitted, f'Unexpected runtime change: {name}'
        assert old.count(b'0.3.3') == 1, f'Expected one version literal: {name}'
        assert old.replace(b'0.3.3', b'0.3.4') == new, f'Non-version change: {name}'
    files.append({'path': name, 'baselineSha256': sha(old), 'currentSha256': sha(new),
                  'identical': identical, 'onlyVersionLiteralChanged': not identical})
assert {f['path'] for f in files if not f['identical']} == permitted
timing = json.loads((ROOT / 'evidence/retained-api-timing-summary.json').read_text())
assert {Path(name).name: sha(data) for name, data in baseline.items()
        if name.startswith('dist/esm/') and name.endswith('.js')} == timing['candidate_runtime_sha256']
oracle = json.loads((ROOT / 'evidence/upstream-0.2.2-provenance.json').read_text())
for name in ['test/fixtures/oracle.json', 'test/fixtures/manifest.json']:
    assert sha((ROOT / name).read_bytes()) == oracle['copied_files_sha256'][name], name
result = {'passed': True, 'version': VERSION, 'baselineVersion': '0.3.3',
          'baselineArtifact': str(BASELINE), 'baselineArtifactSha256': BASELINE_SHA256,
          'files': files, 'runtimeBehaviorUnchangedExceptVersionMetadata': True,
          'oracleCorpusAndManifestUnchanged': True, 'historicalTimedRuntimeMatchesBaseline': True,
          'newPerformanceMeasurements': False}
(EVIDENCE / 'runtime-preservation.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'runtimeFiles': len(files),
                  'identicalFiles': sum(f['identical'] for f in files), 'versionOnlyFiles': len(permitted)}))
