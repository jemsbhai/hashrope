"""Bind passing source and installed-package results to the final inspected tarball."""
from pathlib import Path
import datetime,hashlib,json
ROOT=Path(__file__).resolve().parent.parent
EVIDENCE=ROOT/'evidence'
def read(name):return json.loads((EVIDENCE/name).read_text(encoding='utf-8-sig'))
version=json.loads((ROOT/'package.json').read_text())['version']
artifact=ROOT/'artifacts'/f'hashrope-{version}.tgz'
digest=hashlib.sha256(artifact.read_bytes()).hexdigest()
inspection=read('tarball-inspection.json')
assert inspection['passed'] and inspection['sha256']==digest
for f in inspection['files']:
    assert hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()==f['sha256'],f['path']
sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for folder,pattern in [('src','*.ts'),('test','*.test.mjs')] for p in (ROOT/folder).glob(pattern)}
suites=[];consumers=[]
for runtime,label in [('v20.11.1','node20'),('v24.22.1','node24')]:
    suite=read(f'suite-{label}.json');c=suite['counts']
    assert suite['node']==runtime and suite['exitCode']==0
    assert c['tests']==c['pass']==27636 and all(c[k]==0 for k in ['fail','skipped','cancelled','todo'])
    assert suite['sourceSha256']==sources
    suites.append({'node':runtime,'counts':c})
    consumer=read(f'consumer-verification-{runtime}.json')
    assert consumer['passed'] and consumer['archiveSha256']==digest
    assert len(consumer['checks'])==14 and all(c.get('exitCode')==0 for c in consumer['checks'])
    oracles=[c for c in consumer['checks'] if 'installed-tarball published-wheel oracle parity' in c['name']]
    assert len(oracles)==2
    for check in oracles:
        c=check['summary'];assert c['tests']==c['pass']==27565 and c['fail']==c['skipped']==c['cancelled']==0
    exact=[c for c in consumer['checks'] if 'installed exact and profile regressions' in c['name']]
    assert len(exact)==2
    browser=next(c for c in consumer['checks'] if c['name']=='Real headless browser with fresh temporary profile')
    assert browser['passed']
    consumers.append({'node':runtime,'checks':len(consumer['checks']),'installedOracleTestsPerEntry':27565,'installedExactTestsPerEntry':11,'browserVersion':browser['version']})
timing=read('retained-api-timing-summary.json')
compiled={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'dist/esm').glob('*.js')}
assert compiled==timing['candidate_runtime_sha256'],'Final runtime differs from timed candidate'
preservation=read('source-preservation.json');assert preservation['upstream_unchanged'] and preservation['oracle_corpus_unchanged']
arithmetic=read('retained-arithmetic-audit.json');assert arithmetic['passed'] and all(f['unchanged'] for f in arithmetic['functions'])
report={'verifiedAt':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':True,'version':version,'pythonOracleVersion':'0.2.2','artifact':{'path':str(artifact),'sha256':digest,'bytes':artifact.stat().st_size,'files':len(inspection['files'])},'librarySuites':suites,'cleanConsumers':consumers,'sourceSha256':sources,'upstreamPreservation':preservation,'retainedArithmeticFunctions':len(arithmetic['functions']),'timedRuntimeMatchesFinalArtifact':True,'publicationAttemptedByIntegration':False,'publicationStatus':'Owned by coordinated release process; this report only validates local artifact'}
(EVIDENCE/'release-summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':True,'version':version,'artifactSha256':digest,'testsPerRuntime':27636,'consumerChecksPerRuntime':14}))
