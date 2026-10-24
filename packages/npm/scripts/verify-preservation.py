"""Check the copied upstream source without writing into the standalone checkout."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
version=json.loads((ROOT/'package.json').read_text())['version']
evidence=ROOT/'evidence'/f'release-{version}'
evidence.mkdir(parents=True,exist_ok=True)
manifest=json.loads((ROOT/'evidence/upstream-0.2.2-provenance.json').read_text(encoding='utf-8-sig'))
upstream=Path(manifest['source'])
changed=[name for name,digest in manifest['copied_files_sha256'].items() if hashlib.sha256((upstream/name).read_bytes()).hexdigest()!=digest]
assert not changed,changed
corpus=hashlib.sha256((ROOT/'test/fixtures/oracle.json').read_bytes()).hexdigest()
result={'upstream_unchanged':True,'checked_files':len(manifest['copied_files_sha256']),'oracle_corpus_unchanged':corpus==manifest['copied_files_sha256']['test/fixtures/oracle.json'],'runtime_source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'src').glob('*.ts'))}}
assert result['oracle_corpus_unchanged']
(evidence/'source-preservation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
