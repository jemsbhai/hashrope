import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { existsSync, mkdtempSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const root = process.argv[2] ? path.resolve(process.argv[2]) : fileURLToPath(new URL('../', import.meta.url));
const version = '0.3.3';
const expectedSha256 = 'a6e279f45a3665c6f5366a317affadf5ee4a54f41e8b6de2c2f7dc55dd29bc1c';
const registryUrl = `https://registry.npmjs.org/hashrope/${version}`;
const sha = (data, algorithm='sha256', encoding='hex') => createHash(algorithm).update(data).digest(encoding);
const dir = mkdtempSync(path.join(root, '.release', 'npm-registry-'));
const report = {verifierSha256:sha(readFileSync(fileURLToPath(import.meta.url))),startedAt:new Date().toISOString(),node:process.version,version,registryUrl,expectedSha256,checks:[],passed:false};
const npmCli = [process.env.npm_execpath,
  path.join(path.dirname(process.execPath),'node_modules/npm/bin/npm-cli.js'),
  path.join(path.dirname(process.execPath),'../lib/node_modules/npm/bin/npm-cli.js')
].find(candidate=>candidate&&candidate.endsWith('npm-cli.js')&&existsSync(candidate));
if(!npmCli)throw new Error('Run with a Node installation that includes npm');
const smoke = `
const assert = (await import('node:assert/strict')).default;
assert.equal(api.VERSION,'0.3.3');
for(const name of ['HashProfile','ropeEqual','ropeLcp','ropeRangeEqual','ropeRangeLcp']) assert.equal(typeof api[name],'function');
const h=new api.PolynomialHash(), other=new api.PolynomialHash(api.MERSENNE_127,257n);
const leaf=values=>new api.Leaf(Uint8Array.from(values),h);
const x=leaf([0,131,42]),y=leaf([1,0,42]);
assert.equal(x.hashVal,y.hashVal);assert.equal(api.ropeEqual(x,y),false);assert.equal(api.ropeLcp(x,y),0n);
assert.equal(api.ropeRangeEqual(x,2,y,2,1),true);
assert.equal(api.ropeEqual(x,new api.Leaf(Uint8Array.of(0,131,42),other)),true);
assert.ok(Object.isFrozen(h.profile));
assert.throws(()=>api.ropeConcat(null,x,other),/Incompatible hash profile/);
assert.throws(()=>api.ropeRangeLcp(x,4,x,4,0),RangeError);
const q=1n<<100n,a=new api.RepeatNode(leaf([1,2]),q,h),b=new api.RepeatNode(leaf([1,2,1,2]),q/2n,h);
assert.equal(api.ropeEqual(a,b),true);
assert.equal(api.ropeLcp(a,b),a.len);
assert.equal(api.ropeRangeLcp(a,1,b,3,a.len-3n),a.len-3n);
let da=leaf([1,2]),db=leaf([1,2]);
for(let i=0;i<100;i++){da=new api.Internal(da,da,h);db=new api.Internal(db,db,h);}
assert.equal(api.ropeLcp(da,db),da.len);
const legacy=new api.PolynomialHash(17,2);
const la=new api.Leaf(Uint8Array.of(0,13),legacy),lb=new api.Internal(new api.Leaf(Uint8Array.of(0),legacy),new api.Leaf(Uint8Array.of(13),legacy),legacy);
assert.notEqual(la.hashVal,lb.hashVal);assert.equal(api.ropeEqual(la,lb),true);
const data=Uint8Array.from({length:65536},(_,i)=>i%251),changed=data.slice();changed[65535]^=1;
assert.equal(api.ropeLcp(new api.Leaf(data,h),new api.Leaf(changed,h)),65535n);
console.log(JSON.stringify({passed:true,version:api.VERSION,collisionRejected:true,profilesChecked:true,rangeAndCompressedOracles:true,wordScan:true}));
`;
try {
  const metadataResponse=await fetch(registryUrl,{headers:{accept:'application/json'}});
  assert.equal(metadataResponse.status,200,`Registry returned HTTP ${metadataResponse.status}; package may not be published yet`);
  const metadata=await metadataResponse.json();
  assert.equal(metadata.name,'hashrope');assert.equal(metadata.version,version);
  const tarballUrl=new URL(metadata.dist.tarball);
  assert.equal(tarballUrl.protocol,'https:');assert.equal(tarballUrl.hostname,'registry.npmjs.org');
  assert.equal(tarballUrl.pathname,`/hashrope/-/hashrope-${version}.tgz`);
  let response=await fetch(tarballUrl);
  report.canonicalTarballStatus=response.status;
  if(response.status===404) {
    // A just-published object can be hidden by a negative CDN cache entry.
    // Retry the same anonymous registry artifact with a unique cache key;
    // the pinned SHA256 and registry SHA512 checks remain mandatory below.
    const freshUrl=new URL(tarballUrl);
    freshUrl.searchParams.set('verification',Date.now().toString());
    response=await fetch(freshUrl,{headers:{'cache-control':'no-cache'}});
    report.cacheBypassDownloadUrl=String(freshUrl);
  }
  assert.equal(response.status,200,'Public registry tarball unavailable');
  const bytes=Buffer.from(await response.arrayBuffer());
  const actualSha256=sha(bytes);assert.equal(actualSha256,expectedSha256,'Published tarball differs from tested artifact');
  assert.equal(`sha512-${sha(bytes,'sha512','base64')}`,metadata.dist.integrity);
  const artifact=path.join(dir,`hashrope-${version}.tgz`);writeFileSync(artifact,bytes);
  report.registryArtifact={tarballUrl:String(tarballUrl),sha256:actualSha256,integrity:metadata.dist.integrity,bytes:bytes.length};
  for(const flavor of ['esm','cjs']) {
    const target=path.join(dir,flavor);mkdirSync(target);
    writeFileSync(path.join(target,'package.json'),JSON.stringify({private:true,type:'module'}));
    const install=spawnSync(process.execPath,[npmCli,'install',artifact,'--offline','--ignore-scripts','--no-audit','--no-fund','--package-lock=false'],{cwd:target,encoding:'utf8',timeout:60000});
    assert.equal(install.status,0,install.stdout+install.stderr);
    const prefix=flavor==='esm'?"import * as api from 'hashrope';\n":"import { createRequire } from 'node:module';\nconst api=createRequire(import.meta.url)('hashrope');\n";
    const filename=path.join(target,'verify.mjs');writeFileSync(filename,prefix+smoke);
    const result=spawnSync(process.execPath,[filename],{cwd:target,encoding:'utf8',timeout:30000});
    report.checks.push({flavor,exitCode:result.status,stdout:result.stdout,stderr:result.stderr});
    assert.equal(result.status,0,result.stdout+result.stderr);
  }
  report.passed=true;
} catch(error) {report.error=String(error);throw error;}
finally {
  report.finishedAt=new Date().toISOString();
  writeFileSync(path.join(dir,'verification.json'),JSON.stringify(report,null,2)+'\n');
  writeFileSync(path.join(root,'.release','npm-registry-verification.json'),JSON.stringify(report,null,2)+'\n');
  console.log(JSON.stringify({passed:report.passed,report:path.join(dir,'verification.json'),sha256:report.registryArtifact?.sha256}));
}
