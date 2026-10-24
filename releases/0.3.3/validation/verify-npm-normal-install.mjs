import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {existsSync,mkdirSync,mkdtempSync,readFileSync,readdirSync,writeFileSync} from 'node:fs';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import path from 'node:path';
const root=path.resolve(process.argv[2]??'.');
const sha=b=>createHash('sha256').update(b).digest('hex');
const expected='a6e279f45a3665c6f5366a317affadf5ee4a54f41e8b6de2c2f7dc55dd29bc1c';
const npmCli=[process.env.npm_execpath,path.join(path.dirname(process.execPath),'node_modules/npm/bin/npm-cli.js'),path.join(path.dirname(process.execPath),'../lib/node_modules/npm/bin/npm-cli.js')].find(p=>p&&p.endsWith('npm-cli.js')&&existsSync(p));
assert.ok(npmCli);
const dir=mkdtempSync(path.join(root,'.release','npm-normal-install-'));
const cache=path.join(dir,'empty-cache');mkdirSync(cache);
const userconfig=path.join(dir,'empty.npmrc');writeFileSync(userconfig,'');
writeFileSync(path.join(dir,'package.json'),JSON.stringify({private:true,type:'module'}));
const args=['install','hashrope@0.3.3','--ignore-scripts','--registry=https://registry.npmjs.org/','--cache',cache,'--userconfig',userconfig,'--no-audit','--no-fund','--package-lock=false','--fetch-retries=0','--loglevel=http'];
const report={startedAt:new Date().toISOString(),node:process.version,verifierSha256:sha(readFileSync(fileURLToPath(import.meta.url))),command:['node',npmCli,...args],emptyCacheBeforeInstall:readdirSync(cache).length===0,passed:false};
try{
 const r=spawnSync(process.execPath,[npmCli,...args],{cwd:dir,encoding:'utf8',timeout:120000});
 report.install={exitCode:r.status,stdout:r.stdout,stderr:r.stderr};assert.equal(r.status,0,r.stdout+r.stderr);
 const inspection=JSON.parse(readFileSync(path.join(root,'packages/npm/evidence/tarball-inspection.json'),'utf8'));
 assert.equal(inspection.sha256,expected);
 const hex=Buffer.from(inspection.integrity.slice('sha512-'.length),'base64').toString('hex');
 const cachedTarball=path.join(cache,'_cacache/content-v2/sha512',hex.slice(0,2),hex.slice(2,4),hex.slice(4));
 assert.ok(existsSync(cachedTarball),'Downloaded tarball missing from the initially empty npm cache');
 report.downloadedTarballSha256=sha(readFileSync(cachedTarball));assert.equal(report.downloadedTarballSha256,expected);
 const installed=path.join(dir,'node_modules/hashrope');
 for(const file of inspection.files)assert.equal(sha(readFileSync(path.join(installed,file.path))),file.sha256,file.path);
 report.installedFilesVerified=inspection.files.length;
 const smoke=`import assert from 'node:assert/strict';import * as esm from 'hashrope';import {createRequire} from 'node:module';const cjs=createRequire(import.meta.url)('hashrope');for(const api of [esm,cjs]){assert.equal(api.VERSION,'0.3.3');const h=new api.PolynomialHash();const a=new api.Leaf(Uint8Array.of(0,131),h),b=new api.Leaf(Uint8Array.of(1,0),h);assert.equal(a.hashVal,b.hashVal);assert.equal(api.ropeEqual(a,b),false);assert.equal(api.ropeLcp(a,b),0n);assert.throws(()=>api.ropeConcat(null,a,new api.PolynomialHash(api.MERSENNE_127)),/Incompatible hash profile/);const p=new api.Leaf(Uint8Array.of(1,2),h),x=new api.RepeatNode(p,1n<<100n,h),y=new api.RepeatNode(new api.Leaf(Uint8Array.of(1,2,1,2),h),1n<<99n,h);assert.equal(api.ropeEqual(x,y),true);assert.equal(api.ropeRangeLcp(x,1,y,3,x.len-3n),x.len-3n);}console.log('PASS normal registry install: ESM/CJS exact, collision, profile, compressed range');`;
 writeFileSync(path.join(dir,'verify.mjs'),smoke);const result=spawnSync(process.execPath,['verify.mjs'],{cwd:dir,encoding:'utf8',timeout:30000});
 report.smoke={exitCode:result.status,stdout:result.stdout,stderr:result.stderr};assert.equal(result.status,0,result.stdout+result.stderr);
 report.passed=true;
}catch(error){report.error=String(error);throw error;}
finally{report.finishedAt=new Date().toISOString();writeFileSync(path.join(root,'releases/0.3.3/validation/npm-normal-public-install-latest.json'),JSON.stringify(report,null,2)+'\n');writeFileSync(path.join(dir,'verification.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({passed:report.passed,directory:dir,downloadedTarballSha256:report.downloadedTarballSha256}));}
