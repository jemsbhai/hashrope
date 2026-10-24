"""Generate differential fixtures exclusively from the SHA256-pinned PyPI wheel.

Run with python -I scripts/generate-oracle.py. The isolated interpreter and zip
import check prevent the working checkout (or an installed package) being used.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,random,sys,zipfile
PIN='710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2'
SEED=20261023
ROOT=pathlib.Path(__file__).resolve().parents[1]

def main():
 parser=argparse.ArgumentParser(); parser.add_argument('--wheel',type=pathlib.Path,default=ROOT/'oracle'/'hashrope-0.2.2-py3-none-any.whl'); parser.add_argument('--output',type=pathlib.Path,default=ROOT/'test'/'fixtures')
 args=parser.parse_args(); wheel=args.wheel.resolve(); wheel_bytes=wheel.read_bytes()
 if hashlib.sha256(wheel_bytes).hexdigest()!=PIN: raise SystemExit('wheel SHA256 mismatch: refusing to import')
 sys.path.insert(0,str(wheel))
 import hashrope as hr
 import hashrope.rope as rm
 import hashrope.polynomial_hash as pm
 import hashrope.sliding as sm
 modules={'hashrope':hr,'hashrope.rope':rm,'hashrope.polynomial_hash':pm,'hashrope.sliding':sm}
 origins={k:v.__file__.replace('\\','/') for k,v in modules.items()}
 assert all(str(wheel).replace('\\','/')+'/' in v for v in origins.values()),origins
 assert hr.__version__=='0.2.2'
 rng=random.Random(SEED); rows=[]; specs={}; profiles={
  'm61':{'prime':str(hr.MERSENNE_61),'base':'131'},
  'm127':{'prime':str(hr.MERSENNE_127),'base':'131'},
  'm127b257':{'prime':str(hr.MERSENNE_127),'base':'257'},
  'small':{'prime':'31','base':'3'},
  'nonMersenne':{'prime':'100','base':'3'},
 }
 def emit(kind,profile,args,call):
  row={'id':f'{kind}-{len(rows):06d}','kind':kind,'profile':profile,'args':args}
  try: row['expected']=call()
  except Exception as e: row['error']={'type':type(e).__name__,'message':str(e)}
  rows.append(row)
 def scalar(kind,profile,args,call): emit(kind,profile,args,lambda:str(call()))
 def leaf(data): return ['leaf',bytes(data).hex()]
 def repeat(child,q): return ['repeat',child,str(q)]
 def concat(a,b): return ['concat',a,b]
 def make(spec,h):
  if spec is None:return None
  if spec[0]=='leaf':return hr.rope_from_bytes(bytes.fromhex(spec[1]),h)
  if spec[0]=='directLeaf':return hr.Leaf(bytes.fromhex(spec[1]),h)
  if spec[0]=='repeat':return hr.rope_repeat(make(spec[1],h),int(spec[2]),h)
  if spec[0]=='directRepeat':return hr.RepeatNode(make(spec[1],h),int(spec[2]),h)
  if spec[0]=='concat':return hr.rope_concat(make(spec[1],h),make(spec[2],h),h)
  if spec[0]=='internal':return hr.Internal(make(spec[1],h),make(spec[2],h),h)
  raise AssertionError(spec)
 def shape(node):
  if node is None:return None
  if isinstance(node,hr.Leaf):return ['L',node.data.hex()]
  if isinstance(node,hr.Internal):return ['I',shape(node.left),shape(node.right)]
  return ['R',shape(node.child),str(node.reps)]
 def state(node):
  try:hr.validate_rope(node); valid=True
  except AssertionError:valid=False
  return {'len':str(hr.rope_len(node)),'hash':str(hr.rope_hash(node)),'weight':str(0 if node is None else node.weight),'height':str(hr.rope_height(node)),'bytes':hr.rope_to_bytes(node).hex() if hr.rope_len(node)<=8192 else None,'valid':valid,'shapeSha256':hashlib.sha256(json.dumps(shape(node),separators=(',',':')).encode()).hexdigest()}
 def pair(node,pos,h):return [state(n) for n in hr.rope_split(node,pos,h)]
 # Arithmetic independently covers supported profiles and signed/non-Mersenne quirks.
 for pid,profile in profiles.items():
  p=int(profile['prime']); h=hr.PolynomialHash(p,int(profile['base']))
  data_cases=[b'',b'\0',b'\xff',bytes(range(256)),bytes([0,131]),bytes([1,0]),'A\u00e9\U0001f642Z'.encode()]+[rng.randbytes(rng.randrange(1,70)) for _ in range(24)]
  for data in data_cases:scalar('hash',pid,[data.hex()],lambda data=data:h.hash(data))
  for n in [0,1,2,61,127,256,2**53+1,2**100+123,-1,-2,-127,-(2**70+1)]:scalar('power',pid,[str(n)],lambda n=n:h.power(n))
  for q in [0,1,2,3,4,5,7,12,32,131,2**53+1,2**100+7,-1,-2,-3,-8,-19]:
   for a in [0,1,2,p-1,p+1,-1,-257]:scalar('phi',pid,list(map(str,[q,a,p])),lambda q=q,a=a:hr.phi(q,a,p))
  operands=[0,1,p-1,p,p+1,p*p,2**300+17,-1,-p,-2**300-33]+[rng.randrange(-2**190,2**190) for _ in range(12)]
  for a in operands:scalar('mod',pid,list(map(str,[a,p])),lambda a=a:hr.mersenne_mod(a,p))
  for a in operands:
   for b in [0,1,-7,p+3]:scalar('mul',pid,list(map(str,[a,b,p])),lambda a=a,b=b:hr.mersenne_mul(a,b,p))
  for _ in range(36):
   a,b=rng.randrange(-p*2,p*2),rng.randrange(-p*2,p*2); n=rng.choice([-2,-1,0,1,7,100,2**70+1]); q=rng.choice([-7,-1,0,1,2,17,2**80+1]); d=rng.choice([-3,-1,0,1,3,31,2**62+1]); total=rng.choice([-9,-1,0,1,5,23,2**65+3])
   scalar('concatHash',pid,list(map(str,[a,n,b])),lambda a=a,n=n,b=b:h.hash_concat(a,n,b))
   scalar('repeatHash',pid,list(map(str,[a,d,q])),lambda a=a,d=d,q=q:h.hash_repeat(a,d,q))
   scalar('overlapHash',pid,list(map(str,[a,d,total,b])),lambda a=a,d=d,total=total,b=b:h.hash_overlap(a,d,total,b))
 # Named rope descriptions are independent test inputs, never JS outputs.
 specs.update(empty=None,emptyBytes=leaf(b''),one=leaf(b'a'),bytes=leaf(bytes(range(24))),unicode=leaf('A\u00e9\U0001f642Z'.encode()),collisionA=leaf([0,131]),collisionB=leaf([1,0]),joined=concat(leaf(b'hello '),leaf(b'world')),repeat=repeat(leaf(b'ab'),7),nested=repeat(repeat(leaf(b'xyz'),3),4),repeatInternal=repeat(concat(leaf(b'abc'),leaf(b'defgh')),5),huge=repeat(leaf(b'ab'),2**100+17),hugeNested=repeat(repeat(concat(leaf(b'ab'),leaf(b'cd')),2**60+1),2**53+3))
 for direction in ['append','prepend']:
  s=None
  for i in range(1,65):
   l=leaf(bytes([i%256])*((i%4)+1));s=concat(s,l) if direction=='append' else concat(l,s)
   if i in [3,7,17,32,64]:specs[f'{direction}{i}']=s
 for i in range(18):
  s=None
  for j in range(rng.randint(2,8)):
   part=repeat(leaf(rng.randbytes(rng.randint(1,5))),rng.randint(1,14))
   s=concat(s,part) if rng.randrange(2) else concat(part,s)
  specs[f'random{i}']=s
 specs['repeatJoinLeft']=concat(repeat(concat(leaf(b'a'),leaf(b'bc')),127),leaf(b'!'))
 specs['repeatJoinRight']=concat(leaf(b'!'),repeat(concat(leaf(b'a'),leaf(b'bc')),127))
 specs['unbalanced']=['internal',leaf(b'a'),repeat(leaf(b'z'),12)]
 for pid in ['m61','m127','m127b257','small']:
  profile=profiles[pid];h=hr.PolynomialHash(int(profile['prime']),int(profile['base']))
  for name,spec in specs.items():
   node=make(spec,h); n=hr.rope_len(node)
   emit('rope',pid,[name],lambda node=node:state(node))
   if name=='unbalanced':
    emit('validate',pid,[name],lambda node=node:hr.validate_rope(node));continue
   positions=sorted(set([-100,-1,0,1,n//2,max(0,n-1),n,n+1,2**100+123]+list(range(min(n,25)+1))))
   for pos in positions:emit('split',pid,[name,str(pos)],lambda node=node,pos=pos:pair(node,pos,h))
   ranges={(0,0),(0,n),(0,n+1),(-1,1),(-2,3),(-20,4),(n,1),(n+5,3),(1,-1),(2,-3),(0,-1),(n//2,0)}
   if n<400:
    for _ in range(40):
     start=rng.randrange(n+1);length=rng.randrange(n-start+1);ranges.add((start,length))
    for start in range(min(n,9)+1):
     for length in range(min(n-start,9)+1):ranges.add((start,length))
   else:ranges.update({(1,n-2),(n//2,100),(n-10,10),(2**53+1,1000)})
   for start,length in sorted(ranges):scalar('substr',pid,[name,str(start),str(length)],lambda node=node,start=start,length=length:hr.rope_substr_hash(node,start,length,h))
 # Complete observable streaming snapshots for literal, overlapping and nonoverlap edits.
 def stream_state(sw):return {'hPrefix':str(sw.h_prefix),'lPrefix':str(sw.l_prefix),'pos':str(sw.pos),'windowLen':str(sw.window_len),'W':str(sw.W),'dMax':str(sw.d_max),'mMax':str(sw.m_max),'currentHash':str(sw.current_hash()),'finalHash':str(sw.final_hash()),'window':state(sw.r_window)}
 streams={}
 for pid in ['m61','m127','m127b257','small']:
  profile=profiles[pid];sw=hr.SlidingWindow(16,7,int(profile['prime']),int(profile['base']));steps=[]
  for i in range(130):
   if sw.window_len==0 or rng.randrange(3)==0:op=['bytes',rng.randbytes(rng.randrange(0,40)).hex()];sw.append_bytes(bytes.fromhex(op[1]))
   else:op=['copy',str(rng.randint(1,sw.window_len)),str(rng.randrange(0,65))];sw.append_copy(int(op[1]),int(op[2]))
   steps.append(op);emit('sliding',pid,[pid,i],lambda sw=sw:stream_state(sw))
  streams[pid]={'dMax':'16','mMax':'7','steps':steps}
 # Edge traces retain reference malformed-copy behavior; snapshots after errors confirm no changes.
 for name,steps in {
  'threshold':[['bytes',b'abcd'.hex()],['bytes',b'ef'.hex()],['bytes',b'g'.hex()],['copy','2','10']],
  'quirks':[['bytes',b'abc'.hex()],['copy','1','-1'],['copy','-1','-2'],['copy','0','0']],
  'errors':[['copy','1','0'],['copy','0','1'],['bytes',b'abcdefghi'.hex()],['copy','5','0'],['copy','10','1'],['copy','0','3']],
  'huge':[['bytes',b'ab'.hex()],['copy','2',str(2**100+3)]],
 }.items():
  sw=hr.SlidingWindow(4,2);trace='edge-'+name
  for i,op in enumerate(steps):
   error=None
   try:
    if op[0]=='bytes':sw.append_bytes(bytes.fromhex(op[1]))
    else:sw.append_copy(int(op[1]),int(op[2]))
   except Exception as e:error={'type':type(e).__name__,'message':str(e)}
   emit('slidingEdge','m61',[trace,i],lambda sw=sw,error=error:{'operationError':error,'state':stream_state(sw)})
  streams[trace]={'dMax':'4','mMax':'2','steps':steps}
 # Explicit constructor and helper errors, without invoking nonterminating p<=0 folding.
 for p,b in [(31,1),(31,31),(31,-2),(0,2),(-7,2)]:emit('constructorError','m61',[str(p),str(b)],lambda p=p,b=b:hr.PolynomialHash(p,b))
 for d,m in [(0,0),(0,1),(1,0),(-1,10),(1,-5)]:emit('windowConstructorError','m61',[str(d),str(m)],lambda d=d,m=m:hr.SlidingWindow(d,m))
 error_specs={'emptyLeaf':['directLeaf',''],'nullInternal':['internal',None,leaf(b'a')],'nullRepeat':['directRepeat',None,'2'],'repeatOne':['directRepeat',leaf(b'a'),'1'],'repeatNegative':['directRepeat',leaf(b'a'),'-4'],'repeatHelperNegative':repeat(leaf(b'a'),-1),'repeatHelperNull':repeat(None,2)}
 specs.update(error_specs)
 for name,spec in error_specs.items():emit('nodeError','m61',[name],lambda spec=spec:state(make(spec,hr.PolynomialHash())))
 corpus={'schemaVersion':1,'seed':SEED,'profiles':profiles,'ropes':specs,'streams':streams,'rows':rows}
 data=(json.dumps(corpus,ensure_ascii=True,separators=(',',':'))+'\n').encode()
 with zipfile.ZipFile(wheel) as archive:
  source_hashes={name:hashlib.sha256(archive.read(name)).hexdigest() for name in archive.namelist() if name.startswith('hashrope/') and name.endswith('.py')}
  license_bytes=archive.read('hashrope-0.2.2.dist-info/licenses/LICENSE')
 members={}
 for cls in [hr.PolynomialHash,hr.Leaf,hr.Internal,hr.RepeatNode,hr.SlidingWindow]:
  members[cls.__name__]=sorted(k for k in vars(cls) if not k.startswith('_'))
 manifest={'schemaVersion':1,'pythonVersion':hr.__version__,'wheel':{'filename':wheel.name,'sha256':PIN,'source':'https://pypi.org/project/hashrope/0.2.2/'},'importOrigins':origins,'sourceSha256':source_hashes,'pythonExports':hr.__all__,'classPublicMembers':members,'balanceConstants':{'ALPHA_NUM':rm.ALPHA_NUM,'ALPHA_DEN':rm.ALPHA_DEN},'seed':SEED,'corpus':{'filename':'oracle.json','sha256':hashlib.sha256(data).hexdigest(),'records':len(rows),'kinds':dict(sorted(collections.Counter(r['kind'] for r in rows).items())),'errors':sum('error' in r for r in rows)},'exclusions':['Nonpositive-modulus folding can fail to terminate and is intentionally not invoked.','Python-specific duck types, object hashes/repr, reflection and resource exhaustion are not numerical oracle contracts.'],'licenseSha256':hashlib.sha256(license_bytes).hexdigest()}
 args.output.mkdir(parents=True,exist_ok=True);(args.output/'oracle.json').write_bytes(data);(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'records':len(rows),'errors':manifest['corpus']['errors'],'kinds':manifest['corpus']['kinds'],'sha256':manifest['corpus']['sha256']}))
if __name__=='__main__':main()
