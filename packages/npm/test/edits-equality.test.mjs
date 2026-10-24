import test from 'node:test';
import assert from 'node:assert/strict';
import { PolynomialHash, MERSENNE_61, MERSENNE_127 } from '../dist/esm/polynomial.js';
import { Leaf, Internal, RepeatNode, ropeLen, ropeHash, ropeToBytes, validateRope } from '../dist/esm/nodes.js';
import { ropeConcat, ropeSplit, ropeRepeat, ropeSubstrHash, ropeInsert, ropeDelete, contentEquals } from '../dist/esm/rope.js';
const b=s=>new TextEncoder().encode(s),h=new PolynomialHash(),leaf=s=>new Leaf(b(s),h);

test('persistent insert/delete clamp endpoints and reject unsafe or negative dimensions',()=>{
 const a=leaf('abc'); const n=ropeInsert(a,1n,leaf('xy'),h); assert.deepEqual(ropeToBytes(n),b('axybc'));assert.deepEqual(ropeToBytes(a),b('abc'));
 assert.deepEqual(ropeToBytes(ropeDelete(n,1n,2n,h)),b('abc'));
 assert.deepEqual(ropeToBytes(ropeInsert(a,-1n,leaf('z'),h)),b('zabc'));
 assert.deepEqual(ropeToBytes(ropeInsert(a,100n,leaf('z'),h)),b('abcz'));
 assert.deepEqual(ropeToBytes(ropeDelete(a,-1n,2n,h)),b('c'));
 assert.equal(ropeDelete(a,0n,100n,h),null);assert.equal(ropeDelete(a,100n,1n,h),a);
 assert.equal(ropeInsert(null,5n,a,h),a);assert.equal(ropeDelete(null,0n,3n,h),null);
 assert.throws(()=>ropeDelete(a,0n,-1n,h),RangeError);
 assert.throws(()=>ropeInsert(a,Number.MAX_SAFE_INTEGER+1,null,h),RangeError);
 const huge=ropeRepeat(a,1n<<100n,h), edited=ropeInsert(huge,1n<<99n,leaf('z'),h);
 assert.equal(ropeLen(edited),ropeLen(huge)+1n); assert.equal(ropeLen(ropeDelete(edited,1n<<99n,1n,h)),ropeLen(huge));validateRope(edited);
});

test('exact equality resists collisions and accepts equal bytes across shapes and profiles',()=>{
 for(const prime of [MERSENNE_61,MERSENNE_127]) {
  const hash=new PolynomialHash(prime,131n),a=new Leaf(Uint8Array.of(0,131),hash),c=new Leaf(Uint8Array.of(1,0),hash);
  assert.equal(a.hashVal,c.hashVal);assert.equal(contentEquals(a,c),false);
 }
 assert.equal(contentEquals(null,null),true);assert.equal(contentEquals(null,leaf('a')),false);
 assert.equal(contentEquals(leaf('abc'),new Internal(leaf('a'),leaf('bc'),h)),true);
 assert.equal(contentEquals(leaf('ababab'),new RepeatNode(leaf('ab'),3n,h)),true);
 const other=new PolynomialHash(MERSENNE_127,257n); assert.equal(contentEquals(leaf('abc'),new Leaf(b('abc'),other)),true);
 assert.equal(contentEquals(leaf('abcdef'),leaf('abcdeg')),false);
 const huge=ropeRepeat(leaf('ab'),1n<<100n,h); assert.equal(contentEquals(huge,huge),true);
 assert.equal(contentEquals(huge,new RepeatNode(leaf('ab'),1n<<100n,h)),true);
});

test('seeded persistent edit traces match independent byte and polynomial models',()=>{
 let seed=0x10232026; const rnd=n=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed%n};
 for(const [prime,base] of [[MERSENNE_61,131n],[MERSENNE_127,257n],[31n,3n]]) {
  const hash=new PolynomialHash(prime,base); let root=null,data=new Uint8Array(); const snapshots=[];
  for(let step=0;step<240;step++) {
   const at=rnd(data.length+1);
   if(rnd(3)<2) {
    const value=Uint8Array.from({length:rnd(8)+1},()=>rnd(256));
    root=ropeInsert(root,at,new Leaf(value,hash),hash);data=Uint8Array.from([...data.subarray(0,at),...value,...data.subarray(at)]);
   } else { const len=rnd(12);root=ropeDelete(root,at,len,hash);data=Uint8Array.from([...data.subarray(0,at),...data.subarray(at+len)]); }
   let expected=0n; for(const byte of data) expected=(expected*base+BigInt(byte+1))%prime;
   assert.equal(ropeHash(root),expected);assert.deepEqual(ropeToBytes(root),data);validateRope(root);
   const start=rnd(data.length+1),length=rnd(data.length-start+1);let sub=0n;for(const byte of data.subarray(start,start+length))sub=(sub*base+BigInt(byte+1))%prime;
   assert.equal(ropeSubstrHash(root,start,length,hash),sub);
   if(step%40===0)snapshots.push([root,new Uint8Array(data)]);
  }
  for(const [old,bytes] of snapshots)assert.deepEqual(ropeToBytes(old),bytes);
 }
});
