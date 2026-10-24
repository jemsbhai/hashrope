import test from 'node:test';
import assert from 'node:assert/strict';
import { PolynomialHash, MERSENNE_61, MERSENNE_127 } from '../dist/esm/polynomial.js';
import { Leaf, Internal, RepeatNode, ropeHash, ropeLen, ropeToBytes, validateRope } from '../dist/esm/nodes.js';
import { ropeConcat, ropeSplit, ropeRepeat, ropeSubstrHash, ropeInsert, ropeDelete, contentEquals } from '../dist/esm/rope.js';
const b=s=>new TextEncoder().encode(s), h=new PolynomialHash(), leaf=s=>new Leaf(b(s),h);

test('repeat public zero/one/negative and huge exact integer contracts',()=>{
 const a=leaf('ab'); assert.equal(ropeRepeat(a,0n,h),null); assert.equal(ropeRepeat(a,1n,h),a);
 assert.equal(ropeRepeat(null,0n,h),null); assert.equal(ropeRepeat(null,1n,h),null);
 assert.throws(()=>ropeRepeat(null,2n,h),TypeError); assert.throws(()=>ropeRepeat(a,-1n,h),RangeError);
 assert.throws(()=>ropeRepeat(a,Number.MAX_SAFE_INTEGER+1,h),RangeError);
 const huge=ropeRepeat(a,1n<<100n,h); assert.equal(huge.len,2n<<100n); assert.equal(huge.child,a);
});

test('substring every valid range across leaf/internal/repeat and UTF8 byte boundaries',()=>{
 const variants=[leaf('abcabcabc'),new RepeatNode(leaf('abc'),3n,h),new Internal(leaf('abc'),new RepeatNode(leaf('abc'),2n,h),h),leaf('aé🙂z')];
 for(const n of variants) { const data=ropeToBytes(n); for(let s=0;s<=data.length;s++) for(let len=0;len<=data.length-s;len++) assert.equal(ropeSubstrHash(n,s,len,h),h.hash(data.subarray(s,s+len))); }
 assert.equal(ropeSubstrHash(null,123n,5n,h),0n);
 const huge=new RepeatNode(leaf('abc'),1n<<100n,h);
 assert.equal(ropeSubstrHash(huge,(1n<<99n)*3n+1n,8n,h),h.hash(b('bcabcabc')));
 assert.throws(()=>ropeSubstrHash(huge,Number.MAX_SAFE_INTEGER+1,2,h),RangeError);
});

test('substring preserves shape-dependent out-of-domain release behavior',()=>{
 const a=leaf('ababab'),r=new RepeatNode(leaf('ab'),3n,h);
 assert.equal(ropeSubstrHash(a,6n,1n,h),0n); assert.equal(ropeSubstrHash(r,6n,1n,h),98n);
 assert.equal(ropeSubstrHash(a,-2n,2n,h),0n); assert.equal(ropeSubstrHash(r,-2n,2n,h),12937n);
 assert.equal(ropeSubstrHash(a,0n,-1n,h),h.hash(b('ababa')));
});

test('full coverage uses stored metadata and partial traversal does not rebuild nodes',()=>{
 class Counting extends PolynomialHash { hashes=0; concats=0; hash(data){this.hashes++;return super.hash(data)} hashConcat(...args){this.concats++;return super.hashConcat(...args)} }
 const c=new Counting(),a=new Leaf(b('abcdef'),c),n=new RepeatNode(a,1n<<100n,c);
 c.hashes=0;c.concats=0; assert.equal(ropeSubstrHash(n,0n,n.len,c),n.hashVal);assert.equal(c.hashes,0);assert.equal(c.concats,0);
 const old=n.child; assert.equal(ropeSubstrHash(n,2n,8n,c),c.hash(b('cdefabcd')));assert.equal(n.child,old);
});
