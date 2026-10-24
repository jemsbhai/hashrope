import test from 'node:test';
import assert from 'node:assert/strict';
import { PolynomialHash, MERSENNE_61, MERSENNE_127 } from '../dist/esm/polynomial.js';
import { Leaf, Internal, RepeatNode, ropeHash, ropeLen, ropeToBytes, validateRope, structuralEquals } from '../dist/esm/nodes.js';
import { ropeConcat, ropeSplit, ropeRepeat, ropeSubstrHash, ropeInsert, ropeDelete, contentEquals } from '../dist/esm/rope.js';
const bytes = s => new TextEncoder().encode(s);
const text = n => new TextDecoder().decode(ropeToBytes(n));
const h = new PolynomialHash();
const leaf = s => new Leaf(bytes(s), h);

test('concat identities, persistence, metadata and byte split endpoints', () => {
 const a=leaf('hello '), b=leaf('world'), n=ropeConcat(a,b,h);
 assert.equal(ropeConcat(null,a,h),a); assert.equal(ropeConcat(a,null,h),a);
 assert.equal(text(n),'hello world'); assert.equal(ropeHash(n),h.hash(bytes('hello world')));
 assert.equal(n.left,a); assert.equal(n.right,b);
 assert.deepEqual(ropeSplit(n,6n,h),[a,b]);
 assert.deepEqual(ropeSplit(n,-100n,h),[null,n]); assert.deepEqual(ropeSplit(n,100n,h),[n,null]);
 const [l,r]=ropeSplit(n,3n,h); assert.equal(text(l),'hel'); assert.equal(text(r),'lo world');
 assert.equal(text(n),'hello world'); validateRope(l); validateRope(r);
});

test('balanced append/prepend and unequal byte-sized leaves',()=>{
 let n=null; let expected='';
 for(let i=0;i<350;i++) { const s=String.fromCharCode(33+i%90).repeat(i%7+1); n=i%2?ropeConcat(n,leaf(s),h):ropeConcat(leaf(s),n,h); expected=i%2?expected+s:s+expected; validateRope(n); assert.equal(ropeHash(n),h.hash(bytes(expected))); }
 assert.equal(text(n),expected);
 const large=leaf('x'.repeat(10000)); const joined=ropeConcat(large,leaf('z'),h); assert.equal(joined.weight,2n); validateRope(joined);
});

test('joins and splits repeated children without expanding logical length',()=>{
 const p=new Internal(leaf('a'),leaf('bc'),h), huge=new RepeatNode(p,1n<<100n,h);
 for(const q of [0n,1n,2n,3n,100n,(1n<<99n)*3n+1n,huge.len-1n,huge.len]) {
  const [l,r]=ropeSplit(huge,q,h); assert.equal(ropeLen(l),q); assert.equal(ropeLen(r),huge.len-q);
  assert.equal(h.hashConcat(ropeHash(l),ropeLen(r),ropeHash(r)),huge.hashVal); validateRope(l);validateRope(r);
 }
 for(const count of [2n,3n,7n,31n,100n]) {
  const repeated=new RepeatNode(p,count,h);
  for(const n of [ropeConcat(leaf('z'),repeated,h),ropeConcat(repeated,leaf('z'),h)]) validateRope(n);
  for(let i=0n;i<=repeated.len;i++) { const [l,r]=ropeSplit(repeated,i,h); assert.equal(text(l)+text(r),'abc'.repeat(Number(count))); validateRope(l);validateRope(r); }
 }
});
