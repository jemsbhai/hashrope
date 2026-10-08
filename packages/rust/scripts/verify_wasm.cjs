const fs = require('node:fs');
const assert = require('node:assert/strict');

const path = process.argv[2];
assert.ok(path, 'Usage: node scripts/verify_wasm.cjs target/numeric_probe.wasm');
const { exports: wasm } = new WebAssembly.Instance(
  new WebAssembly.Module(fs.readFileSync(path)), {}
);
// Values from Python pow(131, n, 2**61-1), not a candidate-to-candidate comparison.
assert.equal(wasm.power(4294967296n), 98942608713749192n);
assert.equal(wasm.power(4294967551n), 1154367712108772148n);
assert.equal(wasm.power(18446744073709551615n), 523719968513806577n);
assert.equal(wasm.small_modulus(), 1n);
assert.equal(wasm.max_rope_hash(), 1849259105152891911n);
assert.throws(() => wasm.unaddressable_bytes(), WebAssembly.RuntimeError);
assert.throws(() => wasm.overflowing_repeat(), WebAssembly.RuntimeError);
console.log('PASS: wasm32 full exponents, small modulus, maximum compressed rope, checked capacity');
