//! Runtime regression probe. Compile as a wasm32 cdylib; see VALIDATION.md.
#![allow(dead_code)]
extern crate alloc;
#[path = "../src/polynomial_hash.rs"]
mod polynomial_hash;
#[path = "../src/rope.rs"]
mod rope;

#[no_mangle]
pub extern "C" fn power(n: u64) -> u64 {
    polynomial_hash::PolynomialHash::default_hash().power(n)
}

#[no_mangle]
pub extern "C" fn small_modulus() -> u64 {
    polynomial_hash::mersenne_mod(256, 3)
}

#[no_mangle]
pub extern "C" fn max_rope_hash() -> u64 {
    let mut arena = rope::Arena::new();
    let leaf = arena.from_bytes(b"a");
    let left = arena.repeat(leaf, u64::MAX / 2);
    let right = arena.repeat(leaf, u64::MAX / 2 + 1);
    let node = arena.concat(left, right);
    arena.validate(node);
    assert_eq!(arena.len(node), u64::MAX);
    arena.substr_hash(node, 0, u64::MAX)
}

#[no_mangle]
pub extern "C" fn unaddressable_bytes() {
    let mut arena = rope::Arena::new();
    let leaf = arena.from_bytes(b"a");
    let node = arena.repeat(leaf, 1u64 << 32);
    arena.to_bytes(node);
}

#[no_mangle]
pub extern "C" fn overflowing_repeat() {
    let mut arena = rope::Arena::new();
    let leaf = arena.from_bytes(b"ab");
    arena.repeat(leaf, u64::MAX);
}
