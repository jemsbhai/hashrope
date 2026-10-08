//! # hashrope
//!
//! A BB[2/7] weight-balanced binary tree augmented with polynomial hash
//! metadata at every node.
//!
//! Zero dependencies by default. `no_std + alloc` compatible, including
//! the optional `bigint` feature for Python-compatible arbitrary precision.
//! The default API uses bit-folding Mersenne modular reduction; the optional
//! arbitrary-precision API uses `BigUint` modular arithmetic.
//! Nodes are stored in a contiguous arena for cache-friendly access.

#![no_std]
extern crate alloc;

/// The 127-bit Mersenne prime. Use the optional `bigint` module for hashing
/// with this modulus; the existing `u64` API retains its 61-bit default.
pub const MERSENNE_127: u128 = (1u128 << 127) - 1;

pub mod polynomial_hash;
pub mod rope;
pub mod sliding;

#[cfg(feature = "bigint")]
pub mod bigint;

pub use polynomial_hash::{mersenne_mod, mersenne_mul, phi, PolynomialHash, MERSENNE_61};
pub use rope::{
    rope_concat, rope_from_bytes, rope_hash, rope_height, rope_len, rope_repeat, rope_split,
    rope_substr_hash, rope_to_bytes, validate_rope, Arena, Node, NodeId, NodeInner, LAZY_SENTINEL,
};
pub use sliding::SlidingWindow;
