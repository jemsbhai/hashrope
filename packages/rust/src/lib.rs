//! # hashrope
//!
//! A BB[2/7] weight-balanced binary tree augmented with polynomial hash
//! metadata at every node.
//!
//! Zero dependencies. `no_std + alloc` compatible.
//! All arithmetic uses Mersenne prime modular reduction.

#![no_std]
extern crate alloc;

pub mod polynomial_hash;
pub mod rope;
pub mod sliding;

pub use polynomial_hash::{mersenne_mod, mersenne_mul, phi, PolynomialHash, MERSENNE_61};
pub use rope::{
    rope_concat, rope_from_bytes, rope_hash, rope_len, rope_repeat, rope_split,
    rope_substr_hash, rope_to_bytes, validate_rope, Node,
};
pub use sliding::SlidingWindow;
