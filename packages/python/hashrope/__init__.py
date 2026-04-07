"""
hashrope -- Hash rope data structure with algebraic content addressing.

A balanced binary tree (BB[2/7]) augmented with polynomial hash metadata
at every node. Supports O(log w) concat, split, and substring hashing,
plus O(log q) repetition encoding via RepeatNode and geometric accumulator.

All arithmetic is performed modulo a Mersenne prime using bit-shift
reduction (no expensive division).
"""

__version__ = "0.2.0"

from hashrope.polynomial_hash import (
    PolynomialHash,
    phi,
    mersenne_mod,
    mersenne_mul,
    MERSENNE_61,
    MERSENNE_127,
)
from hashrope.rope import (
    Leaf,
    Internal,
    RepeatNode,
    Node,
    rope_concat,
    rope_split,
    rope_repeat,
    rope_substr_hash,
    rope_len,
    rope_hash,
    rope_height,
    rope_from_bytes,
    rope_to_bytes,
    validate_rope,
)
from hashrope.sliding import SlidingWindow

__all__ = [
    # Polynomial hash
    "PolynomialHash",
    "phi",
    "mersenne_mod",
    "mersenne_mul",
    "MERSENNE_61",
    "MERSENNE_127",
    # Rope nodes
    "Leaf",
    "Internal",
    "RepeatNode",
    "Node",
    # Rope operations
    "rope_concat",
    "rope_split",
    "rope_repeat",
    "rope_substr_hash",
    "rope_len",
    "rope_hash",
    "rope_height",
    "rope_from_bytes",
    "rope_to_bytes",
    "validate_rope",
    # Sliding window
    "SlidingWindow",
]
