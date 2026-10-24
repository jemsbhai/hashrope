"""
Hash rope data structure for compressed-domain hashing.

Implements Part III of the framework:
- Definition 6: Leaf, Internal, RepeatNode
- Definition 7: Invariants I1-I9 (BB[2/7] weight balance)
- Lemma 4: Hash correctness (structural induction)
- Lemma 5: Height bound O(log w)
- Lemma 6: Rotations preserve hash invariant
- Lemma 7: RepeatNode splittability
- Lemma 8: Join complexity O(k · |h_L - h_R|)
- Lemma 8': RepeatNode bisection for rebalancing
- Theorem 6: Join (full specification)
- Theorem 7: Split in O(k · log w)
- Theorem 8: Repeat via RepeatNode in O(k · log q), O(1) space
- Theorem 9: Allocation-free SubstrHash in O(k · log w)
- Theorem 10: Concat = Join

All nodes are immutable (Invariant I9). Operations return new trees
with structural sharing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import gcd
from typing import Optional

from hashrope.polynomial_hash import (
    PolynomialHash,
    HashProfile,
    phi,
    mersenne_mod,
    mersenne_mul,
    MERSENNE_61,
)

# ---------------------------------------------------------------------------
# Balance parameter (Definition 7, Invariant I8)
# ---------------------------------------------------------------------------

# BB[α] with α = 2/7 (Nievergelt-Reingold condition).
# A node is balanced if each child's weight is at least α * parent weight.
ALPHA_NUM = 2
ALPHA_DEN = 7

# Type alias
Node = Optional["Leaf | Internal | RepeatNode"]


def _require_profile(node: Node, h: PolynomialHash) -> None:
    """Check provenance in O(1); constructors preserve it throughout the tree."""
    if node is not None and node.profile is not h.profile and node.profile != h.profile:
        raise ValueError("Incompatible hash profile: rebuild the rope with the requested hasher")


# ---------------------------------------------------------------------------
# Definition 6: Node types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Leaf:
    """
    Leaf node storing a byte sequence (Definition 6a).

    Invariant I1: hash_val = H(data).
    Weight is always 1.
    """

    data: bytes
    hash_val: int
    len: int
    weight: int
    profile: HashProfile = field(compare=False, repr=False)

    def __init__(self, data: bytes, h: PolynomialHash) -> None:
        if len(data) == 0:
            raise ValueError("Leaf cannot be empty")
        object.__setattr__(self, "data", bytes(data))
        object.__setattr__(self, "len", len(data))
        object.__setattr__(self, "hash_val", h.hash(data))
        object.__setattr__(self, "weight", 1)
        object.__setattr__(self, "profile", h.profile)


@dataclass(frozen=True, slots=True)
class Internal:
    """
    Internal node joining two subtrees (Definition 6b).

    Invariant I2: len = left.len + right.len
    Invariant I3: hash_val = left.hash_val · x^(right.len) + right.hash_val
    Invariant I4: weight = left.weight + right.weight
    """

    left: Leaf | Internal | RepeatNode
    right: Leaf | Internal | RepeatNode
    hash_val: int
    len: int
    weight: int
    profile: HashProfile = field(compare=False, repr=False)

    def __init__(
        self,
        left: Leaf | Internal | RepeatNode,
        right: Leaf | Internal | RepeatNode,
        h: PolynomialHash,
    ) -> None:
        _require_profile(left, h)
        _require_profile(right, h)
        object.__setattr__(self, "left", left)
        object.__setattr__(self, "right", right)
        object.__setattr__(self, "len", left.len + right.len)
        object.__setattr__(self, "weight", left.weight + right.weight)
        # I3: hash = left.hash · x^(right.len) + right.hash
        hv = h.hash_concat(left.hash_val, right.len, right.hash_val)
        object.__setattr__(self, "hash_val", hv)
        object.__setattr__(self, "profile", h.profile)


@dataclass(frozen=True, slots=True)
class RepeatNode:
    """
    Repeat node representing child^reps (Definition 6c).

    Invariant I5: len = child.len * reps
    Invariant I6: hash_val = child.hash_val · Φ(reps, x^(child.len))
    Invariant I7: weight = child.weight * reps
    """

    child: Leaf | Internal | RepeatNode
    reps: int
    hash_val: int
    len: int
    weight: int
    profile: HashProfile = field(compare=False, repr=False)

    def __init__(
        self,
        child: Leaf | Internal | RepeatNode,
        reps: int,
        h: PolynomialHash,
    ) -> None:
        if reps < 2:
            raise ValueError(f"RepeatNode reps must be >= 2, got {reps}")
        _require_profile(child, h)
        object.__setattr__(self, "child", child)
        object.__setattr__(self, "reps", reps)
        object.__setattr__(self, "len", child.len * reps)
        object.__setattr__(self, "weight", child.weight * reps)
        # I6: hash = child.hash · Φ(reps, x^d)
        x_d = h.power(child.len)
        phi_val = phi(reps, x_d, h.prime)
        hv = mersenne_mul(child.hash_val, phi_val, h.prime)
        object.__setattr__(self, "hash_val", hv)
        object.__setattr__(self, "profile", h.profile)


# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------


def rope_len(node: Node) -> int:
    """Total byte length of the string represented by node."""
    return 0 if node is None else node.len


def rope_hash(node: Node) -> int:
    """Fingerprint of the represented string; equality of hashes is not proof of equality."""
    return 0 if node is None else node.hash_val


def rope_equal(left: Node, right: Node) -> bool:
    """Exact byte equality, independent of tree shape and hash profile.

    Shared immutable ranges and proven-equal periodic spans are skipped.
    Unshared content may require a linear scan. No expanded rope is allocated.
    """
    if left is right:
        return True
    if rope_len(left) != rope_len(right):
        return False
    # Also independent of custom hash implementations and historically accepted
    # non-Mersenne moduli: exact comparison does not consult hash metadata.
    return rope_lcp(left, right) == rope_len(left)


def _normalize_range(
    node: Leaf | Internal | RepeatNode, start: int, length: int,
) -> tuple[Leaf | Internal | RepeatNode, int]:
    """Descend through a containing child, keeping a virtual byte range."""
    while True:
        if isinstance(node, Internal):
            boundary = node.left.len
            if start >= boundary:
                start -= boundary
                node = node.right
            elif start + length <= boundary:
                node = node.left
            else:
                return node, start
        elif isinstance(node, RepeatNode):
            start %= node.child.len
            if start + length <= node.child.len:
                node = node.child
            else:
                return node, start
        else:
            return node, start


def _leaf_lcp(left: bytes, a: int, right: bytes, b: int, length: int) -> int:
    """Compare in native byte chunks, with bounded temporary slice storage."""
    if left is right and a == b:
        return length
    matched = 0
    while matched < length:
        size = min(65536, length - matched)
        l_chunk = left[a + matched:a + matched + size]
        r_chunk = right[b + matched:b + matched + size]
        if l_chunk != r_chunk:
            # Keep early mismatches cheap. For a late mismatch, native byte
            # comparisons avoid a Python loop over the entire matching prefix.
            if l_chunk[0] != r_chunk[0]:
                return matched
            # Small offsets should not trigger probes over large slices.
            # The bounded scalar loop is used only when this prefix differs.
            if l_chunk[:64] != r_chunk[:64]:
                for i in range(1, min(64, size)):
                    if l_chunk[i] != r_chunk[i]:
                        return matched + i
            lo, hi = 64, size
            # The first mismatch is in [lo, hi). Equal left halves are skipped
            # exactly; unequal ones contain the first mismatch. The intervals
            # halve, so all probes together scan/copy O(size) bytes, not
            # O(size * log(size)). Temporary slices remain bounded by 64 KiB.
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if l_chunk[lo:mid] == r_chunk[lo:mid]:
                    lo = mid
                else:
                    hi = mid
            return matched + lo
        matched += size
    return length


def rope_lcp(left: Node, right: Node) -> int:
    """Exact longest common prefix length in bytes, across any hash profiles.

    Walks virtual ranges without splitting or flattening the ropes. Identical
    immutable ranges are skipped. Two periodic spans need comparison of only
    p + q - gcd(p, q) bytes to prove equality throughout (Fine-Wilf), where p
    and q are their periods. No fingerprints are used as equality evidence.

    Independent equal inputs may require reading the entire common prefix,
    plus traversal, cache and integer-arithmetic costs. There is no general
    logarithmic bound for independently built inputs.
    Auxiliary storage is the traversal stack, a per-call set of verified
    non-leaf ranges (for shared DAGs), and bounded leaf buffers.
    """
    length = min(rope_len(left), rope_len(right))
    if left is right or length == 0:
        return length
    return _range_lcp(left, 0, right, 0, length)


def _check_range(node: Node, start: int, length: int) -> None:
    if start < 0 or length < 0 or start > rope_len(node) or length > rope_len(node) - start:
        raise ValueError("Comparison range is outside the rope")


def rope_range_equal(
    left: Node, left_start: int, right: Node, right_start: int, length: int,
) -> bool:
    """Exactly compare two byte ranges without splitting or rehashing.

    Hash profiles may differ. Negative or out-of-bounds ranges raise ValueError,
    including empty ranges whose start exceeds the rope's length.
    """
    return rope_range_lcp(left, left_start, right, right_start, length) == length


def rope_range_lcp(
    left: Node, left_start: int, right: Node, right_start: int, length: int,
) -> int:
    """Exact LCP of two length-byte ranges, relative to their starts.

    Bounds are checked even for empty or identical ranges. The comparison has
    the same sharing, repetition, and worst-case costs as rope_lcp.
    """
    _check_range(left, left_start, length)
    _check_range(right, right_start, length)
    if length == 0 or (left is right and left_start == right_start):
        return length
    return _range_lcp(left, left_start, right, right_start, length)


def _range_lcp(
    left: Node, left_start: int, right: Node, right_start: int, length: int,
) -> int:

    # Integer entries skip a periodic suffix only AFTER its witness compares
    # equal. Range entries are processed left-to-right, so matched is the
    # absolute offset of the first mismatch if any comparison fails.
    pending = [(left, left_start, right, right_start, length)]
    verified = set()
    matched = 0
    while pending:
        task = pending.pop()
        if isinstance(task, int):
            matched += task
            continue
        if len(task) == 1:
            verified.add(task[0])
            continue
        a_node, a_start, b_node, b_start, size = task
        a_node, a_start = _normalize_range(a_node, a_start, size)
        b_node, b_start = _normalize_range(b_node, b_start, size)
        if a_node is b_node and a_start == b_start:
            matched += size
            continue
        if isinstance(a_node, Leaf) and isinstance(b_node, Leaf):
            common = _leaf_lcp(a_node.data, a_start, b_node.data, b_start, size)
            matched += common
            if common != size:
                return matched
            continue

        key = (id(a_node), a_start, id(b_node), b_start, size)
        if key in verified:
            matched += size
            continue
        # Commit equality only after every child task has succeeded. This
        # keeps independently constructed self-concatenation DAGs compressed.
        pending.append((key,))

        if isinstance(a_node, RepeatNode) and isinstance(b_node, RepeatNode):
            p, q = a_node.child.len, b_node.child.len
            witness = p + q - gcd(p, q)
            if witness < size:
                pending.append(size - witness)
                pending.append((a_node, a_start, b_node, b_start, witness))
                continue

        # Prefer internal boundaries: cutting a Repeat at every copy before
        # exposing the other tree's large shared/repeated child would expand it.
        cut = size
        if isinstance(a_node, Internal):
            cut = min(cut, a_node.left.len - a_start)
        if isinstance(b_node, Internal):
            cut = min(cut, b_node.left.len - b_start)
        if cut == size:
            if isinstance(a_node, RepeatNode):
                cut = min(cut, a_node.child.len - a_start)
            if isinstance(b_node, RepeatNode):
                cut = min(cut, b_node.child.len - b_start)
        pending.append((a_node, a_start + cut, b_node, b_start + cut, size - cut))
        pending.append((a_node, a_start, b_node, b_start, cut))
    return matched


def rope_height(node: Node) -> int:
    """Height of the rope tree. Leaves have height 0, None has height 0."""
    if node is None:
        return 0
    if isinstance(node, Leaf):
        return 0
    if isinstance(node, RepeatNode):
        return 1 + rope_height(node.child)
    # Internal
    return 1 + max(rope_height(node.left), rope_height(node.right))


def _weight(node: Node) -> int:
    return 0 if node is None else node.weight


def _height(node: Node) -> int:
    """Approximate height for join decisions."""
    if node is None:
        return 0
    if isinstance(node, Leaf):
        return 1
    if isinstance(node, RepeatNode):
        return _height(node.child) + 1
    # Internal
    return 1 + max(_height(node.left), _height(node.right))


# ---------------------------------------------------------------------------
# Balance helpers (Invariant I8)
# ---------------------------------------------------------------------------


def _is_balanced(left: Node, right: Node) -> bool:
    """Check BB[2/7] condition for a potential Internal(left, right)."""
    total = _weight(left) + _weight(right)
    if total <= 2:
        return True
    wl = _weight(left)
    # α ≤ wl/total ≤ 1-α  ⟺  α·total ≤ wl ≤ (1-α)·total
    # With α = 2/7: 2·total ≤ 7·wl and 7·wl ≤ 5·total
    return ALPHA_NUM * total <= ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total


def _decompose(
    node: Leaf | Internal | RepeatNode, h: PolynomialHash
) -> tuple[Leaf | Internal | RepeatNode, Leaf | Internal | RepeatNode]:
    """
    Decompose a node into two children for balancing (Lemma 8').

    Internal → its two children.
    RepeatNode → split by reps halving (not byte position).
    Leaf → should never be called on a Leaf.
    """
    if isinstance(node, Internal):
        return node.left, node.right
    if isinstance(node, RepeatNode):
        half = node.reps // 2
        left = _make_repeat(node.child, half, h)
        right = _make_repeat(node.child, node.reps - half, h)
        return left, right
    raise TypeError(f"Cannot decompose Leaf in balance")


def _is_balanced_wt(wl: int, wr: int) -> bool:
    """Check BB[2/7] condition using weights only."""
    total = wl + wr
    if total <= 2:
        return True
    return ALPHA_NUM * total <= ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total


def _balance(
    left: Leaf | Internal | RepeatNode,
    right: Leaf | Internal | RepeatNode,
    h: PolynomialHash,
) -> Leaf | Internal | RepeatNode:
    """
    Adams-style balance with checked rotations (replaces _rebalance).

    Handles RepeatNode via _decompose. Single rotation is used only when
    BOTH the inner node AND the outer pairing are balanced; otherwise
    double rotation decomposes further (bounded by subtree height).

    Termination: M = max(weight(left), weight(right)) strictly decreases
    at every recursive call. Proven analytically (worst-case contraction
    ratio 223/245 ≈ 0.918) and verified computationally over 9.2M
    combinations.
    """
    wl = _weight(left)
    wr = _weight(right)
    total = wl + wr

    if total <= 2:
        return Internal(left, right, h)

    if _is_balanced_wt(wl, wr):
        return Internal(left, right, h)

    if ALPHA_NUM * total > ALPHA_DEN * wl:
        # Left too light (right too heavy) → rotate left
        rl, rr = _decompose(right, h)
        wrl = _weight(rl)
        wrr = _weight(rr)
        # Single rotation: only if BOTH inner and outer are balanced
        if _is_balanced_wt(wl, wrl) and _is_balanced_wt(wl + wrl, wrr):
            new_left = Internal(left, rl, h)
            return Internal(new_left, rr, h)
        else:
            # Double rotation: decompose rl, balance all parts
            rll, rlr = _decompose(rl, h)
            new_left = _balance(left, rll, h)
            new_right = _balance(rlr, rr, h)
            return _balance(new_left, new_right, h)
    else:
        # Left too heavy (right too light) → rotate right
        ll, lr = _decompose(left, h)
        wll = _weight(ll)
        wlr = _weight(lr)
        # Single rotation: only if BOTH inner and outer are balanced
        if _is_balanced_wt(wlr, wr) and _is_balanced_wt(wll, wlr + wr):
            new_right = Internal(lr, right, h)
            return Internal(ll, new_right, h)
        else:
            # Double rotation: decompose lr, balance all parts
            lrl, lrr = _decompose(lr, h)
            new_left = _balance(ll, lrl, h)
            new_right = _balance(lrr, right, h)
            return _balance(new_left, new_right, h)


# ---------------------------------------------------------------------------
# Theorem 6 / Lemma 8: Join (Concat)
# ---------------------------------------------------------------------------


def rope_concat(
    left: Node, right: Node, h: PolynomialHash
) -> Node:
    """
    Join two ropes (Theorem 6, Theorem 10).

    Returns a valid hash rope representing left's string ‖ right's string.
    Time: O(k · |h_L - h_R|).
    """
    _require_profile(left, h)
    _require_profile(right, h)
    if left is None:
        return right
    if right is None:
        return left
    return _join(left, right, h)


def _join(
    left: Leaf | Internal | RepeatNode,
    right: Leaf | Internal | RepeatNode,
    h: PolynomialHash,
) -> Leaf | Internal | RepeatNode:
    """Core join with rebalancing."""
    if _is_balanced(left, right):
        return Internal(left, right, h)

    wl, wr = _weight(left), _weight(right)

    if wl > wr:
        # Left is heavier — descend right spine of left
        if isinstance(left, Internal):
            new_right = _join(left.right, right, h)
            return _balance(left.left, new_right, h)
        elif isinstance(left, RepeatNode):
            # Split by reps count (not byte position) for clean O(log q) recursion
            half = left.reps // 2
            ll = _make_repeat(left.child, half, h)
            lr = _make_repeat(left.child, left.reps - half, h)
            new_right = _join(lr, right, h)
            return _join(ll, new_right, h)
        else:
            # Leaf — just make Internal
            return Internal(left, right, h)
    else:
        # Right is heavier — descend left spine of right
        if isinstance(right, Internal):
            new_left = _join(left, right.left, h)
            return _balance(new_left, right.right, h)
        elif isinstance(right, RepeatNode):
            # Split by reps count (not byte position) for clean O(log q) recursion
            half = right.reps // 2
            rl = _make_repeat(right.child, half, h)
            rr = _make_repeat(right.child, right.reps - half, h)
            new_left = _join(left, rl, h)
            return _join(new_left, rr, h)
        else:
            return Internal(left, right, h)


# ---------------------------------------------------------------------------
# Theorem 7: Split
# ---------------------------------------------------------------------------


def rope_split(
    node: Node, pos: int, h: PolynomialHash
) -> tuple[Node, Node]:
    """
    Split rope at byte position pos (Theorem 7).

    Returns (left, right) where left has pos bytes, right has the rest.
    Time: O(k · log w).

    pos is a byte offset. If the rope stores UTF-8-encoded text, a byte
    offset may fall inside a multibyte character's encoding: the split is
    still exact over bytes (concatenating the halves reproduces the
    original bytes, and all hashes are over bytes), but each half on its
    own may not be valid UTF-8. Character-boundary alignment is the
    caller's responsibility. See the byte-index contract in the README.
    """
    _require_profile(node, h)
    if node is None:
        return None, None
    if pos <= 0:
        return None, node
    if pos >= node.len:
        return node, None
    return _split(node, pos, h)


def _split(
    node: Leaf | Internal | RepeatNode,
    pos: int,
    h: PolynomialHash,
) -> tuple[Node, Node]:
    """Core split. 0 < pos < node.len guaranteed."""

    if isinstance(node, Leaf):
        left = Leaf(node.data[:pos], h)
        right = Leaf(node.data[pos:], h)
        return left, right

    if isinstance(node, Internal):
        ll = node.left.len
        if pos == ll:
            return node.left, node.right
        if pos < ll:
            l1, l2 = _split(node.left, pos, h)
            return l1, rope_concat(l2, node.right, h)
        else:
            r1, r2 = _split(node.right, pos - ll, h)
            return rope_concat(node.left, r1, h), r2

    if isinstance(node, RepeatNode):
        return _split_repeat(node, pos, h)

    raise TypeError(f"Unknown node type: {type(node)}")


def _split_repeat(
    node: RepeatNode, pos: int, h: PolynomialHash
) -> tuple[Node, Node]:
    """
    Split a RepeatNode at byte position pos (Lemma 7).

    Let d = child.len, m = pos // d, r = pos % d.
    """
    d = node.child.len
    m, r = divmod(pos, d)
    q = node.reps

    if r == 0:
        # Split falls on repetition boundary
        left = _make_repeat(node.child, m, h)
        right = _make_repeat(node.child, q - m, h)
        return left, right
    else:
        # Split falls within a repetition
        child_left, child_right = _split(node.child, r, h)

        # Left = Repeat(child, m) ‖ child_left
        left_rep = _make_repeat(node.child, m, h)
        left = rope_concat(left_rep, child_left, h)

        # Right = child_right ‖ Repeat(child, q - m - 1)
        right_rep = _make_repeat(node.child, q - m - 1, h)
        right = rope_concat(child_right, right_rep, h)

        return left, right


def _make_repeat(
    child: Leaf | Internal | RepeatNode, reps: int, h: PolynomialHash
) -> Node:
    """Create a Repeat with proper handling of reps 0 and 1."""
    if reps == 0:
        return None
    if reps == 1:
        return child
    return RepeatNode(child, reps, h)


# ---------------------------------------------------------------------------
# Theorem 8: Repeat
# ---------------------------------------------------------------------------


def rope_repeat(
    node: Leaf | Internal | RepeatNode,
    q: int,
    h: PolynomialHash,
) -> Node:
    """
    Create a RepeatNode representing node^q (Theorem 8).

    Time: O(k · log q) for the Φ computation. O(1) new nodes.
    """
    _require_profile(node, h)
    if q == 0:
        return None
    if q == 1:
        return node
    return RepeatNode(node, q, h)


# ---------------------------------------------------------------------------
# Theorem 9: SubstrHash (allocation-free)
# ---------------------------------------------------------------------------


def rope_substr_hash(
    node: Node,
    start: int,
    length: int,
    h: PolynomialHash,
) -> int:
    """
    Compute H(S[start..start+length-1]) without allocating nodes (Theorem 9).

    Time: O(k · log w). Space: O(log w) stack.

    start and length are byte offsets/counts (see the byte-index contract
    in the README).
    """
    _require_profile(node, h)
    if node is None or length == 0:
        return 0
    return _hash_range(node, start, length, h)


def _hash_range(
    node: Leaf | Internal | RepeatNode,
    start: int,
    length: int,
    h: PolynomialHash,
) -> int:
    """Recursive hash range query."""
    p = h.prime

    # Theorem 9 fast path: a fully covered node answers from stored
    # metadata (Invariants I1/I3/I6) in O(1). This realizes the
    # O(k · log w) bound — the canonical decomposition of any range
    # touches O(log w) fully covered nodes plus at most two partial
    # leaves at the ends.
    if start == 0 and length == node.len:
        return node.hash_val

    if isinstance(node, Leaf):
        return h.hash(node.data[start:start + length])

    if isinstance(node, Internal):
        ll = node.left.len
        if start + length <= ll:
            # Entirely in left
            return _hash_range(node.left, start, length, h)
        if start >= ll:
            # Entirely in right
            return _hash_range(node.right, start - ll, length, h)
        # Spanning
        l_len = ll - start
        r_len = length - l_len
        h_l = _hash_range(node.left, start, l_len, h)
        h_r = _hash_range(node.right, 0, r_len, h)
        return h.hash_concat(h_l, r_len, h_r)

    if isinstance(node, RepeatNode):
        d = node.child.len
        # Map into the repeating pattern
        first_copy = start // d
        start_in_copy = start % d

        end = start + length - 1
        last_copy = end // d
        end_in_copy = end % d

        if first_copy == last_copy:
            # Entirely within one copy
            return _hash_range(node.child, start_in_copy, length, h)

        # Three parts: tail of first copy, full copies, head of last copy
        # 1. Tail of first copy: child[start_in_copy .. d-1]
        tail_len = d - start_in_copy
        h_tail = _hash_range(node.child, start_in_copy, tail_len, h)

        # 2. Full copies in the middle
        full_copies = last_copy - first_copy - 1
        if full_copies > 0:
            x_d = h.power(d)
            phi_val = phi(full_copies, x_d, p)
            h_full = mersenne_mul(node.child.hash_val, phi_val, p)
        else:
            h_full = 0
        full_len = full_copies * d

        # 3. Head of last copy: child[0 .. end_in_copy]
        head_len = end_in_copy + 1
        h_head = _hash_range(node.child, 0, head_len, h)

        # Combine: tail ‖ full ‖ head
        # H(tail ‖ full ‖ head) = H(tail) · x^(full_len + head_len) + H(full) · x^(head_len) + H(head)
        result = mersenne_mod(
            h.power(full_len + head_len) * h_tail
            + h.power(head_len) * h_full
            + h_head,
            p,
        )
        return result

    raise TypeError(f"Unknown node type: {type(node)}")


# ---------------------------------------------------------------------------
# Reconstruction: rope → bytes
# ---------------------------------------------------------------------------


def rope_to_bytes(node: Node) -> bytes:
    """Reconstruct the byte string from a rope (for testing)."""
    if node is None:
        return b""
    parts: list[bytes] = []
    _collect_bytes(node, parts)
    return b"".join(parts)


def _collect_bytes(node: Leaf | Internal | RepeatNode, parts: list[bytes]) -> None:
    if isinstance(node, Leaf):
        parts.append(node.data)
    elif isinstance(node, Internal):
        _collect_bytes(node.left, parts)
        _collect_bytes(node.right, parts)
    elif isinstance(node, RepeatNode):
        for _ in range(node.reps):
            _collect_bytes(node.child, parts)


# ---------------------------------------------------------------------------
# Construction helper
# ---------------------------------------------------------------------------


def rope_from_bytes(data: bytes, h: PolynomialHash) -> Node:
    """Create a rope from a byte string. Returns None for empty.

    All rope operations index into these bytes by byte offset (see the
    byte-index contract in the README).
    """
    if len(data) == 0:
        return None
    return Leaf(data, h)


# ---------------------------------------------------------------------------
# Validation (for testing)
# ---------------------------------------------------------------------------


def validate_rope(node: Node, h: PolynomialHash | None = None) -> None:
    """
    Validate structure and profile consistency; does not recompute fingerprints.

    Raises AssertionError on structural violations and ValueError if the supplied
    hasher has an incompatible profile. This is not a content equality test.
    """
    if node is None:
        return
    if h is not None:
        _require_profile(node, h)

    if isinstance(node, Leaf):
        assert node.weight == 1, f"Leaf weight must be 1, got {node.weight}"
        assert node.len == len(node.data), "Leaf len mismatch"
        return

    if isinstance(node, Internal):
        assert node.profile == node.left.profile == node.right.profile, "Internal profile mismatch"
        # I2
        assert node.len == node.left.len + node.right.len, "Internal len mismatch"
        # I4
        assert node.weight == node.left.weight + node.right.weight, "Internal weight mismatch"
        # I8: balance
        total = node.weight
        wl = node.left.weight
        if total > 2:
            assert ALPHA_NUM * total <= ALPHA_DEN * wl, (
                f"Left child too light: {wl}/{total}"
            )
            assert ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total, (
                f"Left child too heavy: {wl}/{total}"
            )
        # Recurse
        validate_rope(node.left, h)
        validate_rope(node.right, h)
        return

    if isinstance(node, RepeatNode):
        assert node.profile == node.child.profile, "RepeatNode profile mismatch"
        # I5
        assert node.len == node.child.len * node.reps, "RepeatNode len mismatch"
        # I7
        assert node.weight == node.child.weight * node.reps, "RepeatNode weight mismatch"
        assert node.reps >= 2, f"RepeatNode reps must be >= 2, got {node.reps}"
        validate_rope(node.child, h)
        return

    raise TypeError(f"Unknown node type: {type(node)}")
