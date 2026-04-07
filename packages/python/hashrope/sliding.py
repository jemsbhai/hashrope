"""
Sliding window hash rope for bounded-memory streaming hash computation.

Adapted from UHC's SlidingRopeState (Definition 10, Theorem 11).
Provides a generic sliding window over a hash rope without any
compression-format dependencies.

The window maintains a prefix hash and a bounded rope, evicting old
data as new data arrives. At any point, the full-sequence hash
H(T[0..pos-1]) can be computed in O(1) from the prefix hash and
the window rope hash.

Invariant I_slide (holds after every operation):
    (1) h_prefix = H(T[0 .. l_prefix - 1])
    (2) r_window represents T[l_prefix .. pos - 1]
    (3) rope_len(r_window) = pos - l_prefix <= W
    (4) H(T[0..pos-1]) = h_prefix * x^(pos - l_prefix) + rope_hash(r_window)
"""

from __future__ import annotations

from hashrope.polynomial_hash import PolynomialHash, MERSENNE_61
from hashrope.rope import (
    Leaf,
    Node,
    rope_concat,
    rope_split,
    rope_repeat,
    rope_len,
    rope_hash,
)


class SlidingWindow:
    """
    Bounded-memory sliding window over a hash rope.

    Maintains a window of at most W = d_max + m_max bytes as a rope,
    plus a prefix hash covering all evicted bytes. Supports two
    operations:

    - append_bytes(data): append literal bytes
    - append_copy(offset, length): append a copy from within the window

    At any point, current_hash() returns H(T[0..pos-1]) in O(1).

    Parameters
    ----------
    d_max : int
        Maximum copy-back distance. The window retains at least d_max
        bytes after eviction so that future copies can find their source.
    m_max : int
        Maximum copy length. Together with d_max, determines the
        maximum window extent W = d_max + m_max.
    prime : int
        Mersenne prime for hash arithmetic. Default: 2^61 - 1.
    base : int
        Hash base. Default: 131.
    """

    __slots__ = (
        "h", "d_max", "m_max", "W",
        "h_prefix", "l_prefix", "r_window", "pos",
    )

    def __init__(
        self,
        d_max: int = 32768,
        m_max: int = 258,
        prime: int = MERSENNE_61,
        base: int = 131,
    ) -> None:
        if d_max < 1:
            raise ValueError(f"d_max must be >= 1, got {d_max}")
        if m_max < 1:
            raise ValueError(f"m_max must be >= 1, got {m_max}")

        self.h = PolynomialHash(prime=prime, base=base)
        self.d_max = d_max
        self.m_max = m_max
        self.W = d_max + m_max

        # State (Definition 10)
        self.h_prefix: int = 0       # H(T[0..l_prefix-1])
        self.l_prefix: int = 0       # evicted byte count
        self.r_window: Node = None   # rope for T[l_prefix..pos-1]
        self.pos: int = 0            # total bytes ingested

    @property
    def window_len(self) -> int:
        """Current window size in bytes."""
        return rope_len(self.r_window)

    def current_hash(self) -> int:
        """
        Compute H(T[0..pos-1]) from sliding state via Theorem 1:
            H(T) = H_prefix * x^(window_len) + H(window)
        """
        wl = self.window_len
        wh = rope_hash(self.r_window) if self.r_window else 0
        return self.h.hash_concat(self.h_prefix, wl, wh)

    def final_hash(self) -> int:
        """
        FinalHash() -- Theorem 11.

        H(T[0..N-1]) = H_prefix * x^(R_window.len) + R_window.hash
        """
        return self.current_hash()

    def append_bytes(self, data: bytes) -> None:
        """
        Append literal bytes to the stream.

        Parameters
        ----------
        data : bytes
            One or more bytes to append.
        """
        if not data:
            return
        leaf = Leaf(data, self.h)
        self.r_window = rope_concat(self.r_window, leaf, self.h)
        self.pos += len(data)
        self._evict()

    def append_copy(self, offset: int, length: int) -> None:
        """
        Append a copy-from-window reference.

        Copies `length` bytes starting `offset` bytes back from the
        current end of the window. Handles both non-overlapping
        (offset >= length) and overlapping (offset < length) cases.

        Parameters
        ----------
        offset : int
            Distance back from the current end of the window.
            Must be <= pos and <= window_len.
        length : int
            Number of bytes to copy. For overlapping copies
            (offset < length), the pattern of length `offset`
            is repeated as needed.

        Raises
        ------
        ValueError
            If offset exceeds decoded length or window size.
        """
        win_len = self.window_len

        if offset > self.pos:
            raise ValueError(
                f"Invalid copy: offset {offset} exceeds "
                f"decoded length {self.pos}"
            )

        if offset > win_len:
            raise ValueError(
                f"Copy offset {offset} exceeds window size {win_len}. "
                f"Increase d_max (currently {self.d_max}) to at least {offset}."
            )

        start = win_len - offset

        if offset >= length:
            # Non-overlapping
            _, tmp = rope_split(self.r_window, start, self.h)
            source, _ = rope_split(tmp, length, self.h)
            self.r_window = rope_concat(self.r_window, source, self.h)
        else:
            # Overlapping: extract pattern of length `offset`, repeat
            _, tmp = rope_split(self.r_window, start, self.h)
            pattern, _ = rope_split(tmp, offset, self.h)

            q, r = divmod(length, offset)
            rep: Node = rope_repeat(pattern, q, self.h) if q >= 1 else None
            if r > 0:
                partial, _ = rope_split(pattern, r, self.h)
                rep = rope_concat(rep, partial, self.h)
            self.r_window = rope_concat(self.r_window, rep, self.h)

        self.pos += length
        self._evict()

    def _evict(self) -> None:
        """
        Evict excess bytes from the window.

        After eviction: rope_len(r_window) <= d_max, ensuring that
        future copies of distance <= d_max find their source.

        Correctness (Lemma 12): The prefix hash update follows from
        Theorem 1 (concatenation composability).
        """
        win_len = rope_len(self.r_window)
        if win_len <= self.W:
            return

        # Evict down to d_max bytes remaining
        excess = win_len - self.d_max
        r_old, r_keep = rope_split(self.r_window, excess, self.h)

        # Fold evicted bytes into h_prefix via Theorem 1
        old_len = rope_len(r_old)
        old_hash = rope_hash(r_old)
        self.h_prefix = self.h.hash_concat(self.h_prefix, old_len, old_hash)

        self.l_prefix += excess
        self.r_window = r_keep
