"""Burrows-Wheeler Transform and exact matching utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

from .validation import normalize_dna

SENTINEL = "$"


@dataclass(frozen=True)
class BWTIndex:
    """FM-index style structure for exact pattern matching."""

    text: str
    bwt: str
    suffix_array: Tuple[int, ...]
    first_occurrence: Dict[str, int]
    checkpoints: Dict[str, Tuple[int, ...]]

    @classmethod
    def build(cls, sequence: str) -> "BWTIndex":
        """Build an exact-match index for a DNA sequence."""
        normalized = normalize_dna(sequence)
        text = _with_unique_sentinel(normalized)
        bwt, suffix_array = bwt_transform(text, already_terminated=True)
        alphabet = sorted(set(bwt))

        seen = set()
        first_occurrence = {}
        for index, char in enumerate(sorted(bwt)):
            if char not in seen:
                first_occurrence[char] = index
                seen.add(char)

        counts = {char: [0] for char in alphabet}
        running = {char: 0 for char in alphabet}
        for char in bwt:
            running[char] += 1
            for symbol in alphabet:
                counts[symbol].append(running[symbol])

        checkpoints = {char: tuple(values) for char, values in counts.items()}
        return cls(text, bwt, tuple(suffix_array), first_occurrence, checkpoints)

    def exact_match(self, pattern: str) -> List[int]:
        """Return sorted zero-based start positions for an exact DNA pattern."""
        normalized = normalize_dna(pattern, allow_empty=False)
        top = 0
        bottom = len(self.bwt)

        for symbol in reversed(normalized):
            if symbol not in self.first_occurrence:
                return []
            counts = self.checkpoints[symbol]
            top = self.first_occurrence[symbol] + counts[top]
            bottom = self.first_occurrence[symbol] + counts[bottom]
            if top >= bottom:
                return []

        sentinel_position = len(self.text) - 1
        positions = [pos for pos in self.suffix_array[top:bottom] if pos != sentinel_position]
        return sorted(positions)


def bwt_transform(sequence: str, *, already_terminated: bool = False) -> Tuple[str, List[int]]:
    """Return the BWT string and suffix array for a sequence.

    By default a unique ``$`` sentinel is appended. Pass ``already_terminated=True``
    only when the caller has already appended a single trailing sentinel.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    text = sequence if already_terminated else _with_unique_sentinel(normalize_dna(sequence))
    if already_terminated:
        _validate_terminated_text(text)

    suffix_array = sorted(range(len(text)), key=lambda index: text[index:])
    transformed = "".join(text[index - 1] for index in suffix_array)
    return transformed, suffix_array


def inverse_bwt(transformed: str) -> str:
    """Invert a BWT string and return the original sequence without the sentinel."""
    if not isinstance(transformed, str):
        raise TypeError("transformed must be a string")
    if transformed.count(SENTINEL) != 1:
        raise ValueError("BWT input must contain exactly one '$' sentinel")

    table = [""] * len(transformed)
    for _ in transformed:
        table = sorted(char + row for char, row in zip(transformed, table))

    for row in table:
        if row.endswith(SENTINEL):
            return row[:-1]
    raise ValueError("could not invert malformed BWT input")


def find_exact_matches(sequence: str, patterns: Iterable[str]) -> Dict[str, List[int]]:
    """Find exact zero-based pattern positions in a DNA sequence."""
    index = BWTIndex.build(sequence)
    return {pattern: index.exact_match(pattern) for pattern in patterns}


def _with_unique_sentinel(sequence: str) -> str:
    if SENTINEL in sequence:
        raise ValueError("sequence must not contain '$'; it is reserved as the BWT sentinel")
    return f"{sequence}{SENTINEL}"


def _validate_terminated_text(text: str) -> None:
    if text.count(SENTINEL) != 1 or not text.endswith(SENTINEL):
        raise ValueError("terminated BWT text must end with exactly one '$' sentinel")
    normalize_dna(text[:-1])
