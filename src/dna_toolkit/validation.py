"""Input normalization and validation helpers."""

from __future__ import annotations

DNA_ALPHABET = frozenset("ACGTN")


def normalize_dna(sequence: str, *, allow_empty: bool = True) -> str:
    """Return an uppercase DNA sequence after validating its alphabet."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    normalized = "".join(sequence.split()).upper()
    validate_dna(normalized, allow_empty=allow_empty)
    return normalized


def validate_dna(sequence: str, *, allow_empty: bool = True) -> None:
    """Validate that a sequence contains only supported DNA symbols."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    if not allow_empty and not sequence:
        raise ValueError("sequence must not be empty")

    invalid = sorted(set(sequence.upper()) - DNA_ALPHABET)
    if invalid:
        invalid_symbols = ", ".join(repr(symbol) for symbol in invalid)
        raise ValueError(f"sequence contains unsupported DNA symbol(s): {invalid_symbols}")

