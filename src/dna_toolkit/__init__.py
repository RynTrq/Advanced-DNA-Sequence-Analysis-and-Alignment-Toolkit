"""Production-ready DNA sequence analysis primitives."""

from .alignment import (
    AlignmentResult,
    edit_distance_with_gap_penalties,
    hybrid_alignment,
    smith_waterman,
    smith_waterman_with_compression,
)
from .bwt import BWTIndex, bwt_transform, find_exact_matches, inverse_bwt
from .validation import DNA_ALPHABET, normalize_dna, validate_dna

__all__ = [
    "AlignmentResult",
    "BWTIndex",
    "DNA_ALPHABET",
    "bwt_transform",
    "edit_distance_with_gap_penalties",
    "find_exact_matches",
    "hybrid_alignment",
    "inverse_bwt",
    "normalize_dna",
    "smith_waterman",
    "smith_waterman_with_compression",
    "validate_dna",
]

