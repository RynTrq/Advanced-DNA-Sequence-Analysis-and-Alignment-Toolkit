"""Sequence alignment and edit-distance algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from math import inf
from typing import List, Optional, Tuple

from .bwt import bwt_transform, inverse_bwt
from .validation import normalize_dna

MATCH_STATE = 0
INSERT_STATE = 1
DELETE_STATE = 2


@dataclass(frozen=True)
class AlignmentResult:
    """Alignment result with score and local coordinates."""

    score: int
    aligned_seq1: str
    aligned_seq2: str
    start_seq1: int
    end_seq1: int
    start_seq2: int
    end_seq2: int


def edit_distance_with_gap_penalties(
    seq1: str,
    seq2: str,
    consecutive_gap_cost: int,
    isolated_gap_cost: int,
    *,
    mismatch_cost: int = 1,
) -> int:
    """Compute edit distance where opening and extending a gap have different costs."""
    left = normalize_dna(seq1)
    right = normalize_dna(seq2)
    _ensure_non_negative("consecutive_gap_cost", consecutive_gap_cost)
    _ensure_non_negative("isolated_gap_cost", isolated_gap_cost)
    _ensure_non_negative("mismatch_cost", mismatch_cost)

    n, m = len(left), len(right)
    match = [[inf] * (m + 1) for _ in range(n + 1)]
    insert = [[inf] * (m + 1) for _ in range(n + 1)]
    delete = [[inf] * (m + 1) for _ in range(n + 1)]
    match[0][0] = 0

    for j in range(1, m + 1):
        insert[0][j] = isolated_gap_cost if j == 1 else insert[0][j - 1] + consecutive_gap_cost
    for i in range(1, n + 1):
        delete[i][0] = isolated_gap_cost if i == 1 else delete[i - 1][0] + consecutive_gap_cost

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitution = 0 if left[i - 1] == right[j - 1] else mismatch_cost
            previous_best = min(match[i - 1][j - 1], insert[i - 1][j - 1], delete[i - 1][j - 1])
            match[i][j] = previous_best + substitution

            insert[i][j] = min(
                match[i][j - 1] + isolated_gap_cost,
                insert[i][j - 1] + consecutive_gap_cost,
                delete[i][j - 1] + isolated_gap_cost,
            )
            delete[i][j] = min(
                match[i - 1][j] + isolated_gap_cost,
                insert[i - 1][j] + isolated_gap_cost,
                delete[i - 1][j] + consecutive_gap_cost,
            )

    return int(min(match[n][m], insert[n][m], delete[n][m]))


def smith_waterman(
    seq1: str,
    seq2: str,
    *,
    match_score: int = 3,
    mismatch_score: int = -1,
    gap_open: int = -3,
    gap_extend: int = -2,
) -> AlignmentResult:
    """Run local Smith-Waterman alignment with affine gap penalties."""
    left = normalize_dna(seq1)
    right = normalize_dna(seq2)
    _ensure_positive("match_score", match_score)
    _ensure_not_positive("mismatch_score", mismatch_score)
    _ensure_not_positive("gap_open", gap_open)
    _ensure_not_positive("gap_extend", gap_extend)

    n, m = len(left), len(right)
    scores = [[[0, 0, 0] for _ in range(m + 1)] for _ in range(n + 1)]
    traceback: List[List[List[Optional[Tuple[int, int, int]]]]] = [
        [[None, None, None] for _ in range(m + 1)] for _ in range(n + 1)
    ]
    best_score = 0
    best_position = (0, 0, MATCH_STATE)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitution = match_score if left[i - 1] == right[j - 1] else mismatch_score
            match_candidates = [
                (0, None),
                (scores[i - 1][j - 1][MATCH_STATE] + substitution, (i - 1, j - 1, MATCH_STATE)),
                (scores[i - 1][j - 1][INSERT_STATE] + substitution, (i - 1, j - 1, INSERT_STATE)),
                (scores[i - 1][j - 1][DELETE_STATE] + substitution, (i - 1, j - 1, DELETE_STATE)),
            ]
            scores[i][j][MATCH_STATE], traceback[i][j][MATCH_STATE] = max(match_candidates, key=lambda item: item[0])

            insert_candidates = [
                (0, None),
                (scores[i][j - 1][MATCH_STATE] + gap_open, (i, j - 1, MATCH_STATE)),
                (scores[i][j - 1][INSERT_STATE] + gap_extend, (i, j - 1, INSERT_STATE)),
                (scores[i][j - 1][DELETE_STATE] + gap_open, (i, j - 1, DELETE_STATE)),
            ]
            scores[i][j][INSERT_STATE], traceback[i][j][INSERT_STATE] = max(insert_candidates, key=lambda item: item[0])

            delete_candidates = [
                (0, None),
                (scores[i - 1][j][MATCH_STATE] + gap_open, (i - 1, j, MATCH_STATE)),
                (scores[i - 1][j][INSERT_STATE] + gap_open, (i - 1, j, INSERT_STATE)),
                (scores[i - 1][j][DELETE_STATE] + gap_extend, (i - 1, j, DELETE_STATE)),
            ]
            scores[i][j][DELETE_STATE], traceback[i][j][DELETE_STATE] = max(delete_candidates, key=lambda item: item[0])

            for state in (MATCH_STATE, INSERT_STATE, DELETE_STATE):
                if scores[i][j][state] > best_score:
                    best_score = scores[i][j][state]
                    best_position = (i, j, state)

    return _trace_local_alignment(left, right, scores, traceback, best_position, best_score)


def hybrid_alignment(
    seq1: str,
    seq2: str,
    gap_penalty_consecutive: int,
    gap_penalty_isolated: int,
    match_gain: int,
    mismatch_loss: int,
) -> Tuple[int, str, str]:
    """Backward-compatible wrapper returning the notebook's original tuple shape."""
    result = smith_waterman(
        seq1,
        seq2,
        match_score=match_gain,
        mismatch_score=mismatch_loss,
        gap_open=gap_penalty_isolated,
        gap_extend=gap_penalty_consecutive,
    )
    return result.score, result.aligned_seq1, result.aligned_seq2


def smith_waterman_with_compression(seq1: str, seq2: str) -> Tuple[int, str, str]:
    """Round-trip sequence one through BWT, then align the recovered sequence.

    The BWT is a reversible preprocessing representation, not a biologically
    meaningful sequence to align directly. This function verifies the round trip
    and then aligns the original biological sequence.
    """
    original = normalize_dna(seq1)
    transformed, _ = bwt_transform(original)
    recovered = inverse_bwt(transformed)
    if recovered != original:
        raise ValueError("BWT round trip failed; refusing to align corrupted sequence")

    result = smith_waterman(recovered, seq2)
    return result.score, result.aligned_seq1, result.aligned_seq2


def _trace_local_alignment(
    seq1: str,
    seq2: str,
    scores: List[List[List[int]]],
    traceback: List[List[List[Optional[Tuple[int, int, int]]]]],
    best_position: Tuple[int, int, int],
    best_score: int,
) -> AlignmentResult:
    aligned_seq1: List[str] = []
    aligned_seq2: List[str] = []
    i, j, state = best_position
    end_seq1, end_seq2 = i, j

    while i > 0 and j > 0 and scores[i][j][state] > 0:
        previous = traceback[i][j][state]
        if previous is None:
            break
        prev_i, prev_j, prev_state = previous
        if state == MATCH_STATE:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
        elif state == INSERT_STATE:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
        else:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
        i, j, state = prev_i, prev_j, prev_state

    aligned_seq1.reverse()
    aligned_seq2.reverse()
    return AlignmentResult(
        score=int(best_score),
        aligned_seq1="".join(aligned_seq1),
        aligned_seq2="".join(aligned_seq2),
        start_seq1=i,
        end_seq1=end_seq1,
        start_seq2=j,
        end_seq2=end_seq2,
    )


def _ensure_non_negative(name: str, value: int) -> None:
    if value < 0:
        raise ValueError(f"{name} must be non-negative")


def _ensure_positive(name: str, value: int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def _ensure_not_positive(name: str, value: int) -> None:
    if value > 0:
        raise ValueError(f"{name} must be zero or negative")

