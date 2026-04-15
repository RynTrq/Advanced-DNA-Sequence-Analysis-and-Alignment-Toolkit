"""Command-line interface for the DNA toolkit."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Sequence

from .alignment import edit_distance_with_gap_penalties, smith_waterman, smith_waterman_with_compression
from .bwt import BWTIndex, bwt_transform, find_exact_matches, inverse_bwt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DNA sequence analysis and alignment toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    edit = subparsers.add_parser("edit-distance", help="Compute edit distance with affine gap costs")
    edit.add_argument("seq1")
    edit.add_argument("seq2")
    edit.add_argument("--gap-open", type=int, default=3, help="Cost for the first gap in a run")
    edit.add_argument("--gap-extend", type=int, default=2, help="Cost for each consecutive gap")
    edit.add_argument("--mismatch", type=int, default=1, help="Mismatch substitution cost")

    align = subparsers.add_parser("align", help="Run local Smith-Waterman alignment")
    align.add_argument("seq1")
    align.add_argument("seq2")
    align.add_argument("--match", type=int, default=3)
    align.add_argument("--mismatch", type=int, default=-1)
    align.add_argument("--gap-open", type=int, default=-3)
    align.add_argument("--gap-extend", type=int, default=-2)

    bwt = subparsers.add_parser("bwt", help="Build a BWT representation")
    bwt.add_argument("sequence")

    inverse = subparsers.add_parser("inverse-bwt", help="Invert a BWT string")
    inverse.add_argument("transformed")

    match = subparsers.add_parser("match", help="Find exact pattern matches using a BWT index")
    match.add_argument("sequence")
    match.add_argument("patterns", nargs="+")

    compressed = subparsers.add_parser("compressed-align", help="BWT round-trip followed by alignment")
    compressed.add_argument("seq1")
    compressed.add_argument("seq2")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "edit-distance":
        distance = edit_distance_with_gap_penalties(
            args.seq1,
            args.seq2,
            consecutive_gap_cost=args.gap_extend,
            isolated_gap_cost=args.gap_open,
            mismatch_cost=args.mismatch,
        )
        print(distance)
        return 0

    if args.command == "align":
        result = smith_waterman(
            args.seq1,
            args.seq2,
            match_score=args.match,
            mismatch_score=args.mismatch,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
        )
        print(json.dumps(asdict(result), indent=2))
        return 0

    if args.command == "bwt":
        transformed, suffix_array = bwt_transform(args.sequence)
        print(json.dumps({"bwt": transformed, "suffix_array": suffix_array}, indent=2))
        return 0

    if args.command == "inverse-bwt":
        print(inverse_bwt(args.transformed))
        return 0

    if args.command == "match":
        print(json.dumps(find_exact_matches(args.sequence, args.patterns), indent=2))
        return 0

    if args.command == "compressed-align":
        score, aligned_seq1, aligned_seq2 = smith_waterman_with_compression(args.seq1, args.seq2)
        print(json.dumps({"score": score, "aligned_seq1": aligned_seq1, "aligned_seq2": aligned_seq2}, indent=2))
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

