import unittest

from dna_toolkit import (
    edit_distance_with_gap_penalties,
    hybrid_alignment,
    smith_waterman,
    smith_waterman_with_compression,
)


class AlignmentTests(unittest.TestCase):
    def test_edit_distance_prefers_gap_extension_for_gap_runs(self):
        self.assertEqual(
            edit_distance_with_gap_penalties("AAAA", "AA", consecutive_gap_cost=1, isolated_gap_cost=3),
            4,
        )

    def test_edit_distance_handles_empty_sequences(self):
        self.assertEqual(edit_distance_with_gap_penalties("", "ACG", 2, 5), 9)
        self.assertEqual(edit_distance_with_gap_penalties("", "", 2, 5), 0)

    def test_edit_distance_rejects_negative_costs(self):
        with self.assertRaises(ValueError):
            edit_distance_with_gap_penalties("A", "G", -1, 2)

    def test_smith_waterman_traces_from_best_scoring_cell(self):
        result = smith_waterman("TTACGTAA", "GGACGTCC")
        self.assertEqual(result.score, 12)
        self.assertEqual(result.aligned_seq1, "ACGT")
        self.assertEqual(result.aligned_seq2, "ACGT")
        self.assertEqual((result.start_seq1, result.end_seq1), (2, 6))

    def test_smith_waterman_supports_affine_gap_alignment(self):
        result = smith_waterman("ACGTT", "ACTT", match_score=3, mismatch_score=-4, gap_open=-2, gap_extend=-1)
        self.assertEqual(result.aligned_seq1, "ACGTT")
        self.assertEqual(result.aligned_seq2, "AC-TT")
        self.assertEqual(result.score, 10)

    def test_hybrid_alignment_preserves_original_tuple_api(self):
        score, aligned_seq1, aligned_seq2 = hybrid_alignment("ACGT", "ACCT", -2, -3, 3, -1)
        self.assertEqual(score, 8)
        self.assertEqual(aligned_seq1, "ACGT")
        self.assertEqual(aligned_seq2, "ACCT")

    def test_compressed_alignment_round_trips_before_alignment(self):
        score, aligned_seq1, aligned_seq2 = smith_waterman_with_compression("TTACGTAA", "GGACGTCC")
        self.assertEqual(score, 12)
        self.assertEqual(aligned_seq1, "ACGT")
        self.assertEqual(aligned_seq2, "ACGT")


if __name__ == "__main__":
    unittest.main()

