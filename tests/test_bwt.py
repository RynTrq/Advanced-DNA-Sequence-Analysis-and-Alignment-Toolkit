import unittest

from dna_toolkit import BWTIndex, bwt_transform, find_exact_matches, inverse_bwt


class BWTTests(unittest.TestCase):
    def test_bwt_transform_uses_stable_suffix_array_for_repeats(self):
        transformed, suffix_array = bwt_transform("ATAT")
        self.assertEqual(transformed, "TT$AA")
        self.assertEqual(suffix_array, [4, 2, 0, 3, 1])

    def test_inverse_bwt_round_trip(self):
        transformed, _ = bwt_transform("GATTACA")
        self.assertEqual(inverse_bwt(transformed), "GATTACA")

    def test_exact_match_returns_sorted_zero_based_positions(self):
        matches = find_exact_matches("ATATAT", ["ATA", "TAT", "GG"])
        self.assertEqual(matches["ATA"], [0, 2])
        self.assertEqual(matches["TAT"], [1, 3])
        self.assertEqual(matches["GG"], [])

    def test_bwt_index_rejects_empty_patterns(self):
        index = BWTIndex.build("ACGT")
        with self.assertRaises(ValueError):
            index.exact_match("")

    def test_bwt_rejects_reserved_sentinel_in_input(self):
        with self.assertRaises(ValueError):
            bwt_transform("AC$GT")


if __name__ == "__main__":
    unittest.main()

