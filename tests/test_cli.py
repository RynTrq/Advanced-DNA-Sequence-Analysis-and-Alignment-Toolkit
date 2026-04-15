import io
import json
import unittest
from contextlib import redirect_stdout

from dna_toolkit.cli import main


class CLITests(unittest.TestCase):
    def test_edit_distance_command_prints_distance(self):
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["edit-distance", "AAAA", "AA", "--gap-open", "3", "--gap-extend", "1"])
        self.assertEqual(exit_code, 0)
        self.assertEqual(stdout.getvalue().strip(), "4")

    def test_align_command_prints_json_result(self):
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = main(["align", "TTACGTAA", "GGACGTCC"])
        self.assertEqual(exit_code, 0)
        payload = json.loads(stdout.getvalue())
        self.assertEqual(payload["score"], 12)
        self.assertEqual(payload["aligned_seq1"], "ACGT")


if __name__ == "__main__":
    unittest.main()

