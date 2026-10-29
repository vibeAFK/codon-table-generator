"""Unit tests for the codon table generator. Run with: python -m unittest -v"""

import unittest

from codon_table import (CODON_TABLE, validate_sequence, split_codons,
                         translate, mark_codons, get_leftover)


def run(dna):
    """Helper: translate a clean DNA string and return the marked results."""
    return mark_codons(translate(split_codons(dna)))


class TestCodonTable(unittest.TestCase):
    def test_sample_input(self):
        results = run("ATGGCTTTTTAA")
        self.assertEqual([r[1] for r in results], ["ATG", "GCT", "TTT", "TAA"])
        self.assertEqual([r[2] for r in results],
                         ["Methionine", "Alanine", "Phenylalanine", "Stop"])
        self.assertEqual(results[0][3], "Start codon")
        self.assertEqual(results[3][3], "Stop codon")

    def test_lowercase_input(self):
        self.assertEqual(validate_sequence("atggct"), ("ATGGCT", "DNA"))

    def test_invalid_character(self):
        with self.assertRaises(ValueError):
            validate_sequence("ATGXCT")

    def test_empty_input(self):
        with self.assertRaises(ValueError):
            validate_sequence("   ")

    def test_mixed_t_and_u(self):
        with self.assertRaises(ValueError):
            validate_sequence("ATGUCT")

    def test_incomplete_codon(self):
        self.assertEqual(split_codons("ATGGC"), ["ATG"])
        self.assertEqual(get_leftover("ATGGC"), "GC")
        self.assertEqual(get_leftover("ATGGCT"), "")

    def test_rna_input(self):
        seq, mode = validate_sequence("AUGGCUUUUUAA")
        self.assertEqual(mode, "RNA")
        results = run(seq.replace("U", "T"))
        self.assertEqual([r[2] for r in results],
                         ["Methionine", "Alanine", "Phenylalanine", "Stop"])

    def test_64_codons(self):
        self.assertEqual(len(CODON_TABLE), 64)

    def test_all_stop_codons(self):
        for codon in ("TAA", "TAG", "TGA"):
            self.assertEqual(CODON_TABLE[codon], "Stop")


if __name__ == "__main__":
    unittest.main()
