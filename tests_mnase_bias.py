from mnase_bias import *
import unittest

class Test_nucleotide_comp(unittest.TestCase):

    def test_map_kmers_to_positions(self):
        fasta = "tests/map_kmers_to_positions_input.fasta"
        expected = {"AT": [("ENSMUS1", 1), ("ENSMUS1", 5)], "TG": [("ENSMUS1", 2), ("ENSMUS1", 6), ("ENSMUS2", 5)],
                    "GC": [("ENSMUS1", 3)], "CA": [("ENSMUS1", 4)], "AG": [("ENSMUS2", 1)],
                    "GG": [("ENSMUS2", 2), ("ENSMUS2", 3)], "GT": [("ENSMUS2", 4)]}
        observed = map_kmers_to_positions(fasta)
        self.assertEqual(expected, observed)

    def test_map_kmers_to_positions_4mer(self):
        fasta = "tests/map_kmers_to_positions_input.fasta"
        expected = {"ATGC": [("ENSMUS1", 2)], "TGCA": [("ENSMUS1", 3)], "GCAT": [("ENSMUS1", 4)],
                    "CATG": [("ENSMUS1", 5)], "AGGG": [("ENSMUS2", 2)], "GGGT": [("ENSMUS2", 3)],
                    "GGTG": [("ENSMUS2", 4)]}
        observed = map_kmers_to_positions(fasta, k=4, focal_pos=2)
        self.assertEqual(expected, observed)

    def test_transpose_to_new_coords_plus(self):
        transcripts_dict = {}
        transcripts_dict["trans1"] = ["chr1", "3", "14", "trans1", ".", "+"]
        transcripts_dict["trans2"] = ["chr2", "1", "7", "trans2", ".", "-"]
        random_pos = ("trans1", 5)
        length = 4
        expected = "chr1", 8, 12, "+"
        observed = transpose_to_new_coords(transcripts_dict, random_pos, length)
        self.assertEqual(expected, observed)

    def test_transpose_to_new_coords_minus(self):
        transcripts_dict = {}
        transcripts_dict["trans1"] = ["chr1", "3", "8", "trans1", ".", "+"]
        transcripts_dict["trans2"] = ["chr2", "1", "7", "trans2", ".", "-"]
        random_pos = ("trans2", 2)
        length = 3
        expected = "chr2", 2, 5, "-"
        observed = transpose_to_new_coords(transcripts_dict, random_pos, length)
        self.assertEqual(expected, observed)