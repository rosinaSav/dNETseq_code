from splice_distance import *
import numpy as np
import unittest

class Test_splice_distance(unittest.TestCase):

    def test_make_dist_mat(self):
        distances = {"intron1": [4, 2, 4, 4, 3, 6, 8, 7, -1, -2, -3, -4],
                     "intron2": [1, 3, 2, 3, 5],
                     "intron3": [2, -1]}
        lengths = {"intron1": 9,
                     "intron2": 7,
                     "intron3": 4}
        max_dist = 6
        expected = (["intron1", "intron2", "intron3"], np.mat([[0, 0, 1, 1, 3, 0], [0, 1, 1, 2, 0, 1], [0, 0, 1, 0, np.nan, np.nan]]))
        observed = make_dist_mat(distances, max_dist, lengths)
        self.assertEqual(expected[0], observed[0])
        self.assertEqual(str(expected[1]), str(observed[1]))

    def test_write_exon_starts(self):
        valid_junctions = ["ENST1.0", "ENST1.1", "ENST2.0", "ENST2.1"]
        exons = {}
        exons["ENST1"] = [['X', 'FlyBase', 'CDS', 2, 5, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 16, 19, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";']]
        exons["ENST2"] = [['X', 'FlyBase', 'CDS', 16, 19, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 2, 5, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";']]
        file_name = "tests/write_exon_starts_observed.bed"
        limit = 2
        write_exon_starts(valid_junctions, file_name, exons, limit, add_chr=False, from_end=False)
        expected = rw.read_as_string("tests/write_exon_starts_expected.bed")
        observed = rw.read_as_string(file_name)
        self.assertEqual(expected, observed)

    def test_write_intron_starts(self):
        valid_junctions = ["ENST1.0", "ENST1.1", "ENST2.0", "ENST2.1"]
        exons = {}
        exons["ENST1"] = [['X', 'FlyBase', 'CDS', 2, 5, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 16, 19, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";']]
        exons["ENST2"] = [['X', 'FlyBase', 'CDS', 16, 19, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 2, 5, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";']]
        file_name = "tests/write_intron_starts_observed.bed"
        limit = 2
        write_intron_starts(valid_junctions, file_name, exons, limit, add_chr=False)
        expected = rw.read_as_string("tests/write_intron_starts_expected.bed")
        observed = rw.read_as_string(file_name)
        self.assertEqual(expected, observed)

    def test_write_si_pos(self):
        valid_junctions = ["ENST1.0", "ENST1.1", "ENST2.0", "ENST2.1"]
        exons = {}
        exons["ENST1"] = [['X', 'FlyBase', 'CDS', 2, 5, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";'],
                          ['X', 'FlyBase', 'CDS', 16, 19, '.', '+', '0','gene_id "FBgn1"; gene_symbol "Nep3"; transcript_id "FBtr0307555"; transcript_symbol "Nep3-RC";']]
        exons["ENST2"] = [['X', 'FlyBase', 'CDS', 16, 19, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 9, 11, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";'],
                          ['X', 'FlyBase', 'CDS', 2, 5, '.', '-', '0','gene_id "FBgn2"; gene_symbol "Nep4"; transcript_id "FBtr0307555"; transcript_symbol "Nep4-RC";']]
        file_name = "tests/write_si_pos_observed.bed"
        write_si_pos(valid_junctions, file_name, exons, add_chr=False)
        expected = rw.read_as_string("tests/write_si_pos_expected.bed")
        observed = rw.read_as_string(file_name)
        self.assertEqual(expected, observed)

    def test_write_read_lengths(self):
        input_file = "tests/write_read_lengths_input.bed"
        expected_file = "tests/write_read_lengths_expected.txt"
        with open(expected_file) as file:
            expected = file.readlines()
        observed_file = "tests/write_read_lengths_observed.txt"
        write_read_lengths(input_file, observed_file)
        with open(observed_file) as file:
            observed = file.readlines()
        self.assertEqual(expected, observed)

