from nucleotide_comp import *
import numpy as np
import read_and_write as rw
import unittest

class Test_nucleotide_comp(unittest.TestCase):

    def test_calc_nt_freqs1(self):
        motifs = ["AATGGA", "GAGGAG", "AGAGTAA"]
        expected = {"A": 9/19, "T": 2/19, "C": 0/19, "G": 8/19}
        observed = calc_nt_freqs(motifs, 1)
        self.assertEqual(observed,expected)

    def test_calc_nt_freqs2(self):
        motifs = ["AATGGA", "GAGGAG", "AGAGTAA"]
        expected = {"AA": 2/16, "AT": 1/16, "AC": 0/16, "AG": 4/16,
                    "TA": 1/16, "TT": 0/16, "TC": 0/16, "TG": 1/16,
                    "CA": 0/16, "CT": 0/16, "CC": 0/16, "CG": 0/16,
                    "GA": 4/16, "GT": 1/16, "GC": 0/16, "GG": 2/16}
        observed = calc_nt_freqs(motifs, 2)
        self.assertEqual(observed,expected)

    def test_calc_nt_freqs3(self):
        motifs = ["AATGGA", "GAGGAG", "AGAGTAA"]
        expected = {}
        observed = calc_nt_freqs(motifs, 3)
        self.assertEqual(observed["AAT"], 1/13)
        self.assertEqual(observed["CCC"], 0/13)
        self.assertEqual(observed["GGA"], 2/13)
        self.assertEqual(observed["AGG"], 1/13)
        self.assertEqual(observed["GTA"], 1/13)
        self.assertEqual(observed["CGT"], 0/13)
        self.assertEqual(observed["GAG"], 3/13)
        self.assertEqual(observed["TAA"], 1/13)

    def test_get_4fold_deg(self):
        sequence = "CCTCCTAAT"
        expected = [2, 5]
        observed = get_4fold_deg(sequence)
        self.assertEqual(expected, observed)

    def test_get_GC(self):
        sequence = "ATTCGGACNNTTTACG"
        expected = 6/14
        observed = get_GC(sequence)
        self.assertEqual(observed, expected)

    def test_get_GC_2(self):
        sequence = "ATTCGGACNNTTTACG"
        expected = 8/14
        observed = get_GC(sequence, alternative_bases = ["A", "T"])
        self.assertEqual(observed, expected)

    def test_get_GC_3(self):
        sequence = "AttCGGACNNTttACg"
        expected = 8/14
        observed = get_GC(sequence, alternative_bases = ["A", "T"])
        self.assertEqual(observed, expected)

    def test_get_GC4(self):
        names, sequences = rw.read_fasta("tests/test_get_GC4_input.fasta")
        expected = [1.0, 1.0, 1.0, 0.5, 2/3]
        phases = [2, 0, 0, 1, 0]
        observed = [get_GC4(sequences[i], phases[i]) for i in range(len(sequences))]
        self.assertEqual(observed, expected)

    def test_get_ss_strength_upstream_5(self):
        genome = "tests/new_test.fa"
        exons = {"FBtr1": [["1", "Flybase", "exon", 6, 9, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 29, 32, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"]],
                "FBtr2": [["1", "Flybase", "exon", 29, 32, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"]]}
        observed = get_ss_strength(exons, genome, upstream = True, five = True, exonic = 2, intronic = 7)
        expected = {}
        expected["FBtr1.0"] = -38.41
        expected["FBtr1.1"] = -13.74
        expected["FBtr1.2"] = -32.26
        expected["FBtr2.0"] = -1.42
        expected["FBtr2.1"] = -25.66
        self.assertEqual(expected, observed)

    def test_get_ss_strength_upstream_3(self):
        exons = {"FBtr1": [["1", "Flybase", "exon", 6, 9, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 29, 32, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"]],
                "FBtr2": [["1", "Flybase", "exon", 29, 32, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"]]}
        genome = "tests/new_test.fa"
        observed = get_ss_strength(exons, genome, upstream = True, five = False, exonic = 16, intronic = 7)
        expected = {}
        expected["FBtr1.0"] = -0.49
        expected["FBtr1.1"] = -23.16
        expected["FBtr1.2"] = -13.78
        expected["FBtr2.0"] = -16.98
        expected["FBtr2.1"] = -11.01
        self.assertEqual(expected, observed)

    def test_get_ss_strength_downstream_5(self):
        exons = {"FBtr1": [["1", "Flybase", "exon", 6, 9, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 29, 32, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"]],
                "FBtr2": [["1", "Flybase", "exon", 29, 32, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 18, ".", "-", ".", "gene_id \"FBgn2\"; gene_symbol \"Nep3\"; transcript_id \"FBtr2\"; transcript_symbol \"Nep3-RA\";"]]}
        genome = "tests/new_test.fa"
        observed = get_ss_strength(exons, genome, upstream = False, five = True, exonic = 2, intronic = 7)
        expected = {}
        expected["FBtr1.0"] = -13.74
        expected["FBtr1.1"] = -32.26
        expected["FBtr2.0"] = -25.66
        self.assertEqual(expected, observed)

    def test_get_ss_strength_downstream_3(self):
        exons = {"FBtr1": [["1", "Flybase", "exon", 6, 9, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 15, 19, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 22, 24, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"],
                                ["1", "Flybase", "exon", 29, 32, ".", "+", ".", "gene_id \"FBgn1\"; gene_symbol \"Nep3\"; transcript_id \"FBtr1\"; transcript_symbol \"Nep3-RA\";"]]}
        genome = "tests/new_test.fa"
        observed = get_ss_strength(exons, genome, upstream = False, five = False, exonic = 16, intronic = 7)
        expected = {}
        expected["FBtr1.0"] = -23.16
        expected["FBtr1.1"] = -13.78
        self.assertEqual(expected, observed)

    def test_kmer_enrichment_test(self):
        probs = [1/24, 1/6, 1/3]
        expected = [(0.125, 8, 0.1199), (0.5, 2, 0.4213), (1, 1, 1)]
        observed = [kmer_enrichment_test(1, probs[i], 3) for i in range(len(probs))]
        for pos, test in enumerate(observed):
            curr_exp = expected[pos]
            [self.assertAlmostEqual(curr_exp[i], test[i], places=4) for i in range(len(test))]

    def test_make_PPM(self):
        occ_mat = np.array([["A", "G", "G", "A"], ["A", "A", "G", "A"], ["A", "T", "G", "T"], ["A", "C", "N", "T"], ["N", "N", "N", "N"]])
        expected = np.array([[1.0, 0.0, 0.0, 0.0], [0.25, 0.25, 0.25, 0.25], [0.0, 0.0, 0.0, 1.0], [0.5, 0.5, 0.0, 0.0]])
        bases = ["A", "T", "C", "G"]
        observed = make_PPM(occ_mat, bases)
        expected = np.ndarray.tolist(expected)
        observed = np.ndarray.tolist(observed)
        self.assertEqual(expected, observed)

    def test_markov_coding(self):
        seq_frags = ["ATTATA", "TTTAAA", "AAATTT"]
        bases = ["A", "T"]
        expected_freqs = {"AA": 2/6, "AT": 2/6, "TA": 0, "TT": 2/6}
        expected_trans_probs = {"AA": {"A": {0: 0, 1: 0.5, 2: 1}, "T": {0: 1, 1: 0.5, 2: 0}},
                                "AT": {"A": {0: 0.5, 1: 0, 2: 1/2}, "T": {0: 0.5, 1: 1, 2: 1/2}},
                                "TA": {"A": {0: 0.5, 1: 0.5, 2: 0.5}, "T": {0: 0.5, 1: 0.5, 2: 0.5}},
                                "TT": {"A": {0: 1, 1: 0.5, 2: 0}, "T": {0: 0, 1: 0.5, 2: 1}}}
        observed = markov(seq_frags, bases, 2, coding=True)
        expected = expected_freqs, expected_trans_probs
        self.assertEqual(expected, observed)

    def test_markov_coding_dints(self):
        seq_frags = ["ATTATA", "TTTAAA", "AAATTT"]
        bases = ["A", "T"]
        expected_freqs = {"A": 2/3, "T": 1/3}
        expected_trans_probs = {"A": {"A": {0: 0, 1: 0.5, 2: 1}, "T": {0: 1, 1: 0.5, 2: 0}},
                                "T": {"A": {0: 1, 1: 0, 2: 0.25}, "T": {0: 0, 1: 1, 2: 0.75}}}
        observed = markov(seq_frags, bases, 1, coding=True)
        expected = expected_freqs, expected_trans_probs
        self.assertEqual(expected, observed)

    def test_markov_nc(self):
        seq_frags = ["ATTATA", "TTTAAA", "AAATTT"]
        bases = ["A", "T"]
        expected_freqs = {"AA": 4/15, "AT": 1/5, "TA": 1/5, "TT": 1/3}
        expected_trans_probs = {"AA": {"A": 2/3, "T": 1/3}, "AT": {"A": 1/3, "T": 2/3}, "TA": {"A": 1/2, "T": 1/2}, "TT": {"A": 1/2, "T": 1/2}}
        observed = markov(seq_frags, bases, 2, coding=False)
        expected = expected_freqs, expected_trans_probs
        self.assertEqual(expected, observed)

    def test_markov_prob_coding(self):
        seq_frags = ["ATTATA", "TTTAAA", "AAATTT"]
        dint_freqs = {"AA": 2/6, "AT": 2/6, "TA": 0, "TT": 2/6}
        trans_probs = {"AA": {"A": {0: 0, 1: 0.5, 2: 1}, "T": {0: 1, 1: 0.5, 2: 0}},
                                "AT": {"A": {0: 0.5, 1: 0, 2: 1/2}, "T": {0: 0.5, 1: 1, 2: 1/2}},
                                "TA": {"A": {0: 0.5, 1: 0.5, 2: 0.5}, "T": {0: 0.5, 1: 0.5, 2: 0.5}},
                                "TT": {"A": {0: 1, 1: 0.5, 2: 0}, "T": {0: 0, 1: 0.5, 2: 1}}}
        expected = [1/24, 1/6, 1/3]
        observed = [markov_prob(i, dint_freqs, trans_probs, 2, coding=True) for i in seq_frags]
        [self.assertAlmostEqual(expected[i], observed[i]) for i in range(len(observed))]

    def test_markov_prob_nc(self):
        seq_frags = ["ATTATA", "TTTAAA", "AAATTT"]
        dint_freqs = {"AA": 4/15, "AT": 1/5, "TA": 1/5, "TT": 1/3}
        trans_probs = {"AA": {"A": 2/3, "T": 1/3}, "AT": {"A": 1/3, "T": 2/3}, "TA": {"A": 1/2, "T": 1/2}, "TT": {"A": 1/2, "T": 1/2}}
        expected = [1/90, 1/36, 8/405]
        observed = [markov_prob(i, dint_freqs, trans_probs, 2, coding=False) for i in seq_frags]
        [self.assertAlmostEqual(expected[i], observed[i]) for i in range(len(observed))]

    def test_PPM_diff(self):
        PPM1 = np.array([[0.1, 0.0, 0.4, 0.5], [0.0, 0.0, 0.0, 0.1], [0.2, 0.3, 0.2, 0.3]])
        PPM2 = np.array([[0, 0.3, 0.5, 0.2], [0.0, 0.0, 0.0, 0.1], [0.3, 0.2, 0.2, 0.3]])
        expected = 0.22
        observed = PPM_diff(PPM1, PPM2)
        self.assertAlmostEqual(expected, observed)