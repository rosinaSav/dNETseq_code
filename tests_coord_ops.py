from coord_ops import *
from pyfaidx import Fasta
import read_and_write as rw
import unittest

class Test_coord_ops(unittest.TestCase):

    def test_extend_intervals(self):
        input_file = "tests/extend_intervals_input.bed"
        expected = "tests/extend_intervals_expected.bed"
        observed = "tests/extend_intervals_observed.bed"
        left_shift = 4
        right_shift = 5
        extend_intervals(input_file, observed, left_shift, right_shift)
        expected = rw.read_as_string(expected)
        observed = rw.read_as_string(observed)
        self.assertEqual(expected, observed)

    def test_extend_intervals_three_prime(self):
        input_file = "tests/extend_intervals_input.bed"
        expected = "tests/extend_intervals_three_prime_expected.bed"
        observed = "tests/extend_intervals_three_prime_observed.bed"
        left_shift = 4
        right_shift = 2
        extend_intervals(input_file, observed, left_shift, right_shift, three_prime=True)
        expected = rw.read_as_string(expected)
        observed = rw.read_as_string(observed)
        self.assertEqual(expected, observed)

    def test_get_exon_and_intron_coords_plus_with_ups_intron(self):
        exons = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]
        expected = [(12, 20), (20, 38), (38, 70)]
        observed = get_exon_and_intron_coords(exons, with_ups_intron=True)
        self.assertEqual(expected, observed)

    def test_get_exon_and_intron_coords_minus_with_ups_intron(self):
        exons = [["2L", "FlyBase", "exon", 22, 23, ".", "-", ".",
                        "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 18, 19, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 13, 14, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 4, 7, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"]]
        expected = [(3, 12), (12, 17), (17, 21), (21, 23)]
        observed = get_exon_and_intron_coords(exons, with_ups_intron=True)
        self.assertEqual(expected, observed)

    def test_get_exon_and_intron_coords_plus(self):
        exons = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]
        expected = [(12, 20), (20, 29), (29, 38), (38, 59), (59, 70)]
        observed = get_exon_and_intron_coords(exons)
        self.assertEqual(expected, observed)

    def test_get_exon_and_intron_coords_minus(self):
        exons = [["2L", "FlyBase", "exon", 22, 23, ".", "-", ".",
                        "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 18, 19, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 13, 14, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 4, 7, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"]]
        expected = [(3, 7), (7, 12), (12, 14), (14, 17), (17, 19), (19, 21), (21, 23)]
        observed = get_exon_and_intron_coords(exons)
        self.assertEqual(expected, observed)

    def test_get_exon_and_intron_coords_singleexon(self):
        exons = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]
        expected = [(12, 20)]
        observed = get_exon_and_intron_coords(exons)
        self.assertEqual(expected, observed)

    def test_get_exon_number(self):
        input_dict = {}
        input_dict["FBgn1"] = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        input_dict["FBgn2"] = [["2L", "FlyBase", "exon", 14, 350, ".", "+", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]

        input_dict["FBgn3"] = [["2L", "FlyBase", "exon", 14, 350, ".", "+", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]
        valid_junctions = ["FBgn3.+.1", "FBgn1.+.0", "FBgn1.+.2"]
        expected = {"FBgn3.+.1": 1, "FBgn1.+.0": 3, "FBgn1.+.2": 3}
        observed = get_exon_number(input_dict, valid_junctions)
        self.assertEqual(expected, observed)


    def test_extract_3ss(self):
        exons = {}
        exons["FBgn1"] = [["2L", "FlyBase", "exon", 1, 4, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 9, 10, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 12, 14, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 17, 18, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 21, 23, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        exons["FBgn2"] = [["2L", "FlyBase", "exon", 22, 23, ".", "-", ".",
                        "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 18, 19, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 13, 14, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 4, 7, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"]]

        output_file = "tests/extract_3ss_observed.bed"
        expected = rw.read_as_string("tests/extract_3ss_expected.bed")
        extract_3ss(exons, output_file)
        observed = rw.read_as_string(output_file)
        self.assertEqual(expected, observed)

    def test_extract_exon_junctions(self):
        exons = {}
        exons["FBgn1"] = [["2L", "FlyBase", "exon", 1, 4, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 9, 10, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 12, 14, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 17, 18, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 21, 23, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        exons["FBgn2"] = [["2L", "FlyBase", "exon", 22, 23, ".", "-", ".",
                        "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 18, 19, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 13, 14, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 4, 7, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"]]

        output_file = "tests/extract_exon_junctions_observed.bed"
        with open("tests/extract_exon_junctions_expected.bed") as file:
            expected = file.read()
        extract_exon_junctions(exons, output_file)
        with open(output_file) as file:
            observed = file.read()
        self.assertEqual(expected, observed)

    def test_get_exon_rank(self):
        exons = {}
        exons["FBgn1"] = [["2L", "FlyBase", "exon", 1, 4, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 9, 10, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 12, 14, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 17, 18, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 21, 23, ".", "+", ".",
                           "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        exons["FBgn2"] = [["2L", "FlyBase", "exon", 22, 23, ".", "-", ".",
                        "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 18, 19, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 13, 14, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"],
                        ["2L", "FlyBase", "exon", 4, 7, ".", "-", ".",
                       "gene_id \"FBgn2\"; gene_symbol \"drl\"; transcript_id \"FBtr3\"; transcript_symbol \"drl-RA\";"]]
        exon_starts = {}
        exon_starts["FBgn2.0"] = ["chr2L", "17", "19", "FBgn2.0", ".", "-"]
        exon_starts["FBgn2.1"] = ["chr2L", "12", "14", "FBgn2.1", ".", "-"]
        exon_starts["FBgn1.1"] = ["chr2L", "11", "14", "FBgn1.1", ".", "+"]
        exon_starts["FBgn1.3"] = ["chr2L", "20", "23", "FBgn1.1", ".", "+"]
        expected_starts = {"FBgn2.0": 1, "FBgn2.1": 2, "FBgn1.1": 2, "FBgn1.3": 4}
        expected_ends = {"FBgn2.0": 2, "FBgn2.1": 1, "FBgn1.1": 2, "FBgn1.3": 0}
        expected = expected_starts, expected_ends
        observed = get_exon_rank(exons, exon_starts)
        self.assertEqual(expected, observed)

    def test_get_flanking_intron_sizes(self):
        exons = rw.read_gtf("tests/get_flanking_intron_sizes_input.gtf", "exon", gene = False)
        expected = {}
        expected["ENSMUST1"] = {"upstream": [None, 4, 4], "downstream": [4, 4, None]}
        expected["ENSMUST8"] = None
        expected["ENSMUST4"] = {"upstream": [None, 1, 4, 2], "downstream": [1, 4, 2, None]}
        observed = get_flanking_intron_sizes(exons)
        self.assertEqual(observed,expected)

    def test_get_introns(self):
        input_dict = {}
        input_dict["FBgn1"] = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        input_dict["FBgn2"] = [["2L", "FlyBase", "exon", 14, 350, ".", "+", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]

        expected = {}
        expected["FBgn1"] = [["2L", "FlyBase", "intron", 21, 29, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "intron", 39, 59, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        observed = get_introns(input_dict)
        self.assertEqual(expected, observed)

    def test_get_introns_antisense(self):
        input_dict = {}
        input_dict["FBgn1"] = [["2L", "FlyBase", "exon", 60, 70, ".", "-", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "-", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 13, 20, ".", "-", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        input_dict["FBgn2"] = [["2L", "FlyBase", "exon", 14, 350, ".", "-", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]

        expected = {}
        expected["FBgn1"] = [["2L", "FlyBase", "intron", 39, 59, ".", "-", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "intron", 21, 29, ".", "-", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        observed = get_introns(input_dict)
        self.assertEqual(expected, observed)

    def test_get_lengths(self):
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 14, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        valid_junctions = ["FBtr0077982.0", "FBtr0078168.0", "FBtr0078168.1"]
        expected = {"FBtr0077982.0": 6, "FBtr0078168.0": 5, "FBtr0078168.1": 4}
        observed = get_lengths(exons, valid_junctions)
        self.assertEqual(expected, observed)

    def test_get_lengths_intronic(self):
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 14, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        valid_junctions = ["FBtr0077982.0", "FBtr0078168.0", "FBtr0078168.1"]
        expected = {"FBtr0077982.0": 3, "FBtr0078168.0": 3, "FBtr0078168.1": 3}
        observed = get_lengths(exons, valid_junctions, intronic = True)
        self.assertEqual(expected, observed)

    def test_get_sequence(self):
        expected = "GGATGGTGTGGTAG"
        coords = ["2", "havana", "exon", 14, 27, ".", "+", ".", "stuff"]
        with Fasta("tests/test_genome.fa") as pf_genome:
            observed = get_sequence(coords, pf_genome)
        self.assertEqual(observed, expected)

    def test_get_sequence_concat(self):
        expected = "GGATGGTGTGGTAGGGATGGTGTGGTAG"
        coords = [["2", "havana", "exon", 14, 27, ".", "+", ".", "stuff"],
                  ["2", "havana", "exon", 14, 27, ".", "+", ".", "stuff"]]
        with Fasta("tests/test_genome.fa") as pf_genome:
            observed = get_sequence(coords, pf_genome)
        self.assertEqual(observed, expected)

    def test_get_sequence_reverse(self):
        expected = "CTACCACACCATCC"
        coords = ["2", "havana", "exon", 14, 27, ".", "-", ".", "stuff"]
        with Fasta("tests/test_genome.fa") as pf_genome:
            observed = get_sequence(coords, pf_genome, impose_strand = True)
        self.assertEqual(observed, expected)

    def test_get_sequence_reverse_concat(self):
        expected = "CTACCACACCATCCCTAC"
        coords = [["2", "havana", "exon", 14, 27, ".", "-", ".", "stuff"],
                  ["2", "havana", "exon", 24, 27, ".", "-", ".", "stuff"]]
        with Fasta("tests/test_genome.fa") as pf_genome:
            observed = get_sequence(coords, pf_genome, impose_strand = True)
        self.assertEqual(observed, expected)

    def test_get_sequence_bed(self):
        expected = "GGATGGTGTGGTAG"
        coords = ["chr2", 13, 27, "ENSMUST7", ".", "+"]
        with Fasta("tests/test_genome.fa") as pf_genome:
            observed = get_sequence(coords, pf_genome, bed_input = True, strip_chr = True)
        self.assertEqual(observed,expected)

    def test_get_transcripts(self):
        gtf = "tests/get_transcripts_input.gtf"
        observed_file = "tests/get_transcripts_observed.bed"
        expected_file = "tests/get_transcripts_expected.bed"
        get_transcripts(gtf, observed_file)
        expected = rw.read_as_string(expected_file)
        observed = rw.read_as_string(observed_file)
        self.assertEqual(expected, observed)

    def test_get_upstream_intron_size(self):
        exons = rw.read_gtf("tests/get_upstream_intron_size_input.gtf", "exon", gene = False)
        exon_ranks = {"ENSMUST1.0": 0, "ENSMUST1.1": 1, "ENSMUST4.0": 0, "ENSMUST4.2": 2}
        expected = {"ENSMUST1.0": None, "ENSMUST1.1": 4, "ENSMUST4.0": None, "ENSMUST4.2": 4}
        observed = get_upstream_intron_size(exons, exon_ranks)
        self.assertEqual(expected, observed)
        
    def test_intersect_bed_default(self):
        A_file = "tests/intersect_bed_default_input_A.bed"
        B_file = "tests/intersect_bed_default_input_B.bed"
        expected_file = "tests/intersect_bed_default_expected.bed"
        observed_file = "tests/intersect_bed_default_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_default_bedops(self):
        A_file = "tests/intersect_bed_default_bedops_input_A.bed"
        B_file = "tests/intersect_bed_default_bedops_input_B.bed"
        expected_file = "tests/intersect_bed_default_bedops_expected.bed"
        observed_file = "tests/intersect_bed_default_bedops_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, use_bedops = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_no_dups(self):
        A_file = "tests/intersect_bed_no_dups_input_A.bed"
        B_file = "tests/intersect_bed_no_dups_input_B.bed"
        expected_file = "tests/intersect_bed_no_dups_expected.bed"
        observed_file = "tests/intersect_bed_no_dups_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count(self):
        A_file = "tests/intersect_bed_hit_count_input_A.bed"
        B_file = "tests/intersect_bed_hit_count_input_B.bed"
        expected_file = "tests/intersect_bed_hit_count_expected.bed"
        observed_file = "tests/intersect_bed_hit_count_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count_unsorted(self):
        A_file = "tests/intersect_bed_hit_count_unsorted_input_A.bed"
        B_file = "tests/intersect_bed_hit_count_unsorted_input_B.bed"
        expected_file = "tests/intersect_bed_hit_count_unsorted_expected.bed"
        observed_file = "tests/intersect_bed_hit_count_unsorted_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_overlap(self):
        A_file = "tests/intersect_bed_overlap_input_A.bed"
        B_file = "tests/intersect_bed_overlap_input_B.bed"
        expected_file = "tests/intersect_bed_overlap_expected.bed"
        observed_file = "tests/intersect_bed_overlap_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, overlap = 0.5)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand(self):
        A_file = "tests/intersect_bed_force_strand_input_A.bed"
        B_file = "tests/intersect_bed_force_strand_input_B.bed"
        expected_file = "tests/intersect_bed_force_strand_expected.bed"
        observed_file = "tests/intersect_bed_force_strand_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand_hit_count(self):
        A_file = "tests/intersect_bed_force_strand_hit_count_input_A.bed"
        B_file = "tests/intersect_bed_force_strand_hit_count_input_B.bed"
        expected_file = "tests/intersect_bed_force_strand_hit_count_expected.bed"
        observed_file = "tests/intersect_bed_force_strand_hit_count_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True, hit_count = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_write_both(self):
        A_file = "tests/intersect_write_both_input_A.bed"
        B_file = "tests/intersect_write_both_input_B.bed"
        expected_file = "tests/intersect_bed_write_both_expected.bed"
        observed_file = "tests/intersect_bed_write_both_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, write_both = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_intersect(self):
        A_file = "tests/intersect_bed_intersect_input_A.bed"
        B_file = "tests/intersect_bed_intersect_input_B.bed"
        expected_file = "tests/intersect_bed_intersect_expected.bed"
        observed_file = "tests/intersect_bed_intersect_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, intersect = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_intersect_bedops(self):
        A_file = "tests/intersect_bed_intersect_bedops_input_A.bed"
        B_file = "tests/intersect_bed_intersect_bedops_input_B.bed"
        expected_file = "tests/intersect_bed_intersect_bedops_expected.bed"
        observed_file = "tests/intersect_bed_intersect_bedops_observed.bed"
        hk.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, use_bedops = True, intersect = True)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_merge_bed(self):
        in_bed = "tests/merge_bed_input.bed"
        out_bed = "tests/merge_bed_output.bed"
        expected_bed = "tests/merge_bed_expected.bed"
        distance = 2
        merge_bed(in_bed, out_bed, distance)
        expected = rw.read_as_string(expected_bed)
        observed = rw.read_as_string(out_bed)
        self.assertEqual(expected, observed)

    def test_parse_3ss(self):
        input_file = "tests/parse_3ss_input.bed"
        expected = {}
        expected["FBgn1.0"] = [8, 4]
        expected["FBgn1.1"] = [11, 1]
        expected["FBgn1.2"] = [16, 2]
        expected["FBgn1.3"] = [20, 2]
        expected["FBgn2.0"] = [18, 2]
        expected["FBgn2.1"] = [13, 3]
        expected["FBgn2.2"] = [6, 5]
        observed = parse_3ss(input_file)
        self.assertEqual(expected, observed)

    def test_parse_exon_junctions(self):
        input_file = "tests/parse_exon_junctions_input.bed"
        expected = {}
        expected["FBgn1.0.5_2L_3"] = 8
        expected["FBgn1.1.5_2L_9"] = 11
        expected["FBgn1.2.5_2L_13"] = 16
        expected["FBgn1.3.5_2L_17"] = 20
        expected["FBgn2.2.5_2L_21"] = 18
        expected["FBgn2.1.5_2L_16"] = 14
        expected["FBgn2.0.5_2L_12"] = 6
        observed = parse_exon_junctions(input_file)
        self.assertEqual(expected, observed)

    def test_parse_exon_junctions2(self):
        input_file = "tests/parse_exon_junctions_input2.bed"
        expected = {}
        expected["FBgn1.0.5_2L_3"] = 8
        expected["FBgn1.1.5_2L_9"] = 11
        expected["FBgn1.2.5_2L_13"] = 16
        expected["FBgn1.3.5_2L_17"] = 20
        expected["FBgn2.2.5_2L_21"] = 18
        expected["FBgn2.1.5_2L_16"] = 14
        expected["FBgn2.0.5_2L_12"] = 6
        observed = parse_exon_junctions(input_file)
        self.assertEqual(expected, observed)

    def test_peak_pos_in_exon(self):
        exon_starts_file = "tests/peak_pos_in_exon_input_exon_starts.bed"
        peaks_file = "tests/peak_pos_in_exon_input_peaks.bed"
        expected_all = {"FBtr0077982.0": [0, 1, 2, 4],
                        "FBtr0078168.0": [0, 1],
                        "FBtr0078168.2": [3, 4],
                        "FBtr0078168.3": [],
                        "FBtr0077982.1": [],
                        "FBtr0077982.2": [],
                        "FBtr0078168.1": []}
        expected_centres = {"FBtr0077982.0": [0, 1, 4],
                    "FBtr0078168.0": [0],
                    "FBtr0078168.2": [3],
                        "FBtr0078168.3": [],
                        "FBtr0077982.1": [],
                        "FBtr0077982.2": [],
                        "FBtr0078168.1": []}
        expected = expected_all, expected_centres
        observed = peak_pos_in_exon(exon_starts_file, peaks_file)
        self.assertEqual(expected, observed)

    def test_peak_pos_in_exon_from_end(self):
        exon_starts_file = "tests/peak_pos_in_exon_from_end_input_exon_ends.bed"
        peaks_file = "tests/peak_pos_in_exon_from_end_input_peaks.bed"
        expected_all = {"FBtr0077982.0": [3, 5, 6],
                    "FBtr0078168.0": [1, 2],
                    "FBtr0077982.2": [1, 2, 3, 4],
                        "FBtr0077982.1": [],
                        "FBtr0078168.1": [],
                        "FBtr0078168.2": [],
                        "FBtr0078168.3": []}
        expected_centres = {"FBtr0077982.0": [5, 3],
                    "FBtr0078168.0": [1],
                    "FBtr0077982.2": [2],
                        "FBtr0077982.1": [],
                        "FBtr0078168.1": [],
                        "FBtr0078168.2": [],
                        "FBtr0078168.3": []}
        expected = expected_all, expected_centres
        observed = peak_pos_in_exon(exon_starts_file, peaks_file, from_end=True)
        self.assertEqual(expected, observed)

    def test_snr_bed(self):
        inbed = "tests/snr_bed_input.bed"
        outbed = "tests/snr_bed_observed.bed"
        snr_bed(inbed, outbed)
        expected = rw.read_as_string("tests/snr_bed_expected.bed")
        observed = rw.read_as_string(outbed)
        self.assertEqual(expected, observed)

    def test_snr_bed_fiveprime(self):
        inbed = "tests/snr_bed_fiveprime_input.bed"
        outbed = "tests/snr_bed_fiveprime_observed.bed"
        snr_bed(inbed, outbed, five_prime_most=True)
        expected = rw.read_as_string("tests/snr_bed_fiveprime_expected.bed")
        observed = rw.read_as_string(outbed)
        self.assertEqual(expected, observed)

    def test_sort_bed(self):
        infile = "tests/sort_bed_input.bed"
        expected_file = "tests/sort_bed_expected.bed"
        observed_file = "tests/sort_bed_observed.bed"
        hk.remove_file(observed_file)
        sort_bed(infile, observed_file)
        expected = rw.read_many_fields(expected_file, "\t")
        observed = rw.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_trim_sequence_phase0(self):
        sequences = ["GAGTATTCTGAG","GAGTATTCTGAGC","GAGTATTCTGAGCC"]
        expected = ["GAGTATTCTGAG" for i in range(3)]
        observed = [trim_sequence(sequences[i], 0) for i in range(3)]
        self.assertEqual(observed,expected)

    def test_trim_sequence_phase2(self):
        sequences = ["TTGAGTATTCTGAG","TTGAGTATTCTGAGG","TTGAGTATTCTGAGGC"]
        expected = ["GAGTATTCTGAG" for i in range(3)]
        observed = [trim_sequence(sequences[i], 2) for i in range(3)]
        self.assertEqual(observed,expected)

    def test_trim_sequence_phase1(self):
        sequences = ["TGAGTATTCTGAG","TGAGTATTCTGAGG","TGAGTATTCTGAGGC"]
        expected = ["GAGTATTCTGAG" for i in range(3)]
        observed = [trim_sequence(sequences[i], 1) for i in range(3)]
        self.assertEqual(observed,expected)

    def test_write_intron_start_peak_from_exons(self):
        exons = {"FBtr0077958": [['2L', 'FlyBase', 'CDS', 1143371, 1143410, '.', '+', '0', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1145444, 1146472, '.', '+', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1146539, 1146825, '.', '+', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";']],
                 "FBtr0077959": [['2L', 'FlyBase', 'CDS', 1146539, 1146825, '.', '-', '0', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1145444, 1146472, '.', '-', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1143371, 1143410, '.', '-', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";']]}
        outbed = "tests/write_intron_start_peak_from_exons_observed.bed"
        expected = "tests/write_intron_start_peak_from_exons_expected.bed"
        write_intron_start_peak_from_exons(exons, outbed, add_chr=False, alt_start=3, alt_end=6)
        observed = rw.read_as_string(outbed)
        expected = rw.read_as_string(expected)
        self.assertEqual(expected, observed)

    def test_write_si_pos_from_exons(self):
        exons = {"FBtr0077958": [['2L', 'FlyBase', 'CDS', 1143371, 1143410, '.', '+', '0', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1145444, 1146472, '.', '+', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1146539, 1146825, '.', '+', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077958"; transcript_symbol "MFS3-RA";']],
                 "FBtr0077959": [['2L', 'FlyBase', 'CDS', 1143371, 1143410, '.', '-', '0', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1145444, 1146472, '.', '-', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1146539, 1146825, '.', '-', '2', 'gene_id "FBgn0031307"; gene_symbol "MFS3"; transcript_id "FBtr0077959"; transcript_symbol "MFS3-RA";']]}
        outbed = "tests/write_si_pos_from_exons_observed.bed"
        expected = "tests/write_si_pos_from_exons_expected.bed"
        write_si_pos_from_exons(exons, outbed, add_chr=False)
        observed = rw.read_as_string(outbed)
        expected = rw.read_as_string(expected)
        self.assertEqual(expected, observed)