from read_and_write import *
import unittest

class Test_read_and_write(unittest.TestCase):

    def test_read_fasta(self):
        expected = {"1:17-25(+)": "CATAGACA", "2:0-12(+)": "GTCCCCCCCCAA", "1:21-30(-)": "AAATATGTC"}
        observed = read_fasta("tests/read_fasta_input.fasta", return_dict = True)
        self.assertEqual(expected, observed)

    def test_read_gtf_exons_gene(self):
        expected = {}
        expected["FBgn1"] = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]

        expected["FBgn2"] = [["2L", "FlyBase", "exon", 14, 350, ".", "+", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]
        observed = read_gtf("tests/read_gtf_input_file.gtf", "exon", gene = True)
        self.assertEqual(expected, observed)

    def test_read_gtf_start_codons_gene(self):
        expected = {}
        expected["FBgn3"] = [["2L", "FlyBase", "start_codon", 729, 731, ".", "+", "0",
                              "gene_id \"FBgn3\"; gene_symbol \"vrl\"; transcript_id \"FBtr55\"; transcript_symbol \"vrl-RA\";"]]
        observed = read_gtf("tests/read_gtf_input_file.gtf", "start_codon", gene = True)
        self.assertEqual(expected, observed)

    def test_read_gtf_exons_trans(self):
        expected = {}
        expected["FBtr33"] = [["2L", "FlyBase", "exon", 13, 20, ".", "+", ".",
                              "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 30, 38, ".", "+", ".",
                             "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"],
                            ["2L", "FlyBase", "exon", 60, 70, ".", "+", ".",
                            "gene_id \"FBgn1\"; gene_symbol \"drl\"; transcript_id \"FBtr33\"; transcript_symbol \"drl-RA\";"]]
        expected["FBtr44"] = [["2L", "FlyBase", "exon", 14, 350, ".", "+", ".",
                              "gene_id \"FBgn2\"; gene_symbol \"trl\"; transcript_id \"FBtr44\"; transcript_symbol \"trl-RA\";"]]
        observed = read_gtf("tests/read_gtf_input_file.gtf", "exon", gene = False)
        self.assertEqual(expected, observed)

    def test_read_gtf_start_codons_trans(self):
        expected = {}
        expected["FBtr55"] = [["2L", "FlyBase", "start_codon", 729, 731, ".", "+", "0",
                              "gene_id \"FBgn3\"; gene_symbol \"vrl\"; transcript_id \"FBtr55\"; transcript_symbol \"vrl-RA\";"]]
        observed = read_gtf("tests/read_gtf_input_file.gtf", "start_codon", gene = False)
        self.assertEqual(expected, observed)

    def test_read_many_fields(self):
        expected = [["hello", "5", "10"],["goodbye", "4456", "897"]]
        observed = read_many_fields("tests/read_many_fields_input_file.txt", "\t")
        self.assertEqual(expected, observed)