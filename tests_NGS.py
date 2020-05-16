from NGS import *
import read_and_write as rw
import unittest

class Test_NGS(unittest.TestCase):

    def test_analyze_cigar_S_plus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '2M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = "S"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar2_S_minus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = "S"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar3_U_plus(self):
        bed_line = ['chr2L', '6', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '4M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = "U"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar4_U_minus(self):
        bed_line = ['chr2L', '2', '7',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '5M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = "U"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar5_doesnt_reach_intron_plus(self):
        bed_line = ['chr2L', '8', '12',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '4M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar6_doesnt_reach_intron_minus(self):
        bed_line = ['chr2L', '2', '5',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar7_wrong_intron_length_plus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '2M4N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar8_wrong_intron_length_minus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M9N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar9_various_char_plus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '1D1M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = "S"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar10_S_various_char_minus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '1M1D1M3N1I2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = "S"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar11_U_various_char_plus(self):
        bed_line = ['chr2L', '6', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '2M1D1M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = "U"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar12_U_various_char_minus(self):
        bed_line = ['chr2L', '2', '7',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M2I3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = "U"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar13_wrong_intron_pos_plus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '4M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar14_wrong_intron_pos_minus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M3N1D1I133M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar15_minus_various_char(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M3N120I1D1M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = "S"
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar16_upstream_overhang_plus(self):
        bed_line = ['chr2L', '4', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '1M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar17_upstream_overhang_minus(self):
        bed_line = ['chr2L', '3', '9',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M3N1M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar18_downstream_overhang_plus(self):
        bed_line = ['chr2L', '3', '9',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '2M3N1M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar19_downstream_overhang_minus(self):
        bed_line = ['chr2L', '4', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '1M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar20_U_overhang_plus(self):
        bed_line = ['chr2L', '6', '9',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_analyze_cigar21_U_overhang_minus(self):
        bed_line = ['chr2L', '4', '7',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = None
        observed = analyze_cigar(bed_line, 2)
        self.assertEqual(expected, observed)

    def test_check_intermediate_read_plus_True(self):
        bed_line = ['chr2L', '2', '13',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '3M3N5M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = True
        observed = check_intermediate_read(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_intermediate_read_minus_True(self):
        bed_line = ['chr2L', '8', '19',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '5M3N3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '12', '13', 'FBtr0078168.4', '3', '-', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = True
        observed = check_intermediate_read(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_intermediate_read_plus_False(self):
        bed_line = ['chr2L', '2', '14',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '3M3N5M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = False
        observed = check_intermediate_read(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_intermediate_read_minus_False(self):
        bed_line = ['chr2L', '9', '19',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '5M3N3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '12', '13', 'FBtr0078168.4', '3', '-', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = False
        observed = check_intermediate_read(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_position_in_exon_plus_True(self):
        bed_line = ['chr2L', '2', '12',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '3M3N5M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.0', '3', '+', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = True
        observed = check_position_in_exon(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_position_in_exon_minus_True(self):
        bed_line = ['chr2L', '11', '19',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '5M3N3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '12', '13', 'FBtr0078168.0', '3', '-', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = True
        observed = check_position_in_exon(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_position_in_exon_plus_False(self):
        bed_line = ['chr2L', '2', '18',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '3M3N5M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.0', '3', '+', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 2, 5, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 17, 20, '.', '+', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = False
        observed = check_position_in_exon(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_check_position_in_exon_minus_False(self):
        bed_line = ['chr2L', '2', '19',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '5M3N3M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '12', '13', 'FBtr0078168.0', '3', '-', '1']
        exons = {'FBtr0077982': [['2L', 'FlyBase', 'CDS', 1163000, 1163088, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 1161927, 1162146, '.', '-', '1', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']],
                 'FBtr0078168': [['2L', 'FlyBase', 'CDS', 17, 20, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 9, 13, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";'],
                                 ['2L', 'FlyBase', 'CDS', 2, 5, '.', '-', '0', 'gene_id "FBgn0031313"; gene_symbol "CG5080"; transcript_id "FBtr0077982"; transcript_symbol "CG5080-RA";']]}
        expected = False
        observed = check_position_in_exon(bed_line, exons)
        self.assertEqual(expected, observed)

    def test_density_per_transcript(self):
        exon_file = "tests/density_per_transcript_input_exons.gtf"
        polII_bed = "tests/density_per_transcript_input_polII.bed"
        observed = "tests/density_per_transcript_observed.txt"
        expected = "tests/density_per_transcript_expected.txt"
        density_per_transcript(exon_file, polII_bed, observed)
        expected = rw.read_as_string(expected)
        observed = rw.read_as_string(observed)
        self.assertEqual(expected, observed)

    def test_density_per_transcript_bed(self):
        exon_file_bed = "tests/density_per_transcript_bed_input_exons.bed"
        polII_bed = "tests/density_per_transcript_input_polII.bed"
        observed = "tests/density_per_transcript_bed_observed.txt"
        expected = "tests/density_per_transcript_expected.txt"
        density_per_transcript(exon_file_bed, polII_bed, observed, bed_input=True)
        expected = rw.read_as_string(expected)
        observed = rw.read_as_string(observed)
        self.assertEqual(expected, observed)

    def test_get_splice_dist_plus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '+',
                    '16', '2M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '8', '9', 'FBtr0078168.4', '3', '+', '1']
        expected = 2
        observed = get_splice_dist(bed_line)
        self.assertEqual(expected, observed)

    def test_get_splice_dist_minus(self):
        bed_line = ['chr2L', '3', '10',
                    'ST-E00144:677:HKCCLCCXY:5:2121:17178:4456', '255', '-',
                    '16', '2M3N2M', '*', '0', '0',
                    'ACAGATAGGCGACCTACATTTTCATA', 'JJJJJJJJJJJJJJJJJJJJJJ',
                    'NH:i:1', 'HI:i:1', 'AS:i:87', 'nM:i:0', 'chr2L',
                    '4', '5', 'FBtr0078168.4', '3', '-', '1']
        expected = 2
        observed = get_splice_dist(bed_line)
        self.assertEqual(expected, observed)

    def test_get_spliced_reads_no_filters(self):
        reads_file = "tests/get_spliced_reads_input_reads.bed"
        output_file_spliced = "tests/get_spliced_reads_observed_spliced.bed"
        output_file_unspliced = "tests/get_spliced_reads_observed_unspliced.bed"
        exon_junctions_file = "tests/get_spliced_reads_input_junctions.bed"
        expected_file_spliced = "tests/get_spliced_reads_expected_spliced.bed"
        expected_file_unspliced = "tests/get_spliced_reads_expected_unspliced.bed"
        get_spliced_reads(reads_file, exon_junctions_file, output_file_spliced, output_file_unspliced, exons=False, overhang=0)
        expected_spliced = rw.read_as_string(expected_file_spliced)
        expected_unspliced = rw.read_as_string(expected_file_unspliced)
        expected = expected_spliced + expected_unspliced
        observed_spliced = rw.read_as_string(output_file_spliced)
        observed_unspliced = rw.read_as_string(output_file_unspliced)
        observed = observed_spliced + observed_unspliced
        self.assertEqual(expected, observed)

    def test_get_spliced_reads_exon_filters(self):
        reads_file = "tests/get_spliced_reads_input_reads.bed"
        output_file_spliced = "tests/get_spliced_reads_exon_filters_observed_spliced.bed"
        output_file_unspliced = "tests/get_spliced_reads_exon_filters_observed_unspliced.bed"
        exon_junctions_file = "tests/get_spliced_reads_input_junctions.bed"
        expected_file_spliced = "tests/get_spliced_reads_exon_filters_expected_spliced.bed"
        expected_file_unspliced = "tests/get_spliced_reads_expected_unspliced.bed"
        exons = {}
        exons["FBtr0078168"] = [["chr2L", "havana", "exon", 2, 6, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 11, 13, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 17, 19, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 22, 26, "FBtr0078168", ".", "+"]]
        exons["FBtr0078169"] = [["chr2L", "havana", "exon", 22, 26, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 17, 19, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 11, 13, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 2, 6, "FBtr0078169", ".", "-"]]
        get_spliced_reads(reads_file, exon_junctions_file, output_file_spliced, output_file_unspliced, exons=exons, overhang=0)
        expected_spliced = rw.read_as_string(expected_file_spliced)
        expected_unspliced = rw.read_as_string(expected_file_unspliced)
        expected = expected_spliced + expected_unspliced
        observed_spliced = rw.read_as_string(output_file_spliced)
        observed_unspliced = rw.read_as_string(output_file_unspliced)
        observed = observed_spliced + observed_unspliced
        self.assertEqual(expected, observed)

    def test_get_spliced_reads_exon_filters_exon_start_window(self):
        reads_file = "tests/get_spliced_reads_exon_start_window_input_reads.bed"
        output_file_spliced = "tests/get_spliced_reads_exon_filters_exon_start_window_observed_spliced.bed"
        output_file_unspliced = "tests/get_spliced_reads_exon_filters_exon_start_window_observed_unspliced.bed"
        exon_junctions_file = "tests/get_spliced_reads_input_junctions.bed"
        expected_file_spliced = "tests/get_spliced_reads_exon_filters_exon_start_window_expected_spliced.bed"
        expected_file_unspliced = "tests/get_spliced_reads_exon_start_window_expected_unspliced.bed"
        exons = {}
        exons["FBtr0078168"] = [["chr2L", "havana", "exon", 2, 6, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 11, 13, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 17, 19, "FBtr0078168", ".", "+"],
                                ["chr2L", "havana", "exon", 22, 26, "FBtr0078168", ".", "+"]]
        exons["FBtr0078169"] = [["chr2L", "havana", "exon", 22, 26, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 17, 19, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 11, 13, "FBtr0078169", ".", "-"],
                                ["chr2L", "havana", "exon", 2, 6, "FBtr0078169", ".", "-"]]
        get_spliced_reads(reads_file, exon_junctions_file, output_file_spliced, output_file_unspliced, exons=exons, overhang=0, filter_start=1, filter_end=3)
        expected_spliced = rw.read_as_string(expected_file_spliced)
        expected_unspliced = rw.read_as_string(expected_file_unspliced)
        expected = expected_spliced + expected_unspliced
        observed_spliced = rw.read_as_string(output_file_spliced)
        observed_unspliced = rw.read_as_string(output_file_unspliced)
        observed = observed_spliced + observed_unspliced
        self.assertEqual(expected, observed)





