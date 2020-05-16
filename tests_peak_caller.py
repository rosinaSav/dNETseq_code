from peak_caller import *
import read_and_write as rw
import unittest

class Test_peak_caller(unittest.TestCase):

    def test_filter_peaks(self):
        in_peak_bed = "tests/filter_peaks_input_peaks.bed"
        read_bed = "tests/filter_peaks_input_reads.bed"
        min_reads_per_peak = 5
        out_peak_bed = "tests/filter_peaks_observed.bed"
        min_peak_length = 4
        stats_file = "tests/filter_peaks_observed_stats.txt"
        counts_file = "tests/filter_peaks_input_counts.txt"
        filter_peaks(in_peak_bed, read_bed, counts_file, out_peak_bed, min_reads_per_peak, min_peak_length, stats_file)
        expected = rw.read_as_string("tests/filter_peaks_expected.bed")
        observed = rw.read_as_string(out_peak_bed)
        self.assertEqual(expected, observed)
        expected = rw.read_as_string("tests/filter_peaks_expected_stats.txt")
        observed = rw.read_as_string(stats_file)
        self.assertEqual(expected, observed)

    def test_get_reads_per_pos(self):
        reads_file = "tests/get_reads_per_pos_input_reads.bed"
        transcripts_file = "tests/get_reads_per_pos_input_transcripts.bed"
        expected = {}
        expected["chr2L.+.FBtr0307554"] = {"reads": {2: 2, 3: 4, 8: 1, 13: 2,
                                             15: 1, 24: 2, 25: 2,
                                             26: 2, 27: 3, 28: 1, 29: 1,
                                             30: 2, 31: 2, 32: 3, 33: 2},
                                           "coords": (2, 36)}
        expected["chr2L.-.FBtr0307600"] = {"reads": {34: 1, 28: 1, 11: 1, 13: 2, 14: 2, 15: 1, 18: 1, 22: 1, 25: 1},
                                           "coords": (5, 36)}
        observed = get_reads_per_pos(reads_file, transcripts_file)
        self.assertEqual(expected, observed)

    def test_make_read_count_array(self):
        exons = {"reads": {0: 2, 1: 4, 6: 1, 13: 1, 22: 2, 23: 2, 24: 2,
                           25: 3, 26: 1, 27: 1, 28: 2, 29: 2, 30: 3, 31: 2},
                 "coords": (0, 34)}
        expected = np.array([2, 4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 3, 1,
                             1, 2, 2, 3, 2, 0, 0])
        observed = make_read_count_array(exons)
        self.assertTrue((expected == observed).all())

    def test_make_read_count_array_minus(self):
        exons = {"reads": {32: 1, 26: 1, 9: 1, 11: 1, 12: 2, 13: 1, 16: 1},
                 "coords": (3, 34)}
        expected = np.array([0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0])
        observed = make_read_count_array(exons)
        self.assertTrue((expected == observed).all())

    def test_make_reads_dict_from_array(self):
        read_counts = np.array([0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0])
        coords = (3, 34)
        expected = {32: 1, 26: 1, 9: 1, 11: 1, 12: 2, 13: 1, 16: 1}
        observed = make_reads_dict_from_array(read_counts, coords)
        self.assertEqual(expected, observed)

    def test_sig_pos(self):
        coords = (10, 36)
        windowed_counts = np.array([0.5, 1/3, 1/3, 0, 1/3, 1/3, 1/3, 1.02,
                             0, 1/3, 1/3, 0.5,
                             1.5, 4/3, 1, 1/3, 1/3, 1/3, 4, 0, 0, 0, 0, 0, 0, 0])
        threshold = 1.01
        expected = [17, 22, 23, 28]
        observed = sig_pos(coords, windowed_counts, threshold)
        self.assertEqual(expected, observed)

    def test_window_read_count_array_plus(self):
        read_counts = np.array([0, 1, 0, 0, 0, 1, 0, 0,
                                0, 0, 1, 0,
                                1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        window_size = 3
        expected = np.array([1/3, 1/3, 1/3, 0, 1/3, 1/3, 1/3, 0,
                             0, 1/3, 1/3, 2/3,
                             1, 4/3, 1, 1/3, 1/3, 1/3, 1/3, 0, 0, 0, 0, 0, 0, 0])
        observed = window_read_count_array(read_counts, window_size, False)
        expected = list(expected)
        observed = list(observed)
        [self.assertAlmostEqual(expected[i], observed[i]) for i in range(len(expected))]

    def test_window_read_count_array_minus(self):
        read_counts = np.array([0, 1, 0,
                                0, 0, 1, 0, 0, 0,
                                0, 1, 0, 1, 2, 1, 0, 0, 1, 0, 0, 0,
                                0, 0, 0])
        window_size = 3
        expected = np.array([1/3, 1/3, 1/3,
                             0, 1/3, 1/3, 1/3, 0, 0,
                             1/3, 1/3, 2/3, 1, 4/3, 1,
                             1/3, 1/3, 1/3, 1/3, 0, 0,
                             0, 0, 0])
        observed = window_read_count_array(read_counts, window_size, False)
        expected = list(expected)
        observed = list(observed)
        [self.assertAlmostEqual(expected[i], observed[i]) for i in range(len(expected))]

    def test_window_read_count_array_no_slide(self):
        read_counts = np.array([0, 1, 0, 0, 0, 1, 0, 0,
                                0, 0, 1, 0,
                                1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        window_size = 3
        expected = np.array([1/3, 1/3, 0, 0, 0, 1/3, 1/3, 1/3,
                             1/3, 1/3, 1/3, 1,
                             1, 1, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 0, 0, 0, 0, 0, 0])
        observed = window_read_count_array(read_counts, window_size, True)
        expected = list(expected)
        observed = list(observed)
        [self.assertAlmostEqual(expected[i], observed[i]) for i in range(len(expected))]