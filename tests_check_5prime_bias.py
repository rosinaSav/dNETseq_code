from check_5prime_bias import *
import numpy as np
import unittest

class Test_nucleotide_comp(unittest.TestCase):

    def test_extract_true_and_control_string(self):
        fasta_name = "tests/extract_true_and_control_string_input.fasta"
        true_indices = (0, 5)
        control_indices = (8, 13)
        expected1 = np.array([["A", "T", "G", "C", "G"], ["A", "C", "C", "G", "T"], ["T", "G", "G", "G", "C"]])
        expected2 = np.array([["G", "C", "C", "A", "C"], ["T", "C", "A", "G", "A"], ["A", "C", "A", "C", "A"]])
        observed = extract_true_and_control_string(fasta_name, true_indices, control_indices)
        observed = (np.ndarray.tolist(observed[0]), np.ndarray.tolist(observed[1]))
        expected = (np.ndarray.tolist(expected1), np.ndarray.tolist(expected2))
        self.assertEqual(expected, observed)

