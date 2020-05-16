from make_stats_for_splice_sites import *
import numpy as np
import read_and_write as rw
import unittest

class Test_make_stats_for_splice_sites(unittest.TestCase):

    def test_add_to_array(self):
        out_array = np.array([["trans", "number1", "number2"], ["trans1", 5, 6], ["transA", 12, 4], ["trans0", 0, 188]])
        new_dict = {"transA": 2, "trans0": 1, "trans1": 99}
        header = "new_number"
        expected = np.array([["trans", "number1", "number2", "new_number"], ["trans1", 5, 6, 99], ["transA", 12, 4, 2], ["trans0", 0, 188, 1]])
        observed = add_to_array(out_array, new_dict, header)
        self.assertTrue((expected == observed).all())