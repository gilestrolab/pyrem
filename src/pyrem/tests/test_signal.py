__author__ = 'quentin'

import unittest

from pyrem.signal import Signal
import numpy as np

class TestSignal(unittest.TestCase):

    def test_behaves_like_array(self):
        seq = [1,4,5,3,2,1,6,3,2]
        sgn = Signal(seq, 200)
        self.assertEqual(np.sum(seq), np.sum(sgn))


