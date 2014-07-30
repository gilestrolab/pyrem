__author__ = 'quentin'

import unittest
from pyrem.signal.polygram import *
from pyrem.signal.signal import *
import numpy as np
import os
import tempfile


class TestPolygram(unittest.TestCase):

    np.random.seed(1)
    rw1 = np.cumsum(np.random.normal(0,1,2000))
    probs = np.random.random_sample(1000)
    vals = (np.random.random_sample(1000) * 4 +1).astype(np.int)

    def test_init_sanity(self):
        an = Annotation(self.vals,fs=10, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})
        #Polygram([an, c1])


        an = Annotation(self.vals,fs=1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=100,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})


        Polygram([an[:10], c1[:1000]])
        Polygram([an[:10], c1[:901]])

