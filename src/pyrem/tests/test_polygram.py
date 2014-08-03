__author__ = 'quentin'

import unittest
import os
import tempfile

from pyrem.polygram import *
from pyrem.time_series import *


class TestPolygram(unittest.TestCase):

    np.random.seed(1)
    rw1 = np.cumsum(np.random.normal(0,1,2000))
    probs = np.random.random_sample(1000)
    vals = (np.random.random_sample(1000) * 4 +1).astype(np.int)

    def test_init_sanity(self):
        # an = Annotation(self.vals,fs=10, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        # c1 = Signal(self.rw1, fs=10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})
        #Polygram([an, c1])


        an = Annotation(self.vals,fs=1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=100,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        # perfect match:
        Polygram([an[:10], c1[:1000]])

        # tolerable match
        Polygram([an[:10], c1[:1099]])
        # intolerable match raises a value error
        self.assertRaises(ValueError, lambda: Polygram([an[:10], c1[:1100]]))

        # channel with the minimum sampling rate is the reference all channels should have equal of longer durationq
        self.assertRaises(ValueError, lambda: Polygram([an[:10], c1[:999]]))

        # channel with the same names (homonnymous) are forbidden !
        c2 = Signal(self.rw1, fs=100,type="eeg", name="bar", metadata={"animal":"joe", "treatment":18})
        self.assertRaises(ValueError, lambda: Polygram([an[:10], c2[:1000]]))



    def test_merging(self):
        an = Annotation(self.vals,fs=1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=100,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})
        c2 = Signal(self.rw1, fs=10,type="eeg", name="blah", metadata={"animal":"joe", "treatment":18})
        c3 = Signal(self.rw1, fs=50,type="eeg", name="argh", metadata={"animal":"joe", "treatment":18})

        pol = Polygram([an[:19], c1[:1999]])

        pol2 = Polygram([c2[:199], c3[:999]])

        pol.merge(pol2)
        pol.merge(c3)
        self.assertRaises(ValueError, lambda : pol.merge(c1))






    def test_slicing(self):
        an = Annotation(self.vals,fs=1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=100,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        pol = Polygram([an[:19], c1[:1999]])

        pol2 = pol["1s":"2s"]
        # sig = pol2[1]
        #sig += 1
        res = pol2["foo"]
        ans = c1["1s":"2s"]

        self.assertTrue(np.allclose(res, ans))

        self.assertTrue(pol2["foo"] is pol2["foo"])
        self.assertRaises(ValueError, lambda :pol2["DUMMY_NAME"])


    def test_copy(self):
        an = Annotation(self.vals,fs=50, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=100,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        pol = Polygram([an, c1])
        pol2 = pol.copy()
        sig = pol2["foo"]
        print sig
        sig += 1


        self.assertFalse(np.allclose(pol["foo"], sig ))
        self.assertFalse(np.allclose(c1, pol2["foo"]))
        self.assertTrue(c1 is pol["foo"])

    def test_windowing(self):
        an = Annotation(self.vals,fs=.1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        pol = Polygram([an[:10], c1[:1001]])

        for c, p in  pol.iter_window(1,1):
            self.assertEqual(dict([(c.name,c.size) for c in p.channels]),  {"bar":1, "foo": 10})

    def test_save_load(self):

        file, path = tempfile.mkstemp(suffix=".pkl")

        an = Annotation(self.vals,fs=.1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        pol = Polygram([an[:10], c1[:1001]])


        pol.save(path, 9)
        import joblib as pkl
        b = pkl.load(path)
        try:
            os.remove(path)
        except Exception as e:
            print e
            pass
        self.assertEqual(repr(pol[0]), repr(b[0]))
        self.assertEqual(repr(pol[1]), repr(b[1]))
        #self.assertTrue(compare_annots(a,b))
    def test_dummy(self):
        an = Annotation(self.vals,fs=.1, observation_probabilities=self.probs, type="vigilance_state", name="bar", metadata={"animal":"joe", "treatment":18})
        c1 = Signal(self.rw1, fs=10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

        print Polygram([an[:10], c1[:1001]])
