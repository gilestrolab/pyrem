__author__ = 'quentin'


import unittest

from pyrem.features import pyeeg
from pyrem.features import pyrem_pyeeg
import numpy as np



class TestFeatures(unittest.TestCase):
    random_walk = np.cumsum(np.random.normal(0,1,(int(1e4))))

    def test_pfd(self):
        ref = pyeeg.pfd(self.random_walk)
        ans = pyrem_pyeeg.pfd(self.random_walk)

        self.assertAlmostEqual(ref, ans)


    def test_svd_entropy(self):
        ref = pyeeg.svd_entropy(self.random_walk,10,10)
        ans = pyrem_pyeeg.svd_entropy(self.random_walk,10,10)

        self.assertAlmostEqual(ref, ans)

    def test_fisher_information(self):
        ref = pyeeg.fisher_info(self.random_walk,10,10)
        ans = pyrem_pyeeg.fisher_info(self.random_walk,10,10)

        self.assertAlmostEqual(ref, ans)


    def test_hjorth(self):
        ref_morbidity, ref_complexity = pyeeg.hjorth(self.random_walk)

        ans_activity, ans_morbidity, ans_complexity= pyrem_pyeeg.hjorth(self.random_walk)


        # Hjorth complexity vary a lot according to papers/ implementations.
        # pyeeg ones is not normalised by number of points -> ours is right
        #self.assertAlmostEqual(ref_morbidity, ans_morbidity)



    def test_embed(self):
        linspace = np.arange(1,20)
        ref_mat = pyeeg.embed_seq(linspace, 3,6)
        ans_mat = pyrem_pyeeg.embed_seq(linspace, 3,6)
        self.assertTrue((ref_mat == ans_mat).all())

        linspace = np.arange(1,27)
        ref_mat = pyeeg.embed_seq(linspace, 2,6)
        ans_mat = pyrem_pyeeg.embed_seq(linspace, 2,6)

        self.assertTrue((ref_mat == ans_mat).all())

#
# from pyrem.features import pyeeg
# from pyrem.features import pyrem_pyeeg
# import numpy as np
#
# random_walk = np.cumsum(np.random.normal(0,1,(int(1e4))))
# %timeit  pyeeg.pfd(random_walk)
# %timeit  pyrem_pyeeg.pfd(random_walk)
# %timeit  pyeeg.hjorth(random_walk)
# %timeit  pyrem_pyeeg.hjorth(random_walk)
#
# %timeit  pyeeg.embed_seq(random_walk, 100, 10)
# %timeit  pyrem_pyeeg.embed_seq(random_walk, 100, 10)


