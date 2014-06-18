__author__ = 'quentin'

__author__ = 'quentin'

import numpy as np
import scipy.stats as stats
import pywt

from feature_base import FeatureGroup


class WaveletsFeaturesBase(FeatureGroup):
    prefix = "wavelets"
    _max_level = 5
    _wavelet =  None
    def __init__(self):
        if self._wavelet is None:
            raise NotImplementedError
        self.prefix = self.prefix + "." + self._wavelet
    def _make_feature_vec(self,signal):
        coeffs = pywt.wavedec(signal, self._wavelet, level=self._max_level)

        coeff_names = ["_cA_%i" % (len(coeffs) -1 ) ]
        for n in range(len(coeffs) -1, 0, -1):
             coeff_names.append("_cD_%i" % (n))

        out = dict()

        for c, n in zip(coeffs, coeff_names):
            c =  np.abs(c)
            out["mean" + n] = np.mean(c)
            out["sd"+ n] = np.std(c)
            out["median"+ n] = np.median(c)
            out["skew"+ n] = stats.skew(c)
            out["kurtosis"+ n] = stats.kurtosis(c)

        return out

class WaveletsFeaturesDB1(WaveletsFeaturesBase): _wavelet =  "db1"
class WaveletsFeaturesDB2(WaveletsFeaturesBase): _wavelet =  "db2"
class WaveletsFeaturesDB3(WaveletsFeaturesBase): _wavelet =  "db3"
class WaveletsFeaturesDB4(WaveletsFeaturesBase): _wavelet =  "db4"



