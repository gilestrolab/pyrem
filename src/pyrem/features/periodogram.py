__author__ = 'quentin'

import numpy as np
import scipy.stats as stats
from scipy import signal as sig
from pyrem.features.feature_base import FeatureGroup


class PeriodFeatures(FeatureGroup):
    prefix = "welch"


    def _make_feature_vec(self,signal):

        freqs, pow = sig.welch(signal, signal.sampling_freq)

        pow_f = freqs * pow / np.sum(pow)

        out = dict()

        out["mode"] = freqs[np.argmax(pow)]
        out["mean"] = np.mean(pow_f) * freqs.size
        out["sd"] = np.std(pow_f) * freqs.size
        out["mean"] = np.mean(pow_f) * freqs.size
        out["median"] = np.median(pow_f)  * freqs.size
        out["skew"] = stats.skew(pow_f)  * freqs.size
        out["kurtosis"] = stats.kurtosis(pow_f)  * freqs.size

        return out
