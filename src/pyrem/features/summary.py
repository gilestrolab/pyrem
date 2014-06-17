__author__ = 'quentin'

import numpy as np
import scipy.stats as stats

from feature_base import FeatureGroup


class PowerFeatures(FeatureGroup):
    prefix = "power"
    def _make_feature_vec(self,signal):
        out  = dict()
        powers = signal ** 2
        out["mean"] = np.mean(powers)
        out["median"] = np.median(powers)
        out["skew"] = stats.skew(powers)
        out["kurtosis"] = stats.kurtosis(powers)
        return out



