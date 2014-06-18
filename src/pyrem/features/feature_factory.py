__author__ = 'quentin'

import pandas as pd
from summary import PowerFeatures
from wavelets import *

class FeatureFactory(object):
    def __init__(self):
        self.feature_groups = [
            PowerFeatures(),
            # WaveletsFeaturesDB1(),
            # WaveletsFeaturesDB2(),
            # WaveletsFeaturesDB3(),
            WaveletsFeaturesDB4()
        ]
    def __call__(self, t, signal):
        dfs = [ group(t,signal) for group in self.feature_groups]
        return pd.concat(dfs, axis=1)