__author__ = 'quentin'

import pandas as pd
from summary import PowerFeatures
from periodogram import *
from entropy import *
from non_linear import *
from hjorth import *
import multiprocessing as mp

__FEATURE_GROUPS = [
        PowerFeatures(),
        EntropyFeatures(),
        NonLinearFeatures(),
        Hjorth(),
        # WaveletsFeaturesDB1(),
        # WaveletsFeaturesDB2(),
        # WaveletsFeaturesDB3(),
        PeriodFeatures()
    ]

def _make_features(t_signal):
    t, signal = t_signal
    dfs = [ group(t,signal) for group in __FEATURE_GROUPS]
    return pd.concat(dfs, axis=1)


class FeatureFactory(object):


    def __call__(self, signal, t=np.NaN):
        return _make_features((t,signal))

    def make_features_for_epochs(self,signal, length, lag, processes=1):
        time_signals = [(t,s) for t, s in signal.embed_seq(length, lag)]
        if processes == 1:
            dfs = map(_make_features, time_signals)
        else:
            pool = mp.Pool(processes)
            dfs = pool.map(_make_features, time_signals)

        features = pd.concat(dfs)

        return features