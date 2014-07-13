__author__ = 'quentin'

import pandas as pd
import numpy as np

class FeatureFactory(object):

    def __init__(self, feature_groups):
        self._feature_group = feature_groups

    def _make_features(self,signal):

        dfs = [ group(signal) for group in self._feature_group]
        return pd.concat(dfs, axis=1)


    def __call__(self, signal, t=np.NaN):
        return self._make_features(signal)

    def make_features_for_epochs(self,signal, length, lag, processes=1):
        rows = []
        for t, s in signal.embed_seq(length, lag):
            for c in s.channels():
                row = self._make_features(c)

                row["channel"] = [c.channel_types[0]]
                row.index = [t]
                rows.append(row)
        features = pd.concat(rows)
        return features