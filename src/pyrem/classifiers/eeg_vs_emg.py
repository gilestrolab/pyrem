__author__ = 'quentin'




import pandas as pd
from pyrem.features.summary import PowerFeatures
from pyrem.features.periodogram import PeriodFeatures


__FEATURE_GROUPS = [
        PowerFeatures(),
        PeriodFeatures()
    ]


def _make_features(channel):
    dfs = [ group(channel) for group in __FEATURE_GROUPS]
    return pd.concat(dfs, axis=1)



class EEGsvEMG(object):
    __FEATURE_GROUPS = [
            PowerFeatures(),
            PeriodFeatures()
        ]


    def __init__(self):
        self._features = pd.DataFrame()
    def make_features(self, data_iterator, data, processes=1):
        all_dfs = []
        for  channel  in data_iterator:
            groups_dfs = [group(channel, data) for group in self.__FEATURE_GROUPS]
            all_dfs.append(pd.concat(groups_dfs, axis=1))

        features = pd.concat(all_dfs)

        return  features

    def train(self, labels, data_iterator, data):
        new_features_df = self.make_features(data_iterator, data)
        new_features_df.index = labels
        self._features = pd.concat([self._features, new_features_df])
        #todo append = cleanner



    def predict(self):
        pass