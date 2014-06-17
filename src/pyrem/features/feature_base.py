__author__ = 'quentin'

import pandas as pd

class FeatureGroup(pd.DataFrame):
    prefix = None
    def __init__(self, signal):
        feature_dict = self._make_feature_vec(signal)

        super(FeatureGroup,self).__init__(data=feature_dict, index=[0])

        if len(feature_dict) > 1 and self.prefix is None:
            raise Exception("More than one features in this group. You need a prefix to identify this group")

        if self.prefix:
            self.columns = [self.prefix + "_" + c for c in self.columns]

    def _make_feature_vec(self,signal):
        raise NotImplementedError