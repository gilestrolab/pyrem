__author__ = 'quentin'

import pandas as pd




class FeatureGroup(object):
    prefix = None
    def __call__(self, t, signal):
        feature_dict = self._make_feature_vec(signal)

        data_frame = pd.DataFrame(feature_dict, index=[t])


        if len(feature_dict) > 1 and self.prefix is None:
            raise Exception("More than one features in this group. You need a prefix to identify this group")

        if self.prefix:
            data_frame .columns = [self.prefix + "_" + c for c in data_frame .columns]
        return data_frame

    def _make_feature_vec(self,signal):
        raise NotImplementedError