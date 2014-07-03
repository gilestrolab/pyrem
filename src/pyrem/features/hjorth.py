__author__ = 'quentin'

import numpy as np
from pyrem.features.feature_base import FeatureGroup


class Hjorth(FeatureGroup):
    prefix = "hjorth"



    def _make_feature_vec(self,signal):


        first_deriv = np.diff(signal)
        second_deriv = np.diff(signal,2)

        var_zero = np.var(signal)
        var_d1 = np.var(first_deriv)
        var_d2 = np.var(second_deriv)

        out = dict()
        out["activity"] = var_zero ** 2.0
        out["morbidity"] = var_d1 / var_zero
        out["complexity"] = np.abs(var_d2 / var_d1) - out["morbidity"] ** 2.0

        return out
