"""
The goal of this submodule is to provide a flexible interface to compute arbitrary features on each channel and epoch (temporal slices) of a multivariate time series (Polygraph).
Features are grouped in families of several features (e.g. Power Features may contain mean power, variance of power, ...).
Feature factory computes features for arbitrary feature families and group them in a data.frame
"""

__author__ = 'quentin'


import pandas as pd
from univariate import *
import scipy.stats as stats
import numpy as np
import pywt
from scipy import signal as sig


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

class FeatureFamilyBase(object):
    """
    A feature family object is a process returning a vector of features upon analysis of some data.
    Features are returned as a pandas DataFrame object, with column names for features. Each feature name is prefixed by the name of the
    feature family. This is an abstract class designed to be derived by:

    """
    prefix = None
    def __call__(self, data):
        """
        Compute one vector of features from polygraph.

        :param data : A polygraph
        :type data : :class:`~pyrem.signal.polygraph.Polygraph`
        :return: a one-row dataframe
        :rtype: :class:`~pandas.DataFrame`
        """
        if data.nchannels != 1:
            raise NotImplementedError("Only data with one channel can be analysed")
        feature_dict = self._make_feature_vec(data)

        data_frame = pd.DataFrame(feature_dict, index=[None])

        if len(feature_dict) > 1 and self.prefix is None:
            raise Exception("More than one features in this group. You need a prefix to identify this group")

        if self.prefix:
            data_frame.columns = [self.prefix + "_" + c for c in data_frame .columns]
        return data_frame

    def _make_feature_vec(self,data):
        raise NotImplementedError


class PowerFeatures(FeatureFamilyBase):
    prefix = "power"
    def _make_feature_vec(self, channel):
        data = channel.data.flatten()
        out  = dict()

        powers = data ** 2
        out["mean"] = np.mean(powers)
        out["sd"] = np.std(powers)
        out["median"] = np.median(powers)
        out["skew"] = stats.skew(powers)
        out["kurtosis"] = stats.kurtosis(powers)

        return out


class NonLinearFeatures(FeatureFamilyBase):
    prefix = "nl"
    def _make_feature_vec(self,channel):
        data = channel.data.flatten()
        out = dict()
        out["hurst"] = hurst(data)
        return out


class HjorthFeatures(FeatureFamilyBase):
    prefix = "hjorth"
    def _make_feature_vec(self, channel):
        data = channel.data.flatten()
        a,m,c = hjorth(data)
        out = {"activity":a, "morbidity":m, "complexity":c}
        return out


class EntropyFeatures(FeatureFamilyBase):
    prefix = "entropy"

    def _make_feature_vec(self,channel):
        data = channel.data.flatten()
        out = dict()
        out["spectral"] = spectral_entropy(data,np.arange(0,50), channel.sampling_freq) # fixme magic number here
        out["svd"] = svd_entropy(data, 3,3) # fixme magic number here
        return out

class WaveletsFeaturesBase(FeatureFamilyBase):
    prefix = "wavelets"
    _max_level = 5
    _wavelet =  None
    def __init__(self):
        if self._wavelet is None:
            raise NotImplementedError
        self.prefix = self.prefix + "." + self._wavelet

    def _make_feature_vec(self,channel):
        data = channel.data.flatten()
        coeffs = pywt.wavedec(data, self._wavelet, level=self._max_level)
        coeff_names = ["_cA_%i" % (len(coeffs) -1 ) ]
        for n in range(len(coeffs) -1, 0, -1):
             coeff_names.append("_cD_%i" % (n))
        out = dict()
        for c, n in zip(coeffs, coeff_names):
            c =  np.abs(c)
            out["mean" + n] = np.mean(c)
            out["sd"+ n] = np.std(c)
            # out["median"+ n] = np.median(c)
            # out["skew"+ n] = stats.skew(c)
            # out["kurtosis"+ n] = stats.kurtosis(c)
        return out

class WaveletsFeaturesDB1(WaveletsFeaturesBase): _wavelet =  "db1"
class WaveletsFeaturesDB2(WaveletsFeaturesBase): _wavelet =  "db2"
class WaveletsFeaturesDB3(WaveletsFeaturesBase): _wavelet =  "db3"
class WaveletsFeaturesDB4(WaveletsFeaturesBase): _wavelet =  "db4"


class PeriodFeatures(FeatureFamilyBase):
    prefix = "welch"

    def _make_feature_vec(self, channel):
        data = channel.data.flatten()

        freqs, pow = sig.welch(data, channel.sampling_freq)

        pow_f = freqs * pow/np.sum(pow)

        out = dict()

        out["mode"] = freqs[np.argmax(pow)]
        out["mean"] = np.mean(pow_f) * freqs.size
        out["sd"] = np.std(pow_f) * freqs.size
        out["mean"] = np.mean(pow_f) * freqs.size
        out["median"] = np.median(pow_f)  * freqs.size
        out["skew"] = stats.skew(pow_f)  * freqs.size
        out["kurtosis"] = stats.kurtosis(pow_f)  * freqs.size

        return out