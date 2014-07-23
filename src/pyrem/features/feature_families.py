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
        dfs = [ group.make_vector(signal) for group in self._feature_group]
        return pd.concat(dfs, axis=1)


    def make_vector(self, signal, t=np.NaN):
        return self._make_features(signal)

    def major_annotations(self, annotations):
        a = annotations.flatten()
        r = np.real(a)
        uniqs = np.unique(r)
        i = np.imag(a)
        probs = []
        for u in uniqs:
            eqs = (r == u)
            probs.append(np.sum(i[eqs]))

        probs = np.array(probs)
        probs /= np.sum(i)


        max_prob_idx =  np.argmax(probs)

        return uniqs[max_prob_idx] + probs[max_prob_idx]* 1j


    def make_features_for_epochs(self, data, length, lag, add_major_annotations=False, processes=1):
        r"""
        Compute features, for all channels and all epochs, of a polygraph.

        :param data: A polygraph
        :type data: :class:`~pyrem.signal.polygraph.Polygraph`

        :param length: the length of the epoch in seconds (see :meth:`~pyrem.signal.polygraph.Polygraph.embed_seq`)
        :type length: float
        :param lag: the lag of the epoch in seconds (see :meth:`~pyrem.signal.polygraph.Polygraph.embed_seq`)
        :type lag: float
        :return: a dataframe
        :rtype: :class:`~pandas.DataFrame`
        """
        rows = []
        for t, s in data.embed_seq(length, lag):
            if add_major_annotations:
                majors = [self.major_annotations(a) for a in s.annotations()]




            for c in s.channels():
                row = self._make_features(c)

                row["channel"] = [c.channel_types[0]]
                if add_major_annotations:
                    for ann, maj  in zip(s.annotation_types, majors):
                        row[ann+"_value"] = np.real(maj)
                        row[ann+"_prob"] = np.imag(maj)


                row.index = [t]
                rows.append(row)

        features = pd.concat(rows)
        return features

class FeatureFamilyBase(object):
    r"""
    A feature family object is a process returning a vector of features upon analysis of some data.
    Features are returned as a pandas DataFrame object, with column names for features. Each feature name is prefixed by the name of the
    feature family. This is an abstract class designed to be derived by:

    1. Defining a ``prefix`` attribute. It will add the name of the family to the name of the features.
    2. Overriding the ``_make_feature_vec`` method. It should return a dictionary of scalars, each being a feature.

    """
    prefix = None
    def make_vector(self, data):
        """
        Compute one vector of features from polygraph.

        :param data: A polygraph
        :type data: :class:`~pyrem.signal.polygraph.Polygraph`
        :return: a one-row dataframe
        :rtype: :class:`~pandas.DataFrame`
        """
        if data.n_channels != 1:
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
        out = dict()

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
        #out["dfa"] = dfa(data)
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
        #out["spectral"] = spectral_entropy(data,np.arange(0,50), channel.sampling_freq) # fixme magic number here
        out["svd"] = svd_entropy(data, 3,3) # fixme magic number here
        out["fisher"] = fisher_info(data, 3,3)
        #out["apent"] = ap_entropy(data, 2,5000)
        out["sample_2_1000"] = samp_entropy(data, 2, 1000)
        out["sample_2_2000"] = samp_entropy(data, 2, 2000)
        out["sample_2_5000"] = samp_entropy(data, 2, 5000)
        out["sample_3_1000"] = samp_entropy(data, 3, 1000)
        out["sample_3_2000"] = samp_entropy(data, 3, 2000)
        out["sample_3_5000"] = samp_entropy(data, 3, 5000)
        out["sample_4_1000"] = samp_entropy(data, 4, 1000)
        out["sample_4_2000"] = samp_entropy(data, 4, 2000)
        out["sample_4_5000"] = samp_entropy(data, 4, 5000)

        return out

class MSEFeatures(FeatureFamilyBase):
    prefix = "mse"

    def _make_feature_vec(self,channel):
        data = channel.data.flatten()
        out = dict()
        #out["spectral"] = spectral_entropy(data,np.arange(0,50), channel.sampling_freq) # fixme magic number here
        r = 0.15
        for i in range(1,4):

            name = "%i_%i" %(i, int(r*100))
            out[name] = samp_entropy(data,2,r,i)
            # print i, int(r*100), out[name]

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
        # out["entropy"] = stats.entropy(pow_f)  * freqs.size #TODO

        return out