__author__ = 'quentin'

import numpy as np
import scipy.stats as stats
from scipy import signal as sig
from pyrem.features.feature_base import FeatureGroup


class EntropyFeatures(FeatureGroup):
    prefix = "entropy"


    def spectral_entropy(self, signal, Band):
        """
        Code adapted from Forrest Sheng Bao's pyeeg
        """
        X = signal
        Fs = signal.sampling_freq
        C = np.fft.fft(X)
        C = abs(C)
        Power = np.zeros(len(Band)-1);
        for Freq_Index in xrange(0,len(Band)-1):
            Freq = float(Band[Freq_Index])                                        ## Xin Liu
            Next_Freq = float(Band[Freq_Index+1])
            Power[Freq_Index] = sum(C[np.floor(Freq / Fs*len(X)):np.floor(Next_Freq / Fs*len(X))])


        Power_Ratio = Power/float(np.sum(Power))

        Spectral_Entropy = 0
        for i in xrange(0, len(Power_Ratio) - 1):
            Spectral_Entropy += Power_Ratio[i] * np.log(Power_Ratio[i])
        Spectral_Entropy /= np.log(len(Power_Ratio))    # to save time, minus one is omitted
        return -1 * Spectral_Entropy



    def svd_entropy(self, signal, lag_s):
        """
        Code adapted from Forrest Sheng Bao's pyeeg
        """

        lag = lag_s * signal.sampling_freq
        mat = np.array([s for _,s in signal.embed_seq(lag*2,lag)])
        W = np.linalg.svd(mat, compute_uv = 0)
        W /= sum(W) # normalize singular values

        return -1*sum(W * np.log(W))

    #todo approximate entropy !

    def _make_feature_vec(self,signal):

        out = dict()
        out["spectral"] = self.spectral_entropy(signal,np.arange(0,50)) # fixme magic number here
        out["svd"] = self.svd_entropy(signal, 2) # fixme magic number here
        return out
