__author__ = 'quentin'


import numpy as np
from pyrem.features.feature_base import FeatureGroup


class NonLinearFeatures(FeatureGroup):
    prefix = "nl"

    def hurst_exponent(self, signal):
        """
        Code adapted from Forrest Sheng Bao's pyeeg
        """
        X = signal
        N = len(X)

        T = np.array([float(i) for i in xrange(1,N+1)])
        Y = np.cumsum(X)
        Ave_T = Y/T

        S_T = np.zeros((N))
        R_T = np.zeros((N))
        for i in xrange(N):
            S_T[i] = np.std(X[:i+1])
            X_T = Y - T * Ave_T[i]
            R_T[i] = max(X_T[:i + 1]) - min(X_T[:i + 1])

        R_S = R_T / S_T
        R_S = np.log(R_S)
        n = np.log(T).reshape(N, 1)
        H = np.linalg.lstsq(n[1:], R_S[1:])[0]
        return H[0]


    def hurst(self, signal):
        """
        from:
        http://drtomstarke.com/index.php/calculation-of-the-hurst-exponent-to-test-for-trend-and-mean-reversion/
        """
        tau = []; lagvec = []

        #  Step through the different lags
        for lag in range(2,20):

        #  produce price difference with lag
            pp = np.subtract(signal[lag:],signal[:-lag])

        #  Write the different lags into a vector
            lagvec.append(lag)

        #  Calculate the variance of the difference vector
            tau.append(np.std(pp))

        #  linear fit to double-log graph (gives power)
        m = np.polyfit(np.log10(lagvec),np.log10(tau),1)

        # calculate hurst
        hurst = m[0]

        return hurst



    def _make_feature_vec(self,signal):

        out = dict()

        #out["hurstExp"] = self.hurst_exponent(signal)
        out["hurst"] = self.hurst(signal)


        return out
