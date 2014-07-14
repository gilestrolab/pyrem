__author__ = 'quentin'
import numpy as np

def _embed_seq(X,tau,de):

    N =len(X)

    if de * tau > N:
        raise ValueError("Cannot build such a matrix, because D * Tau > N")

    if tau<1:
        raise ValueError("Tau has to be at least 1")


    Y=np.zeros((de, N - (de - 1) * tau))

    for i in range(de):
        Y[i] = X[i *tau : i*tau + Y.shape[1] ]

    return Y.T

def _make_cmp(X, M, R, in_range_i, in_range_j, ap_ent=True):
     #Then we make Cmp
    N = len(X)

    Emp = _embed_seq(X, 1, M + 1)
    inrange_cmp = np.abs(Emp[in_range_i,-1] - Emp[in_range_j,-1]) <= R

    in_range_cmp_i = in_range_i[inrange_cmp]


    if not ap_ent:
        Cmp = np.bincount(in_range_cmp_i, minlength=N-M-1)
        return Cmp.astype(np.float)
    else:
        Cmp = np.bincount(in_range_cmp_i, minlength=N-M)
        in_range_cmp_j = in_range_j[inrange_cmp]
        Cmp += np.bincount(in_range_cmp_j, minlength=N-M)

        return Cmp.astype(np.float)

def _make_cm(X,M,R,ap_ent=True):
    N = len(X)

    # we pregenerate all indices
    if ap_ent:
        i_idx,j_idx  = np.triu_indices(N - M)
    else:
        i_idx,j_idx  = np.triu_indices(N - M, k=1)
    #i_idx,j_idx = i_idx.astype(np.uint16), j_idx.astype(np.uint16)

    # We start by making Cm
    Em = _embed_seq(X, 1, M)
    dif =  np.abs(Em[i_idx] - Em[j_idx])
    max_dist = np.max(dif, 1)
    inrange_cm = max_dist <= R


    in_range_i = i_idx[inrange_cm]
    in_range_j = j_idx[inrange_cm]

    if not ap_ent:
        Cm = np.bincount(in_range_i, minlength=N-M-1)
        return Cm.astype(np.float), in_range_i, in_range_j
    else:
        Cm = np.bincount(in_range_i, minlength=N-M+1)
        Cm += np.bincount(in_range_j, minlength=N-M+1)

        inrange_last = np.max(np.abs(Em[:-1] - Em[-1]),1) <= R
        Cm[inrange_last] += 1
        # all matches + self match
        Cm[-1] += np.sum(inrange_last) + 1

        return Cm.astype(np.float), in_range_i, in_range_j

def pfd(a):
    r"""
    Compute Petrosian Fractal Dimension of a time series [PET95]_.

    It is defined by:

    .. math::

        \frac{log(N)}{log(N) + log(\frac{N}{N+0.4N_{\delta}})}


    Where:

    :math:`N` is the length of the time series, and

    :math:`N_{\delta}` is the number of sign changes.


    .. [PET95]  A. Petrosian, Kolmogorov complexity of finite sequences and recognition of different preictal EEG patterns, in ,
        Proceedings of the Eighth IEEE Symposium on Computer-Based Medical Systems, 1995, 1995, pp. 212-217.



    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :return: the Petrosian Fractal Dimension; a scalar.
    :rtype: float
    """
    diff = np.diff(a)
    # x[i] * x[i-1] for i in t0 -> tmax
    prod = diff[1:-1] * diff[0:-2]

    # Number of sign changes in derivative of the signal
    N_delta = np.sum(prod < 0)
    n = len(a)
    return np.log(n)/(np.log(n)+np.log(n/(n+0.4*N_delta)))

def hjorth(X):
    r"""
    Compute Hjorth parameters [HJO70]_.


    .. math::

        Activity = m_0 = \sigma_{a}^2

    .. math::

        Morbidity = m_2 = \sigma_{d}/ \sigma_{a}

    .. math::
        Morbidity = m_4 =  \frac{\sigma_{dd}/ \sigma_{d}}{m_2}


    Where:

    :math:`\sigma_{x}^2` is the mean power of a signal :math:`x`. That is, its variance, if it's mean is zero.

    :math:`a`, :math:`d` and :math:`dd` represent the original signal, its first and second derivatives, respectively.


    .. [HJO70] B. Hjorth, EEG analysis based on time domain properties,
        Electroencephalography and Clinical Neurophysiology, vol. 29, no. 3, pp. 306-310, Sep. 1970.

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :return: activity, complexity and morbidity
    :rtype: tuple(float, float, float)
    """

    first_deriv = np.diff(X)
    second_deriv = np.diff(X,2)

    var_zero = np.mean(X ** 2)
    var_d1 = np.mean(first_deriv ** 2)
    var_d2 = np.mean(second_deriv ** 2)


    activity = var_zero
    morbidity = np.sqrt(var_d1 / var_zero)
    complexity = np.sqrt(var_d2 / var_d1) / morbidity

    return activity, morbidity, complexity

def svd_entropy(a, tau, de):
    r"""
    Compute the Single Value Decomposition entropy of a signal with embedding dimension "de" and delay "tau" [PYEEG]_.
    The result differs from PyEEG implementation because :math:`log_2` is used, according to the definition in the paper.

    .. [PYEEG] F. S. Bao, X. Liu, and C. Zhang, PyEEG: An Open Source Python Module for EEG/MEG Feature Extraction,
    Computational Intelligence and Neuroscience, vol. 2011, p. e406391, Mar. 2011.

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :param tau: the delay
    :type tau: int
    :param de: the embedding dimension
    :type de: int
    :return: the SVD entropy, a scalar
    :rtype: float
    """

    mat =  _embed_seq(a, tau, de)
    W = np.linalg.svd(mat, compute_uv = False)
    W /= sum(W) # normalize singular values
    return -1*sum(W * np.log2(W))

def fisher_info(a, tau, de):
    r"""
    Compute the Fisher information of a signal with embedding dimension "de" and delay "tau" [PYEEG]_.
    Vectorised version of the PyEEG function.

    .. [PYEEG] F. S. Bao, X. Liu, and C. Zhang, PyEEG: An Open Source Python Module for EEG/MEG Feature Extraction,
    Computational Intelligence and Neuroscience, vol. 2011, p. e406391, Mar. 2011.

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :param tau: the delay
    :type tau: int
    :param de: the embedding dimension
    :type de: int
    :return: the Fisher information, a scalar
    :rtype: float
    """

    mat =  _embed_seq(a, tau, de)
    W = np.linalg.svd(mat, compute_uv = False)
    W /= sum(W) # normalize singular values
    FI_v = (W[1:] - W[:-1]) **2 / W[:-1]

    return np.sum(FI_v)

def ap_entropy(a, m, R):
    r"""
    Compute the approximate entropy of a signal with embedding dimension "de" and delay "tau" [PYEEG]_.
    Vectorised version of the PyEEG function. Faster than PyEEG, but still critically slow.

    .. [PYEEG] F. S. Bao, X. Liu, and C. Zhang, PyEEG: An Open Source Python Module for EEG/MEG Feature Extraction,
    Computational Intelligence and Neuroscience, vol. 2011, p. e406391, Mar. 2011.

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :param m: the scale
    :type m: int
    :param R: The tolerance
    :type R: float`
    :return: the approximate entropy, a scalar
    :rtype: float
    """

    N = len(a)
    Cm, in_range_i, in_range_j = _make_cm(a,m,R)

    Cmp = _make_cmp(a, m, R, in_range_i, in_range_j)

    Cm /= float((N - m +1 ))
    Cmp /= float(N - m)

    Phi_m, Phi_mp = np.sum(np.log(Cm)),  np.sum(np.log(Cmp))
    Ap_En = (Phi_m - Phi_mp) / (N - m)
    return Ap_En

def samp_entropy(a, m, R):
    r"""
    Compute the sample entropy of a signal with embedding dimension "de" and delay "tau" [PYEEG]_.
    Vectorised version of the PyEEG function. Faster than PyEEG, but still critically slow.

    .. [PYEEG] F. S. Bao, X. Liu, and C. Zhang, PyEEG: An Open Source Python Module for EEG/MEG Feature Extraction,
    Computational Intelligence and Neuroscience, vol. 2011, p. e406391, Mar. 2011.

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :param m: the scale
    :type m: int
    :param R: The tolerance
    :type R: float`
    :return: the approximate entropy, a scalar
    :rtype: float
    """


    Cm, in_range_i, in_range_j = _make_cm(a,m,R, False)
    Cmp = _make_cmp(a, m, R, in_range_i, in_range_j, False)

    Samp_En = np.log(sum(Cm)/sum(Cmp))

    return Samp_En


def spectral_entropy(a, sampling_freq, bands=None):

    r"""
    Compute spectral entropy of a  signal with respect to frequency bands.
    The power spectrum is computed through fft. Then, it is normalised and assimilated to a probability density function.
    The entropy of the signal :math:`x` can be expressed by:

    .. math::

        H(x) =  -\sum_{f=0}^{f = f_s/2} PSD(f) log_2[PSD(f)]

    Where:

    :math:`PSD` is the normalised power spectrum (Power Spectrum Density), and

    :math:`f_s` is the sampling frequency

    :param a: a one dimensional array representing a time series
    :type a: np.ndarray
    :param sampling_freq: the sampling frequency
    :type sampling_freq:  float
    :param bands: a list of numbers delimiting the bins of the frequency bands. If None the entropy is computed over the whole range of the DFT (from 0 to :math:`f_s/2`)
    :return: the spectral entropy; a scalar
    """



    psd = np.abs(np.fft.rfft(a))**2
    psd /= np.sum(psd) # psd as a pdf (normalised to one)

    if bands is None:
        power_per_band= psd[psd>0]
    else:
        freqs = np.fft.rfftfreq(a.size, 1/float(sampling_freq))
        bands = np.asarray(bands)

        freq_limits_low = np.concatenate([[0.0],bands])
        freq_limits_up = np.concatenate([bands, [np.Inf]])

        power_per_band = [np.sum(psd[np.bitwise_and(freqs >= low, freqs<up)])
                for low,up in zip(freq_limits_low, freq_limits_up)]

        power_per_band= power_per_band[ power_per_band > 0]

    return - np.sum(power_per_band * np.log2(power_per_band))


def dfa(X, Ave = None, L = None, sampling= 1):
    X = np.array(X)
    if Ave is None:
        Ave = np.mean(X)
    Y = np.cumsum(X)
    Y -= Ave
    if not L:
        max_power = np.int(np.log2(len(X)))-4
        L = X.size / 2 ** np.arange(4,max_power)
    F = np.zeros(len(L)) # F(n) of different given box length n

    for i,n in enumerate(L):
        sampled = 0
        for j in xrange(0,len(X) -n ,n):

            if np.random.rand() < sampling:
                F[i] += np.polyfit(np.arange(j,j+n), Y[j:j+n],1, full=True)[1]
                sampled += 1
        if sampled > 0:
            F[i] /= float(sampled)

    LF = np.array([(l,f) for l,f in zip(L,F) if l>0]).T

    F = np.sqrt(LF[1])

    Alpha = np.polyfit(np.log(LF[0]), np.log(F),1)[0]
    return Alpha

def hurst(signal):
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



#
#
# def hfd_new(X, Kmax):
#     """ Compute Higuchi Fractal Dimension of a time series X, kmax
#      is an HFD parameter
#     """
#     L = []
#     x = []
#     N = len(X)
#     for k in range(1,Kmax):
#         Lk = []
#         for m in range(k):
#             max = int((N-m)/k)
#
#             diffs = X[m+k:len(X)-k+1:k]
#             #print np.sum(np.abs(diffs))
#             # print len(ar[m: :k])
#
#
#             print k,m,len(diffs)
#
#             Lmk = np.sum(np.abs(diffs))*(N - 1)/np.floor((N - m) / float(k)) / k
#             #print k,m, np.sum(np.abs(diffs)) / float(len(diffs))
#             Lk.append(Lmk)
#
#
#
#         L.append(np.log(np.mean(Lk)))
#         x.append([np.log(float(1) / k), 1])
#
#     # print L,x
#     (p, r1, r2, s)=np.linalg.lstsq(x, L)
#     return p[0]
#
#



def hfd(X, Kmax):
    """ Compute Higuchi Fractal Dimension of a time series X, kmax
     is an HFD parameter
    """
    L = []
    x = []
    N = len(X)
    for k in xrange(1,Kmax):
        Lk = []
        for m in xrange(0,k):
            Lmk = 0
            test = []

            for i in xrange(1,int(np.floor((N-m)/k))):
                test.append(X[m+i*k])
                Lmk += abs(X[m+i*k] - X[m+i*k-k])
            #print k,m, Lmk
            print k,m,len(test)
            # print k,m,Lmk
            Lmk = Lmk*(N - 1)/np.floor((N - m) / float(k)) / k
            Lk.append(Lmk)


        L.append(np.log(np.mean(Lk)))
        x.append([np.log(float(1) / k), 1])
    (p, r1, r2, s)=np.linalg.lstsq(x, L)
    return p[0]

#
# ar = np.arange(100)
# hfd(ar, 10)
# print "==================="
# hfd_new(ar, 10)

# def hfd_exple(pow=10):
#     start = 1000
#     N = np.int64(2 **pow)
#     i=start
#     out = np.zeros((N))
#     while i < N + start:
#         out[i-start] = np.sum(np.random.normal(0,1,i))
#         i += 1
#     return out
#
# a = hfd_exple(15)


# hfd(a, 4)






