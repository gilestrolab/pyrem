
__author__ = 'quentin'

import pywt
import numpy as np

from pyrem.polygram import Polygram


def decompose_signal(signal, levels=(1,2,3,4,5), wavelet="db4", resample_before=None, mode="per", keep_a=True):
    """
    A wrapper around :func:`~pywt.wavedec`. It performs discrete wavelet decomposition and return the coefficients as a :class:`pyrem.polygram.Polygram`.
    It allows to select the wavelet coefficients to keep and can perform preliminary resampling of the signal.
    In the resulting polygram, the names of the coefficients will be suffixed by an an identifier describing their respective levels (i.e. cD_1, cD_2, ..., cD_N, cA_N).
    The sampling frequency of coefficients will also be automatically computed.

    >>> import numpy as np
    >>> import pyrem as pr

    >>> noise = np.random.normal(size=int(1e6))
    >>> sig = pr.time_series.Signal(noise, 256.0,
    >>>      name="channel_1")
    >>>      name="channel_1")
    >>> pol = pr.wavelet_decomposition decompose_signal(sig)
    >>> pol
        Polygram
        ----
        Duration:	1:05:06.250000 (HH:mm:ss)
        N signals:	6
        N annotations:	0
        Metadata:
                None
        ----
        Channel information:
                     Name  Type fs(Hz)        Duration
        0  channel_1_cA_5  None    8.0  1:05:06.250000
        1  channel_1_cD_1  None  128.0  1:05:06.250000
        2  channel_1_cD_2  None   64.0  1:05:06.250000
        3  channel_1_cD_3  None   32.0  1:05:06.250000
        4  channel_1_cD_4  None   16.0  1:05:06.250000
        5  channel_1_cD_5  None    8.0  1:05:06.250000

    :param signal: the time series to be decomposed
    :type signal: :class:`pyrem.time_series.Signal`

    :param levels: the levels to keep (e.g. [1,2,3,5]) will not return the coeficient cD_4
    :type levels: list([int])
    :param wavelet: the type of wavelet (see :func:`~pywt.wavedec`)
    :type wavelet: str
    :param resample_before: the sampling frequency at which to resample the signal before the decomposition
    :type resample_before: float

    :param mode: (see :func:`~pywt.wavedec`)
    :type mode: str

    :param keep_a: whether to keep the last coefficient (cA_N)

    :return: A polygram with all the requested coefficients
    :rtype: :class:`pyrem.polygram.Polygram`
    """
    max_level = max(levels)

    if resample_before:
        tmp_signal = signal.resample(resample_before)
    else:
        tmp_signal = signal

    # this is where the magic happens
    coeffs = pywt.wavedec(tmp_signal, wavelet, level=max_level, mode=mode)


    coeff_names = ["_cA_%i" % (len(coeffs) -1 )]

    coeff_levels = [(len(coeffs) -1 )]
    for n in range(len(coeffs) -1, 0, -1):
        coeff_levels.append(n)
        coeff_names.append("_cD_%i" % (n))


    coeff_fs = tmp_signal.fs / 2.0 ** np.arange(max_level, 0,-1)
    coeff_fs = np.concatenate([[coeff_fs[0]],coeff_fs])
    coeff_names = [signal.name + n for n in coeff_names]
    signals = []
    for l, fs, n, sig in zip(coeff_levels, coeff_fs, coeff_names, coeffs):
        if l in levels:
            signals.append(signal._copy_attrs_to_array(sig, fs=fs, name=n))
    #trim channels that are too long
    if not keep_a:
        signals = signals[1:]

    min_duration = min([s.duration for s in signals])
    signals = [s[:min_duration] for s in signals]

    return Polygram(signals)
