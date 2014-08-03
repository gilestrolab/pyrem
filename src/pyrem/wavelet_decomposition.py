__author__ = 'quentin'

import pywt
import numpy as np

from pyrem.polygram import Polygram


def decompose_signal(signal, levels=(1,2,3,4,5), wavelet="db4", resample_before=None, mode="per", keep_a=True):
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
