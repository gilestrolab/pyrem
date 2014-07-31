__author__ = 'quentin'

import pywt
import numpy as np
from pyrem.signal.signal import Signal
from pyrem.signal.polygram import Polygram

def decompose_signal(signal, wavelet="db4", levels=[1,2,3,4,5]):
    max_level = max(levels)

    coeffs = pywt.wavedec(signal, wavelet, level=max_level)


    coeff_names = ["_cA_%i" % (len(coeffs) -1 )]
    coeff_levels = [(len(coeffs) -1 )]
    for n in range(len(coeffs) -1, 0, -1):
        coeff_levels.append(n)
        coeff_names.append("_cD_%i" % (n))


    coeff_fs = signal.fs / 2.0 ** np.arange(max_level, 0,-1)
    coeff_fs = np.concatenate([[coeff_fs[0]],coeff_fs])

    signals = []
    for l, fs, n, sig in zip(coeff_levels, coeff_fs, coeff_names, coeffs):
        if l in levels:
            signals.append(Signal(sig,fs, name=n))
    #trim channels that are too long
    min_duration = min([s.duration for s in signals])
    signals = [s[:min_duration] for s in signals]

    return Polygram(signals)
