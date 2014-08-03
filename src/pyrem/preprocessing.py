__author__ = 'quentin'

import numpy as np
import scipy.signal as sig

def _butter_bandpass(lowcut, highcut, fs, order=5,typ="band"):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = sig.butter(order, [low, high], btype=typ)
    return b, a


def _butter_bandpass_filter(data, lowcut, highcut, fs, order=5, typ="band"):
    b, a = _butter_bandpass(lowcut, highcut, fs, order=order, typ=typ)
    y = sig.lfilter(b, a, data)
    return y



def preprocess_eegs(pol, low_cut=0.5, high_cut=60, stop_cut=50.0, order=6):

    stop_low = stop_cut - 1
    stop_high = stop_cut + 1

    # "EEG" string must be in the name of the channel
    #eegs_channels = [i for i,c in enumerate(pol.channel_types)]

    out_data = []
    for i,c in enumerate(pol.channel_types):
        data = pol[i].data.flatten()
        if "EEG" in c.upper():
            data  = _butter_bandpass_filter(data  , low_cut, high_cut, pol.sampling_freq, order=6)
            data  = _butter_bandpass_filter(data  , stop_low,stop_high, pol.sampling_freq, order=6,typ="bandstop")
        out_data.append(data)

    return  pol.copy(np.array(out_data).T)
