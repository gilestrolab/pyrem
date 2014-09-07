from pyrem.polygram import Polygram

__author__ = 'quentin'

#todo edf, csv



import scipy.io as scio
import pandas as pd
import numpy as np
from pyrem.time_series import Signal, Annotation
import joblib

from pyrem.time_series import BiologicalTimeSeries, Annotation, Signal


def polygram_from_pkl(filename):
    return joblib.load(filename)

def signal_from_pkl(filename):
    return joblib.load(filename)

#
# def signal_from_csv(file_name, sampling_freq):
#     data = pd.read_csv(file_name, engine="c", header=None, dtype=np.float32)
#     return Signal(data, sampling_freq)


def _annotation_from_spike_txt(filename, doubt_chars):
    """
    Hacky parser needed because spike does not export epochs of exactly 5.000 s!!!!!!!

    :param filename:
    :return:
    """
    df = pd.read_csv(filename, skiprows=16, sep="\t", names=["t", "x1", "x2", "x3","x4", "x5","y"], header=None, na_values="nan")
    ts = df["y"]
    ts.index = pd.to_datetime(df["t"] * 1e9)
    annotations  = ts.dropna()
    annotations = annotations.resample("5s", how="first")
    annotations = annotations.fillna(method="ffill")


    np_ord = np.vectorize(lambda x : ord(x.upper()))

    annot_values = np_ord(np.array(annotations)).flatten()

    annot_values = annot_values.astype(np.uint8)

    annot_probas = [0 if a in doubt_chars else 1 for a in annot_values]

    return  Annotation(annot_values, 1/5.0, annot_probas, name="vigilance_state", type="vigilance_state")



def polygram_from_spike_matlab_file(signal_filename, annotation_filename, fs, annotation_fs, channel_names, channel_types, doubt_chars,resample_signals, metadata={}):

    """
    This function loads a matlab file exported by spike to
    as a polygraph.

    :param signal_filename: the matlab file name
    :return: a polygram
    """

    an = _annotation_from_spike_txt(annotation_filename, doubt_chars)


    type_for_name = dict(zip(channel_names, channel_types))
    matl = scio.loadmat(signal_filename, squeeze_me=True, struct_as_record=False)

    data_channels = {}
    # annotation_channels = {}

    for k in matl.keys():
        # exclude metadata such as "__global__", "__version__" ...
        if not k.startswith("__"):
            obj = matl[k]
            channel_number = int(k.split("_")[-1][2:])
            if "values" in dir(obj):
                channel_id  = channel_names[channel_number-1]
                data_channels[channel_id] = obj.values
            elif "text" in dir(obj):
                pass
                # annotation_channels["Stage_%i" % (channel_number-1)] = obj.text
    del matl

    crop_at = np.min([i.size for _,i in data_channels.items()])

    for k,a in data_channels.items():
        data_channels[k] = a[:crop_at]

    signals = [Signal(data,fs, name=name, type=type_for_name[name] ) for name,data in data_channels.items()]
    del data_channels


    if resample_signals:
        signals = [ s.resample(resample_signals) for s in signals]



    #signals = [s[:an.duration]for s in signals]



    signals.append(an)
    return Polygram(signals)

