from pyrem.polygram import Polygram

__author__ = 'quentin'

#todo edf, csv



import scipy.io as scio
# from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
from pyrem.time_series import Signal, Annotation
import  joblib

def polygram_from_pkl(filename):
    return joblib.load(filename)

def polygram_from_spike_matlab_file(filename, fs, annotation_fs, channel_names, channel_types, doubt_chars,resample_signals, metadata={}):

    """
    This function loads a matlab file exported by spike to
    as a polygraph.

    :param file:the matlab file name
    :return: a polygram
    """

    type_for_name = dict(zip(channel_names, channel_types))
    matl = scio.loadmat(filename, squeeze_me=True, struct_as_record=False)

    data_channels = {}
    annotation_channels = {}

    for k in matl.keys():
        # exclude metadata such as "__global__", "__version__" ...
        if not k.startswith("__"):
            obj = matl[k]
            channel_number = int(k.split("_")[-1][2:])
            if "values" in dir(obj):
                channel_id  = channel_names[channel_number-1]
                data_channels[channel_id] = obj.values
            elif "text" in dir(obj):
                annotation_channels["Stage_%i" % (channel_number-1)] = obj.text
    del matl

    crop_at = np.min([i.size for _,i in data_channels.items()])
    print "crop_at =", crop_at, "max size was:" , np.max([i.size for _,i in data_channels.items()])
    for k,a in data_channels.items():
        data_channels[k] = a[:crop_at]

    signals = [Signal(data,fs, name=name, type=type_for_name[name] ) for name,data in data_channels.items()]
    del data_channels


    if resample_signals:
        signals = [ s.resample(resample_signals) for s in signals]

    annotations = pd.DataFrame(annotation_channels)
    del annotation_channels
    # trick to remove every second annotation (they are unpredictably redundant)
    annotations = annotations[1:annotations.shape[0]:2]
    #annotations2 = annotations[0:annotations.shape[0]:2]

    # check we have removed the right (mainly empty) annotations
    # print annotations.shape[0] *5.0 / 60.0 / 60.0
    # print annotations2.shape[0] *5.0 / 60.0 / 60.0




    annotations[annotations[annotations.columns[0]] == ""] = "?"


    np_ord = np.vectorize(lambda x : ord(x.upper()))

    annot_values = np_ord(np.array(annotations)).flatten()
    #
    annot_values = annot_values.astype(np.uint8)

    # last annotation is truncated

    annot_values = annot_values[:annot_values.size -1]


    annot_probas = [0 if a in doubt_chars else 1 for a in annot_values]
    print "000000000000", np.sum(np.array(annot_probas) == 0)
    an = Annotation(annot_values, annotation_fs, annot_probas, name="vigilance_state", type="vigilance_state")


    signals = [s[:an.duration]for s in signals]

    print an.duration
    # print signals[0].duration.total_seconds(), an.duration.total_seconds()
    # print 200  - 200*(signals[0].duration.total_seconds() - an.duration.total_seconds() )/signals[0].duration.total_seconds()

    signals.append(an)
    return Polygram(signals)

    #

    #
    #
    #
    # x = np.linspace(0,data.shape[0],annot_values.shape[0])
    #
    # inter_f = interp1d(x, annot_values, "nearest", axis=0)
    #
    # annot_values = inter_f(np.arange(data.shape[0]))
    #
    # metadata["input_file"] = filename
    #
    #
    # # polygraph =  pr.Polygraph(data, sampling_rate,
    # #                           annot_values, df.columns,
    # #                           annotation_types=["vigil"], metadata=metadata)
    #

#    return polygraph






