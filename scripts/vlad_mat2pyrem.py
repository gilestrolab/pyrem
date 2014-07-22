__author__ = 'quentin'


import scipy.io as scio
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
import pyrem as pr

VIGILANCE_MAP = {"w1":"D", "r3":"D", "nr2":"D" ,
                 "w":"W", "mt":"M", "r":"R", "nr":"N"}

ANNOTATION_VALUES = {
    "w": ord("W") + 1j,
    "W": ord("W") + 1j,
    "n": ord("N") + 1j,
    "nr": ord("N") + 1j,
    "N": ord("N") + 1j,
    "r": ord("R") + 1j,
    "R": ord("R") + 1j,
    "m": ord("M") + 1j,
    "mt": ord("M") + 1j,
    "M": ord("M") + 1j,
    "d": 0 + 0j,
    "w1": 0 + 0j,
    "nr2": 0 + 0j,
    "r3": 0 + 0j,
    "D": 0 + 0j,
    "?": 0 + 0j }

MATLAB_FILE_SUFFIX = {
    "EMG":"EMG",
    "f":"Frontal_EEG",
    "o":"Occipital_EEG",
}

ANNOTATION_SUFFIX = "o-VSspec"

def make_annotations(data_dic):
    max = 0
    for k in ANNOTATION_VALUES.keys():
        try:
            tmp_max = np.max(data_dic[k])

            if tmp_max > max:
                max = tmp_max
        except Exception as e:
            #print e
            pass
    annotations = np.zeros(max, dtype=np.complex64) + ANNOTATION_VALUES["d"]


    for k,item in ANNOTATION_VALUES.items():
        try:
            data = data_dic[k]
            if len(data) > 0:
                annotations[data -1 ] = item
        except Exception as e:
            pass


    return np.reshape(annotations, (annotations.size,1))


def load_from_mat(file_name_prefix ,sampling_rate=256, metadata={}):
    """

    :param file:the matlab file name
    :return: a polygraph
    """


    file_name_annotations = file_name_prefix + "-" + ANNOTATION_SUFFIX + ".mat"

    matl = scio.loadmat(file_name_annotations, squeeze_me=True, struct_as_record=False)
    annot_values = make_annotations(matl)
    signals = []
    good_keys = []
    for k, item in MATLAB_FILE_SUFFIX.items():
        file_name = file_name_prefix + "-" + k + ".mat"
        matl = scio.loadmat(file_name, squeeze_me=True, struct_as_record=False)
        signals.append(matl["resampled_sig"].astype(np.float32))
        good_keys.append(item)
    data = np.array(signals).T

    #
    x = np.linspace(0,data.shape[0],annot_values.shape[0])
    inter_f = interp1d(x, annot_values, "nearest", axis=0)
    annot_values = inter_f(np.arange(data.shape[0]))

    metadata["input_file"] = file_name_prefix

    polygraph =  pr.Polygraph(data,sampling_rate,annot_values,good_keys,annotation_types=["vigilance"], metadata=metadata)

    return polygraph

#
# data = load_from_mat("/stk/VM_SHARE_FOLDER/test2.mat")
#
# for k,i in matl.items():
#     try:
#         print np.min(i), k
#     except:
#         pass



