__author__ = 'quentin'


import scipy.io as scio
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
import pyrem as pr

DOUBT_CHARS = {ord("?"), ord("D")}
def load_from_mat(file_name,sampling_rate=200, metadata={}):
    """
    This function loads a matlab file exported by spike to
    as a polygraph.

    :param file:the matlab file name
    :return: a polygraph
    """


    matl = scio.loadmat(file_name, squeeze_me=True, struct_as_record=False)
    data_channels = {}
    annotation_channels = {}

    for k in matl.keys():
        # exclude metadata such as "__global__", "__version__" ...
        if not k.startswith("__"):

            obj = matl[k]
            channel_id = int(k.split("_")[-1][2:])
            if "values" in dir(obj):
                data_channels[channel_id] = obj.values
            elif "text" in dir(obj):
                annotation_channels[channel_id] = obj.text

    data = np.array(pd.DataFrame(data_channels))
    annotations = pd.DataFrame(annotation_channels)

    annotations = annotations[annotations[annotations.columns[0]] != ""]
    #return annotations
    np_ord = np.vectorize(lambda x : ord(x.upper()))

    annot_values = np_ord(np.array(annotations))
    # return annot_values


    annot_values = annot_values.astype(np.complex64)
    annot_values += 1.0j

    for d in DOUBT_CHARS:
        print np.sum(np.real(annot_values) == d)
        annot_values[np.real(annot_values) == d] = 0.0+0.0j

    x = np.linspace(0,data.shape[0],annot_values.shape[0])
    inter_f = interp1d(x, annot_values, "nearest", axis=0)
    annot_values = inter_f(np.arange(data.shape[0]))

    metadata["input_file"] = file_name

    polygraph =  pr.Polygraph(data,sampling_rate,annot_values,data_channels.keys(),annotation_types=None, metadata=metadata)


    return polygraph

data = load_from_mat("/stk/VM_SHARE_FOLDER/test2.mat")


for t,d in data.embed_seq(5,1):
    for a in d.annotations():
        r = np.real(a)
        uniqs = np.unique(r)
        i = np.imag(a)
        probs = []
        for u in uniqs:
            eqs = (r == u)
            probs.append(np.sum(i[eqs]))

        probs = np.array(probs)
        probs /= np.sum(i)
        # if len(uniqs) >1:
        #     print uniqs
        #     print probs





