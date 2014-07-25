__author__ = 'quentin'


import scipy.io as scio
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
import pyrem as pr
import glob
import os

DOUBT_CHARS = {ord("?"), ord("D")}
CHANNEL_ID_MAP = [
                "EEG_parietal_cereb",
                "EEG_parietal_frontal",
                "EMG_1",
                "EMG_2"
                  ]
INPUT_PATTERN = "/data/pyrem/Ellys/mats/*.mat"
# INPUT_PATTERN = "/data/pyrem/Ellys/mats/GFP*A*.mat"
OUT_DIR = "/data/pyrem/Ellys/pkls"

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
            channel_number = int(k.split("_")[-1][2:])


            if "values" in dir(obj):
                channel_id  = CHANNEL_ID_MAP[channel_number-1]
                data_channels[channel_id] = obj.values
            elif "text" in dir(obj):

                annotation_channels["Stage_%i" % (channel_number-1)] = obj.text

    crop_at = np.min([i.size for _,i in data_channels.items()])

    for k,a in data_channels.items():
        data_channels[k] = a[:crop_at]



    df = pd.DataFrame(data_channels)

    data = np.array(df)
    annotations = pd.DataFrame(annotation_channels)
    # trick to remove every second annotation (they are unpredictibly rudundant)
    annotations = annotations[1:annotations.shape[0]:2]
    annotations2 = annotations[0:annotations.shape[0]:2]

    # check we have removed the right (mainly empty) annotations
    print annotations.shape[0] *5.0 / 60.0 / 60.0
    print annotations2.shape[0] *5.0 / 60.0 / 60.0

    annotations = annotations[annotations[annotations.columns[0]] != ""]

    np_ord = np.vectorize(lambda x : ord(x.upper()))

    annot_values = np_ord(np.array(annotations))



    annot_values = annot_values.astype(np.complex64)
    annot_values += 1.0j

    for d in DOUBT_CHARS:
        annot_values[np.real(annot_values) == d] = 0.0+0.0j

    x = np.linspace(0,data.shape[0],annot_values.shape[0])

    inter_f = interp1d(x, annot_values, "nearest", axis=0)

    annot_values = inter_f(np.arange(data.shape[0]))

    metadata["input_file"] = file_name

    polygraph =  pr.Polygraph(data, sampling_rate,
                              annot_values, df.columns,
                              annotation_types=["vigil"], metadata=metadata)


    return polygraph

# data = load_from_mat("/data/pyrem/Ellys/mats/GFP_2.mat")
# #
if __name__== "__main__":
    files = glob.glob(INPUT_PATTERN)
    for f in sorted(files):

        pol = load_from_mat(f)
        new_file_name = os.path.basename(f).split(".")[0] + ".pkl"
        out_path = os.path.join(OUT_DIR,new_file_name)
        print "Converting: " + f + "\nInto: "+ out_path
        pol.save(out_path)



    #
    #
    #
    # for t,d in data.embed_seq(5,1):
    #     for a in d.annotations():
    #         r = np.real(a)
    #         uniqs = np.unique(r)
    #         i = np.imag(a)
    #         probs = []
    #         for u in uniqs:
    #             eqs = (r == u)
    #             probs.append(np.sum(i[eqs]))
    #
    #         probs = np.array(probs)
    #         probs /= np.sum(i)
            # if len(uniqs) >1:
            #     print uniqs
            #     print probs





