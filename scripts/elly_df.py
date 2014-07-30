__author__ = 'quentin'

import glob
import os
import pandas as pd
import pyrem as pr
from multiprocessing import Pool
import numpy as np


DATA_FILE_PATTERN= "/data/pyrem/Ellys/pkls/*.pkl"
SAMPLING_RATE = 200
OUT_CSV = "/data/pyrem/Ellys/all_features.csv"
LAG_WINDOW = 5
N_PROCESSES = 6

# DATA_FILE_PATTERN= "/data/pyrem/Ellys/pkls/GFP_*A*.pkl"
# OUT_CSV = "/tmp/all_features.csv"




def features_one_file(f):
    file_name = os.path.basename(f).split(".")[0]
    treatment, animal = file_name.split("_")

    pol = pr.polygraph_from_pkl(f)
    pol = pol.normalise()
    pol = pr.preprocess_eegs(pol)
    print "processing " + f

    tmp_df = feature_factory.make_features_for_epochs(pol,10,LAG_WINDOW, add_major_annotations=True)

    tmp_df["animal"] = animal
    tmp_df["treatment"] = treatment


    return tmp_df

if __name__ == "__main__":

    files = glob.glob(DATA_FILE_PATTERN)

    feature_factory = pr.features.FeatureFactory([
        pr.features.PeriodFeatures(),
        pr.features.PowerFeatures(),
        pr.features.NonLinearFeatures(),
        # pr.features.EntropyFeatures(),
        pr.features.HjorthFeatures(),
        pr.features.WaveletsFeaturesDB4(),
        pr.features.MSEFeatures(),


    ])
    if N_PROCESSES > 1 :
        p = Pool(N_PROCESSES)
        dfs = p.map(features_one_file, sorted(files))
    else:

        dfs = map(features_one_file, sorted(files))

    out_df = pd.concat(dfs)

    out_df.to_csv(OUT_CSV, float_format="%e")
