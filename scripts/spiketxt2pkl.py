__author__ = 'quentin'

import sys
import numpy as np
import pandas as pd
import pyrem as pr

CHANNELS = ['e','e','m','m']
METADATA = {"organism":"Mus musculus", "sex":"M", "author":"Valentina Ferretti"}
SAMPLING_FREQ = 200

if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    df = pd.read_csv(input, sep="\t")
    df = np.array(df.drop("Time",axis=1).dropna())
    met = METADATA
    met["raw_file_name"] = input
    sig = pr.Signal(df,SAMPLING_FREQ, metadata=METADATA, channel_types=CHANNELS, normalise=False)
    sig.save(output)
