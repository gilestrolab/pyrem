__author__ = 'quentin'
import glob
import os
import pandas as pd

from pyrem.feature_families import *
from pyrem.wavelet_decomposition import decompose_signal
from pyrem.io import polygram_from_pkl
from scipy.signal import welch
from multiprocessing import Pool
import numpy as np



DATA_FILE_PATTERN= "/data/pyrem/Ellys/pkls/GFP_*.pkl"

WINDOW_SIZE = 5
WINDOW_LAG = 100.0


def make_periodogram(a):
    f, Pxx_den = welch(a, a.fs, nperseg=512)

    den = np.log10(Pxx_den)
    return den

import pylab as pl
def data_for_one_file(file, channel_name, dfs):
    pol = polygram_from_pkl(file)

    for t, w in pol.iter_window(WINDOW_SIZE, WINDOW_LAG):

        eeg = w[channel_name]
        ann = w["vigilance_state"]
        if ann.probas.all() >0:
            y = ann.values[0]
            periodo = make_periodogram(eeg)

            try:
                dfs[y].append(periodo)
            except KeyError:
                dfs[y]=[periodo]
    return dfs

files = glob.glob(DATA_FILE_PATTERN)
dfs = {}
for f in sorted(files):
    print f
    dfs = data_for_one_file(f, "EEG_parietal_frontal", dfs)
    # dfs = data_for_one_file(f, "EMG_REF", dfs)

out = {}
cols = {78:"r", 82:"m", 87:"b"}
for y,ar in dfs.items():

    q1,q2,q3  = np.percentile(ar,[25,50,75], 0)
    #stds = np.std(np.array(ar),0 )
    #out[y] =stds

    print y,len(ar)
    x = np.arange(0, q1.size)

    pl.fill_between(x,q3,q1,color=cols[y], alpha=0.15)
    pl.plot(q2,color=cols[y])
#pl.show()
#
# hlines = 160.0 / (2.0 **np.arange(1,7))
hlinesold = 256/ (2.0 **np.arange(1,7))
# pl.vlines(hlines, -100, +100)
pl.vlines(hlinesold, -100, +100, linestyle="dashed")
pl.show()



# axarr[0].plot(x, y)
# axarr[0].set_title('Sharing X axis')
# axarr[1].scatter(x, y)
posit= {78:0, 82:1, 87:2}
for y,ar in dfs.items():
    npar = 1/-np.array(ar)
    npar = npar[:,0:25]

    for a in npar:
        pl.plot(a,color=cols[y], alpha=0.1)

#
# vlines = 160.0 / (2.0 **np.arange(1,7))
# pl.vlines(vlines, -100, +100)

pl.show()



outl = []
ys = []
for y,ar in dfs.items():
    npar = 1/-np.array(ar)
    npar = npar[:,0:25]
    sums = np.sum(npar,1)
    npar = npar/ sums[:,None]


    for a in npar:
        r = np.arange(0, a.size)
        m = np.sum(a * r)

        h = -np.sum(a * np.log2(a))
        order = np.argsort(-a)
        ranks = order.argsort()
        #print y, order
        ys.append(y)
        outl.append(
           # {
           #      "y":y,
           #      "mod":np.argmax(a),
           #      "h":h,
           #      "mean":m}
            order
        )


out_df = pd.DataFrame(np.array(outl))
out_df["y"] = ys
gr = out_df.groupby("y")

for y,a in zip(ys, outl):
    pl.plot(a,color=cols[y], alpha=0.1)
