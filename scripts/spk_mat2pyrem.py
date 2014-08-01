__author__ = 'quentin'
from pyrem.signal.io import polygram_from_spike_matlab_file
from multiprocessing import pool
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
# INPUT_PATTERN = "/data/pyrem/Ellys/mats/GFP_F*.mat"
OUT_DIR = "/data/pyrem/Ellys/pkls"


RESAMPLE_AT = 256
DEFAULT_FS  = 200
FSS = {
"GFP_sleep":199.92,
"GFP_A":200.11,
"GFP_D":199.98,
"GFP_F":199.9,
"GFP_1":199.98,
"TelC_C":199.99,
}

EXCLUDED = {
"TelC_E"

}
N_PROCESS = 3

def save_one_pol(f):

    basename = os.path.basename(f).split(".")[0]
    if basename in EXCLUDED:
        return

    new_file_name = basename + ".pkl"
    out_path = os.path.join(OUT_DIR,new_file_name)

    print "=======================\nConverting: " + f + "\nInto: "+ out_path
    try:

        fs = FSS[basename]
        print "applying custom fs:", fs

    except KeyError as k:

        fs = DEFAULT_FS
    try:
        pol = polygram_from_spike_matlab_file(f, fs, 1/5.0,CHANNEL_ID_MAP,DOUBT_CHARS, resample_signals=RESAMPLE_AT)



        pol.save(out_path)
        print "Saved"
    except Exception as e:
        print "ERRORRRRRRRRRR!", e
        pass

if __name__== "__main__":
    files = glob.glob(INPUT_PATTERN)
    clust = pool.Pool(N_PROCESS)

    clust.map(save_one_pol, sorted(files))
