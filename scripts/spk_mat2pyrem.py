__author__ = 'quentin'
from pyrem.signal.io import polygram_from_spike_matlab_file
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




if __name__== "__main__":
    files = glob.glob(INPUT_PATTERN)
    for f in sorted(files):
        new_file_name = os.path.basename(f).split(".")[0] + ".pkl"
        out_path = os.path.join(OUT_DIR,new_file_name)
        print "=======================\nConverting: " + f + "\nInto: "+ out_path
        try:
            pol = polygram_from_spike_matlab_file(f, 200.0, 1/5.0,CHANNEL_ID_MAP,DOUBT_CHARS)



            pol.save(out_path)
            print "Saved"
        except Exception as e:
            print e
            pass