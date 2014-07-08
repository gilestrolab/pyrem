
import glob
from pyrem.classifiers.eeg_vs_emg import EEGsvEMG
from pyrem.signal.signal import signal_from_csv
DATA_FILE_PATTERN= "/stk/eeg_vs_emg_data/*.txt"
SAMPLING_FREQ = 200.0

if __name__ == "__main__":
    files = glob.glob(DATA_FILE_PATTERN)
    classif = EEGsvEMG()
    for f in sorted(files):
        print "Processing: " + f
        signal = signal_from_csv(f, SAMPLING_FREQ)
        labels = ["a", "b", "c", "d"]

        for t,subsignal in signal.embed_seq(30,1):
            print t
            classif.train(labels, subsignal.signal_iter(), subsignal)
        #todel

        classif._features.to_csv("/home/quentin/Desktop/eeg_emg_test.csv", float_format="%e")





