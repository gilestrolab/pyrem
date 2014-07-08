from pyrem.classifiers.eeg_vs_emg import EEGsvEMG
DATA_FILE_PATTERN= "/stk/eeg_vs_emg_data/*.pkl"




if __name__ == "__main__":

    classif = EEGsvEMG(DATA_FILE_PATTERN)
    classif.save("EEGsvEMG.pkl")

    # files = glob.glob(DATA_FILE_PATTERN)
    # for f in sorted(files):
    #     print "Predicting: " + f
    #
    #     signal = load_signal(f)
    #     a = classif.predict_samples(signal)
    #




