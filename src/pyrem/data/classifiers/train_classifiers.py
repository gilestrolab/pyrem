from pyrem.classifiers.eeg_vs_emg import EEGsvEMG
import pyrem as pr
import glob
DATA_FILE_PATTERN= "/data/pyrem/eeg_vs_emg_data/*.pkl"

from sklearn.externals import joblib

if __name__ == "__main__":

    classif = EEGsvEMG()
    classif.train_from_polygraph_file_list(DATA_FILE_PATTERN)
    classif.save("EEGsvEMG.pkl")

    # # prediction verification:

    # files = glob.glob(DATA_FILE_PATTERN)
    # for f in sorted(files):
    #     print "Predicting: " + f
    #
    #     polyg = pr.polygraph_from_pkl(f)
    #     a, proba = classif.predict(polyg)
    #     print a, proba
    # #




