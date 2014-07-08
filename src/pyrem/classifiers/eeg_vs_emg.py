__author__ = 'quentin'


import numpy as np
import pandas as pd
import glob
import cPickle
from sklearn.preprocessing import LabelEncoder
from scipy.stats import mode
from pyrem.features.summary import PowerFeatures
from pyrem.features.periodogram import PeriodFeatures
from pyrem.signal.signal import load_signal
from sklearn.ensemble import RandomForestClassifier

class EEGsvEMG(object):
    __FEATURE_GROUPS = [
            PowerFeatures(),
            PeriodFeatures()
        ]


    def __init__(self, data_file_pattern):
        self._features = pd.DataFrame()
        self.classifier = RandomForestClassifier(oob_score=True, n_estimators=500, n_jobs=8)
        self.label_encoder = LabelEncoder()
        self.n_chunks_for_predict = 25
        files = glob.glob(data_file_pattern)
        for f in sorted(files):
            print "Processing: " + f

            signal = load_signal(f)
            labels = signal.channel_types

            for t,subsignal in signal.embed_seq(30,1):
                self.train(labels, subsignal.signal_iter(), subsignal)

        self.complete_training()


    def make_features(self, data_iterator, data, processes=1):
        all_dfs = []
        for  channel  in data_iterator:
            groups_dfs = [group(channel, data) for group in self.__FEATURE_GROUPS]
            all_dfs.append(pd.concat(groups_dfs, axis=1))

        features = pd.concat(all_dfs)

        return  features

    def save(self, filename):
        cPickle.dump(self, Popen(filename, "w"), cPickle.HIGHEST_PROTOCOL)

    def train(self, labels, data_iterator, data):
        new_features_df = self.make_features(data_iterator, data)
        new_features_df.index = labels
        self._features = pd.concat([self._features, new_features_df])
        #todo append = cleanner

    def complete_training(self):
        x = self._features
        y = self.label_encoder.fit_transform(self._features.index)
        # drop features to save memory ;)
        self._features = None
        self.classifier.fit(x,y)

    def predict (self,  data_iterator, data):
        x = self.make_features(data_iterator, data)
        y = self.classifier.predict(x)
        return y
        #return self.label_encoder.inverse_transform(y)

    def predict_samples(self, signal):
        predictions = []
        n_chunks = float(signal.ntimepoints) / signal.sampling_freq / 30.0
        p =  n_chunks / self.n_chunks_for_predict
        if p < 1:
            p = 1
        else:
            p = int(p)
        for i,(t,subsignal) in enumerate(signal.embed_seq(30,1)):
            if i % p == 0:
                predictions.append(self.predict(subsignal.signal_iter(), subsignal))
        all_predicts = np.array(predictions)
        mo, n = mode(all_predicts)
        prop = n / all_predicts.shape[0]

        annots = []

        for p,m in zip(prop.flatten(), mo.flatten()):
            if p > 0.75: #fixme magic number
                annots.append(self.label_encoder.inverse_transform(int(m)))
            else:
                annots.append(None)

        return annots