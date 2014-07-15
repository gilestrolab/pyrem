__author__ = 'quentin'


import glob
from sklearn.externals import joblib as pkl
from sklearn.preprocessing import LabelEncoder
from pyrem.features.feature_families import *

from pyrem.signal.polygraph import polygraph_from_pkl
from sklearn.ensemble import RandomForestClassifier

class EEGsvEMG(object):
    __FEATURE_GROUPS = [
            PowerFeatures(),
            HjorthFeatures()
        ]


    def __init__(self):
        self._features = pd.DataFrame()
        self.classifier = RandomForestClassifier(oob_score=True, n_estimators=500, n_jobs=8, max_depth=25)
        self.label_encoder = LabelEncoder()
        self._is_complete = False
        self.n_chunks_for_predict = 25
        self._EPOCH_LENGTH = 30 # seconds

    def train_from_polygraph_file_list(self, data_file_pattern):
        files = glob.glob(data_file_pattern)
        for f in sorted(files):
            print "Training with: " + f

            polygraph = polygraph_from_pkl(f)

            self.train(polygraph)

        self.complete_training()


    def make_features(self,epoch):
        groups_dfs = [group(epoch) for group in self.__FEATURE_GROUPS]
        features = pd.concat(groups_dfs, axis=1)
        return  features

    def save(self, filename):
        if not self._is_complete:
            raise Exception("Trainning not complete. Call `complete_training()' method")
        pkl.dump(self, filename, compress=9)

    def train(self, polygraph):

        for lab, channel in zip(polygraph.channel_types, polygraph.channels()):
            for _, epoch in channel.embed_seq(self._EPOCH_LENGTH,7): # fixme magic number
                feature_row = self.make_features(epoch)
                feature_row.index = [lab]

                self._features = pd.concat([self._features, feature_row])

        #todo append = cleanner

    def complete_training(self):
        x = self._features
        y = self.label_encoder.fit_transform(self._features.index)
        # drop features to save memory ;)
        self._features = None
        self.classifier.fit(x,y)
        self._is_complete = True

    def predict (self,  polygraph):
        idxs = []
        features = []
        for channel_idx, channel in enumerate(polygraph.channels()):
            for _, epoch in channel.embed_seq(self._EPOCH_LENGTH,10): # fixme magic number
                features.append(self.make_features(epoch))
                idxs.append(channel_idx)

        features = pd.concat(features)
        features.index = idxs

        y_pred = self.classifier.predict(features)
        y_pred = self.label_encoder.inverse_transform(y_pred)

        df = pd.DataFrame({"channel":idxs, "label":y_pred})

        a = pd.pivot_table(df, index="channel", columns="label", aggfunc=len).fillna(0)
        labels_pr = a.apply(lambda x: pd.Series({"pr": x[np.argmax(x)]/ float(np.sum(x)), "label":np.argmax(x)}), 1)
        labels = labels_pr.apply(lambda x : x["label"] if x["pr"] > 0.75 else "NaN" ,1)


        return list(labels), list(labels_pr["pr"])
