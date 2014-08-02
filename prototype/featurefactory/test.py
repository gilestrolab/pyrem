from pyrem.signal.io import polygram_from_pkl
from pyrem.features.univariate import *
from pyrem.signal.wavelet_decomposition import decompose_signal
from pyrem.signal.signal import Signal,Annotation
import pandas as pd
from scipy import stats
from datetime import timedelta

class FeatureFamilyBase(object):
    r"""
    A feature family object is a process returning a vector of features upon analysis of some data.
    Features are returned as a pandas DataFrame object, with column names for features. Each feature name is prefixed by the name of the
    feature family. This is an abstract class designed to be derived by:

    1. Defining a ``prefix`` attribute. It will add the name of the family to the name of the features.
    2. Overriding the ``_make_feature_vec`` method. It should return a dictionary of scalars, each being a feature.

    """
    prefix = None
    def make_vector(self, signal):
        """
        Compute one vector of features from polygraph.

        :param data: A signal
        :type data: :class:`~pyrem.signal.polygraph.Polygraph`
        :return: a one-row dataframe
        :rtype: :class:`~pandas.DataFrame`
        """
        if not self._check_channel_type(signal):
            return

        feature_dict = self._make_feature_vec(signal)

        data_frame = pd.DataFrame(feature_dict, index=[None])

        if len(feature_dict) > 1 and self.prefix is None:
            raise Exception("More than one features in this group. You need a prefix to identify this group")

        if self.prefix:
            data_frame.columns = [signal.name +"."+self.prefix + "." + c for c in data_frame.columns]
        return data_frame

    def _check_channel_type(self,data):
        return NotImplementedError

    def _make_feature_vec(self,data):
        raise NotImplementedError

class SignalFeatureBase(FeatureFamilyBase):
    def _check_channel_type(self,channel):
        return isinstance(channel, Signal)

class AnnotationFeatureBase(FeatureFamilyBase):
    def _check_channel_type(self,channel):
        return isinstance(channel, Annotation)

class VigilState(AnnotationFeatureBase):
    prefix = "vigil"
    def _make_feature_vec(self, channel):

        r = channel.values
        uniqs = np.unique(r)
        i = channel.probas
        probs = []
        for u in uniqs:
            eqs = (r == u)
            probs.append(np.sum(i[eqs]))

        probs = np.array(probs)
        probs /= np.sum(i)
        max_prob_idx =  np.argmax(probs)

        out = dict()
        out["value"] = uniqs[max_prob_idx]
        out["proba"] = probs[max_prob_idx]
        if out["value"] == 63 and  out["proba"] > 0:
            print "r", r
            print "i", i
            print "probs", probs
            raise Exception
        return out

class PowerFeatures(SignalFeatureBase):
    prefix = "power"
    def _make_feature_vec(self, channel):

        out = dict()

        powers = channel ** 2
        out["mean"] = np.mean(powers)
        out["sd"] = np.std(powers)
        out["median"] = np.median(powers)
        out["skew"] = stats.skew(powers)
        out["kurt"] = stats.kurtosis(powers)

        return out

#
class NonLinearFeatures(SignalFeatureBase):
    prefix = "nl"
    def _make_feature_vec(self,channel):
        out = dict()
        out["hurst"] = hurst(channel)
        #out["dfa"] = dfa(data)
        return out
#
#
class HjorthFeatures(SignalFeatureBase):
    prefix = "hjorth"
    def _make_feature_vec(self, channel):
        a,m,c = hjorth(channel)
        out = {"activity":a, "morbidity":m, "complexity":c}
        return out
#


pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/GFP_A.pkl")

eegs = decompose_signal(pol["EEG_parietal_cereb"], levels=[3,4,5])

emgs = decompose_signal(pol["EMG_1"],[1,2,3,4],keep_a=False)

pol2 = eegs.merge(emgs)
pol2 = pol2.merge(pol["vigilance_state"])
print np.sum(pol2["vigilance_state"].probas == 0)
print np.sum(pol["vigilance_state"].probas == 0)


pol2 = pol2.map_signal_channels(lambda c : (c - np.mean(c))/ np.std(c))

raise Exception

print np.mean(pol2[0])
#pol2.show()

feature_factory = [
                    # PowerFeatures(),
                    # HjorthFeatures(),
                    # NonLinearFeatures(),
                    VigilState(),]
all_rows = []
try:
    for t, w in pol2.iter_window(20,1):
        dfs = []
        for c in w.channels:
            for f in feature_factory:
                feature_vec = f.make_vector(c)
                if not feature_vec is None:
                    dfs.append(feature_vec)
        # print int(100 * t/ pol2.duration.total_seconds()), "%"
        print t

        print c.values
        print c.probas
        print c
        row = pd.concat(dfs, axis=1)
        row.index = [t]
        all_rows.append(row)
    df = pd.concat(all_rows)
except Exception as e:
    print e
    pol2.show()

print df
OUT_CSV = "/tmp/my_csv.out"

df.to_csv(OUT_CSV, float_format="%e")


# v = np.array(df["vigilance_state.vigil.value"])
# p = np.array(df["vigilance_state.vigil.proba"])
# v = np.array(df["EEG_parietal_cereb_cD_3.power.median"])
#
# sig = Signal(v, 1/20.0,  name="new_sig")
# pol = pol.merge(sig)
#
# pol.show()



    #print



