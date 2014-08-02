__author__ = 'quentin'



import datetime
import numpy as np
import joblib as pkl
from scipy import signal
import pandas as pd
from pyrem.signal.utils import str_to_time
from scikits.samplerate import resample

SIGNALY_DPI = 328
SIGNAL_FIGSIZE = (30, 5)
MAX_POINTS_AMPLITUDE_PLOT = 1000




def signal_from_csv(file_name, sampling_freq):
    data = pd.read_csv(file_name, engine="c", header=None, dtype=np.float32)
    return BiologicalTimeSeries(data, sampling_freq)


def _normalise(mat):
    out = (mat - np.mean(mat,0)) / np.std(mat,0)
    return out


def signal_from_pkl(filename):
    return pkl.load(filename)


class BiologicalTimeSeries(np.ndarray):
    def __new__(cls, data, fs, type=None, name=None, metadata=None):

        obj = np.array(data).view(cls)
        # if len(obj.shape) != 1:
        #     raise ValueError("A Signal object can only be build from a 1D array")

        # add the new attribute to the created instance
        obj.__type = type
        obj.__fs = float(fs)
        obj.__metadata = metadata
        obj.__name = name

        # Finally, we must return the newly created object:
        return obj


    def __reduce__(self):
        state = list(np.ndarray.__reduce__(self))
        new_state = list(state[-1])
        new_state.append(self.__type)
        new_state.append(self.__fs)
        new_state.append(self.__metadata)
        new_state.append(self.__name)
        state[-1] = tuple(new_state)

        return tuple(state)

    def __setstate__(self, state):
        list_state = list(state)
        self.__name = list_state.pop()
        self.__metadata = list_state.pop()
        self.__fs = list_state.pop()
        self.__type = list_state.pop()
        return np.ndarray.__setstate__(self,tuple(list_state))
    #
    def save(self, filename, compression_level=5):
        pkl.dump(self,filename,compression_level)

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return

        self.__type = getattr(obj, 'type', None)
        self.__fs = getattr(obj, 'fs', None)
        self.__metadata = getattr(obj, 'metadata', None)
        self.__name = getattr(obj, 'name', None)

    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)

    @property
    def fs(self):
        return self.__fs

    @property
    def type(self):
        return self.__type


    @property
    def metadata(self):
        return self.__metadata

    def rename(self, name):
        self.__name =  name

    @property
    def name(self):
        return self.__name

    @property
    def duration(self):
        return self._time_from_idx(self.size)
#

    def _idx_from_time(self, time):
        return  int(time.total_seconds() * self.fs)


    def _time_from_idx(self, idx):
        start = datetime.datetime.fromtimestamp(0)
        end = datetime.datetime.fromtimestamp(float(idx) / self.fs)
        return  end - start
#
    def copy(self):
        return self._copy_attrs_to_array(self)


    def _copy_attrs_to_array(self, a, **kwargs):
        dic = {"type":self.type,
        "fs":self.fs,
        "metadata":self.metadata,
        "name":self.name}
        dic =  dict(dic.items() + kwargs.items())
        return self.__new__(type(self),a, **dic)

    def resample(self, new_fs):
        raise NotImplementedError

    def __repr__(self):
        if self.metadata:

            metadata = "\n".join(["\t\t%s:\t%s" % (k, str(v)) for k,v in self.metadata.items()])
        else:
            metadata = "\t\tNone"

        out = ["\n" + type(self).__name__ + "\n",
               "Name:\t%s" % (self.name),
               "Duration:\t%s (HH:mm:ss)" % (str(self.duration)),
               "Sampling freq:\t%f Hz" % (self.fs),
               "Type:\t%s" % (self.type),
               "N points:\t%i" % (self.size),
               "Metadata:\n%s" % (metadata),
               ]

        return "\n".join(out)

    def _get_time_slice( self, key ):

        if not key.start:
            start_idx = 0
        elif isinstance(key.start, str):
            start = str_to_time(key.start)
            start_idx = self._idx_from_time(start)
        elif isinstance(key.start, datetime.timedelta):
            start_idx = self._idx_from_time(key.start)

        else:
            raise NotImplementedError()

        if not key.stop:
            stop_idx = self.size
        elif isinstance(key.stop, str):
            stop = str_to_time(key.stop)

            stop_idx= self._idx_from_time(stop)

        elif isinstance(key.stop, datetime.timedelta):
            stop_idx = self._idx_from_time(key.stop)
        else:
            raise NotImplementedError()

        if start_idx > stop_idx:
            raise Exception("The starting time(%s), MUST be before the end time(%s)" % (str(start),str(stop)))

        # important detail: we never want an empty subchannel.
        # fixme it would be clever to implement this in the integer indexing engine as well...?
        if start_idx == stop_idx:
            stop_idx +=1
        return self[start_idx: stop_idx]

    def __getitem__( self, key ) :

        if isinstance( key, slice ):

            return self._get_time_slice(key)
        else:
            return np.ndarray.__getitem__(self,key)

    def iter_window(self, length, lag):
        """
        Iterate through an array by successive overlapping slices.
        Also returns the center of the slice

        :param lag: the ratio of overlap (1= no overlap, 0= completely overlapped)
        :param length:of the epoch (in second)
        :return: (centre_of_window, sub_signal)
        """
        if lag<=0:
            raise Exception("lag has to be  greater than one")

        n_points = int(self.fs * length)


        lag_in_points = int(n_points * lag)

        for i in np.arange(0, self.size - n_points, lag_in_points):
            out = self[i:i+n_points]
            #onset = out._time_from_idx(i)
            centre = ( i + float(out.size)/2.0) / self.fs
            if out.size < n_points:
                return
            yield centre , out


class Signal(BiologicalTimeSeries):
    #dummry init for pycharm completion
    def __init__(self,data, fs, **kwargs):
        pass
    def __new__(cls,data, fs, **kwargs):
        try:
            obj = np.array(data).view(cls).astype(np.float32)
        except:
            raise ValueError("Data could not be understood as an array of float32")

        if len(obj.shape) != 1:
            raise ValueError("A Signal object can only be build from a 1D array")

        return BiologicalTimeSeries.__new__(cls, data, fs, **kwargs)

    def resample(self, target_fs, mode="sinc_best"):
        # num = target_fs * self.size / self.fs
        ratio = target_fs / self.fs

        out = resample(self,ratio,mode)
        real_fs = out.size * self.fs / float(self.size)
        out = self._copy_attrs_to_array(out, fs=real_fs)
        return out


class Annotation(BiologicalTimeSeries):
    # dummy for doc and pycharm
    def __init__(self,data, fs, observation_probabilities=None, **kwargs):
        pass
    def __new__(cls,data, fs, observation_probabilities=None, **kwargs):
        try:
            chars = np.array(data,dtype=np.uint8)
        except:
            raise ValueError("could not understand input as an array of uint8")

        if len(chars.shape) != 1:
            raise ValueError("An Annotation object can only be build from a 1D array")

        if observation_probabilities is None:
            observation_probabilities = np.array([1.0] * data.size, dtype=np.float32)

        data = np.recarray((chars.size,),
                          dtype=[("values", np.uint8), ("probas", np.float32)])
        data["values"] = chars
        data["probas"] = observation_probabilities

        return BiologicalTimeSeries.__new__(cls, data, fs, **kwargs)

    def resample(self, new_fs):
        raise NotImplementedError #fixme

    @property
    def values(self):
        return self["values"]

    @property
    def probas(self):
        return self["probas"]
