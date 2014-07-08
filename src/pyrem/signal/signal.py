__author__ = 'quentin'

import numpy as np
from datetime import timedelta
import pandas as pd

import datetime

__author__ = 'quentin'

import numpy as np
import cPickle
from scipy.interpolate import interp1d
import pandas as pd
# for plotting signals in ipn:
SIGNALY_DPI = 328
SIGNAL_FIGSIZE = (30, 5)


def signal_from_csv(file_name, sampling_freq):
    data = pd.read_csv(file_name, engine="c", header=None, dtype=np.float32)
    return Signal(data, sampling_freq)


def _normalise(mat):
    out = (mat - np.mean(mat,0)) / np.std(mat,0)
    return out


def load_signal(filename):
    return cPickle.load(open(filename))


class Signal(np.recarray):
    def __new__(cls, input, sampling_freq, normalise=True, channel_types=[],metadata=None):
        if not isinstance(input,np.recarray):
            data = np.asarray(input)
            # force array to be 2d

            if data.ndim != 2:
                data = np.reshape(data, (1,-1))

            if normalise:
                data = _normalise(data)
            types = [("c%i" %(i),np.float) for i in range(data.shape[1])]
            obj = np.recarray((data.shape[0],),types).view(cls)
            for i,(c,t) in enumerate(types):
                obj[c] = data[:,i]
        else:
            obj = np.copy(input).view(cls)
        # add the new attribute to the created instance
        obj.__normalised = normalise
        obj.__sampling_freq = float(sampling_freq)
        obj.__metadata = metadata
        obj.__channel_types = channel_types

        # Finally, we must return the newly created object:
        return obj

    def __deepcopy__(self):
        return Signal(self, self.sampling_freq, self.normalised)

    def __reduce__(self):
        state = list(np.ndarray.__reduce__(self))
        new_state = list(state[-1])
        new_state.append(self.__normalised)
        new_state.append(self.__sampling_freq)
        new_state.append(self.__metadata)
        new_state.append(self.__channel_types)
        state[-1] = tuple(new_state)

        return tuple(state)

    def __setstate__(self, state):
        list_state = list(state)
        self.__channel_types = list_state.pop()
        self.__metadata = list_state.pop()
        self.__sampling_freq = list_state.pop()
        self.__normalised = list_state.pop()

        return np.ndarray.__setstate__(self,tuple(list_state))

    def save(self, filename):
        with open(filename, "w") as f:
            cPickle.dump(self,f, cPickle.HIGHEST_PROTOCOL)

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return

        self.__normalised = getattr(obj, 'normalised', None)
        self.__sampling_freq = getattr(obj, 'sampling_freq', None)
        self.__metadata = getattr(obj, 'metadata', None)
        self.__channel_types = getattr(obj, 'channel_types', None)

    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)

    @property
    def nchannels(self):
        return len(self.dtype)

    @property
    def nsignals(self):
        return len([True for n in self.dtype.names if n.startswith('c')])

    @property
    def metadata(self):
        return self.__metadata
    @property
    def channel_types(self):
        return self.__channel_types

    @property
    def sampling_freq(self):
        return self.__sampling_freq

    @property
    def normalised(self):
        return self.__normalised

    @property
    def ntimepoints(self):
        return self.shape[0]

    @property
    def duration(self):
        return self._time_from_idx(float(self.ntimepoints))


    def _time_from_idx(self, idx):
        start = datetime.datetime.fromtimestamp(0)
        end = datetime.datetime.fromtimestamp(idx / self.sampling_freq)
        return  str(end - start)


    def signal_iter(self, return_names=False):
         for c_name in self.dtype.names:
             if c_name.startswith('c'):
                 if return_names:
                    yield  c_name, self[c_name]
                 else:
                     yield  self[c_name]

    def annot_iter(self, return_names=False):
        for c_name in self.dtype.names:
            if c_name.startswith('a'):
                if return_names:
                    yield  c_name, self[c_name]
                else:
                    yield  self[c_name]

    def channel_iter(self):
        for s in self.signal_iter():
            yield s

        for a in self.annot_iter():
            yield a

    def resample(self, new_sampling_freq):

        new_step = self.sampling_freq / float(new_sampling_freq)
        new_t = np.arange(0, self.ntimepoints, new_step)
        new_t = new_t[new_t <= self.ntimepoints -1]
        old_t= np.arange(0, self.ntimepoints)

        out = np.recarray((new_t.size,),self.dtype)

        for name, channel in self.channel_iter():

            if channel.dtype==float:
                kind="linear"
            else:
                kind="nearest"

            f = interp1d(old_t, channel, assume_sorted=True, kind=kind)
            out[name] = f(new_t)
        return Signal(out, new_sampling_freq)


    def embed_seq(self, length, lag):
        """
        Iterate through an array by successive overlapping slices.
        Also returns the center of the slice

        :param lag: the ratio of overlap (1= no overlap, 0= completely overlapped)
        :param length:of the epoch (in second)
        :return: a signal
        """

        if lag<=0 or lag>1:
            raise Exception("lag has to be between 0 and 1")

        n_points = int(self.sampling_freq * length)

        lag_in_points = int(n_points * lag)

        for i in np.arange(0, self.ntimepoints - n_points, lag_in_points):
            out = self[i:i+n_points]
            centre = ( i + float(out.ntimepoints)/2.0) / self.sampling_freq
            if out.ntimepoints < n_points:
                return
            yield centre, out
#
#
    def _create_fig(self, *args, **kwargs):
        from matplotlib import pyplot as plt

        title = "Duration = %s; at = %fHz" % (self.duration , self.sampling_freq)
        
        f, axarr = plt.subplots(self.nsignals , sharex=True, sharey=True)

        axarr[0].set_title(title)
        for i,(name,s) in enumerate(self.signal_iter()):
            axarr[i].plot(s, *args, **kwargs)
        out = plt

        location, _ = plt.xticks()
        plt.xticks(location, [self._time_from_idx(l) for l in location], rotation=45)

        return out

    def plot(self, *args, **kwargs):
        """
        Plots the signal using :mod:`matplotlib.pyplot`.

        :param args: arguments to pass to :func:`~matplotlib.pyplot.plot`
        :param kwargs: keyword arguments to pass to :func:`~matplotlib.pyplot.plot`
        :return: the result of the :func:`matplotlib.pyplot.plot` function.
        """
        return self._create_plot(*args, **kwargs)

    def _repr_png_(self):
        from IPython.core.pylabtools import print_figure
        from matplotlib import pyplot as plt

        fig = self._create_fig()
        data = print_figure(fig, 'png')
        plt.close(fig)
        return data

    @property
    def png(self):
        from IPython.display import Image
        return Image(self._repr_png_(), embed=True)


