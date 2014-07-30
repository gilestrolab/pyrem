__author__ = 'quentin'



import datetime

__author__ = 'quentin'

import numpy as np
import joblib as pkl
from scipy import signal

# from scipy.interpolate import interp1d
import pandas as pd
# for plotting signals in ipn:
SIGNALY_DPI = 328
SIGNAL_FIGSIZE = (30, 5)
MAX_POINTS_AMPLITUDE_PLOT = 1000

import re
from datetime import timedelta

REGEX = re.compile(r'((?P<hours>\d+?)h)?((?P<minutes>\d+?)m)?((?P<seconds>\d+?)s)?((?P<milliseconds>\d+?)w)?')



def str_to_time(str):
    parts = REGEX.match(str)

    if not parts or len(parts.group()) <1:
        raise ValueError("The string `%s' does not contain time information.\nThe syntax is, for instance, 1h99m2s" %(str))
    parts = parts.groupdict()
    time_params = {}
    for (name, param) in parts.iteritems():
        if param:
            time_params[name] = int(param)
    return timedelta(**time_params)





def signal_from_csv(file_name, sampling_freq):
    data = pd.read_csv(file_name, engine="c", header=None, dtype=np.float32)
    return Signal(data, sampling_freq)


def _normalise(mat):
    out = (mat - np.mean(mat,0)) / np.std(mat,0)
    return out


def signal_from_pkl(filename):
    return pkl.load(filename)


class Signal(np.ndarray):
    def __new__(cls, input, fs, type=None, name=None, metadata=None):

        obj = np.array(input,dtype=np.float32).view(Signal)
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


    @property
    def name(self):
        return self.__name

    @property
    def duration(self):
        return self._time_from_idx(self.size)
#

    def _idx_from_time(self, time):
        return  int(time.microseconds * 1e-6 * self.fs)


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
        return Signal(a, **dic)

    def resample(self, new_fs, type="fourier"):
        if type == "fourier":
            num = new_fs * self.size / self.fs
            out = signal.resample(self, int(num))
            out = self._copy_attrs_to_array(out, fs=new_fs)
            return out
        else:
            raise NotImplementedError

    def __repr__(self):
        metadata = "\n".join(["\t\t%s:\t%s" % (k, str(v)) for k,v in self.metadata.items()])

        out = ["\n" + type(self).__name__ + "\n",
               "Name:\t%s" % (self.name),
               "Duration:\t%s (HH:mm:ss)" % (str(self.duration)),
               "Sampling freq:\t%f Hz" % (self.fs),
               "Type:\t%s" % (self.type),
               "N points:\t%i" % (self.size),
               "metadata:\n%s" % (metadata),
               ]

        return "\n".join(out)
        #
        # new_step = self.sampling_freq / float(new_sampling_freq)
        # new_t = np.arange(0, self.ntimepoints, new_step)
        # new_t = new_t[new_t <= self.ntimepoints -1]
        # old_t= np.arange(0, self.ntimepoints)
        #
        # out = np.recarray((new_t.size,),self.dtype)
        #
        # for name, channel in self.channel_iter():
        #
        #     if channel.dtype==float:
        #         kind="linear"
        #     else:
        #         kind="nearest"
        #
        #     f = interp1d(old_t, channel, assume_sorted=True, kind=kind)
        #     out[name] = f(new_t)
        # return Signal(out, new_sampling_freq)


    # def __getslice__(self, start, stop) :
    #
    #     return self.__getitem__(slice(start, stop))

    def _get_time_slice( self, key ):

        if not key.start:
            start_idx = 0
        elif isinstance(key.start, str):
            start = str_to_time(key.start)
            print start
            start_idx = self._idx_from_time(start)

        else:
            raise NotImplementedError()

        if not key.stop:
            stop_idx = self.size
        elif isinstance(key.stop, str):
            stop = str_to_time(key.stop)

            stop_idx= self._idx_from_time(stop)

        else:
            raise NotImplementedError()


        if start_idx > stop_idx:
            raise Exception("The starting time(%s), MUST be before the end time(%s)" % (str(start),str(stop)))


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
    def _represent_long_time_series(self,max_points=MAX_POINTS_AMPLITUDE_PLOT):
        from matplotlib import pyplot as plt
        winsize_npoints = float(self.size) / float(max_points)
        secs = winsize_npoints / self.fs
        print secs
        mins, maxes, means, sds,xs = [],[],[],[],[]
        for c, w in  self.iter_window(secs,1):
            mins.append(np.min(w))
            maxes.append(np.max(w))
            means.append(np.mean(w))
            sds.append(np.std(w))
            xs.append(c)

        means = np.array(means)
        mean_plus_sd = means +sds
        mean_minus_sd = means - sds
        plt.figure()


        plt.fill_between(xs,mins,maxes, facecolor=(0,0,1,0.6),edgecolor=(0,0,0,0.2), antialiased=True)
        plt.fill_between(xs,mean_minus_sd, mean_plus_sd, facecolor=(1,0.5,0,0.9),edgecolor=(0,0,0,0), antialiased=True)
        plt.plot(xs,means,"-", linewidth=1, color='k')
        plt.show()



#
# #
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

#     def plot(self, *args, **kwargs):
#         """
#         Plots the signal using :mod:`matplotlib.pyplot`.
#
#         :param args: arguments to pass to :func:`~matplotlib.pyplot.plot`
#         :param kwargs: keyword arguments to pass to :func:`~matplotlib.pyplot.plot`
#         :return: the result of the :func:`matplotlib.pyplot.plot` function.
#         """
#         return self._create_plot(*args, **kwargs)
#
#     def _repr_png_(self):
#         from IPython.core.pylabtools import print_figure
#         from matplotlib import pyplot as plt
#
#         fig = self._create_fig()
#         data = print_figure(fig, 'png')
#         plt.close(fig)
#         return data
#
#     @property
#     def png(self):
#         from IPython.display import Image
#         return Image(self._repr_png_(), embed=True)
#
#

