"""

"""
import datetime

__author__ = 'quentin'

import numpy as np

SIGNALY_DPI = 328
SIGNAL_FIGSIZE = (30, 5)

def _normalise(obj):
    out = (obj - np.mean(obj)) / np.std(obj)
    return out


class Signal(np.ndarray):
    def __new__(cls, input, sampling_freq, normalised=True):

        obj = np.asarray(input).view(cls)
        if normalised:"
            obj = _normalise(obj)
        # add the new attribute to the created instance
        obj.__normalised = normalised
        obj.__sampling_freq = float(sampling_freq)

        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return

        self.__normalised = getattr(obj, 'normalised', None)
        self.__sampling_freq = getattr(obj, 'sampling_freq', None)



    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)



    @property
    def sampling_freq(self):
        return self.__sampling_freq

    @property
    def normalised(self):
        return self.__normalised

    @property
    def duration(self):
        return self._time_from_idx(float(self.size))


    def _time_from_idx(self, idx):
        start = datetime.datetime.fromtimestamp(0)
        end = datetime.datetime.fromtimestamp(idx / self.sampling_freq)
        return  str(end - start)


    def resample(self, new_sampling_freq):
        # new_size = self.size * float(new_sampling_freq) /self.sampling_freq
        new_step = self.sampling_freq / float(new_sampling_freq)
        new_t = np.arange(0, self.size, new_step)
        new_t = new_t[new_t <= self.size -1]
        old_t= np.arange(0, self.size)
        return Signal(np.interp(new_t, old_t, self), new_sampling_freq)

    def embed_seq(self, length, lag):
        """
        Iterate through an array by successive (typically overlapping) slices.
        Also returns the center of the slice

        :param lag:
        :param length:
        :return:
        """

        if lag<1:
            raise Exception("lag has to be at least 1")


        margin = (lag-1)/2

        for i in np.arange(0, self.size - length, lag):
            out = self[i:i+length]
            centre = i + float(out.size)/2.0
            if out.size < length:
                return
            yield centre, out


    def _create_plot(self, *args, **kwargs):
        from matplotlib import pyplot as plt
        title = "Duration = %s; at = %fHz" % (self.duration , self.sampling_freq)

        out = plt.plot(self, *args, **kwargs)
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
        fig = plt.figure(figsize=SIGNAL_FIGSIZE, dpi=SIGNALY_DPI)
        ax = self._create_plot()
        data = print_figure(fig, 'png')
        plt.close(fig)
        return data

    @property
    def png(self):
        from IPython.display import Image
        return Image(self._repr_png_(), embed=True)


