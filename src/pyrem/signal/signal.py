"""

"""
from matplotlib.pyplot import title

__author__ = 'quentin'

import numpy as np

SIGNALY_DPI = 328
SIGNAL_FIGSIZE = (30, 5)

class Signal(np.ndarray):
    def __new__(cls, input, sampling_freq):

        obj = np.asarray(input).view(cls)
        # add the new attribute to the created instance
        obj.__sampling_freq = float(sampling_freq)
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.__sampling_freq = getattr(obj, 'sampling_freq', None)

    @property
    def sampling_freq(self):
        return float(self.__sampling_freq)


    def _time_from_idx(self, idx, format="seconds"):
        multipliers = {"seconds":1.0, "minutes":60.0, "hours":60*60.0}
        return  idx * self.sampling_freq * multipliers[format]


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

        length = self.size
        time_points = self._time_from_idx(np.arange(0,length,1))

        title = "Length = %is; sampling frequency = %f" % (self._time_from_idx(length), self.sampling_freq)
        return plt.plot( time_points, self, *args, **kwargs)

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