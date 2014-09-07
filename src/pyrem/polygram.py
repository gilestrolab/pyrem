r"""
===========================
Polygram
===========================


This module provides :class:`~pyrem.polygram.Polygram`; a container for biological time series
such as :class:`~pyrem.time_series.Signal` and :class:`~pyrem.time_series.Annotation`.
In this respect, it is inspired from pandas :class:`~pandas.TimeSeries` and :class:`~pandas.DataFrame`.
You can think about it as a dataframe where each column is a signal, or an annotation, and each row a time point.

The originality of :`~pyrem.polygram.Polygram` is to be able to deal with **heterogeneous (between signals) sampling rates**.
It contains time series with the same approximate duration, but different number of points.
This is typical when dealing with physiological time series because different variable will be
recorded at different sampling rates (see for instance, the [EDF]_ data format).
Another situation it which this could  be useful, is when performing a wavelet decomposition of a signal.
Indeed, one would obtain a set of time series (coefficients) of the same duration, but with different sampling rates (i.e. :math:`fs_{D_N} = 2fs_{D_{N+1}}`).

Systematically resampling signals, and annotations, to the maximal sampling rate is not  trivial, and would impact
significantly computational efficiency.


.. [EDF] B. Kemp and J. Olivan, "European data format 'plus' (EDF+), an EDF alike standard format for the exchange of physiological data,"
        Clinical Neurophysiology, vol. 114, no. 9, pp. 1755-1761, Sep. 2003.



----------------------------
Constructing a Polygram
----------------------------

First, let us create a couple of :class:`~pyrem.time_series.BiolgicalTimeSeries`:

>>> import numpy as np
>>> from pyrem.time_series import Annotation, Signal
>>> from pyrem.polygram import Polygram
>>>
>>> # create an Annotation with 1000 random values, sampled at 1.0Hz
>>> probs = np.random.random_sample(1000)
>>> vals = (np.random.random_sample(1000) * 4 +1).astype(np.int)
>>> annot = Annotation(vals,fs=1.0, observation_probabilities=probs, type="vigilance_state", name="state")
>>>
>>> # now a random walk signal of 100000 points at 100.0Hz
>>> rw = np.cumsum(np.random.normal(0,1,100000))
>>> sig = Signal(rw, fs=100.0,type="eeg", name="eeg1")
>>>
>>> # Once we have our time series, we can just do:
>>> pol = Polygram([annot, sig])
>>> #printing the object shows the characteristic of each channels
>>> pol
Polygram
-----
Duration:	0:16:40 (HH:mm:ss)
N signals:	1
N annotations:	1
Metadata:
		None
----
Channel information:
    Name             Type fs(Hz) Duration
0   eeg1              eeg  100.0  0:16:40
1  state  vigilance_state    1.0  0:16:40

.. note ::
    **Slightly different durations are allowed**

    The constructor will raise an error if the provided channels do not have the same duration:

    >>> Polygram([annot[:"11m"], sig[:"10m"]])
    ValueError
    'Channels must have approximately the same length.
    The durations of the input channels are:['0:10:00', '0:11:00']'

    However, in practice, it is almost impossible to obtain discrete signal of the exact same duration.
    Imagine, for instance that you have a first signal of 14 points at 3Hz (~ 4.667s), and a second signal of 5 points at 1Hz (5.0s).
    In this case, it is impossible to have exactly 14/3s of signal form a 1Hz signal.
    This could be represented by:

    >>> 0123456789abcd-   #  3Hz => one symbol/point
    >>> AAABBBCCCDDDEEE   # 1Hz => one LETTER/point
    >>> AAABBBCCCDDD---   # 1Hz => one LETTER/point

    Here, neither the second nor the third signal match, exactly, the duration of the first, but bot are approximately the same duration as the first.

    **A Polygram will tolerate this sort of mismatch if and only if all pairs of channels are within one period of the time series  with the channel longest period.**



------------------
Accessing channels
------------------
Often, you will want to extract a channel by name:

>>> pol.channel_names
['eeg1', 'state']
>>> pol['eeg1']
Signal
----
Name:	eeg1
Duration:	0:16:40 (HH:mm:ss)
Sampling freq:	100.000000 Hz
Type:	eeg
N points:	100000
Metadata:
		None

>>> # this is equivalent to
>>> pol[0]

You can also iterate through channels:

>>> [c.size for c in pol.channels]
[100000, 1000]

----------------------------
With strings and time-deltas
----------------------------

Because time series are potentially at different sampling rates, it makes no sense to index a polygram by range of integers:

>>> #does NOT work
>>> # pol[10:20]

Instead, time string and :class:`datetime.timedelta` can be used for extracting a sub_polygram:

>>> pol["1m":"2m"]

Indexing rules are similar to :mod:`~pyrem.time_series`

.. note ::
    **Indexing does NOT deep copy**

    When getting an epoch (temporal slice), of a polygram, the channel in the new polygram are *views* to the underlying data of the original channel.
    Like for numpy arrays, *modifying the data in a sub-polygram will modify the parent polygram*. To avoid this behaviour, one can call :func:`~pyrem.polygram.Polygram.copy`

----------------------------
Epoch iteration
----------------------------

If you want to extract features for each epoch and each channel, you may want to use th
:func:`~pyrem.polygram.Polygram.iter_window` iterator.
It works like the :func:`~pyrem.time_series.BiologicalTimeSeries.iter_window`

"""

from datetime import timedelta

__author__ = 'quentin'

import numpy as np
import joblib as pkl
import pandas as pd
from pyrem.time_series import Signal, Annotation
from pyrem.visualization import PolygramDisplay


class Polygram(object):

    def __init__(self, channels, metadata=None):
        """
        :param channels: a list of time series with approximately the same duration
        :type channels: list(:class:`~pyrem.time_series.BiolgicalTimeSeries`)
        :param metadata:    a dictionary of additional information (e.g. experimental variables)
        :type metadata: dict
        """


        annotations_channels = [c for c in channels if isinstance(c, Annotation)]
        signal_channels = [c for c in channels if isinstance(c, Signal)]

        annotations_channels = sorted(annotations_channels, key= lambda x: x.name)
        signal_channels = sorted(signal_channels, key= lambda x: x.name)

        self._channels = signal_channels + annotations_channels
        self._metadata= metadata

        durations = [c.duration for c in self.channels]

        max_duration =  max(durations)
        fs_max_duration =  channels[np.argmax(durations)].fs

        for c in self.channels:
            if not self._test_duration( max_duration, fs_max_duration, c.duration, c.fs):
                raise ValueError("Channels must have approximately the same length."
                                 "\nThe durations of the input channels are:\n%s", (str([str(c.duration) for c in self.channels ])))

        duplicate_names  = set([x for x in self.channel_names if self.channel_names.count(x) > 1])
        if len(duplicate_names) > 0:
            raise ValueError("Channels CANNOT have the same name. Duplicated names:\n %s"
                             % (str("\n".join(duplicate_names))))


    def _test_duration(self, max_duration, fs_max_duration, duration, fs):
        delta = max_duration - duration
        smallest_fs = min([fs_max_duration, fs])
        longest_period = 1.0 / smallest_fs

        if delta.total_seconds() >= longest_period:
            return False

        return True

    def copy(self):
        """
        Deep copy of an Polygram

        :return: a new Polygram with the same values
        :rtype: :class:`~pyrem.polygram.Polygram`
        """
        new_channels = [c.copy() for c in self.channels]
        return Polygram(new_channels, self.metadata)

    def save(self, filename, compression_level=5):
        """
        Efficiently save a Polygram using joblib

        :param filename: the output file name
        :type filename: str
        :param compression_level: an integer between 1 and 9. More is better, but slower. 5 is generally a good compromise
        :type compression_level: int
        """
        pkl.dump(self,filename,compression_level)


    def _get_time_slice( self, key ):
        if isinstance(key.start, str) and isinstance(key.start, str):
            channel_views = [c[key.start: key.stop] for c in self.channels]
            return Polygram(channel_views, metadata=self.metadata)
        else:
            raise NotImplementedError


    def __getitem__( self, key ) :

        if isinstance( key, list):
            key = set(key)
        if isinstance( key, set):
            return Polygram([self.__getitem__(k) for k in key], metadata=self.metadata)

        if isinstance( key, slice ):
            return self._get_time_slice(key)
        elif isinstance( key, int):
            return self._channels[key]

        elif isinstance( key, str):
            if key not in self.channel_names:
                raise ValueError("`%s' is not a valid channel names.\nChannel names are:\n%s"
                                 % (key,self.channel_names))
            return self[self.channel_names.index(key)]
        else:
            print key
            raise NotImplementedError
    def merge(self, obj, trim_channel=True):
        """
        Adds channels from a polygram to another polygram, or append a time series to a polygram

        :param obj: either a polygram or a time series to be added
        :type obj:  :class:`~pyrem.polygram.Polygram` or :class:`~pyrem.time_series.BiologicalTimeSeries`
        :param trim_channel: whether the new channel(s), if they have a longer duration, would be shortened to match the existing polygram.
        :type trim_channel: bool

        """

        if isinstance(obj, list):
            out = self.copy()
            for c in obj:
                out = out.merge(c)
            return out
        elif isinstance(obj, Polygram):
            print obj
            out = self.copy()
            print "!!!!", out.duration, self.duration
            for c in obj.channels:
                out = out.merge(c)
            return out

        channel  = obj
        if channel.duration <= self.duration or (not trim_channel):
            appended_channel = channel
        else:
            appended_channel = channel[:self.duration]

        new_channels = self._channels + [appended_channel]

        return Polygram(new_channels, self.metadata)


    def iter_window(self, length, lag):
    #     """
    #     Iterate through an array by successive overlapping slices.
    #     Also returns the center of the slice
    #
    #     :param lag: the ratio of overlap (1= no overlap, 0= completely overlapped)
    #     :param length:of the epoch (in second)
    #     :return: (centre_of_window, sub_signal)
    #     """
        if lag<=0:
            raise Exception("lag has to be  greater than one")

        min_duration = min([ c.duration for c in self.channels])
        length_td = timedelta(seconds=length  * lag)
        lag_td = length_td
        start = timedelta()

        while start < (min_duration - length_td):
            stop = start + length_td
            center = ((stop + start) /2).total_seconds()
            yield center, Polygram([c[start:stop] for c in self.channels], self.metadata)
            start += lag_td


    def __repr__(self):
        if self.metadata:
            metadata = "\n".join(["\t\t%s:\t%s" % (k, str(v)) for k,v in self.metadata.items()])
        else:
            metadata = "\t\tNone"


        out = ["\n" + type(self).__name__ + "\n",
               "Duration:\t%s (HH:mm:ss)" % (str(self.duration)),
               "N signals:\t%i" % (self.n_signals),
               "N annotations:\t%i" % (self.n_annotations),
               "Metadata:\n%s" % (metadata),
               ]

        headers = ["Name", "Type", "fs(Hz)", "Duration"]
        detail = []
        for c in self.channels:
            row = [c.name, c.type, c.fs, c.duration]
            row = [str(r) for r in row]
            detail.append(row)

        detail = str(pd.DataFrame(detail, columns=headers, index=[i for i in range(len(detail))]))

        out = "\n".join(out)
        out = out + "\n\nChannel information:\n" + detail
        return out

    def show(self):
        """
        Interactively displays a polygram using matplotlib. **Very urresponsive and prototypical at the minute**
        """
        PolygramDisplay(self).show()

    @property
    def channels(self):
        """
        An iterator through the all the channels
        """
        for c in self._channels:
            yield c

    @property
    def signal_channels(self):
        """
        An iterator through the all the *signal* channels
        """
        for c in self._channels:
            if isinstance(c,Signal):
                yield c

    @property
    def annotation_channels(self):
        """
        An iterator through the all the *annotation* channels
        """
        for c in self._channels:
            if isinstance(c,Annotation):
                yield c

    def map_signal_channels(self, fun):
        """
        Applies a function to all signal channels and returns a new Polygram with modified channels

        An example of how to normalise all signal channels

        >>> pol_norm =  pol.map_signal_channels(
        >>>         lambda x: (x - np.mean(x))/np.std(x))
        >>> np.mean(pol[0])
        >>> np.mean(pol_norm[0])


        :param fun: a function to be applied
        :type fun: callable
        :return: a new polygram
        :rtype:  :class:`~pyrem.polygram.Polygram`
        """

        if not callable(fun):
            raise ValueError("fun must be a function!")
        out = self.copy()
        for i, c in enumerate(out._channels):
            if isinstance(c,Signal):
                out._channels[i] = fun(c)
        return out


    @property
    def n_channels(self):
        """
        :return: The total number of channels
        :rtype: int
        """
        return len(self._channels)

    @property
    def n_signals(self):
        """
        :return: The total number of *signal* channels
        :rtype: int
        """
        return len([_ for _ in self.signal_channels])

    @property
    def n_annotations(self):
        """
        :return: The total number of *annotation* channels
        :rtype: int
        """
        return len([_ for _ in self.annotation_channels])

    @property
    def metadata(self):
        """
        :return: the metadata of this polygram
        :rtype: dict
        """
        return self._metadata

    @property
    def channel_names(self):
        """
        :return: The list of channel names
        :rtype: list(str)
        """
        return [c.name for c in self.channels]

    @property
    def channel_types(self):
        """
        :return: the types of all channels
        :rtype: list(str)
        """
        return [c.type for c in self.channels]
    @property
    def duration(self):
        """
        :return: The duration total of the polygram. That is the duration of the channel with the longest duration
        :rtype: :class:`datetime.timedelta`
        """
        return max([c.duration for c in self.channels])
