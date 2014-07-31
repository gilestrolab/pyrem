from datetime import timedelta

__author__ = 'quentin'

import numpy as np
import joblib as pkl


class Polygram(object):

    def __init__(self, channels, metadata=None):
        self._channels = channels
        self._metadata= metadata

        durations = [c.duration for c in self.channels]

        max_duration =  max(durations)
        fs_max_duration =  channels[np.argmax(durations)].fs

        for c in self.channels:
            if not self._test_duration( max_duration, fs_max_duration, c.duration, c.fs):
                raise ValueError("Channels must have approximately the same length")

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
        new_channels = [c.copy() for c in self.channels]
        return Polygram(new_channels, self.metadata)

    def save(self, filename, compression_level=5):
        pkl.dump(self,filename,compression_level)


    def _get_time_slice( self, key ):
        if isinstance(key.start, str) and isinstance(key.start, str):
            channel_views = [c[key.start: key.stop] for c in self.channels]
            return Polygram(channel_views, metadata=self.metadata)
        else:
            raise NotImplementedError


    def __getitem__( self, key ) :
        if isinstance( key, slice ):
            return self._get_time_slice(key)
        elif isinstance( key, int):
            return self.channels[key]

        elif isinstance( key, str):
            if key not in self.channel_names:
                raise ValueError("`%s' is not a valid channel names.\nChannel names are:\n%s"
                                 % (key,self.channel_names))
            return self[self.channel_names.index(key)]
        else:
            print key
            raise NotImplementedError

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
        length_td = timedelta(seconds=length)
        lag_td = length_td * lag
        start = timedelta()

        while start < (min_duration - length_td):
            stop = start + length_td
            center = ((stop + start) /2).total_seconds()
            yield center, Polygram([c[start:stop] for c in self.channels], self.metadata)
            start += lag_td


    @property
    def channels(self):
        return self._channels

    @property
    def metadata(self):
        return self._metadata

    @property
    def channel_names(self):
        return [c.name for c in self.channels]

    @property
    def channel_types(self):
        return [c.type for c in self.channels]
