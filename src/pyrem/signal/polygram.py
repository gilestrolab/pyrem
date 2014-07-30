__author__ = 'quentin'

import numpy as np
class Polygram(object):

    def __init__(self, channels, metadata=None):
        self._channels = channels
        self._metadata= metadata
        for i in self.channels:
            print type(i)

        durations = [c.duration for c in self.channels]

        max_duration =  max(durations)
        fs_max_duration =  channels[np.argmax(durations)].fs

        for c in channels:

            print self._test_duration( max_duration, fs_max_duration, c.duration, c.fs)


    def _test_duration(self, max_duration, fs_max_duration, duration, fs):
        delta = max_duration - duration
        smallest_fs = min([fs_max_duration, fs])
        longest_period = 1.0 / smallest_fs

        if delta.total_seconds() >= longest_period:
            return False

        return True


        #epsilon =






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
