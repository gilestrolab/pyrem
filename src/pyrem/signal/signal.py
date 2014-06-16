"""

"""
__author__ = 'quentin'

import numpy as np

class Signal(np.ndarray):
    def __init__(self, seq, sampling_freq):

        #super(Signal, self).__init__(seq)

        self.__sampling_freq = sampling_freq

    @property
    def sampling_freq(self):
        return self.__sampling_freq