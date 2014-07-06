__author__ = 'quentin'


from pyrem.signal.signal import Signal
#from pyrem.signal.annotations import Annotation
import numpy as np
import cPickle
import json


def parse_signal(filename, sampling_rate=200, resampling=256):
    arr = np.genfromtxt(filename, delimiter=',')
    if arr.ndim ==1:
        signals = [Signal(arr,sampling_rate)]
    else:
        signals = [Signal(s,sampling_rate) for s in arr.T]

    signals = [s.resample(resampling) for s in signals]

    multiplex = Multiplex(signals)

    return multiplex

class Multiplex(list):

    def __init__(self,signals):
        super(Multiplex, self).__init__(signals)
        self._annotations = []

    def save(self, filename):
        with open(filename, "w") as f:
            cPickle.dump(self,f, cPickle.HIGHEST_PROTOCOL)
    @property
    def annotations(self):
        return self._annotations

    def annotations_as_json(self):
        return json.dumps(self._annotations)



    # def save(self, filename):
    #     np.savez_compressed(filename, *self)