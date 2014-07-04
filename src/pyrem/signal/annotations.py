__author__ = 'quentin'


class Annotation(object):

    """
    Each annotation of a multiplexed signal is a data structure with a start, an end, a confidence level, a "type" and optionally, a value.
    Two annotations with the same type will not be allowed to overlap.
    """

    def __init__(self, start, end, confidence, type, value=True):
        """
        :param start: starting time in seconds
        :param end: ending time in second
        :param confidence: the confidence level (from 0 to 1)
        :param type: The name of the type of annotation (e.g. "spindle")
        :param value: The value for this annotation
        """
        self._start, self._end, self._confidence, self._type, self._value = start, end, confidence, type, value

    def to_dict(self):
        out = {"start":self.start,
               "end":self.end,
               "confidence":self.confidence,
               "type":self.type,
               "value":self.value}
        return out

    @property
    def start (self):
        return self._start
    @property
    def end (self):
        return self._end
    @property
    def confidence(self):
        return self._confidence
    @property
    def type(self):
        return self._type
    @property
    def value(self):
        return self._value
