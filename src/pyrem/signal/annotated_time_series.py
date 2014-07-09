__author__ = 'quentin'

import numpy as np
import datetime

class AnnotatedTimeSeries(object):
    def __init__(self,
                 data,
                 sampling_rate,
                 annotations=None,
                 channel_types=None,
                 metadata=None
                ):
        self.data = np.asarray(data)

        if annotations is None:
            self.annotations = None
        else:
            #todo check size should be equal to the number of rows
            self.annotations = np.asarray(annotations)
        if channel_types is None:
            channel_types = ["NaN"] * self.nchannels
        else:
            if len(channel_types) != self.nchannels:
                raise Exception("the number of channels does not match the number of elements in channel types")

        if not metadata:
            self._metadata = dict()
        else:
            self._metadata = metadata

        self._sampling_freq = sampling_rate
        self._channel_types= channel_types

    @property
    def sampling_freq(self):
        return self._sampling_freq

    @property
    def channel_types(self):
        return self._channel_types

    @property
    def ntimepoints(self):
        return self.data.shape[0]
    @property
    def metadata(self):
        return self._metadata

    def set_data(self, new_data):
        self.data = new_data

    @property
    def nchannels(self):
        return self.data.shape[1]


    def _soft_copy(self,new_data,new_annotations=None, new_channel_types=None):
        if new_annotations is None:
            annotations = self.annotations
        else:
            annotations = new_annotations

        if new_channel_types is None:
            new_channel_types = self.channel_types
        else:
            new_channel_types = new_channel_types

        return AnnotatedTimeSeries(new_data, self.sampling_freq, annotations = annotations, channel_types=new_channel_types, metadata = self.metadata)


    def channels(self):
        for i in range(self.nchannels):
            yield self[i]

    def __getitem__( self, key ) :

        # slice -> time chunk
        if isinstance( key, slice ) :

            sub_data = self.data[key]
            if self.annotations is None:
                sub_annotations =None
            else:
                sub_annotations = self.annotations[key]


            return self._soft_copy(sub_data, sub_annotations)

        else:
            # todo channel name indexation
            # if isinstance( key, str) :
            #     self.channel_names

            if isinstance( key, int ) :

                if key >= self.data.shape[1] or key < 0 :
                    raise IndexError, "The index (%d) is out of range."%key


                sub_data = self.data[:,key].reshape((self.ntimepoints,1))
                new_channel_types = [self.channel_types[key]]
                return self._soft_copy(sub_data, new_channel_types = new_channel_types)

            else:
                raise TypeError, "Invalid argument type."


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


    @property
    def duration(self):
        return self._time_from_idx(float(self.ntimepoints))


    def _time_from_idx(self, idx):
        start = datetime.datetime.fromtimestamp(0)
        end = datetime.datetime.fromtimestamp(idx / self.sampling_freq)
        return  end - start


    def __repr__(self):
        metadata = "\n".join(["\t\t%s:\t%s" % (k, str(v)) for k,v in self.metadata.items()])

        out = ["\n" + type(self).__name__ + "\n",
               "N channels:\t%i" % (self.nchannels),
               "duration:\t%s (HH:mm:ss)" % (str(self.duration)),
               "sampling freq:\t%f Hz" % (self.sampling_freq),
               "Channel types:\t%s" % (str(self.channel_types)),

               "metadata:\n%s" % (metadata),
               ]
        return "\n".join(out)




a = AnnotatedTimeSeries(np.reshape(np.arange(0,1200),(300,4)),10, channel_types=["EEG","EMG","NaN","EEG"],metadata={"a":1,"b":"yoyo"})
#print a.data
b = a[0:2]
#print b.data

c = a[1]
print c