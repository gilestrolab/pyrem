__author__ = 'quentin'

import datetime
from sklearn.externals import joblib as pkl
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

def polygraph_from_pkl(filename):
    return pkl.load(filename)


def polygraph_from_csv(file_name, sampling_freq, comment= "#"):
    #todo
    data = pd.read_csv(file_name, engine="c", header=None, dtype=np.float32, comment=comment)
    return Polygraph(data, sampling_freq)


class Polygraph(object):

    def __init__(self,
                 data,
                 sampling_rate,
                 annotations=None,
                 channel_types=None,
                 annotation_types=None,
                 metadata=None,
                ):
        self.data = np.asarray(data)

        #todo force 2d array for annot and data
        #self.data.reshape((self.ntimepoints,1))

        if annotations is None:
            self._annotations = None
        else:
            self._annotations = np.asarray(annotations)
            if annotations.shape[0] != self.ntimepoints:

                raise Exception("The length of the provided annotations does not match the length of the data:\n %i != %i" %(
                                                annotations.shape[0],self.ntimepoints))

            if annotation_types is None:
                annotation_types= ["NaN"] * self.n_annotations

            if len(annotation_types) != annotations.shape[1]:
                print len(annotation_types), annotations
                raise Exception("the number of annotations does not match the number of elements in annotation types")


        if channel_types is None:
            channel_types = ["NaN"] * self.n_channels
        else:
            if len(channel_types) != self.n_channels:
                raise Exception("the number of channels does not match the number of elements in channel types")

        if not metadata:
            self._metadata = dict()
        else:
            self._metadata = metadata

        self._sampling_freq = sampling_rate
        self._channel_types= channel_types
        self._annotation_types= annotation_types


    def apply_channels(self, function, *args, **kwargs):
        new_data = []

        for c in self.channels():
            new_data.append(function(c.data.flatten(),*args, **kwargs))

        new_data = np.array(new_data).T
        return self._soft_copy(new_data)

    def normalise(self):
        return self.apply_channels( lambda x : (x - np.mean(x)) / np.std(x))

    @property
    def sampling_freq(self):
        return self._sampling_freq

    @property
    def channel_types(self):
        return self._channel_types
    @property
    def annotation_types(self):
        return self._annotation_types

    @property
    def ntimepoints(self):
        return self.data.shape[0]

    @property
    def metadata(self):
        return self._metadata

    def set_data(self, new_data):
        self.data = new_data

    @property
    def n_channels(self):
        return self.data.shape[1]

    @property
    def annotation_data(self):
        return self._annotations

    @property
    def n_annotations(self):
        return self._annotations.shape[1]




    def _soft_copy(self,new_data,new_annotations=None, new_channel_types=None):
        if new_annotations is None:
            annotations = self._annotations
        else:
            annotations = new_annotations

        if new_channel_types is None:
            new_channel_types = self.channel_types
        else:
            new_channel_types = new_channel_types

        return Polygraph(new_data, self.sampling_freq, annotations = annotations,
                                   channel_types = new_channel_types,
                                   annotation_types = self.annotation_types,
                                   metadata = self.metadata)



    def channels(self):
        for i in range(self.n_channels):
            yield self[i]

    def annotations(self):
        for i in range(self.n_annotations):
            yield self._annotations[:,i]


    def __getitem__( self, key ) :

        # slice -> time chunk
        if isinstance( key, slice ) :

            if isinstance(key.start, int):
                sub_data = self.data[key]
                if self._annotations is None:
                    sub_annotations = None
                else:
                    sub_annotations = self._annotations[key]
            #time string slices ;)

            else:
                raise NotImplementedError("TODO")
                pass


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

        :param lag: the ratio of overlap (e.g. 1= no overlap, 0= completely overlapped, 10= 9 * length between end of e and start of e+1)
        :param length:of the epoch (in second)
        :return: a signal
        """

        if lag<=0:
            raise Exception("lag has to be  greater than one")

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

    # def _idx_from_time(self, time):
    #     ref = datetime.datetime.fromtimestamp(0)
    #
    #     end = datetime.datetime.fromtimestamp(idx / self.sampling_freq)
    #     return  end - start


    def resample(self, new_sampling_freq):

        new_step = self.sampling_freq / float(new_sampling_freq)
        new_t = np.arange(0, self.ntimepoints, new_step)
        new_t = new_t[new_t <= self.ntimepoints -1]
        old_t= np.arange(0, self.ntimepoints)



        new_data = []
        for c in self.channels():
            f = interp1d(old_t, c.data.flatten(), assume_sorted=True, kind="linear")
            new_data.append(f(new_t))

        new_data = np.array(new_data).T

        if self._annotations is None:
            new_annotations = None
        else:
            new_annotations = []
            for c in self.annotations():
                f = interp1d(old_t, c.flatten(), assume_sorted=True, kind="nearest")
                new_annotations.append(f(new_t))

            new_annotations= np.array(new_annotations).T


        return Polygraph(new_data,  new_sampling_freq, channel_types=self.channel_types,
                                   annotations=new_annotations, metadata=self.metadata, annotation_types=self.annotation_types)


    def __repr__(self):
        metadata = "\n".join(["\t\t%s:\t%s" % (k, str(v)) for k,v in self.metadata.items()])

        out = ["\n" + type(self).__name__ + "\n",
               "N channels:\t%i" % (self.n_channels),
               "duration:\t%s (HH:mm:ss)" % (str(self.duration)),
               "sampling freq:\t%f Hz" % (self.sampling_freq),
               "Channel types:\t%s" % (str(self.channel_types)),
               "N points:\t%i" % (self.ntimepoints),
               "metadata:\n%s" % (metadata),
               ]

        return "\n".join(out)

    def summary(self):
        first = self[0:10]
        last = self[-10:]
        channels_f = pd.DataFrame(first.data, columns=self.channel_types)
        channels_l = pd.DataFrame(last.data, index=range(self.ntimepoints - 10, self.ntimepoints),
                                                columns=self.channel_types)

        annotation_f = pd.DataFrame(first.annotation_data, columns=self.annotation_types)
        annotation_l = pd.DataFrame(last.annotation_data, index=range(self.ntimepoints - 10, self.ntimepoints),
                                                columns=self.annotation_types)

        out = pd.concat([pd.concat([channels_f, channels_l]), pd.concat([annotation_f, annotation_l])], axis=1)

        return out


#
    def _create_fig(self, *args, **kwargs):
        from matplotlib import pyplot as plt

        title = "Duration = %s; at = %fHz" % (str(self.duration) , self.sampling_freq)

        f, axarr = plt.subplots(self.n_channels , sharex=True)

        axarr[0].set_title(title)
        for i, (c,ty) in enumerate(zip(self.channels(), self.channel_types)):
            axarr[i].plot(c.data.flatten(), *args, **kwargs)
            axarr[i].set_ylabel('Channel #%i\n (%s)' % (i + 1, ty))
        out = plt

        location, _ = plt.xticks()
        plt.xticks(location, [self._time_from_idx(l) for l in location], rotation=45)

        return out

    def show(self, *args, **kwargs):
        self.plot( *args, **kwargs).show()
    def plot(self, *args, **kwargs):
        """
        Plots the signal using :mod:`matplotlib.pyplot`.

        :param args: arguments to pass to :func:`~matplotlib.pyplot.plot`
        :param kwargs: keyword arguments to pass to :func:`~matplotlib.pyplot.plot`
        :return: the result of the :func:`matplotlib.pyplot.plot` function.
        """
        return self._create_fig(*args, **kwargs)


    def _repr_png_(self):
        from IPython.core.pylabtools import print_figure
        from matplotlib import pyplot as plt

        fig = self._create_fig()
        data = print_figure(fig, 'png')
        plt.close(fig)
        return data

    @property
    def png(self):
        from IPython.display import Image
        return Image(self._repr_png_(), embed=True)

    def save(self, filename, compression_level=5):
        pkl.dump(self,filename, compress=compression_level)


# annotations = np.random.normal(10,10,(100000,2)) + 1j
# a = AnnotatedTimeSeries( np.random.normal(10,100,(100000,4)),256, annotations = annotations, channel_types=["EEG","EMG","NaN","EEG"], annotation_types=["a", "b"] ,metadata={"a":1,"b":"yoyo"})
