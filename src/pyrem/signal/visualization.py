__author__ = 'quentin'
from datetime import timedelta
from scipy import stats

from scipy.ndimage.interpolation import zoom
import numpy as np
import pylab as pl
from pyrem.signal.signal import Signal, Annotation
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as plt


class PolygramDisplay(object):
    def __init__(self, polygram, max_point_amplitude_plot=1000):

        self.polygram = polygram
        self.max_point_amplitude_plot = max_point_amplitude_plot
        self.fig, self.axarr = pl.subplots(self.polygram.n_channels, sharex=True)

        self.fig.subplots_adjust(hspace=0)


        self._redraw(None, init=True)
        self._redraw(None)
        #ax2.plot(np.linspace(0,a.duration.total_seconds(),a.size),a)
        #self.ax_update(ax2)

        pl.show()


    def _redraw(self, _, init=False):

        for ax, sig in zip(self.axarr, self.polygram.channels):
            if not init:
                ax.clear()
                #ax.set_autoscale_on(False) # Otherwise, infinite loop
                ax.autoscale(enable=False, axis='x')
                ax.autoscale(enable=True, axis='y')
                ax.callbacks.connect('xlim_changed', self._redraw)

            if isinstance(sig, Signal):
                self._plot_signal_on_ax(sig, ax, init)
            elif isinstance(sig, Annotation):
                self._plot_annotation_on_ax(sig, ax,init)
            else:
                raise ValueError("The time series is a %s" % str(type(sig)))
            #pl.setp([ax.get_xticklabels()], visible=False)
            axis_title = "%s\n(@%sHz)" % (sig.name, str(round(sig.fs,3)))
            ax.set_ylabel(axis_title)

    def _plot_annotation_on_ax(self, signal, ax, autoscale=False, colourmap="flag"):

        if autoscale:
            xstart = 0
            xdelta = signal.duration.total_seconds()
        else:

            xstart,ystart,xdelta,ydelta = ax.viewLim.bounds

        if xstart <0:
            start_time = timedelta()
        else:
            start_time = timedelta(seconds=xstart)

        stop_time = timedelta(seconds=xdelta) +  timedelta(seconds=xstart)
        sub_sig = signal[start_time:stop_time]
        xs =np.linspace(0, sub_sig.duration.total_seconds() ,sub_sig.size) + start_time.total_seconds()
        ys = sub_sig.values
        probs = sub_sig.probas

        ys = ys.reshape((1,ys.size))

        zoom_f = float(self.max_point_amplitude_plot)/ sub_sig.size

        ys = zoom(ys,[1, zoom_f], order=0)


        ax.imshow(ys, extent=[np.min(xs), np.max(xs), 1.5, -0.5], aspect="auto",
                  cmap=colourmap, vmin=0, vmax=255, origin='lower')


        ax.plot(xs,probs,"-", color="k", linewidth=3)
        ax.plot(xs,probs,"-", color="y", linewidth=1,alpha=0.5)


        jet = cm = pl.get_cmap(colourmap)
        cNorm  = colors.Normalize(vmin=0, vmax=255)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        states = np.unique(ys)

        boxes = [pl.Rectangle((0, 0), 1, 1, fc=scalarMap.to_rgba(col)) for col in states]
        labels = [chr(s) for  s in states]
        pl.legend(boxes,labels, loc='lower right')

        n_labels = 8 #fixme magic number

        if len(xs) > n_labels:
            trimming = int(float(len(xs)) / float(n_labels))
            xs_trimmed = np.round(xs[::trimming])
        else:
            xs_trimmed = xs

        time_strings = [str(timedelta(seconds=s)) for s in xs_trimmed]


        ax.set_xticks(xs_trimmed)
        ax.set_xticklabels(time_strings, rotation=70)

        return

    def _plot_signal_on_ax(self, signal, ax, autoscale=False):

        if autoscale:
            xstart = 0
            xdelta = signal.duration.total_seconds()
        else:

            xstart,ystart,xdelta,ydelta = ax.viewLim.bounds
        n_viewed_points = xdelta * signal.fs

        if n_viewed_points < self.max_point_amplitude_plot*5:
            if xstart <0:
                start_time = timedelta()
            else:
                start_time = timedelta(seconds=xstart)

            stop_time = timedelta(seconds=xdelta) +  timedelta(seconds=xstart)
            sub_sig = signal[start_time:stop_time]
            xs =np.linspace(0, sub_sig.duration.total_seconds() ,sub_sig.size) + start_time.total_seconds()
            ax.plot(xs, sub_sig,"-", linewidth=1, color=(0,0,1,0.5))

            return



        winsize_npoints = float(n_viewed_points) / float(self.max_point_amplitude_plot)
        secs = winsize_npoints / signal.fs


        start = int(xstart  * signal.fs)
        if start <0:
                start=0
        stop = int((xstart + xdelta) * signal.fs)
        sub_sig = signal[start:stop]


        mins, maxes, means, sds,xs = [],[],[],[],[]
        for c, w in  sub_sig.iter_window(secs,1):

            mins.append(np.min(w))
            maxes.append(np.max(w))
            means.append(np.mean(w))
            sds.append(np.std(w))
            xs.append(c+start / sub_sig.fs)

        means = np.array(means)
        mean_plus_sd = means +sds
        mean_minus_sd = means - sds


        ax.fill_between(xs,mins,maxes, facecolor=(0,0,1,0.6),edgecolor=(0,0,0,0.2), antialiased=True)
        ax.fill_between(xs,mins,maxes, facecolor=(0,0,1,0.6),edgecolor=(0,0,0,0.2), antialiased=True)
        ax.fill_between(xs,mean_minus_sd, mean_plus_sd, facecolor=(1,0.5,0,0.9),edgecolor=(0,0,0,0), antialiased=True)
        ax.plot(xs,means,"-", linewidth=1, color='k')

        n_labels = 8 #fixme magic number

        if len(xs) > n_labels:
            trimming = int(float(len(xs)) / float(n_labels))
            xs_trimmed = np.round(xs[::trimming])
        else:
            xs_trimmed = xs

        time_strings = [str(timedelta(seconds=s)) for s in xs_trimmed]


        ax.set_xticks(xs_trimmed)
        ax.set_xticklabels(time_strings, rotation=70)