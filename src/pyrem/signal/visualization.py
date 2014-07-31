__author__ = 'quentin'
from datetime import timedelta
import numpy as np
import pylab as pl



class PolygramDisplay(object):
    def __init__(self, polygram, max_point_amplitude_plot=1000):

        self.polygram = polygram
        self.max_point_amplitude_plot = max_point_amplitude_plot
        self.fig, self.axarr = pl.subplots(self.polygram.n_channels, sharex=True)

        for ax, sig in zip(self.axarr, self.polygram.channels):
            self._plot_signal_on_ax(sig, ax,True)

        self._redraw(None)
        #ax2.plot(np.linspace(0,a.duration.total_seconds(),a.size),a)
        #self.ax_update(ax2)

        pl.show()


    def _redraw(self, _):
        for ax, sig in zip(self.axarr, self.polygram.channels):
            ax.clear()
            ax.callbacks.connect('xlim_changed', self._redraw)
            #ax.set_autoscale_on(False) # Otherwise, infinite loop
            ax.autoscale(enable=False, axis='x')
            ax.autoscale(enable=True, axis='y')

            self._plot_signal_on_ax(sig, ax)



    def _plot_annotation_on_ax(self, signal, ax, autoscale=False):
        pass

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

        time_strings = [str(timedelta(seconds=s)) for s in xs]
        if len(time_strings) > n_labels:
            trimming = int(float(len(time_strings)) / float(n_labels))
            xs = xs[::trimming]
            time_strings = time_strings[::trimming]
