from pyrem.signal.signal  import *
import numpy as np
import pylab as pl
MAX_POINTS_AMPLITUDE_PLOT = 1000
np.random.seed(1)
v = np.random.normal(0,1,(int(1e6)))

a = Signal(v, 10,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})

class SignalDisplay(object):
    def __init__(self, signal, h=500, w=500):
        self.height = h
        self.width = w
        self.signal = signal

    
    def __call__(self,max_points=MAX_POINTS_AMPLITUDE_PLOT):
        pass
        
    def ax_update(self, ax):
        ax.clear()
        ax.callbacks.connect('xlim_changed', md.ax_update)
        #~ ax.callbacks.connect('ylim_changed', md.ax_update)

        ax.set_autoscale_on(False) # Otherwise, infinite loop
        
        #Get the range for the new area
        xstart,ystart,xdelta,ydelta = ax.viewLim.bounds
        n_viewed_points = xdelta * self.signal.fs
        print n_viewed_points
        if n_viewed_points < MAX_POINTS_AMPLITUDE_PLOT*5:
            start = int(xstart  * self.signal.fs)
            if start <0:
                start=0
            stop = int((xstart + xdelta) * self.signal.fs)
            sub_sig = self.signal[start:stop]
            xs =np.linspace(start, sub_sig.duration.total_seconds() ,sub_sig.size)
            ax.plot(xs,sub_sig,"-", linewidth=1, color='k')
            return
            
        max_points=MAX_POINTS_AMPLITUDE_PLOT 
        
        winsize_npoints = float(n_viewed_points) / float(max_points)
        secs = winsize_npoints / self.signal.fs
        print secs
        
        start = int(xstart  * self.signal.fs)
        if start <0:
                start=0
        stop = int((xstart + xdelta) * self.signal.fs)
        sub_sig = self.signal[start:stop]
        #xxs =np.linspace(start, sub_sig.duration.total_seconds() ,sub_sig.size)
            
        mins, maxes, means, sds,xs = [],[],[],[],[]
        for c, w in  sub_sig.iter_window(secs,1):
            mins.append(np.min(w))
            maxes.append(np.max(w))
            means.append(np.mean(w))
            sds.append(np.std(w))
            xs.append(c+start / sub_sig.fs)
        
        print xs[0]
        means = np.array(means)
        mean_plus_sd = means +sds
        mean_minus_sd = means - sds
        
        
        ax.fill_between(xs,mins,maxes, facecolor=(0,0,1,0.6),edgecolor=(0,0,0,0.2), antialiased=True)
        ax.fill_between(xs,mins,maxes, facecolor=(0,0,1,0.6),edgecolor=(0,0,0,0.2), antialiased=True)
        ax.fill_between(xs,mean_minus_sd, mean_plus_sd, facecolor=(1,0.5,0,0.9),edgecolor=(0,0,0,0), antialiased=True)
        ax.plot(xs,means,"-", linewidth=2, color='k')

        n_labels = 8 #fixme magic number

        time_strings = [str(timedelta(seconds=s)) for s in xs]
        if len(time_strings) > n_labels:
            trimming = int(float(len(time_strings)) / float(n_labels))
            xs = xs[::trimming]
            time_strings = time_strings[::trimming]


        #~ ax.set_xticks(xs, time_strings, rotation=45)
        
        #Get the number of points from the number of pixels in the window
        dims = ax.axesPatch.get_window_extent().bounds
        self.width = int(dims[2] + 0.5)
        self.height = int(dims[2] + 0.5)
        return
        
md = SignalDisplay(a)
md()
fig1, ax2 = pl.subplots(1, 1)
ax2.plot(np.linspace(0,a.duration.total_seconds(),a.size),a)
md.ax_update(ax2)

pl.show()
   
