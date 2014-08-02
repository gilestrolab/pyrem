from pyrem.signal.signal import Signal
import numpy as np
import pylab as pl
from scipy.ndimage.filters import gaussian_filter
from joblib import Memory

from tempfile import mkdtemp
cachedir = mkdtemp()

memory = Memory(cachedir=cachedir, verbose=0)

np.random.seed(1)
a = np.random.normal(0,1,(int(2.4e7)))
a = a * np.linspace(1,2,a.size) + np.linspace(0,15,a.size) + np.sin(np.linspace(0,12,a.size)) 

sub_sig = Signal(a, 256,type="eeg", name="foo", metadata={"animal":"joe", "treatment":18})


def make_hist(w, rang):
    return np.histogram(w,bins=w.size/10, range=rang)[0]

@memory.cache
def make_img(sub_sig):
    s = sub_sig.duration.total_seconds() / 512.0
    rang = (np.min(sub_sig), np.max(sub_sig)) 
    l = []
    ts = []
    xs = []
    for c, w in  sub_sig.iter_window(s,1):
        ts.append(c)
        xs.append(np.mean(w))    
        l.append(make_hist(w,rang))
    img = np.log10(np.array(l).T +1)
    img = gaussian_filter(img, 3)
    return img, rang, ts, xs
    
def hist(sub_sig):
    img, rang, ts, xs = make_img(sub_sig)
    #~ pl.imshow(img, extent=[0, 1, rang[0], rang[1]], aspect="auto", origin='lower')
    pl.imshow(img, extent=[min(ts), max(ts), rang[0], rang[1]], aspect="auto", origin='lower', cmap="gray")
    pl.plot(ts, xs, '-')
    #~ pl.show()
    
hist(sub_sig)
pl.show()
hist(sub_sig[0:int(1e6)])
hist(sub_sig[0:int(1e5)])

