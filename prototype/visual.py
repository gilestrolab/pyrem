from pyrem.signal.signal  import *
from pyrem.signal.polygram import *
from pyrem.signal.visualization import PolygramDisplay


import numpy as np

np.random.seed(1)
rw1 = np.cumsum(np.random.normal(0,5,int(1e5)))
rw2 = np.random.normal(0,2,int(1e5))
rw2 = rw2 * np.linspace(1,3, rw2.size)
rw3 = np.random.normal(0,2,int(1e5/2))
rw3 = rw3 / np.linspace(1,3, rw3.size)
probs = np.random.random_sample(100)
vals = (np.random.random_sample(100) * 4 +1).astype(np.int)

an = Annotation(vals,fs=0.256, observation_probabilities=probs, type="vigilance_state", name="an", metadata={"animal":"joe", "treatment":18})
c1 = Signal(rw1, fs=256,type="eeg", name="c1", metadata={"animal":"joe", "treatment":18})
c2 = Signal(rw2, fs=256,type="eeg", name="c2", metadata={"animal":"joe", "treatment":18})
c3 = Signal(rw3, fs=256/2,type="eeg", name="c3", metadata={"animal":"joe", "treatment":18})

pol = Polygram([c1,c2, c3, an])


md = PolygramDisplay(pol)
