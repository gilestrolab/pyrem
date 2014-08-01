__author__ = 'quentin'



from pyrem.signal.wavelet_decomposition import decompose_signal
from pyrem.signal.io import polygram_from_pkl
from pyrem.signal.visualization import PolygramDisplay
import numpy as np

pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/GFP_A.pkl")

print pol.channel_names
pol2 = pol.map_signal_channels(lambda mat: (mat - np.mean(mat,0)) / np.std(mat,0))
print pol2.channel_names

PolygramDisplay(pol2)