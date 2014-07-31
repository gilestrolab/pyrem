__author__ = 'quentin'

from pyrem.signal.wavelet_decomposition import decompose_signal
from pyrem.signal.io import polygram_from_pkl
from pyrem.signal.visualization import PolygramDisplay

pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/TelC_2.pkl")

sig = pol["EEG_parietal_frontal"].resample(256.0)

pol2 = decompose_signal(sig,levels=[3,4,5])

PolygramDisplay(pol2)