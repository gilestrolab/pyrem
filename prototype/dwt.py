__author__ = 'quentin'

from pyrem.signal.wavelet_decomposition import decompose_signal
from pyrem.signal.io import polygram_from_pkl
from pyrem.signal.visualization import PolygramDisplay

pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/GFP_A.pkl")

print pol.channel_names
sig = pol["EEG_parietal_frontal"]

print "applying ressampling + DWTD..."
pol2 = decompose_signal(sig,levels=[3,4,5])
for i in pol2.signal_channels:
    print i.name, i.duration, i.fs

i = sig
print i.name, i.duration, i.fs

pol2 = pol2.merge(sig)
# pol2.append_channel(pol["EMG_1"])
# pol2.append_channel(pol["EMG_2"])
pol2 = pol2.merge(pol["vigilance_state"])

print "Displaying..."
pol2.show()