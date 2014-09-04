__author__ = 'quentin'


from pyrem.io import polygram_from_pkl
from pyrem.wavelet_decomposition import decompose_signal
from pyrem.time_series import Annotation
from pyrem.time_series import Signal
from pyrem.polygram import Polygram

import pylab as pl
import pandas as pd
import numpy as np

# DATA_FILE_PATTERN=


#df = pd.read_csv("/tmp/telc4_res.csv")
#an = Annotation(df["pred"], 0.2, df["conf_preds"], name="prediction")

#pol1  = Polygram([an])

pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/TelC_4.pkl")

#pol = pol1.merge(pol)

z = decompose_signal(pol["EEG_parietal_frontal"], [1,2,3,4,5,6])
z = z.merge(pol["EEG_parietal_frontal"])
z["16h45m":"16h55m"].show()

dwtp = decompose_signal(pol["EEG_parietal_frontal"], [6], keep_a=False)

t = dwtp[0] ** 2

N=25
tt = np.log10(np.convolve(t, np.ones((N,))/N,"same"))
t = Signal(tt,t.fs,name="Power in cD_6")
print t.duration

# dwt_expl = dwtp["13h41m30s":"14h44m05s"]
dwt_expl = Polygram([t])

pol_small = pol
# pol_small = pol["13h41m30s":"14h44m05s"]

dwt_expl = dwt_expl.merge(pol_small["vigilance_state"])
dwt_expl = dwt_expl.merge(pol_small["EMG_REF"])
dwt_expl = dwt_expl.merge(pol_small["EEG_parietal_frontal"])
dwt_expl = dwt_expl.merge(pol_small["prediction"])

# dwt_expl.show()
# dwt_expl["2h55m":"3h20m"].show()

# wake_expl = pol["3h19m34s":"3h19m39s"]
# nrem_expl = pol["3h18m40s":"3h18m45s"]









#
# def plot_one_trace(refs,str):
#     xs =[ p[str] for p in refs]
#
#     f, (ax1, ax2, ax3) = pl.subplots(3, sharex=True, sharey=True)
#     ax1.plot(xs[0])
#     ax2.plot(xs[1])
#     ax3.plot(xs[2])
#
#     f.subplots_adjust(hspace=0)
#     pl.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#     return f

# from matplotlib.backends.backend_pdf import PdfPages
#
# refs = [nrem_expl, rem_expl, wake_expl]
# pp = PdfPages("/tmp/eeg_emg_expls.pdf")
#
# pp.savefig(plot_one_trace(refs,str="EEG_parietal_frontal"))
# pp.savefig(plot_one_trace(refs,str="EMG_REF"))
# pp.close()

