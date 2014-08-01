__author__ = 'quentin'



from pyrem.signal.io import polygram_from_pkl
from pyrem.signal.visualization import PolygramDisplay
from pyrem.signal.polygram import Polygram




pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/GFP_sleep.pkl")
c1 = pol[0]
c2 = pol[0].resample(256.0)
c1.rename("boom")

print c2

pol = Polygram([c1,c2])
md = PolygramDisplay(pol)