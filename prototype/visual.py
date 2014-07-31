
from pyrem.signal.io import polygram_from_pkl
from pyrem.signal.visualization import PolygramDisplay
from pyrem.signal.polygram import Polygram




pol = polygram_from_pkl("/data/pyrem/Ellys/pkls/TelC_2.pkl")
#pol = Polygram(pol.channels[0:2])
md = PolygramDisplay(pol)
