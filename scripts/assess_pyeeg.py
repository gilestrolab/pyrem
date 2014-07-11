__author__ = 'quentin'


"""
Issues with pyeeg:
General: no use of vectoristation and high level numpy functions such as np.diff (instead of manual implementation of numerical derivatives)

Specific:
    * Hjorth parameters are not normalises by the length of the time series.
    * HFD is for Higuchi, not Hjorth Fractal Dimention.
"""