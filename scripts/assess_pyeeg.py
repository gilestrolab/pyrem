__author__ = 'quentin'


"""
Issues with pyeeg:
General: no use of vectoristation and high level numpy functions such as np.diff (instead of manual implementation of numerical derivatives)

Specific:
    * Hjorth parameters are not normalises by the length of the time series.

    * HFD is for Higuchi, not Hjorth Fractal Dimention.

    * PFD returns:
        log10(n)/(log10(n)+log10(n/n+0.4*N_delta))
        From Petrosian 1995:
        log10(n)/(
            log10(n) + log10(n/  (   n+0.4*N_delta  ) )
            )
    * SVD ent and Spectral entropy:
        should be log2 not log
    *spectral entropy. the power spectrum should be the **square** of the absolute fft
"""