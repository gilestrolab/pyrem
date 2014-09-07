"""
Pyrem is a python package to help feature extraction from physiological time series such as electroencephalograms(EEGs) and such.

It provides:

1. Data structures for time series (:mod:`pyrem.time_series`) based on numpy arrays. This extends functionality of numpy arrays by:

    a. Providing extra attributes such as sampling frequency and metadata
    b. Allow new functionality such as time-string indexing (e.g. :code:`signal["28m":"2h22.4s"]`)

2. A data structure for annotations (:class:`pyrem.time_series.Annotation`) based on numpy recarrays. This allows to describe time series of discrete states linked to a confidence/probability of observation of each states.

3. A data structure for collection of time series (:mod:`pyrem.polygram`). It features:

    a. Support for heterogeneous sampling rates between time series.
    b. An iterator through arbitrary epochs (temporal slices)
    c. Slicing using time strings and seamless merging between polygrams

4. Implementations of algorithms often used in analysis of EEGs (see :mod:`pyrem.univariate`). They are essentially  faster and curated reimplementation of the function available in PyEEG_.

5. Utilities to load, save and visualize the data (still in development).

6. Wrappers around samplerate_ and pywt_ libraries, to efficiently resample, and compute discrete wavelet decomposition on signals.

.. _PyEEG: http://www.hindawi.com/journals/cin/2011/406391/
.. _samplerate: http://www.ar.media.kyoto-u.ac.jp/members/david/softwares/samplerate/sphinx/index.html
.. _pywt: http://www.pybytes.com/pywavelets/

"""


import univariate
import time_series
from time_series import Signal, Annotation

import polygram
from polygram import Polygram

import feature_families


import wavelet_decomposition
