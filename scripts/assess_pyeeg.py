__author__ = 'quentin'


"""
Issues with pyeeg:
General: no use of vectoristation and high level numpy functions such as np.diff (instead of manual implementation of numerical derivatives)

Specific:
    * Hjorth parameters are not normalised by the length of the time series.

    * HFD is for Higuchi, not Hjorth Fractal Dimension.
                .. math::

            \frac{N-1}{[\frac{N-m}{k}]}(in pyeeg) instead of  \frac{N-1}{[\frac{N-m}{k} ]\dot{} k}


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

import pyeeg
import pyrem.univariate as univ
import numpy
from timeit import Timer
import pandas as pd
import numpy as np
import sys


MIN_EPOCH_N = 256 * 5
MAX_EPOCH_N = 256 * 30
EPOCH_STEP = 256 * 5
N_REPLICATES = 5

SPECT_ENT_BANDS = 2 ** np.arange(0,8)/2

fun_to_test = [
                  {"times":100,"name":"hfd", "is_original":True,"fun": lambda x: pyeeg.hfd(x,2**3)},
                  {"times":100,"name":"hfd", "is_original":False,"fun": lambda x: univ.hfd(x,2**3)},
                  {"times":100,"name":"hjorth", "is_original":True,"fun": lambda x: pyeeg.hjorth(x)},
                  {"times":100,"name":"hjorth", "is_original":False,"fun": lambda x: univ.hjorth(x)},
                  {"times":100,"name":"pfd", "is_original":True, "fun":lambda x: pyeeg.pfd(x)},
                  {"times":100,"name":"pfd", "is_original":False, "fun":lambda x: pyeeg.pfd(x)},
                  {"times":2,"name":"samp_ent", "is_original":True, "fun":lambda x: pyeeg.samp_entropy(x,2,1.5)},
                  {"times":10,"name":"samp_ent", "is_original":False, "fun":lambda x: univ.samp_entropy(x,2,1.5,relative_r=False)},
                  {"times":2,"name":"ap_ent", "is_original":True, "fun":lambda x: pyeeg.ap_entropy(x,2,1.5)},
                  {"times":10,"name":"ap_ent", "is_original":False, "fun":lambda x: univ.ap_entropy(x,2,1.5)},
                  {"times":10,"name":"svd_ent", "is_original":True, "fun":lambda x: pyeeg.svd_entropy(x,2,3)},
                  {"times":100,"name":"svd_ent", "is_original":False, "fun":lambda x: univ.svd_entropy(x,2,3)},
                  {"times":10,"name":"fisher_info", "is_original":True, "fun":lambda x: pyeeg.fisher_info(x,2,3)},
                  {"times":100, "name":"fisher_info", "is_original":False, "fun":lambda x: univ.fisher_info(x,2,3)},
                  {"times":100,"name":"spectral_entropy", "is_original":True, "fun":lambda x: pyeeg.spectral_entropy(x,SPECT_ENT_BANDS,256)},
                  {"times":100, "name":"spectral_entropy", "is_original":False, "fun":lambda x: univ.spectral_entropy(x,256, SPECT_ENT_BANDS)},

    ]


def make_one_rep():
    ldfs = []
    for n in range(MIN_EPOCH_N, MAX_EPOCH_N + 1, EPOCH_STEP):
        a = numpy.random.normal(size=n)
        for fun in fun_to_test:
            f = lambda : fun["fun"](a)
            t=Timer(f)
            numb = fun["times"]
            dt = t.timeit(number=numb)
            fun["log10_dt"] = np.log10(dt/float(numb))
            fun["n"] = n

        df = pd.DataFrame(fun_to_test)
        ldfs.append(df)
        print n
    return pd.concat(ldfs)

if __name__ == "__main__":
    out_file_name = sys.argv[1]

    df = pd.concat([make_one_rep() for i in range(N_REPLICATES )])
    df.to_csv(out_file_name, float_format="%e")





r"""
rm(list=ls())

library(ggplot2)
library(nlme)

make_coeffs <- function(df, is_ori){

    out <- data.frame(t(sapply(lmList(log10_dt ~ n | name  , subset(df, is_original == is_ori)), function(c)c$coefficients)))
    print(out)
    colnames(out) <- c("b","a")

    out$orig <- out$b + out$a * 256*5
    out$b <- NULL
    out$a <- round(1/out$a)
    out$orig <- round(10^(out$orig + 3),3)
    colnames(out) <- sprintf("%s_%s", is_ori, c("a","tn1000"))
    out
}
df <- read.csv("./assess_pyeeg_out_2.csv")
df$n <- df$n - 1280
means_df <- aggregate(log10_dt ~ n* name* is_original, df, mean)
ggplot(means_df, aes(x=n, y=log10_dt, colour=name, linetype=is_original)) + geom_line() + geom_point()
means_df <- aggregate(log10_dt ~ n* name* is_original, df, mean)

original <- make_coeffs(df, "True")
pr <- make_coeffs(df, "False")
cbind(original,pr)



summary(lmList(log10_dt ~ n * is_original | name ,df))
"""