__author__ = 'quentin'

cimport numpy as np
import numpy as np

DTYPE = np.float
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float_t DTYPE_t
# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.
#


def embed_seq(X,Tau,D):

    N =len(X)

    if D * Tau > N:
        raise ValueError("Cannot build such a matrix, because D * Tau > N")

    if Tau<1:
        raise ValueError("Tau has to be at least 1")


    Y=np.zeros((D, N - (D - 1) * Tau))

    for i in range(D):
        Y[i] = X[i *Tau : i*Tau + Y.shape[1] ]

    return Y.T

def in_range(np.ndarray Template, np.ndarray Scroll, float  Distance):
    cdef int N = len(Template)
    cdef int i

    for i in range(0,  N):
            if abs(Template[i] - Scroll[i]) > Distance:
                 return False
    return True



def ap_entropy(np.ndarray X, int M, float R):

    cdef int N = len(X)

    cdef np.ndarray  Em = embed_seq(X, 1, M)
    cdef np.ndarray Emp = embed_seq(X, 1, M + 1)


    cdef np.ndarray Cm = np.zeros(N - M + 1, dtype=DTYPE)
    cdef np.ndarray Cmp = np.zeros(N - M, dtype=DTYPE)

#     # in case there is 0 after counting. Log(0) is undefined.
    cdef int i,j
    for i in xrange(0, N - M):
        for j in xrange(i, N - M): # start from i, self-match counts in ApEn
#             print i,j
#             print np.array([Em[i], Em[j]])

            if in_range(Em[i], Em[j], R):
                Cm[i] += 1                                                                                        ### Xin Liu
                Cm[j] += 1
                if np.abs(Emp[i][-1] - Emp[j][-1]) <= R: # check last one
                    Cmp[i] += 1
                    Cmp[j] += 1
        if in_range(Em[i], Em[N-M], R):
            Cm[i] += 1
            Cm[N-M] += 1
        # try to count Cm[j] and Cmp[j] as well here

#        if max(abs(Em[N-M]-Em[N-M])) <= R: # index from 0, so N-M+1 is N-M  v 0.01b_r1
#    if in_range(Em[i], Em[N - M], R):  # for Cm, there is one more iteration than Cmp
#            Cm[N - M] += 1 # cross-matches on Cm[N - M]

    Cm[N - M] += 1 # Cm[N - M] self-matches
#    import code;code.interact(local=locals())
    Cm /= (N - M +1 )
    Cmp /= ( N - M )
#    import code;code.interact(local=locals())
    Phi_m, Phi_mp = sum(np.log(Cm)),  sum(np.log(Cmp))

    Ap_En = (Phi_m - Phi_mp) / (N - M)

    return Ap_En
