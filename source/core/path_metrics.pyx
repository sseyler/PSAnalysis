import numpy as np
cimport numpy as np
from libc.math cimport fmin, fmax, sqrt
cimport cython
from numpy import amin, amax


def vsqnorm(v, axis=None):
    return (v*v).sum(axis=axis)


def rmsdMatrix(P, Q, axis=(1,2)):
    return np.array([vsqnorm(p - Q, axis=axis) for p in P])


def hausdorff(d, N):
    return ( max( amax(amin(d, axis=0)), amax(amin(d, axis=1)) ) / N  )**0.5


def frechet(d, N):
    Np, Nq = d.shape
    cd = -numpy.ones((Np, Nq))

    def cD(i, j):
        if cd[i,j] != -1 : return cd[i,j]
        if i > 0:
            if j > 0: cd[i,j] = max( min(cD(i-1,j),cD(i,j-1),cD(i-1,j-1)), d[i,j] )
            else:     cd[i,j] = max( cD(i-1,0), d[i,0] )
        elif j > 0:   cd[i,j] = max( cD(0,j-1), d[0,j] )
        else:         cd[i,j] = d[0,0]
        return        cd[i,j]

    return ( cD(Np-1, Nq-1) / N )**0.5


################################################################################
# RMSD Matrix
#---------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
cdef float[:,::1] rmsdMatrix_c(float[:,::1] P, float[:,::1] Q):
    cdef np.intp_t lenP = P.shape[0]
    cdef np.intp_t lenQ = Q.shape[0]
    cdef np.intp_t i, j, k
    cdef float s, diff
    cdef float[:,::1] d = np.empty((lenP, lenQ), dtype='float32')
    for i in xrange(lenP):
        for j in xrange(lenQ):
            s = 0.0
            for k in xrange(P.shape[1]):
                diff = P[i,k] - Q[j,k]
                s += diff*diff
            d[i,j] = s
    return d
#===============================================================================


################################################################################
# Hausdorff
#---------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def hausdorff(float[:,::1] P, float[:,::1] Q, np.intp_t N):
    assert (int(P.shape[1]) == int(3*N) and P.shape[1] == Q.shape[1])
    cdef float[:,::1] d = rmsdMatrix_c(P, Q)

    return sqrt( fmax( amax(amin(d, axis=0)), amax(amin(d, axis=1)) ) / N  )
#===============================================================================


################################################################################
# Frechet
#---------------------------------------
ctypedef float (*cD_ptr)(float[:,::1], float[:,::1], np.intp_t, np.intp_t)

cdef float fmin3(float a, float b, float c):
    return fmin(fmin(a,b),c)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef float cD(float[:,::1] d, float[:,::1] cd, np.intp_t i, np.intp_t j):
    cdef np.intp_t im1 = i-1
    cdef np.intp_t jm1 = j-1
    if cd[i,j] != -1 : return cd[i,j]
    if i > 0:
        if j > 0: cd[i,j] = fmax( fmin3(cD(d,cd,i,jm1),cD(d,cd,im1,jm1),cD(d,cd,im1,j)), d[i,j] )
        else:     cd[i,0] = fmax( cD(d,cd,im1,0), d[i,0] )
    elif j > 0:   cd[0,j] = fmax( cD(d,cd,0,jm1), d[0,j] )
    else:         cd[0,0] = d[0,0]
    return        cd[i,j]

@cython.boundscheck(False)
@cython.wraparound(False)
def frechet(float[:,::1] P, float[:,::1] Q, np.intp_t N):
    cdef np.intp_t lenP = P.shape[0]
    cdef np.intp_t lenQ = Q.shape[0]
    assert (int(P.shape[1]) == int(3*N) and P.shape[1] == Q.shape[1])
    cdef float[:,::1] d = rmsdMatrix_c(P, Q)
    cdef float[:,::1] cd = -np.ones((lenP, lenQ), dtype='float32')
    cdef cD_ptr couplingDistance = &cD

    return sqrt( couplingDistance(d, cd, lenP-1, lenQ-1) / N )
#===============================================================================