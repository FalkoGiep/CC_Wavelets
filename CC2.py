"""
Falko Giepmans, 05-08-2020

This is a code to calculate two-term connection coefficients for all(?) finite
supported discrete wavelet families on an unbounded interval.

References:
Fukuda, N., On the Wavelet-Galerkin Method with Deslauriers-Dubuc interpolating
    Scaling functions, 2013
Goedecker, S., Wavelets and their application for the solution of partial 
    differential equations in physics, 2009
"""

import itertools
import warnings
import scipy

import numpy as np

from   scipy.special import factorial


# =============================================================================
# 
# =============================================================================

def moment(a, i, j):
    """
    The moment M^j_i for all(?) scaling functions (Goedecker, 2009)
    
    INPUT:
    a       :  np.array         Wavelet filter
    i       :  int
    j       :  int
    
    OUTPUT:
    M^j_i   :  float
    """
    N = a.size
    return (i-(int( N/2)-1))**j 


def twoterm_connection_coefficients( a, d ):
    """
    Calculates the two-term connection coefficients CC^(0,d) of a wavelet.
    
    INPUT:
    a       :  np.array         Wavelet filter
    d       :  int              Order of derivative
    
    OUTPUT:
    CC      :  np.array         non-zero connection coefficients in order
    """
    
    """
    Eigen-value problem due to the scaling equations using the auto-correlation
    of the wavelet filter (Fukuda, 2013):
    """
    a_c = np.correlate( a, a, mode = "full")
    N_c = len(a_c)
    N   = N_c - 2
    T = np.zeros((N,N))
    for i,j in itertools.product(range(N), repeat=2): 
            if  -1 < j - 2*i + N < N_c:
                T[i,j] = a_c[ j - 2*i + N ]
    
    T -= 2**(1-d)*np.eye(N)
    b = np.zeros([N])   
    
    """
    Since the eigenvector is determined up to a constant, we alse need a 
    normalization equation ( Goedecker, 2009):
    """
    M = np.zeros([1, N])
    for i in range(0,N):
        M[0,i] += moment(a, i, d) 
    A = np.vstack([T,M])
    b = np.hstack([b, [factorial(d)]])
    
    """
    A least squares algorithm is used to solve the over-determined system.
    One can also use np.linalg.lstsq with rcond = None. In my experience 
    however, np.linalg.lstsq does not always return residuals correctly.
    """
    CC, residuals, rank, singular_values  = scipy.linalg.lstsq(A, b)
    
    if abs( residuals ) >= 10**-30:
        msg = 'Residue of lstsq algorithm is {:.2e}!'.format(residuals)
        warnings.warn(msg)
        
    return CC

