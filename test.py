"""
Falko Giepmans, 05-08-2020

This script contains several functions to validate the algorithm for computing
three-term connection coefficients.
"""

from CC3 import threeterm_connection_coefficients as cc3
from CC2 import twoterm_connection_coefficients as cc2

import itertools
import pywt

import numpy as np

from cc3LIB import *
# =============================================================================
# 
# =============================================================================

def cc3dict( a, d1,d2,d3 ):
    """
    Gives the three-term connection coefficients for a wavelet in a dictonary.

    Parameters
    ----------
    a : np.array
        Wavelet filter.
    d1,d2,d3 : Int
        Order of derivative


    Returns
    -------
    ccdict : dictonary
        Dictonary of all non-zero connection coefficients with (l,m) as keys

    """
    idx, Tind, CC = cc3(a,d1,d2,d3)
    ccdict = {}
    for ii in range(0,len(Tind)):
            ccdict[Tind[ii]] = CC[ii]
    return ccdict

def CheckWithData( Data, a, d1, d2, d3 ):
    """
    Checks whether the computed connection coefficients match a given set of
    connection coefficients.

    Parameters
    ----------
    Data : Dictonary
        Dictonary of the known connection coefficient values.
    a : np.array
        Wavelet filter.
    d1,d2,d3 : Int
        Order of derivative

    Returns
    -------
    None

    """
    
    # Dictonary of the computed connection coefficients
    CC3_d1d2d3 = cc3dict( a, d1, d2, d3 )
    
    for key in Data.keys():
        a = Data[key]
        b = CC3_d1d2d3[key] if key in CC3_d1d2d3.keys() else 0
        
        print( '{}  \t {:.2e}'.format( key, ( a - b)/( b if abs(b) > 10**-10 else 1 ) ) )
    


def CheckSum( a, d1, d2 ):
    """
    Checks whether the sum over m of CC^(d1,d2,0)_(l,m) = CC^(d1,d2)_(l),
    where the latter is the two-term connection coefficients.
    
    INPUT:
    a       :  np.array         Wavelet filter
    d1,d2   :  int              Order of derivative
    
    OUTPUT:
    Error   :  np.array         
    """
    N        = a.size
    
    """
    The two-term connection coefficients and three-term connection coefficients:
    """
    CC2_d1d2  = (-1)**d1*cc2( a, d1+d2 )
    CC3_d1d20 = cc3dict( a, d1, d2, 0 )
    
    sumrow = []
    """
    Summation over the three-term connection coefficients:
    """
    for l in range(-(N-2), (N-2)+1):
        summation = 0
        for m in range(-(N-2), (N-2)+1):
            summation  += CC3_d1d20[(l, m)] if abs(l - m) < N-1 else 0
        sumrow.append( summation )
    
    sumrow = np.array( sumrow ) 
    
    
    Error  = abs( ( sumrow - CC2_d1d2 ) )/ ( sumrow + 10**-10 )
    print( "The two-term connection coefficients are:" )
    print( np.array_str( CC2_d1d2, precision = 2) )
    print( "Relative error for the check-sum procedure is:" )
    print( np.array_str( Error, precision = 2) )
    return Error

if __name__ == "__main__":
    """
    Pywt uses a different normalization and therefore we need to multiply the 
    arrays with 2**0.5 for the algorithm to work properly.
    """
    db2   = np.array( pywt.Wavelet( 'db2' ).filter_bank[2] ) * 2**0.5
    Error = CheckSum( db2, 1, 0 )
    
    db3   = np.array( pywt.Wavelet( 'db3' ).filter_bank[2] ) * 2**0.5
    CheckWithData( db3_100, db3, 1, 0, 0 )
    
