# -*- coding: utf-8 -*-
import numpy as np
def score_spectra( X, method='DEuc' ):
    """Returns the correaltion coefficient matrix of the variables in X.

    Returns the correlation coefficient matrix of the varibles in X.
    Allows correlation between derivates.
    
    Parameters
    ----------
    X : array_like
        A 2-D array containing multiple variables and observations. 
        Each row of X representes a variable, and each column a single observation of all those variables.
    method : str, optional
        Correlation method, valid entries: Cor, DCor, Euc and DEuc
        Cor: Correlation coefficient
        DCor: Squared first difference correlation coefficient
        Euc: Euclidean cosine
        DEuc: Squared first difference Euclidean cosine
    
    Returns
    ----------
    Rs: ndarray
        The correaltion coefficient matrix of the variables
    """
    Rs = np.zeros( ( len( X ), len( X ) ) )
    for c1, A1 in enumerate( X ):
        dA1 = np.diff( A1 )
        mA1 = A1 - np.mean( A1 )
        mdA1 = dA1 - np.mean( dA1 )
        for c2, A2 in enumerate( X ):
            dA2 = np.diff( A2 )
            mA2 = A2 - np.mean( A2 )
            mdA2 = dA2 - np.mean( dA2 )
            if method=='Cor':
                Rs[c1][c2] = ( np.sum( ( mA1 ) * ( mA2 ) )**2 ) / ( np.sum( mA1**2 ) * np.sum( mA2**2 ) )
            elif method=='DCor':
                Rs[c1][c2] = ( ( np.sum( ( mdA1 ) * ( mdA2 ) ) )**2 ) / ( np.sum( mdA1**2 ) * np.sum( mdA2**2 ) )
            elif method=='Euc':
                Rs[c1][c2] = ( np.sum( A1 * A2 )**2 ) / ( np.sum( A1**2 ) * np.sum( A2**2 ) )
            elif method=='DEuc':
                Rs[c1][c2] = ( np.sum( dA1 * dA2 )**2 ) / ( np.sum( dA1**2 ) * np.sum( dA2**2 ) )
    return Rs