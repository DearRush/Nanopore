# -*- coding: utf-8 -*-
import scipy
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from NDFs import *

def cluster_xdrNDF( X, n_bins=500 ):
    nx = np.linspace( 0, 1, n_bins + 1 )
    ny = np.linspace( 0, 10, n_bins + 1 )
    
    H, _, _ = np.histogram2d( X[:,0], X[:,1], 
                             bins=(nx, ny) )
    H = H.T
    nx = np.linspace( 0, 100, n_bins )
    ny = np.linspace( 0, 10, n_bins )
    Ires_H = np.sum( H, axis=0 )
    
    POPT_Ires = []
    POPT_beta = []
    h = np.copy( Ires_H )
    Peaks = scipy.signal.find_peaks( savgol_filter( Ires_H, 51, 2 ), height=np.max(Ires_H)/20 )[0]
    Px = []
    n_x = int( 5 * n_bins / 100 )
    for c, p, peak in sorted( zip( Ires_H[ Peaks ], nx[ Peaks ], Peaks ) ):
        try:
            bounds = ( ( 0, 0, 0 ), ( np.inf, 100, 10 ) )
            popt, _ = curve_fit( rNDF, nx[peak-n_x:peak+n_x], h[peak-n_x:peak+n_x], p0=[c,p,1], bounds=bounds)
            if rNDF( p, *popt ) > np.max(Ires_H)/20:
                h -= rNDF( nx, *popt )
                Px.append( popt[1] )
                POPT_Ires.append( popt )
                
                w = rNDF( nx, *popt )
                weight = []
                for i in range( H.shape[0] ):
                    weight.append( w )
                weight = np.array( weight )
                hy = np.average( H, axis=1, weights=weight)
                
                bounds = ( ( 0, 0, 0 ), ( np.inf, np.inf, np.inf ) )
                popt2, _ = curve_fit( rNDF, ny, hy, p0=[1,ny[ np.argmax( hy ) ],0.1])
                POPT_beta.append( popt2 )
        except:
            pass
    p0 = []
    b0 = ()
    b1 = ()
    for i, j in zip( POPT_Ires, POPT_beta ):
        p0.append( np.array( ( i[0], i[1], j[1], i[2], j[2], 0, 0 ) ) )
        b0 += ( 0, 0, 0, 0, 0, -1, 0 )
        b1 += ( np.inf, 100, np.inf, 10, np.inf, 1, np.inf)
    
    return nx, ny, H, np.asarray( p0 ), ( b0, b1 )

def cluster_me2dNDF( X, n_bins=500 ):
    nx, ny, H, p0, b = cluster_xdrNDF( X, n_bins=n_bins )
    popt, pcov = curve_fit( me2dNDF, ( nx, ny ), H.ravel(), p0=p0.ravel(), bounds=b )
    n_cols = 7
    n_rows = int( len( popt ) / 7 )
    P = popt.reshape( n_rows, n_cols )
    return nx, ny, H, P


















