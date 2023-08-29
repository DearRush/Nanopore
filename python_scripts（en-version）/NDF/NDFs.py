# -*- coding: utf-8 -*-
from scipy import asarray as ar, exp
import numpy as np

# rNDF: regular Normal Distribution Function
def rNDF( x, a, x0, sigma ):
    return a*exp(-(x-x0)**2/(2*sigma**2))

# sNDF: super Normal Distribution Function
def sNDF( x, A, x0, sigma, B, C ):
    E = -( ( x - x0 )**2 / ( 2 * sigma**2 ) )**B
    return ( A*exp( E ) ) + C

# e2dNDF: elliptical 2-dimensional Normal Distribution Function
def e2dNDF( X, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = np.meshgrid(X[0], X[1])
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

# me2dNDF: multiple elliptical 2-dimensional Normal Distribution Function
def me2dNDF( X, *P):
    x, y = np.meshgrid(X[0], X[1])
    n_cols = 7
    n_rows = int( len( P ) / 7 )
    Hx = np.zeros( ( len( x ), len( y ) ) )
    for p in np.asarray( P ).reshape( n_rows, n_cols ):
        Hx += e2dNDF( X, *p ).reshape(len(x), len(y))
    return Hx.ravel()

