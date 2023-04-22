# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit
import sys
notebook_path="D:/prp"
sys.path.insert(0, notebook_path + "/python_scripts")
import nanolyse as nl
import os
from NDFs import sNDF

from ipywidgets import IntProgress
from IPython.display import display
import time


def get_events( fname, levels ):
    signal, sampling_period = nl.loadAxon( fname )
    signal_no_filt, _ = nl.loadAxon( fname )
    Fs = int( 1 / sampling_period )
    if len( signal ) > 1:
        N_offset = 1.8
    else:
        N_offset = 0
    for c, s in enumerate( signal ):
        signal[c] = nl.filter_gaussian( signal[c][int(Fs*N_offset)::], sampling_period, Fs )
        l0, l1, L0_start, L0_end, L1_start, L1_end = nl.thresholdsearch( signal, sampling_period, levels, trace=c )
        l0 = []
        l1 = []
        # Reconconstruct the original events from the signal, to not lose information
        for l0_start, l0_end, l1_start, l1_end in zip( L0_start, L0_end, L1_start, L1_end ):
            l0.append( signal_no_filt[c][int(Fs*N_offset)+l0_start:int(Fs*N_offset)+l0_end] )
            l1.append( signal_no_filt[c][int(Fs*N_offset)+l1_start:int(Fs*N_offset)+l1_end] )
        l0 = np.array( l0 ); l1 = np.array( l1 );
        if c==0:
            L0 = l0
            L1 = l1
        else:
            try:
                L0 = np.concatenate( ( L0, l0 ) )
                L1 = np.concatenate( ( L1, l1 ) )
            except:
                pass
    return L0, L1, sampling_period

def events_sNDF( fname, levels, load_data=True ):
    f = IntProgress(min=0, max=100, description='Fitting parameters:')
    display(f)
    print( fname )
    
    K = False
    if load_data==True:
        try:
            super_fit = list( np.load( os.path.dirname( fname ) + "/super_fit.npy" ) )
            super_fit_cov = list( np.load( os.path.dirname( fname ) + "/super_fit_cov.npy" ) )
            SD_2 = list( np.load( os.path.dirname( fname ) + "/SD_2.npy" ) )
            K = True
        except:
            pass
    if K==False:
        # Get the estimated event locations using RED
        L0, L1, sampling_period = get_events( fname, levels )
        print( len( L1 ) )
        super_fit = []
        super_fit_cov = []
        SD_2 = []
        for i in range( len( L1 ) - 1 ):
            # Safety function, if something crashes (e.g. no good fit) it will go to the except statement.
            try:
                # Combine the event back with it's surrounding baseline
                Y = np.concatenate( ( L0[i], L1[i], L0[i+1] ) )
                
                # Get some range of x
                x = np.linspace( 0, len( Y )*sampling_period, len( Y ) )
                
                # Estimate location of x0
                x0 = ( len( L0[i] ) + len( L1[i] )/2 ) * sampling_period
                
                # Estimate the sigma
                sigma_x0 = ( len( L1[i] )/2 ) * sampling_period
                
                # Parameters [amplitude, x0, sigma_x, Beta, offset]
                # Beta=2 -> Normal distribution
                p0 = ( abs( np.mean(L1[i])-np.mean( np.concatenate( ( L0[i], L0[i+1] ) ) ) ), 
                      x0, sigma_x0, 2, np.mean( np.concatenate( ( L0[i], L0[i+1] ) ) ) )
                
                # Fit a super normal distribution to the event
                popt, pcov = curve_fit( sNDF, x, -1*abs( Y ), p0=p0 )
                
                SD_2.append( np.average( ( Y - sNDF( x, *popt ) ) ** 2, weights=sNDF( x, *popt ) ) )
                
                # Append the resulting parameters
                super_fit.append( popt )
                super_fit_cov.append( pcov )
                
                f.value = int( 100 * ( i / len( L1 ) ) )
            except:
                # Some times fitting doesn't work, just ignore the "event" and move on.
                #print( 'Fitting error on event: ' + str( i ) )
                pass
        try:
            np.save( os.path.dirname( fname ) + "/super_fit", np.array( super_fit ) )
            np.save( os.path.dirname( fname ) + "/super_fit_cov", np.array( super_fit_cov ) )
            np.save( os.path.dirname( fname ) + "/SD_2", np.array( SD_2 ) )
        except:
            #print( 'Error saving: ' + str( fname ) )
            pass
    f.value = 100
    return super_fit, SD_2


def features_sNDF( super_fit ):
    #SF = [ i[0] for i in super_fit ]
    SF = super_fit
    Iex = np.array([ i[0]/abs(i[-1]) for i in SF ])
    DWT = np.array([ 2 * i[2] * np.sqrt( ( 2 * np.log(2)**(1/i[3]) ) ) for i in SF ])
    beta = np.array([ i[3] for i in SF ])
    return Iex, DWT, beta