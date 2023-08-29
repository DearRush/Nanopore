import copy
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy import signal
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from IPython.display import Latex


class NanoMatch:
    def __init__( self, Iex_bins, Iex_edges, min_iex=0.4, max_iex=1, min_dwelltime=5e-4, min_beta=1, max_shift=0.05 ):
        self.Iex_bins = Iex_bins
        self.Iex_edges = Iex_edges
        self.min_iex = min_iex
        self.max_iex = max_iex
        self.min_dwelltime = min_dwelltime
        self.min_beta = min_beta
        self.max_shift = 0.05
    
    def _get_rows( self, iex, beta=None, dwelltime=None ):
        try: return np.where( (~np.isnan(iex)) & (beta>self.min_beta) & (iex<self.max_iex) & (iex>self.min_iex) & (dwelltime>self.min_dwelltime) )
        except: return np.where( (iex<self.max_iex) & (iex>self.min_iex) )
    
    def _get_spectrum( self, iex, beta=None, dwelltime=None ):
        try: iex_spectrum, _  = np.histogram( iex[ self._get_rows( iex, beta, dwelltime ) ], bins=self.Iex_bins )
        except: iex_spectrum, _  = np.histogram( iex, bins=self.Iex_bins )
        return iex_spectrum / sum( iex_spectrum )
    
    def _make_plot( self, protein_names, xlabel='', ylabel='' ):
        fig_size = (10, 7)
        fig = plt.figure( figsize=fig_size )
        gs = plt.GridSpec( 4, 3, top=1.0, wspace=0.5, hspace=0.7)
        
        ax = []
        for c, p in enumerate( protein_names ):
            col_number = int( (c) % 4 )
            row_number = int( ( c - col_number ) / 3 )
            ax.append( fig.add_subplot( gs[ col_number, row_number ] ) )
            ax[-1].set_title( protein_names[c], fontsize=12 )
            ax[-1].set_xlabel( xlabel, fontsize=10 )
            ax[-1].set_ylabel( ylabel, fontsize=10 )
        return fig, ax
    
    def allign( self, char_vars, shift_resolution=200, get_error=False ):
        shift_line = np.linspace( -self.max_shift, self.max_shift, shift_resolution )
        char_vars = list( char_vars )
        
        Iex_spectrum = []
        for iex, iex_SD, dwelltime, beta in char_vars:
            if len( iex )!=0:
                Iex_spectrum.append( self._get_spectrum( iex, beta, dwelltime) )
        
        no_spec = np.array( [ i for i in range( len( Iex_spectrum ) ) ] )
        min_shift = np.argmin([ np.sum( [ ( Iex_spectrum[i] - Iex_spectrum[j] ) ** 2 for j in np.where( no_spec!=i )[0] ] ) for i in no_spec ])
        
        Iex_cor = []
        Iex_SD_cor = []
        Dwelltime_cor = []
        Beta_cor = []
        Error_cor = []
        Shift = []
        for iex, iex_SD, dwelltime, beta in char_vars:
            if len( iex )!=0:
                error = [ np.sum( ( Iex_spectrum[min_shift] - self._get_spectrum( iex+j, beta, dwelltime) ) ** 2 ) for j in shift_line ]
                iex_shift = iex+shift_line[ np.argmin( error ) ]
                keep_rows = self._get_rows( iex_shift, beta, dwelltime )
                Iex_cor.append( iex_shift[ keep_rows ] )
                Iex_SD_cor.append( np.array( iex_SD )[ keep_rows ] )
                Dwelltime_cor.append( dwelltime[ keep_rows ] )
                Beta_cor.append( beta[ keep_rows ] )
                Error_cor.append( error )
                Shift.append( shift_line )
        if get_error==False:
            return Iex_cor, Iex_SD_cor, Dwelltime_cor, Beta_cor
        else:
            return Iex_cor, Iex_SD_cor, Dwelltime_cor, Beta_cor, Error_cor, Shift
    
    def get_features( self, char_vars, allign=False ):
        Iex = []; IexSD2 = []; Dwelltime = []; Beta = [];
        for iex, iex_SD, dwelltime, beta in char_vars:
            if len( iex )!=0:
                keep_rows = self._get_rows( iex, beta, dwelltime )
                Iex.append( iex[keep_rows] )
                IexSD2.append( np.array( iex_SD )[keep_rows] )
                Dwelltime.append( dwelltime[keep_rows] )
                Beta.append( beta[keep_rows] )
        if allign==True:
            Iex, IexSD2, Dwelltime, Beta = self.allign( char_vars )
        return Iex, IexSD2, Dwelltime, Beta
    
    def get_spectra_all( self, results, allign=True ):
        protein_names = []
        iex_spectra = []
        for c, p in enumerate( np.unique( results.loc[:,'Protein'] ) ):
            data = results.where( results.loc[:,'Protein'] == p ).dropna()
            char_vars = list( zip( data['Iex'], data['Iex SD2'], data['dwelltime'], data['beta'] ) )
            Iex, IexSD2, Dwelltime, Beta = self.get_features( char_vars, allign=allign )
            Iex_spectrum = []
            for iex in Iex:
                Iex_spectrum.append( self._get_spectrum( iex ) )
            iex_spectra.append( Iex_spectrum )
            protein_names.append( p )
        return protein_names, iex_spectra
    
    def plot_shift_error( self, results ):
        protein_names, iex_spectra = self.get_spectra_all( results )
        fig, ax = self._make_plot( protein_names, xlabel=r'$\Delta$I$_{ex}$\%', ylabel=r'Error (arb. units)' )
        for c, p in enumerate( np.unique( results.loc[:,'Protein'] ) ):
            data = results.where( results.loc[:,'Protein'] == p ).dropna()
            char_vars = list( zip( data['Iex'], data['Iex SD2'], data['dwelltime'], data['beta'] ) )
            _, _, _, _, Error, Shift = self.allign( char_vars, get_error=True )
            for shift, error in zip( Shift, Error ):
                ax[c].plot( shift*100, error )
            ax[c].set_xlim([ -self.max_shift*100, self.max_shift*100 ])
            #ax[c].set_ylim([ 5e-4, 1e-1 ])
        return fig
    
    def plot_Iex_dwelltime( self, results, allign=True ):
        protein_names, iex_spectra = self.get_spectra_all( results )
        fig, ax = self._make_plot( protein_names, xlabel=r'I$_{ex}$\%', ylabel=r'Dwelltime (s)' )
        for c, p in enumerate( np.unique( results.loc[:,'Protein'] ) ):
            data = results.where( results.loc[:,'Protein'] == p ).dropna()
            char_vars = list( zip( data['Iex'], data['Iex SD2'], data['dwelltime'], data['beta'] ) )
            Iex, IexSD2, Dwelltime, Beta = self.get_features( char_vars, allign=allign )
            for iex, dwelltime in zip( Iex, Dwelltime ):
                ax[c].semilogy( iex*100, dwelltime, '.k', markersize=0.05 )
            ax[c].set_xlim([ 40, 100 ])
            ax[c].set_ylim([ 5e-4, 1e-1 ])
        return fig
    
    def plot_Iex_spectra( self, results, allign=True ):
        protein_names, iex_spectra = self.get_spectra_all( results, allign )
        fig, ax = self._make_plot( protein_names, xlabel=r'I$_{ex}$\%', ylabel=r'P(X) \%' )
        for c, X in enumerate( iex_spectra ):
            ax[c].plot( self.Iex_edges*100, 100 * np.mean( X, axis=0 ), 'k' )
            ax[c].plot( self.Iex_edges*100, 100 * ( np.mean( X, axis=0 ) + np.std( X, axis=0 ) ), ':r' )
            ax[c].plot( self.Iex_edges*100, 100 * ( np.mean( X, axis=0 ) - np.std( X, axis=0 ) ), ':r' )
            ax[c].set_xlim([40,100])
            ax[c].set_ylim([0,max(130 * np.mean( X, axis=0 ))])
        return fig
    
