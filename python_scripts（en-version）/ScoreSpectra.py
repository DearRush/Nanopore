import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.patches as patches

class ScoreSpectra:
    def __init__( self ):
        pass
    
    def _get_features( self, A ):
        dA = np.diff( A )
        mA = A - np.mean( A )
        mdA = dA - np.mean( dA )
        return dA, mA, mdA
    
    def correlate( self, A1, A2, method='Cor'):
        dA1, mA1, mdA1 = self._get_features( A1 )
        dA2, mA2, mdA2 = self._get_features( A2 )
        if method=='Cor':
            return ( np.sum( ( mA1 ) * ( mA2 ) )**2 ) / ( np.sum( mA1**2 ) * np.sum( mA2**2 ) )
        elif method=='DCor':
            return ( ( np.sum( ( mdA1 ) * ( mdA2 ) ) )**2 ) / ( np.sum( mdA1**2 ) * np.sum( mdA2**2 ) )
        elif method=='Euc':
            return ( np.sum( A1 * A2 )**2 ) / ( np.sum( A1**2 ) * np.sum( A2**2 ) )
        elif method=='DEuc':
            return ( np.sum( dA1 * dA2 )**2 ) / ( np.sum( dA1**2 ) * np.sum( dA2**2 ) )
    
    def leave_one_out( self, X, method='DEuc' ):
        Rs = np.zeros( ( len( X ), len( X ) ) )
        for c0, x in enumerate( X ):
            sample_scores = []
            for c1, sample in enumerate( x ):
                scores = []; Database = [];
                for c2, x0 in enumerate( X ):
                    if c2==c0: Database.append( np.mean( [ x0[c3] for c3 in range( len( x0 ) ) if c3!=c1 and c2==c0 ], axis=0 ) )
                    else: Database.append( np.mean( x0, axis=0 ) )
                for database in Database:
                    scores.append( self.correlate( sample, database, method=method ) )
                sample_scores.append( 100 * np.array( scores ) / sum( scores ) )
            for c1, j in enumerate( np.mean( sample_scores, axis=0 ) ):
                Rs[c0][c1] = j
        return Rs
    
    def score_spectra( self, X, method='DEuc' ):
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
            for c2, A2 in enumerate( X ):
                Rs[c1][c2] = self.correlate( A1, A2, method=method)
        return Rs
    
    def plot_Rs( self, protein_names, Rs, title=''):
        fig_size = (14, 5)
        fig = plt.figure( figsize=fig_size )
        gs = plt.GridSpec( 1, 2, top=0.9, wspace=0.4, hspace=0)
        ax1 = fig.add_subplot( gs[0,1] )
        im = ax1.imshow(Rs, cmap="OrRd")
        cbar = ax1.figure.colorbar(im, ax=ax1)
        ax1.set_xticks(np.arange(len(protein_names)))
        ax1.set_yticks(np.arange(len(protein_names)))
        ax1.set_xticklabels(protein_names, rotation='vertical')
        ax1.set_yticklabels(protein_names)
        ax1.set_title( title, fontsize=16)
        ax1.text(11, -1.2, r'Identical signal', fontsize=10, weight='bold', ha='left', va='center')
        ax1.text(11, 9.8, r'Dissimilar signal', fontsize=10, weight='bold', ha='left', va='center')
        ax1.text(12, 3.8, r'$P(X)$', fontsize=16, weight='bold', ha='center', va='center', rotation=270)

        rec = []
        for c, i in enumerate(Rs):
            idx = np.argmax(i)
            rec.append( patches.Rectangle((c-0.5, idx-0.5), 1, 1, linewidth=3, edgecolor='k', facecolor='none') )
            ax1.add_patch(rec[-1])

        return fig