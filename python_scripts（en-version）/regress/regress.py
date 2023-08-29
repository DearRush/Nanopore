# -*- coding: utf-8 -*-
import numpy as np
class regress:
    def __init__( self, x, y, order, ex_0=False ):
        
        self.ex_0 = ex_0
        x = np.asarray( x, dtype=np.float64 )[::np.newaxis]
        y = np.asarray( y, dtype=np.float64 )[::np.newaxis]
        if self.ex_0==False:
            X = np.expand_dims( x ** 0, axis=1)
        for i in range( order ):
            try:
                X = np.append( X, np.expand_dims( x ** ( i + 1 ), axis=1), axis=1 )
            except:
                X = np.expand_dims( x ** ( i + 1 ), axis=1)
        self.b = np.dot( np.dot( np.linalg.inv( np.dot( np.transpose( X ), X ) ), np.transpose( X ) ), y )
        self.H = np.dot( np.dot( X, np.linalg.inv( np.dot( np.transpose( X ), X ) ) ), np.transpose( X ) )
        
        # Calculate the errors
        Y = np.array([ i for i in y ])
        e_hat = Y - self.predict( x )
        self.SStot = np.dot( np.transpose( Y - np.mean( Y ) ), Y - np.mean( Y ) )
        self.SSres = np.dot( np.transpose( e_hat ), e_hat ) / ( len( x ) - order )
        self.Rsquared  = 1 - ( self.SSres / self.SStot )
        DFe = len( x ) - len( self.b ) - 1
        DFt = len( x ) - 1
        self.AdjustRsquared = 1 - ( ( self.SSres / DFe ) / ( self.SStot / DFt ) )
        self.cov_b = np.diag( np.dot( self.SSres, np.linalg.inv( np.dot( np.transpose( X ), X ) ) ) )
        self.sigma_b = np.sqrt( self.cov_b )
        
    def predict( self, x ):
        x = np.asarray( x, dtype=np.float64 )[::np.newaxis]
        l = len( self.b )
        if self.ex_0==False:
            X = np.expand_dims( x ** 0, axis=1)
            l -= 1
        for i in range( l ):
            try:
                X = np.append( X, np.expand_dims( x ** ( i + 1 ), axis=1), axis=1 )
            except:
                X = np.expand_dims( x ** ( i + 1 ), axis=1)
        
        Y = np.dot( X, self.b )
        return Y