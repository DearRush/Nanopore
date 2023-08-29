import logging
import numpy as np
import pandas as pd
import nanolyse as nl
notebook_path="D:/prp"
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def collect_data( self, main_index, offset=1.8, Fs=5000, dwelltime=4e-4, sNDF=False, reload=False ):
    self.main_index = main_index
    self.main_index_cols = [ i for i in self.main_index.columns ]
    self.offset = offset
    for i in range( len( self.main_index ) ):
        results_fname = notebook_path + '/data/' + self.main_index[self.main_index_cols[1]][i] + "/results.pkl"
        try:
            # In this try block, we try to load the results file generated.
            # An error will we raised if reload == False, or if the results don't excist yet.
            # Uppon error, the block will be skipped.
            if reload == False:
                raise Exception('Reloading data ' + self.main_index[self.main_index_cols[1]][i])
            df = pd.read_pickle( results_fname )
            logger.debug( 'Succes while loading: ' + results_fname )
        except:
            # The results data could not be loaded, so calcualte everything from scratch.
            # See the above block for more information.
            idx = notebook_path + '/data/' + self.main_index[self.main_index_cols[1]][i] + '/index.csv'
            I0 = self.main_index[self.main_index_cols[2]][i]
            SD = self.main_index[self.main_index_cols[3]][i]
            levels = ( False, I0, SD*3 )
            self.index = pd.read_csv( idx )
            self.index_cols = [ i for i in self.index.columns ]
            
            
            BLANK_signal, DATA_signal = self.load_signals( i )
            BLANK_signal = nl.filter_gaussian( BLANK_signal, self.sampling_period, Fs )
            DATA_signal = nl.filter_gaussian( DATA_signal, self.sampling_period, Fs )
            BLANK_events_THS = nl.thresholdsearch( [BLANK_signal], self.sampling_period, levels, dwelltime=dwelltime )
            DATA_events_THS = nl.thresholdsearch( [DATA_signal], self.sampling_period, levels, dwelltime=dwelltime )
            
            if sNDF==True: # Slow precise method
                BLANK_fit, BLANK_fit_cov, BLANK_Iex_SD_2 = nl.fit_sNDF( BLANK_events_THS, self.sampling_period )
                DATA_fit, DATA_fit_cov, DATA_Iex_SD_2 = nl.fit_sNDF( DATA_events_THS, self.sampling_period )
                BLANK_Iex, BLANK_dwelltime, BLANK_beta = nl.features_sNDF( BLANK_fit )
                DATA_Iex, DATA_dwelltime, DATA_beta = nl.features_sNDF( DATA_fit )
                BLANK_beta = np.sqrt( np.array( [ j[3] for j in BLANK_fit ] ) )
                DATA_beta = np.sqrt( np.array( [ j[3] for j in DATA_fit ] ) )
            else:
                BLANK_Iex, BLANK_Iex_SD_2, BLANK_dwelltime = nl.get_features_THS( BLANK_events_THS, self.sampling_period )
                DATA_Iex, DATA_Iex_SD_2, DATA_dwelltime = nl.get_features_THS( DATA_events_THS, self.sampling_period )
                DATA_beta = False
                BLANK_beta = False
                BLANK_fit = np.array([levels, self.sampling_period])
                DATA_fit = np.array([levels, self.sampling_period])
            
            df = pd.DataFrame({'Folder':self.main_index[self.main_index_cols[1]][i],
                               'Protein':self.main_index[self.main_index_cols[0]][i], 
                               'Blank Signal':[BLANK_signal], 
                               'Signal':[DATA_signal],
                               'Iex BLANK':[BLANK_Iex],
                               'Iex SD2 BLANK':[BLANK_Iex_SD_2], 
                               'dwelltime BLANK':[BLANK_dwelltime],
                               'beta BLANK':[BLANK_beta],
                               'fit params BLANK':[BLANK_fit],
                               'Iex':[DATA_Iex],
                               'Iex SD2':[DATA_Iex_SD_2], 
                               'dwelltime':[DATA_dwelltime], 
                               'beta':[DATA_beta], 
                               'fit params':[DATA_fit]}) 
            df.to_pickle( results_fname )