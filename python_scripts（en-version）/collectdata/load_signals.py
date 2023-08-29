import logging
import numpy as np
import nanolyse as nl
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_signals( self, i ):
    BLANK_signal = np.array([]); DATA_signal = np.array([]);
    for blank, fname in zip( self.index[self.index_cols[0]], self.index[self.index_cols[1]] ):
        try:
            fname_load = self.notebook_path + '/data/' + self.main_index[self.main_index_cols[1]][i] + fname + '.abf'
            signal, self.sampling_period = nl.loadAxon( self.notebook_path + '/data/' + self.main_index[self.main_index_cols[1]][i] + fname + '.abf' )
            for s in signal:
                if blank==True:
                    BLANK_signal = np.concatenate( ( BLANK_signal, s[int(self.offset/self.sampling_period)::] ) )
                else:
                    DATA_signal = np.concatenate( ( DATA_signal, s[int(self.offset/self.sampling_period)::] ) )
            logger.debug('Succes: ' + fname_load )
        except:
            try:
                logger.error('Unable to process data: ' + fname_load )
            except:
                logger.error('UNKNOWN ERROR WHILE LOADING DATA')
    return BLANK_signal, DATA_signal