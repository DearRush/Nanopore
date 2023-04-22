import numpy as np
import pandas as pd
import nanolyse as nl

def collect_data( notebook_path, main_index, offset=1.8, Fs=5000, dwelltime=4e-4, sNDF=False ):
    try:
        unpickled_results = pd.read_pickle(notebook_path + "/data/results.pkl" )
        unpickled_results_cols = [ i for i in unpickled_results.columns ]
        analysed_folders = [ unpickled_results[unpickled_results_cols[0]][i] for i in range( len( unpickled_results ) ) ]
    except:
        pass
    main_index_cols = [ i for i in main_index.columns ]
    results = pd.DataFrame( columns=np.array(['Folder',
                                              'Protein', 
                                              'Blank Signal', 
                                              'Signal',
                                              'Iex BLANK', 
                                              'Iex SD2 BLANK',
                                              'dwelltime BLANK',
                                              'beta BLANK',
                                              'fit params BLANK', 
                                              'Iex', 
                                              'Iex SD2',
                                              'dwelltime', 
                                              'beta',
                                              'fit params']) )
    results_cols = [ i for i in results.columns ]
    K = False
    for i in range( len( main_index ) ):
        try:
            if main_index[main_index_cols[1]][i] in analysed_folders:
                K = True
        except:
            pass
        if K == False:
            idx = notebook_path + '/data/' + main_index[main_index_cols[1]][i] + '/index.csv'
            I0 = main_index[main_index_cols[2]][i]
            SD = main_index[main_index_cols[3]][i]
            levels = ( False, I0, SD*3 )
            index = pd.read_csv( idx )
            index_cols = [ i for i in index.columns ]
            BLANK_signal = np.array([])
            DATA_signal = np.array([])
            for blank, fname in zip( index[index_cols[0]], index[index_cols[1]] ):
                signal, sampling_period = nl.loadAxon( notebook_path + '/data/' + main_index[main_index_cols[1]][i] + fname + '.abf' )
                for s in signal:
                    if blank==True:
                        BLANK_signal = np.concatenate( ( BLANK_signal, s[int(offset/sampling_period)::] ) )
                    else:
                        DATA_signal = np.concatenate( ( DATA_signal, s[int(offset/sampling_period)::] ) )
            
            BLANK_signal = nl.filter_gaussian( BLANK_signal, sampling_period, Fs )
            DATA_signal = nl.filter_gaussian( DATA_signal, sampling_period, Fs )
            BLANK_events_THS = nl.thresholdsearch( [BLANK_signal], sampling_period, levels, dwelltime=dwelltime )
            DATA_events_THS = nl.thresholdsearch( [DATA_signal], sampling_period, levels, dwelltime=dwelltime )
            DATA_beta = False
            BLANK_beta = False
            if sNDF==True: # Slow precise method
                BLANK_fit, BLANK_fit_cov, BLANK_Iex_SD_2 = nl.fit_sNDF( BLANK_events_THS, sampling_period )
                DATA_fit, DATA_fit_cov, DATA_Iex_SD_2 = nl.fit_sNDF( DATA_events_THS, sampling_period )
                BLANK_Iex, BLANK_dwelltime, BLANK_beta = nl.features_sNDF( BLANK_fit )
                DATA_Iex, DATA_dwelltime, DATA_beta = nl.features_sNDF( DATA_fit )
                BLANK_beta = np.sqrt( np.array( [ j[3] for j in BLANK_fit ] ) )
                DATA_beta = np.sqrt( np.array( [ j[3] for j in DATA_fit ] ) )
            else: # Quick method
                BLANK_Iex, BLANK_Iex_SD_2, BLANK_dwelltime = nl.get_features_THS( BLANK_events_THS, sampling_period )
                DATA_Iex, DATA_Iex_SD_2, DATA_dwelltime = nl.get_features_THS( DATA_events_THS, sampling_period )
                BLANK_beta = np.zeros( len( BLANK_Iex ) )
                DATA_beta = np.zeros( len( DATA_Iex ) )
                BLANK_fit = np.array([levels, sampling_period])
                DATA_fit = np.array([levels, sampling_period])
            df = pd.DataFrame({'Folder':main_index[main_index_cols[1]][i],
                               'Protein':main_index[main_index_cols[0]][i], 
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
            results = results.append(df, ignore_index = True)
        K = False
    try:
        results = pd.concat( [ unpickled_results, results ] )
    except:
        pass
    results.to_pickle( notebook_path + "/data/results.pkl" )
    return results