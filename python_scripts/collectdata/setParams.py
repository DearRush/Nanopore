import logging
import numpy as np
import nanolyse as nl
import pandas as pd

def setParams( self, notebook_path, results_cols ):
    self.notebook_path = notebook_path
    self.results = pd.DataFrame( columns=np.array(results_cols) )
    self.results_cols = [ i for i in self.results.columns ]