# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:19:43 2021

@author: ctorti
"""

from importlib import reload
import io_tools.fetch_params_and_data
reload(io_tools.fetch_params_and_data)
from io_tools.fetch_params_and_data import Fetcher

class FetcherTester:
    """
    Test fetch_params_and_data.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    params : Fetcher Object, optional
        A Fetcher object. The default is None.
    
    Returns
    -------
    Fetcher Object
        A Fetcher Object containing various parameters for the specified runID.
    
    Notes
    -----
    To run unit test:
        from testing.test_fetch_params_and_data import FetcherTester
        
        testResult = FetcherTester(runID)
    """
    
    def __init__(self, runID, params=None):
        cfgDir = r'C:\Code\WP1.3_multiple_modalities\Python_code\configs'
        
        print('Fetching parameters (locally) and data (from XNAT) without pre-',
              'existing XNAT session...\n')
        params = Fetcher(cfgDir, runID)
        
        #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
        #params = Fetcher(cfgDir, params.xnatSession)
        
        self.params = params