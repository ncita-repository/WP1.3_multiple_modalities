# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:19:43 2021

@author: ctorti
"""

from importlib import reload
import io_tools.fetch_config
reload(io_tools.fetch_config)
from io_tools.fetch_config import ConfigFetcher

class ConfigFetcherTester:
    """
    Test ConfigFetcher in fetch_config.py
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    
    Returns
    -------
    ConfigFetcherTester Object
        Containing various parameters for the specified runID.
    
    Notes
    -----
    To run unit test:
        from testing.test_fetch_config import ConfigFetcherTester
        
        testResult = ConfigFetcherTester(runID)
    """
    
    def __init__(self, cfgDir, runID):
        
        #print('Fetching parameters from local drive...\n')
        
        # Instanstantiate a ConfigFetcher object:
        cfgObj = ConfigFetcher(cfgDir=cfgDir, runID=runID)
        
        #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
        #params = Fetcher(cfgDir, params.xnatSession)
        
        # Get the config settings:
        cfgObj.get_config()
        
        self.cfgDir = cfgObj.cfgDir
        self.cfgDict = cfgObj.cfgDict