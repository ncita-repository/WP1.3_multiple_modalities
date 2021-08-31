# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:19:43 2021

@author: ctorti
"""

from importlib import reload
import testing.test_fetch_config
reload(testing.test_fetch_config)
import io_tools.download_data
reload(io_tools.download_data)

from testing.test_fetch_config import ConfigFetcherTester
from io_tools.download_data import DataDownloader

class DataDownloaderTester:
    """
    Test DataDownloader in download_data.py
    
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
    
    def __init__(self, cfgDir, runID):
        # Run ConfigFetcherTester:
        cfgObj = ConfigFetcherTester(cfgDir=cfgDir, runID=runID)
        
        # Instantiate a DataDownloader object:
        params = DataDownloader(cfgDict=cfgObj.cfgDict)
        
        # Download data and create pathsDict:
        print('Downloading data from XNAT...\n')
        params.download_and_get_pathsDict()
        
        #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
        #params = Fetcher(cfgDir, params.xnatSession)
        
        self.cfgDict = params.cfgDict
        self.pathsDict = params.pathsDict
        self.xnatSession = params.xnatSession
        self.params = params
        
    def test_fetch_params_and_data_210812(runID):
        """
        Test fetch_params_and_data.
        
        Parameters
        ----------
        runID : str
            The ID that determines the main configuration parameters to use for the
            run (i.e. the key in the base level dictionary contained in 
            main_params.json).
        
        Returns
        -------
        Fetcher Object
            A Fetcher Object containing various parameters for the specified runID.
        
        Notes
        -----
        To run unit test:
            from testing.test_fetch_params_and_data import test_fetch_params_and_data
            
            params = test_fetch_params_and_data(runID)
        """
        
        cfgDir = r'C:\Code\WP1.3_multiple_modalities\Python_code\configs'
        
        print('Fetching parameters (locally) and data (from XNAT) without pre-',
              'existing XNAT session...\n')
        params = Fetcher(cfgDir, runID)
        
        #print('\n\n\nFetching parameters and data using pre-existing XNAT session...\n')
        #params = Fetcher(cfgDir, params.xnatSession)
        
        return params