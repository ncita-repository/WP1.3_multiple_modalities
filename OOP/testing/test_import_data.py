# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:38:46 2021

@author: ctorti
"""

from importlib import reload

#import testing.test_fetch_config
#reload(testing.test_fetch_config)

import testing.test_download_data
reload(testing.test_download_data)

import io_tools.import_data
reload(io_tools.import_data)

#from testing.test_fetch_config import ConfigFetcherTester
from testing.test_download_data import DataDownloaderTester
from io_tools.import_data import DataImporter


class DataImporterTester:
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Fetcher Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    runID : str
        The ID that determines the main configuration parameters to use for the
        run (i.e. the key in the base level dictionary contained in 
        main_params.json).
    params : ConfigFetcher Object, optional
        A ConfigFetcher object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_import_data import DataImporterTester
        
        testResult = ImporterTester(runID, params)
    
    To integrated test:
        from testing.test_import_data import DataImporterTester
        
        testResult = ImporterTester(runID)
    """
    
    def __init__(self, cfgDir, runID, params=None):
        
        if params == None:
            # Run DataDownloaderTester to get the config parameters and 
            # download the data from XNAT:
            testResult = DataDownloaderTester(cfgDir=cfgDir, runID=runID)
            params = testResult.params
            
            # Or run the internal operations within DataDownloaderTester:
            """
            # Run ConfigFetcherTester:
            cfgObj = ConfigFetcherTester(cfgDir=cfgDir, runID=runID)
            
            # Instantiate a DataDownloader object:
            params = DataDownloader(cfgDict=cfgObj.cfgDict)
            """
        
        # Initialise a DataImporter objects for source and target:
        #self.DataImporter(params=params, srcORtrg='src')
        srcDataset = DataImporter(params=params, srcORtrg='src')
        #self.DataImporter(params=params, srcORtrg='trg')
        trgDataset = DataImporter(params=params, srcORtrg='trg')
        
        # Import source and target data:
        print('Importing source dataset...\n')
        srcDataset.import_data(params)
        print('Importing target dataset...\n')
        trgDataset.import_data(params)
        
        self.params = params
        self.srcDataset = srcDataset
        self.trgDataset = trgDataset

def test_import_data_210812(runID, params=None):
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Fetcher Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(runID, params)
    
    To integrated test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(runID)
    """
        
    if params == None:
        # Fetch the parameters:
        print('Fetching parameters (locally) and data (from XNAT)...\n')
        params = test_fetch_params_and_data(runID)
    
    # Import data for source and target:
    print('Importing source dataset...\n')
    srcDataset = Importer(params=params, srcORtrg='src')
    
    print('Importing target dataset...\n')
    trgDataset = Importer(params=params, srcORtrg='trg')
    
    return srcDataset, trgDataset

def test_import_data_210809(params=None):
    """
    Run unit or integrated test of import_data.
    
    To run unit test, provide Params Object. To run integrated test, leave 
    input argument blank.
    
    Parameters
    ----------
    params : Params Object, optional
        A Params object. The default is None.
        
    Returns
    -------
    srcDataset : Dataset Object
        Dataset Object for Source DICOM series.
    trgDataset : Dataset Object
        Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data(Params)
    
    To integrated test of import_data:
        from testing.test_import_data import test_import_data
        
        srcDataset, trgDataset = test_import_data()
    """
    
    if params == None:
        # Fetch the parameters:
        params = test_fetch_params_and_data()
    
    # Import data for source and target:
    srcDataset = Dataset(xnatParams=params.srcXnatParams,
                         genParams=params.genParams)
    
    trgDataset = Dataset(xnatParams=params.trgXnatParams,
                         genParams=params.genParams)
    
    return srcDataset, trgDataset