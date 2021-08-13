# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:52:32 2021

@author: ctorti
"""


from importlib import reload

#import testing.test_import_data
#reload(testing.test_import_data)
#from testing.test_import_data import FetcherTester

import testing.test_import_data
reload(testing.test_import_data)
from testing.test_import_data import ImporterTester

import propagate_tools.propagate
reload(propagate_tools.propagate)
from propagate_tools.propagate import Propagator

class PropagatorTester:
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object containing various parameters. The default is None.
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
        
    Returns
    -------
    newDataset : Dataset Object
        New Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_propagate_data import PropagatorTester
        
        testResult = PropagatorTester(
            runID, params, srcDataset, trgDataset
            )
    
    To integrated test:
        from testing.test_propagate_data import PropagatorTester
        
        testResult = PropagatorTester(runID)
    """
    
    def __init__(self, runID, params=None, srcDataset=None, trgDataset=None):
        if (params == None or srcDataset == None or trgDataset == None):
            # Import data for source:
            print('Running ImporterTester...\n')
            testResult = ImporterTester(runID=runID)
            
            params = testResult.params
            srcDataset = testResult.srcDataset
            trgDataset = testResult.trgDataset
        
        self.params = params
        self.srcDataset = srcDataset
        self.trgDataset = trgDataset
        
        newDataset = Propagator(params, srcDataset, trgDataset)
        
        self.newDataset = newDataset

def test_propagator_210812(runID, params=None, srcDataset=None, trgDataset=None):
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    params : Fetcher Object, optional
        A Fetcher object containing various parameters. The default is None.
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
        
    Returns
    -------
    newDataset : Dataset Object
        New Dataset Object for Target DICOM series.
    
    Notes
    -----
    To run unit test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(
            runID, params, srcDataset, trgDataset
            )
    
    To integrated test:
        from testing.test_propagate_data import test_propagator
        
        newDataset = test_propagator(runID)
    """
    
    if params == None:
        # Fetch the parameters:
        print('Fetching parameters (locally) and data (from XNAT)...\n')
        params = test_fetch_params_and_data(runID)
    
    if srcDataset == None:
        # Import data for source:
        print('Importing source dataset...\n')
        srcDataset = Importer(params=params, srcORtrg='src')
    
    if trgDataset == None:
        # Import data for target:
        print('Importing target dataset...\n')
        trgDataset = Importer(params=params, srcORtrg='trg')
    
    return srcDataset, trgDataset