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

import io_tools.import_dro
reload(io_tools.import_dro)

import image_tools.propagate
reload(image_tools.propagate)

import dro_tools.create_dro
reload(dro_tools.create_dro)

import seg_tools.create_seg
reload(seg_tools.create_seg)

import plotting_tools.plotting
reload(plotting_tools.plotting)

from testing.test_import_data import DataImporterTester
from io_tools.import_dro import DroImporter
from image_tools.propagate import Propagator
from dro_tools.create_dro import DroCreator
from seg_tools.create_seg import create_seg

from plotting_tools.plotting import plot_pixarrs_from_list_of_segs_v3

class PropagatorTester:
    # TODO update docstrings
    """
    Run unit or integrated test of propagate.
    
    To run unit test, provide Importer Objects for source and target. To run 
    integrated test, leave optional input arguments blank.
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
        
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
    
    def __init__(
            self, cfgDir, runID, srcDataset=None, trgDataset=None, params=None
            ):
        if (params == None or srcDataset == None or trgDataset == None):
            # Import data for source:
            print('Running ImporterTester...\n')
            testResult = DataImporterTester(cfgDir=cfgDir, runID=runID)
            
            params = testResult.params
            srcDataset = testResult.srcDataset
            trgDataset = testResult.trgDataset
        
        self.params = params
        self.srcDataset = srcDataset
        self.trgDataset = trgDataset
        useCaseToApply = params.cfgDict['useCaseToApply']
        
        
        #if '5' in params.cfgDict['useCaseToApply']:
        #    """
        #    Image registration will be required unless a suitable DRO can be
        #    found on XNAT.  Search XNAT for an appropriate DRO:
        #    """
        # The following was under the if statement:
        # Instantiate a droData Object:
        droData = DroImporter()
        droData.fetch_dro(params)
        self.droData = droData
        dro = droData.dro
        #print(f'type(droData.dro) = {type(droData.dro)}\n')
        
        #print(f'dir(srcDataset) = {dir(srcDataset)}\n')
        
        # Instantiate a Propagator object:
        newDataset = Propagator(srcDataset, trgDataset, params)
        
        # Propagate/copy the source ROI Collection:
        newDataset.propagate(srcDataset, trgDataset, params, dro)
        
        self.newDataset = newDataset
        
        """ Creating a new ROI Collection and DRO (if applicable) is not
        required but adding it to more fully test (23/08/2021). """
        newRoicol = create_seg(srcDataset, trgDataset, newDataset, params)
        
        self.newDataset.roicol = newRoicol
        
        """ Note create_seg will need to be wrapped in another function along
        with the yet-to-be-refactored create_rts. """
        # TODO refactor create_rts and wrap it and create_seg into a new func (create_new_roicol)
        
        if '5' in useCaseToApply:
            # Instantiate a DroCreator Object:
            newDroData = DroCreator(newDataset, params)
            
            # Create the new DRO:
            newDroData.create_dro(srcDataset, trgDataset, newDataset, params)
        
            self.newDroData = newDroData
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        p2c = cfgDict['p2c']
        useCaseToApply = cfgDict['useCaseToApply']
        forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        
        if roicolMod == 'RTSTRUCT':
            exportDir = cfgDict['rtsPlotsExportDir']
        else:
            exportDir = cfgDict['segPlotsExportDir']
        
        listOfSegs = [srcDataset.roicol, newRoicol]
        listOfDcmDirs = [srcDataset.dicomDir, newDataset.dicomDir]
        listOfPlotTitles = ['Source', 'New Target']
        
        plot_pixarrs_from_list_of_segs_v3(
        listOfSegs, listOfDcmDirs, listOfPlotTitles, exportPlot=exportPlot,
            exportDir=exportDir, runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, p2c=p2c
        )

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