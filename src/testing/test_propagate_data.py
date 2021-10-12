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

import io_tools.propagate
reload(io_tools.propagate)

import dro_tools.create_dro
reload(dro_tools.create_dro)

import dicom_tools.create_seg
reload(dicom_tools.create_seg)

import plotting_tools.general
reload(plotting_tools.general)

import plotting_tools.res_reg_results
reload(plotting_tools.res_reg_results)

import io_tools.exports
reload(io_tools.exports)

import os
import time
from pathlib import Path
import SimpleITK as sitk
from testing.test_import_data import DataImporterTester
from io_tools.import_dro import DroImporter
from io_tools.propagate import Propagator
from dro_tools.create_dro import DroCreator
from dicom_tools.create_seg import create_seg
#from plotting_tools.plotting import plot_pixarrs_from_list_of_segs_v3
from plotting_tools.general import (
    plot_pixarrBySeg, plot_pixarrs_from_list_of_segs_and_images
    )
from plotting_tools.res_reg_results import(
    plot_metricValues_v_iters, compare_res_results 
    )
from io_tools.exports import (export_im, export_list_to_txt, export_newRoicol)

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
        """ 
        useCaseToApply isn't an attribute until the Propagator class is
        implemented
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        """
        
        
        #if '5' in params.cfgDict['useCaseToApply']:
        #    """
        #    Image registration will be required unless a suitable DRO can be
        #    found on XNAT.  Search XNAT for an appropriate DRO:
        #    """
        # The following was under the if statement:
        
        if params.cfgDict['useDroForTx']:
            # Instantiate a droData Object:
            droData = DroImporter()
            droData.fetch_dro(params)
            self.droData = droData
            dro = droData.dro
            #print(f'type(droData.dro) = {type(droData.dro)}\n')
        else:
            dro = None
        
        #print(f'dir(srcDataset) = {dir(srcDataset)}\n')
        
        # Instantiate a Propagator object:
        newDataset = Propagator(srcDataset, trgDataset, params)
        
        # Propagate/copy the source ROI Collection:
        newDataset.execute(srcDataset, trgDataset, params, dro)
        
        self.newDataset = newDataset
        # Update params:
        params = self.params
        
        useCaseToApply = params.cfgDict['useCaseToApply']
        roicolMod = params.cfgDict['roicolMod']
        
        """ Creating a new ROI Collection and DRO (if applicable) is not
        required but adding it to more fully test (23/08/2021). """
        newRoicol = create_seg(srcDataset, trgDataset, newDataset, params)
        
        self.newDataset.roicol = newRoicol
        
        """
        Export the new ROI Collection - this is now done in a separate
        class - see dicom_tools.create_roicol.py
        """
        #self.export_newRoiCol()
        
        
        """ Create and export a new DRO (if applicable) - this is done in a 
        separate class so test outside of this Propagator-specific test...
        
        if '5' in useCaseToApply:
            print('\n*** Commenting out create & export DRO from',
                  'test_propagate_data.py\n')
            
            # Instantiate a DroCreator Object:
            newDroData = DroCreator(newDataset, params)
            
            # Create the new DRO:
            newDroData.create_dro(srcDataset, trgDataset, newDataset, params)
            
            # Export the new DRO:
            newDroData.export_dro(params)
            
            self.newDroData = newDroData
        """
    
        """ Export label images, tranforms and transform parameters """
        #self.export_labims()
        #self.export_tx_and_params()
        
        """ Plot results """
        #self.plot_metric_v_iters()
        #self.plot_res_results()
        #self.plot_roi_over_dicom_ims()