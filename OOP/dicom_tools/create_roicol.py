# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 19:07:06 2021

@author: ctorti
"""

from importlib import reload
import dro_tools.create_spa_dro
import dro_tools.create_def_dro
import dro_tools.matrices
import general_tools.general
reload(general_tools.general)
reload(dro_tools.create_spa_dro)
reload(dro_tools.create_def_dro)
reload(dro_tools.matrices)
import dicom_tools.create_rts
reload(dicom_tools.create_rts)
import dicom_tools.create_seg
reload(dicom_tools.create_seg)
import dicom_tools.error_check_rts
reload(dicom_tools.error_check_rts)
import dicom_tools.error_check_seg
reload(dicom_tools.error_check_seg)

import os
from pathlib import Path
import time
#from datetime import datetime
#from pydicom import dcmread
#from pydicom.uid import generate_uid
#import SimpleITK as sitk
#import numpy as np
from dicom_tools.create_seg import create_seg
from dicom_tools.create_rts import create_rts
from dicom_tools.error_check_rts import error_check_rts
from dicom_tools.error_check_seg import error_check_seg

#from general_tools.general import (
#    reduce_list_of_str_floats_to_16, generate_reg_fname
#    )


class RoicolCreator:
    # TODO update docstrings
    """
    This class creates a RoicolCreator Object to create a DICOM-RTSTRUCT (RTS)
    or DICOM-SEG (SEG).
    
    Parameters
    ----------
    srcDataset : Importer Object
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    self.roicol : Pydicom Object
        The new RTS or SEG object.
        
    Note
    ----
    
    """
    
    """
    def __init__(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        
        Initialise the object.
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        trgDataset : Importer Object
            Object containing various data for the target DICOM series.
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.srcDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        self.trgDataset : Importer Object
            Object containing various data for the target DICOM series.
        self.newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        self.params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
    """
    
    def create_roicol(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        """
        Create a new ROI Collection (RTS/SEG).
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        trgDataset : Importer Object
            Object containing various data for the target DICOM series.
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.roicol : Pydicom Object
            New ROI Collection (RTS/SEG).
        """
        
        roicolMod = params.cfgDict['roicolMod']
        
        if roicolMod == 'RTSTRUCT':
            self.roicol = create_rts(
                srcDataset, trgDataset, newDataset, params
                )
        else:
            self.roicol = create_seg(
                srcDataset, trgDataset, newDataset, params
                )
    
    def error_check_roicol(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        """
        Error check the new ROI Collection (RTS/SEG).
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        trgDataset : Importer Object
            Object containing various data for the target DICOM series.
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.logList : list of strs
            A list of strings that describe any errors that are found.  If no 
            errors are found an empty list ([]) will be returned.
        self.Nerrors : int
            The number of errors found.
        """
        
        roicolMod = params.cfgDict['roicolMod']
        
        if roicolMod == 'RTSTRUCT':
            self.logList, self.Nerrors = error_check_rts(
                self.roicol, trgDataset, params
                )
        else:
            self.logList, self.Nerrors = error_check_seg(
                self.roicol, trgDataset, params
                )
    
    def export_roicol(self, params, fname=''):
        """
        Export ROI Collection (RTS/SEG) to disk.  
        
        Parameters
        ----------
        self.roicol : Pydicom object
            New ROI Collection (RTS/SEG) to be exported.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        fname : str, optional
            File name to assign. The default value is ''. If empty the file
            name will be of the form:
                "{runID}_{roicolMod}_{currentDateTime}.dcm"
            where roicolMod = RTSTRUCT or SEG
        
        Returns
        -------
        self.roicolFpath : str
            Full path of the exported new RTS/SEG file.
        """
        
        cfgDict = params.cfgDict
        roicol = self.roicol
        
        runID = cfgDict['runID']
        #srcRoicolFpath = cfgDict['srcRoicolFpath']
        roicolMod = cfgDict['roicolMod']
        #useDroForTx = cfgDict['useDroForTx']
        
        if roicolMod == 'RTSTRUCT': 
            exportDir = cfgDict['rtsExportDir']
        else:
            exportDir = cfgDict['segExportDir']
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        if fname == '':
            #if useDroForTx:
            #    fname = f'{runID}_DRO_newRoicol_{currentDateTime}.dcm'
            #else:
            #    fname = f'{runID}_newRoicol_{currentDateTime}.dcm' 
            fname = f'{runID}_{roicolMod}_{currentDateTime}.dcm' 
        
        if not os.path.isdir(exportDir):
            #os.mkdir(ExportDir)
            Path(exportDir).mkdir(parents=True)
        
        fpath = os.path.join(exportDir, fname)
        
        roicol.save_as(fpath)
            
        print(f'New {roicolMod} exported to:\n {fpath}\n')
        
        self.roicolFpath = fpath
        
    