# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 19:07:06 2021

@author: ctorti
"""

from importlib import reload
import dicom_tools.create_rts
reload(dicom_tools.create_rts)
import dicom_tools.create_seg
reload(dicom_tools.create_seg)
import dicom_tools.error_check_rts
reload(dicom_tools.error_check_rts)
import dicom_tools.error_check_seg
reload(dicom_tools.error_check_seg)
import xnat_tools.im_sessions_exps
reload(xnat_tools.im_sessions_exps)
import xnat_tools.im_assessors
reload(xnat_tools.im_assessors)
import plotting_tools.general
reload(plotting_tools.general)

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
from xnat_tools.im_sessions_exps import get_exp_id_from_label
from xnat_tools.im_assessors import upload_im_asr
from plotting_tools.general import (
    plot_pixarrs_from_list_of_segs_and_dicom_ims, 
    plot_contours_from_list_of_rtss_and_dicom_ims
    )

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
            
            # datetime when RTS file was generated without fractional seconds:
            dateTime = roicol.StructureSetDate + '_' +\
                str(round(float(roicol.StructureSetTime)))
        else:
            exportDir = cfgDict['segExportDir']
            
            # datetime when the SEG file was generated without fractional s:
            dateTime = roicol.ContentDate + '_' +\
                str(round(float(roicol.ContentTime)))
        
        #currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        if fname == '':
            #if useDroForTx:
            #    fname = f'{runID}_DRO_newRoicol_{currentDateTime}.dcm'
            #else:
            #    fname = f'{runID}_newRoicol_{currentDateTime}.dcm'
            #fname = f'{runID}_{roicolMod}_{currentDateTime}.dcm'
            fname = f'{runID}_{roicolMod}_{dateTime}.dcm'
        
        if not os.path.isdir(exportDir):
            #os.mkdir(ExportDir)
            Path(exportDir).mkdir(parents=True)
        
        fpath = os.path.join(exportDir, fname)
        
        roicol.save_as(fpath)
            
        print(f'New {roicolMod} exported to:\n {fpath}\n')
        
        self.roicolFpath = fpath
      
    def upload_roicol(self, params, collLab=''):
        """
        Upload ROI Collection to XNAT.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        collLab : str, optional
            Label to apply to the ROI Collection.  If not provided the file 
            name from self.roicolFpath will be used. The default value is ''.
        
        Returns
        -------
        None.
        """
        
        roicolFpath = self.roicolFpath
        
        xnatSession = params.xnatSession
        cfgDict = params.cfgDict
        
        url = cfgDict['url']
        projID = cfgDict['projID']
        subjLab = cfgDict['subjLab']
        expLab = cfgDict['trgExpLab']
        
        expID = get_exp_id_from_label(
            url=url, 
            proj_id=projID,
            subj_label=subjLab,
            exp_label=expLab,
            session=params.xnatSession
        )

        xnatSession = upload_im_asr(
            roicol_fpath=roicolFpath, 
            url=url,
            proj_id=projID, 
            session_id=expID, 
            coll_label=collLab,
            session=xnatSession
        )
        
    def plot_roi_over_dicoms(
            self, srcDataset, trgDataset, newDataset, params
            ):
        """ 
        Plot copied/propagated ROI overlaid on DICOM images 
        """
        
        print('* Plotting ROI over DICOMs..\n')
        
        cfgDict = params.cfgDict
        runID = cfgDict['runID']
        roicolMod = cfgDict['roicolMod']
        exportPlot = cfgDict['exportPlots']
        #p2c = cfgDict['p2c']
        p2c = False
        useCaseToApply = cfgDict['useCaseToApply']
        forceReg = cfgDict['forceReg']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        initMethod = cfgDict['initMethod']
        resInterp = cfgDict['resInterp']
        #trgIm = trgDataset.image
        #resIm = newDataset.image # resampled source or source-registered-to-target 
        
        #resExportDir = cfgDict['resPlotsExportDir']
        
        if roicolMod == 'RTSTRUCT':
            roiExportDir = cfgDict['rtsPlotsExportDir']
        else:
            roiExportDir = cfgDict['segPlotsExportDir']
        
        # Prepare plot title for new dataset:
        method = 'Src ROI '
        if useCaseToApply in ['1', '2a']:
            method += 'NRPC to Trg over Trg DICOMs '
        elif useCaseToApply == '2b':
            method += 'RPC to Trg over Trg DICOMs '
        elif useCaseToApply in ['3a', '4a']:
            method += f'NRPP to Trg over Trg DICOMs \n(res, {resInterp})'
        elif useCaseToApply in ['3b', '4b']:
            method += f'RPP to Trg over Trg DICOMs \n(res, {resInterp})'
        elif useCaseToApply == '5a':
            if useDroForTx:
                method += 'NRPP to Trg over Trg DICOMs ' +\
                     f'\n({regTxName} DRO, {resInterp})'
            else:
                method += 'NRPP to Trg over Trg DICOMs ' +\
                    f'\n({regTxName} reg, {initMethod}, {resInterp})'
        elif useCaseToApply == '5b':
            if useDroForTx:
                method += 'RPP to Trg over Trg DICOMs ' +\
                    f'\n({regTxName} DRO, {resInterp})'
            else:
                method += 'RPP to Trg over Trg DICOMs ' +\
                    f'\n({regTxName} reg, {initMethod}, {resInterp})'
        
        #print(method)
        
        listOfRoicol = [srcDataset.roicol, self.roicol]
        """
        Note: newDataset does not have a dicomDir since DICOMs are not
        generated for the resampled/registered source dataset, but for the
        purposes of pulling metadata (sizes, spacings, IPPs, directions, etc)
        relevant to the resampled/registered dataset, getting it from the
        target dataset is equivalent.
        """
        listOfDicomDirs = [srcDataset.dicomDir, trgDataset.dicomDir]
        
        listOfImages = [srcDataset.dcmIm, newDataset.dcmIm]
        listOfPlotTitles = ['Src ROI over Src DICOMs', method]
        
        """
        plot_pixarrs_from_list_of_segs_v3(
        listOfSegs, listOfDicomDirs, listOfPlotTitles, exportPlot=exportPlot,
            exportDir=exportDir, runID=runID, useCaseToApply=useCaseToApply,
            forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
            initMethod=initMethod, resInterp=resInterp, p2c=p2c
        )
        """
        
        # Prepare filename for exported plot:
        
        currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        if useDroForTx:
            fname = f'{runID}_DRO'
        else:
            fname = f'{runID}_'
        fname += method.replace(' ', '_').replace(',', '')
        fname = fname.replace('(', '').replace(')', '').replace('\n', '')
        fname += f'_{currentDateTime}'
        
        if roicolMod == 'RTSTRUCT':
            plot_contours_from_list_of_rtss_and_dicom_ims(
                listOfRoicol, listOfImages, listOfDicomDirs, listOfPlotTitles, 
                exportPlot=exportPlot, exportDir=roiExportDir,
                runID=runID, useCaseToApply=useCaseToApply,
                forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
                initMethod=initMethod, resInterp=resInterp, 
                fname=fname, fontSize=12, p2c=p2c
                #fname='', p2c=p2c
            )
        else:
            plot_pixarrs_from_list_of_segs_and_dicom_ims(
                listOfRoicol, listOfImages, listOfDicomDirs, listOfPlotTitles, 
                exportPlot=exportPlot, exportDir=roiExportDir,
                runID=runID, useCaseToApply=useCaseToApply,
                forceReg=forceReg, useDroForTx=useDroForTx, regTxName=regTxName,
                initMethod=initMethod, resInterp=resInterp, 
                fname=fname, fontSize=12, p2c=p2c
                #fname='', p2c=p2c
            )