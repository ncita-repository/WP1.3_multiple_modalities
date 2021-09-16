# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:58:18 2021

@author: ctorti
"""


import os
from pathlib import Path
import time
from io_tools.exports import export_dict_to_json
from io_tools.imports import import_dict_from_json


def create_config_files(cfgDir, password=""):
    """
    Generate a JSON files for all run configurations.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the JSON files.
    password : str, optional
        Password for XNAT connection. Default value is "".
    
    Returns
    -------
    None.
    
    Notes
    -----
    XNAT config settings:
        If the values for the keys "url", "user" or "password" are left empty,    
        the user will be prompted to input it/them.
    
    Resampling / Registration / Transformation settings:
        If the values for the keys "srcFidsFpath" and "trgFidsFpath" are not
        empty, the files containing the fiducials will be used to initialise
        image registration.  If no fiducials are available leave empty.
        
        By default image registration will only be performed if the source and 
        target images have different FrameOfReferenceUIDs.  However it is possible
        for differing ImagePositionPatient and/or ImageOrientationPatient to exist be
        between two images with the same FoR (perhaps down to human error). 
        
        If results of an ROI copy between two datasets with the same FoR arise
        try setting the parameter 'forceReg' to True (it is False by default) to
        force registration to be performed despite equal FoRs.
        
        Image registration can be initialised using one of three methods:
            1. 'landmarks'
            2. 'moments' (or 'centerofgravity')
            3. 'geometry' (or 'geometricalcenter')
            
        If 'landmarks' is chosen but srcFidsFpath and trgFidsFpath are empty or
        not formatted properly, a 'moments'-based initialisation will be performed.
        
        If useDroForTx = False, the transformation parameters will be determined
        based on image registration rather than use of the parameters present in
        the DRO.
        
        Suitable limits for the maximum iterations ("maxIters") are 512 for 
        rigid/affine, and 100 for bspline registrations.
        
        The following three possible values are allowed for the key "resInterp"
        (resampling interpolator):
            1. 'NearestNeighbor'
            2. 'LabelGaussian'
            3. 'BlurThenLinear'
    """
    
    # Current working directory:
    #cwd = os.getcwd()
    
    # Default download directory for XNAT data:
    ddd = os.path.join(Path.home(), "Downloads", "xnat_downloads")
    
    # Current date:
    currentDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    """ 
    Select one of the following suggested download directories for XNAT data:
    """
    #rootExportDir = ddd
    rootExportDir = os.path.join(ddd, currentDate)
    #rootExportDir = os.path.join(ddd, '2021-07-12') # change accordingly
    
    """
    Define the export directories for the (new) target ROI Collections, 
    DRO, plots, transforms, binary label maps, logs and configuration files:
    """
    rtsExportDir = os.path.join(rootExportDir, 'new_RTS')
    segExportDir = os.path.join(rootExportDir, 'new_SEG')
    
    droExportDir = os.path.join(rootExportDir, 'new_DRO')
    
    rtsPlotsExportDir = os.path.join(rootExportDir, 'plots_RTS')
    segPlotsExportDir = os.path.join(rootExportDir, 'plots_SEG')
    
    # Location where resampled/transformed/registered plots will be saved:
    resPlotsExportDir = os.path.join(rootExportDir, 'plots_res')
    
    txExportDir = os.path.join(rootExportDir, 'transforms')
    imExportDir = os.path.join(rootExportDir, 'images')
    labimExportDir = os.path.join(rootExportDir, 'label_images')
    
    logsExportDir = os.path.join(rootExportDir, 'logs')
    
    #cfgExportDir = os.path.join(rootExportDir, 'configs')
    
    """
    Define the directory path of the sample DROs and fiducials:
    """
    sampleDroDir = r'C:\Code\WP1.3_multiple_modalities\sample_DROs'
    fidsDir = r'C:\Code\WP1.3_multiple_modalities\OOP\fiducials'
    
    """
    Chose whether or not to export the new ROI Collection (i.e. RTS or SEG),
    DRO, transforms, label images, plots and logs:
    """
    #exportRoicol = False
    exportRoicol = True
    exportIm = False # 3D DICOM images
    #exportDro = False
    exportDro = True
    exportTx = True
    exportLabim = True
    #exportPlots = False
    exportPlots = True
    #exportLogs = False
    exportLogs = True
    
    
    """
    The parameter addToRoicolLabel allows for adding text to the 
    StructureSetLab (SSL) of the new DICOM-RTSTRUCT or SeriesDescription (SD)
    of the new DICOM-SEG file, 
    e.g. New SSL/SD = SSD/SD copied from the Source RTS/SEG + addToRoicolLabel
    The default value of addToRoicolLab will the empty string, which will
    result in the addition of the timestamp of the file creation,
    e.g. "{copied SSL/SD} 20210323 101748".
    """
    addToRoicolLab = '' # the current timestamp will be added
    
    """ 
    When searching XNAT for the source ROI Collection that matches the required
    metadata multiple hits may arise. Chose how to proceed from the following
    3 options:
        1. Select the oldest ROI Collection (whichSrcRoicol = 'oldest')
        2. Select the newest ROI Collection (whichSrcRoicol = 'newest')
        3. Allow the user to decide (whichSrcRoicol = 'user')
    """
    whichSrcRoicol = 'oldest'
    #whichSrcRoicol = 'newest'
    #whichSrcRoicol = 'user'
    
    
    """
    Choose whether or not to use a valid DRO to transform (propagate) the 
    source ROI Collection to the target domain or to perform image registration
    irrespective of a valid DRO:
    
    Update: Rather than defining the value here it is defined in the definition
    of the dictionaries (see below).
    """
    #useDroForTx = True # use transform in DRO, by-passing image registration
    #useDroForTx = False # apply image registration even if DRO is available
    
    """
    Chose whether or not to upload the (new) target DRO to XNAT:
    """
    #uploadDro = False
    uploadDro = True # see Notes
    
    """
    Chose whether or not to overwrite an existing DRO to XNAT:
    """
    overwriteDro = False
    #overwriteDro = True # see Notes
    
    """
    Chose whether or not to print results to console (e.g. for debugging):
    """
    #p2c = False
    p2c = True
    
    # Initialise dictionary (of dictionaries) to store all configurations:
    cfg = {}
    
    # Datasets to propagate a SEG using affine registration of soft tissue:
    runID = 'Soft_tissue_seg_aff'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': password,
        'projID' : 'Soft-tissue-Sarcoma', 
        'subjLab' : 'STS_004',
        'srcExpLab' : 'STS_004_PETCT',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'CT_Ser2_tumour',
        #'srcRoicolName' : 'Left_thigh_tumour',
        'srcRoiName' : None,
        'roicolMod' : 'SEG',
        'trgExpLab' : 'STS_004_MR',
        'trgScanID' : '501',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : os.path.join(fidsDir, r'STS_004_CT_fiducials.txt'),
        'trgFidsFpath' : os.path.join(fidsDir, r'STS_004_MR_fiducials.txt'), 
        'forceReg' : False,
        'regTxName' : 'affine',
        'initMethod' : 'landmarks',
        'maxIters' : 512,
        #'useDroForTx' : True,
        'useDroForTx' : False,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1),
        'sampleDroDir' : sampleDroDir,
        'rootExportDir' : rootExportDir,
        #'cfgExportDir' : cfgExportDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir,
        'exportRoicol' : exportRoicol,
        'exportDro' : exportDro,
        'exportTx' : exportTx,
        'exportIm' : exportIm,
        'exportLabim' : exportLabim,
        'exportPlots' : exportPlots,
        'exportLogs' : exportLogs,
        'whichSrcRoicol' : whichSrcRoicol,
        'addToRoicolLab' : addToRoicolLab,
        'uploadDro' : uploadDro,
        'overwriteDro' : overwriteDro,
        'p2c' : p2c
        }
    
    # Datasets to propagate a RTS using affine registration of soft tissue 
    # (start with copy of 'Soft_tissue_seg_aff' since nearly identical):
    runID = 'Soft_tissue_rts_aff'
    cfg[runID] = dict(cfg['Soft_tissue_seg_aff'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    #cfg[runID]['srcRoicolName'] = 'CT_Ser2_tumour'
    
    # Datasets to propagate SEG using bspine registration of soft tissue:
    runID = 'Soft_tissue_seg_def'
    cfg[runID] = dict(cfg['Soft_tissue_seg_aff'])
    cfg[runID]['runID'] = runID
    cfg[runID]['regTxName'] = 'bspline'
    cfg[runID]['maxIters'] = 100
    
    # Datasets to propagate RTS using bspine registration of soft tissue:
    runID = 'Soft_tissue_rts_def'
    cfg[runID] = dict(cfg['Soft_tissue_seg_aff'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    
    
    
    """
    Following are datasets for test runs using the NCITA_TEST collection (brain
    images) to either propagate ("relationship-preserving copy") or copy 
    ("direct copy") SEG and RTSTRUCT ROI Collections using affine registrations:
    
    https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g/edit#
    """
    
    # Datasets to propagate RTS within the same DICOM series:
    runID = 'NCITA_test_RR1'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': password,
        'projID' : 'NCITA_TEST', 
        'subjLab' : 'TCGA-BB-A5HY',
        'srcExpLab' : 'Session2',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Left eye CT',
        'srcRoiName' : None,
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'Session2',
        'trgScanID' : '4',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : None,
        'trgFidsFpath' : None, 
        'forceReg' : False,
        'regTxName' : 'affine',
        #'initMethod' : 'moments', # 07/09/21
        'initMethod' : 'geometry', # 07/09/21
        'maxIters' : 512,
        'useDroForTx' : False,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1),
        'sampleDroDir' : sampleDroDir,
        'rootExportDir' : rootExportDir,
        #'cfgExportDir' : cfgExportDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir,
        'exportRoicol' : exportRoicol,
        'exportDro' : exportDro,
        'exportTx' : exportTx,
        'exportIm' : exportIm,
        'exportLabim' : exportLabim,
        'exportPlots' : exportPlots,
        'exportLogs' : exportLogs,
        'whichSrcRoicol' : whichSrcRoicol,
        'addToRoicolLab' : addToRoicolLab,
        'uploadDro' : uploadDro,
        'overwriteDro' : overwriteDro,
        'p2c' : p2c
        }
    
    # Datasets to propagate RTS across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_RR2'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session4'
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    cfg[runID]['trgExpLab'] = 'Session4'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR3'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_RR4'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '11'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR5'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '11'
    cfg[runID]['srcRoicolName'] = 'Ventricles'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR7'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '12'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR8'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '13'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR9'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR10'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session3'
    cfg[runID]['trgScanID'] = '6'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR11'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session3'
    cfg[runID]['trgScanID'] = '6'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR12'
    cfg[runID] = dict(cfg['NCITA_test_RR11'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session5'
    cfg[runID]['srcScanID'] = '12'
    cfg[runID]['srcRoicolName'] = 'Left eye'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR13'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session1'
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR14'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session1'
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR15'
    cfg[runID] = dict(cfg['NCITA_test_RR12'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session1'
    cfg[runID]['trgScanID'] = '5'
    
    # Dataset to propagate RTS to an existing RTS:
    """ 
    Note: This was added to the original list of test configurations.
    """
    runID = 'NCITA_test_RR16'
    cfg[runID] = dict(cfg['NCITA_test_RR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    # Datasets to copy RTS within the same DICOM series:
    """ Use case 1 """
    runID = 'NCITA_test_RD1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 5
    cfg[runID]['srcRoicolName'] = 'Right eye CT (single slice)'
    cfg[runID]['trgScanID'] = '2'
    cfg[runID]['trgSlcNum'] = 6
    
    # Datasets to copy RTS within the same DICOM series:
    """ Use case 2a """
    runID = 'NCITA_test_RD2'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '4'
    
    # Datasets to copy RTS across two DICOM series with different IOPs:
    """ Use case 3a
    Note: Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_RD3'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 65
    
    # Datasets to copy RTS from one DICOM series to another with equal slice
    # thickness:
    """ Use case 5a
    Note: Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_RD4'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session6'
    cfg[runID]['trgScanID'] = '3'
    cfg[runID]['trgSlcNum'] = 8
    
    # Datasets to copy RTS from one DICOM series to another with smaller
    # slice thickness:
    """ Use case 5a
    Note: Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_RD5'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 78
    
    # Datasets to copy RTS across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    """ Use case 5a """
    runID = 'NCITA_test_RD6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 18
    cfg[runID]['trgSlcNum'] = 21
    
    # Datasets to propagate SEG within the same DICOM series:
    runID = 'NCITA_test_SR1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_SR2'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR3'
    cfg[runID] = dict(cfg['NCITA_test_RR3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_SR4'
    cfg[runID] = dict(cfg['NCITA_test_RR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR5'
    cfg[runID] = dict(cfg['NCITA_test_RR5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR6'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR7'
    cfg[runID] = dict(cfg['NCITA_test_RR7'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR8'
    cfg[runID] = dict(cfg['NCITA_test_RR8'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR9'
    cfg[runID] = dict(cfg['NCITA_test_RR9'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR10'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR11'
    cfg[runID] = dict(cfg['NCITA_test_RR11'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR12'
    cfg[runID] = dict(cfg['NCITA_test_RR12'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR13'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR14'
    cfg[runID] = dict(cfg['NCITA_test_RR14'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR15'
    cfg[runID] = dict(cfg['NCITA_test_RR15'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    """
    Note: This was added to the original list of test configurations.
    """
    runID = 'NCITA_test_SR16'
    cfg[runID] = dict(cfg['NCITA_test_RR16'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD1'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD2'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with different IOPs:
    """
    Note: Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_SD3'
    cfg[runID] = dict(cfg['NCITA_test_RD3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with equal slice
    # thickness:
    """
    Note: Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_SD4'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with smaller
    # slice thickness:
    """
    Note: Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_SD5'
    cfg[runID] = dict(cfg['NCITA_test_RD5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    runID = 'NCITA_test_SD6'
    cfg[runID] = dict(cfg['NCITA_test_RD6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # TODO Verify the above NCITA_TEST configs and comments
    
    #export_dict_to_json(dictionary=cfg, filename='main_params.json',
    #                    exportDir=cfgDir)
    
    # Export each sub-dictionary to a unique JSON file:
    for runID in cfg.keys():
        
        export_dict_to_json(
            dictionary=cfg[runID], filename=f'{runID}.json', exportDir=cfgDir
            )

def create_main_params_file_OLD(cfgDir):
    """
    Generate a JSON file containing parameters for defined configurations of
    different runs.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the config file.
    
    Returns
    -------
    None.
    
    Notes
    -----
    XNAT config settings:
        If the values for the keys "url", "user" or "password" are left empty,    
        the user will be prompted to input it/them.
    
    Resampling / Registration / Transformation settings:
        If the values for the keys "srcFidsFpath" and "trgFidsFpath" are not
        empty, the files containing the fiducials will be used to initialise
        image registration.  If no fiducials are available leave empty.
        
        By default image registration will only be performed if the source and 
        target images have different FrameOfReferenceUIDs.  However it is possible
        for differing ImagePositionPatient and/or ImageOrientationPatient to exist be
        between two images with the same FoR (perhaps down to human error). 
        
        If results of an ROI copy between two datasets with the same FoR arise
        try setting the parameter 'forceReg' to True (it is False by default) to
        force registration to be performed despite equal FoRs.
        
        Image registration can be initialised using one of three methods:
            1. 'landmarks'
            2. 'moments' (or 'centerofgravity')
            3. 'geometry' (or 'geometricalcenter')
            
        If 'landmarks' is chosen but srcFidsFpath and trgFidsFpath are empty or
        not formatted properly, a 'moments'-based initialisation will be performed.
        
        If useDroForTx = False, the transformation parameters will be determined
        based on image registration rather than use of the parameters present in
        the DRO.
        
        Suitable limits for the maximum iterations ("maxIters") are 512 for 
        rigid/affine, and 100 for bspline registrations.
        
        The following three possible values are allowed for the key "resInterp"
        (resampling interpolator):
            1. 'NearestNeighbor'
            2. 'LabelGaussian'
            3. 'BlurThenLinear'
    """
    
    cfg = {}
    
    # Datasets to propagate a SEG using affine registration of soft tissue:
    cfg['Soft_tissue_seg_aff'] = {
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': "",
        'projID' : 'Soft-tissue-Sarcoma', 
        'subjLab' : 'STS_004',
        'srcExpLab' : 'STS_004_PETCT',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'CT_Ser2_tumour',
        #'srcRoicolName' : 'Left_thigh_tumour',
        'srcRoiName' : None,
        'roicolMod' : 'SEG',
        'trgExpLab' : 'STS_004_MR',
        'trgScanID' : '501',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : os.path.join(cwd, r'STS_004_CT_fiducials.txt'),
        'trgFidsFpath' : os.path.join(cwd, r'STS_004_MR_fiducials.txt'), 
        'forceReg' : False,
        'regTxName' : 'affine',
        'initMethod' : 'landmarks',
        'maxIters' : 512,
        'useDroForTx' : True,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1)
        }
    
    # Datasets to propagate a RTS using affine registration of soft tissue 
    # (start with copy of 'Soft_tissue_seg_aff' since nearly identical):
    cfg['Soft_tissue_rts_aff'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_rts_aff']['roicolMod'] = 'RTSTRUCT'
    #cfg['Soft_tissue_rts_aff']['srcRoicolName'] = 'CT_Ser2_tumour'
    
    # Datasets to propagate SEG using bspine registration of soft tissue:
    cfg['Soft_tissue_seg_def'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_seg_def']['regTxName'] = 'bspline'
    cfg['Soft_tissue_seg_def']['maxIters'] = 100
    
    # Datasets to propagate RTS using bspine registration of soft tissue:
    cfg['Soft_tissue_rts_def'] = dict(cfg['Soft_tissue_seg_aff'])
    cfg['Soft_tissue_rts_def']['roicolMod'] = 'RTSTRUCT'
    
    
    
    """
    Following are datasets for test runs using the NCITA_TEST collection (brain
    images) to either propagate ("relationship-preserving copy") or copy 
    ("direct copy") SEG and RTSTRUCT ROI Collections using affine registrations:
    
    https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g/edit#
    """
    
    # Datasets to propagate RTS within the same DICOM series:
    runID = 'NCITA_test_RR1'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin", 
        'password': "",
        'projID' : 'NCITA_TEST', 
        'subjLab' : 'TCGA-BB-A5HY',
        'srcExpLab' : 'Session2',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Left eye CT',
        'srcRoiName' : None,
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'Session2',
        'trgScanID' : '4',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : None,
        'trgFidsFpath' : None, 
        'forceReg' : False,
        'regTxName' : 'affine',
        'initMethod' : 'moments',
        'maxIters' : 512,
        'useDroForTx' : False,
        'applyPreResBlur' : False,
        'preResVar' : (1,1,1),
        'resInterp' : 'BlurThenLinear',
        'applyPostResBlur' : True,
        'postResVar' : (1,1,1)
        }
    
    # Datasets to propagate RTS across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_RR2'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session4'
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    cfg[runID]['trgExpLab'] = 'Session4'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR3'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_RR4'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '11'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR5'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '11'
    cfg[runID]['srcRoicolName'] = 'Ventricles'
    cfg[runID]['trgScanID'] = '10'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR7'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '12'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR8'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '13'
    
    # Datasets to propagate RTS from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_RR9'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '14'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR10'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session3'
    cfg[runID]['trgScanID'] = '6'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR11'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '13'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR12'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session5'
    cfg[runID]['srcScanID'] = '12'
    cfg[runID]['srcRoicolName'] = 'Left eye'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR13'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session1'
    cfg[runID]['trgScanID'] = '5'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR14'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    
    # Datasets to propagate RTS across two DICOM series with different IOPs:
    runID = 'NCITA_test_RR15'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session5'
    cfg[runID]['srcScanID'] = '12'
    cfg[runID]['srcRoicolName'] = 'Left eye'
    
    # Dataset to propagate RTS to an existing RTS:
    """ 
    Note:
        This was not one of the original test configurations so unknown whether
    this will work.
    """
    runID = 'NCITA_test_RR16'
    cfg[runID] = dict(cfg['NCITA_test_RR14'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session4'
    cfg[runID]['trgScanID'] = '11'
    cfg[runID]['trgRoicolName'] = 'Ventricles'
    
    # Datasets to copy RTS within the same DICOM series:
    runID = 'NCITA_test_RD1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 5
    cfg[runID]['srcRoicolName'] = 'Right eye CT (single slice)'
    cfg[runID]['trgScanID'] = '2'
    cfg[runID]['trgSlcNum'] = 6
    
    # Datasets to copy RTS within the same DICOM series:
    runID = 'NCITA_test_RD2'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '4'
    
    # Datasets to copy RTS across two DICOM series with different IOPs:
    """
    Note:
        Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_RD3'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 65
    
    # Datasets to copy RTS from one DICOM series to another with equal slice
    # thickness:
    """
    Note:
        Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_RD4'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgExpLab'] = 'Session6'
    cfg[runID]['trgScanID'] = '3'
    cfg[runID]['trgSlcNum'] = 8
    
    # Datasets to copy RTS from one DICOM series to another with smaller
    # slice thickness:
    """
    Note:
        Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_RD5'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '5'
    cfg[runID]['trgSlcNum'] = 78
    
    # Datasets to copy RTS across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    runID = 'NCITA_test_RD6'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = 18
    cfg[runID]['trgSlcNum'] = 21
    
    # Datasets to propagate SEG within the same DICOM series:
    runID = 'NCITA_test_SR1'
    cfg[runID] = dict(cfg['NCITA_test_RR1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with the same slice
    # thickness:
    runID = 'NCITA_test_SR2'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR3'
    cfg[runID] = dict(cfg['NCITA_test_RR3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with smaller
    # slice thickness:
    runID = 'NCITA_test_SR4'
    cfg[runID] = dict(cfg['NCITA_test_RR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR5'
    cfg[runID] = dict(cfg['NCITA_test_RR5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR6'
    cfg[runID] = dict(cfg['NCITA_test_RR6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR7'
    cfg[runID] = dict(cfg['NCITA_test_RR7'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR8'
    cfg[runID] = dict(cfg['NCITA_test_RR8'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG from one DICOM series to another with greater
    # slice thickness:
    runID = 'NCITA_test_SR9'
    cfg[runID] = dict(cfg['NCITA_test_RR9'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR10'
    cfg[runID] = dict(cfg['NCITA_test_RR10'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR11'
    cfg[runID] = dict(cfg['NCITA_test_RR11'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR12'
    cfg[runID] = dict(cfg['NCITA_test_RR12'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR13'
    cfg[runID] = dict(cfg['NCITA_test_RR13'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR14'
    cfg[runID] = dict(cfg['NCITA_test_RR14'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to propagate SEG across two DICOM series with different IOPs:
    runID = 'NCITA_test_SR15'
    cfg[runID] = dict(cfg['NCITA_test_RR15'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD1'
    cfg[runID] = dict(cfg['NCITA_test_RD1'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG within the same DICOM series:
    runID = 'NCITA_test_SD2'
    cfg[runID] = dict(cfg['NCITA_test_RD2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with different IOPs:
    """
    Note:
        Source slice 5 maps to Target slices 55-64.
    """
    runID = 'NCITA_test_SD3'
    cfg[runID] = dict(cfg['NCITA_test_RD3'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with equal slice
    # thickness:
    """
    Note:
        Source slice 5 maps to Target slice 7.
    """
    runID = 'NCITA_test_SD4'
    cfg[runID] = dict(cfg['NCITA_test_RD4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG from one DICOM series to another with smaller
    # slice thickness:
    """
    Note:
        Source slice 5 maps to Target slices 67-77.
    """
    runID = 'NCITA_test_SD5'
    cfg[runID] = dict(cfg['NCITA_test_RD5'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Datasets to copy SEG across two DICOM series with equal slice
    # thickness, demonstrating the displacement within the same series:
    runID = 'NCITA_test_SD6'
    cfg[runID] = dict(cfg['NCITA_test_RD6'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # TODO Verify the above NCITA_TEST configs and comments
    
    export_dict_to_json(dictionary=cfg, filename='main_params.json',
                        exportDir=cfgDir)
    
def add_password_to_main_params_file_OLD(cfgDir, password):
    
    fpath = os.path.join(cfgDir, r'main_params.json')
    
    cfg = import_dict_from_json(fpath)
    
    for runID in cfg.keys():
        cfg[runID]['password'] = password
    
    export_dict_to_json(dictionary=cfg, filename='main_params.json',
                        exportDir=cfgDir)
    
def create_general_params_file_OLD(cfgDir):
    """
    Create general parameters.
    
    It is assumed that the general parameters will not likely be changed 
    between runs with different datasets and/or resampling / registering / 
    transforming settings. The file general_params.json will be exported to a 
    directory named "config" within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the config file.
    
    Returns
    -------
    None.
    
    Notes
    -----
    If uploadDro = True, the DRO will be uploaded to XNAT as a subject 
    assessor.
    
    If overwriteDro = True (and uploadDro = True) and there exists a DRO
    (as a subject assessor) in XNAT, the existing DRO will be over-written.
    """
    
    
    """ 
    Select one of the following suggested download directories for XNAT data:
    """
    #rootExportDir = ddd
    rootExportDir = os.path.join(ddd, currentDate)
    #rootExportDir = os.path.join(ddd, '2021-07-12') # change accordingly
    
    """
    Define the export directories for the (new) target ROI Collections, 
    DRO, plots, logs and configuration files:
    """
    rtsExportDir = os.path.join(rootExportDir, 'new_RTS')
    segExportDir = os.path.join(rootExportDir, 'new_SEG')
    
    droExportDir = os.path.join(rootExportDir, 'new_DRO')
    
    rtsPlotsExportDir = os.path.join(rootExportDir, 'plots_RTS')
    segPlotsExportDir = os.path.join(rootExportDir, 'plots_SEG')
    
    logsExportDir = os.path.join(rootExportDir, 'logs')
    
    #cfgExportDir = os.path.join(rootExportDir, 'configs')
    
    """
    Chose whether or not to export the (new) target ROI Collection (i.e.
    RTS or SEG), DRO, plots and logs:
    """
    #exportRoicol = False
    exportRoicol = True
    
    #exportDro = False
    exportDro = True
    
    #exportPlots = False
    exportPlots = True
    #exportLogs = False
    exportLogs = True
    
    
    """
    Chose whether or not to upload the (new) target DRO to XNAT:
    """
    #uploadDro = False
    uploadDro = True # see Notes
    
    """
    Chose whether or not to overwrite an existing DRO to XNAT:
    """
    overwriteDro = False
    #overwriteDro = True # see Notes
    
    """
    Chose whether or not to print results to console (e.g. for debugging):
    """
    #p2c = False
    p2c = True
    
    cfg = {'rootExportDir' : rootExportDir,
           #'cfgExportDir' : cfgExportDir,
           'rtsExportDir' : rtsExportDir,
           'segExportDir' : segExportDir,
           'droExportDir' : droExportDir,
           'rtsPlotsExportDir' : rtsPlotsExportDir,
           'segPlotsExportDir' : segPlotsExportDir,
           'logsExportDir' : logsExportDir,
           'exportTrgRoicol' : exportTrgRoicol,
           'exportDro' : exportDro,
           'exportPlots' : exportPlots,
           'exportLogs' : exportLogs,
           'uploadDro' : uploadDro,
           'overwriteDro' : overwriteDro,
           'p2c' : p2c}
    
    export_dict_to_json(dictionary=cfg, filename='general_params.json',
                        exportDir=cfgDir)

def create_xnat_cfg_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a configuration file for XNAT connections.
    
    The file xnat_config.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Note
    ----
    If the values for keys "user" and "password" are left empty, i.e.
    
    cfg = {"url" : "http://10.1.1.20", 
           "user": "", 
           "password": ""}
    
    the user will be prompted to enter the credentials.
    """
    
    cfg = {"url" : "http://10.1.1.20", 
           "username": "admin", 
           "password": "admin"}
    
    export_dict_to_json(dictionary=cfg, filename='xnat_config.json',
                        exportDir=cfgDir)

def create_src_xnat_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a parameters file for XNAT parameters for Source.
    
    The file xnat_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    roicol => ROI Collection (i.e. RTS or SEG)
    
    projID = XNAT project ID
    subjLab = XNAT subject label
    expLab = XNAT experiment label
    scanID = XNAT scan ID
    slcNum = Slice number within a colection of DICOM image stack
    roicolMod = Modality of the ROI Collection
    roicolName = XNAT name of the ROI Collection
    roiName = DICOM ROIName or SegmentLabel of an ROI Collection
    """
    
    cfg = {'projID' : 'Soft-tissue-Sarcoma', 
           'subjLab' : 'STS_004',
           'expLab' : 'STS_004_PETCT',
           'scanID' : '2',
           'slcNum' : None,
           'roicolMod' : 'RTSTRUCT',
           'roicolName' : 'CT_Ser2_tumour',
           'roiName' : None}
    
    export_dict_to_json(dictionary=cfg, filename='src_xnat_params.json',
                        exportDir=cfgDir)

def create_trg_xnat_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create a parameters file for XNAT parameters for Target.
    
    The file xnat_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    roicol => ROI Collection (i.e. RTS or SEG)
    
    projID = XNAT project ID
    subjLab = XNAT subject label
    expLab = XNAT experiment label
    scanID = XNAT scan ID
    slcNum = Slice number within a colection of DICOM image stack
    roicolMod = Modality of the ROI Collection
    roicolName = XNAT name of the ROI Collection
    roiName = DICOM ROIName or SegmentLabel of an ROI Collection
    """
    
    cfg = {'projID' : 'Soft-tissue-Sarcoma', 
           'subjLab' : 'STS_004',
           'expLab' : 'STS_004_MR',
           'scanID' : '501',
           'slcNum' : None,
           'roicolMod' : 'RTSTRUCT',
           'roicolName' : None,
           'roiName' : None}
           #'addToTrgRoicolName' : f' from {srcRoicolName}'}
    
    export_dict_to_json(dictionary=cfg, filename='trg_xnat_params.json',
                        exportDir=cfgDir)

def create_res_params_OLD(cfgDir):
    """
    OBSOLETE (09/08/21)
    
    Create resampling / registration / transformation-related parameters.
    
    The file res_params.json will be exported to a directory named "config"
    within the current working directory.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory containing config files.
    
    Returns
    -------
    None.
    
    Notes
    -----
    By default image registration will only be performed if the source and 
    target images have different FrameOfReferenceUIDs.  However it is possible
    for differing ImagePositionPatient and/or ImageOrientationPatient to exist be
    between two images with the same FoR (perhaps down to human error). 
    
    If results of an ROI copy between two datasets with the same FoR arise
    try setting the parameter 'forceReg' to True (it is False by default) to
    force registration to be performed despite equal FoRs.
    
    Image registration can be initialised using one of three methods:
        1. 'landmarks'
        2. 'moments' (or 'centerofgravity')
        3. 'geometry' (or 'geometricalcenter')
        
    If 'landmarks' is chosen but srcFidsFpath and trgFidsFpath are empty or
    not formatted properly, a 'moments'-based initialisation will be performed.
    
    If useDroForTx = False, the transformation parameters will be determined
    based on image registration rather than use of the parameters present in
    the DRO.
    """
    
    
    """ 
    Define the filepaths to the source and target fiducials txt files (if
    applicable):
    """
    #srcFidsFpath = '' # no fiducials
    #trgFidsFpath = '' # no fiducials
    srcFidsFpath = os.path.join(cwd, r'STS_004_CT_fiducials.txt')
    trgFidsFpath = os.path.join(cwd, r'STS_004_MR_fiducials.txt')
    
    
    """ 
    Define registration parameters:
    """
    forceReg = False
    #forceReg = True # see Notes
    
    #regTxName = 'rigid'
    regTxName = 'affine'
    #regTxName = 'bspline'
    
    initMethod = 'landmarks' # see Notes
    #initMethod = 'moments'
    #initMethod = 'geometry'
    
    maxIters = 512 # suitable for rigid/affine
    #maxIters = 100 # suitable for bspline
    
    
    """ 
    Define resampling / transforming parameters:
    """
    useDroForTx = True
    #useDroForTx = False # see Notes
    
    applyPreResBlur = False
    #applyPreResBlur = True
    
    preResVar = (1,1,1) # variance
    
    #resInterp = 'NearestNeighbor'
    #resInterp = 'LabelGaussian'
    resInterp = 'BlurThenLinear'
    
    #applyPostResBlur = False
    applyPostResBlur = True
    
    postResVar = (1,1,1) # variance
    
    cfg = {'srcFidsFpath' : srcFidsFpath,
           'trgFidsFpath' : trgFidsFpath, 
           'forceReg' : forceReg,
           'regTxName' : regTxName,
           'initMethod' : initMethod,
           'maxIters' : maxIters,
           'useDroForTx' : useDroForTx,
           'applyPreResBlur' : applyPreResBlur,
           'preResVar' : preResVar,
           'resInterp' : resInterp,
           'applyPostResBlur' : applyPostResBlur,
           'postResVar' : postResVar}
    
    export_dict_to_json(dictionary=cfg, filename='reg_res_tx_params.json',
                        exportDir=cfgDir)