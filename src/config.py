# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:58:18 2021

@author: ctorti
"""


import os
#from pathlib import Path
import time
import argparse
from io_tools.exports import export_dict_to_json
#from io_tools.imports import import_dict_from_json


def create_config_files(cfgDir):
    """
    Generate a JSON files for all run configurations.
    
    Parameters
    ----------
    cfgDir : str
        The path to the directory to export the JSON files.
    
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
    
    # Current date:
    currentDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    # Try to get current working directory from an environment variable - this 
    # will be updated by running update_cfgDict() in io_tools.fetch_config.py:
    """
    cwd = os.getenv('workdir')
    print(f'os.getenv (in config.py) = {cwd}')
    if cwd == None:
        cwd = os.getcwd()
    print(f'cwd (in config.py) = {cwd}')
    """
    cwd = "" # this will be updated later
    #cwd = "src"  # 03/11/21
    
    # Default download directory for XNAT data:
    #downloadDir = os.path.join(Path.home(), "Downloads", "xnat_downloads")
    #downloadDir = os.path.join(cwd, "xnat_downloads", currentDate)
    # Default download directory for XNAT data within cwd:
    downloadDir = os.path.join(r"xnat_downloads", currentDate)
    
    # Root input and output directories:
    #rootInputDir = os.path.join(cwd, r'inputs')
    #rootOutputDir = os.path.join(cwd, r'outputs')
    # Inputs and outputs directories within cwd:
    inputsDir = r'inputs' # 03/11/21
    outputsDir = r'outputs' # 03/11/21
    #inputsDir = os.path.join(cwd, r'inputs') # 03/11/21
    #outputsDir = os.path.join(cwd, r'outputs') # 03/11/21
    
    # Directories containing sample DROs and fiducials:
    #sampleDroDir = os.path.join(rootInputDir, r'sample_DROs')
    #fidsDir = os.path.join(rootInputDir, r'fiducials')
    # Directories containing sample DROs and fiducials within cwd:
    sampleDroDir = os.path.join(inputsDir, r'sample_DROs')
    #fidsDir = os.path.join(outputsDir, r'fiducials') # 03/11/21
    fidsDir = os.path.join(inputsDir, r'fiducials') # 03/11/21
    
    # Directories for the (new) target ROI Collections, DRO, plots, transforms,
    # binary label maps, logs and configuration files:
    rtsExportDir = os.path.join(outputsDir, r'new_RTS')
    segExportDir = os.path.join(outputsDir, r'new_SEG')
    droExportDir = os.path.join(outputsDir, r'new_DRO')
    rtsPlotsExportDir = os.path.join(outputsDir, 'plots_RTS')
    segPlotsExportDir = os.path.join(outputsDir, 'plots_SEG')
    # resampled/transformed/registered plots
    resPlotsExportDir = os.path.join(outputsDir, r'plots_res')
    txExportDir = os.path.join(outputsDir, r'transforms')
    imExportDir = os.path.join(outputsDir, r'images')
    labimExportDir = os.path.join(outputsDir, r'label_images')
    logsExportDir = os.path.join(outputsDir, r'logs')
    
    
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
    p2c = False
    #p2c = True
    
    # Initialise dictionary (of dictionaries) to store all configurations:
    cfg = {}
    
    # Datasets to propagate a SEG using affine registration of soft tissue:
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin",
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
        'p2c' : p2c,
        'cwd' : cwd,
        'downloadDir' : downloadDir,
        'inputsDir' : inputsDir,
        'outputsDir' : outputsDir,
        'sampleDroDir' : sampleDroDir,
        'fidsDir' : fidsDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir
        }
    
    # Datasets to propagate a RTS using affine registration of soft tissue 
    # (start with copy of 'Soft_tissue_seg_aff' since nearly identical):
    runID = 'Soft_tissue_RTS_CT_to_MR_aff_using_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    #cfg[runID]['srcRoicolName'] = 'CT_Ser2_tumour'
    
    # Datasets to propagate SEG using bspine registration of soft tissue:
    runID = 'Soft_tissue_SEG_CT_to_MR_def_using_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['regTxName'] = 'bspline'
    cfg[runID]['maxIters'] = 100
    
    # Datasets to propagate RTS using bspine registration of soft tissue:
    runID = 'Soft_tissue_RTS_CT_to_MR_def_using_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_def_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    
    """ 
    As above but without using fiducials.
    
    Note 01/10/21:
        These runs fail to execute to completion since the registrations go
        uttery awry, resulting in empty resampled label images.
    """
    
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_no_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcFidsFpath'] = ''
    cfg[runID]['trgFidsFpath'] = ''
    cfg[runID]['initMethod'] = 'geometry'
    #cfg[runID]['initMethod'] = 'moments'
    
    runID = 'Soft_tissue_RTS_CT_to_MR_aff_no_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    
    runID = 'Soft_tissue_SEG_CT_to_MR_def_no_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['regTxName'] = 'bspline'
    cfg[runID]['maxIters'] = 100
    
    runID = 'Soft_tissue_RTS_CT_to_MR_def_no_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_def_no_fiducials'])
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
        'p2c' : p2c,
        'cwd' : cwd,
        'downloadDir' : downloadDir,
        'inputsDir' : inputsDir,
        'outputsDir' : outputsDir,
        'sampleDroDir' : sampleDroDir,
        'fidsDir' : fidsDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir
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
    
    # Datasets to propagate RTS across two DICOM series with the same slice
    # thicknes:
    """ 
    Force registration to account for different patient position despite same
    FORuid.
    """
    runID = 'NCITA_test_RR2_force_reg'
    cfg[runID] = dict(cfg['NCITA_test_RR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    """ 
    Use a pre-archived DRO rather than image registration to account for 
    different patient position despite same FORuid.
    """
    runID = 'NCITA_test_RR2_use_dro'
    cfg[runID] = dict(cfg['NCITA_test_RR2_force_reg'])
    cfg[runID]['runID'] = runID
    cfg[runID]['useDroForTx'] = True
    
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
    
    # Datasets to propagate RTS from one DICOM series to another with smaller
    # slice thickness:
    """ 
    Force registration to account for different patient position despite same
    FORuid.
    """
    runID = 'NCITA_test_RR4_force_reg'
    cfg[runID] = dict(cfg['NCITA_test_RR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    """ 
    Use a pre-archived DRO rather than image registration to account for 
    different patient position despite same FORuid.
    """
    runID = 'NCITA_test_RR4_use_dro'
    cfg[runID] = dict(cfg['NCITA_test_RR4_force_reg'])
    cfg[runID]['runID'] = runID
    cfg[runID]['useDroForTx'] = True
    
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
    
    # Datasets to propagate SEG across two DICOM series with the same slice
    # thickness:
    """ 
    Force registration to account for different patient position despite same
    FORuid.
    """
    runID = 'NCITA_test_SR2_force_reg'
    cfg[runID] = dict(cfg['NCITA_test_SR2'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    """ 
    Use a pre-archived DRO rather than image registration to account for 
    different patient position despite same FORuid.
    """
    runID = 'NCITA_test_SR2_use_dro'
    cfg[runID] = dict(cfg['NCITA_test_SR2_force_reg'])
    cfg[runID]['runID'] = runID
    cfg[runID]['useDroForTx'] = True
    
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
    
    # Datasets to propagate SEG from one DICOM series to another with smaller
    # slice thickness:
    """ 
    Force registration to account for different patient position despite same
    FORuid.
    """
    runID = 'NCITA_test_SR4_force_reg'
    cfg[runID] = dict(cfg['NCITA_test_SR4'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    """ 
    Use a pre-archived DRO rather than image registration to account for 
    different patient position despite same FORuid.
    """
    runID = 'NCITA_test_SR4_use_dro'
    cfg[runID] = dict(cfg['NCITA_test_SR4_force_reg'])
    cfg[runID]['runID'] = runID
    cfg[runID]['useDroForTx'] = True
    
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
    
    
    """
    Following are datasets for test runs using the ACRIN collection (brain
    images) from TCIA, testing how well the algorithm is able to deal with
    propagating SEG and RTSTRUCT ROI Collections from one imaging perspective 
    (e.g. sagittal, coronal) to another (e.g. axial):
    """
    
    # Datasets to propagate SEG from sagital view to axial:
    runID = 'ACRIN_SEG_MR4_SAG_to_MR12_AX'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.20", 
        'username': "admin",
        'projID' : 'ACRIN', 
        'subjLab' : 'ACRIN-FMISO-Brain-011',
        'srcExpLab' : 'ACRIN-FMISO-Brain-011_MR_4',
        'srcScanID' : '11',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Tumour SAG',
        'srcRoiName' : None,
        'roicolMod' : 'SEG',
        'trgExpLab' : 'ACRIN-FMISO-Brain-011_MR_12',
        'trgScanID' : '3',
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
        'p2c' : p2c,
        'cwd' : cwd,
        'downloadDir' : downloadDir,
        'inputsDir' : inputsDir,
        'outputsDir' : outputsDir,
        'sampleDroDir' : sampleDroDir,
        'fidsDir' : fidsDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir
        }
    
    # Datasets to propagate RTS from coronal view to axial:
    runID = 'ACRIN_RTS_MR4_COR_to_MR12_AX'
    cfg[runID] = dict(cfg['ACRIN_SEG_MR4_SAG_to_MR12_AX'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '10'
    cfg[runID]['srcRoicolName'] = 'Tumour COR'
    cfg[runID]['roicolMod'] = 'RTSTRUCT'
    cfg[runID]['trgScanID'] = '8'
    
    
    """
    Following are datasets for test runs using the Chest Imaging with Clinical 
    and Genomic Correlates Representing a Rural COVID-19 Positive Population 
    collection from TCIA, testing how well the algorithm is able to deal with
    propagating SEG and RTSTRUCT ROI Collections from one imaging perspective 
    (e.g. sagittal, coronal) to another (e.g. axial) within the same imaging
    session:
    """
    
    # Datasets to propagate SEG from axial view to coronal:
    runID = 'COVID_RTS_4_AX_to_80642_COR'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.21", 
        'username': "admin",
        'projID' : 'COVID_19_AR', 
        'subjLab' : 'COVID-19-AR-16406488',
        'srcExpLab' : 'CT-PE-CHEST-63916',
        'srcScanID' : '4',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Lungs',
        'srcRoiName' : None,
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'CT-PE-CHEST-63916',
        'trgScanID' : '80642',
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
        'p2c' : p2c,
        'cwd' : cwd,
        'downloadDir' : downloadDir,
        'inputsDir' : inputsDir,
        'outputsDir' : outputsDir,
        'sampleDroDir' : sampleDroDir,
        'fidsDir' : fidsDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir
        }
    
    # Datasets to propagate RTS from axial view to sagittal:
    runID = 'COVID_RTS_4_AX_to_80643_SAG'
    cfg[runID] = dict(cfg['COVID_RTS_4_AX_to_80642_COR'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '80643'

    # Datasets to propagate SEG from coronal view to sagittal:
    runID = 'COVID_SEG_80642_COR_to_80643_SAG'
    cfg[runID] = dict(cfg['COVID_RTS_4_AX_to_80643_SAG'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcScanID'] = '80642'
    cfg[runID]['roicolMod'] = 'SEG'
    cfg[runID]['trgScanID'] = '80643'
    
    # Datasets to propagate SEG from coronal view to axial:
    runID = 'COVID_SEG_80642_COR_to_4_AX'
    cfg[runID] = dict(cfg['COVID_SEG_80642_COR_to_80643_SAG'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '4'
    
    
    
    """
    Following are datasets for test runs using the Chest Imaging with Clinical 
    and Genomic Correlates Representing a Rural COVID-19 Positive Population 
    collection from TCIA, testing how well the algorithm is able to deal with
    propagating SEG and RTSTRUCT ROI Collections from one imaging perspective 
    (e.g. sagittal, coronal) to another (e.g. axial) across different imaging
    sessions:
    """
    
    # Datasets to propagate SEG from axial view to coronal:
    runID = 'COVID_RTS_5_AX_to_80344_COR'
    cfg[runID] = {
        'runID' : runID,
        # XNAT config settings:
        'url' : "http://10.1.1.21", 
        'username': "admin", 
        'projID' : 'COVID_19_AR', 
        'subjLab' : 'COVID-19-AR-16406502',
        'srcExpLab' : 'CT-PE-CHEST-28055',
        'srcScanID' : '5',
        'srcSlcNum' : None,
        'srcRoicolName' : 'Lungs axial',
        'srcRoiName' : None,
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'CT-CHEST-WO-CONTRAST-99968',
        'trgScanID' : '80344',
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
        'p2c' : p2c,
        'cwd' : cwd,
        'downloadDir' : downloadDir,
        'inputsDir' : inputsDir,
        'outputsDir' : outputsDir,
        'sampleDroDir' : sampleDroDir,
        'fidsDir' : fidsDir,
        'rtsExportDir' : rtsExportDir,
        'segExportDir' : segExportDir,
        'droExportDir' : droExportDir,
        'txExportDir' : txExportDir,
        'imExportDir' : imExportDir,
        'labimExportDir' : labimExportDir,
        'logsExportDir' : logsExportDir,
        'rtsPlotsExportDir' : rtsPlotsExportDir,
        'segPlotsExportDir' : segPlotsExportDir,
        'resPlotsExportDir' : resPlotsExportDir
        }
    
    # Datasets to propagate RTS from axial view to sagittal:
    runID = 'COVID_RTS_5_AX_to_80345_SAG'
    cfg[runID] = dict(cfg['COVID_RTS_5_AX_to_80344_COR'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '80345'

    # Datasets to propagate SEG from coronal view to sagittal:
    runID = 'COVID_SEG_80344_COR_to_8041_SAG'
    cfg[runID] = dict(cfg['COVID_RTS_5_AX_to_80344_COR'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'CT-CHEST-WO-CONTRAST-99968'
    cfg[runID]['srcScanID'] = '80344'
    cfg[runID]['srcRoicolName'] = 'Lungs coronal'
    cfg[runID]['roicolMod'] = 'SEG'
    cfg[runID]['trgExpLab'] = 'CT-PE-CHEST-28055'
    cfg[runID]['trgScanID'] = '8041'
    
    # Datasets to propagate SEG from coronal view to axial:
    runID = 'COVID_SEG_80344_COR_to_8040_COR'
    cfg[runID] = dict(cfg['COVID_SEG_80344_COR_to_8041_SAG'])
    cfg[runID]['runID'] = runID
    cfg[runID]['trgScanID'] = '8040'
    
    #export_dict_to_json(dictionary=cfg, filename='main_params.json',
    #                    exportDir=cfgDir)
    
    # Export each sub-dictionary to a unique JSON file:
    for runID in cfg.keys():
        
        export_dict_to_json(
            dictionary=cfg[runID], filename=f'{runID}.json', exportDir=cfgDir
            )
        
    print(f"\nConfiguration files have been exported to {cfgDir}")



if __name__ == '__main__':
    """
    Run config.py as a script.
    
    Example usage in a console:
    
    cd C:\Code\WP1.3_multiple_modalities\src
    python config.py configs
    """
    
    parser = argparse.ArgumentParser(
        description='Arguments for create_config_files()'
        )
    
    parser.add_argument(
        "cfgDir", 
        help="The directory to export config files"
        )
    
    args = parser.parse_args()
    
    # Run create_config_files():
    create_config_files(args.cfgDir)