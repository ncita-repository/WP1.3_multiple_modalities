# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:58:18 2021

@author: ctorti
"""


import os
from pathlib import Path
import time
from io_tools.exports import export_dict_to_json
#from io_tools.imports import import_dict_from_json


# Current working directory:
cwd = os.getcwd()

# Default download directory for XNAT data:
ddd = os.path.join(Path.home(), "Downloads", "xnat_downloads")

# Current date:
currentDate = time.strftime("%Y-%m-%d", time.gmtime())

def create_general_params(cfgDir):
    """
    Create general parameters.
    
    The file general_params.json will be exported to a directory named "config"
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
    #exportTrgRoicol = False
    exportTrgRoicol = True
    
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

def create_xnat_cfg(cfgDir):
    """
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

def create_src_xnat_params(cfgDir):
    """
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

def create_trg_xnat_params(cfgDir):
    """
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

def create_res_params(cfgDir):
    """
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