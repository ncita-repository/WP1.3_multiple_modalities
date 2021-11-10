# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:07:09 2021

@author: ctorti
"""



"""
This is a tool to define directories and filepaths that will be assumed to be
provided by the front-end (OHIF-Viewer) and passed to the back-end.  The
purpose is to create a JSON file of a dictionary of dictionaries that define
a number of different use case scenarios.

For now only two principle sets of data will be used, but with a number of 
sub-uses that covers a broader range of possibilities.

Refer to the "ROI Copying Design Spec" for further details:
https://docs.google.com/document/d/1VT7M3-9OvC2dxDR3gvGUFs9OPuy1lDrdGWtOObQNP7M/edit#heading=h.moz6ya56wvsq

and to "Testing of RTS/SEG copying algorithm with NCITA_TEST data":
https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g/edit#

The two datasets that link to the NCITA_TEST runs are RR3 and RR11.  Since the
front end will deal with "simple"/"direct" copies of entities (contours or
segmentations), the back-end will only need to deal with cases that require
resampling (i.e. different voxel resolutions) or image registration (different
frames of reference) between the source ("src") and target ("trg") datasets.  
Hence only tests with the prefix "RR" (Rts Relationship-preserving) or "SR"
(Seg Relationship-preserving) will need to be considered within the full tests
covered in the "Testing of RTS/SEG copying algorithm with NCITA_TEST data" 
document.

The keys of the main dictionary ("runID"s) will be of the form:

"{mod}R{N}_{object}"

where mod = "R" or "S" for RTSTRUCT or SEG
      N = "3" or "11" (for the time being)
      object = "contour" or "segmentation" or "roi" or "segment" or "roicol"

Hence all three scales are covered: a single contour or segmentation, a
collection of contours or segmentations within an entire ROI or segment, and
a collection of ROIs or segments within an ROI Collection.

In the stand-alone version of the code, the DataDownloader class is used to
download data from XNAT using XNAT REST APIs (add_dirs_and_fpaths()), and the 
paths of the downloaded data is added to cfgDict.  Then the DataImporter class
is used to import the data based on the paths in cfgDict.

Since downloading will be by-passed for the OHIF integration, the data has been
put in src/inputs/scans within sub-directories with names that take the form:
    {srcExpLab}_Scan{srcScanID} for source data, and
    {trgExpLab}_Scan{trgScanID} for target data
    
    e.g. Session2_Scan2

Likewise the RTSTRUCTs and SEGs are in src/inputs/roicols with filenames of the
form:
    {mod}_YYYYMMDD_HHMMSS_{srcExpLab}_{roicolName}
    
    e.g. RTSTRUCT_20210131_214838_Session4_Nasal_cavity
    
The paths will be added manually here.
"""



import os
#from pathlib import Path
#import time
#import argparse
from io_tools.exports import export_dict_to_json
from io_tools.imports import import_dict_from_json



def create_cfgDict_json():
    """
    Create a dictionary of dictionaries storing paths to files and directories
    that will mimic that required of the front-end to be passed to the back-end.
    
    Returns
    -------
    None.
    
    Note
    ----
    The int values for srcSlcNum are zero-indexed.
    """
    
    # Try to get the current working directory from an environmnt variable:
    cwd = os.getenv('workdir')
    #print(f'os.getenv (in fetch_config.py) = {cwd}')
    if cwd == None:
        cwd = os.getcwd()
    #print(f'cwd (in fetch_config.py) = {cwd}')
    
    globalVars_fpath = os.path.join(cwd, 'global_variables.json')
    
    # Get the global variables:
    globalVars = import_dict_from_json(globalVars_fpath)
    
    cfgDict = {}
    
    """ Runs that use image resampling """
    
    """ RTSTRUCTs """
    
    """ 
    Propagate the contour on slice 7 from ROI 'Left eye CT"' in ROI Collection
    'Left eye CT' from 'Session2', scan 2 to scan 5.
    """
    runID = 'RR3_contour'
    cfgDict[runID] = {
        'runID' : runID,
        'srcExpLab' : 'Session2',
        'srcScanID' : '2',
        'srcSlcNum' : 7,
        'srcRoicolName' : 'Left eye CT',
        'srcRoiName' : 'Left eye CT',
        'roicolMod' : 'RTSTRUCT',
        'trgExpLab' : 'Session2',
        'trgScanID' : '4',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFpath' : None,
        'trgFidsFpath' : None, 
        'forceReg' : globalVars['forceReg'],
        'useDroForTx' : globalVars['useDroForTx'],
        'regTxName' : globalVars['regTxName'],
        'initMethod' : globalVars['initMethod'],
        'maxIters' : globalVars['maxIters'],
        'applyPreResBlur' : globalVars['applyPreResBlur'],
        'preResVar' : globalVars['preResVar'],
        'resInterp' : globalVars['resInterp'],
        'applyPostResBlur' : globalVars['applyPostResBlur'],
        'postResVar' : globalVars['postResVar'],
        'addToRoicolLab' : globalVars['addToRoicolLab'],
        'exportRoicol' : globalVars['exportRoicol'],
        'exportDro' : globalVars['exportDro'],
        'exportTx' : globalVars['exportTx'],
        'exportIm' : globalVars['exportIm'],
        'exportLabim' : globalVars['exportLabim'],
        'exportPlots' : globalVars['exportPlots'],
        'exportLogs' : globalVars['exportLogs'],
        #'whichSrcRoicol' : whichSrcRoicol,
        #'addToRoicolLab' : addToRoicolLab,
        'uploadDro' : globalVars['uploadDro'],
        'overwriteDro' : globalVars['overwriteDro'],
        'p2c' : globalVars['p2c'],
        #'cwd' : globalVars['cwd'],
        'cwd' : cwd,
        #'downloadDir' : downloadDir,
        'inputsDir' : globalVars['inputsDir'],
        'outputsDir' : globalVars['outputsDir'],
        'sampleDroDir' : globalVars['sampleDroDir'],
        'fidsDir' : globalVars['fidsDir'],
        'srcDicomDir' : os.path.join(
            cwd, 'inputs', 'scans', 'Session2_Scan2'
            ),
        'trgDicomDir' : os.path.join(
            cwd, 'inputs', 'scans', 'Session2_Scan5'
            ),
        'srcRoicolFpath' : os.path.join(
            cwd, 'inputs', 'roicols', 
            'RTSTRUCT_20210201_000732_Session2_Left_eye_CT.dcm'
            ),
        'trgRoicolFpath' : None,
        'droDir' : os.path.join(cwd, 'inputs', 'dro'),
        'rtsExportDir' : globalVars['rtsExportDir'],
        'segExportDir' : globalVars['segExportDir'],
        'droExportDir' : globalVars['droExportDir'],
        'txExportDir' : globalVars['txExportDir'],
        'imExportDir' : globalVars['imExportDir'],
        'labimExportDir' : globalVars['labimExportDir'],
        'logsExportDir' : globalVars['logsExportDir'],
        'rtsPlotsExportDir' : globalVars['rtsPlotsExportDir'],
        'segPlotsExportDir' : globalVars['segPlotsExportDir'],
        'resPlotsExportDir' : globalVars['resPlotsExportDir']
        }
    
    """ 
    Propagate the entire ROI 'Left eye CT' in ROI Collection 'Left eye CT'
    from 'Session2', scan 2 to scan 5.
    """
    runID = 'RR3_roi'
    cfgDict[runID] = dict(cfgDict['RR3_contour']) # copy existing and modify:
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcSlcNum'] = None
    
    """ 
    Propagate the entire ROI Collection 'Left eye CT' from 'Session2', scan 2 
    to scan 5.
    
    Note:
        This should yield the same results as 'RR3_roi' since there is only 1
        ROI in this ROI Collection.
    """
    runID = 'RR3_roicol'
    cfgDict[runID] = dict(cfgDict['RR3_roi'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcRoiName'] = None
    
    
    """ SEGs """
    
    """ 
    Propagate the segmentation on slice 7 from segment 'Left eye CT' in ROI 
    Collection 'Left eye CT' from 'Session2', scan 2 to scan 5.
    """
    runID = 'SR3_segmentation'
    cfgDict[runID] = dict(cfgDict['RR3_contour'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['roicolMod'] = 'SEG'
    cfgDict[runID]['srcRoicolFpath'] = os.path.join(
        cwd, 'inputs', 'roicols',
        'SEG_20210213_140215_Session2_Left_eye_CT.dcm'
        )
    
    """ 
    Propagate the entire segment 'Left eye CT' in ROI Collection 'Left eye CT'
    from 'Session2', scan 2 to scan 5.
    """
    runID = 'SR3_segment'
    cfgDict[runID] = dict(cfgDict['SR3_segmentation'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcSlcNum'] = None
    
    """ 
    Propagate the entire ROI Collection 'Left eye CT' from 'Session2', scan 2 
    to scan 5.
    
    Note:
        This should yield the same results as 'RR3_segment' since there is only
        1 ROI in this ROI Collection.
    """
    runID = 'SR3_roicol'
    cfgDict[runID] = dict(cfgDict['SR3_segment'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcRoiName'] = None
                        
    
    
    """ Runs that use image registration """
    
    """ RTSTRUCTs """
    
    """ 
    Propagate the contour on slice 18 from ROI 'Nasal cavity (right)' in 
    ROI Collection 'Nasal cavity' from 'Session4', scan 13 to 'Session3', 
    scan 6.
    """
    runID = 'RR11_contour'
    cfgDict[runID] = dict(cfgDict['RR3_contour'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcExpLab'] = 'Session4'
    cfgDict[runID]['srcScanID'] = '13'
    cfgDict[runID]['srcSlcNum'] = 18
    cfgDict[runID]['srcRoicolName'] = 'Nasal cavity'
    cfgDict[runID]['srcRoiName'] = 'Nasal cavity (right)'
    cfgDict[runID]['trgExpLab'] = 'Session3'
    cfgDict[runID]['trgScanID'] = '6'
    cfgDict[runID]['srcDicomDir'] = os.path.join(
        cwd, 'inputs', 'scans', 'Session4_Scan13'
        )
    cfgDict[runID]['trgDicomDir'] = os.path.join(
        cwd, 'inputs', 'scans', 'Session3_Scan6'
        )
    cfgDict[runID]['srcRoicolFpath'] = os.path.join(
        cwd, 'inputs', 'roicols',
        'RTSTRUCT_20210131_214838_Session4_Nasal_cavity.dcm'
        )
    
    """ 
    Propagate the entire ROI 'Nasal cavity (right)' in ROI Collection 
    'Nasal cavity' from 'Session4', scan 13 to 'Session3', scan 6.
    """
    runID = 'RR11_roi'
    cfgDict[runID] = dict(cfgDict['RR11_contour'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcSlcNum'] = None
    
    """ 
    Propagate the entire ROI Collection 'Nasal cavity' from 'Session4', scan 13
    to 'Session3', scan 6.
    """
    runID = 'RR11_roicol'
    cfgDict[runID] = dict(cfgDict['RR11_roi'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcRoiName'] = None
    
    
    """ SEGs """
    
    """ 
    Propagate the segmentation on slice 18 from ROI 'Nasal cavity (right)' in 
    ROI Collection 'Nasal cavity' from 'Session4', scan 13 to 'Session3', 
    scan 6.
    """
    runID = 'SR11_segmentation'
    cfgDict[runID] = dict(cfgDict['RR11_contour'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['roicolMod'] = 'SEG'
    cfgDict[runID]['srcRoicolFpath'] = os.path.join(
        cwd, 'inputs', 'roicols',
        'SEG_20210213_152014_Session4_Nasal_cavity.dcm'
        )
    
    """ 
    Propagate the entire segment 'Nasal cavity (right)' in ROI Collection 
    'Nasal cavity' from 'Session4', scan 13 to 'Session3', scan 6.
    """
    runID = 'SR11_segment'
    cfgDict[runID] = dict(cfgDict['SR11_segmentation'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcSlcNum'] = None
    
    """ 
    Propagate the entire ROI Collection 'Nasal cavity' from 'Session4', scan 13
    to 'Session3', scan 6.
    """
    runID = 'SR11_roicol'
    cfgDict[runID] = dict(cfgDict['SR11_segment'])
    cfgDict[runID]['runID'] = runID
    cfgDict[runID]['srcRoiName'] = None
    
    
    # Export the dictionary to a JSON file:
    filename = 'cfgDict_ohif.json'
    export_dict_to_json(
        dictionary=cfgDict, filename=filename, exportDir=cwd
        )
    
    filepath = os.path.join(cwd, filename)
    print(f"\nRun configs file has been exported to {filepath}")



if __name__ == '__main__':
    """
    Run create_cfgDict_json.py as a script.
    
    Example usage in a console:
    
    cd C:\Code\WP1.3_multiple_modalities\src
    python create_cfgDict_json.py
    """
    
    # Run create_cfgDict_json():
    create_cfgDict_json()