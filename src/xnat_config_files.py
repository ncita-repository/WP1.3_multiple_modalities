# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:58:18 2021

@author: ctorti
"""


import os
#from pathlib import Path
#import time
#import argparse
from io_tools.exports import export_dict_to_json
#from io_tools.imports import import_dict_from_json


def create_xnat_config_files():
    """
    Generate a JSON files for various XNAT configurations.
    
    Files will be exported to src/xnat_configs/
    
    Parameters
    ----------
    None.
    
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
    #currentDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    # Get current working directory:
    cwd = os.getcwd()
    
    # Directory where XNAT config files will be exported to:
    cfgDir = os.path.join(cwd, 'xnat_configs')
    
    # Initialise dictionary (of dictionaries) to store all configurations:
    cfg = {}
    
    
    
    
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
        # XNAT parameters:
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
        'trgRoiName' : None
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
    Following are tests from the same collection but covering only two pairs of
    DICOM scans - one that requires resampling, and the other image 
    registration to propagate the entity.
    
    For each pair, 3 entities are considered - the contour/segmentation, an
    entire ROI/segment, and an entire ROI Collection.
    """
    
    """ Runs that use image resampling """
    
    """ RTSTRUCTs """
    
    # Propagate the contour on slice 7 from ROI 'Left eye CT"' in ROI Collection
    # 'Left eye CT' from 'Session2', scan 2 to scan 5.
    runID = 'RR3_contour'
    cfg[runID] = {
        'runID' : runID,
        # XNAT parameters:
        'url' : 'http://10.1.1.20', 
        'username': 'admin',
        'projID' : 'NCITA_TEST', 
        'subjLab' : 'TCGA-BB-A5HY',
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
        'trgRoiName' : None
        }
    
    # Propagate the entire ROI 'Left eye CT' in ROI Collection 'Left eye CT'
    # from 'Session2', scan 2 to scan 5.
    runID = 'RR3_roi'
    cfg[runID] = dict(cfg['RR3_contour']) # copy existing and modify:
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = None
    
    # Propagate the entire ROI Collection 'Left eye CT' from 'Session2', scan 2 
    # to scan 5.
    """
    Note:
        This should yield the same results as 'RR3_roi' since there is only 1
        ROI in this ROI Collection.
    """
    runID = 'RR3_roicol'
    cfg[runID] = dict(cfg['RR3_roi'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcRoiName'] = None
    
    
    """ SEGs """
    
    # Propagate the segmentation on slice 7 from segment 'Left eye CT' in ROI 
    # Collection 'Left eye CT' from 'Session2', scan 2 to scan 5.
    runID = 'SR3_segmentation'
    cfg[runID] = dict(cfg['RR3_contour'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
     
    # Propagate the entire segment 'Left eye CT' in ROI Collection 'Left eye CT'
    # from 'Session2', scan 2 to scan 5.
    runID = 'SR3_segment'
    cfg[runID] = dict(cfg['SR3_segmentation'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = None
    
    # Propagate the entire ROI Collection 'Left eye CT' from 'Session2', scan 2 
    # to scan 5.
    """
    Note:
        This should yield the same results as 'RR3_segment' since there is only
        1 ROI in this ROI Collection.
    """
    runID = 'SR3_roicol'
    cfg[runID] = dict(cfg['SR3_segment'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcRoiName'] = None
                        
    
    
    """ Runs that use image registration """
    
    """ RTSTRUCTs """
    
    # Propagate the contour on slice 18 from ROI 'Nasal cavity (right)' in 
    # ROI Collection 'Nasal cavity' from 'Session4', scan 13 to 'Session3', 
    # scan 6.
    runID = 'RR11_contour'
    cfg[runID] = dict(cfg['RR3_contour'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcExpLab'] = 'Session4'
    cfg[runID]['srcScanID'] = '13'
    cfg[runID]['srcSlcNum'] = 18
    cfg[runID]['srcRoicolName'] = 'Nasal cavity'
    cfg[runID]['srcRoiName'] = 'Nasal cavity (right)'
    cfg[runID]['trgExpLab'] = 'Session3'
    cfg[runID]['trgScanID'] = '6'
    
    # Propagate the entire ROI 'Nasal cavity (right)' in ROI Collection 
    # 'Nasal cavity' from 'Session4', scan 13 to 'Session3', scan 6.
    runID = 'RR11_roi'
    cfg[runID] = dict(cfg['RR11_contour'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = None
    
    # Propagate the entire ROI Collection 'Nasal cavity' from 'Session4', scan 13
    # to 'Session3', scan 6.
    runID = 'RR11_roicol'
    cfg[runID] = dict(cfg['RR11_roi'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcRoiName'] = None
    
    
    """ SEGs """
    
    # Propagate the segmentation on slice 18 from ROI 'Nasal cavity (right)' in 
    # ROI Collection 'Nasal cavity' from 'Session4', scan 13 to 'Session3', 
    # scan 6.
    runID = 'SR11_segmentation'
    cfg[runID] = dict(cfg['RR11_contour'])
    cfg[runID]['runID'] = runID
    cfg[runID]['roicolMod'] = 'SEG'
    
    # Propagate the entire segment 'Nasal cavity (right)' in ROI Collection 
    # 'Nasal cavity' from 'Session4', scan 13 to 'Session3', scan 6.
    runID = 'SR11_segment'
    cfg[runID] = dict(cfg['SR11_segmentation'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcSlcNum'] = None
     
    # Propagate the entire ROI Collection 'Nasal cavity' from 'Session4', scan 13
    # to 'Session3', scan 6.
    runID = 'SR11_roicol'
    cfg[runID] = dict(cfg['SR11_segment'])
    cfg[runID]['runID'] = runID
    cfg[runID]['srcRoiName'] = None
    
    
    """
    Following are datasets for test runs using the ACRIN collection (brain
    images) from TCIA, testing how well the algorithm is able to deal with
    propagating SEG and RTSTRUCT ROI Collections from one imaging perspective 
    (e.g. sagittal, coronal) to another (e.g. axial).
    
    Data source:
    https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948305#4613f7d5a2514c2195445afde60ce068
    
    Collection ID = ACRIN-FMISO-Brain
    Subject ID = ACRIN-FMISO-Brain-011
    """
    
    # Datasets to propagate SEG from sagital view to axial:
    runID = 'ACRIN_SEG_MR4_SAG_to_MR12_AX'
    cfg[runID] = {
        'runID' : runID,
        # XNAT parameters:
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
        'trgRoiName' : None
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
    session.
    
    Data source:
    https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=70226443
    
    Collection ID = COVID-19-AR
    Subject ID = COVID-19-AR-16406488
    """
    
    # Datasets to propagate SEG from axial view to coronal:
    runID = 'COVID_RTS_4_AX_to_80642_COR'
    cfg[runID] = {
        'runID' : runID,
        # XNAT parameters:
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
        'trgRoiName' : None
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
    sessions.
    
    Data source:
    https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=70226443
    
    Collection ID = COVID-19-AR
    Subject ID = COVID-19-AR-16406488
    """
    
    # Datasets to propagate SEG from axial view to coronal:
    runID = 'COVID_RTS_5_AX_to_80344_COR'
    cfg[runID] = {
        'runID' : runID,
        # XNAT parameters:
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
        'trgRoiName' : None
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
    
    
    
    """
    Datasets that involve affine or bspline registration with or without
    fiducials.
    
    Data source:
    https://wiki.cancerimagingarchive.net/display/Public/Soft-tissue-Sarcoma
    
    Collection ID = Soft-tissue-Sarcoma
    Subject ID = STS_004
    """
    
    # Datasets to propagate a SEG using affine registration of soft tissue:
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'
    cfg[runID] = {
        'runID' : runID,
        # XNAT parameters:
        'url' : "http://10.1.1.20", 
        'username': "admin",
        'projID' : 'Soft-tissue-Sarcoma', 
        'subjLab' : 'STS_004',
        'srcExpLab' : 'STS_004_PETCT',
        'srcScanID' : '2',
        'srcSlcNum' : None,
        'srcRoicolName' : 'CT_Ser2_tumour',
        'srcRoiName' : None,
        'roicolMod' : 'SEG',
        'trgExpLab' : 'STS_004_MR',
        'trgScanID' : '501',
        'trgSlcNum' : None,
        'trgRoicolName' : None,
        'trgRoiName' : None,
        # Resampling / Registration / Transformation settings:
        'srcFidsFname' : 'STS_004_CT_fiducials.txt',
        'trgFidsFname' : 'STS_004_MR_fiducials.txt',
        'initMethod' : 'landmarks'
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
    As above but forcing registration.
    """
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_using_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_RTS_CT_to_MR_aff_using_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_RTS_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_SEG_CT_to_MR_def_using_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_def_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_RTS_CT_to_MR_def_using_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_RTS_CT_to_MR_def_using_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    
    """ 
    As above but without using fiducials.
    
    Note 01/10/21:
        These runs fail to execute to completion since the registrations go
        uttery awry, resulting in empty resampled label images.
    """
    
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_no_fiducials'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_using_fiducials'])
    cfg[runID]['runID'] = runID
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
    As above but forcing registration.
    """
    
    runID = 'Soft_tissue_SEG_CT_to_MR_aff_no_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_aff_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_RTS_CT_to_MR_aff_no_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_RTS_CT_to_MR_aff_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_SEG_CT_to_MR_def_no_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_SEG_CT_to_MR_def_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    runID = 'Soft_tissue_RTS_CT_to_MR_def_no_fiducials_force_reg'
    cfg[runID] = dict(cfg['Soft_tissue_RTS_CT_to_MR_def_no_fiducials'])
    cfg[runID]['runID'] = runID
    cfg[runID]['forceReg'] = True
    
    
    
    # Export each sub-dictionary to a unique JSON file:
    for runID in cfg.keys():
        export_dict_to_json(
            dictionary=cfg[runID], filename=f'{runID}.json', exportDir=cfgDir
            )
        
    print(f"\nConfiguration files have been exported to {cfgDir}")


if __name__ == '__main__':
    """
    Run create_xnat_config_files.py as a script.
    
    Example usage in a console:
    
    python xnat_config_files.py
    """
    
    create_xnat_config_files()