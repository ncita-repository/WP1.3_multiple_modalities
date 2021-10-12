# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 10:39:57 2021

@author: ctorti
"""

from importlib import reload

import dicom_tools.dcm_metadata
reload(dicom_tools.dcm_metadata)
import image_tools.attrs_info
reload(image_tools.attrs_info)
#import image_tools.imports
#reload(image_tools.imports)
import general_tools.general
reload(general_tools.general)
import general_tools.geometry
reload(general_tools.geometry)

from dicom_tools.dcm_metadata import get_dcm_uids#, get_roicol_labels
from image_tools.attrs_info import get_im_attrs
#from image_tools.imports import import_im
from io_tools.imports import import_im
from general_tools.general import are_items_equal_to_within_eps
from general_tools.geometry import get_im_extent


def are_inputs_valid(mainParams, trgRoiNames):
    """
    Check if the combination of input parameters are valid.

    Parameters
    ----------
    mainParams : dict
        Dictionary of main parameters for run.
    trgRoiNames : list of strs or None
        List of ROI Names / segment labels of the Target ROI Collection 
        (RTS/SEG), or None if there is no Target ROI Collection.
    
    Returns
    -------
    valid : bool
        True if the combinations of inputs are valid.
    msg : str
        The error message. msg should be '' if valid = True, and non-empty
        otherwise.
    """
    
    mod = mainParams['roicolMod']
    srcRoiName = mainParams['srcRoiName']
    srcSlcNum = mainParams['srcSlcNum']
    
    trgRoiName = mainParams['trgRoiName']
    trgSlcNum = mainParams['trgSlcNum']
    
    """ Useful variables:
    RorS = 'ROI' or 'segment',  
    CorS = 'contours' or 'segmentations', and
    NorL = 'name' or 'label'
    """
    if mod == 'RTSTRUCT':
        RorS = 'ROI'
        CorS = 'contours'
        NorL = 'name'
    else:
        RorS = 'segment'
        CorS = 'segmentations'
        NorL = 'label'
    
    valid = False
    
    msg = ''
    
    if not srcRoiName:
        if srcSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied. \nBut srcSlcNum (= "\
                + f"{srcSlcNum}) is specified (i.e. != None). "
        
        if trgRoiName:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied.\nBut trgRoiName (= "\
                + f"{trgRoiName}) is specified(i.e. != None). "
        
        if trgSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All {CorS} in all {RorS}s "\
                + f"in the Source {mod} are to be copied. \nBut trgSlcNum (= "\
                + f"{trgSlcNum}) is specified (i.e. != None). "
    else:
        if not srcSlcNum and trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => All {CorS} in the {RorS} with"\
                + f" {NorL} '{srcRoiName}' are to be copied.\nBut trgSlcNum "\
                + f"(= {trgSlcNum}) is specified. "
        
        if srcSlcNum and not trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => The {CorS}(s) for the "\
                + f"specified DICOM slice in the {RorS} with {NorL} "\
                + f"'{srcRoiName}' are to be copied.\nBut trgSlcNum (= "\
                + f"{trgSlcNum}) is not specified. "
        
        if trgRoiNames != None:
            T = len(trgRoiNames)
            
            if T > 1 and not trgRoiName:
                msg += f"srcRoiName (= {srcRoiName}) has been specified but "\
                    + f"trgRoiName (= {trgRoiName}) has not, and there are "\
                    + f"{T} {RorS}s in the Target ROI Collection. "
    
    if msg:
        msg += "This is not a valid combination of inputs."
        #raise Exception(msg)
    else:
        valid = True
    
    return valid, msg

def which_use_case(config):
    """
    Determine which of 5 Use Cases applies.
    
    Various data is imported based on file paths and directories stored in the
    config dictionary.
    
    Parameters
    ----------
    config : dict
        Dictionary of parameters for run.
    params : DataDownloader Object, optional
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs). The default value is None.
    
    Returns
    -------
    useCaseThatApplies : str
        String that denotes the use case that applies.
    useCaseToApply : str
        String that denotes the use case that will be applied.
    
    Note
    ----
    useCaseToApply will be the same as useCaseThatApplies unless forceReg is 
    True and useCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    srcDcmDir = config['srcDicomDir']
    
    trgDcmDir = config['trgDicomDir']
    trgSlcNum = config['trgSlcNum']
    
    forceReg = config['forceReg']
    
    p2c = config['p2c']
    
    # Get the DICOM UIDs:
    srcStudyUID, srcSeriesUID, srcFORUID,\
    srcSOPUIDs = get_dcm_uids(dicomDir=srcDcmDir)
    
    trgStudyUID, trgSeriesUID, trgFORUID,\
    trgSOPUIDs = get_dcm_uids(dicomDir=trgDcmDir)
        
    if trgSlcNum:
        if trgSlcNum > len(trgSOPUIDs) - 1:
            msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                  + f'the number of Target DICOMs ({len(trgSOPUIDs)}).'
            raise Exception(msg)
    
    # Get the Image Attributes for Source and Target:
    srcSize, srcSpacing, srcST, srcIPPs, srcDirs,\
    warnings = get_im_attrs(dicomDir=srcDcmDir, package='sitk')
    
    trgSize, trgSpacing, trgST, trgIPPs, trgDirs,\
    warnings = get_im_attrs(dicomDir=trgDcmDir, package='sitk')
    
    
    # Are the Source and Target part of the same study, series and do they 
    # have the same SOP UIDs?
    if srcStudyUID == trgStudyUID:
        sameStudy = True
    else:
        sameStudy = False
        
    if srcSeriesUID == trgSeriesUID:
        sameSeries = True
    else:
        sameSeries = False
        
    if srcFORUID == trgFORUID:
        sameFOR = True
    else:
        sameFOR = False
        
    if srcSOPUIDs == trgSOPUIDs:
        sameSOPs = True
    else:
        sameSOPs = False
    
    # Are the Source and Target spacings, directions and origins equal to  
    # within 1e-06?
    sameSpacings = are_items_equal_to_within_eps(srcSpacing, trgSpacing)
    sameDirs = are_items_equal_to_within_eps(srcDirs, trgDirs)
    sameOrigins = are_items_equal_to_within_eps(srcIPPs[0], trgIPPs[0])
    sameST = are_items_equal_to_within_eps([srcST], [trgST])
    
    if srcSize == trgSize:
        sameSize = True
    else:
        sameSize = False
    
    # Check which case follows:
    if sameFOR:
        if sameDirs:
            """ 
            Comment 11/02 #1:
                
            It's not entirely necessary to resample the Source image to the 
            Target domain just because their origins differ. It depends on
            whether the slice in Source coincide with slices in Target.
            
            Even if the slices don't coincide, an approximation can be made as
            follows.  When dealing contours, the in-plane coordinates don't 
            need to be changed, and the destination slice can be determined by 
            matching on the out-of-plane (z) coordinates as best as possible. 
            When dealing with segmentations, the destination slice can be 
            determined by matching on the z components of the IPPs as best as
            possible.  Once the destination slice number is known, the 
            necessary shift in pixels for the in-plane components can be 
            applied accordingly. 
            
            Perhaps the approximation can be justified if the difference 
            between the planes in Source and Target are close to an integer.
            If the difference is large, the Source image should be resampled to
            the Target image.  For now I will treat UseCase 2 data sets as 
            UseCase 3 if their origins are not the same. Further complexity
            (taking into account the error in the above approximation) can be
            considered for future work. 
            """
            if sameSpacings and sameOrigins:
            #if sameSpacings:
                if sameSeries:
                    useCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(trgSlcNum, int):
                        useCaseThatApplies = '2a' # Direct copy
                    else:
                        useCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ 
                Comment 11/02 #2:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. 
                """
                if isinstance(trgSlcNum, int):
                    useCaseThatApplies = '3a' # Direct copy
                else:
                    useCaseThatApplies = '3b' # Relationship-preserving copy
        else:
            if isinstance(trgSlcNum, int):
                useCaseThatApplies = '4a' # Direct copy
            else:
                useCaseThatApplies = '4b' # Relationship-preserving copy
    else:
        if isinstance(trgSlcNum, int):
            useCaseThatApplies = '5a' # Direct copy
        else:
            useCaseThatApplies = '5b' # Relationship-preserving copy
    
    if p2c:
        # Print comparisons to the console:
        srcIm = import_im(srcDcmDir)
        trgIm = import_im(trgDcmDir)
        
        # Get the image extents
        srcImExtent = get_im_extent(srcIm)
        trgImExtent = get_im_extent(trgIm)
        
        # Are the image extents equal to within 1e-06?
        sameExtents = are_items_equal_to_within_eps(srcImExtent, trgImExtent)
    
        #CompareImageAttributes(srcDcmDir, trgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running which_use_case():')
        print('\n\n', '-'*120)
        print('The Source and Target are in the/have:')
        
        if sameStudy:
            print('   * same study')
        else:
            print('   * different studies')
        
        if sameFOR:
            print('   * same Frame of Reference')
        else:
            print('   * different Frames of Reference')
            
        if sameSeries:
            print('   * same series')
        else:
            print('   * different series')
            
        if sameSOPs:
            print('   * same DICOMs')
        else:
            print('   * different DICOMs')
        
        if sameSize:
            print(f'   * same image size\n      Source/Target: {srcSize}')
        else:
            print(f'   * different image sizes:\n      Source: {srcSize}',
                  f'\n      Target: {trgSize}')
        
        if sameSpacings:
            print(f'   * same voxel sizes\n      Source/Target: {srcSpacing}')
        else:
            print(f'   * different voxel sizes:\n      Source: {srcSpacing}',
                  f'\n      Target: {trgSpacing}')
        
        if sameST:
            print(f'   * same slice thickness\n      Source/Target: {srcST}')
        else:
            print(f'   * different slice thickness:\n      Source: {srcST}',
                  f'\n      Target: {trgST}')
        
        if sameExtents:
            print(f'   * same image extents\n      Source/Target: {srcImExtent}')
        else:
            print(f'   * different image extents:\n      Source: {srcImExtent}',
                  f'\n      Target: {trgImExtent}')
        
        if sameOrigins:
            print(f'   * same origin\n      Source/Target: {srcIPPs[0]}')
        else:
            print(f'   * different origins:\n      Source: {srcIPPs[0]}',
                  f'\n      Target: {trgIPPs[0]}')
            
        if sameDirs:
            print('   * same patient orientation\n',
                  f'      Source/Target: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'                      {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'                      {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]')
        else:
            print('   * different patient orientations:\n',
                  f'      Source: [{srcDirs[0]}, {srcDirs[1]}, {srcDirs[2]},\n',
                  f'               {srcDirs[3]}, {srcDirs[4]}, {srcDirs[5]},\n',
                  f'               {srcDirs[6]}, {srcDirs[7]}, {srcDirs[8]}]\n',
                  f'      Target: [{trgDirs[0]}, {trgDirs[1]}, {trgDirs[2]},\n',
                  f'               {trgDirs[3]}, {trgDirs[4]}, {trgDirs[5]},\n',
                  f'               {trgDirs[6]}, {trgDirs[7]}, {trgDirs[8]}]')
    
    if True:#p2c:
        #PrintTitle(f'Case {UseCase} applies.')
        if forceReg and useCaseThatApplies in ['3a', '3b', '4a', '4b']:
            useCaseToApply = useCaseThatApplies.replace('3', '5').replace('4', '5')
            
            msg = f'\n*UseCase {useCaseThatApplies} applies but since '\
                  + f'forceReg = {forceReg}, UseCase {useCaseToApply} will '\
                  + 'be applied.'
            print(msg)
        else:
            useCaseToApply = useCaseThatApplies
            
            print(f'\n*UseCase {useCaseThatApplies} applies.')
            
        print('-'*120)
        
    return useCaseThatApplies, useCaseToApply