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
import image_tools.imports
reload(image_tools.imports)
import general_tools.general
reload(general_tools.general)
import general_tools.geometry
reload(general_tools.geometry)

from dicom_tools.dcm_metadata import get_dcm_uids, get_roicol_labels
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

def are_inputs_valid_210810(params, trgRoiNames):
    """
    Check if the combination of input parameters are valid.

    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
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
    
    mod = params.srcXnatParams['roicolMod']
    srcRoiName = params.srcXnatParams['roiName']
    srcSlcNum = params.srcXnatParams['slcNum']
    
    trgRoiName = params.trgXnatParams['roiName']
    trgSlcNum = params.trgXnatParams['slcNum']
    
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

def are_inputs_valid_OLD(srcRoi, srcRoiName, srcSlcNum, trgRoi, trgRoiName, 
                     trgSlcNum, p2c=False):
    """
    Check if the combination of inputs are valid.

    Parameters
    ----------
    srcRoi : Pydicom object
        The Source RTS/SEG object.
    srcRoiName : string
        The Source ROIName or SegmentLabel.
    srcSlcNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If srcSlcNum = None, a relationship-preserving copy will be made.
    trgRoi : Pydicom object or None
        The Target RTS/SEG object if applicable; None otherwise.
    trgRoiName : string or None
        The Target ROIName or SegmentLabel if applicable; None otherwise.
    trgSlcNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If trgSlcNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    p2c : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    IsValid : boolean
        True if the combinations of inputs are valid.
    """
    
    msg = ''
              
    #if srcRoiName == None:
    if not srcRoiName:
        #if srcSlcNum != None:
        if srcSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut srcSlcNum (= {srcSlcNum}) is "\
                   + "specified (i.e. != None). "
        
        #if trgRoiName != None:
        if trgRoiName:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut trgRoiName (= {trgRoiName}) is specified "\
                   + "(i.e. != None). "
        
        #if trgSlcNum != None:
        if trgSlcNum:
            msg += f"srcRoiName = {srcRoiName} => All contours/segmentations "\
                   + "in all ROIs/segments in the Source RTS/SEG are to be "\
                   + f"copied. \nBut trgSlcNum (= {trgSlcNum}) is "\
                   + "specified (i.e. != None). "
        
    else:
        #if srcSlcNum == None and trgSlcNum != None:
        if not srcSlcNum and trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => All contours/segmentations"\
                   + f"in the ROI/segment with name/label '{srcRoiName}' are "\
                   + f"to be copied. \nBut trgSlcNum (= {trgSlcNum}) is "\
                   + "specified. "
        
        #if srcSlcNum != None and trgSlcNum == None:
        if srcSlcNum and not trgSlcNum:
            msg += f"srcSlcNum = {srcSlcNum} => The contour(s)/"\
                   + "segmentation(s) for the specified DICOM slice in the "\
                   + f"ROI/segment with name/label '{srcRoiName}' are to be "\
                   + f"copied. \nBut trgSlcNum (= {trgSlcNum}) is not "\
                   + "specified. "
        
        if trgRoi != None:
            """ Get the names/labels of all ROIs/segments in trgRoi: """
            trgRoiNames = get_roicol_labels(trgRoi)
            
            T = len(trgRoiNames)
            
            #if T > 1 and trgRoiName == None:
            if T > 1 and not trgRoiName:
                msg += f"srcRoiName (= {srcRoiName}) has been specified but "\
                       + f"trgRoiName (= {trgRoiName}) has not been specified"\
                       + f" and there are {T} ROIs/segments in trgRoi. "
        
    if msg:
        msg += "This is not a valid combination of inputs."
            
        raise Exception(msg)
    
    return True

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

def which_use_case_210810(srcXnatParams, trgXnatParams, resParams, genParams):
    """
    Determine which of 5 Use Cases applies.
    
    Parameters
    ----------
    srcXnatParams : dict
        Dictionary of XNAT parameters for Source.
    trgXnatParams : dict
        Dictionary of XNAT parameters for Target.
    resParams : dict
        Dictionary of general parameters.
    genParams : dict
        Dictionary containing resampling, registration and transformation
        related parameters.
    
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
    
    srcDcmDir = srcXnatParams['dicomDir']
    
    trgDcmDir = trgXnatParams['dicomDir']
    trgSlcNum = trgXnatParams['slcNum']
    
    forceReg = resParams['forceReg']
    
    p2c = genParams['p2c']
    
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
            Comment 11/02:
                
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
                Comment 11/02:
                    
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

def which_use_case_Pre_210810(srcSlcNum, trgSlcNum, srcDcmDir, trgDcmDir, forceReg, 
                   p2c=False):
    """
    Determine which of 5 Use Cases applies.
    
    Parameters
    ----------
    srcSlcNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If srcSlcNum = None, a relationship-preserving copy will be made.
    trgSlcNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If trgSlcNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    srcDcmDir : string
        Directory containing the Source DICOMs.
    trgDcmDir : string
        Directory containing the Target DICOMs.
    forceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    p2c : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    useCaseThatApplies : string
        String that denotes the Use Case that applies.
    useCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as useCaseThatApplies unless forceRegistration is True and
        useCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    from dicom_tools.metadata import get_dcm_uids
    from image_tools.attrs_info import get_im_attrs
    from general_tools.general import items_unique_to_within
    
    
    """ Get the DICOM UIDs. """
    srcStudyuid, srcSeriesuid, srcFORuid,\
    srcSOPuids = get_dcm_uids(DicomDir=srcDcmDir)
    
    trgStudyuid, trgSeriesuid, trgFORuid,\
    trgSOPuids = get_dcm_uids(DicomDir=trgDcmDir)
        
    if trgSlcNum:
        if trgSlcNum > len(trgSOPuids) - 1:
            msg = f'The input argument "trgSlcNum" = {trgSlcNum} exceeds '\
                  + f'the number of Target DICOMs ({len(trgSOPuids)}).'
            
            raise Exception(msg)
    
    
    """ Get the Image Attributes for Source and Target using SimpleITK. """
    srcSize, srcSpacing, srcST, srcIPPs, srcDirs,\
    warnings = get_im_attrs(DicomDir=srcDcmDir, Package='sitk')
    
    trgSize, trgSpacing, trgST, trgIPPs, trgDirs,\
    warnings = get_im_attrs(DicomDir=trgDcmDir, Package='sitk')
    
    
    """ Are the Source and Target part of the same study, series and do they 
    have the same SOP UIDs? """
    if srcStudyuid == trgStudyuid:
        sameStudy = True
    else:
        sameStudy = False
        
    if srcSeriesuid == trgSeriesuid:
        SameSeries = True
    else:
        SameSeries = False
        
    if srcFORuid == trgFORuid:
        sameFOR = True
    else:
        sameFOR = False
        
    if srcSOPuids == trgSOPuids:
        sameSOPs = True
    else:
        sameSOPs = False
    
    """ Are the Source and Target spacings, directions and origins equal to  
    within 1e-06? """
    sameSpacings = items_unique_to_within(srcSpacing, trgSpacing)
    sameDirs = items_unique_to_within(srcDirs, trgDirs)
    sameOrigins = items_unique_to_within(srcIPPs[0], trgIPPs[0])
    sameST = items_unique_to_within([srcST], [trgST])
    
    if srcSize == trgSize:
        sameSize = True
    else:
        sameSize = False
    
        
    """ Check which case follows. """
    if sameFOR:
        if sameDirs:
            """ 
            Comment 11/02:
                
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
            considered for future work. """
            if sameSpacings and sameOrigins:
            #if sameSpacings:
                if SameSeries:
                    useCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(trgSlcNum, int):
                        useCaseThatApplies = '2a' # Direct copy
                    else:
                        useCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ Comment 11/02:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. """
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
    
        
    """ Print comparisons to the console. """
    if p2c:
        from ImageTools import import_im, GetImageExtent
        
        srcIm = import_im(srcDcmDir)
        trgIm = import_im(trgDcmDir)
        
        """ Get the image extents. """
        srcImExtent = GetImageExtent(srcIm)
        trgImExtent = GetImageExtent(trgIm)
        
        """ Are the image extents equal to within 1e-06? """
        sameExtents = items_unique_to_within(srcImExtent, trgImExtent)
    
        #CompareImageAttributes(srcDcmDir, trgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running WhichUseCase():')
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
            
        if SameSeries:
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
