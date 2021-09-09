# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:50:51 2021

@author: ctorti
"""

from importlib import reload
import image_tools.attrs_info
reload(image_tools.attrs_info)
import seg_tools.metadata
reload(seg_tools.metadata)
import seg_tools.pixarrs
reload(seg_tools.pixarrs)
import conversion_tools.pixarrs_ims
reload(conversion_tools.pixarrs_ims)

import numpy as np
from copy import deepcopy
from general_tools.console_printing import print_inds_by_roi
from general_tools.console_printing import print_pixarr_shape_by_seg
from image_tools.attrs_info import (
    get_im_attrs, get_im_attrs_from_list_of_dicomDir
    )
from dicom_tools.metadata import get_dcm_uids
from dicom_tools.imports import get_dcm_fpaths
from seg_tools.metadata import get_p2sIndsBySeg
from seg_tools.pixarrs import get_pixarrBySeg
from conversion_tools.pixarrs_ims import imBySeg_to_pixarrBySeg


def get_seg_data_of_interest(
        seg, allF2SindsBySeg, segNums, segLabs, segLab, slcNum, p2c=False
        ):
    # TODO tidy docstrings 16/07
    """   
    Get data of interest from a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        The SEG object.
    allF2SindsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in seg.
    segNums : list of ints
        List of the segment number(s) that match the segment of interest.  
        The list may have a single int.
    segLabs : list of strs
        List of all segment labels. Only used if printing results to console.
    segLab : str or None
        The SEG label of interest, or None (if all segments are of interest).
    slcNum : int or None
        The index within the DICOM stack of the segmentation of interest, or
        None if all segmentations are of interest.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
    f2sIndsBySeg : list of a list of int
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg, i.e. the list is a sub-list 
        of allF2SindsBySeg.
        
    Notes
    -----
    There are 3 possible main use cases:
        1. Copy a single segmentation
        2. Copy an entire segment
        3. Copy all segments
        
    Which main case applies depends on the inputs:
        slcNum
        segLab
    
    If slcNum != None (i.e. if a slice number is defined) a single 
    segmentation will be copied from the segment that matches segLab 
    (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple segments and 
    segLab = None.  Possible outcomes:

        - Raise exception with message "There are multiple segments so the 
        description of the desired segment must be provided"
        - Copy the segmentation to all segments that have a segmentation on 
        slice slcNum (* preferred option?)
        
    If slcNum = None but segLab != None, all segmentations within the
    segment given by srcSegLab will be copied (--> Main Use Case 2).  
    
    If slcNum = None and segLab = None, all segmentations within all
    segments will be copied (--> Main Use Case 3).
    
    If slcNum = None and segLab = None, all segments will be copied.  
    If not, get the segment number (zero-indexed) to be copied, and reduce 
    pixarrBySeg and f2sIndsBySeg to the corresponding item.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_seg_data_of_interest():')
        print('\n\n', '-'*120)
        print(f'   segLab = {segLab}')
        print(f'   slcNum = {slcNum}')
    
    if seg == None:
        return None, None
    
    # Get all pixel arrays by seg:
    """ 16/07: making changes to get_pixarrBySeg... """
    allPixarrBySeg = get_pixarrBySeg(
        seg=seg, f2sIndsBySeg=allF2SindsBySeg, p2c=p2c
        )
    
    #print('\n\nallPixarrBySeg =', allPixarrBySeg)
    #print(f'\n\nallF2SindsBySeg = {allF2SindsBySeg}')
    
    Nsegs = len(allF2SindsBySeg)
    
    # The shape of the first pixel array:
    F, R, C = allPixarrBySeg[0].shape
    
    if p2c:
        print_inds_by_roi(allF2SindsBySeg)
    
    # Initialise variable that indicates if the SEG data was reduced by
    # frame/slice number (FromSliceNum):
    reducedBySlc = False
            
    if slcNum:
        # Limit data to those which belong to the chosen slice number:
        pixarrBySeg = []
        f2sIndsBySeg = []
    
        for s in range(Nsegs):
            # The pixel array and frame-to-slice indices for this segment:
            allPixarr = deepcopy(allPixarrBySeg[s])
            allF2Sinds = deepcopy(allF2SindsBySeg[s])
                
            # Get the frame number(s) that relate to slcNum:
            frameNums = [i for i, x in enumerate(allF2Sinds) if x==slcNum]
            
            if p2c:
                print(f'\nallF2Sinds = {allF2Sinds}')
                print(f'frameNums = {frameNums}')
            
            if len(frameNums) != len(allF2Sinds):
                # Some frames will be rejected:
                reducedBySlc = True
                
                # Initialise pixarr with the reduced number of frames matching
                # len(frameNums):
                pixarr = np.zeros((len(frameNums), R, C), dtype='uint')
                f2sInds = []
                
                if frameNums:
                    # Keep only the indeces and frames that relate to frameNums:
                    for f in range(len(frameNums)):
                        pixarr[f] = allPixarr[frameNums[f]]
                        
                        f2sInds.append(allF2Sinds[frameNums[f]])
                    
                    if p2c:
                        print(f'f2sInds = {f2sInds}')
                    
                    # Append non-empty pixel arrays and non-empty f2sInds:
                    pixarrBySeg.append(pixarr)
                    f2sIndsBySeg.append(f2sInds)
            """ 
            else 
            all frames to remain
            """
        
        #print('\n\npixarrBySeg =', pixarrBySeg)
        
        if reducedBySlc:
            # Replace allPixarrBySeg and allF2SindsBySeg with pixarrBySeg and 
            # f2sIndsBySeg in case further restricting of data is required 
            # below:
            allPixarrBySeg = deepcopy(pixarrBySeg)
            allF2SindsBySeg = deepcopy(f2sIndsBySeg)
        """ 
        else 
            allPixarrBySeg and allF2SindsBySeg remain unchanged
        """
        
        #print('\n\nallPixarrBySeg =', allPixarrBySeg)
        
        if p2c and reducedBySlc:
            print('\n   After limiting data to those that relate to slice',
                  f'number {slcNum}:')#, the f2sIndsBySeg =')
            print_inds_by_roi(f2sIndsBySeg)
            print('   pixarrBySeg = ', pixarrBySeg)
            print_pixarr_shape_by_seg(pixarrBySeg)
    
    # Initialise variable that indicates if the SEG data was reduced by
    # chosen segment label (segLab):
    #reducedBySeg = False
    
    if segLab:
        # Limit data to those which belong to the chosen segment(s):
        pixarrBySeg = []
        f2sIndsBySeg = []
        
        #print(f'\nsegNums = {segNums}')
        #print(f'allF2SindsBySeg = {allF2SindsBySeg}')
        
        if len(segNums) != len(allF2SindsBySeg):
            #reducedBySeg = True
            
            pixarrBySeg = [allPixarrBySeg[s] for s in segNums]
            f2sIndsBySeg = [allF2SindsBySeg[s] for s in segNums]
            
            if p2c:
                # Limit the list of segment names to those that belong to the
                # chosen segment(s):
                segLabs = [segLabs[i] for i in segNums]
            
                print('\n   After limiting data to those whose segment name',
                      f'matches {segLabs}:')
                print_inds_by_roi(f2sIndsBySeg)
                print_pixarr_shape_by_seg(pixarrBySeg)
                
        else:
            pixarrBySeg = deepcopy(allPixarrBySeg)
            f2sIndsBySeg = deepcopy(allF2SindsBySeg)
        
    else:
        pixarrBySeg = deepcopy(allPixarrBySeg)
        f2sIndsBySeg = deepcopy(allF2SindsBySeg)
    
    if p2c:
        print('\n   Final outputs of get_seg_data_of_interest():')
        print_pixarr_shape_by_seg(pixarrBySeg)
        print_inds_by_roi(f2sIndsBySeg)
        print('-'*120)
        
    return pixarrBySeg, f2sIndsBySeg

def raise_error_if_no_seg_data_of_interest(f2sIndsBySeg, segLabs, slcNum):
    """
    Raise error if there is no matching SEG data of interest.

    Parameters
    ----------
    f2sIndsBySeg : list of a list of int
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg, i.e. the list is a sub-list 
        of allF2SindsBySeg.
    segLabs : list of strs
        List (for each segment) of segment label(s).
    slcNum : int or None

    Returns
    -------
    None.

    """
    
    if f2sIndsBySeg == []:
        msg = "There are no segmentations "
        
        if slcNum:
            msg += f"on slice {slcNum} "
              
        msg += f"in any segment. There are {len(segLabs)} segments:"
        
        for i in range(len(segLabs)):
            msg += f"\nSegment {i+1} with label '{segLabs[i]}' has "\
                   + f"{len(f2sIndsBySeg[i])} segmentations that correspond "\
                   + f"to slices {f2sIndsBySeg[i]}" 
        raise Exception(msg)

def get_seg_data_from_list_of_segs(
        listOfSegs, listOfDicomDirs, p2c=False):
    """ 
    Note:
    
    listOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    listOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_seg_data_from_list_of_segs():')
        print('\n\n', '-'*120)
    
    listOfPixarrBySeg = []
    listOfF2SindsBySeg = []
    #listOfSegNums = []
    listOfImSizes = []
    listOfImSpacings = []
    listOfImSlcThicks = []
    listOfImIPPs = []
    listOfImDirections = []
    listOfDicomFpaths = []
    
    for i in range(len(listOfSegs)):
        if p2c:
            print(f'\n\nlistOfSegs[{i}]:')
            
        studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(listOfDicomDirs[i])
        
        if listOfSegs[i]:
            f2sIndsBySeg = get_p2sIndsBySeg(listOfSegs[i], SOPUIDs)
            
            pixarrBySeg = get_pixarrBySeg(listOfSegs[i], f2sIndsBySeg, p2c)
            
            if p2c:      
                print(f'f2sIndsBySeg = {f2sIndsBySeg}') 
        else:
            pixarrBySeg = None
            f2sIndsBySeg = None
            
            if p2c:      
                print('No segmentations for this dataset')
        
        listOfPixarrBySeg.append(pixarrBySeg)
        listOfF2SindsBySeg.append(f2sIndsBySeg)
            
        if p2c:      
            print(f'\nlen(listOfF2SindsBySeg) = {len(listOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
        imSize, imSpacings, imSlcThick, imIPPs, imDirections,\
            warnings = get_im_attrs(listOfDicomDirs[i], p2c=p2c)
        
        listOfImSizes.append(imSize)
        listOfImSpacings.append(imSpacings)
        listOfImSlcThicks.append(imSlcThick)
        listOfImIPPs.append(imIPPs)
        listOfImDirections.append(imDirections)
        listOfDicomFpaths.append(get_dcm_fpaths(listOfDicomDirs[i]))
        
    #listOfDicomFpaths, listOfIPPs, listOfDirections, listOfSpacings\
    #    = GetImAttributesFromListOfDcmDirs(listOfDicomDirs)
    
    if p2c:
        print('-'*120)
        
    return listOfPixarrBySeg, listOfF2SindsBySeg, listOfImSizes,\
        listOfImSpacings, listOfImIPPs, listOfImDirections, listOfDicomFpaths

def get_seg_data_from_list_of_labimByRoi(
        listOfLabimByRoi, listOfDicomDirs, p2c=False):
    """ 
    listOfLabimByRoi is a list (for each dataset, which could be of length 1) 
    of a list (for each segment) of label images (SimpleITK Images).
    
    listOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    """
    
    if p2c:
        print('\n\n' + '-'*120)
        print('---> get_seg_data_from_list_of_labimByRoi():\n')
        #print('-'*120)
    
    listOfPixarrByRoi = []
    listOfF2SindsByRoi = []
    
    for i in range(len(listOfLabimByRoi)):
        labimByRoi = listOfLabimByRoi[i]
        
        if p2c:
            print(f'listOfLabimByRoi[{i}]:\n')
        
        if labimByRoi:
            pixarrByRoi, f2sIndsByRoi = imBySeg_to_pixarrBySeg(
                labimByRoi, p2c
                )
            
            if p2c:      
                print(f'  f2sIndsByRoi = {f2sIndsByRoi}\n')  
        else:
            #pixarrByRoi = None
            #f2sIndsByRoi = None
            pixarrByRoi = []
            f2sIndsByRoi = []
            
            if p2c:      
                print('  No segmentations for this dataset.\n')
        
        listOfPixarrByRoi.append(pixarrByRoi)
        listOfF2SindsByRoi.append(f2sIndsByRoi)
            
        if p2c:      
            print(f'  len(listOfF2SindsByRoi) = {len(listOfF2SindsByRoi)}\n')
            print(f'  listOfF2SindsByRoi = {listOfF2SindsByRoi}\n')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    listOfDicomFpaths, listOfSizes, listOfSpacings, listOfSlcThick, listOfIPPs,\
        listOfDirections = get_im_attrs_from_list_of_dicomDir(listOfDicomDirs)
    
    if p2c:
        print('<--- get_seg_data_from_list_of_labimByRoi()')
        print('-'*120 + '\n')
        
    return listOfPixarrByRoi, listOfF2SindsByRoi, listOfDicomFpaths,\
           listOfIPPs, listOfDirections, listOfSpacings