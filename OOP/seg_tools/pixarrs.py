# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:19:16 2021

@author: ctorti
"""


import numpy as np
from dicom_tools.metadata import get_dcm_uids
from seg_tools.metadata import get_frameNums, get_f2sIndsBySeg
from general_tools.console_printing import (
    print_indsByRoi, print_shape_of_pixarrBySeg
    )

def get_pixarr_in_seg(seg, dicomDir, searchStr):
    """
    Get the frames in a SEG's pixel array that belong to a segment with
    specified label.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    searchStr : str
        All or part of the Segment Label containing the segment of interest.
                       
    Returns
    -------
    pixarrInSeg : Numpy array
        Sub-array of the SEG's pixel array that contains the frames that belong
        to the segment of interest. 
    f2sIndsInSeg : list of ints
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in pixarr.
    """
    
    """ 16/07: Making changes to metadata.py in seg_tools... """
    
    pixarr_all = seg.pixel_array
    
    F_all, R, C = pixarr_all.shape
    
    # Get the frame numbers of the SEG's pixel array that correspond to the 
    # segment of interest, and the corresponding Per-frame Functional Groups
    # Sequence-to-slice indices:
    frameNumsInSeg, f2sIndsInSeg = get_frameNums(seg, searchStr, dicomDir)
    
    F = len(frameNumsInSeg)
    
    pixarrInSeg = np.zeros((F, R, C), dtype='uint')
    
    for i in range(F):  
        pixarrInSeg[i] = pixarr_all[frameNumsInSeg[i]]
    
    return pixarrInSeg, f2sIndsInSeg

def get_pixarrBySeg(seg, f2sIndsBySeg, p2c=False):
    # seg, dicomDir, p2c=False
    # 16/07 need to tidy up docstrings
    """
    Get a list of pixel arrays grouped by segment in a SEG.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
    f2sIndsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
    """
    
    #p2c = params.genParams['p2c']
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    pixarr_all = seg.pixel_array
    
    """ If the seg has a single frame it will have shape (NumOfRows, NumOfCols).
    For compatibility with multi-framed pixel arrays, convert single-framed
    pixel arrays to shape (1, NumOfRows, NumOfCols). """
    
    if len(pixarr_all.shape) < 3:
        R, C = pixarr_all.shape
        
        pixarr_all = np.reshape(pixarr_all, (1, R, C))
    
    F_all, R, C = pixarr_all.shape
    
    Nsegs = len(f2sIndsBySeg)
    
    if p2c:
        print(f'pixarr_all.shape = {pixarr_all.shape}')
        #print('f2sIndsBySeg =')
        print_indsByRoi(f2sIndsBySeg)
        
    pixarrBySeg = []
    
    j = 0 # total frame counter for all frames in pixarr_all
    
    for s in range(Nsegs):
        # The number of frames in pixarr_all that belong to this segment:
        F = len(f2sIndsBySeg[s])
        
        # Initialise the pixel array for this segment:
        pixarr = np.zeros((F, R, C), dtype='uint')
        
        if F > 0:
            for i in range(F):
                pixarr[i] = pixarr_all[j]
                
                j += 1 # increment the total frame counter
        
        pixarrBySeg.append(pixarr)
    
    if p2c:
        print_shape_of_pixarrBySeg(pixarrBySeg)
        print('-'*120)
        
    return pixarrBySeg


""" Old below """





def get_pixarr_in_seg_OLD(seg, dicomDir, searchStr):
    """
    Get the frames in a SEG's pixel array that belong to a segment with
    specified label.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    searchStr : str
        All or part of the Segment Label containing the segment of interest.
                       
    Returns
    -------
    pixarrInSeg : Numpy array
        Sub-array of the SEG's pixel array that contains the frames that belong
        to the segment of interest. 
    f2sIndsInSeg : list of ints
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in pixarr.
    """
    
    pixarr_all = seg.pixel_array
    
    F_all, R, C = pixarr_all.shape
    
    # Get the frame numbers of the SEG's pixel array that correspond to the 
    # segment of interest, and the corresponding Per-frame Functional Groups
    # Sequence-to-slice indices:
    frameNumsInSeg, f2sIndsInSeg = get_frameNums(seg, searchStr, 
                                                           dicomDir)
    
    F = len(frameNumsInSeg)
    
    pixarrInSeg = np.zeros((F, R, C), dtype='uint')
    
    for i in range(F):  
        pixarrInSeg[i] = pixarr_all[frameNumsInSeg[i]]
    
    return pixarrInSeg, f2sIndsInSeg

def get_pixarrBySeg_OLD(seg, dicomDir, p2c=False):
    """
    Get a list of pixel arrays grouped by segment in a SEG.
    
    Parameters
    ----------
    seg : Pydicom object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
    f2sIndsBySeg : list of a list of ints
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in pixarrBySeg.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_pixarr_by_seg():')
        print('\n\n', '-'*120)
        
    pixarr_all = seg.pixel_array
    
    """ If the seg has a single frame it will have shape (NumOfRows, NumOfCols).
    For compatibility with multi-framed pixel arrays, convert single-framed
    pixel arrays to shape (1, NumOfRows, NumOfCols). """
    
    if len(pixarr_all.shape) < 3:
        R, C = pixarr_all.shape
        
        pixarr_all = np.reshape(pixarr_all, (1, R, C))
    
    F_all, R, C = pixarr_all.shape
    
    studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(dicomDir)
    
    f2sIndsBySeg = get_f2sIndsBySeg(seg, SOPUIDs)
    
    Nsegs = len(f2sIndsBySeg)
    
    if p2c:
        print(f'   pixarr_all.shape = {pixarr_all.shape}')
        #print('   f2sIndsBySeg =')
        print_indsByRoi(f2sIndsBySeg)
        
    pixarrBySeg = []
    
    j = 0 # total frame counter for all frames in pixarr_all
    
    for s in range(Nsegs):
        # The number of frames in pixarr_all that belong to this segment:
        F = len(f2sIndsBySeg[s])
        
        # Initialise the pixel array for this segment:
        pixarr = np.zeros((F, R, C), dtype='uint')
        
        if F > 0:
            for i in range(F):
                pixarr[i] = pixarr_all[j]
                
                j += 1 # increment the total frame counter
        
        pixarrBySeg.append(pixarr)
    
    if p2c:
        print_shape_of_pixarrBySeg(pixarrBySeg)
        print('-'*120)
        
    return pixarrBySeg, f2sIndsBySeg
