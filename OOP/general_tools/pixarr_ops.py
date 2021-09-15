# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 12:42:57 2021

@author: ctorti
"""


""" Functions that apply basic operations on pixel arrays. """

import numpy as np
from dicom_tools.metadata import (
    get_dcm_uids, get_roicol_nums, get_roicol_labels
    )
#from seg_tools.metadata import get_f2sIndsBySeg


def mean_frame_in_pixarr(
        pixarr, f2sInds, makeBinary=True, binaryThresh=0.5, p2c=False
        ):
    """
    Perform pixel-by-pixel mean of all frames in a pixel array. Output a 
    binary pixel array if binary=True (or undefined).
    
    Parameters
    ----------
    pixarr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
    f2sInds : List of ints
        A list (for each frame) of the slice numbers that correspond to each 
        frame in pixarr.
    makeBinary : bool, optional (True by default)
        If makeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by binaryThresh.
    binaryThresh : float, optional (0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if makeBinary = True.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    meanPixarr : Numpy array
        Mean pixel array.
    meanF2Sind : int or float
        The integer or fractional frame-to-slice index corresponding to the 
        mean pixel array.  meanF2Sind will be an integer if meanF2Sind % 1 = 0,
        and a fractional float otherwise.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of mean_frame_in_pixarr():')
        print('\n\n', '-'*120)
        
    F, R, C = pixarr.shape
    
    if p2c:
        print(f'\npixarr.shape = {pixarr.shape}')
        #print(f'np.amax(pixarr, axis=0) = {np.amax(pixarr, axis=0)}')
        #print(f'np.max(pixarr, axis=0) = {np.max(pixarr, axis=0)}')
        print(f'Max along each frame = {[np.amax(pixarr[f]) for f in range(pixarr.shape[0])]}')
    
    meanPixarr = np.mean(pixarr, axis=0)
    
    if p2c:
        print(f'Max of meanPixarr after averaging = {np.amax(meanPixarr)}')
            
    if makeBinary:
        #meanPixarr = (meanPixarr >= binaryThresh) * meanPixarr
        meanPixarr = (meanPixarr >= binaryThresh) * 1
        
        if p2c:
            print('Max of meanPixarr after binary thresholding =',
                  f'{np.amax(meanPixarr)}')
    
    meanPixarr = np.reshape(meanPixarr, (1, R, C))
    
    if p2c:
        print(f'Max of meanPixarr after reshape = {np.amax(meanPixarr)}')
    
    meanPixarr = meanPixarr.astype(dtype='uint8')
    
    if p2c:
        print(f'Max of meanPixarr after change to uint8 = {np.amax(meanPixarr)}')
    
    # Find the mean slice index:
    meanF2Sind = sum(f2sInds)/len(f2sInds)
    
    if meanF2Sind % 1 == 0:
        meanF2Sind = int(meanF2Sind)
    
    if p2c:
        print('-'*120)
        
    return meanPixarr, meanF2Sind

def mean_frame_in_pixarrBySeg(
        pixarrBySeg, f2sIndsBySeg, makeBinary=True, binaryThresh=0.5, p2c=False
        ):
    """
    Perform pixel-by-pixel mean of all frames in all pixel arrays in a list
    of pixel arrays. 
    
    No changes will be made to empty or single-framed pixel arrays.  If there 
    are less than 2 frames in any given pixel array the pixel array will be 
    returned unchanged. Likewise for the list of frame-to-slice-indices-by-seg.
    
    Parameters
    ----------
    pixarrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    f2sIndsBySeg : list of a list of ints
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in pixarrBySeg.
    makeBinary : bool, optional (True by default)
        If makeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by binaryThreshold.
    binaryThresh : float, optional (0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if makeBinary = True.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    meanPixarrByRoi : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
    meanF2SindByRoi : list of a list of an int or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the mean pixel array.  meanF2Sind will be 
        an integer if meanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of mean_frame_in_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    meanPixarrBySeg = []
    meanF2SindBySeg = []
    
    for s in range(len(pixarrBySeg)):
        pixarr = pixarrBySeg[s]
        f2sInds = f2sIndsBySeg[s]
        
        # The number of frames in pixarr:
        F = pixarr.shape[0]
        
        if F > 1: 
            if p2c:
                print(f'\npixarrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be averaged.')
            
            meanPixarr, meanF2Sind = mean_frame_in_pixarr(
                pixarr, f2sInds, makeBinary, binaryThresh, p2c
                )
            
            meanPixarrBySeg.append(meanPixarr)
            meanF2SindBySeg.append([meanF2Sind])
            
            if p2c:
                print(f'\npixarr had {F} frames but now has',
                      f'{meanPixarr.shape[0]}.')
        else:
            meanPixarrBySeg.append(pixarr)
            
            meanF2SindBySeg.append(f2sInds)
            
            if p2c:
                print(f'\npixarr has {F} frames so no frame averaging required.')
    
    if p2c:
        print('-'*120)
        
    return meanPixarrBySeg, meanF2SindBySeg

def or_frame_of_pixarrBySeg(pixarrBySeg, f2sIndsBySeg, p2c=False):
    """
    Perform pixel-by-pixel logical OR of all frames in all pixel arrays in a 
    list of pixel arrays. 
    
    No changes will be made to empty or single-framed pixel arrays. If there 
    are less than 2 frames in any given pixel array the pixel array will be 
    returned unchanged. Likewise for the list of frame-to-slice-indices-by-seg.
    
    Parameters
    ----------
    pixarrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    f2sIndsBySeg : list of a list of ints
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in pixarrBySeg.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    orPixarrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
    meanF2SindBySeg : list of a list of an int or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the OR pixel array.  meanF2Sind will be 
        an integer if meanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of or_frame_of_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    # Get the shape the pixel arrays (the first will do):
    F, R, C = pixarrBySeg[0].shape
    
    orPixarrBySeg = []
    meanF2SindBySeg = []
    
    for s in range(len(pixarrBySeg)):
        pixarr = pixarrBySeg[s]
        f2sInds = f2sIndsBySeg[s]
        
        F = pixarr.shape[0] # number of frames in pixarr
        
        if F > 1: 
            if p2c:
                print(f'\npixarrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be logically ORed.')
            
            orPixarr = pixarr.any(axis=0)
            
            # Reshape to a (1, R, C) pixel array:
            orPixarr = np.reshape(orPixarr, (1, R, C))
            
            # Convert from boolean to uint8:
            orPixarr = orPixarr.astype('uint8')
            
            orPixarrBySeg.append(orPixarr)
            
            # Find the mean slice index:
            meanF2Sind = sum(f2sInds)/len(f2sInds)
            
            """
            14/02: 
            It may be useful to return the fractional index but for now
            I'll simply round it to the nearest integer...
            if meanF2Sind % 1 == 0:
                meanF2Sind = int(meanF2Sind)
            """
            meanF2Sind = round(meanF2Sind)
            
            meanF2SindBySeg.append([meanF2Sind])
            
            if p2c:
                print(f'\npixarr had {F} frames but now has,',
                      f'{orPixarr.shape[0]}.')
        else:
            orPixarrBySeg.append(pixarr)
            
            meanF2SindBySeg.append(f2sInds)
            
            if p2c:
                print(f'\npixarr has {F} frames so no logical ORing required.')
    
    if p2c:
        print(f'meanF2SindBySeg = {meanF2SindBySeg}')
        print('-'*120)
        
    return orPixarrBySeg, meanF2SindBySeg

def get_frame_from_pixarr(seg, dicomDir, searchStr, slcNum, p2c):
    """
    Extract a single (2D) frame from a 3D pixel array that matches a given
    segment label.  
    
    Parameters
    ----------
    seg : Pydicom Object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    searchStr : str
        All or part of the Segment label containing the segmentation to be
        copied.
    slcNum : int
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    frame : Numpy array
        The desired frame from pixarr.
    """
    
    segLabels = get_roicol_labels(seg)
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  segNum is equivalent to SegmentNumber in SegmentSequence. """
    segNum = get_roicol_nums(seg, searchStr)
    segNum = segNum[0] # 09/07/21 Not sure about this
    
    studyuid, seriesuid, FORuid, SOPuids = get_dcm_uids(dicomDir)
    
    PFFGStoSliceIndsBySeg = get_PFFGS_to_sliceindsBySeg(seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[segNum]
    
    if not slcNum in PFFGStoSliceIndsInSeg:
        msg = 'There is no segmentation in the segment matching the term '\
              + f'"{searchStr}" for slice {slcNum} in the SEG.'
              
        raise Exception(msg)
    
    # Get the corresponding frameNum:
    
    n = 0 # initialise frame number counter
    
    for i in range(len(segLabels)):
        if i == segNum:
            """ 
            This segment contains the frame to be copied, whose position is
            given by the index of slcNum in 
            SrcPFFGStoSliceIndsBySeg[i] plus any frames that preceeded it 
            (i.e. the frame counter n).
            """
            frameNum = n + PFFGStoSliceIndsBySeg[i].index(slcNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSliceIndsBySeg[i]) 
    
    if p2c:
        print('\n\nResults of get_frame_from_pixarr:')
        print(f'   slcNum = {slcNum} relates to frameNum = {frameNum} in',
              'PixelArray')
    
    pixarr = seg.pixel_array
    
    if len(pixarr.shape) < 3:
        # This is a single-frame pixarr, i.e. shape (R, C). Reshape it to 
        # shape (1, R, C):
        R, C = pixarr.shape
        
        pixarr = np.reshape(pixarr, (1, R, C))
    else:
        # This is a multi-frame pixarr, i.e. shape (AllF, R, C).
        AllF, R, C = pixarr.shape
    
    # Initialise frame:
    frame = np.zeros((1, R, C), dtype='uint')
    
    frame[0] = pixarr[frameNum]
            
    return frame

def get_frames_from_pixarr(seg, dicomDir, searchStr, p2c):
    """
    Extract all frames from a 3D pixel array that match a given segment label.  
    
    Parameters
    ----------
    seg : Pydicom Object
        SEG object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    searchStr : str
        All or part of the Segment Label of the segment containing the 
        segmentation to be copied.
    slcNum : int
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    frames : Numpy array
        The desired frames from pixarr.
    """
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  segNum is equivalent to SegmentNumber in SegmentSequence. """
    segNum = get_roicol_nums(seg, searchStr)
    segNum = segNum[0] # 09/07/21 Not sure about this
    
    studyuid, seriesuid, foruid, SOPuids = get_dcm_uids(dicomDir)
    
    PFFGStoSliceIndsBySeg = get_PFFGS_to_sliceindsBySeg(seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[segNum]
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSliceIndsInSeg)
    
    # Get the number of rows and columns required for frames:
    pixarr = seg.pixel_array
    
    if len(pixarr.shape) < 3:
        # This is a single-frame pixarr, i.e. shape (1, R, C).
        R, C = pixarr.shape
    else:
        # This is a multi-frame pixarr, i.e. shape (AllF, R, C).
        AllF, R, C = pixarr.shape
    
    frames = np.zeros((F, R, C), dtype='uint')
    
    for f in range(F):
        slcNum = PFFGStoSliceIndsInSeg[f]
        
        frame = get_frame_from_pixarr(
            seg, dicomDir, searchStr, slcNum, p2c
            )
        
        frames[f] = frame
    
    if p2c:
        print('\n\n***Result from get_frames_from_pixarr:')
        print(f'   There are {F} frames in the segment matching',
              f'"{searchStr}".')
        print(f'   PFFGStoSliceIndsInSeg = {PFFGStoSliceIndsInSeg}')
        print(f'   pixarr.shape = {pixarr.shape}')
        if len(pixarr.shape) > 2:
            [print(f'   np.amax(pixarr[{i}]) = {np.amax(pixarr[i])}') for i in range(pixarr.shape[0])]
        else:
            print(f'   np.amax(pixarr) = {np.amax(pixarr)}')
        print(f'   frames.shape = {frames.shape}')
        if len(frames.shape) > 2:
            [print(f'   np.amax(frames[{i}]) = {np.amax(frames[i])}') for i in range(frames.shape[0])]
        else:
            print(f'   np.amax(frames) = {np.amax(frames)}')
            
    return frames