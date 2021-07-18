# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 12:42:57 2021

@author: ctorti
"""


""" Functions that apply basic operations on pixel arrays. """

import numpy as np
from dicom_tools.metadata import get_dcm_uids, get_roicol_nums, get_roicol_labels
from seg_tools.metadata import get_PFFGS_to_sliceindsBySeg


def mean_frame_in_pixarr(PixArr, F2Sinds, MakeBinary=True, BinaryThresh=0.5,
                         LogToConsole=False):
    """
    Perform pixel-by-pixel mean of all frames in a pixel array. Output a 
    binary pixel array if binary=True (or undefined).
    
    Parameters
    ----------
    PixArr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
    F2Sinds : List of ints
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.
    MakeBinary : bool, optional (True by default)
        If MakeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by BinaryThreshold.
    BinaryThreshold : float, optional (0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if MakeBinary = True.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    MeanPixArr : Numpy array
        Mean pixel array.
    MeanF2Sind : int or float
        The integer or fractional frame-to-slice index corresponding to the 
        mean pixel array.  MeanF2Sind will be an integer if MeanF2Sind % 1 = 0,
        and a fractional float otherwise.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of mean_frame_in_pixarr():')
        print('\n\n', '-'*120)
        
    F, R, C = PixArr.shape
    
    if LogToConsole:
        print(f'\nPixArr.shape = {PixArr.shape}')
        #print(f'np.amax(PixArr, axis=0) = {np.amax(PixArr, axis=0)}')
        #print(f'np.max(PixArr, axis=0) = {np.max(PixArr, axis=0)}')
        print(f'Max along each frame = {[np.amax(PixArr[f]) for f in range(PixArr.shape[0])]}')
    
    MeanPixArr = np.mean(PixArr, axis=0)
    
    if LogToConsole:
        print(f'Max of MeanPixArr after averaging = {np.amax(MeanPixArr)}')
            
    if MakeBinary:
        #MeanPixArr = (MeanPixArr >= BinaryThresh) * MeanPixArr
        MeanPixArr = (MeanPixArr >= BinaryThresh) * 1
        
        if LogToConsole:
            print('Max of MeanPixArr after binary thresholding =',
                  f'{np.amax(MeanPixArr)}')
    
    MeanPixArr = np.reshape(MeanPixArr, (1, R, C))
    
    if LogToConsole:
        print(f'Max of MeanPixArr after reshape = {np.amax(MeanPixArr)}')
    
    MeanPixArr = MeanPixArr.astype(dtype='uint8')
    
    if LogToConsole:
        print(f'Max of MeanPixArr after change to uint8 = {np.amax(MeanPixArr)}')
    
    """ Find the mean slice index: """
    MeanF2Sind = sum(F2Sinds)/len(F2Sinds)
    
    if MeanF2Sind % 1 == 0:
        MeanF2Sind = int(MeanF2Sind)
    
    if LogToConsole:
        print('-'*120)
        
    return MeanPixArr, MeanF2Sind

def mean_frame_in_pixarrBySeg(PixArrBySeg, F2SindsBySeg, MakeBinary=True, 
                              BinaryThresh=0.5, LogToConsole=False):
    """
    Perform pixel-by-pixel mean of all frames in all pixel arrays in a list
    of pixel arrays. 
    
    No changes will be made to empty or single-framed pixel arrays.  If there 
    are less than 2 frames in any given pixel array the pixel array will be 
    returned unchanged. Likewise for the list of frame-to-slice-indices-by-seg.
    
    Parameters
    ----------
    PixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    F2SindsBySeg : list of a list of ints
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in PixArrBySeg.
    MakeBinary : bool, optional (True by default)
        If MakeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by BinaryThreshold.
    BinaryThreshold : float, optional (0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if MakeBinary = True.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    MeanPixArrByRoi : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
    MeanF2SindByRoi : list of a list of an int or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the mean pixel array.  MeanF2Sind will be 
        an integer if MeanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of mean_frame_in_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    MeanPixArrBySeg = []
    MeanF2SindBySeg = []
    
    for s in range(len(PixArrBySeg)):
        PixArr = PixArrBySeg[s]
        F2Sinds = F2SindsBySeg[s]
        
        """ The number of frames in PixArr """
        F = PixArr.shape[0]
        
        if F > 1: 
            if LogToConsole:
                print(f'\nPixArrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be averaged.')
            
            MeanPixArr, MeanF2Sind\
                = mean_frame_in_pixarr(PixArr, F2Sinds, MakeBinary, 
                                       BinaryThresh, LogToConsole)
            
            MeanPixArrBySeg.append(MeanPixArr)
            MeanF2SindBySeg.append([MeanF2Sind])
            
            if LogToConsole:
                print(f'\nPixArr had {F} frames but now has',
                      f'{MeanPixArr.shape[0]}.')
        else:
            MeanPixArrBySeg.append(PixArr)
            
            MeanF2SindBySeg.append(F2Sinds)
            
            if LogToConsole:
                print(f'\nPixArr has {F} frames so no frame averaging required.')
    
    if LogToConsole:
        print('-'*120)
        
    return MeanPixArrBySeg, MeanF2SindBySeg

def or_frame_of_pixarrBySeg(PixArrBySeg, F2SindsBySeg, LogToConsole=False):
    """
    Perform pixel-by-pixel logical OR of all frames in all pixel arrays in a 
    list of pixel arrays. 
    
    No changes will be made to empty or single-framed pixel arrays. If there 
    are less than 2 frames in any given pixel array the pixel array will be 
    returned unchanged. Likewise for the list of frame-to-slice-indices-by-seg.
    
    Parameters
    ----------
    PixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    F2SindsBySeg : list of a list of ints
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in PixArrBySeg.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    OrPixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
    MeanF2SindBySeg : list of a list of an int or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the OR pixel array.  MeanF2Sind will be 
        an integer if MeanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of or_frame_of_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    # Get the shape of any pixel array:
    F, R, C = PixArrBySeg[0].shape
    
    OrPixArrBySeg = []
    MeanF2SindBySeg = []
    
    for s in range(len(PixArrBySeg)):
        PixArr = PixArrBySeg[s]
        F2Sinds = F2SindsBySeg[s]
        
        """ The number of frames in PixArr """
        F = PixArr.shape[0]
        
        if F > 1: 
            if LogToConsole:
                print(f'\nPixArrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be logically ORed.')
            
            OrPixArr = PixArr.any(axis=0)
            
            """ Reshape to a (1, R, C) pixel array: """
            OrPixArr = np.reshape(OrPixArr, (1, R, C))
            
            """ Convert from boolean to uint8: """
            OrPixArr = OrPixArr.astype('uint8')
            
            OrPixArrBySeg.append(OrPixArr)
            
            """ Find the mean slice index: """
            MeanF2Sind = sum(F2Sinds)/len(F2Sinds)
            
            """
            14/02: It may be useful to return the fractional index but for now
            I'll simply round it to the nearest integer...
            if MeanF2Sind % 1 == 0:
                MeanF2Sind = int(MeanF2Sind)
            """
            MeanF2Sind = round(MeanF2Sind)
            
            MeanF2SindBySeg.append([MeanF2Sind])
            
            if LogToConsole:
                print(f'\nPixArr had {F} frames but now has,',
                      f'{OrPixArr.shape[0]}.')
        else:
            OrPixArrBySeg.append(PixArr)
            
            MeanF2SindBySeg.append(F2Sinds)
            
            if LogToConsole:
                print(f'\nPixArr has {F} frames so no logical ORing required.')
    
    if LogToConsole:
        print(f'MeanF2SindBySeg = {MeanF2SindBySeg}')
        print('-'*120)
        
    return OrPixArrBySeg, MeanF2SindBySeg

def get_frame_from_pixarr(Seg, DicomDir, SearchString, SliceNum, LogToConsole):
    """
    Extract a single (2D) frame from a 3D pixel array that matches a given
    segment label.  
    
    Parameters
    ----------
    Seg : Pydicom Object
        SEG object.
    DicomDir : str
        Directory containing the corresponding DICOMs.
    SearchString : str
        All or part of the Segment label containing the segmentation to be
        copied.
    SliceNum : int
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    Frame : Numpy array
        Frame from PixArr.
    """
    
    SegLabels = get_roicol_labels(Seg)
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = get_roicol_nums(Seg, SearchString)
    SegNum = SegNum[0] # 09/07/21 Not sure about this
    
    Studyuid, Seriesuid, FORuid, SOPuids = get_dcm_uids(DicomDir)
    
    PFFGStoSliceIndsBySeg = get_PFFGS_to_sliceindsBySeg(Seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[SegNum]
    
    if not SliceNum in PFFGStoSliceIndsInSeg:
        msg = 'There is no segmentation in the segment matching the term '\
              + f'"{SearchString}" for slice {SliceNum} in the SEG.'
              
        raise Exception(msg)
    
    # Get the corresponding FrameNum:
    
    n = 0 # initialise frame number counter
    
    for i in range(len(SegLabels)):
        if i == SegNum:
            """ 
            This segment contains the frame to be copied, whose position is
            given by the index of FromSliceNum in 
            SrcPFFGStoSliceIndsBySeg[i] plus any frames that preceeded it 
            (i.e. the frame counter n).
            """
            FrameNum = n + PFFGStoSliceIndsBySeg[i].index(SliceNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSliceIndsBySeg[i]) 
    
    
    if LogToConsole:
        print('\n\nResults of GetFrameFromPixArr:')
        print(f'   SliceNum = {SliceNum} relates to FrameNum = {FrameNum} in',
              'PixelArray')
    
    PixArr = Seg.pixel_array
    
    if len(PixArr.shape) < 3:
        # This is a single-frame PixArr, i.e. shape (R, C). Reshape it to 
        # shape (1, R, C):
        R, C = PixArr.shape
        
        PixArr = np.reshape(PixArr, (1, R, C))
    else:
        # This is a multi-frame PixArr, i.e. shape (AllF, R, C).
        AllF, R, C = PixArr.shape
    
    # Initialise Frame:
    Frame = np.zeros((1, R, C), dtype='uint')
    
    Frame[0] = PixArr[FrameNum]
            
    return Frame

def get_frames_from_pixarr(Seg, DicomDir, SearchString, LogToConsole):
    """
    Extract all frames from a 3D pixel array that match a given segment label.  
    
    Parameters
    ----------
    Seg : Pydicom Object
        SEG object.
    DicomDir : str
        Directory containing the corresponding DICOMs.
    SearchString : str
        All or part of the Segment Label of the segment containing the 
        segmentation to be copied.
    SliceNum : int
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    Frame : Numpy array
        Frame from PixArr.
    """
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = get_roicol_nums(Seg, SearchString)
    SegNum = SegNum[0] # 09/07/21 Not sure about this
    
    Studyuid, Seriesuid, FORuid, SOPuids = get_dcm_uids(DicomDir)
    
    PFFGStoSliceIndsBySeg = get_PFFGS_to_sliceindsBySeg(Seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[SegNum]
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSliceIndsInSeg)
    
    # Get the number of rows and columns required for Frames:
    PixArr = Seg.pixel_array
    
    if len(PixArr.shape) < 3:
        """ This is a single-frame PixArr, i.e. shape (1, R, C). """
        R, C = PixArr.shape
    else:
        """ This is a multi-frame PixArr, i.e. shape (AllF, R, C). """
        AllF, R, C = PixArr.shape
    
    Frames = np.zeros((F, R, C), dtype='uint')
    
    for f in range(F):
        SliceNum = PFFGStoSliceIndsInSeg[f]
        
        Frame = get_frame_from_pixarr(Seg, DicomDir, SearchString, SliceNum, 
                                      LogToConsole)
        
        Frames[f] = Frame
    
    if LogToConsole:
        print('\n\n***Result from GetFramesFromPixArr:')
        print(f'   There are {F} frames in the segment matching',
              f'"{SearchString}".')
        print(f'   PFFGStoSliceIndsInSeg = {PFFGStoSliceIndsInSeg}')
        print(f'   PixArr.shape = {PixArr.shape}')
        if len(PixArr.shape) > 2:
            [print(f'   np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        else:
            print(f'   np.amax(PixArr) = {np.amax(PixArr)}')
        print(f'   Frames.shape = {Frames.shape}')
        if len(Frames.shape) > 2:
            [print(f'   np.amax(Frames[{i}]) = {np.amax(Frames[i])}') for i in range(Frames.shape[0])]
        else:
            print(f'   np.amax(Frames) = {np.amax(Frames)}')
            
    return Frames