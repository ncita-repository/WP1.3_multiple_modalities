# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:54:12 2020

@author: ctorti
"""




# Import packages and functions:
#import os
#import time
#from copy import deepcopy
#import pydicom
#from pydicom.pixel_data_handlers.numpy_handler import pack_bits
#from pydicom.uid import generate_uid
#import numpy as np
#import SimpleITK as sitk




"""
******************************************************************************
******************************************************************************
SEGMENT FUNCTIONS
******************************************************************************
******************************************************************************
"""      


def ProportionOfPixArrInExtent(PixArr, FrameToSliceInds, SrcDicomDir, 
                               TrgDicomDir, LogToConsole=False):
    """
    COMMENT 29/01/2021:
        I was using the same image (i.e. same DICOM directory) to generate the
        pixel array (i.e. PixArr2PtsByContour()) and to generate the vertices
        (i.e. GetImageVertices()).  The former had to be the source DICOMs and
        the latter the target DICOMs.
        
    
    Determine what proportion of voxels that make up a 3D pixel array intersect
    with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    PixArr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    SrcDicomDir : string
        Directory containing the DICOMs that relate to PixArr.
    
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArr2PtsByContour
    from ImageTools import ImportImage, GetImageVertices
    from GeneralTools import IsPointInPolygon
    
    F, R, C = PixArr.shape
    
    # Convert the pixel array to a list (for each frame in PixArr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    PtsByObjByFrame,\
    CntDataByObjByFrame = PixArr2PtsByContour(PixArr, 
                                              FrameToSliceInds, 
                                              SrcDicomDir)
    
    Image = ImportImage(TrgDicomDir)
    
    Vertices = GetImageVertices(Image)
    
    IsInside = []
    
    for f in range(len(PtsByObjByFrame)):
        for o in range(len(PtsByObjByFrame[f])):
            for pt in PtsByObjByFrame[f][o]:
                IsInside.append(IsPointInPolygon(pt, Vertices))
    
    FracProp = sum(IsInside)/len(IsInside)
    
    if LogToConsole:
        print('\n\n\nFunction ProportionOfPixArrInExtent():')
        print(f'   PixArr has {F}x{R}x{C} (FxRxC)')
        print(f'   FrameToSliceInds = {FrameToSliceInds}')
        print(f'   len(PtsByObjByFrame) = {len(PtsByObjByFrame)}')
        print(f'   len(CntDataByObjByFrame) = {len(CntDataByObjByFrame)}')
        for f in range(len(PtsByObjByFrame)):
            print(f'   Frame {f} has {len(PtsByObjByFrame[f])} objects')
            for o in range(len(PtsByObjByFrame[f])):
                print(f'      Object {o} has {len(PtsByObjByFrame[f][o])} points')
                ##print(f'   len(CntDataByFrame[{i}]) = {len(CntDataByFrame[i])}')
        #print(f'   \nPtsByFrame = {PtsByFrame}')
        #print(f'   \nCntDataByFrame = {CntDataByFrame}')
        
        if False:
            print(f'\n   The vertices of the image grid:')
            for i in range(len(Vertices)):
                print(f'   {Vertices[i]}')
            
            print(f'\n   The points from the input indices:')
            for i in range(len(PtsByObjByFrame[f][o])):
                print(f'   {PtsByObjByFrame[f][o][i]}')
                
            print(f'\n   len(IsInside) = {len(IsInside)}')
            print(f'   sum(IsInside) = {sum(IsInside)}')
            
        if FracProp < 1:
            print(f'\n   The vertices of the image grid:')
            for i in range(len(Vertices)):
                print(f'   {Vertices[i]}')
            N = len(IsInside)
            limit = 10
            if N > limit:
                print(f'\n   The first {limit} points that lie outside of extent:')
                #print([f'\n   {PtsByObjByFrame[0][0][i]}' for i in range(5)])
                for i in range(limit):
                    print(f'   {PtsByObjByFrame[0][0][i]}')
            else:
                print(f'\n   The first {N} points that lie outside of extent:')
                #print([f'\n   {PtsByObjByFrame[0][0][i]}' for i in range(N)])
                for i in range(N):
                    print(f'   {PtsByObjByFrame[0][0][i]}')
    
    return FracProp






def ProportionOfSegsInExtent(PixArrBySeg, F2SindsBySeg, SrcDicomDir, 
                             TrgDicomDir, LogToConsole=False):
    """
    COMMENT 29/01/2021:
        I was using the same image (i.e. same DICOM directory) to generate the
        pixel array (i.e. PixArr2PtsByContour()) and to generate the vertices
        (i.e. GetImageVertices()).  The former had to be the source DICOMs and
        the latter the target DICOMs.
        
    
    Determine what proportion of voxels that make up a list of 3D pixel arrays 
    intersect with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    PixArrBySeg : list of Numpy arrays
        A list (for each segment) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    
    F2SindsBySeg : List of list of integers
        A list (for each segment) of a list (for each frame) of the slice 
        numbers that correspond to each frame in each pixel array in 
        PixArrBySeg.   
        
    SrcDicomDir : string
        Directory containing the DICOMs that relate to PixArrBySeg.
    
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArrByRoi2PtsByCntByRoi
    from RtsTools import ProportionOfRoisInExtent
    
    # Convert the list of pixel arrays-by-segment to a list of 
    # points-by-contour-by-ROI:
    PtsByCntByRoi, CntDataByCntByRoi,\
    C2SindsByRoi = PixArrByRoi2PtsByCntByRoi(PixArrBySeg, F2SindsBySeg, 
                                             SrcDicomDir, Thresh=0.5)
    
    FracProp = ProportionOfRoisInExtent(PtsByCntByRoi, TrgDicomDir, 
                                        LogToConsole)
    
    return FracProp








def GetFrameNums(Seg, SearchString, DicomDir):
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    SearchString : string
        All or part of the Segment Label containing the segmentation of
        interest.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
                            
    Outputs:
    *******
    
    FrameNums : list of integer
        List of the frame numbers that correspond to all frames in the SEG's 
        pixel array that belongs to the segment of interest.  If only one frame 
        corresponds to the segment, the returned list will be of length one 
        (e.g. [3]).
        
    PFFGStoSliceInds : list of integers
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    # Get the DimensionIndexValues (DIVs) and group the DIVs by segment.
    """ Note:  The DIVs are integers that start from 1. """
    DIVs = GetDIVs(Seg)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir)
    
    # Get the PerFrameFunctionalGroupsSequence-to-slice indices
    # (PFFGStoSliceInds) and PFFGStoSliceInds grouped by segment
    # (PFFGStoSliceIndsBySeg):
    """ Note:  The PFFGStoSliceInds are integers that start from 0. """
    PFFGStoSliceInds = GetPFFGStoSliceInds(Seg, SOPuids)
    
    PFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=PFFGStoSliceInds, 
                                               DIVs=DIVs)
    
    SegNum = GetRoiNum(Seg, SearchString)
    
    # The PFFGStoSliceInds for the segment of interest:
    PFFGStoSliceInds = PFFGStoSliceIndsBySeg[SegNum]
    
    # The number of segments:
    S = len(PFFGStoSliceIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSliceInds)
    
    n = 0 # initialise frame number counter
    
    FrameNums = []
    
    for s in range(S):
        if s == SegNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of PFFGStoSliceInds[f] in PFFGStoSliceIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                FrameNum = n + PFFGStoSliceIndsBySeg[s].index(PFFGStoSliceInds[f])
                
                FrameNums.append(FrameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSliceIndsBySeg[s])
    
    return FrameNums, PFFGStoSliceInds








def GetRSOPuidsInRIS(Seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence
    of a SEG.
    
    Inputs:
    ******
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    *******
        
    RSOPuids : list of strings
        List of ReferencedSOPInstanceUIDs.
    """
    
    RSOPuids = [] 
    
    sequences = Seg.ReferencedSeriesSequence[0]\
                   .ReferencedInstanceSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids





def GetRIStoSliceInds(Seg, SOPuids):
    """
    Get the slice numbers that correspond to each ReferencedInstanceSequence
    in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        List of SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    RIStoSliceInds : list of integers
        Slice numbers that correspond to each ReferencedInstanceSequence in 
        Seg.
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPuids = GetRSOPuidsInRIS(Seg)
    
    RIStoSliceInds = [] 
    
    for Ruid in RSOPuids:
        # Find the matching index of Ruid in SOPuids:
        RIStoSliceInds.append(SOPuids.index(Ruid))
        
    return RIStoSliceInds





def GetRSOPuidsInPFFGS(Seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the 
    Per-FrameFunctionalGroupsSequence of a SEG.
    
    Inputs:
    ******
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    *******
    
    RSOPuids : list of strings
        List of ReferencedSOPUIDs.
    """
    
    RSOPuids = [] 
    
    sequences = Seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids





def GetPFFGStoSliceInds(Seg, SOPuids):
    """
    Get a list of the slice numbers that correspond to each  
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Inputs:
    ******
        
    Seg : Pydcom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
        
    PFFGStoSliceInds : list of integers
        List (for each frame) of the slice numbers that correspond to each 
        Per-FrameFunctionalGroupsSequence in Seg.
    """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    RSOPuids = GetRSOPuidsInPFFGS(Seg)
    
    PFFGStoSliceInds = [] 
    
    for i in range(len(RSOPuids)):
        if RSOPuids[i] in SOPuids:
            # Find the matching index of RSOPuids[i] in SOPuids:
            PFFGStoSliceInds.append(SOPuids.index(RSOPuids[i]))
        else:
            msg = f'ReferencedSOPInstanceUID[{i}], {RSOPuids[i]}, is not in '\
                  + 'the list of SOPInstanceUIDs.'
            
            raise Exception(msg)
        
    return PFFGStoSliceInds







def GetPFFGStoSliceIndsBySeg(Seg, SOPuids):
    """
    Get a list (per segment) of the slice numbers that correspond to each 
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Inputs:
    ******
        
    Seg : Pydcom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
        
    PFFGStoSliceIndsBySeg : list of a list of integers
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-FrameFunctionalGroupsSequence in Seg.
    """
    
    PFFGStoSliceInds = GetPFFGStoSliceInds(Seg, SOPuids)
    
    DIVs = GetDIVs(Seg)
    
    PFFGStoSliceIndsBySeg = GroupListBySegment(ListToGroup=PFFGStoSliceInds, 
                                               DIVs=DIVs)
    
    return PFFGStoSliceIndsBySeg






    
    
def GetDIVs(Seg):
    """
    Get the DimensionIndexValues in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
        
    Outputs:
    *******
    
    DIVs : list of lists of integers
        List of DimensionIndexValues in Seg.
        
    """
    
    # Initialise the list of DIVs:
    DIVs = [] 
    
    sequences = Seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        DIV = sequence.FrameContentSequence[0]\
                      .DimensionIndexValues
        
        DIVs.append(DIV)
        
    return DIVs






def GroupListBySegment(ListToGroup, DIVs):
    """
    Group a list of items by common segment/ROI using the DimensionIndexValues.
    
    Inputs:
    ******
    
    ListToGroup : list
        The list to be grouped.
        
    DIVs : list of lists of integers
        List of DimensionIndexValues.
        
        
    Outputs:
    *******
    
    GroupedList : list of lists
        ListToGroup grouped by ROI.
        
        
    Note:
    ----
    
    DIVs is a list of lists of the form:
            
    [ [1, 10], [1, 11], [1, 12], ... [2, 6], [2, 7], [2, 8], ... ]
    
    where each list of index values that belong to a common ROI are separated
    into separate lists.
    
    The goal is to group a list (of anything) by ROI.  In the case of the DIVs,
    that is:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ],
      [ [2, 6], [2, 7], [2, 8], ... ],
      ... 
     ]
      
     
    If there's only one ROI, ListToGroup will still be nested in a list to 
    maintain continuity.  In the case of the DIVs with only one ROI:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ]
     ]
    """
    
    # Get the list of ROI numbers:
    RoiNums = []
    
    for DIV in DIVs:
        RoiNums.append(DIV[0])
        
    UniqueRoiNums = list(set(RoiNums))
    
    if len(UniqueRoiNums) == 1:
        GroupedList = [ListToGroup]
        
    else:
        GroupedList = []
        
        for RoiNum in UniqueRoiNums:
            ListThisRoi = [] # initialise
            
            for i in range(len(ListToGroup)):
                if DIVs[i][0] == RoiNum:
                    ListThisRoi.append(ListToGroup[i])
        
            GroupedList.append(ListThisRoi)
             
    return GroupedList











def AddToPFFGS(Seg):
    """
    Append the last item in the Per-FrameFunctionalGroupsSequence of the
    SEG to itself (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    *******
    
    NewSeg : Pydicom object
        Modified ROI object with Per-FrameFunctionalGroupsSequence lengthened
        by one.
    """
    
    from copy import deepcopy
    
    # Use Seg as a template for NewSeg: 
    NewSeg = deepcopy(Seg)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = deepcopy(Seg.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    NewSeg.PerFrameFunctionalGroupsSequence.append(last)
             
    return NewSeg






def AddToSS(Seg):
    """
    Append the last item in the SegmentSequence of the SEG Object to itself 
    (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ******
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    *******
    
    NewSeg : Pydicom object
        Modified Seg with SegmentSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Seg as a template for NewSeg: 
    NewSeg = deepcopy(Seg)
             
    # The last item in Segment Sequence:
    last = deepcopy(Seg.SegmentSequence[-1])
    
    # Append to the Segment Sequence:
    NewSeg.SegmentSequence.append(last)
             
    return NewSeg







def ChangeNumOfSegSequences(Seg, Dicoms, NumOfFrames, LogToConsole=False):
    """
    Adjust the number of sequences in a SEG object based so that they are 
    consistent with the number of DICOMs in the series and the number of
    frames in the pixel array.
         
    
    Inputs:
    ******
                              
    Seg : Pydicom object
        SEG object.
        
    Dicoms : List of Pydicom objects
        List of DICOM objects that relate to the SEG object to be created.
    
    NumOfFrames : integer
        The number of frames in the SEG object's pixel array.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
        
    Seg : Pydicom object
        SEG object with modified number of sequences.
    """
    
    
    """ 
    Modify the number of sequences in ReferencedSeriesSequence. 
    The number of sequences must be equal to the number of DICOMs.
    """
    
    from copy import deepcopy
    
    # The number of DICOMs:
    D = len(Dicoms)
    
    # The number of sequences in the ReferencedSeriesSequence:
    R = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
    
    
    # Add/remove sequences by appending/popping the last sequence N times, 
    # where N is the difference between D and R:
    if D > R:
        for i in range(D - R):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence\
               .append(deepcopy(Seg.ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1]))          
    else:
        for i in range(R - D):
            Seg.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()

    
    """ 
    Modify the number of sequences in SegmentSequence. 
    There must be only one sequence - for the (one) new segment.
    """
    
    S = len(Seg.SegmentSequence)
    
    if S > 1:
        for i in range(S - 1):
            Seg.SegmentSequence.pop()
    
    
    """ 
    Modify the number of sequences in Per-FrameFunctionGroupsSequence. 
    The number of sequences must be equal to NumOfFrames.
    """
    
    # The number of sequences in the Per-FrameFunctionGroupsSequence:
    P = len(Seg.PerFrameFunctionalGroupsSequence)
    
    if NumOfFrames > P:
        for i in range(NumOfFrames - P):
            Seg.PerFrameFunctionalGroupsSequence\
               .append(deepcopy(Seg.PerFrameFunctionalGroupsSequence[-1]))          
    else:
        for i in range(P - NumOfFrames):
            Seg.PerFrameFunctionalGroupsSequence.pop()

    
    
    
    if True:#LogToConsole:       
        R2 = len(Seg.ReferencedSeriesSequence[0].ReferencedInstanceSequence)
        
        S2 = len(Seg.SegmentSequence)
        
        P2 = len(Seg.PerFrameFunctionalGroupsSequence)
        
        print(f'\nThere are {D} DICOMs and there were {R} sequences in',
              f'ReferencedSeriesSequence. Now there are {R2} sequences.')
        
        print(f'\nThere were {S} sequences in SegmentSequence. Now there are',
              f'{S2} sequences.')
        
        print(f'\nThere are {NumOfFrames} frames in the pixel array and there',
              f'were {P} sequences in PerFrameFunctionalGroupsSequence. Now',
              f'there are {P2} sequences.')
        
    
    return Seg








def GetFrameFromPixArr(Seg, DicomDir, SearchString, SliceNum, LogToConsole):
    """
    Extract a single (2D) frame from a 3D pixel array that matches a given
    segment label.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    SearchString : string
        All or part of the Segment label containing the segmentation to be
        copied.
    
    SliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    Frame : Numpy array
        Frame from PixArr.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiLabels
    from DicomTools import GetRoiNum
    
    SegLabels = GetRoiLabels(Seg)
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = GetRoiNum(Seg, SearchString)
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    PFFGStoSliceIndsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
    PFFGStoSliceIndsInSeg = PFFGStoSliceIndsBySeg[SegNum]
    
    if not SliceNum in PFFGStoSliceIndsInSeg:
        msg = f'There is no segmentation in the segment matching the term '\
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








def GetFramesFromPixArr(Seg, DicomDir, SearchString, LogToConsole):
    """
    Extract all frames from a 3D pixel array that match a given segment label.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    SearchString : string
        All or part of the Segment Label of the segment containing the 
        segmentation to be copied.
    
    SliceNum : integer
        Slice index of the Source DICOM stack corresponding to the segmentation
        to be copied (counting from 0).
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    Frame : Numpy array
        Frame from PixArr.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    # Determine which segment contains the segmentation to be copied, and
    # the corresponding frame number in the pixel array:
    """ Note:  SegNum is equivalent to SegmentNumber in SegmentSequence. """
    SegNum = GetRoiNum(Seg, SearchString)
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    PFFGStoSliceIndsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
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
        
        Frame = GetFrameFromPixArr(Seg, DicomDir, SearchString, SliceNum, 
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








def GetPixArrInSeg(Seg, DicomDir, SearchString):
    """
    Get the frames in a SEG's pixel array that belong to a segment with
    specified label.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    SearchString : string
        All or part of the Segment Label containing the segment of interest.
                       
                            
    Outputs:
    *******
    
    PixArrInSeg : Numpy array
        Sub-array of the SEG's pixel array that contains the frames that belong
        to the segment of interest. 
        
    PFFGStoSliceIndsInSeg : List of integers
        List (for each segmentation) of the slice numbers that correspond to 
        each frame in PixArr.
    """
    
    import numpy as np
    
    AllPixArr = Seg.pixel_array
    
    AllF, R, C = AllPixArr.shape
    
    # Get the frame numbers of the SEG's pixel array that correspond to the 
    # segment of interest, and the corresponding Per-frame Functional Groups
    # Sequence-to-slice indices:
    FrameNumsInSeg, PFFGStoSliceIndsInSeg = GetFrameNums(Seg, SearchString, 
                                                         DicomDir)
    
    F = len(FrameNumsInSeg)
    
    PixArrInSeg = np.zeros((F, R, C), dtype='uint')
    
    for i in range(F):  
        PixArrInSeg[i] = AllPixArr[FrameNumsInSeg[i]]
    
        
    return PixArrInSeg, PFFGStoSliceIndsInSeg







def GetPixArrBySeg(Seg, DicomDir, LogToConsole=False):
    """
    Get a list of pixel arrays grouped by segment in a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
        
    F2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in PixArrBySeg.
    """
    
    import numpy as np
    from DicomTools import GetDicomSOPuids
    
    if LogToConsole:
        from GeneralTools import PrintIndsByRoi, PrintPixArrShapeBySeg
        
        print('\n\n', '-'*120)
        print('Results of GetPixArrBySeg():')
        
    AllPixArr = Seg.pixel_array
    
    """ If the Seg has a single frame it will have shape (NumOfRows, NumOfCols).
    For compatibility with multi-framed pixel arrays, convert single-framed
    pixel arrays to shape (1, NumOfRows, NumOfCols). """
    
    if len(AllPixArr.shape) < 3:
        R, C = AllPixArr.shape
        
        AllPixArr = np.reshape(AllPixArr, (1, R, C))
    
    AllF, R, C = AllPixArr.shape
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    F2SindsBySeg = GetPFFGStoSliceIndsBySeg(Seg, SOPuids)
    
    Nsegs = len(F2SindsBySeg)
    
    if LogToConsole:
        print(f'   AllPixArr.shape = {AllPixArr.shape}')
        #print('   F2SindsBySeg =')
        PrintIndsByRoi(F2SindsBySeg)
        
    PixArrBySeg = []
    
    j = 0 # total frame counter for all frames in AllPixArr
    
    for s in range(Nsegs):
        # The number of frames in AllPixArr that belong to this segment:
        F = len(F2SindsBySeg[s])
        
        # Initialise the pixel array for this segment:
        PixArr = np.zeros((F, R, C), dtype='uint')
        
        if F > 0:
            for i in range(F):
                PixArr[i] = AllPixArr[j]
                
                j += 1 # increment the total frame counter
        
        PixArrBySeg.append(PixArr)
    
    if LogToConsole:
        PrintPixArrShapeBySeg(PixArrBySeg)
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg







def GetSegDataOfInterest(Seg, FromSliceNum, FromSegLabel, DicomDir, 
                         LogToConsole=False):
    """
    Get data of interest from a SEG.
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object from an SEG file.
        
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a 
        single segmentation). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
        
    FromSegLabel : string
        All or part of the Source segment Description of the segment containing 
        the segmentation(s) to be copied.
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment.  The pixel array is a sub-array of the (entire)
        SEG's pixel array. 
        
    F2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in PixArrBySeg.
        
        
    Notes:
    *****
    
    There are 3 possible main use cases:
        1. Copy a single segmentation
        2. Copy an entire segment
        3. Copy all segments
        
    Which main case applies depends on the inputs:
        FromSliceNum
        FromSegLabel
    
    If FromSliceNum != None (i.e. if a slice number is defined) a single 
    segmentation will be copied from the segment that matches FromSegLabel 
    (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple segments and 
    FromSegLabel = None.  Possible outcomes:

        - Raise exception with message "There are multiple segments so the 
        description of the desired segment must be provided"
        - Copy the segmentation to all segments that have a segmentation on 
        slice FromSliceNum (* preferred option?)
        
    If FromSliceNum = None but FromSegLabel != None, all segmentations within 
    the segment given by FromSegLabel will be copied (--> Main Use Case 2).  
    
    If FromSliceNum = None and FromSegLabel = None, all segmentations within 
    all segments will be copied (--> Main Use Case 3).
    
    If FromSliceNum = None and FromSegLabel = None, all segments will be copied.  
    If not, get the segment number (zero-indexed) to be copied, and reduce 
    PixArrBySeg and F2SindsBySeg to the corresponding item.
    """
    
    import numpy as np
    from DicomTools import GetRoiNums, GetRoiLabels
    from copy import deepcopy
    
    if LogToConsole:
        from GeneralTools import PrintIndsByRoi, PrintPixArrShapeBySeg
        
        print('\n\n', '-'*120)
        print('Results of GetSegDataOfInterest():')
        print(f'   FromSegLabel = {FromSegLabel}')
        print(f'   FromSliceNum = {FromSliceNum}')
    
    AllPixArrBySeg, AllF2SindsBySeg = GetPixArrBySeg(Seg, DicomDir, 
                                                     LogToConsole)
    
    #print('\n\nAllPixArrBySeg =', AllPixArrBySeg)
    #print(f'\n\nAllF2SindsBySeg = {AllF2SindsBySeg}')
    
    Nsegs = len(AllF2SindsBySeg)
    
    """ Get the segment labels: """
    SegLabels = GetRoiLabels(Roi=Seg)
    
    """ The shape of the first pixel array: """
    F, R, C = AllPixArrBySeg[0].shape
    
    if LogToConsole:
        print(f'   SegLabels = {SegLabels}')
        #print('   AllF2SindsBySeg =')
        PrintIndsByRoi(AllF2SindsBySeg)
    
    """ Initialise variable that indicates if the SEG data was reduced by
    frame/slice number (FromSliceNum): """
    ReducedBySlice = False
            
    if FromSliceNum:
        """ Limit data to those which belong to the chosen slice number. """
        
        PixArrBySeg = []
        F2SindsBySeg = []
    
        for s in range(Nsegs):
            """ The pixel array and frame-to-slice indices for this segment:"""
            AllPixArr = deepcopy(AllPixArrBySeg[s])
            AllF2Sinds = deepcopy(AllF2SindsBySeg[s])
                
            """ Get the frame number(s) that relate to FromSliceNum: """
            FromFrmNums = [i for i, x in enumerate(AllF2Sinds) if x==FromSliceNum]
            
            if LogToConsole:
                print(f'\nAllF2Sinds = {AllF2Sinds}')
                print(f'FromFrmNums = {FromFrmNums}')
            
            if len(FromFrmNums) != len(AllF2Sinds):
                """ Some frames will be rejected: """
                ReducedBySlice = True
                
                """ Initialise PixArr with the reduced number of frames 
                matching len(FromFrmNums): """
                PixArr = np.zeros((len(FromFrmNums), R, C), dtype='uint')
                F2Sinds = []
                
                if FromFrmNums:
                    """ Keep only the indeces and frames that relate to 
                    FromFrmNums. """
                    #PixArr = [AllPixArr[f] for f in FromFrmNums]
                    
                    #""" Set PixArr to the first frame in FromFrmNums: """
                    #PixArr = AllPixArr[0]
                    #
                    #if len(FromFrmNums) > 1:
                    #    """ Add additional frames: """
                    #    for f in range(1, len(FromFrmNums)):
                    #        PixArr = np.vstack((PixArr, AllPixArr[f]))
                    
                    for f in range(len(FromFrmNums)):
                        PixArr[f] = AllPixArr[FromFrmNums[f]]
                        
                        F2Sinds.append(AllF2Sinds[FromFrmNums[f]])
                    
                    if LogToConsole:
                        print(f'F2Sinds = {F2Sinds}')
                    
                    """ Append non-empty pixel arrays and non-empty F2Sinds:"""
                    PixArrBySeg.append(PixArr)
                    F2SindsBySeg.append(F2Sinds)
                
                #else:
                #    #PixArr = []
                #    #PixArr = np.zeros((0, R, C), dtype='uint')
                #    F2Sinds = []
                
                #""" Append empty pixel arrays and empty F2Sinds: """
                #PixArrBySeg.append(PixArr)
                #F2SindsBySeg.append(F2Sinds)
            """ else all frames to remain. """
        
        #print('\n\nPixArrBySeg =', PixArrBySeg)
        
        if ReducedBySlice:
            """ Replace AllPixArrByRoi and AllF2SindsBySeg with PixArrBySeg  
            and F2SindsBySeg in case further restricting of data is required 
            below: """
            AllPixArrBySeg = deepcopy(PixArrBySeg)
            AllF2SindsBySeg = deepcopy(F2SindsBySeg)
        """ else AllPixArrBySeg and AllF2SindsBySeg remain unchanged. """
        
        #print('\n\nAllPixArrBySeg =', AllPixArrBySeg)
        
        if LogToConsole and ReducedBySlice:
            print('\n   After limiting data to those that relate to slice',
                  f'number {FromSliceNum}:')#, the F2SindsBySeg =')
            PrintIndsByRoi(F2SindsBySeg)
            print('   PixArrBySeg = ', PixArrBySeg)
            PrintPixArrShapeBySeg(PixArrBySeg)
    
    
    """ Initialise variable that indicates if the SEG data was reduced by
    chosen segment label (FromSegLabel): """
    ReducedBySeg = False
    
    if FromSegLabel:
        """ Limit data to those which belong to the chosen segment(s). """
        
        PixArrBySeg = []
        F2SindsBySeg = []
        
        """ Get the segment number(s) whose name matches FromSegLabel: """
        FromSegNums = GetRoiNums(Roi=Seg, SearchString=FromSegLabel)
        
        #print(f'\nFromSegNums = {FromSegNums}')
        #print(f'AllF2SindsBySeg = {AllF2SindsBySeg}')
        
        """ 10/02:  Not sure about this if statement... """
        #if len(FromSegNums) != len(AllF2SindsBySeg):
        #    #ReducedBySeg = True
        #
        #    """ Get the names of all ROIs: """
        #    FromSegNames = GetRoiLabels(Roi=Seg)
        #
        #    """ Limit the list of segment descriptions to those that belong to 
        #    the chosen segment(s): """
        #    FromSegNames = [FromSegNames[i] for i in FromSegNums]
        #    
        #    
        #    PixArrBySeg = [AllPixArrBySeg[s] for s in FromSegNums]
        #    F2SindsBySeg = [AllF2SindsBySeg[s] for s in FromSegNums]
        #    
        #    if LogToConsole:
        #        print('\n   After limiting data to those whose segment name',
        #              f'matches {FromSegNames}:')#, the F2SindsBySeg =')
        #        PrintIndsByRoi(F2SindsBySeg)
        #        PrintPixArrShapeBySeg(PixArrBySeg)
        #        
        #else:
        #    PixArrBySeg = deepcopy(AllPixArrBySeg)
        #    F2SindsBySeg = deepcopy(AllF2SindsBySeg)
        
        
    
        """ Limit the list of segment descriptions to those that belong to 
        the chosen segment(s): """
        FromSegLabels = [SegLabels[i] for i in FromSegNums]
        
        PixArrBySeg = [AllPixArrBySeg[s] for s in FromSegNums]
        F2SindsBySeg = [AllF2SindsBySeg[s] for s in FromSegNums]
        
        if len(F2SindsBySeg) != len(AllF2SindsBySeg):
            ReducedBySeg = True
        
        if LogToConsole and ReducedBySeg:
            print('\n   After limiting data to those whose segment name',
                  f'matches {FromSegLabels}:')#, the F2SindsBySeg =')
            PrintIndsByRoi(F2SindsBySeg)
            PrintPixArrShapeBySeg(PixArrBySeg)
    
    if LogToConsole:
        print('\n   Final outputs of GetSegDataOfInterest():')
        PrintPixArrShapeBySeg(PixArrBySeg)
        PrintIndsByRoi(F2SindsBySeg)
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg







def AddCopiedPixArr(OrigPixArr, OrigPFFGStoSliceInds, 
                    PixArrToAdd, PFFGStoSliceIndsToAdd):
    """
    Add frame(s) from a pixel array to an existing pixel array.
       
    
    Inputs:
    ******
    
    OrigPixArr : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames that belong to
        the segment of interest. 
        
    OrigPFFGStoSliceInds : list of integers
        List (for each segmentation) of slice numbers that correspond to each 
        frame in OrigPixArr.
        
    PixArrToAdd : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames that originated
        from a single segmentation.  F will be one if a single segmentation is
        to be copied, but may be greater than one, if, for example, the 
        segmentation to be copied was resampled to a number of segmentations.
        
    PFFGStoSliceIndsToAdd : list of integers
        List (for each segmentation) of slice numbers that correspond to each 
        frame in PixArrToAdd.
        
        
    Outputs:
    *******
    
    NewPixArr : Numpy array
        A FxRxC (Frames x Rows x Columns) array of all frames in 
        OrigPixArrInSeg and PixArrToAddToSeg. Note that F will not necessarily
        be the sum of the number of frames in OrigPixArrInSeg and 
        PixArrToAddToSeg.  Some frames in both pixel arrays may correspond
        to the same slice number, in which case, any duplicate frame(s) in
        OrigPixArrInSeg will be over-written with the frame(s) in
        PixArrToAddToSeg that correspond to the same slice number.
        
    NewPFFGStoSliceInds : list of integers
        Modified list (for each segmentation) of slice numbers that correspond   
        to each frame in NewPixArr.
    """
    
    import numpy as np
    from GeneralTools import UniqueItems
    
    # Combine OrigPFFGStoSliceInds and PFFGStoSliceIndsToAdd:
    NewPFFGStoSliceInds = OrigPFFGStoSliceInds + PFFGStoSliceIndsToAdd
    
    # Remove any duplicate slice numbers:
    NewPFFGStoSliceInds = UniqueItems(NewPFFGStoSliceInds)
    
    F = len(NewPFFGStoSliceInds)
    
    OrigF, R, C = OrigPixArr.shape
    
    NewPixArr = np.zeros((F, R, C), dtype='uint')
    
    for FrameNum in range(F):
        # The corresponding slice number for this frame:
        SliceNum = NewPFFGStoSliceInds[FrameNum]
        
        if SliceNum in PFFGStoSliceIndsToAdd:
            """
            This frame is to be over-written by the corresponding frame in
            PixArrToAdd.
            """
            
            ind = PFFGStoSliceIndsToAdd.index(SliceNum)
            
            NewPixArr[FrameNum] = PixArrToAdd[ind]
            
        else:
            """
            This existing frame (in OrigPixArr) is to be preserved.
            """
            
            ind = OrigPFFGStoSliceInds.index(SliceNum)
            
            NewPixArr[FrameNum] = OrigPixArr[ind]
    
    
    return NewPixArr, NewPFFGStoSliceInds











def CreateSeg(SrcSeg, TrgSeg, TrgPixArrBySeg, TrgF2SindsBySeg, TrgDicomDir, 
              FromSegLabel, AddTxtToSegLabel='', RefAllSOPs=False,
              LogToConsole=False):
    """
    Create an SEG object for the target dataset.
         
    
    Inputs:
    ******
    
    SrcSeg : Pydicom object
        The Source SEG object. If TrgSeg = None, SrcSeg will be used as a 
        template for NewTrgSeg.
    
    TrgSeg : Pydicom object or None
        If TrgSeg != None, the TrgSeg will be modified to create NewTrgSeg.
        If TrgSeg = None, the SrcSeg will be used as a template.

    TrgPixArrBySeg : list of Numpy arrays
        List (for each segment) of the pixel array containing the frames that
        belong to each segment for Target.  The pixel array is a sub-array of 
        the (entire) SEG's pixel array. 
        
    TrgF2SindsBySeg : List of a list of integers
        List (for each segment) of the slice numbers that correspond to each
        frame in each pixel array in TrgPixArrBySeg.
                    
    TrgDicomDir : string
        Directory containing the Target DICOMs.
    
    FromSegLabel : string
        All or part of the Source SegmentLabel of the segment containing the 
        segmentations(s) that were copied.
    
    AddTxtToSegLabel : string (optional; '' by default)
        String to be added to SeriesDescription for the new SEG. 
    
    RefAllSOPs : boolean (optional; False by default)
        If True, all SOPInstanceUIDs in the series will be referenced in
        ReferencedSeriesSequence[0].ReferencedInstanceSequence. If False, only
        the SOPInstanceUIDs that have segmentations will be referenced.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    NewTrgSeg : Pydicom object
        New SEG object.
    """
    
    from pydicom.uid import generate_uid
    import time
    from copy import deepcopy
    import numpy as np
    from GeneralTools import FlattenList, UniqueItems
    from DicomTools import ImportDicoms, GetRoiNums#, GetRoiLabels
    from ImageTools import GetImageAttributes
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of CreateSeg():')
        print(f'   TrgF2SindsBySeg = {TrgF2SindsBySeg}')
        print(f'   len(TrgPixArrBySeg) = {len(TrgPixArrBySeg)}')
            
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    NumOfDicoms = len(TrgDicoms)
    
    NumOfFramesBySeg = []
    for s in range(len(TrgPixArrBySeg)):
        NumOfFramesBySeg.append(TrgPixArrBySeg[s].shape[0])
    
    #print(f'\n\n\nNumOfFramesBySeg = {NumOfFramesBySeg}')
    
    NumOfFrames = sum(NumOfFramesBySeg)
    
    #NumOfSegs = len(TrgF2SindsBySeg)
    
    """ Get the list of SegmentNumber for the segments of interest in SrcSeg. """
    SegNums = GetRoiNums(SrcSeg, FromSegLabel)
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=TrgDicomDir, Package='sitk')
    
    UniqueF2Sinds = UniqueItems(Items=TrgF2SindsBySeg, IgnoreZero=False, 
                                MaintainOrder=True)
    
    FlattenedF2Sinds = FlattenList(TrgF2SindsBySeg)
    
    if LogToConsole:
        print(f'   NumOfDicoms = {NumOfDicoms}')
        print(f'   SegNums = {SegNums}')
        print(f'   NumOfFramesBySeg = {NumOfFramesBySeg}')
        print(f'   NumOfFrames = {NumOfFrames}')
    
    # Use TrgSeg or SrcSeg as a template for NewTrgSeg:
    if TrgSeg:
        NewTrgSeg = deepcopy(TrgSeg)
    else:
        NewTrgSeg = deepcopy(SrcSeg)
    
    # Generate a new SOPInstance UID:
    NewSOPuid = generate_uid()
    NewTrgSeg.SOPInstanceUID = deepcopy(NewSOPuid)
    NewTrgSeg.file_meta.MediaStorageSOPInstanceUID = deepcopy(NewSOPuid)
    
    # Generate a new SeriesInstance UID:
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    """ If TrgSeg != None, some tags will not need to be replaced.
    If TrgSeg = None, use the corresponding values in the first Target 
    DICOM. """
    
    if TrgSeg == None:
        NewTrgSeg.StudyDate = TrgDicoms[0].StudyDate
        NewTrgSeg.SeriesDate = TrgDicoms[0].SeriesDate
        #NewTrgSeg.ContentDate = TrgDicoms[0].ContentDate
        NewTrgSeg.ContentDate = NewDate
        NewTrgSeg.StudyTime = TrgDicoms[0].StudyTime
        NewTrgSeg.SeriesTime = TrgDicoms[0].SeriesTime
        #NewTrgSeg.ContentTime = TrgDicoms[0].ContentTime
        NewTrgSeg.ContentTime = NewTime
        NewTrgSeg.Manufacturer = TrgDicoms[0].Manufacturer
        NewTrgSeg.SeriesDescription += AddTxtToSegLabel
        NewTrgSeg.PatientName = TrgDicoms[0].PatientName
        NewTrgSeg.PatientID = TrgDicoms[0].PatientID
        NewTrgSeg.PatientBirthDate = TrgDicoms[0].PatientBirthDate
        NewTrgSeg.PatientSex = TrgDicoms[0].PatientSex
        NewTrgSeg.PatientAge = TrgDicoms[0].PatientAge
        NewTrgSeg.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
        NewTrgSeg.StudyID = TrgDicoms[0].StudyID
        NewTrgSeg.SeriesNumber = TrgDicoms[0].SeriesNumber
        NewTrgSeg.InstanceNumber = TrgDicoms[0].InstanceNumber
        NewTrgSeg.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
        NewTrgSeg.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
        NewTrgSeg.Rows = TrgDicoms[0].Rows
        NewTrgSeg.Columns = TrgDicoms[0].Columns
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PlaneOrientationSequence[0]\
                 .ImageOrientationPatient = TrgDicoms[0].ImageOrientationPatient
           
        """ 
        Notes:
            
        1. SliceThickness in PixelMeasuresSequence appears to be the z-Spacing 
        rather than SliceThickness from the DICOM metadata. 
        
        2. The character limit of VR DS is 16 characters.
        """
        Zspacing = f"{Spacings[2]}"
        
        if len(Zspacing) > 16:
            Zspacing = Zspacing[:16]
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SliceThickness = deepcopy(Zspacing)
        
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .SpacingBetweenSlices = deepcopy(Zspacing)
    
        NewTrgSeg.SharedFunctionalGroupsSequence[0]\
                 .PixelMeasuresSequence[0]\
                 .PixelSpacing = TrgDicoms[0].PixelSpacing
    
    
    #""" Modify SeriesDescription. <-- this was done above """
    #
    #if len(SegNums) == 1:
    #    """ There is only one segment. The new Target SeriesDescription will be
    #    the Source SeriesDescription with the Source SegmentLabel that was 
    #    copied + AddTxtToSegLabel. """
    #    NewSD = SrcSeg.SeriesDescription + '_' \
    #            + SrcSeg.SegmentSequence[SegNums[0]].SegmentLabel \
    #            + AddTxtToSegLabel
    #            
    #    NewTrgSeg.SeriesDescription = NewSD
    #else:
        
    
    
    """ Modify NumberOfFrames. """
    
    NewTrgSeg.NumberOfFrames = f"{len(FlattenedF2Sinds)}"
    
    
    """ Modify ReferencedSeriesSequence. """
    
    NewTrgSeg.ReferencedSeriesSequence[0]\
             .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    """ Modify ReferencedInstanceSequence in ReferencedSeriesSequence. 
    Note:
        If TrgSeg was used as a template no changes should need to be made. But 
        since it's not computationally expensive to do so, ensure that all 
        TrgDicoms are referenced anyhow. """
    
    """ Loop through each TrgDicom. 
    14/02: Appending sequences to reference every SOP UID is incredibly time
    consuming so instead will only reference the SOP UIDs that have 
    segmentations. 
    
    At 21:20 on 14/02 changed from looping through range(NumOfFrames) to
    range(len(UniqueF2Sinds)), since some slices are repeated.
    
    15/02/21: Variable RefAllSOPs allows for both scenarios.
    """
    
    """ Define the required number of ReferencedInstanceSequence and slice
    indices that correspond to the SOPInstanceUIDs to reference. """
    if RefAllSOPs:
        ReqNumOfRIS = deepcopy(NumOfDicoms)
        ReqIndsToRef = list(range(NumOfDicoms))
    else:
        ReqNumOfRIS = len(UniqueF2Sinds)
        ReqIndsToRef = deepcopy(UniqueF2Sinds)
        
    for i in range(ReqNumOfRIS):#range(len(UniqueF2Sinds)):#range(NumOfFrames):#range(NumOfDicoms):
        #s = FlattenedF2Sinds[i]
        #s = UniqueF2Sinds[i]
        s = ReqIndsToRef[i]
        
        """ The number of sequences in ReferencedSeriesSequence: """
        RSS = deepcopy(NewTrgSeg.ReferencedSeriesSequence[0]\
                                .ReferencedInstanceSequence)
        
        N = len(RSS)
        
        if i > N - 1:
            """ Increase the sequence by one. """
            LastItem = deepcopy(RSS[-1])
            
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.append(LastItem)
        
        """ Update the sequence. """
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPClassUID = TrgDicoms[s].SOPClassUID
        
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
    
    """ Check if there are more sequences in ReferencedInstanceSequence than 
    required. """
    NumOfRIS = len(NewTrgSeg.ReferencedSeriesSequence[0]\
                            .ReferencedInstanceSequence)
    
    if NumOfRIS > ReqNumOfRIS:#len(UniqueF2Sinds):#NumOfFrames:#NumOfDicoms:
        for i in range(NumOfRIS - ReqNumOfRIS):#len(UniqueF2Sinds)):#NumOfFrames):#NumOfDicoms):
            """ Remove the last sequence. """
            NewTrgSeg.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence.pop()
    
    
    
    """ Modify SegmentSequence. """
    
    """ Loop through all SegNums. """
    for s in range(len(SegNums)):
        SegNum = SegNums[s]
        
        """ Replace the s^th SegmentSequence in NewTrgSeg with the SegNum^th 
        SegmentSequence in SrcSeg."""
        NewTrgSeg.SegmentSequence[s] = deepcopy(SrcSeg.SegmentSequence[SegNum])
        
        NewTrgSeg.SegmentSequence[s].SegmentNumber = int(s+1)
        
        NewTrgSeg.SegmentSequence[s]\
                 .SegmentLabel = SrcSeg.SegmentSequence[SegNum].SegmentLabel
        
    """ Check if there are more sequences in SegmentSequence than required. """
    N = len(NewTrgSeg.SegmentSequence)
    
    if N > len(SegNums):
        for i in range(N - len(SegNums)):
            """ Remove the last sequence: """
            NewTrgSeg.SegmentSequence.pop()
    
    
    """ Modify PerFrameFunctionalGroupsSequence. """
    
    """ Get the ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence."""
    RefSOPsInRIS = GetRSOPuidsInRIS(NewTrgSeg)
    
    ReqNumOfPFFGS = len(FlattenedF2Sinds)
    
    i = 0 # total frame counter (for all segments)
    
    for s in range(len(TrgF2SindsBySeg)):
        for f in range(len(TrgF2SindsBySeg[s])):
            PFFGS = deepcopy(NewTrgSeg.PerFrameFunctionalGroupsSequence)
            
            N = len(PFFGS)
            
            if i > N - 1:
                """ Increase the sequence by one. """
                LastItem = deepcopy(PFFGS[-1])
                
                NewTrgSeg.PerFrameFunctionalGroupsSequence.append(LastItem)
            
            """ The DICOM slice number: """
            d = TrgF2SindsBySeg[s][f]
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPClassUID = TrgDicoms[d].SOPClassUID
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .DerivationImageSequence[0]\
                     .SourceImageSequence[0]\
                     .ReferencedSOPInstanceUID = TrgDicoms[d].SOPInstanceUID
            
            """ The index of the ReferencedSOPInstanceUID in 
            ReferencedInstanceSequence is d, since the sequence was overwritten
            in the same order as the list of DICOMs in TrgDicoms. 
            Note:
                While s (segment number) and d (DICOM slice number) are both
                0-indexed, DimensionIndexValues are 1-indexed.
            
            14/02/21: Following change to no longer reference all SOP UIDs 
            instead of mirroring the sequence of ReferencedSOPInstanceUIDs in
            ReferencedInstanceSequence (which matched the order of DICOMs), and
            hence, the second element in DimensionIndexValues was d + 1, now
            the second value will simply be f + 1.
            
            At 21:20 on 14/02:  It's wrong to use f+1 since there will be more
            frames than the number of ReferencedSOPInstanceUIDs if there are
            more than one segmentation on any given frame. So instead need to 
            get the index of the ReferencedSOPInstanceUID in the 
            ReferencedInstanceSequence that matches the d^th DICOM 
            SOPInstanceUID.
            
            Alternatively since UniqueF2Sinds was created maintaining the order
            of indices in TrgF2SindsBySeg, the index (+1) of the UniqueF2Sind 
            that matches slice number d will be the required second element.
            
            15/02: Now the index is the index of slice d within ReqIndsToRef.
            """
            
            #ind = UniqueF2Sinds.index(d)
            #ind = ReqIndsToRef.index(d)
            
            """ Some SEGs are not rendered at all so try more fail-proof
            method. """
            
            SOPuid = TrgDicoms[d].SOPInstanceUID
            
            """ Get the index of SOPuid in the ReferencedSOPInstanceUIDs from
            ReferencedInstanceSequence (RefSOPsInRIS): """
            ind = RefSOPsInRIS.index(SOPuid)
            
        
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .FrameContentSequence[0]\
                     .DimensionIndexValues = [s + 1, ind + 1]# [s + 1, f + 1] # [s + 1, d + 1]
            
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .PlanePositionSequence[0]\
                     .ImagePositionPatient = TrgDicoms[d].ImagePositionPatient
               
            NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                     .SegmentIdentificationSequence[0]\
                     .ReferencedSegmentNumber = s + 1
            
            i += 1 # increment the total frame count
    
    """ Check if there are more sequences in PerFrameFunctionalGroupsSequence 
    than required. """
    NumOfPFFGS = len(NewTrgSeg.PerFrameFunctionalGroupsSequence)
    
    if NumOfPFFGS > ReqNumOfPFFGS:
        for i in range(NumOfPFFGS - ReqNumOfPFFGS):
            """ Remove the last sequence: """
            NewTrgSeg.PerFrameFunctionalGroupsSequence.pop()        
    
    
    """ Modify PixelData. """
    
    PixArr = deepcopy(TrgPixArrBySeg[0])
    
    if len(TrgPixArrBySeg) > 1:
        """ Vertically stack all pixel arrays into a single pixel array. """
        for s in range(1, len(TrgPixArrBySeg)):
            PixArr = np.vstack((PixArr, TrgPixArrBySeg[s]))
    
    """ Ravel (to 1D array) PixArr and pack bits (see
    https://github.com/pydicom/pydicom/issues/1230). """
    packed = pack_bits(PixArr.ravel())
    
    """ Convert PixArr to bytes. """
    #NewTrgSeg.PixelData = PixArr.tobytes() <-- doesn't work
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    
    if LogToConsole:
        print('-'*120)
        
    return NewTrgSeg







def ErrorCheckSeg(Seg, DicomDir, RefAllSOPs=False, LogToConsole=False):
    """
    Check a SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the SEG.  
    
    Inputs:
    ******
    
    Seg : Pydicom object
        SEG object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Seg.
    
    RefAllSOPs : boolean (optional; False by default)
        If True, all SOPInstanceUIDs in the series were referenced in
        ReferencedSeriesSequence[0].ReferencedInstanceSequence. If False, only
        the SOPInstanceUIDs that have segmentations were referenced.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe any errors that are found.  If no 
        errors are found an empty list ([]) will be returned.
        
    Nerrors : integer
        The number of errors found.
    """
    
    import importlib
    import GeneralTools
    importlib.reload(GeneralTools)
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import ImportDicom
    from ImageTools import GetImageAttributes
    from GeneralTools import AreItemsEqualToWithinEpsilon
    from GeneralTools import AreListsEqualToWithinEpsilon
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ErrorCheckSeg():')
        
    # The maximum allowed error between two lists/items:
    epsilon = 1e-05
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='pydicom')
    
    PixArr = Seg.pixel_array
    
    """ The number of frames, rows and columns in the SEG's pixel array: """
    if len(PixArr.shape) > 2:
        NumOfFrames, NumOfRows, NumOfCols = PixArr.shape
    else:
        NumOfRows, NumOfCols = PixArr.shape
        NumOfFrames = 1
    
    #NumOfDicoms = len(SOPuids)
    
    if LogToConsole:
        print(f'NumOfFrames = {NumOfFrames}')
    
    LogList = []
    Nerrors = 0
    
    
    """Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID."""
    MSSOPuid = deepcopy(Seg.file_meta.MediaStorageSOPInstanceUID)
    SOPuid = deepcopy(Seg.SOPInstanceUID)
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} matches '\
              + f'SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not '\
              + f'match SOPInstanceUID {SOPuid}.\n'
        Nerrors = Nerrors + 1 
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether the number of sequences in 
    ReferencedInstanceSequence matches the number of DICOMs.
    14/02/21: Due to the excessive time taken to append sequences there was a
    change to only reference the SOP UIDs that have segmentations. 
    15/02/21: Variable RefAllSOPs allows for both scenarios. """
    
    RIS = deepcopy(Seg.ReferencedSeriesSequence[0]\
                      .ReferencedInstanceSequence)
    
    if RefAllSOPs:
        if len(RIS) == len(SOPuids):
            msg = f'INFO:  The number of SOPInstanceUIDs, {len(SOPuids)},'\
                  + ' matches the number of sequences in '\
                  + f'ReferencedInstanceSequence, {len(RIS)}.\n'
        else:
            msg = f'ERROR:  The number of SOPInstanceUIDs, {len(SOPuids)}, '\
                  + 'does not match the number of sequences in '\
                  + f'ReferencedInstanceSequence, {len(RIS)}.\n'
            LogList.append(msg) 
            Nerrors = Nerrors + 1
        if LogToConsole:
            print(msg)
    
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs do not
    match the SOPInstanceUIDs."""        
    RefSOPuidsInRIS = [RIS[i].ReferencedSOPInstanceUID for i in range(len(RIS))]

    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInRIS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
    if Inds:
        for i in range(len(Inds)):
            uid = RefSOPuidsInRIS[Inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + f'ReferencedInstanceSequence[{Inds[i]}] does not match '\
                  + 'any SOPInstanceUID.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
            if LogToConsole:
                print(msg)
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in '\
                  + f'ReferencedInstanceSequence matches the '\
                  + 'SOPInstanceUIDs.\n'
        LogList.append(msg)
        if LogToConsole:
            print(msg)
                
                
    
    """Determine whether any of the SOPInstanceUIDs are not referenced in the
    ReferencedSOPInstanceUIDs.
    14/02/21: Following the change to not reference every SOP UID this check is
    no longer required: 
    15/02/21: Variable RefAllSOPs allows for both scenarios. """        
    
    if RefAllSOPs:
        IsMatch = [SOPuid in RefSOPuidsInRIS for SOPuid in SOPuids]
        
        #NumOfMatches = IsMatch.count(True)
        
        # Find the indices of any non-matching SOP UIDs:
        Inds = [i for i, x in enumerate(IsMatch) if x==False]
        
        if Inds:        
            for i in range(len(Inds)):
                uid = SOPuids[Inds[i]]
                
                msg = f'ERROR:  SOPInstanceUID {uid} is not referenced in any '\
                      + 'ReferencedSOPInstanceUID in '\
                      + 'ReferencedInstanceSequence.\n'
                LogList.append(msg)
                Nerrors = Nerrors + 1
                if LogToConsole:
                    print(msg)
        else:
            msg = f'INFO:  All SOPInstanceUIDs are referenced in the '\
                      + 'ReferencedSOPInstanceUIDs in '\
                      + 'ReferencedInstanceSequence.\n'        
            LogList.append(msg)
            if LogToConsole:
                print(msg)
    
    
                           
    if Seg.ReferencedSeriesSequence[0].SeriesInstanceUID == Dicom.SeriesInstanceUID:
        msg = f'INFO:  SEG SeriesInstanceUID in ReferencedSeriesSequence, '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, '\
              + 'matches DICOM SeriesInstanceUID, '\
              + f'{Dicom.SeriesInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG SeriesInstanceUID in ReferencedSeriesSequence, '\
              + f'{Seg.ReferencedSeriesSequence[0].SeriesInstanceUID}, does '\
              + 'not match DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
        Nerrors = Nerrors + 1 
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
       
    if Seg.PatientName == Dicom.PatientName:
        msg = f'INFO:  SEG PatientName, {Seg.PatientName}, matches DICOM '\
              + f'PatientName, {Dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  SEG PatientName, {Seg.PatientName}, does not match '\
              + f'DICOM PatientName, {Dicom.PatientName}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.PatientID == Dicom.PatientID:
        msg = f'INFO:  SEG PatientID, {Seg.PatientID}, matches the '\
              + f'DICOM PatientID, {Dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  SEG PatientID, {Seg.PatientID}, does not match the '\
              + f'DICOM PatientID, {Dicom.PatientID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.StudyInstanceUID == Dicom.StudyInstanceUID:
        msg = f'INFO:  SEG StudyInstanceUID, {Seg.StudyInstanceUID}, '\
              + 'matches the DICOM StudyInstanceUID, '\
              + f'{Dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  SEG StudyInstanceUID, {Seg.StudyInstanceUID}, does '\
              + 'not match the DICOM StudyInstanceUID, '\
              + f'{Dicom.StudyInstanceUID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    #if Seg.StudyID == Dicom.StudyID:
    #    msg = f'INFO:  The SEG StudyID {Seg.StudyID} matches the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #else:
    #    msg = f'ERROR:  The SEG StudyID {Seg.StudyID} does not match the'\
    #          + f' DICOM StudyID {Dicom.StudyID}.'
    #
    #    Nerrors = Nerrors + 1
    #    
    #LogList.append(msg)
    #    
    #if LogToConsole:
    #    print(msg)
            
    if Seg.FrameOfReferenceUID == Dicom.FrameOfReferenceUID:
        msg = f'INFO:  SEG FrameOfReferenceUID, {Seg.FrameOfReferenceUID}, '\
              + 'matches DICOM FrameOfReferenceUID, '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  SEG FrameOfReferenceUID, {Seg.FrameOfReferenceUID}, '\
              + 'does not match DICOM FrameOfReferenceUID, '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if int(Seg.NumberOfFrames) == NumOfFrames:
        msg = f'INFO:  NumberOfFrames, {Seg.NumberOfFrames}, matches the '\
              + f'number of frames in PixelData, {NumOfFrames}.\n'
    else:
        msg = f'ERROR:  NumberOfFrames, {Seg.NumberOfFrames}, does not '\
              + f'match the number of frames in PixelData, {NumOfFrames}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.Rows == Dicom.Rows:
        msg = f'INFO:  SEG Rows, {Seg.Rows}, matches DICOM Rows, '\
              + f'{Dicom.Rows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {Seg.Rows}, does not match DICOM Rows, '\
              + f'{Dicom.Rows}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Seg.Columns == Dicom.Columns:
        msg = f'INFO:  SEG Columns, {Seg.Columns}, matches DICOM '\
              + f'Columns, {Dicom.Columns}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {Seg.Columns}, does not match DICOM '\
              + f'Columns, {Dicom.Columns}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """ This is already checked below
    if Seg.Rows == NumOfRows:
        msg = f'INFO:  SEG Rows, {Seg.Rows}, matches the number of rows in '\
              + f'PixelData, {NumOfRows}.\n'
    else:
        msg = f'ERROR:  SEG Rows, {Seg.Rows}, does not match the number of '\
              + f'rows in PixelData, {NumOfRows}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
        
    if Seg.Columns == NumOfCols:
        msg = f'INFO:  SEG Columns, {Seg.Columns}, matches the number of '\
              + f'columns in PixelData, {NumOfCols}.\n'
    else:
        msg = f'ERROR:  SEG Columns, {Seg.Columns}, does not match the '\
              + f'number of columns in PixelData, {NumOfCols}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    """
    
    
    """Keep track of the number of sequences in SegmentSequence."""
    NumSS = len(Seg.SegmentSequence)
    
    """Verify that the SegmentNumber in each SegmentSequence increments from 
    1."""
    SS = deepcopy(Seg.SegmentSequence)
    
    for i in range(NumSS):
        N = int(SS[i].SegmentNumber)
        
        #print(f'\n\nN = {N}, type(N) = {type(N)}')
        
        if N == i + 1:
            msg = f'INFO:  SegmentNumber in SegmentSequence[{i}] is {N}.\n'
            LogList.append(msg)
        else:
            msg = f'ERROR:  SegmentNumber in SegmentSequence[{i}] is {N}.'\
                  + f' It should be {i + 1}.\n'
            Nerrors = Nerrors + 1
            LogList.append(msg)

    if LogToConsole:
        print(msg)
    
    
    """Check various tags in SharedFunctionalGroupsSequence."""
    SegIOP = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PlaneOrientationSequence[0]\
                         .ImageOrientationPatient)
    
    SegIOP = [float(item) for item in SegIOP]
    
    DcmIOP = [float(item) for item in Dicom.ImageOrientationPatient]
    
    if SegIOP == DcmIOP:
        msg = 'INFO:  SEG ImageOrientationPatient in '\
              + f'SharedFunctionalGroupsSequence, {SegIOP} matches '\
              + f'DICOM ImageOrientationPatient, {DcmIOP}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegIOP, DcmIOP, epsilon):
            msg = 'INFO:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence, {SegIOP}, is within '\
                  + f'{epsilon} of DICOM ImageOrientationPatient, {DcmIOP}.\n'
        else:
            msg = 'ERROR:  SEG ImageOrientationPatient in '\
                  + f'SharedFunctionalGroupsSequence, {SegIOP}, is not within'\
                  + f' {epsilon} of DICOM ImageOrientationPatient, {DcmIOP}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
       
    
    SegST = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .SliceThickness)
    
    """ The SliceThickness appears to be the z-Spacing rather than the
    SliceThickness from the DICOM metadata. """
    
    if float(SegST) == Spacings[2]:
        msg = 'INFO:  SEG SliceThickness in SharedFunctionalGroupsSequence,'\
              + f' {SegST}, matches the slice thickness calculated '\
              + 'from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient, {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegST), Spacings[2], epsilon):
            msg = 'INFO:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence, {SegST}, is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM'\
                  + ' ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
        else:
            print(f'\n{float(SegST) - Spacings[2]}')
            
            msg = 'ERROR:  SEG SliceThickness in '\
                  + f'SharedFunctionalGroupsSequence, {SegST}, is not within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    SegSBS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                         .PixelMeasuresSequence[0]\
                         .SpacingBetweenSlices)
    
    if float(SegSBS) == Spacings[2]:
        msg = 'INFO:  SEG SpacingBetweenSlices in '\
              + f'SharedFunctionalGroupsSequence, {SegSBS}, matches the slice'\
              + ' thickness calculated from DICOM ImagePositionPatient and '\
              + f'ImageOrientationPatient, {Spacings[2]}.\n'
    else:
        if AreItemsEqualToWithinEpsilon(float(SegSBS), Spacings[2], epsilon):
            msg = 'INFO:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence, {SegSBS}, is within '\
                  + f'{epsilon} of the slice thickness calculated from DICOM '\
                  + 'ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
        else:
            #print(f'\n{float(SegSBS) - Spacings[2]}')
            msg = 'ERROR:  SEG SpacingBetweenSlices in '\
                  + f'SharedFunctionalGroupsSequence, {SegSBS}, is not within'\
                  + f' {epsilon} of the slice thickness calculated from DICOM'\
                  + ' ImagePositionPatient and ImageOrientationPatient, '\
                  + f'{Spacings[2]}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    
    SegPS = deepcopy(Seg.SharedFunctionalGroupsSequence[0]\
                        .PixelMeasuresSequence[0]\
                        .PixelSpacing)
    
    SegPS = [float(item) for item in SegPS]
    
    DcmPS = [float(item) for item in Dicom.PixelSpacing]
    
    if SegPS == DcmPS:
        msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence, '\
              + f'{SegPS}, matches the DICOM PixelSpacing, {DcmPS}.\n'
    else:
        if AreListsEqualToWithinEpsilon(SegPS, DcmPS, epsilon):
            msg = 'INFO:  SEG PixelSpacing in SharedFunctionalGroupsSequence,'\
                  + f' {SegPS} is within {epsilon} of the DICOM PixelSpacing,'\
                  + f' {DcmPS}.\n'
        else:
            msg = 'ERROR:  SEG PixelSpacing in SharedFunctionalGroupsSequence'\
                  + f', {SegPS}, is not within {epsilon} of the DICOM '\
                  + f'PixelSpacing, {DcmPS}.\n'
            Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Verify that the number of sequences in PerFrameFunctionalGroupsSequence
    is equal to NumberOfFrames."""
    PFFGS = deepcopy(Seg.PerFrameFunctionalGroupsSequence)
    
    if len(PFFGS) == int(Seg.NumberOfFrames):
        msg = f'INFO:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, '\
              + f'matches NumberOfFrames {Seg.NumberOfFrames}.\n'
    else:
        msg = f'ERROR:  The number of sequences in '\
              + f'PerFrameFunctionGroupsSequence, {len(PFFGS)}, does not '\
              + f'match NumberOfFrames {Seg.NumberOfFrames}.\n'
        Nerrors = Nerrors + 1
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether any of the ReferencedSOPInstanceUIDs in
    PerFrameFunctionGroupsSequence do not match the SOPInstanceUIDs."""        
    RefSOPuidsInPFFGS = [PFFGS[i].DerivationImageSequence[0]\
                                 .SourceImageSequence[0]\
                                 .ReferencedSOPInstanceUID for i in range(len(PFFGS))]
    
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInPFFGS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching ReferencedSOPInstanceUID:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
    if Inds:        
        for i in range(len(Inds)):
            uid = RefSOPuidsInPFFGS[Inds[i]]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + f'PerFrameFunctionalGroupsSequence[{Inds[i]}] does not '\
                  + 'match any DICOM SOPInstanceUID.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in '\
              + 'PerFrameFunctionalGroupsSequence match the DICOM '\
              + 'SOPInstanceUIDs.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine whether the DimensionIndexValues are sensible integers,
    and if the indexed ReferencedSOPInstanceUID agrees with the 
    ReferencedSOPInstanceUID within the SourceImageSequence."""
    DIVs = [PFFGS[i].FrameContentSequence[0]\
                    .DimensionIndexValues for i in range(len(PFFGS))]
    
    # Determine whether any of the first elements in the DIVs are not 1  
    # (there should only be one segment, so the first element should always
    # be 1):
    FirstElements = [int(DIVs[i][0]) for i in range(len(DIVs))]
    
    #NumOfMatches = FirstElements.count('1')
    
    ## Find the indices of any non-matching elements:
    #Inds = [i for i, x in enumerate(FirstElements) if x != 1]
    # Find the indices of any elements that exceed the number of sequences in
    # SegmentSequence, NumSS:
    Inds = [i for i, x in enumerate(FirstElements) if x > NumSS]
    
    if Inds:        
        for i in range(len(Inds)):
            DIV = DIVs[Inds[i]]
            
            msg = f'ERROR:  The first element in DimensionalIndexValue[{Inds[i]}],'\
                  + f' {DIV}, exceeds the number of sequences in '\
                  + f'SegmentSequene, {NumSS}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The first element in each DimensionalIndexValue are '\
              + 'consistent with the number of sequences in SegmentSequence, '\
              + f'{NumSS}.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
    # Determine if any of the second elements in the DIVs exceed the number
    # of sequences in ReferencedInstanceSequence:
    SecondElements = [int(DIVs[i][1]) for i in range(len(DIVs))]
    
    # Find the indices of any elements that exceed len(RIS):
    Inds = [i for i, x in enumerate(SecondElements) if x > len(RIS)]

    if Inds:       
        for i in range(len(Inds)):
            DIV = DIVs[Inds[i]]
            
            msg = f'ERROR:  The second element in DimensionalIndexValue[{Inds[i]}],'\
                  + f' {DIV}, exceeds the number of sequences in '\
                  + f'ReferencedInstanceSequence, {len(RIS)}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The second elements in DimensionalIndexValue '\
                  + f'matches the number of sequences in '\
                  + f'ReferencedInstanceSequence.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    #print(f'\n\n\nlen(RefSOPuidsInPFFGS) = {RefSOPuidsInPFFGS}')
    #print(f'len(SecondElements) = {len(SecondElements)}')
    #print(f'DIVs = {DIVs}')
    
    
    """Determine if the ReferencedSOPInstanceUID in the 
    ReferencedInstanceSequence indexed by the second elements in the DIVs 
    match the ReferencedSOPInstanceUID in the SourceImageSequence."""
    n = 0
    for i in range(len(SecondElements)):
        RefSOPinPFFGS = RefSOPuidsInPFFGS[i]
        
        # Need to -1 since i is zero-indexed:
        ind = SecondElements[i] - 1
        
        #print(f'ind = {ind}')
        
        RefSOPinRIS = RefSOPuidsInRIS[ind]
        
        if RefSOPinPFFGS == RefSOPinRIS:
            msg = 'INFO:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence[{i}] matches '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + f'ReferencedInstanceSequence[{DIVs[i]}], as indexed by '\
                  + f'DimensionalIndexValue[{i}], {DIVs[i]}.\n'
        else:
            msg = 'ERROR:  ReferencedSOPInstanceUID referenced in '\
                  + f'SourceImageSequence[{i}] does not match '\
                  + 'ReferencedSOPInstanceUID referenced in '\
                  + f'ReferencedInstanceSequence[{DIVs[i]}], as indexed by '\
                  + f'DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
            n += 1
                  
    if n == 0:
        msg = 'INFO:  The ReferencedSOPInstanceUIDs referenced in '\
                  + f'SourceImageSequence match the '\
                  + 'ReferencedSOPInstanceUIDs referenced in '\
                  + f'ReferencedInstanceSequence, as indexed by '\
                  + f'DimensionalIndexValue.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine if ImagePositionPatient in PerFrameFunctionalGroupsSequence 
    matches ImagePositionPatient of the DICOM as indexed by the DIVs."""
    IPPsInPFFGS = [PFFGS[i].PlanePositionSequence[0]\
                           .ImagePositionPatient for i in range(len(PFFGS))]
    
    # Convert from strings to floats:
    IPPsInPFFGS = [[float(item) for item in IPP] for IPP in IPPsInPFFGS]
    
    c = 0
    for i in range(len(IPPsInPFFGS)):
        if RefSOPuidsInPFFGS[i] in SOPuids:
            ind = SOPuids.index(RefSOPuidsInPFFGS[i])
            
            IPP = IPPs[ind]
            
            #if IPPsInPFFGS[i] != IPP:
            if AreListsEqualToWithinEpsilon(IPPsInPFFGS[i], IPP, epsilon):
                msg = f'INFO: ImagePositionPatient, {IPPsInPFFGS[i]}, '\
                      + f'referenced in PlanePositionSequence[{i}] is within '\
                      + f'{epsilon} of the DICOM ImagePositionPatient, {IPP},'\
                      + f' indexed by the second element, {SecondElements[i]}'\
                      + f', in DimensionalIndexValue[{i}], {DIVs[i]}.\n'
            else:
                msg = f'ERROR: ImagePositionPatient, {IPPsInPFFGS[i]}, '\
                      + f'referenced in PlanePositionSequence[{i+1}] is not '\
                      + f'within {epsilon} of the DICOM ImagePositionPatient,'\
                      + f' {IPP}, indexed by the second element, '\
                      + f'{SecondElements[i]}, in DimensionalIndexValue[{i}],'\
                      + f' {DIVs[i]}.\n'
                LogList.append(msg)
                Nerrors = Nerrors + 1
                c += 1
        
        if c == 0:
            msg = f'INFO: All ImagePositionPatient referenced in '\
                      + f'PlanePositionSequence are within {epsilon} of the '\
                      + f'DICOM ImagePositionPatient indexed by the second '\
                      + f'element in DimensionalIndexValue.\n'
            LogList.append(msg)
            
        #else:
        #    msg = f'.\n'
        #    
        #    Nerrors = Nerrors + 1
                  
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Check that the ReferencedSegmentNumber match the first elements in the
    DimensionIndexValues.
    """
    RSNs = [PFFGS[i].SegmentIdentificationSequence[0]\
                    .ReferencedSegmentNumber for i in range(len(PFFGS))]
    
    Diffs = [FirstElements[i] - RSNs[i] for i in range(len(RSNs))]
    
    Inds = [i for i, x in enumerate(Diffs) if x > 0]
    
    if Inds:       
        for i in range(len(Inds)):
            msg = f'ERROR:  ReferencedSegmentNumber in '\
                  + 'SegmentIdentificationSequence of '\
                  + f'PerFrameFunctionalGroupsSequence[{Inds[i]}], '\
                  + f'{RSNs[Inds[i]]}, does not match the first element in '\
                  + f'DimensionIndexValues[{Inds[i]}], {DIVs[Inds[i]]}.\n'
            LogList.append(msg)
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  Each ReferencedSegmentNumber in '\
              + 'SegmentIdentificationSequence of '\
              + 'PerFrameFunctionalGroupsSequence matches the first element '\
              + 'in DimensionIndexValues.\n'
        LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Determine if the dimensions of the pixel array in PixelData do not
    match the expected dimensions (the number of frames (PixArrF) has 
    already been compared to NumberOfFrames)."""
    if NumOfRows == Seg.Rows:
        msg = f'INFO:  The number of rows in PixelData, {NumOfRows}, '\
              + f'matches the number of SEG Rows, {Seg.Rows}.\n'
        LogList.append(msg)
    else:
        msg = f'ERROR:  The number of rows in PixelData, {NumOfRows}, does '\
              + f'not match the number of SEG Rows, {Seg.Rows}.\n'
        LogList.append(msg)
        Nerrors = Nerrors + 1
    
    
    if LogToConsole:
        print(msg)
            
    if NumOfCols == Seg.Columns:
        msg = f'INFO:  The number of columns in PixelData, {NumOfCols}, '\
              + f'matches the number of SEG Columns, {Seg.Columns}.\n'
        LogList.append(msg)
    else:
        msg = f'ERROR:  The number of columns in PixelData, {NumOfCols}, does'\
              + f' not match the number of SEG Columns, {Seg.Columns}.\n'
        LogList.append(msg)
        Nerrors = Nerrors + 1
    
    if LogToConsole:
        print(msg)
    
    
    print(f'There were {Nerrors} errors found in the SEG.\n')
    
    if LogToConsole:
        print('-'*120)
            
    return LogList, Nerrors







def GetSegDataFromListOfSegs(ListOfSegs, ListOfDicomDirs, SearchString, 
                             LogToConsole=False):
    """ 
    Note:
    
    ListOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    
    ListOfPlotTitles is a list (which could be of length 1) of text for each
    plot title.
    """
    
    import importlib
    import DicomTools
    importlib.reload(DicomTools)
    
    #import numpy as np
    from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    #from DicomTools import GetRoiNums
    from ImageTools import GetImageAttributes
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of GetSegDataFromListOfSegs():')
    
    ListOfPixArrBySeg = []
    ListOfF2SindsBySeg = []
    #ListOfSegNums = []
    ListOfDcmFpaths = []
    ListOfIPPs = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfSegs)):
        if LogToConsole:
            print(f'\n\nListOfSegs[{i}]:')
        
        if ListOfSegs[i]:
            PixArrBySeg, F2SindsBySeg = GetPixArrBySeg(ListOfSegs[i], 
                                                       ListOfDicomDirs[i], 
                                                       LogToConsole)
            
            if LogToConsole:      
                print(f'F2SindsBySeg = {F2SindsBySeg}')
            
        else:
            PixArrBySeg = None
            F2SindsBySeg = None
            
            if LogToConsole:      
                print(f'No segmentations for this dataset')
        
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
        ListOfPixArrBySeg.append(PixArrBySeg)
        ListOfF2SindsBySeg.append(F2SindsBySeg)
        ListOfDcmFpaths.append(DcmFpaths)
        #ListOfOrigins.append(IPPs[0])
        ListOfIPPs.append(IPPs)
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
            
        if LogToConsole:      
            print(f'\nlen(ListOfF2SindsBySeg) = {len(ListOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPixArrBySeg, ListOfF2SindsBySeg, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings