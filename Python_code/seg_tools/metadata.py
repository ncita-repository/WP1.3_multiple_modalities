# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:27:27 2021

@author: ctorti
"""

def get_frame_nums(PFFGStoSlcIndsBySeg, segNum):
    # seg, searchStr, dicomDir
    # this was modified on 16/07
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Parameters
    ----------
    PFFGStoSlcIndsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-frame Functional Groups Sequence.
    segNum : int
        The segment number of interest.
    
    Returns
    -------
    frameNums : list of ints
        List of the frame numbers that correspond to all frames in the SEG's 
        pixel array that belongs to the segment of interest.  If only one frame 
        corresponds to the segment, the returned list will be of length one 
        (e.g. [3]).
    PFFGStoSlcInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    # The PFFGStoSlcInds for the segment of interest:
    PFFGStoSlcInds = PFFGStoSlcIndsBySeg[segNum]
    
    # The number of segments:
    S = len(PFFGStoSlcIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSlcInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of PFFGStoSlcInds[f] in PFFGStoSlcIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + PFFGStoSlcIndsBySeg[s].index(PFFGStoSlcInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSlcIndsBySeg[s])
    
    return frameNums, PFFGStoSlcInds

def get_RSOP_uids_in_RIS(seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence
    of a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    
    Returns
    -------
    RSOPUIDs : list of strs
        List of ReferencedSOPInstanceUIDs.
    """
    
    RSOPUIDs = [] 
    
    sequences = seg.ReferencedSeriesSequence[0]\
                   .ReferencedInstanceSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        RSOPUIDs.append(uid)
        
    return RSOPUIDs

def get_RIS_to_sliceinds(seg, SOPUIDs):
    """
    Get the slice numbers that correspond to each ReferencedInstanceSequence
    in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    SOPUIDs : list of strs
        List of SOPInstanceUIDs of the DICOMs.
        
    Returns
    -------
    RIStoSlcInds : list of ints
        Slice numbers that correspond to each ReferencedInstanceSequence in 
        seg.
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPUIDs = get_RSOP_uids_in_RIS(seg)
    
    RIStoSlcInds = [] 
    
    for Ruid in RSOPUIDs:
        # Find the matching index of Ruid in SOPUIDs:
        RIStoSlcInds.append(SOPUIDs.index(Ruid))
        
    return RIStoSlcInds

def get_RSOP_uids_in_PFFGS(seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the 
    Per-FrameFunctionalGroupsSequence of a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    
    Returns
    -------
    RSOPUIDs : list of strs
        List of ReferencedSOPUIDs.
    """
    
    RSOPUIDs = [] 
    
    sequences = seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPUIDs.append(uid)
        
    return RSOPUIDs

def get_PFFGS_to_sliceinds(seg, SOPUIDs):
    """
    Get a list of the slice numbers that correspond to each  
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    SOPUIDs : list of strs
        SOPInstanceUIDs of the DICOMs.
        
    Returns
    -------
    PFFGStoSlcInds : list of ints
        List (for each frame) of the slice numbers that correspond to each 
        Per-FrameFunctionalGroupsSequence in seg.
    """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    RSOPUIDs = get_RSOP_uids_in_PFFGS(seg)
    
    PFFGStoSlcInds = [] 
    
    for i in range(len(RSOPUIDs)):
        if RSOPUIDs[i] in SOPUIDs:
            # Find the matching index of RSOPUIDs[i] in SOPUIDs:
            PFFGStoSlcInds.append(SOPUIDs.index(RSOPUIDs[i]))
        else:
            msg = f'ReferencedSOPInstanceUID[{i}], {RSOPUIDs[i]}, is not in '\
                  + 'the list of SOPInstanceUIDs.'
            
            raise Exception(msg)
        
    return PFFGStoSlcInds

def get_PFFGS_to_sliceindsBySeg(seg, SOPUIDs):
    """
    Get a list (per segment) of the slice numbers that correspond to each 
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    SOPUIDs : list of strs
        SOPInstanceUIDs of the DICOMs.
    
    Returns
    -------
    PFFGStoSlcIndsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-FrameFunctionalGroupsSequence in seg.
    """
    
    PFFGStoSlcInds = get_PFFGS_to_sliceinds(seg, SOPUIDs)
    
    DIVs = get_DIVs(seg)
    
    PFFGStoSlcIndsBySeg = group_list_by_seg(listToGroup=PFFGStoSlcInds, 
                                            DIVs=DIVs)
    
    return PFFGStoSlcIndsBySeg

def get_DIVs(seg):
    """
    Get the DimensionIndexValues in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    DIVs : list of lists of ints
        List of DimensionIndexValues in seg.
    """
    
    # Initialise the list of DIVs:
    DIVs = [] 
    
    sequences = seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        DIV = sequence.FrameContentSequence[0].DimensionIndexValues
        
        DIVs.append(DIV)
        
    return DIVs

def group_list_by_seg(listToGroup, DIVs):
    """
    Group a list of items by common segment/ROI using the DimensionIndexValues.
    
    Parameters
    ----------
    listToGroup : list
        The list to be grouped.
        
    DIVs : list of lists of ints
        List of DimensionIndexValues.
        
    Returns
    -------
    groupedList : list of lists
        listToGroup grouped by ROI.
        
    Notes
    -----
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
    
    If there's only one ROI, listToGroup will still be nested in a list to 
    maintain continuity.  In the case of the DIVs with only one ROI:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ]
     ]
    """
    
    # Get the list of ROI numbers:
    roiNums = []
    
    for DIV in DIVs:
        roiNums.append(DIV[0])
        
    uniqueRoiNums = list(set(roiNums))
    
    if len(uniqueRoiNums) == 1:
        groupedList = [listToGroup]
    else:
        groupedList = []
        
        for roiNum in uniqueRoiNums:
            listThisRoi = [] # initialise
            
            for i in range(len(listToGroup)):
                if DIVs[i][0] == roiNum:
                    listThisRoi.append(listToGroup[i])
        
            groupedList.append(listThisRoi)
             
    return groupedList





""" Move the two functions below elsewhere: """

def AddToPFFGS(seg):
    """
    Append the last item in the Per-FrameFunctionalGroupsSequence of the
    SEG to itself (i.e. increase the length of the sequence by one). 
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    newSeg : Pydicom Object
        Modified ROI object with Per-FrameFunctionalGroupsSequence lengthened
        by one.
    """
    
    from copy import deepcopy
    
    # Use seg as a template for newSeg: 
    newSeg = deepcopy(seg)
             
    # The last item in Per-frame Functional Groups Sequence:
    last = deepcopy(seg.PerFrameFunctionalGroupsSequence[-1])
    
    # Append to the Per-frame Functional Groups Sequence:
    newSeg.PerFrameFunctionalGroupsSequence.append(last)
             
    return newSeg

def AddToSS(seg):
    """
    Append the last item in the SegmentSequence of the SEG Object to itself 
    (i.e. increase the length of the sequence by one). 
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    newSeg : Pydicom Object
        Modified seg with SegmentSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use seg as a template for newSeg: 
    newSeg = deepcopy(seg)
             
    # The last item in Segment Sequence:
    last = deepcopy(seg.SegmentSequence[-1])
    
    # Append to the Segment Sequence:
    newSeg.SegmentSequence.append(last)
             
    return newSeg



""" Old functions below """


def get_frame_nums_OLD(seg, searchStr, dicomDir):
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Parameters
    ----------
    seg : Pydicom Object
        SEG object.
    searchStr : str
        All or part of the Segment Label containing the segmentation of
        interest.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    frameNums : list of ints
        List of the frame numbers that correspond to all frames in the SEG's 
        pixel array that belongs to the segment of interest.  If only one frame 
        corresponds to the segment, the returned list will be of length one 
        (e.g. [3]).
    PFFGStoSlcInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    from dicom_tools.metadata import get_dcm_uids, get_roicol_nums
    
    # Get the DimensionIndexValues (DIVs) and group the DIVs by segment.
    """ Note:  The DIVs are integers that start from 1. """
    DIVs = get_DIVs(seg)
    
    # Get the DICOM SOP UIDs:
    studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(dicomDir)
    
    # Get the PerFrameFunctionalGroupsSequence-to-slice indices
    # (PFFGStoSlcInds) and PFFGStoSlcInds grouped by segment
    # (PFFGStoSlcIndsBySeg):
    """ Note:  The PFFGStoSlcInds are integers that start from 0. """
    PFFGStoSlcInds = get_PFFGS_to_sliceinds(seg, SOPUIDs)
    
    PFFGStoSlcIndsBySeg = group_list_by_seg(listToGroup=PFFGStoSlcInds, 
                                              DIVs=DIVs)
    
    segNum = get_roicol_nums(seg, searchStr)
    
    # The PFFGStoSlcInds for the segment of interest:
    PFFGStoSlcInds = PFFGStoSlcIndsBySeg[segNum]
    
    # The number of segments:
    S = len(PFFGStoSlcIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(PFFGStoSlcInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of PFFGStoSlcInds[f] in PFFGStoSlcIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + PFFGStoSlcIndsBySeg[s].index(PFFGStoSlcInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(PFFGStoSlcIndsBySeg[s])
    
    return frameNums, PFFGStoSlcInds