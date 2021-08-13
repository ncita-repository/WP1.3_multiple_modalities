# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:27:27 2021

@author: ctorti
"""

def get_frameNums(p2sIndsBySeg, segNum):
    # seg, searchStr, dicomDir
    # this was modified on 16/07
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Parameters
    ----------
    p2sIndsBySeg : list of a list of ints
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
    p2sInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    # The p2sInds for the segment of interest:
    p2sInds = p2sIndsBySeg[segNum]
    
    # The number of segments:
    S = len(p2sIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(p2sInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of p2sInds[f] in p2sIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + p2sIndsBySeg[s].index(p2sInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(p2sIndsBySeg[s])
    
    return frameNums, p2sInds

def get_rsopuids_in_ris(seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the ReferencedInstanceSequence
    of a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    
    Returns
    -------
    rsopuids : list of strs
        List of ReferencedSOPInstanceUIDs.
    """
    
    rsopuids = [] 
    
    sequences = seg.ReferencedSeriesSequence[0]\
                   .ReferencedInstanceSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        rsopuids.append(uid)
        
    return rsopuids

def get_r2sInds(seg, sopuids):
    """
    Get the slice numbers that correspond to each ReferencedInstanceSequence
    in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    sopuids : list of strs
        List of SOPInstanceUIDs of the DICOMs.
        
    Returns
    -------
    r2sInds : list of ints
        Slice numbers that correspond to each ReferencedInstanceSequence in 
        seg.
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    rsopuids = get_rsopuids_in_ris(seg)
    
    r2sInds = [] 
    
    for uid in rsopuids:
        # Find the matching index of uid in sopuids:
        r2sInds.append(sopuids.index(uid))
        
    return r2sInds

def get_rsopuids_in_pffgs(seg):
    """
    Get the list of ReferencedSOPInstanceUIDs in the 
    Per-FrameFunctionalGroupsSequence of a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    
    Returns
    -------
    rsopuids : list of strs
        List of Referencedsopuids.
    """
    
    rsopuids = [] 
    
    sequences = seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        rsopuids.append(uid)
        
    return rsopuids

def get_p2sInds(seg, sopuids):
    """
    Get a list of the slice numbers that correspond to each  
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    sopuids : list of strs
        SOPInstanceUIDs of the DICOMs.
        
    Returns
    -------
    p2sInds : list of ints
        List (for each frame) of the slice numbers that correspond to each 
        Per-FrameFunctionalGroupsSequence in seg.
    """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    rsopuids = get_rsopuids_in_pffgs(seg)
    
    p2sInds = [] 
    
    for i in range(len(rsopuids)):
        if rsopuids[i] in sopuids:
            # Find the matching index of rsopuids[i] in sopuids:
            p2sInds.append(sopuids.index(rsopuids[i]))
        else:
            msg = f'ReferencedSOPInstanceUID[{i}], {rsopuids[i]}, is not in '\
                  + 'the list of SOPInstanceUIDs.'
            
            raise Exception(msg)
        
    return p2sInds

def get_p2sIndsBySeg(seg, sopuids):
    """
    Get a list (per segment) of the slice numbers that correspond to each 
    Per-FrameFunctionalGroupsSequence in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
    sopuids : list of strs
        SOPInstanceUIDs of the DICOMs.
    
    Returns
    -------
    p2sIndsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-FrameFunctionalGroupsSequence in seg.
    """
    
    p2sInds = get_p2sInds(seg, sopuids)
    
    divs = get_divs(seg)
    
    p2sIndsBySeg = group_list_by_seg(listToGroup=p2sInds, 
                                           divs=divs)
    
    return p2sIndsBySeg

def get_divs(seg):
    """
    Get the DimensionIndexValues in a SEG.
    
    Parameters
    ----------
    seg : Pydicom Object
        ROI object from a SEG file.
        
    Returns
    -------
    divs : list of lists of ints
        List of DimensionIndexValues in seg.
    """
    
    # Initialise the list of divs:
    divs = [] 
    
    sequences = seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        div = sequence.FrameContentSequence[0].DimensionIndexValues
        
        divs.append(div)
        
    return divs

def group_list_by_seg(listToGroup, divs):
    """
    Group a list of items by common segment/ROI using the DimensionIndexValues.
    
    Parameters
    ----------
    listToGroup : list
        The list to be grouped.
        
    divs : list of lists of ints
        List of DimensionIndexValues.
        
    Returns
    -------
    groupedList : list of lists
        listToGroup grouped by ROI.
        
    Notes
    -----
    divs is a list of lists of the form:
            
    [ [1, 10], [1, 11], [1, 12], ... [2, 6], [2, 7], [2, 8], ... ]
    
    where each list of index values that belong to a common ROI are separated
    into separate lists.
    
    The goal is to group a list (of anything) by ROI.  In the case of the divs,
    that is:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ],
      [ [2, 6], [2, 7], [2, 8], ... ],
      ... 
     ]
    
    If there's only one ROI, listToGroup will still be nested in a list to 
    maintain continuity.  In the case of the divs with only one ROI:
        
    [ [ [1, 10], [1, 11], [1, 12], ... ]
     ]
    """
    
    # Get the list of ROI numbers:
    roiNums = []
    
    for div in divs:
        roiNums.append(div[0])
        
    uniqueRoiNums = list(set(roiNums))
    
    if len(uniqueRoiNums) == 1:
        groupedList = [listToGroup]
    else:
        groupedList = []
        
        for roiNum in uniqueRoiNums:
            listThisRoi = [] # initialise
            
            for i in range(len(listToGroup)):
                if divs[i][0] == roiNum:
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
    p2sInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    from dicom_tools.metadata import get_dcm_uids, get_roicol_nums
    
    # Get the DimensionIndexValues (divs) and group the divs by segment.
    """ Note:  The divs are integers that start from 1. """
    divs = get_divs(seg)
    
    # Get the DICOM SOP UIDs:
    studyUID, seriesUID, FORUID, sopuids = get_dcm_uids(dicomDir)
    
    # Get the PerFrameFunctionalGroupsSequence-to-slice indices
    # (p2sInds) and p2sInds grouped by segment
    # (p2sIndsBySeg):
    """ Note:  The p2sInds are integers that start from 0. """
    p2sInds = get_p2sInds(seg, sopuids)
    
    p2sIndsBySeg = group_list_by_seg(listToGroup=p2sInds, 
                                           divs=divs)
    
    segNum = get_roicol_nums(seg, searchStr)
    
    # The p2sInds for the segment of interest:
    p2sInds = p2sIndsBySeg[segNum]
    
    # The number of segments:
    S = len(p2sIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(p2sInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of p2sInds[f] in p2sIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + p2sIndsBySeg[s].index(p2sInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(p2sIndsBySeg[s])
    
    return frameNums, p2sInds