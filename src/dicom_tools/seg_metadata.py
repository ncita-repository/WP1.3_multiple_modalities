# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 17:27:27 2021

@author: ctorti
"""

def get_frameNums(f2sIndsBySeg, segNum):
    """
    Get the frame number(s) in a SEG's pixel array that matches the segment 
    label provided.
    
    Parameters
    ----------
    f2sIndsBySeg : list of a list of ints
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
    f2sInds : list of ints
        List of the slice numbers that correspond to each Per-frame Functional 
        Groups Sequence in the segment of interest. 
    """
    
    # The f2sInds for the segment of interest:
    f2sInds = f2sIndsBySeg[segNum]
    
    # The number of segments:
    S = len(f2sIndsBySeg)
    
    # The number of frames in the segment of interest:
    F = len(f2sInds)
    
    n = 0 # initialise frame number counter
    
    frameNums = []
    
    for s in range(S):
        if s == segNum:
            """ 
            This segment contains the frame(s) of interest, which are given by
            the index of f2sInds[f] in f2sIndsBySeg[s] plus 
            any frames that preceeded it (i.e. the frame counter n).
            """
            
            for f in range(F):
                frameNum = n + f2sIndsBySeg[s].index(f2sInds[f])
                
                frameNums.append(frameNum)
            
        else:
            # Add to the frame counter the number of frames in this segment:
            n += len(f2sIndsBySeg[s])
    
    return frameNums, f2sInds

def get_RSOPuids_in_RIS(seg):
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
    rsopuids = get_RSOPuids_in_RIS(seg)
    
    r2sInds = [] 
    
    for uid in rsopuids:
        # Find the matching index of uid in sopuids:
        r2sInds.append(sopuids.index(uid))
        
    return r2sInds

def get_RSOPuids_in_PFFGS(seg):
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

def get_f2sInds(seg, sopuids):
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
    f2sInds : list of ints
        List (for each frame) of the slice numbers that correspond to each 
        Per-FrameFunctionalGroupsSequence in seg.
    """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    rsopuids = get_RSOPuids_in_PFFGS(seg)
    
    f2sInds = [] 
    
    for i in range(len(rsopuids)):
        if rsopuids[i] in sopuids:
            # Find the matching index of rsopuids[i] in sopuids:
            f2sInds.append(sopuids.index(rsopuids[i]))
        else:
            msg = f'ReferencedSOPInstanceUID[{i}], {rsopuids[i]}, is not in '\
                  + 'the list of SOPInstanceUIDs.'
            
            raise Exception(msg)
        
    return f2sInds

def get_f2sIndsBySeg(seg, sopuids):
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
    f2sIndsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of the slice numbers 
        that correspond to each Per-FrameFunctionalGroupsSequence in seg.
    """
    
    f2sInds = get_f2sInds(seg, sopuids)
    
    divs = get_DIVs(seg)
    
    f2sIndsBySeg = group_list_by_seg(listToGroup=f2sInds, 
                                           divs=divs)
    
    return f2sIndsBySeg

def get_DIVs(seg):
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