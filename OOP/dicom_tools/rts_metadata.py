# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 16:21:16 2021

@author: ctorti
"""

#from copy import deepcopy
from conversion_tools.inds_pts_cntdata import (
    cntdata_to_pts, ptsByCntByRoi_to_cntdataByCntByRoi
    )
from general_tools.console_printing import (
    print_indsByRoi, print_ptsByCntByRoi
)
    

def get_RSOPuidsByRoi(rts):
    """
    Get the list of ReferencedSOPInstance UIDs for each ContourSequence
    for each ROI in a RT-STRUCT.
    
    Parameters
    ----------
    rts : Pydicom object
        Pydicom representation of an RTSTRUCT file.
    
    Returns
    -------
    RSOPuidsByRoi : list of list of strs
        List (for each ROI) of a list (for each ContourSequence) of 
        ReferencedSOPInstanceUIDs.
    
    Notes
    -----
    ContourSequence = rts.ROIContourSequence[0]\
                         .ContourSequence
    
    and
    
    ReferencedSOPInstanceUID = rts.ROIContourSequence[0]\
                                  .ContourSequence[i]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourSequence.
    
    If there is only one ROI Contour Sequence RSOPuids will be a list of
    length 1 of a list of strings, e.g.
        
    [['uid1', 'uid2', ...]]
    """
    
    RSOPuidsByRoi = []
    
    roiCntSeqs = rts.ROIContourSequence
    
    for r in range(len(roiCntSeqs)):
        cntSeqs = roiCntSeqs[r].ContourSequence
        
        # Initialise the list of ReferencedSOPInstanceUIDs for this ROI:
        RSOPuids = []
        
        for seq in cntSeqs:
            uid = seq.ContourImageSequence[0].ReferencedSOPInstanceUID
            
            RSOPuids.append(uid)
        
        RSOPuidsByRoi.append(RSOPuids)
        
    return RSOPuidsByRoi

def get_c2sIndsByRoi(rts, sopuids):
    """
    Get the slice numbers that correspond to each ContourSequence for all ROIs
    in an RTSTRUCT.
    
    Parameters
    ----------
    rts : Pydicom object
        Pydicom representation of an RTSTRUCT file.
    sopuids : list of strs
        List of SOPInstanceUIDs of the DICOM series that the RTSTRUCT relates
        to.
    
    Returns
    -------
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each sequence in ContourSequence.
    c2sInds : list of ints
        List (for each contour) of slice numbers that correspond to each 
        sequence in ContourSequence (i.e. c2sIndsByRoi flattened).
    """
    
    # Get the ReferencedSOPInstanceUIDs from the ContourSequence by ROI:
    RSOPuidsByRoi = get_RSOPuidsByRoi(rts)
    
    c2sInds = []
    c2sIndsByRoi = []
    
    # Loop through each list of ReferencedSOPInstanceUIDs:
    for i in range(len(RSOPuidsByRoi)):
        inds = []
    
        for RefUid in RSOPuidsByRoi[i]:
            # Find the matching index of RefUid in sopuids:
            inds.append(sopuids.index(RefUid))
            
        c2sInds.extend(inds)
        c2sIndsByRoi.append(inds)
        
    return c2sIndsByRoi, c2sInds

def get_ptsByCntByRoi(rts, sopuids, p2c=False):
    """  
    Get a list of the physical points grouped by contour grouped by ROI in an 
    RTS.
    
    Parameters
    ----------
    rts : Pydicom object
        Pydicom representation of an RTSTRUCT file.
    sopuids : list of strs
        List of SOPInstanceUIDs of the DICOM series that the RTSTRUCT relates
        to.
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
    
    Returns
    -------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    cntdataByCntByRoi : list of a list of a list of strs
        List (for each ROI) of a list (for all contours) of a flat list of 
        coordinates in ptsByCntByRoi converted from floats to strings.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour.
    c2sInds : list of ints
        List (for each contour) of slice numbers that correspond to each 
        contour.
    """
    
    #from copy import deepcopy
    #from DicomTools import GetDicomSOPuids
    #from ConversionTools import ContourData2Points
    #from GeneralTools import PrintIndsByRoi
    
    # Get the ContourSequence-to-slice indices by ROI:
    c2sIndsByRoi, c2sInds = get_c2sIndsByRoi(rts, sopuids)
    
    numRois = len(c2sIndsByRoi)
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        #print(f'   c2sIndsByRoi = {c2sIndsByRoi}')
        print_indsByRoi(c2sIndsByRoi)
        print(f'numRois = {numRois}')
    
    ptsByCntByRoi = []
    
    for r in range(numRois):
        # Number of contours in this ROI:
        numCnts = len(c2sIndsByRoi[r])
        
        # Get a list of the ContourSequences in this ROI:
        #sequences = deepcopy(rts.ROIContourSequence[r].ContourSequence)
        sequences = list(rts.ROIContourSequence[r].ContourSequence)
        
        # Get a list of all ReferencedSOPInstanceUIDs from ContourSequences:
        #RefSopUids = [sequence.ContourImageSequence[0]\
        #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        if False:#p2c:
            print(f'      r = {r}')
            print(f'      numCnts = {numCnts}')
            print(f'      number of Contour Sequences = {len(sequences)}')
        
        ptsByCnt = []
        
        # Loop through each contour sequence:
        for c in range(len(sequences)):
            # Get the indices of all matching RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first 
            # match, and there may be more than one!
            #inds = [i for i, e in enumerate(RefSopUids) if e==SopUid]
            
            # Iterate for each index in inds:
            #for ind in inds:
            
            #ContourSequence = deepcopy(ContourSequences[c])
            
            #ContourData = [float(item) for item in ContourSequence.ContourData]
            
            #cntData = deepcopy(sequence.ContourData)
            cntData = list(sequences[c].ContourData)
            
            pts = cntdata_to_pts(cntData)
            
            ptsByCnt.append(pts)
            
            if False:#p2c:
                print(f'         c = {c}')
                print(f'         len(pts) = {len(pts)}')
        
        ptsByCntByRoi.append(ptsByCnt)
    
    cntdataByCntByRoi = ptsByCntByRoi_to_cntdataByCntByRoi(
        ptsByCntByRoi
        )
    
    if p2c:
        #R = len(ptsByCntByRoi)
        #print(f'\n   len(ptsByCntByRoi) = {R}')
        #for r in range(R):
        #    print(f'      len(ptsByCntByRoi[{r}]) = {len(ptsByCntByRoi[r])}')
        print_ptsByCntByRoi(ptsByCntByRoi)
        print('-'*120)

    return ptsByCntByRoi, cntdataByCntByRoi, c2sIndsByRoi, c2sInds
