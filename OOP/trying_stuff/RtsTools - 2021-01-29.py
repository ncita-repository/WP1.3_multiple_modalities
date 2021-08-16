# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:55:17 2020

@author: ctorti
"""




# Import packages and functions:
#import os
#import time
#from copy import deepcopy
#from pydicom import dcmread
#from pydicom.uid import generate_uid
#from DicomTools import GetRoiLabels
#from DicomTools import GetDicomSOPuids
#from ImageTools import GetImageAttributes
#from ConversionTools import Point2Index
#from ConversionTools import Indices2Mask


"""
******************************************************************************
******************************************************************************
DICOM RTSTRUCT FUNCTIONS
******************************************************************************
******************************************************************************
"""


def GetContourNumByRoi(Rts):
    """
    COMMENT 20/01/21: 
        This function and the description below is misleading.  What it
        actually outputs is a list (for each ROI) of a list (for each contour)
        of the ROI number.  It comes as little surprise that this function 
        isn't used anywhere.
    
    
    Get the list of the Referenced ROI number for each ContourSequence in a
    RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    ConNumsByRoi - List of list of integers
        The Referenced ROI numbers for all Contour Sequences grouped by ROI
        number.  The integers are from 0, unlike ROI Number which begin from 1. 
        
        
    Note:
    ----
        If there is only one ROI Contour Sequence, ContourNumByRoi will be a 
        list of length 1, e.g.
        
        [[0, 0, ...]]
    """
    
    ConNumsByRoi = []
    
    Rois = Rts.ROIContourSequence
    
    for i in range(len(Rois)):
        # The number of contours in this sequence:
        n = len(Rois[i].ContourSequence)
        
        ConNumsByRoi.append([i]*n)
        
    return ConNumsByRoi






def GetNumOfContoursByRoi(Rts):
    """
    Get the number of contours in each ContourSequence in a RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    NumOfConsByRoi - List of integers
        The number of contours in each ContourSequence (counting from 0). 
        
    Note:
        If there is only one ROI Contour Sequence, NumOfConsByRoi will be a 
        list of length 1, e.g. [3].
    """
    
    NumOfConsByRoi = []
    
    Rois = Rts.ROIContourSequence
    
    for i in range(len(Rois)):
        # The number of contours in this sequence:
        n = len(Rois[i].ContourSequence)
        
        NumOfConsByRoi.append(n)
        
    return NumOfConsByRoi







def GroupListByRoi(ListToGroup, NumOfConsByRoi):
    """
    Group a list of items by common contour ROI.
    
    Inputs:
    ******
    
    ListToGroup : list
        List to be grouped.
        
    NumOfConsByRoi : list integers 
        The number of contours in each ContourSequence (counting from 0).
        
        
    Outputs:
    *******
    
    GroupedList : list of lists 
        The items in ListToGroup grouped by ROI.
    """
    
    GroupedList = []
    
    n = 0
    
    for roi in range(len(NumOfConsByRoi)):
        # The number of contours in this ROI:
        C = NumOfConsByRoi[roi]
        
        ListThisRoi = []
        
        for c in range(C):
            ListThisRoi.append(ListToGroup[n])
            
            n += 1
            
        GroupedList.append(ListThisRoi)
        
    return GroupedList








def GetRSOPuidsInCIS(Rts):
    """
    Get the list of referenced SOP Instance UIDs in the ContourImageSequence
    of a RTS ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
    
        
    Outputs:
    *******
    
    RSOPuids : list of strings
        List of Referenced SOP UIDs.
    """
    
    RSOPuids = [] 
    
    sequences = Rts.ReferencedFrameOfReferenceSequence[0]\
                   .RTReferencedStudySequence[0]\
                   .RTReferencedSeriesSequence[0]\
                   .ContourImageSequence
    
    for sequence in sequences:
        uid = sequence.ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids






def GetCIStoSliceInds(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each ContourImageSequence in a
    RTSTRUCT ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    CIStoSliceInds : list of integers
        List of slice numbers that correspond to each sequence in
        ContourImageSequence.
    """
    
    # Get the ReferencedSOPInstanceUIDs from ContourImageSequence:
    RSOPuids = GetRSOPuidsInCIS(Rts)
    
    CIStoSliceInds = [] 
    
    for RefUid in RSOPuids:
        # Find the matching index of RefUid in SOPuids:
        CIStoSliceInds.append(SOPuids.index(RefUid))
        
    return CIStoSliceInds






def GetRSOPuidsInCSByRoi(Rts):
    """
    Get the list of ReferencedSOPInstance UIDs for each ContourSequence
    of a RTS ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    RSOPuidsByRoi : list of list of strings
        List (for each ROI) of a list (for each contour) of 
        ReferencedSOPInstanceUIDs.
        
        
    Notes:
    -----
    
    If there is only one ROI Contour Sequence RSOPuids will be a list of
    length 1 of a list of strings, e.g.
        
    [['uid1', 'uid2', ...]]
    """
    
    RSOPuidsByRoi = []
    
    Rois = Rts.ROIContourSequence
    
    for roi in Rois:
        Sequences = roi.ContourSequence
        
        # Initialise the list of ReferencedSOPInstanceUIDs for this ROI:
        RSOPuids = []
        
        for sequence in Sequences:
            uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
            
            RSOPuids.append(uid)
        
        RSOPuidsByRoi.append(RSOPuids)
        
    return RSOPuidsByRoi






def GetCStoSliceIndsByRoi(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each ContourSequence for all ROIs
    in a RTSTRUCT ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    SOPuids : list of strings
        List of SOPInstanceUIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    CStoSliceIndsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each sequence in ContourSequence.
    """
    
    # Get the ReferencedSOPInstanceUIDs from the ContourSequence by ROI:
    RSOPuidsByRoi = GetRSOPuidsInCSByRoi(Rts)
    
    CStoSliceIndsByRoi = [] 
    
    # Loop through each list of ReferencedSOPInstanceUIDs:
    for i in range(len(RSOPuidsByRoi)):
        inds = []
    
        for RefUid in RSOPuidsByRoi[i]:
            # Find the matching index of RefUid in SOPuids:
            inds.append(SOPuids.index(RefUid))
            
        CStoSliceIndsByRoi.append(inds)
        
    return CStoSliceIndsByRoi









def AddToCIS(Rts):
    """
    Append the last item in the Contour Image Sequence.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    NewRts : Pydicom object
        Modified ROI object with ContourImageSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Rts as a template for NewRts: 
    NewRts = deepcopy(Rts)
        
    # The last item in ContourImageSequence:
    last = deepcopy(Rts.ReferencedFrameOfReferenceSequence[0]\
                       .RTReferencedStudySequence[0]\
                       .RTReferencedSeriesSequence[0]\
                       .ContourImageSequence[-1])
    
    # Append to ContourImageSequence:
    NewRts.ReferencedFrameOfReferenceSequence[0]\
          .RTReferencedStudySequence[0]\
          .RTReferencedSeriesSequence[0]\
          .ContourImageSequence\
          .append(last)
             
    return NewRts









def AddToCS(Rts):
    """
    Append the last item in the ContourSequence.
    
    Inputs:
    *******
        
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    NewRts : Pydicom object
        Modified ROI object with ContourSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Rts as a template for NewRts: 
    NewRts = deepcopy(Rts)
        
    # The last item in ContourSequence:
    last = deepcopy(Rts.ROIContourSequence[0]\
                       .ContourSequence[-1])
    
    # Append to ContourSequence:
    NewRts.ROIContourSequence[0]\
          .ContourSequence\
          .append(last)
             
    return NewRts









def GetPtsByCntByRoi(Rts, DicomDir, LogToConsole=False):
    """
    Note 20/01/21:
        The name of this function was previously GetPtsByRoi.
        
        
    Get a list of the physical points grouped by contour grouped by ROI in an 
    RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each list in PtsByRoi. 
    """
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from ConversionTools import ContourData2Points
    
    # Get the DICOM SOPInstanceUIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourImageSequence-to-slice indices:
    #CIStoSliceInds = GetCIStoSliceInds(Rts, SopUids)
    
    # Get the ContourSequence-to-slice indices by ROI:
    C2SindsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    Nrois = len(C2SindsByRoi)
    
    if LogToConsole:
        print(f'\nC2SindsByRoi = {C2SindsByRoi}')
        print(f'Nrois = {Nrois}')
    
    PtsByCntByRoi = []
    
    for r in range(Nrois):
        # Number of contours in this ROI:
        NumC2Sinds = len(C2SindsByRoi[r])
        
        # Get a list of the ContourSequences in this ROI:
        CS = deepcopy(Rts.ROIContourSequence[r].ContourSequence)
        
        NumCS = len(CS)
        
        # Get a list of all ReferencedSOPInstanceUIDs from ContourSequences:
        #RefSopUids = [sequence.ContourImageSequence[0]\
        #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        if LogToConsole:
            print(f'   r = {r}')
            print(f'   NumC2Sinds = {NumC2Sinds}')
            print(f'   NumCS = {NumCS}')
        
        
        PtsByCnt = []
        
        
        # Loop through each contour sequence:
        for c in range(NumCS):
            
            # Get the indices of all matching RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first 
            # match, and there may be more than one!
            #inds = [i for i, e in enumerate(RefSopUids) if e==SopUid]
            
            # Iterate for each index in inds:
            #for ind in inds:
            
            #ContourSequence = deepcopy(ContourSequences[c])
            
            #ContourData = [float(item) for item in ContourSequence.ContourData]
            
            CntData = deepcopy(CS[c].ContourData)
            
            Pts = ContourData2Points(CntData)
            
            PtsByCnt.append(Pts)
            
            if LogToConsole:
                print(f'      c = {c}')
                print(f'      len(Pts) = {len(Pts)}')
            
            
        PtsByCntByRoi.append(PtsByCnt)

    return PtsByCntByRoi, C2SindsByRoi











def GetPtsByCntForRoi(Rts, DicomDir, SearchString='', LogToConsole=False):
    """
    Note 20/01/21:
        The name of this function was previously GetPtsInRoi.
        
        
    Get a list of the physical points grouped by contour for an ROI that 
    matches a search string in an RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    SearchString : string (optional; '' by default)
        All or part of the ROI name of the ROI of interest.  If '', the first
        ROI will be used.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PtsByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates. 
        
    C2Sinds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        list of points in PointsByContour.
    """
    
    from DicomTools import GetRoiNum
    
    #PtsByRoi, CStoSliceIndsByRoi = GetPtsByRoi(Rts, DicomDir)
    PtsByCntByRoi, C2SindsByRoi = GetPtsByCntByRoi(Rts, DicomDir)
    
    if SearchString:
        RoiNum = GetRoiNum(Rts, SearchString)
    else:
        RoiNum = 0
    
    #PointsByContour = PtsByRoi[RoiNum]
    PtsByCnt = PtsByCntByRoi[RoiNum]
    
    C2Sinds = C2SindsByRoi[RoiNum]
    
    if LogToConsole:
        print('\n\nResults of GetPtsInRoi:')
        print(f'   len(PtsByCntByRoi) = {len(PtsByCntByRoi)}')
        print(f'   len(C2SindsByRoi) = {len(C2SindsByRoi)}')
        print(f'   len(PtsByCnt) = {len(PtsByCnt)}')
        print(f'   len(C2Sinds) = {len(C2Sinds)}')
        print('end of GetPtsByRoi.\n')
        
    return PtsByCnt, C2Sinds










def GetPtsByCntForSliceNum(Rts, SliceNum, DicomDir, SearchString='', 
                           LogToConsole=False):
    """
    NOTE 20/01/21:
        This is a revision of GetPtsInContour, which had a misleading name, 
        since the function returns the points grouped by contour for a slice of
        interest, i.e. the result may include multiple contours.
        
        This function is quite time consuming to run.
        
        
    Get a list of the physical points grouped by all contours that relate to a
    specified slice number and ROI that matches a search string in an RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
        
    SliceNum : integer (0-indexed)
        The slice number of interest.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    SearchString : string (optional; '' by default)
        All or part of the ROI name of the ROI containing the contour(s) to be
        copied.  If '', the first ROI will be used.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
                            
    Outputs:
    *******
    
    PtsByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates that define all contours for the slice of
        interest. 
    """
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    from ConversionTools import ContourData2Points
    
    RoiCS = deepcopy(Rts.ROIContourSequence)
    
    NumOfRois = len(RoiCS)
    
    if SearchString:
        RoiNum = GetRoiNum(Rts, SearchString)
        
        if RoiNum > NumOfRois - 1:
            msg = f'There are only {NumOfRois} ROIs in the RTS.\n'\
                  + f'RoiNum = {RoiNum} > {NumOfRois}.'
                  
            raise Exception(msg)
    else:
        RoiNum = 0
    
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourSequence-to-slice indices by ROI:
    C2SindsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    C2Sinds = C2SindsByRoi[RoiNum]
    
    #print(f'\nCStoSliceInds = {CStoSliceInds}')
    
    if not SliceNum in C2Sinds:
        raise Exception(f'There is no contour on slice {SliceNum} in the RTS.')
        
    # The contour number(s) of interest:
    #ContourNum = CStoSliceInds.index(SliceNum) # 20/01/21 this yields the 
    # first match only- ignores any others
    ContourNums = [i for i, e in enumerate(C2Sinds) if e==SliceNum]
    
    #print(f'\nContourNums = {ContourNums}')
    
    # The ContourSequences in the ROI of interest:
    #ContourSequence = deepcopy(Rts.ROIContourSequence[RoiNum].ContourSequence)
    CS = deepcopy(RoiCS[RoiNum].ContourSequence)
    
    # The number of contours in this ROI:
    NumCS = len(CS)
    
    #print(f'\nNumOfContours = {NumOfContours}')
    
    if len(ContourNums) > NumCS:
        msg = f'There are only {NumCS} contours in the RTS but there '\
              + f'are {len(ContourNums)} contour numbers of interest.'
                        
        raise Exception(msg)
        
    CntData = []
    
    for c in range(len(ContourNums)):
        CntData.append(CS[ContourNums[c]].ContourData)
        
    #print(f'\nContourData = {ContourData}')
    
    PtsByCnt = []
    
    for c in range(len(CntData)):
        PtsByCnt.append(ContourData2Points(CntData[c]))
    
    #CStoSliceInd = [SliceNum]
    
    if LogToConsole:
        print(f'\n\nResults of GetPtsByCntForSliceNum():')
        print(f'   NumOfRois = {NumOfRois}')
        print(f'   C2SindsByRoi = {C2SindsByRoi}')
        print(f'   RoiNum = {RoiNum}')
        print(f'   C2Sinds = C2SindsByRoi[RoiNum={RoiNum}] =',
              f'{C2Sinds}')
        print(f'   ContourNums = {ContourNums} (the contour numbers that lie',
              f'on SliceNum {SliceNum})')
        print(f'   len(CntData) = {len(CntData)} (= no. of contours',
              f'that exist for SliceNum {SliceNum})')
        Npts = []
        for c in range(len(CntData)):
            Npts.append(len(CntData[c])//3)
        print(f'   The no. of points in each contour in ContourData is {Npts}')
        print(f'   len(PtsByCnt) = {len(PtsByCnt)} (= no. of contours',
              f'that exist for SliceNum {SliceNum})')
        Npts = []
        for c in range(len(PtsByCnt)):
            Npts.append(len(PtsByCnt[c]))
        print(f'   The no. of points in each contour in PtsByCnt is {Npts}')
        #numofpts = []
        #for c in range(NumOfContours):
        #    contourdata = deepcopy(ContourSequence[c].ContourData)
        #
        #    pts = ContourData2Points(contourdata)
        #    
        #    numofpts.append(len(pts))
        #print(f'   Number of pts in each CS = {numofpts}')
        print('end of results.')
    
    return PtsByCnt








def GetPtsByCntBySliceNum(Rts, DicomDir, SearchString='', LogToConsole=False):
    """
    Get a list of the physical points grouped by all contours grouped by slice 
    number for all slices referenced in ContourImageSequence, and whose ROI 
    matches a search string in an RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    SearchString : string (optional; '' by default)
        All or part of the ROI name of the ROI containing the contour(s) to be
        copied.  If '', the first ROI will be used.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
                            
    Outputs:
    *******
    
    PtsByCntBySlice : list of a list of a list of floats
        List (for each slice) of a list (for each contour) of a list (for each 
        point) of a list (for each dimension) of coordinates that define all 
        contours for all slices that have contours. 
    """
    
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    #PtsByCntByRoi, C2SindsByRoi = GetPtsByCntByRoi(Rts, DicomDir, LogToConsole)
    
    # Get the DICOM SOPInstanceUIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourSequence-to-slice indices by ROI:
    C2SindsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    if SearchString:
        RoiNum = GetRoiNum(Rts, SearchString)
    else:
        RoiNum = 0
    
    C2Sinds = C2SindsByRoi[RoiNum]
    
    PtsByCntBySlice = []
    
    """Since any given slice number can be repeated in C2Sinds, rather than
    iterating through each C2Sinds, iterate over its (unique) set."""
    #for SliceNum in C2Sinds:
    for SliceNum in list(set(C2Sinds)):
        # Get the points-by-contour for this SliceNum:
        PtsByCnt = GetPtsByCntForSliceNum(Rts, SliceNum, DicomDir, 
                                          SearchString='', LogToConsole=False)
              
        PtsByCntBySlice.append(PtsByCnt)
    
    
    if LogToConsole:
        print(f'\n\nResults of GetPtsByCntBySliceNum():')
        print(f'   NumOfRois = {len(C2SindsByRoi)}')
        print(f'   C2SindsByRoi = {C2SindsByRoi}')
        print(f'   RoiNum = {RoiNum}')
        print(f'   C2Sinds = C2SindsByRoi[RoiNum={RoiNum}] =',
              f'{C2Sinds}')
        NcntsBySlice = []
        NptsByCntBySlice = []
        for s in range(len(PtsByCntBySlice)):
            NcntsBySlice.append(len(PtsByCntBySlice[s]))
            NptsByCnt = []
            for c in range(len(PtsByCntBySlice[s])):
                NptsByCnt.append(len(PtsByCntBySlice[s][c]))
            NptsByCntBySlice.append(NptsByCnt)
        print('   The no. of contours in each slice in PtsByCntBySlice is',
              f'{NcntsBySlice}')
        print('   The no. of points in each contour in each slice in',
              f'PtsByCntBySlice is {NptsByCntBySlice}')
        print('end of results.')
    
    return PtsByCntBySlice











def AddCopiedPts(OrigPtsByCnt, PtsToAdd, CStoSliceInds, SliceNum):
    """
    Add points in a contour to an existing list of points-by-contour.
       
    
    Inputs:
    ******
    
    OrigPtsByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates. 
        
    PtsToAdd : list of a list of floats
        List (for each point) of a list (for each dimension) of coordinates. 
        
    CStoSliceInds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        sequence in ContourSequence.
        
    SliceNum : integer
        The slice number (counting from 0) that corresponds to the contour of
        interest.
    
    
    Outputs:
    *******
    
    NewCntDataByCnt : list of a list of floats
        Modfied list (for each contour) of a flat list of [x, y, z] coordinates
        of the polygons that define each contour including PtsToAdd.
        
    NewPtsByCnt : list of a list of a list of floats
        A list (for each contour) of a list (for each point) of a list (for 
        each dimension) of the polygons that define each contour including 
        PtsToAdd.
    """
    
    from ConversionTools import Points2ContourData
    
    if SliceNum in CStoSliceInds:
        """
        PtsToAdd will overwrite existing points for slice SliceNum in 
        OrigPtsByCnt.
        """
        
        from copy import deepcopy
        
        NewPtsByCnt = deepcopy(OrigPtsByCnt)
        
        ind = CStoSliceInds.index(SliceNum)
        
        NewPtsByCnt[ind] = PtsToAdd
        
    else:
        """
        PtsToAdd will be added to OrigPtsByCnt.
        """
        
        from GeneralTools import AppendItemToListAndSort
        
        NewCStoSliceInds = AppendItemToListAndSort(OrigList=CStoSliceInds, 
                                                   ItemToAppend=SliceNum)
        
        NewPtsByCnt = []
        
        for i in range(len(NewCStoSliceInds)):
            s = NewCStoSliceInds[i]
            
            if s == SliceNum:
                NewPtsByCnt.append(PtsToAdd)
            else:
                ind = CStoSliceInds.index(s)
                
                NewPtsByCnt.append(OrigPtsByCnt[ind])
                
    
    # Create a list of flat ContourData for each list of contour points in 
    # NewPtsByCnt:
    NewCntDataByCnt = []
    
    for Points in NewPtsByCnt:
        NewCntDataByCnt.append(Points2ContourData(Points))
    
    
    return NewCntDataByCnt, NewPtsByCnt






    
    
def AddCopiedPtsByCnt(OrigPtsByCnt, OrigCStoSliceInds, 
                      PtsToAddByCnt, CStoSliceIndsToAddByCnt,
                      LogToConsole=False):
    """
    Add points in a contour to an existing list of points-by-contour.
       
    
    Inputs:
    ******
    
    OrigPtsByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates in the original list of points-by-contour. 
        
    OrigCStoSliceInds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        sequence in ContourSequence in the original list of points-by-contour.
        
    PtsToAddByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates in the list of points-by-contour to be added. 
        
    CStoSliceIndsToAddByCnt : list of integers
        List (for each contour) of slice numbers that correspond to each 
        contour in PtsToAddByCnt.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    NewCntDataByCnt : list of a list of floats
        Modfied list (for each contour) of a flat list of [x, y, z] coordinates
        of the polygons that define each contour including PtsToAdd.
        
    NewPtsByCnt : list of a list of a list of floats
        Modified list (for each contour) of a list (for each point) of a list 
        (for each dimension) of the polygons that define each contour including 
        PtsToAdd.
        
    NewCStoSliceInds : list of integers
        Modified list (for each contour) of slice numbers that correspond to  
        each sequence in ContourSequence in the original list of points-by-
        contour.
    """
    
    from GeneralTools import UniqueItems
    from ConversionTools import Points2ContourData
    
    if LogToConsole:
        print('\n\nRunning AddCopiedPtsByCnt:')
        print(f'   OrigCStoSliceInds = {OrigCStoSliceInds}')
        print(f'   CStoSliceIndsToAddByCnt = {CStoSliceIndsToAddByCnt}')
    
    # Combine OrigCStoSliceInds and CStoSliceIndsToAddByCnt:
    NewCStoSliceInds = OrigCStoSliceInds + CStoSliceIndsToAddByCnt
    
    #print(f'\nNewCStoSliceInds = {NewCStoSliceInds}')
    
    # Remove any duplicate slice numbers:
    NewCStoSliceInds = UniqueItems(NewCStoSliceInds)
    
    if LogToConsole:
        print('   NewCStoSliceInds after removing duplicates =',
              f'{NewCStoSliceInds}')
    
    #C = len(NewCStoSliceInds)
    
    #print(f'\nOrigPtsByCnt:\n\n', OrigPtsByCnt)
    
    if LogToConsole:
        print(f'   len(PtsToAddByCnt) = {len(PtsToAddByCnt)}')
        #print(f'\nPtsToAddByCnt:\n\n', PtsToAddByCnt)
    
    NewCntDataByCnt = []
    NewPtsByCnt = []
    
    for SliceNum in NewCStoSliceInds:
        if LogToConsole:
            print(f'   SliceNum = {SliceNum}')
        
        if SliceNum in CStoSliceIndsToAddByCnt:
            """
            This contour is to be over-written by the corresponding points in
            PtsToAddByCnt.
            """
            
            ind = CStoSliceIndsToAddByCnt.index(SliceNum)
            
            if LogToConsole:
                print(f'      This contour is to be over-written with contour',
                      f'{ind} in PtsToAddByCnt')
                
                print(f'      len(PtsToAddByCnt[{ind}]) =',
                      f'{len(PtsToAddByCnt[ind])}')
            
            NewPtsByCnt.append(PtsToAddByCnt[ind])
            
            NewCntDataByCnt.append(Points2ContourData(PtsToAddByCnt[ind]))
            
        else:
            """
            This existing contour is to be preserved.
            """
            
            ind = OrigCStoSliceInds.index(SliceNum)
            
            if LogToConsole:
                print('      This contour is to be preserved (contour {ind}',
                      'from OrigPtsByCnt)')
                
                print(f'      len(OrigPtsByCnt[{ind}]) =',
                      f'{len(OrigPtsByCnt[ind])}')
            
            NewPtsByCnt.append(OrigPtsByCnt[ind])
            
            NewCntDataByCnt.append(Points2ContourData(OrigPtsByCnt[ind]))
        
        if LogToConsole:
            print(f'      len(NewPtsByCnt[-1]) = {len(NewPtsByCnt[-1])}')
    
    if LogToConsole:
        print('end of AddCopiedPtsByCnt.')
    
    return NewCntDataByCnt, NewPtsByCnt, NewCStoSliceInds










def GetIndsByRoi(Rts, DicomDir):
    """
    Get the indices of the physical points in all contours for all ROIs.
    
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
                       
                            
    Outputs:
    *******
    
    IndsByCntByRoi : list of list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each 
        dimension) of indices of the corresponding physical points.
    """
    
    #from ImageTools import GetImageAttributes
    from ImageTools import ImportImage
    from ConversionTools import ContourToIndices
    
    #PtsByRoi = GetPtsByRoi(Rts, DicomDir)
    PtsByCntByRoi = GetPtsByCntByRoi(Rts, DicomDir)
    
    #Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    RefIm = ImportImage(DicomDir)
    
    Nrois = len(PtsByCntByRoi)
    
    IndsByCntByRoi = []
    
    for r in range(Nrois):
        IndsByCnt = []
        
        # The list of points in the contours for this ROI:
        PtsByCnt = PtsByCntByRoi[r]
        
        for Pts in PtsByCnt:
            #Inds = ContourToIndices(Pts, IPPs[0], Dirs, Spacings)
            Inds = ContourToIndices(Pts, RefIm)
                
            IndsByCnt.append(Inds)
            
        IndsByCntByRoi.append(IndsByCnt)
        
    
    return IndsByCntByRoi










def GetMaskFromContoursForSliceNum(Rts, SliceNum, DicomDir, RefImage, 
                                   SearchString='', LogToConsole=False):
    """
    NOTE 20/01/21:
        This function is modified from GetMaskFromRts, but has been renamed to
        provide a clearer idea of its purpose.
        
        
    Get a 2D mask (pixel array) from all contours in a RTS that correspond to a
    slice number of interest.
       
    
    Inputs:
    ******
    
    Rts : Pydicom object
        DICOM-RTSTRUCT object.
        
    SliceNum : integer (0-indexed)
        Slice index of the DICOM stack corresponding to the contour(s) of
        interest.
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
        
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
    SearchString : string (optional; '' by default)
        All or part of the ROI name of the ROI containing the contour to be
        copied.  If '', the first ROI will be used.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
        
    Outputs:
    *******
    
    Mask : Numpy data array
        A 2D mask of the contour(s) specified by FromRoiLabel and FromSliceNum.
        
    
    Note:
    ----
        At present multiple contours on the slice of interest will be converted
        to a single 2D mask.  It might be useful to generate a distinct 2D mask
        for each contour (i.e. a 3D mask for slices that contain multiple
        contours).
    """
    
    #import importlib
    #import ConversionTools
    #importlib.reload(ConversionTools)
    
    #from ImageTools import GetImageAttributes
    from ConversionTools import Point2Index
    from ConversionTools import Indices2Mask
    
    #Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    # Get the contour points (in the Patient Coordinate System):
    #Points, CStoSliceInds = GetPtsInContour(Rts, SearchString, SliceNum, 
    #                                        DicomDir)
    #Points = GetPtsInContour(Rts, SearchString, SliceNum, DicomDir)
    PtsByCnt = GetPtsByCntForSliceNum(Rts, SliceNum, DicomDir, SearchString, 
                                      LogToConsole)
    
    #print(f'\nlen(PtsByCnt) = {len(PtsByCnt)}')
    #print(f'\nPtsByCnt = {PtsByCnt}')
    
    # Convert the points to indices (i.e. the Image Coordinate System):
    Inds = []
            
    for c in range(len(PtsByCnt)):
        #print(f'\nlen(PtsByCnt[{c}]) = {len(PtsByCnt[c])}')
        
        for point in PtsByCnt[c]:
            inds = Point2Index(point, RefImage)
            
            # Convert fractional indices to integers:
            inds = [int(round(item)) for item in inds]
            
            Inds.append(inds)
    
    Mask = Indices2Mask(Inds, RefImage)
    
    return Mask











def InitialiseRts(RtsTemplate, RoiNum, CStoSliceInds, CntDataByObjByFrame,
                  DicomDir, NamePrefix='', LogToConsole=False):
    """
    UPDATE 19/01/21:
        Since there can be more than one contour in any given frame this
        version allows for more than one ContourSequence in ROIContourSequence.
        This requires an additional input: either PtsByObjByFrame or
        CntDataByObjByFrame to gauge the required number of ContourSequences
        based on the number of objects.
        
        
        
    Initialise a RTS object based on an existing RTS object (RtsTemplate).
         
    
    Inputs:
    ******
                              
    RtsTemplate : Pydicom object
        RTS object that will be used as a template to create the new object.
    
    RoiNum : integer
        The index of the ROI that contains the contour to be copied (counting
        from 0).
    
    CIStoSliceInds : List of integers
        List (for each sequence) of slice numbers that correspond to each 
        ContourSequence.
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    Rts : Pydicom object
        New RTS object.
    """
    
    from pydicom.uid import generate_uid
    import time
    from copy import deepcopy
    from DicomTools import ImportDicom
    from ConversionTools import NumOfListsAtDepthTwo
    
    # The ROIName of the ROI containing the contour to be copied:
    RoiName = RtsTemplate.StructureSetROISequence[RoiNum].ROIName
    
    
    # Use RtsTemplate as a template for NewRts:
    Rts = deepcopy(RtsTemplate)
    
    Dicom = ImportDicom(DicomDir)
    
    # Generate a new SOPInstanceUID:
    Rts.SOPInstanceUID = generate_uid()
    
    Rts.file_meta.MediaStorageSOPInstanceUID = deepcopy(Rts.SOPInstanceUID) # 18/12
    
    Rts.StudyDate = Dicom.StudyDate
    Rts.StudyTime = Dicom.StudyTime
    Rts.SeriesDate = Dicom.SeriesDate
    Rts.PatientName = Dicom.PatientName
    Rts.PatientID = Dicom.PatientID
    Rts.PatientBirthDate = Dicom.PatientBirthDate
    Rts.PatientSex = Dicom.PatientSex
    Rts.StudyInstanceUID = Dicom.StudyInstanceUID
    
    # Generate a new SeriesInstanceUID:
    Rts.SeriesInstanceUID = generate_uid()
    
    Rts.StudyID = Dicom.StudyID
    Rts.FrameOfReferenceUID = Dicom.FrameOfReferenceUID
    
    Rts.ReferencedFrameOfReferenceSequence[0]\
       .FrameOfReferenceUID = Dicom.FrameOfReferenceUID # 14/12
       
    Rts.ReferencedFrameOfReferenceSequence[0]\
       .RTReferencedStudySequence[0]\
       .ReferencedSOPInstanceUID = Dicom.StudyInstanceUID # 14/12
    
    Rts.ReferencedFrameOfReferenceSequence[0]\
       .RTReferencedStudySequence[0]\
       .RTReferencedSeriesSequence[0]\
       .SeriesInstanceUID = Dicom.SeriesInstanceUID # 14/12
       
    """I thought the following was missing but it's done in ModifyRts():"""
    #Rts.StructureSetROISequence[0]\
    #   .ReferencedFrameOfReferenceUID = Dicom.FrameOfReferenceUID # 25/01/21 
       
    
    # Start timing:
    #times = []
    #times.append(time.time())
    
    """ 
    Modify the number of sequences in ContourImageSequence (CIS). 
    The number of sequences CIS is equal to the length of CStoSliceInds.
    """
    # The original and required number of sequences in ContourImageSequence:
    OrigCIS = len(Rts.ReferencedFrameOfReferenceSequence[0]\
                     .RTReferencedStudySequence[0]\
                     .RTReferencedSeriesSequence[0]\
                     .ContourImageSequence)
    
    ReqCIS = len(CStoSliceInds)
    
    # Add/remove sequences by appending/popping the last sequence N times, 
    # where N is the difference between ReqCIS and OrigCIS:
    if ReqCIS > OrigCIS:
        for i in range(ReqCIS - OrigCIS):
            Rts.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0]\
               .RTReferencedSeriesSequence[0]\
               .ContourImageSequence\
               .append(deepcopy(Rts.ReferencedFrameOfReferenceSequence[0]\
                                   .RTReferencedStudySequence[0]\
                                   .RTReferencedSeriesSequence[0]\
                                   .ContourImageSequence[-1]))          
    else:
        for i in range(OrigCIS - ReqCIS):
            Rts.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0]\
               .RTReferencedSeriesSequence[0]\
               .ContourImageSequence.pop()

    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'ContourImageSequence.')
        
    """ 
    Modify the number of sequences in StructureSetROISequence. 
    There must be only one.
    """
    
    OrigS = len(Rts.StructureSetROISequence)
    
    ReqS = 1
    
    if ReqS > OrigS:
        for i in range(ReqS - OrigS):
            Rts.StructureSetROISequence\
               .append(deepcopy(Rts.ReferencedFrameOfReferenceSequence[0]\
                                   .RTReferencedStudySequence[0]\
                                   .RTReferencedSeriesSequence[0]\
                                   .ContourImageSequence[-1]))
            
    else:
        for i in range(OrigS - ReqS):
            Rts.StructureSetROISequence.pop()
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'StructureSetROISequence.')
        
        
    """ 
    Modify the number of sequences in ROIContourSequence. 
    There must be only one.
    """
    
    OrigR = len(Rts.ROIContourSequence)
    
    ReqR = 1
    
    if ReqR > OrigR:
        for i in range(ReqR - OrigR):
            Rts.ROIContourSequence\
               .append(deepcopy(Rts.ReferencedFrameOfReferenceSequence[0]\
                                   .RTReferencedStudySequence[0]\
                                   .RTReferencedSeriesSequence[0]\
                                   .ContourImageSequence[-1]))
            
    else:
        for i in range(OrigR - ReqR):
            Rts.ROIContourSequence.pop()
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'ROIContourSequence.')
        
        
    """ 
    Modify the number of sequences in ContourSequence.  
    The number of sequences must be equal to the total number of objects (i.e.
    contours) in CntDataByObjByFrame.
    """
    
    # The number of sequences in the last ContourSequence:
    OrigCS = len(Rts.ROIContourSequence[-1].ContourSequence)
    
    ReqCS = NumOfListsAtDepthTwo(CntDataByObjByFrame)
    
    if ReqCS > OrigCS:
        for i in range(ReqCS - OrigCS):
            Rts.ROIContourSequence[-1].ContourSequence\
               .append(deepcopy(Rts.ROIContourSequence[-1].ContourSequence[-1]))          
    else:
        for i in range(OrigCS - ReqCS):
            Rts.ROIContourSequence[-1].ContourSequence.pop()

    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'ContourSequence.')
        
        
    """ 
    Modify the number of sequences in RTROIObservationsSequence. 
    There must be only one.
    """
    
    OrigO = len(Rts.RTROIObservationsSequence)
    
    ReqO = 1
    
    if ReqO > OrigO:
        for i in range(ReqO - OrigO):
            Rts.RTROIObservationsSequence\
               .append(deepcopy(Rts.RTROIObservationsSequence[-1]))          
    else:
        for i in range(OrigO - ReqO):
            Rts.RTROIObservationsSequence.pop()
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'RTROIObservationsSequence.')
        
    
    """ 
    Modify other tags.
    """
    #if 'RTS' in NamePrefix:
    #    SSL = NamePrefix + '_from_' + Rts.StructureSetLabel
    #else:
    #    SSL = 'RTS_' + NamePrefix + '_from_' + Rts.StructureSetLabel
        
    #Rts.StructureSetLabel = SSL
    Rts.StructureSetLabel = NamePrefix
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    Rts.StructureSetDate = NewDate
    Rts.StructureSetTime = NewTime
    
    # Modify the ROI Name to RoiName:
    Rts.StructureSetROISequence[0].ROIName = RoiName
    
    
    
    """
    Output some results to the console.
    """
    if LogToConsole:
        NewCIS = len(Rts.ReferencedFrameOfReferenceSequence[0]\
                        .RTReferencedStudySequence[0]\
                        .RTReferencedSeriesSequence[0]\
                        .ContourImageSequence)
        
        NewS = len(Rts.StructureSetROISequence)
        
        NewR = len(Rts.ROIContourSequence)
        
        NewCS = len(Rts.ROIContourSequence[-1].ContourSequence)
        
        NewO = len(Rts.RTROIObservationsSequence)
        
        print(f'\nThere were {OrigCIS} sequences in ContourImageSequence. ',
              f'Now there are {NewCIS} sequences.')
        
        print(f'\nThere were {OrigS} sequences in StructureSetROISequence.',
              f'Now there are {NewS} sequences.')
        
        print(f'\nThere were {OrigR} sequences in ROIContourSequence.',
              f'Now there are {NewR} sequences.')
        
        print(f'\nThere were {OrigCS} sequences in ContourSequence.',
              f' Now there are {NewCS} sequences.')
        
        print(f'\nThere were {OrigO} sequences in RTROIObservationsSequence.',
              f'Now there are {NewO} sequences.')

    return Rts









def ModifyRts(Rts, CntDataByObjByFrame, PtsByObjByFrame, CStoSliceInds, 
              DicomDir, LogToConsole=False):
    """
    COMMENT 19/01/21:
        Need to update the "Inputs" section below following change to input
        arguments and change to PixArr2PtsByContour.
        
        
    Modify various tag values of the new RTS ROI so that they are consistent
    with the corresponding tag values of the DICOMs (e.g. SOPClassUID, 
    SOPInstanceUID), and so that the ROIContourSequence values are consistent
    with those of the original RTS ROI for existing contours, and with the
    contour that has been copied from the original RTS ROI.
         
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
    ******
                              
    Rts : Pydicom object
        The RTS ROI object to modify.
    
    CntDataByCnt : list of a list of strings
        A list (for each contour) of a flat list of [x, y, z] coordinates
        of the polygons that define each contour as strings.
        
    PtsByCnt : list of a list of a list of floats
        A list (for each contour) of a list (for each point) of a list (for 
        each dimension) of the polygons that define each contour.
    
    CStoSliceInds : List of integers
        List (for each contour) of slice numbers that correspond to each 
        list of contour data in CntDataByCnt.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    Rts : Pydicom object
        Modified RTS ROI object.
    """
    
    import time
    from DicomTools import ImportDicoms
    
    Dicoms = ImportDicoms(DicomDir)
    
    if False:
        CIS = Rts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence
               
        print(f'\nlen(CIS) = {len(CIS)}')
        
        print(f'\nlen(CStoSliceInds) = {len(CStoSliceInds)}')
      
    if LogToConsole:
        print('\n\nRunning ModifyRts:')
        print(f'   CStoSliceInds = {CStoSliceInds}')
        print(f'   len(Dicoms) = {len(Dicoms)}')
    
    
    ObjNum = 0 # initial value of the object/contour number
    
    # Start timing:
    times = []
    times.append(time.time())
    
    #for i in range(len(CStoSliceInds)):
    for f in range(len(CntDataByObjByFrame)):
        # The DICOM slice index for this frame:
        s = CStoSliceInds[f]
        
        if LogToConsole:
            print(f'\n   f = {f}, CStoSliceInds[{f}] = {CStoSliceInds[f]}')
            #print(f'      DICOM slice number = {s}')
        
        # Modify the ReferencedSOPInstanceUID in ContourImageSequence:
        Rts.ReferencedFrameOfReferenceSequence[0]\
           .RTReferencedStudySequence[0]\
           .RTReferencedSeriesSequence[0]\
           .ContourImageSequence[f]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        ## Modify the ReferencedSOPInstanceUID in ROIContourSequence: (20/01/21)
        #Rts.ROIContourSequence[0]\
        #   .ContourSequence[i]\
        #   .ContourImageSequence[0]\
        #   .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # The contour points for this frame:
        PointsByObj = PtsByObjByFrame[f]
        CntDataByObj = CntDataByObjByFrame[f]
        
        for o in range(len(PointsByObj)):
            # Modify the ReferencedSOPInstanceUID in ROIContourSequence: (20/01/21) 
            Rts.ROIContourSequence[0]\
               .ContourSequence[ObjNum]\
               .ContourImageSequence[0]\
               .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
           
            # Modify NumberOfContourPoints:
            Rts.ROIContourSequence[0]\
               .ContourSequence[ObjNum]\
               .NumberOfContourPoints = f"{len(PointsByObj[o])}"
            
            # The contour points as a flat list of coordinates as strings:
            Rts.ROIContourSequence[0]\
               .ContourSequence[ObjNum]\
               .ContourData = CntDataByObj[o]
            
            if LogToConsole:
                print(f'      len(PointsByObj[{o}]) = {len(PointsByObj[o])}')
            
            # Increment the object/contour number:
            ObjNum = ObjNum + 1
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to modify ContourData.')
        
        
    # Modify the StructureSetROISequence:
    Rts.StructureSetROISequence[0].RoiNumber = f"1"
        
    Rts.StructureSetROISequence[0]\
       .ReferencedFrameOfReferenceUID = Dicoms[0].FrameOfReferenceUID
       
    # Modify ROIContourSequence:
    Rts.ObservationNumber = f"{1}"
        
    Rts.ReferencedROINumber = f"{1}"
    
    if LogToConsole:
        print('end of ModifyRts.')
    
    return Rts







def ErrorCheckRts(Rts, DicomDir, LogToConsole=False):
    """
    Check a RTS for errors in dependencies based on provided directory of
    the DICOMs that relate to the RTS.  
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Rts.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe info or any errors that are found.
    
    Nerrors : integer
        The number of errors found.
    """        
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import ImportDicom
    from ImageTools import GetImageAttributes
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    Size, Spacings, ST,\
    IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    
    
    LogList = []
    Nerrors = 0
    
    """Determine whether MediaStorageSOPInstanceUID matches SOPInstanceUID."""
    MSSOPuid = deepcopy(Rts.file_meta.MediaStorageSOPInstanceUID)
    SOPuid = deepcopy(Rts.SOPInstanceUID)
    
    if MSSOPuid == SOPuid:
        msg = f'INFO:  MediaStorageSOPInstanceUID {MSSOPuid} does not '\
              + f'match SOPInstanceUID {SOPuid}.\n'
    else:
        msg = f'ERROR:  MediaStorageSOPInstanceUID {MSSOPuid} does not '\
              + f'match SOPInstanceUID {SOPuid}.\n'
              
        Nerrors = Nerrors + 1
    
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    if Rts.PatientName == Dicom.PatientName:
        msg = f'INFO:  RTS PatientName {Rts.PatientName} matches DICOM '\
              + f'PatientName {Dicom.PatientName}.\n'
    else:
        msg = f'ERROR:  RTS PatientName {Rts.PatientName} does not match '\
              + f'DICOM PatientName {Dicom.PatientName}.\n'
        
        Nerrors = Nerrors + 1
    
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
        
    if Rts.PatientID == Dicom.PatientID:
        msg = f'INFO:  RTS PatientID {Rts.PatientID} matches '\
              + f'DICOM PatientID {Dicom.PatientID}.\n'
    else:
        msg = f'ERROR:  RTS PatientID {Rts.PatientID} does not match '\
              + f'DICOM PatientID {Dicom.PatientID}.\n'
              
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
        
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    if Rts.StudyInstanceUID == Dicom.StudyInstanceUID:
        msg = f'INFO:  RTS StudyInstanceUID {Rts.StudyInstanceUID} '\
              + f'matches DICOM StudyInstanceUID {Dicom.StudyInstanceUID}.\n'
    else:
        msg = f'ERROR:  RTS StudyInstanceUID {Rts.StudyInstanceUID} does not '\
              + f'match DICOM StudyInstanceUID {Dicom.StudyInstanceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    #if Rts.StudyID == Dicom.StudyID:
    #    msg = f'INFO:  The RTS StudyID {Rts.StudyID} matches '\
    #          + f'the DICOM StudyID {Dicom.StudyID}.\n'
    #else:
    #    msg = f'ERROR:  The RTS StudyID {Rts.StudyID} does not match '\
    #          + f'the DICOM StudyID {Dicom.StudyID}.\n'
    #
    #    Nerrors = Nerrors + 1
    #    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #    
    #if LogToConsole:
    #    print(msg)
            
    if Rts.FrameOfReferenceUID == Dicom.FrameOfReferenceUID:
        msg = f'INFO:  RTS FrameOfReferenceUID {Rts.FrameOfReferenceUID}'\
              + '  matches the DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = f'ERROR:  RTS FrameOfReferenceUID {Rts.FrameOfReferenceUID}'\
              + '  does not match the DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
      
    
    #RFORS = deepcopy(Rts.ReferencedFrameOfReferenceSequence[0])
    RFORS = deepcopy(Rts.ReferencedFrameOfReferenceSequence)
    
    #FORuidInRFORS = deepcopy(RFORS.FrameOfReferenceUID)
    FORuidInRFORS = deepcopy(RFORS[0].FrameOfReferenceUID)
    
    if FORuidInRFORS == Dicom.FrameOfReferenceUID:
        msg = 'INFO:  FrameOfReferenceUID in '\
              + f'ReferencedFrameOfReferenceSequence {FORuidInRFORS} '\
              + 'matches DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = 'ERROR:  FrameOfReferenceUID in '\
              + f'ReferencedFrameOfReferenceSequence {FORuidInRFORS} does'\
              + ' not match DICOM FrameOfReferenceUID '\
              + f'{Dicom.FrameOfReferenceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Verify that the ReferencedSOPInstanceUID in 
    ReferencedFrameOfReferenceSequence matches the RTS StudyInstanceUID."""
    #RefSOPuid = RFORS.RTReferencedStudySequence[0].ReferencedSOPInstanceUID
    RefSOPuid = RFORS[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID
    
    if RefSOPuid == Rts.StudyInstanceUID:
        msg = 'INFO:  ReferencedSOPInstanceUID in '\
              + f'ReferencedFrameOfReferenceSequence {RefSOPuid} matches '\
              + f'the RTS StudyInstanceUID {Rts.StudyInstanceUID}.\n'
    else:
        msg = 'ERROR:  ReferencedSOPInstanceUID in '\
              + f'ReferencedFrameOfReferenceSequence {RefSOPuid} does not'\
              + f' match the RTS StudyInstanceUID {Rts.StudyInstanceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    """Verify that the SeriesInstanceUID in ReferencedFrameOfReferenceSequence
    matches the DICOM SeriesInstanceUID."""
    #RtsSeriesuid = deepcopy(RFORS.RTReferencedStudySequence[0]\
    #                             .RTReferencedSeriesSequence[0]\
    #                             .SeriesInstanceUID)
    RtsSeriesuid = deepcopy(RFORS[0].RTReferencedStudySequence[0]\
                                    .RTReferencedSeriesSequence[0]\
                                    .SeriesInstanceUID)
    
    if RtsSeriesuid == Dicom.SeriesInstanceUID:
        msg = 'INFO:  SeriesInstanceUID in '\
              + f'ReferencedFrameOfReferenceSequence {RtsSeriesuid} '\
              + 'matches DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
    else:
        msg = 'ERROR:  SeriesInstanceUID in '\
              + f'ReferencedFrameOfReferenceSequence {RtsSeriesuid} does '\
              + 'not match DICOM SeriesInstanceUID '\
              + f'{Dicom.SeriesInstanceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    
    """Verify that the ReferencedSOPInstanceUIDs in 
    ReferencedFrameOfReferenceSequence match the DICOM SOPInstanceUIDs."""
    #CIS = deepcopy(RFORS.RTReferencedStudySequence[0]\
    #                    .RTReferencedSeriesSequence[0]\
    #                    .ContourImageSequence)
    CIS = deepcopy(RFORS[0].RTReferencedStudySequence[0]\
                           .RTReferencedSeriesSequence[0]\
                           .ContourImageSequence)
    
    RefSOPuidsInRFORS = [CIS[i].ReferencedSOPInstanceUID for i in range(len(CIS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [RefSOPuid in SOPuids for RefSOPuid in RefSOPuidsInRFORS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if Inds:
        NonMatchingUids = [RefSOPuidsInRFORS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in '\
                  + 'ContourImageSequence {i+1} does not match any DICOM '\
                  + 'SOPInstanceUID.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in '\
              + 'ContourImageSequence match the DICOM '\
              + 'SOPInstanceUIDs.\n'
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
            
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
            
            
    """Verify that the StructureSetROISequence has only 1 sequence."""
    N = len(Rts.StructureSetROISequence)
    
    if N == 1:
        msg = f'INFO:  There are {N} sequences in StructureSetROISequence.'\
              + ' There should be 1.\n'
    else:
        msg = f'ERROR:  There are {N} sequences in StructureSetROISequence.'\
              + ' There should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
    
    
    """Verify that the ROINumber in StructureSetROISequence is 1."""
    N = int(deepcopy(Rts.StructureSetROISequence[0].ROINumber))
    
    if N == 1:
        msg = f'INFO:  ROINumber in StructureSetROISequence is {N}.  It '\
              + 'should be 1.\n'
    else:
        msg = f'ERROR:  ROINumber in StructureSetROISequence is {N}.  It '\
              + 'should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    
    """Verify that the ReferencedFrameOfReferenceUID in StructureSetROISequence
    matches the DICOM FrameOfReferenceUID."""
    RefFORuidInSSRS = deepcopy(Rts.StructureSetROISequence[0]\
                                  .ReferencedFrameOfReferenceUID)
    
    if RefFORuidInSSRS == Dicom.FrameOfReferenceUID:
        msg = 'INFO:  ReferencedFrameOfReferenceUID in '\
              + f'StructureSetROISequence {RefFORuidInSSRS} matches DICOM '\
              + f'FrameOfReferenceUID {Dicom.FrameOfReferenceUID}.\n'
    else:
        msg = 'ERROR:  ReferencedFrameOfReferenceUID in '\
              + f'StructureSetROISequence {RefFORuidInSSRS} does not match '\
              + f'DICOM FrameOfReferenceUID {Dicom.FrameOfReferenceUID}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
            
    
    """Verify that the ROIContourSequence has only 1 sequence."""
    N = len(Rts.ROIContourSequence)
    
    if N == 1:
        msg = f'INFO:  There are {N} sequences in ROIContourSequence. There '\
              + 'should be 1.\n'
    else:
        msg = f'ERROR:  There are {N} sequences in ROIContourSequence. There '\
              + 'should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)

    
    
    """Verify that the number of sequences in ContourSequence matches the
    number of sequences in ContourImageSequence.
    
    19/01/21:
        There's no reason why the above need be true.  The number of sequences
        in ContourSequence must be equal to the number of contours/objects in 
        TrgCntDataByObjByFrame.  Commenting the code below."""
    
    CS = deepcopy(Rts.ROIContourSequence[0].ContourSequence)
    
    #print(f'\nlen(CS) = {len(CS)}, len(CIS) = {len(CIS)}')
    """
    if len(CS) == len(CIS):
        msg = 'INFO:  The number of sequences in ContourSequence '\
              + f'{len(CS)} matches the number of sequences in '\
              + f'ContourImageSequence {len(CIS)}.\n'
    else:
        msg = 'ERROR:  The number of sequences in ContourSequence '\
              + f'{len(CS)} does not match the number of sequences in '\
              + f'ContourImageSequence {len(CIS)}.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    """
    
    
    """Verify that the ReferencedSOPInstanceUIDs in ContourSequence
    match the ReferencedSOPInstanceUIDs in ContourImageSequence."""
    RefSOPuidsInRCS = [CS[i].ContourImageSequence[0].ReferencedSOPInstanceUID for i in range(len(CS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [uid in RefSOPuidsInRFORS for uid in RefSOPuidsInRCS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if Inds:
        NonMatchingUids = [RefSOPuidsInRCS[i] for i in Inds]
        
        for i in range(len(NonMatchingUids)):
            uid = NonMatchingUids[i]
            
            msg = f'ERROR:  ReferencedSOPInstanceUID {uid} in ContourSequence'\
                  + f'{i+1} does not match any ReferencedSOPInstanceUID in '\
                  + 'ContourImageSequence.\n'
                  
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ReferencedSOPInstanceUIDs in ContourSequence'\
                  + f' match the ReferencedSOPInstanceUIDs in '\
                  + 'ContourImageSequence.\n'
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
                
    
    
    """Verify that the ContourNumbers in each ContourSequence range from 1 to
    len(CS)."""
    ExpectedNums = list(range(1, len(CS)+1))
    
    ContourNums = [int(CS[i].ContourNumber) for i in range(len(CS))]
    
    # Determine whether the contour numbers match the expected contour numbers:
    IsMatch = [ContourNums[i] == ExpectedNums[i] for i in range(len(CS))]
    #IsMatch = [ContourNums[i] == i+1 for i in range(len(CS))]
    
    # Find the indices of any non-matching contour numbers:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if Inds:        
        for i in range(len(Inds)):
            cn = ContourNums[Inds[i]]
            en = ExpectedNums[Inds[i]]
            
            msg = f'ERROR:  ContourNumber {i+1} is {cn} but is expected to be'\
                  + f'{en}.\n'
                  
            Nerrors = Nerrors + 1
    else:
        msg = f'INFO:  The ContourNumbers are as expected.\n'
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    
    """Verify that the number of elements in ContourData is 3x 
    NumberOfContourPoints for each ContourSequence."""
    NCP = [int(CS[i].NumberOfContourPoints) for i in range(len(CS))]
    
    NCD = [len(CS[i].ContourData) for i in range(len(CS))]
    
    IsMatch = [NCD[i] == 3*NCP[i] for i in range(len(CS))]
    
    # Find the indices of any non-zero modulos:
    Inds = [i for i, x in enumerate(IsMatch) if False]
    
    if Inds:        
        for i in range(len(Inds)):
            ncp = NCP[Inds[i]]
            
            ncd = NCD[Inds[i]]
            
            msg = f'ERROR:  The number of elements in ContourData {ncd} is '\
                  + f'not 3x NumberOfContourPoints {ncp} for '\
                  + f'ContourSequence {i+1}.\n'
            
            Nerrors = Nerrors + 1
    else:
        msg = 'INFO:  The number of elements in ContourData are 3x '\
              + 'NumberOfContourPoints for ContourSequence.\n'
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    
    """Verify that the ReferencedROINumber in ROIContourSequence is 1."""
    N = int(deepcopy(Rts.ROIContourSequence[0].ReferencedROINumber))
    
    if N == 1:
        msg = f'INFO:  ReferencedROINumber in ROIContourSequence is {N}. '\
              + 'It should be 1.\n'
    else:
        msg = f'ERROR:  ReferencedROINumber in ROIContourSequence is {N}. '\
              + 'It should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    
    """Verify that the ROIObservationsSequence has only 1 sequence."""
    N = len(Rts.RTROIObservationsSequence)
    
    if N == 1:
        msg = f'INFO:  There are {N} sequences in RTROIObservationsSequence.'\
              + ' There should be 1.\n'
    else:
        msg = f'ERROR:  There are {N} sequences in RTROIObservationsSequence.'\
              + ' There should be 1.\n'
              
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)


    
    """Verify that the ObservationNumber is 1."""
    N = int(deepcopy(Rts.RTROIObservationsSequence[0].ObservationNumber))
    
    if N == 1:
        msg = 'INFO:  ObservationNumber in RTROIObservationsSequence is '\
              + f'{N}. It should be 1.\n'
    else:
        msg = 'ERROR:  ObservationNumber in RTROIObservationsSequence is '\
              + f'{N}. It should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
            
            
            
    """Verify that the ReferencedROINumber is 1."""
    N = int(deepcopy(Rts.RTROIObservationsSequence[0].ReferencedROINumber))
    
    if N == 1:
        msg = 'INFO:  ReferencedROINumber in RTROIObservationsSequence is '\
              + f'{N}. It should be 1.\n'
    else:
        msg = 'ERROR:  ReferencedROINumber in RTROIObservationsSequence is '\
              + f'{N}. It should be 1.\n'
        
        Nerrors = Nerrors + 1
    
    #msg = msg + f'Nerrors = {Nerrors}.\n'
    
    LogList.append(msg)
    
    if LogToConsole:
        print(msg)
    
    
    print(f'\nThere were {Nerrors} errors found in the RTS.')
    
    return LogList, Nerrors









def GetIndsOfRefSOPsFromCSToCIS(Rts):
    """
    Get the indices of the ReferencedSOPInstanceUIDs in the 
    ContourImageSequences for each ReferencedSOPInstanceUID in the 
    ContourSequences.
    
    Returns -1 if a ReferencedSOPInstanceUID in ContourSequences is not found
    in ContourImageSequence.
    """
    
    CIS = Rts.ReferencedFrameOfReferenceSequence[0]\
             .RTReferencedStudySequence[0]\
             .RTReferencedSeriesSequence[0]\
             .ContourImageSequence
               
    RefSOPsInCIS = [CIS[i].ReferencedSOPInstanceUID for i in range(len(CIS))]
    
    
    
    CS = Rts.ROIContourSequence[0]\
            .ContourSequence
    
    RefSOPsInCS = [CS[i].ContourImageSequence[0]\
                        .ReferencedSOPInstanceUID for i in range(len(CS))]
       
    Inds = []
    
    for i in range(len(CS)):
        if RefSOPsInCS[i] in RefSOPsInCIS:
            Inds.append(RefSOPsInCIS.index(RefSOPsInCS[i]))
        else:
            Inds.append(-1)
    
    return Inds







def GetContoursFromListOfRtss(ListOfRtss, ListOfDicomDirs, SearchString='', 
                              LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    from DicomTools import GetRoiNum
    from DicomTools import GetDicomFpaths
    from ImageTools import GetImageAttributes
    
    ListOfPtsByCnt = []
    ListOfC2Sinds = []
    #ListOfPtsByCntBySlice = []
    ListOfRoiNums = []
    ListOfDcmFpaths = []
    ListOfOrigins = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfRtss)):
        if LogToConsole:
            print(f'\nListOfRtss[{i}]:')
        
        """
        11/12/20:  
            There should only be one ROI so perhaps this should be checked
            and raise exception if the number of ROIs is greater than one.
            And if not, RoiNum will always be = 0 (no need to GetRoiNum).
        """
        
        if ListOfRtss[i]:
            RoiNum = GetRoiNum(ListOfRtss[i], SearchString)
            
            PtsByCnt, C2Sinds = GetPtsByCntForRoi(ListOfRtss[i], 
                                                  ListOfDicomDirs[i], 
                                                  SearchString,
                                                  LogToConsole)
            
            #PtsByCntBySlice = GetPtsByCntBySliceNum(ListOfRtss[i], 
            #                                        ListOfDicomDirs[i], 
            #                                        SearchString,
            #                                        LogToConsole)
            
            if LogToConsole:      
                print(f'\nlen(PtsByCnt) in "{SearchString}" = {len(PtsByCnt)}')
                print(f'C2Sinds = {C2Sinds}')
                [print(f'len(PtsByCnt{i}) = {len(PtsByCnt[i])}') for i in range(len(PtsByCnt))]
                #print(f'\nlen(PtsByCntBySlice) in "{SearchString}" =',
                #      f'{len(PtsByCntBySlice)}')
                #[print(f'len(PtsByCntBySlice{i}) = {len(PtsByCntBySlice[i])}') for i in range(len(PtsByCntBySlice))]
            
        else:
            RoiNum = None
            PtsByCnt = None
            C2Sinds = None
            #PtsByCntBySlice = None
            
            if LogToConsole:      
                print(f'No contours for this dataset')
        
        
        ListOfRoiNums.append(RoiNum)
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(ListOfDicomDirs[i])
        
        #ListOfPtsByCntBySlice.append(PtsByCntBySlice)
        ListOfPtsByCnt.append(PtsByCnt)
        ListOfC2Sinds.append(C2Sinds)
        ListOfDcmFpaths.append(DcmFpaths)
        ListOfOrigins.append(IPPs[0])
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
        
        if LogToConsole:      
            #print(f'len(Contours) in "{SearchString}" = {len(Contours)}')
            print(f'\nlen(ListOfPtsByCnt) = {len(ListOfPtsByCnt)}')
            #print(f'C2Sinds = {C2Sinds}')
            #[print(f'len(Contours{i}) = {len(Contours[i])}') for i in range(len(Contours))]
            #print(f'\nlen(ListOfPtsByCntBySlice) = {len(ListOfPtsByCntBySlice)}')
        
        
        #SOPuids = GetDicomSOPuids(ListOfDicomDirs[i])
        
        #CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(ListOfRtss[i], SOPuids)
        
        #ListOfC2Sinds.append(CStoSliceIndsByRoi[RoiNum])
        
        #if LogToConsole:
        #    print(f'ListOfC2Sinds[{i}] = {ListOfC2Sinds[i]}')
    
            
    #return ListOfPtsByCntBySlice, ListOfPtsByCnt, ListOfC2Sinds, ListOfRoiNums,\
    #       ListOfDcmFpaths, ListOfOrigins, ListOfDirections, ListOfSpacings
    return ListOfPtsByCnt, ListOfC2Sinds, ListOfRoiNums,\
           ListOfDcmFpaths, ListOfOrigins, ListOfDirections, ListOfSpacings





def CentroidOfContour(Contour):
    """
    Compute the centroid of a list of points (or indices).  
    
    Inputs:
    ******
    
    Contour : list of list of floats/integers
        List of a list of points or indices.
    
    
    Outputs:
    *******
    
    Centroid : list of floats/integers
        Centroid of Contour.
    """
    
    import numpy as np
    
    Contour = np.array(Contour)
    
    return Contour.mean(axis=0)






def LabelContoursByCentroid(Contours): # incomplete
    """
    
    This function is incomplete 
    
    Given a list of points (or indices) assign an integer label to each contour 
    using the contour centroids as a best guess of which ROI they belong to.  
    
    Inputs:
    ******
    
    Contours : list of list of a list of floats/integers
        A list of a list of points/indices for each contour.
    
    
    Outputs:
    *******
    
    RoiNums : list of integers
        An integer label assigning each contour to a ROI.
    """
    
    from GeneralTools import Distance
    
    Centroids = []
    
    for Contour in Contours:
        Centroids.append(CentroidOfContour(Contour))
        
        
    RoiNums = []
    RoiCentroids = []
    
    # Assign the first contour the number 0, and the first RoiCentroid 
    # Centroids[0]:
    RoiNums.append(0)
    RoiCentroids.append(Centroids[0])
    
    
    # Loop through Centroids starting from the 2nd one:
    for c in range(1, len(Centroids)):
        distances = []
        
        for r in range(len(RoiCentroids)):
            # Get the distance between the centroids but ignore the z-component:
            distance = Distance(Centroids[c][:2], RoiCentroids[r][:2])
            
            distances.append(distance)
            
            
        # Find index of minimum distance:
        ind = distances.index(min(distances))
        
        #if ind in RoiNums:
            
            
        
        
    return RoiNums