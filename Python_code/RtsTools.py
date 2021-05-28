# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:55:17 2020

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
DICOM RTSTRUCT FUNCTIONS
******************************************************************************
******************************************************************************
"""


def GetRtsDataOfInterest(RtsFpath, SliceNum, RoiName, DicomDir, 
                         LogToConsole=False):
    """
    Get data of interest from an RTS.
    
    Inputs:
    ******
    
    RtsFpath : string
        Filepath of Source RTSTRUCT file.
        
    SliceNum : integer or None
        The slice indeces within the DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SliceNum = None, a relationship-preserving copy will be made.
        
    RoiName : string
        All or part of the ROIName of the ROI containing the contour(s) to be
        copied.
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates to be copied from 
        the RTS.
        
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour to be copied from the RTS.
        
        
    Notes:
    *****
    
    There are 3 possible main use cases:
        1. Copy a single contour
        2. Copy an entire ROI
        3. Copy all ROIs
        
    Which main case applies depends on the inputs:
        SliceNum
        RoiName
    
    If SliceNum != None (i.e. if a slice number is defined) a single contour
    will be copied from the ROI that matches RoiName (--> Main Use Case 1).  
    
    Need to decide what to do if there are multiple ROIs and 
    RoiName = None.  Possible outcomes:

        - Raise exception with message "There are multiple ROIs so the ROI
        label must be provided"
        - Copy the contour to all ROIs that have a contour on slice 
        SliceNum (* preferred option?)
        
    If SliceNum = None but RoiName != None, all contours within the ROI given
    by FromRoiName will be copied (--> Main Use Case 2).  
    
    If SliceNum = None and RoiName = None, all contours within all ROIs will be
    copied (--> Main Use Case 3).
    
    
    If SliceNum = None and RoiName = None, all ROIs will be copied.  
    If not, get the ROI number (zero-indexed) to be copied, and reduce 
    PtsByCntByRoi and C2SindsByRoi to the corresponding item.
    """
    
    #import DicomTools
    #import importlib
    #importlib.reload(DicomTools)
    from DicomTools import GetRoiNums, GetRoiLabels
    from copy import deepcopy
    from pydicom import dcmread
    
    Rts = dcmread(RtsFpath)
    
    # Set LogToConsole for functions within this function:
    LtoC = LogToConsole 
    #LtoC = False
    
    AllPtsByCntByRoi, AllC2SindsByRoi = GetPtsByCntByRoi(Rts, DicomDir, LtoC)
    
    Nrois = len(AllC2SindsByRoi) # number of ROIs
    
    if LogToConsole:
        from GeneralTools import PrintIndsByRoi, PrintPtsByCntByRoi
        
        print('\n\n', '-'*120)
        print('Results of GetRtsDataOfInterest():')
        print('   AllC2SindsByRoi =')
        PrintIndsByRoi(AllC2SindsByRoi)
        
    
    ReducedBySlice = False
            
    if SliceNum:
        """ Limit data to those which belong to the chosen slice number. """
        
        PtsByCntByRoi = []
        C2SindsByRoi = []
    
        for r in range(Nrois):
            """ All points-by-contour and contour-to-slice indices for this 
            ROI: """
            AllPtsByCnt = deepcopy(AllPtsByCntByRoi[r])
            AllC2Sinds = deepcopy(AllC2SindsByRoi[r])
                
            """ Get the contour number(s) that relate to FromSliceNum: """
            FromCntNums = [i for i, e in enumerate(AllC2Sinds) if e==SliceNum]
            
            #print(f'\nFromCntNums = {FromCntNums}')
            #print(f'AllC2Sinds = {AllC2Sinds}')
            
            if len(FromCntNums) != len(AllC2Sinds):
                """ Some contours will be rejected. """
                ReducedBySlice = True
            
                if FromCntNums:
                    """ Keep only the indeces and points that relate to 
                    FromCntNums: """
                    PtsByCnt = [AllPtsByCnt[c] for c in FromCntNums]
                    C2Sinds = [AllC2Sinds[c] for c in FromCntNums]
                else:
                    PtsByCnt = []
                    C2Sinds = []
                
                PtsByCntByRoi.append(PtsByCnt)
                C2SindsByRoi.append(C2Sinds)
            #else:
                # All contours to remain.
                
        
        if ReducedBySlice:
            """ Replace AllPtsByCntByRoi and AllC2SindsByRoi with PtsByCntByRoi  
            and C2SindsByRoi in case further restricting of data is required 
            below: """
            AllPtsByCntByRoi = deepcopy(PtsByCntByRoi)
            AllC2SindsByRoi = deepcopy(C2SindsByRoi)
        #else:
            # AllPtsByCntByRoi and AllC2SindsByRoi remain unchanged.
        
        if LogToConsole and ReducedBySlice:
            print('\n   After limiting data to those that relate to slice',
                  f'number {SliceNum}, the C2SindsByRoi =')
            PrintIndsByRoi(C2SindsByRoi)
    
    
    #ReducedByRoi = False
    
    if RoiName:
        """ Limit data to those which belong to the chosen ROI(s). """
        
        PtsByCntByRoi = []
        C2SindsByRoi = []
        
        """ Get the ROI number(s) whose name matches RoiName: """
        RoiNums = GetRoiNums(Roi=Rts, SearchString=RoiName)
        
        #print(f'\nRoiNums = {RoiNums}')
        #print(f'AllC2SindsByRoi = {AllC2SindsByRoi}')
        
        if len(RoiNums) != len(AllC2SindsByRoi):
            #ReducedByRoi = True
            
            PtsByCntByRoi = [AllPtsByCntByRoi[r] for r in RoiNums]
            C2SindsByRoi = [AllC2SindsByRoi[r] for r in RoiNums]
            
            if LogToConsole:
                """ Get the names of all ROIs: """
                AllRoiNames = GetRoiLabels(Roi=Rts)
            
                """ Limit the list of ROI names to those that belong to the chosen 
                ROIs(s): """
                RoiNames = [AllRoiNames[i] for i in RoiNums]
            
                print('\n   After limiting data to those whose ROI name',
                      f'matches {RoiNames}, the C2SindsByRoi =')
                PrintIndsByRoi(C2SindsByRoi)
                
        else:
            PtsByCntByRoi = deepcopy(AllPtsByCntByRoi)
            C2SindsByRoi = deepcopy(AllC2SindsByRoi)
    else:
            PtsByCntByRoi = deepcopy(AllPtsByCntByRoi)
            C2SindsByRoi = deepcopy(AllC2SindsByRoi)
    
    if LogToConsole:
        print('\n   Final outputs of GetRtsDataOfInterest():')
        PrintPtsByCntByRoi(PtsByCntByRoi)
        PrintIndsByRoi(C2SindsByRoi)
        print('-'*120)
        
    return PtsByCntByRoi, C2SindsByRoi








def ProportionOfContourInExtent(Points, TrgDicomDir, LogToConsole=False):
    """
    COMMENT 26/01/2021:
        I was previously mostly getting 100% proportions and the odd 0%.  Then
        most or all results were 0%.  After recent revisions I get 100% for
        a case where the contour points/segmentation pixels are clearly not 
        within the Target image domain.
        
    
    Determine what proportion of voxels that make up a 3D pixel array intersect
    with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    Points : list of a list of floats
        A list (for each point) of a list (for each dimension) of coordinates, 
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of the points in Contour 
        that intersect the extent of RefImage.
    """
    
    from ImageTools import ImportImage, GetImageVertices
    from GeneralTools import IsPointInPolygon
    
    Image = ImportImage(TrgDicomDir)
    
    Vertices = GetImageVertices(Image)
    
    IsInside = []
    
    for pt in Points:
        IsInside.append(IsPointInPolygon(pt, Vertices))
    
    FracProp = sum(IsInside)/len(IsInside)
    
    
    if LogToConsole:
        print('\n\n\nFunction ProportionOfContourInExtent():')
        print(f'   Contour has {len(Points)} points')
        
        print(f'\n   The vertices of the image grid:')
        for i in range(len(Vertices)):
            print(f'   {Vertices[i]}')
        N = len(IsInside)
        limit = 10
        if N > limit:
            print(f'\n   The first {limit} points:')
            for i in range(limit):
                print(f'   {Points[i]}')
        else:
            print(f'\n   The first {N} points:')
            for i in range(N):
                print(f'   {Points[i]}')
                
            
        if FracProp < 1:
            print(f'\n   The vertices of the image grid:')
            for i in range(len(Vertices)):
                print(f'   {Vertices[i]}')
            N = len(IsInside)
            limit = 3
            if N > limit:
                print(f'\n   The first {limit} points that lie outside of extent:')
                for i in range(limit):
                    print(f'   {Points[i]}')
            else:
                print(f'\n   The first {N} points that lie outside of extent:')
                for i in range(N):
                    print(f'   {Points[i]}')
                    
    return FracProp






def ProportionOfRoisInExtent(PtsByCntByRoi, TrgDicomDir, LogToConsole=False):
    """
    COMMENT 26/01/2021:
        I was previously mostly getting 100% proportions and the odd 0%.  Then
        most or all results were 0%.  After recent revisions I get 100% for
        a case where the contour points/segmentation pixels are clearly not 
        within the Target image domain.
        
    
    Determine what proportion of the points in a list of points by contour by
    ROI intersect with the 3D volume of a 3D image.
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
    TrgDicomDir : string
        Directory containing the DICOMs that relate to the image grid/extent
        within which the test is for.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    FracProp : float
        The fractional proportion (normalised to 1) of the points in Contour 
        that intersect the extent of RefImage.
    """
    
    from ImageTools import ImportImage, GetImageVertices
    from GeneralTools import IsPointInPolygon
    
    
    Image = ImportImage(TrgDicomDir)
    
    Vertices = GetImageVertices(Image)
    
    IsInside = []
    
    # Iterate through ROIs:
    for r in range(len(PtsByCntByRoi)):
        # Iterate through contours:
        for c in range(len(PtsByCntByRoi[r])):
            # Iterate through points:
            for p in range(len(PtsByCntByRoi[r][c])):
                IsInside.append(IsPointInPolygon(PtsByCntByRoi[r][c][p], 
                                                 Vertices))
    
    FracProp = sum(IsInside)/len(IsInside)
                    
    return FracProp







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
    Get all ReferencedSOPInstanceUID in ContourImageSequence of a RT-STRUCT.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
    
        
    Outputs:
    *******
    
    RSOPuids : list of strings
        List of Referenced SOP UIDs.
    
    
    Notes:
    *****
    
    ContourImageSequence = Rts.ReferencedFrameOfReferenceSequence[0]\
                              .RTReferencedSeriesSequence[0]\
                              .ContourImageSequence
    
    and
    
    ReferencedSOPInstanceUID = Rts.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[i]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourImageSequence.
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





def GetRSOPuidsInCS(Rts):
    """
    Get all ReferencedSOPInstanceUID in ContourSequence of a RT-STRUCT.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
    
        
    Outputs:
    *******
    
    RSOPuids : list of strings
        List of Referenced SOP UIDs.
    
    
    Notes:
    *****
    
    ContourSequence = Rts.ROIContourSequence[0]\
                         .ContourSequence
    
    and
    
    ReferencedSOPInstanceUID = Rts.ROIContourSequence[0]\
                                  .ContourSequence[i]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourSequence.
    """
    
    RSOPuids = [] 
    
    sequences = Rts.ROIContourSequence[0]\
                   .ContourSequence
    
    for sequence in sequences:
        uid = sequence.ContourImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids






def GetCIStoSliceInds(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each ReferencedSOPInstanceUID in
    each ContourImageSequence in a RTSTRUCT ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    CIS2Sinds : list of integers
        List of slice numbers that correspond to each sequence in
        ContourImageSequence.
        
    
    Notes:
    *****
    
    ContourImageSequence = Rts.ReferencedFrameOfReferenceSequence[0]\
                              .RTReferencedSeriesSequence[0]\
                              .ContourImageSequence
    
    and
    
    ReferencedSOPInstanceUID = Rts.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[i]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourImageSequence.
    """
    
    # Get the ReferencedSOPInstanceUIDs from ContourImageSequence:
    RSOPuids = GetRSOPuidsInCIS(Rts)
    
    CIS2Sinds = [] 
    
    for RefUid in RSOPuids:
        # Find the matching index of RefUid in SOPuids:
        CIS2Sinds.append(SOPuids.index(RefUid))
        
    return CIS2Sinds





def GetCStoSliceInds(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each ReferencedSOPInstanceUID in
    each ContourSequence in a RTSTRUCT ROI.
    
    WARNING:
    *******
        If there's more than one ROI in Rts it will only return the CS2Sinds
        for the first ROI.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from a RTS file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    *******
    
    CS2Sinds : list of integers
        List of slice numbers that correspond to each sequence in
        ContourImageSequence.
    
    
    Notes:
    *****
    
    ContourSequence = Rts.ROIContourSequence[0]\
                         .ContourSequence
    
    and
    
    ReferencedSOPInstanceUID = Rts.ROIContourSequence[0]\
                                  .ContourSequence[i]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourSequence.
    """
    
    # Get the ReferencedSOPInstanceUIDs from ContourSequence:
    RSOPuids = GetRSOPuidsInCS(Rts)
    
    CS2Sinds = [] 
    
    for RefUid in RSOPuids:
        # Find the matching index of RefUid in SOPuids:
        CS2Sinds.append(SOPuids.index(RefUid))
        
    return CS2Sinds







def GetRSOPuidsInCSByRoi(Rts):
    """
    Get the list of ReferencedSOPInstance UIDs for each ContourSequence
    for each ROI in a RT-STRUCT.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    *******
    
    RSOPuidsByRoi : list of list of strings
        List (for each ROI) of a list (for each ContourSequence) of 
        ReferencedSOPInstanceUIDs.
        
        
    Notes:
    -----
    
    
    Notes:
    *****
    
    ContourSequence = Rts.ROIContourSequence[0]\
                         .ContourSequence
    
    and
    
    ReferencedSOPInstanceUID = Rts.ROIContourSequence[0]\
                                  .ContourSequence[i]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPInstanceUID
    
    for the i^th ContourSequence.
    
    If there is only one ROI Contour Sequence RSOPuids will be a list of
    length 1 of a list of strings, e.g.
        
    [['uid1', 'uid2', ...]]
    """
    
    RSOPuidsByRoi = []
    
    ROICS = Rts.ROIContourSequence
    
    for r in range(len(ROICS)):
        Sequences = ROICS[r].ContourSequence
        
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
    from GeneralTools import PrintIndsByRoi
    
    # Get the DICOM SOPInstanceUIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourImageSequence-to-slice indices:
    #CIStoSliceInds = GetCIStoSliceInds(Rts, SopUids)
    
    # Get the ContourSequence-to-slice indices by ROI:
    C2SindsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    Nrois = len(C2SindsByRoi)
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Results of GetPtsByCntByRoi():')
        #print(f'   C2SindsByRoi = {C2SindsByRoi}')
        PrintIndsByRoi(C2SindsByRoi)
        print(f'   Nrois = {Nrois}')
    
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
        
        if False:#LogToConsole:
            print(f'      r = {r}')
            print(f'      NumC2Sinds = {NumC2Sinds}')
            print(f'      NumCS = {NumCS}')
        
        
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
            
            if False:#LogToConsole:
                print(f'         c = {c}')
                print(f'         len(Pts) = {len(Pts)}')
            
            
        PtsByCntByRoi.append(PtsByCnt)
    
    if LogToConsole:
        R = len(PtsByCntByRoi)
        
        print(f'\n   len(PtsByCntByRoi) = {R}')
        
        for r in range(R):
            print(f'      len(PtsByCntByRoi[{r}]) = {len(PtsByCntByRoi[r])}')
        
        print('-'*120)

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
        Modified list (for each contour) of a flat list of [x, y, z] coordinates
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
        A 2D mask of the contour(s) specified by FromRoiName and FromSliceNum.
        
    
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








def CreateRts(SrcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi, TrgPtsByCntByRoi, 
              TrgC2SindsByRoi, TrgDicomDir, #SrcRoiName, 
              AddTxtToSSLabel='', LogToConsole=False):
    """
    Create an RTS object for the target dataset. 
    
    Inputs:
    ******
                              
    SrcRtsFpath : string
        Filepath of the Source RTS file.
    
    TrgRtsFpath : string or None
        If TrgRtsFpath != None, the Target RTS will be modified to create 
        NewTrgRts.  If TrgRtsFpath = None, the Source RTS will be used as a 
        template.
        
    TrgCntDataByCntByRoi : list of a list of a list of strings
        A list (for each ROI) of a list (for each contour) of a flat list of 
        [x, y, z] coordinates (as strings) of the polygons that define each 
        Target contour.
        
    TrgPtsByCntByRoi : list of a list of a list of a list of floats
        A list (for each ROI) of a list (for each contour) of a list (for each 
        point) of a list (for each dimension) of the polygons that define each 
        Target contour.
    
    TrgC2SindsByRoi : list of a list of integers
        A list (for each ROI) of a list (for each contour) of slice numbers 
        that correspond to each Target contour.
    
    TrgDicomDir : string
        Directory containing the Target DICOMs.
    
    AddTxtToSSLabel : string (optional; '' by default)
        Text to be added to StructureSetLabel for the new RTS.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ******
        
    NewTrgRts : Pydicom object
        New RTS object for the Target dataset.
    
    
    Notes:
    *****
    
    Performance of adding/removing sequences using append()/pop():
        
    * Took 329.9 ms to add 10000 sequences to ContourImageSequence 
    (33.0 us/sequence)
    * Took 165.4 ms to remove 10000 sequences to ContourImageSequence 
    (16.5 us/sequence)
    
    * Took 177.7 ms to add 10000 sequences to ContourSequence 
    (17.8 us/sequence)
    * Took 81.7 ms to remove 10000 sequences to ContourSequence 
    (8.2 us/sequence)
    
    19/02:
        After avoiding the use of deepcopy to increase sequences I've had to
        re-instate them since the tag values were being duplicated despite
        calls to modify them.  For Test RR3, creating the new RTS with 45 
        contours on 45 slices previously took 2 s to execute without using 
        deepcopy.  Using deepcopy it takes 11 s.
        
        This has also lengthened the time it takes to error check the RTS. The
        same RTS object formerly took 9 s to error check.  Now it takes 31 s.
        It's not clear why error checking should take longer.
        
        Test RR4 previously took 1 s to create the RTS and 5 s to error check
        it.  Now it takes 115 s and 343 s.
    """
    
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #import time
    import datetime
    from copy import deepcopy
    from DicomTools import ImportDicoms#, GetRoiNums, GetRoiLabels
    from GeneralTools import UniqueItems
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of CreateRts():')
        print(f'   TrgC2SindsByRoi = {TrgC2SindsByRoi}')
        print(f'   len(TrgPtsByCntByRoi) = {len(TrgPtsByCntByRoi)}')
        for r in range(len(TrgPtsByCntByRoi)):
            print(f'   len(TrgPtsByCntByRoi[{r}]) = {len(TrgPtsByCntByRoi[r])}')
        print(f'   len(TrgCntDataByCntByRoi) = {len(TrgCntDataByCntByRoi)}')
        for r in range(len(TrgCntDataByCntByRoi)):
            print(f'   len(TrgCntDataByCntByRoi[{r}] = {len(TrgCntDataByCntByRoi[r])}')
    
    SrcRts = dcmread(SrcRtsFpath)
    
    #print(f'\n\n\nTrgRtsFpath = {TrgRtsFpath}')
    
    """ Use TrgRts or SrcRts as a template for NewTrgRts. """
    if TrgRtsFpath:
        NewTrgRts = dcmread(TrgRtsFpath)
    else:
        NewTrgRts = deepcopy(SrcRts)
    
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    UniqueTrgC2Sinds = UniqueItems(Items=TrgC2SindsByRoi, IgnoreZero=False, 
                                   MaintainOrder=True)
    
    
    """ Generate a new SOPInstanceUID. """
    NewSOPuid = generate_uid()
    NewTrgRts.SOPInstanceUID = deepcopy(NewSOPuid)
    NewTrgRts.file_meta.MediaStorageSOPInstanceUID = deepcopy(NewSOPuid)
    
    """ Generate a new SeriesInstanceUID. """
    NewTrgRts.SeriesInstanceUID = generate_uid()
    
    #NewDate = time.strftime("%Y%m%d", time.gmtime())
    #NewTime = time.strftime("%H%M%S", time.gmtime())
    
    TimeNow = datetime.datetime.now()
    NewDate = TimeNow.strftime('%Y%m%d')
    NewTime = TimeNow.strftime('%H%M%S.%f')
    
    """ If TrgRtsFpath != None, some tags will not need to be replaced.
    If TrgRtsFpath = None, use the corresponding values in the first Target 
    DICOM. """
    
    if TrgRtsFpath == None:
        NewTrgRts.StudyDate = TrgDicoms[0].StudyDate
        NewTrgRts.StudyTime = TrgDicoms[0].StudyTime
        NewTrgRts.Manufacturer = TrgDicoms[0].Manufacturer
        NewTrgRts.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
        #NewTrgRts.SeriesDate = Dicoms[0].SeriesDate
        NewTrgRts.PatientName = TrgDicoms[0].PatientName
        NewTrgRts.PatientID = TrgDicoms[0].PatientID
        NewTrgRts.PatientBirthDate = TrgDicoms[0].PatientBirthDate
        NewTrgRts.PatientSex = TrgDicoms[0].PatientSex
        NewTrgRts.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
        NewTrgRts.StudyID = TrgDicoms[0].StudyID
        NewTrgRts.SeriesNumber = TrgDicoms[0].SeriesNumber
        NewTrgRts.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
        NewTrgRts.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
        
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
          
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .ReferencedSOPInstanceUID = TrgDicoms[0].StudyInstanceUID
        
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
       

    NewTrgRts.StructureSetLabel += AddTxtToSSLabel
    """ Modify StructureSetDate and StructureSetTime to the present: """
    NewTrgRts.StructureSetDate = NewDate
    NewTrgRts.StructureSetTime = NewTime
        
        
    # Start timing:
    #times = []
    #times.append(time.time())
    
    """ Modify ContourImageSequence. """
    
    """ Loop through each index in UniqueC2Sinds: """
    for i in range(len(UniqueTrgC2Sinds)):
        """ The DICOM slice number: """
        s = UniqueTrgC2Sinds[i]
        
        #print(f'\nUniqueC2Sinds[{i}] = {UniqueC2Sinds[i]}')
        #print(f'SOP UID is {Dicoms[s].SOPInstanceUID}')
        
        N = len(NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                         .RTReferencedStudySequence[0]\
                         .RTReferencedSeriesSequence[0]\
                         .ContourImageSequence)
        
        if i > N - 1:
            """ Increase the sequence by one. """
            #NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
            #         .RTReferencedStudySequence[0]\
            #         .RTReferencedSeriesSequence[0]\
            #         .ContourImageSequence.append(NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
            #                                               .RTReferencedStudySequence[0]\
            #                                               .RTReferencedSeriesSequence[0]\
            #                                               .ContourImageSequence[-1])
                     
            LastItem = deepcopy(NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                                         .RTReferencedStudySequence[0]\
                                         .RTReferencedSeriesSequence[0]\
                                         .ContourImageSequence[-1])
            
            NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                     .RTReferencedStudySequence[0]\
                     .RTReferencedSeriesSequence[0]\
                     .ContourImageSequence.append(LastItem)
        
        """ Update the sequence: """
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPClassUID = TrgDicoms[s].SOPClassUID
           
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
           
    
    """ Check if there are more sequences in ContourImageSequence than 
    required: """
    N = len(NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                     .RTReferencedStudySequence[0]\
                     .RTReferencedSeriesSequence[0]\
                     .ContourImageSequence)
    
    if N > len(UniqueTrgC2Sinds):
        for i in range(N - len(UniqueTrgC2Sinds)):
            """ Remove the last sequence: """
            NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                     .RTReferencedStudySequence[0]\
                     .RTReferencedSeriesSequence[0]\
                     .ContourImageSequence.pop()
    
            
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to modify 'ContourImageSequence.')
    
    
    """ Modify StructureSetROISequence. """
    
    """ Loop through each TrgC2SindsByRoi: """
    for r in range(len(TrgC2SindsByRoi)):
        #""" The number of contours in this ROI: """
        #C = len(TrgC2SindsByRoi[r])
        
        N = len(NewTrgRts.StructureSetROISequence)
        
        if r > N - 1:
            """ Increase the sequence by one. """
            
            LastItem = deepcopy(NewTrgRts.StructureSetROISequence[-1])
            
            NewTrgRts.StructureSetROISequence.append(LastItem)
    
    
        """ Update the sequence: """
        NewTrgRts.StructureSetROISequence[r].ROINumber = f"{r+1}"
        
        NewTrgRts.StructureSetROISequence[r]\
                 .ReferencedFrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
           
        NewTrgRts.StructureSetROISequence[r]\
                 .ROIName = NewTrgRts.StructureSetROISequence[r]\
                                     .ROIName + AddTxtToSSLabel
    
    
    """ Check if there are more sequences in StructureSetROISequence than 
    required: """
    N = len(NewTrgRts.StructureSetROISequence)
    
    if N > len(TrgC2SindsByRoi):
        for i in range(N - len(TrgC2SindsByRoi)):
            """ Remove the last sequence: """
            NewTrgRts.StructureSetROISequence.pop()
               
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to modify 'StructureSetROISequence.')
        
        
    """ 
    Modify ROIContourSequence. 
    """
    
    #print(f'\nTrgC2SindsByRoi = {TrgC2SindsByRoi}')
    
    
    """ Loop through each ROI in TrgC2SindsByRoi: """
    for r in range(len(TrgC2SindsByRoi)):
        
        n = 0 # total contour counter for this ROI
        
        """ The number of ROIContourSequences: """
        N = len(NewTrgRts.ROIContourSequence)
        
        if r > N - 1:
            """ Increase the sequence by one. """
            LastItem = deepcopy(NewTrgRts.ROIContourSequence[-1])
            
            NewTrgRts.ROIContourSequence.append(LastItem)
        
        """ Loop through each contour in TrgC2SindsByRoi[r]: """
        for c in range(len(TrgC2SindsByRoi[r])):
            #print(f'\n\n\nn = {n}')
            
            """ The DICOM slice number: """
            s = TrgC2SindsByRoi[r][c]
            
            """ The number of ContourSequences: """
            N = len(NewTrgRts.ROIContourSequence[r].ContourSequence)
        
            if n > N - 1:
                """ Increase the sequence by one. """
                #NewTrgRts.ROIContourSequence[r]\
                #         .ContourSequence.append(NewTrgRts.ROIContourSequence[r]\
                #                                          .ContourSequence[-1])
                         
                LastItem = deepcopy(NewTrgRts.ROIContourSequence[r]\
                                             .ContourSequence[-1])
                
                NewTrgRts.ROIContourSequence[r]\
                         .ContourSequence.append(LastItem)
                   
        
            #print(f'\nRts.ROIContourSequence[{r}].ContourSequence[{i}].ContourImageSequence[0].ReferencedSOPInstanceUID =',
            #      f'{Rts.ROIContourSequence[r].ContourSequence[i].ContourImageSequence[0].ReferencedSOPInstanceUID}')
            #print(f'\nDicoms[{s}].SOPInstanceUID = {Dicoms[s].SOPInstanceUID}')
            
            """ Update the sequence: """
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[n]\
                     .ContourImageSequence[0]\
                     .ReferencedSOPClassUID = TrgDicoms[s].SOPClassUID
               
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[n]\
                     .ContourImageSequence[0]\
                     .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
            
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[n]\
                     .NumberOfContourPoints = f"{len(TrgPtsByCntByRoi[r][c])}"
             
            #print(f'\nn = {n}')
            
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[n]\
                     .ContourNumber = f"{n+1}"
                     
            #print(f'NewTrgRts.ROIContourSequence[{r}].ContourSequence[{n}].ContourNumber =',
            #      f'{NewTrgRts.ROIContourSequence[r].ContourSequence[n].ContourNumber}')
               
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[n]\
                     .ContourData = TrgCntDataByCntByRoi[r][c]
            
            NewTrgRts.ROIContourSequence[r]\
                     .ReferencedROINumber = f"{r+1}"
            
            n += 1 # increment the total contour counter
        
        
    
        """ Check if there are more sequences in ContourSequence than 
        required: """
        N = len(NewTrgRts.ROIContourSequence[r].ContourSequence)
        
        if N > len(TrgC2SindsByRoi[r]):
            #print(f'There are {N} sequences in ROIContourSequence[{r}].ContourSequence',
            #      f'but only {len(TrgC2SindsByRoi[r])} contours in this ROI.')
            for i in range(N - len(TrgC2SindsByRoi[r])):
                """ Remove the last sequence: """
                NewTrgRts.ROIContourSequence[r].ContourSequence.pop()
    
    
    #print('\nCheck ContourNumbers:')
    #for r in range(len(NewTrgRts.ROIContourSequence)):
    #    for n in range(len(NewTrgRts.ROIContourSequence[r].ContourSequence)):
    #        print(f'\nn = {n}')
    #                 
    #        print(f'NewTrgRts.ROIContourSequence[{r}].ContourSequence[{n}].ContourNumber =',
    #              f'{NewTrgRts.ROIContourSequence[r].ContourSequence[n].ContourNumber}')
    
    """ Check if there are more sequences in ROIContourSequence than 
    required: """
    N = len(NewTrgRts.ROIContourSequence)
    
    if N > len(TrgC2SindsByRoi):
        for i in range(N - len(TrgC2SindsByRoi)):
            """ Remove the last sequence: """
            NewTrgRts.ROIContourSequence.pop()
        
        
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to modify 'ROIContourSequence.')
    
        
        
    """ 
    Modify RTROIObservationsSequence. 
    """
    
    """ Loop through each TrgC2SindsByRoi: """
    for r in range(len(TrgC2SindsByRoi)):
        #""" The number of contours in this ROI: """
        #C = len(TrgC2SindsByRoi[r])
        
        N = len(NewTrgRts.RTROIObservationsSequence)
        
        if r > N - 1:
            """ Increase the sequence by one. """
            LastItem = deepcopy(NewTrgRts.RTROIObservationsSequence[-1])
            
            NewTrgRts.RTROIObservationsSequence.append(LastItem)
    
        """ Update the sequence: """
        NewTrgRts.RTROIObservationsSequence[r].ObservationNumber = f"{r+1}"
        
        NewTrgRts.RTROIObservationsSequence[r].ReferencedROINumber = f"{r+1}"
           
    
    """ Check if there are more sequences in RTROIObservationsSequence than 
    required: """
    N = len(NewTrgRts.RTROIObservationsSequence)
    
    if N > len(TrgC2SindsByRoi):
        for i in range(N - len(TrgC2SindsByRoi)):
            """ Remove the last sequence: """
            NewTrgRts.RTROIObservationsSequence.pop()
            
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if True:#LogToConsole:
    #    print(f'\nTook {Dtime} s to add/remove sequences from',
    #          'RTROIObservationsSequence.')
    
    if LogToConsole:
        print('-'*120)
    
    return NewTrgRts








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
    
    
    Notes:
    *****
    
    19/02:
        After avoiding the use of deepcopy to increase sequences I've had to
        re-instate them since the tag values were being duplicated despite
        calls to modify them.  For Test RR3, creating the new RTS with 45 
        contours on 45 slices previously took 2 s to execute without using 
        deepcopy.  Using deepcopy it takes 11 s.
        
        This has also lengthened the time it takes to error check the RTS. The
        same RTS object formerly took 9 s to error check.  Now it takes 31 s.
        It's not clear why error checking should take longer.
        
        Test RR4 previously took 1 s to create the RTS and 5 s to error check
        it.  Now it takes 115 s and 343 s.
    """        
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids, ImportDicom
    from ImageTools import GetImageAttributes
    #from GeneralTools import PrintTitle
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    Dicom = ImportDicom(DicomDir)
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    
    LogList = []
    Nerrors = 0
    
    if LogToConsole:
        print('\n\n', '-'*120)
        #PrintTitle('Error checking RTS..')
        print('Results of running ErrorCheckRts():')
    
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
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
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
            
            
            
    #"""Verify that the StructureSetROISequence has only 1 sequence."""
    #N = len(Rts.StructureSetROISequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in StructureSetROISequence.'\
    #          + ' There should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in StructureSetROISequence.'\
    #          + ' There should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)
            
    
    
    #"""Verify that the ROINumber in StructureSetROISequence is 1."""
    #N = int(deepcopy(Rts.StructureSetROISequence[0].ROINumber))
    #
    #if N == 1:
    #    msg = f'INFO:  ROINumber in StructureSetROISequence is {N}.  It '\
    #          + 'should be 1.\n'
    #else:
    #    msg = f'ERROR:  ROINumber in StructureSetROISequence is {N}.  It '\
    #          + 'should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)
    
    
    
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
            
            
    
    #"""Verify that the ROIContourSequence has only 1 sequence."""
    #N = len(Rts.ROIContourSequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in ROIContourSequence. There '\
    #          + 'should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in ROIContourSequence. There '\
    #          + 'should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)

    
    
    
    
    """Verify that the ReferencedSOPInstanceUIDs in ROIContourSequence
    match the ReferencedSOPInstanceUIDs in ContourImageSequence."""
    CS = deepcopy(Rts.ROIContourSequence[0].ContourSequence)
    
    RefSOPuidsInRCS = [CS[i].ContourImageSequence[0].ReferencedSOPInstanceUID for i in range(len(CS))]
    
    # Determine whether the Ref SOP UIDs match the SOP UIDs:
    IsMatch = [uid in RefSOPuidsInRFORS for uid in RefSOPuidsInRCS]
    
    #NumOfMatches = IsMatch.count(True)
    
    # Find the indices of any non-matching Ref SOP UIDs:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
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
    
    #print(f'\n\n\nContourNums = {ContourNums}')
    
    # Determine whether the contour numbers match the expected contour numbers:
    IsMatch = [ContourNums[i] == ExpectedNums[i] for i in range(len(CS))]
    #IsMatch = [ContourNums[i] == i+1 for i in range(len(CS))]
    
    # Find the indices of any non-matching contour numbers:
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
    if Inds:        
        for i in range(len(Inds)):
            cn = ContourNums[Inds[i]]
            en = ExpectedNums[Inds[i]]
            
            msg = f'ERROR:  ContourNumber {i+1} is {cn} but is expected to be'\
                  + f' {en}.\n'
            LogList.append(msg)      
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
    Inds = [i for i, x in enumerate(IsMatch) if x==False]
    
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
    
    
    
    #"""Verify that the ReferencedROINumber in ROIContourSequence is 1."""
    #N = int(deepcopy(Rts.ROIContourSequence[0].ReferencedROINumber))
    #
    #if N == 1:
    #    msg = f'INFO:  ReferencedROINumber in ROIContourSequence is {N}. '\
    #          + 'It should be 1.\n'
    #else:
    #    msg = f'ERROR:  ReferencedROINumber in ROIContourSequence is {N}. '\
    #          + 'It should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)
    
    
    
    #"""Verify that the ROIObservationsSequence has only 1 sequence."""
    #N = len(Rts.RTROIObservationsSequence)
    #
    #if N == 1:
    #    msg = f'INFO:  There are {N} sequences in RTROIObservationsSequence.'\
    #          + ' There should be 1.\n'
    #else:
    #    msg = f'ERROR:  There are {N} sequences in RTROIObservationsSequence.'\
    #          + ' There should be 1.\n'
    #          
    #    Nerrors = Nerrors + 1
    
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)


    
    #"""Verify that the ObservationNumber is 1."""
    #N = int(deepcopy(Rts.RTROIObservationsSequence[0].ObservationNumber))
    #
    #if N == 1:
    #    msg = 'INFO:  ObservationNumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #else:
    #    msg = 'ERROR:  ObservationNumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)
            
            
            
    #"""Verify that the ReferencedROINumber is 1."""
    #N = int(deepcopy(Rts.RTROIObservationsSequence[0].ReferencedROINumber))
    #
    #if N == 1:
    #    msg = 'INFO:  ReferencedROINumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #else:
    #    msg = 'ERROR:  ReferencedROINumber in RTROIObservationsSequence is '\
    #          + f'{N}. It should be 1.\n'
    #    
    #    Nerrors = Nerrors + 1
    #
    ##msg = msg + f'Nerrors = {Nerrors}.\n'
    #
    #LogList.append(msg)
    #
    #if LogToConsole:
    #    print(msg)
    
    
    #print(f'There were {Nerrors} errors found in the RTS.\n')
    
    if LogToConsole:
        print('-'*120)
    
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








def GetRtsDataFromListOfRtss(ListOfRtss, ListOfDicomDirs, LogToConsole=False):
    """ 
    Note:
    
    ListOfRtss is a list (which could be of length 1) of RTS (Pydicom) objects.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each RTS.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #from DicomTools import GetRoiNum
    from DicomTools import GetDicomFpaths
    from ImageTools import GetImageAttributes
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of GetRtsDataFromListOfRtss():')
    
    ListOfPtsByCntByRoi = []
    ListOfC2SindsByRoi = []
    #ListOfPtsByCntBySlice = []
    #ListOfRoiNums = []
    ListOfDcmFpaths = []
    #ListOfOrigins = []
    ListOfIPPs = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfRtss)):
        if LogToConsole:
            print(f'\nListOfRtss[{i}]:')
        
        if ListOfRtss[i]:
            
            PtsByCntByRoi, C2SindsByRoi = GetPtsByCntByRoi(ListOfRtss[i], 
                                                           ListOfDicomDirs[i], 
                                                           LogToConsole)
            
            if LogToConsole:      
                #print(f'\nlen(PtsByCntByRoi) in "{SearchString}" = ',
                #      f'{len(PtsByCntByRoi)}')
                print(f'C2SindsByRoi = {C2SindsByRoi}')
                [print(f'len(PtsByCntByRoi{i}) = {len(PtsByCntByRoi[i])}') for i in range(len(PtsByCntByRoi))]
                
        else:
            PtsByCntByRoi = None
            C2SindsByRoi = None
            
            if LogToConsole:      
                print(f'No contours for this dataset')
        
        
        #ListOfRoiNums.append(RoiNum)
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
        #ListOfPtsByCntBySlice.append(PtsByCntBySlice)
        ListOfPtsByCntByRoi.append(PtsByCntByRoi)
        ListOfC2SindsByRoi.append(C2SindsByRoi)
        ListOfDcmFpaths.append(DcmFpaths)
        #ListOfOrigins.append(IPPs[0])
        ListOfIPPs.append(IPPs)
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
        
        if LogToConsole:      
            #print(f'len(Contours) in "{SearchString}" = {len(Contours)}')
            print(f'\nlen(ListOfPtsByCntByRoi) = {len(ListOfPtsByCntByRoi)}')
            #print(f'C2Sinds = {C2Sinds}')
            #[print(f'len(Contours{i}) = {len(Contours[i])}') for i in range(len(Contours))]
            #print(f'\nlen(ListOfPtsByCntBySlice) = {len(ListOfPtsByCntBySlice)}')
        
        
        #SOPuids = GetDicomSOPuids(ListOfDicomDirs[i])
        
        #CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(ListOfRtss[i], SOPuids)
        
        #ListOfC2Sinds.append(CStoSliceIndsByRoi[RoiNum])
        
        #if LogToConsole:
        #    print(f'ListOfC2Sinds[{i}] = {ListOfC2Sinds[i]}')
    
    if LogToConsole:
        print('-'*120)
        
    return ListOfPtsByCntByRoi, ListOfC2SindsByRoi, ListOfDcmFpaths,\
           ListOfIPPs, ListOfDirections, ListOfSpacings
    #return ListOfPtsByCntByRoi, ListOfC2SindsByRoi, ListOfDcmFpaths





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






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / ROI / RTSTRUCT
******************************************************************************
******************************************************************************
"""

def CopyRts(SrcRtsFpath, SrcSliceNum, SrcRoiName, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgRtsFpath=None, TrgRoiName=None, TrgSliceNum=None, 
            ResInterp='BlurThenLinear', PreResVariance=(1,1,1),
            ApplyPostResBlur=False, PostResVariance=(1,1,1), 
            ForceReg=False, SelxOrSitk='Sitk', Transform='affine', 
            MaxIters='512', InitMethod='centerofgravity', 
            SrcFidsFpath='moving_fiducials.txt', 
            TrgFidsFpath='target_fiducials.txt', TxInterp='NearestNeighbor', 
            ApplyPostTxBin=True, ApplyPostTxBlur=True, PostTxVariance=(1,1,1),  
            TxMatrix=None, GridDims=None, GridRes=None, VectGridData=None,
            TxtToAddToTrgRoiName='', LogToConsole=False,
            DictOfInputs={}, ListOfInputs=[], ListOfTimings=[]):
    """
    Note 01/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific contour
            2. All contours in a specific ROI
            3. All contours in all ROIs
        
     
    Inputs:
    ******
    
    SrcRtsFpath : string
        Filepath of the Source RTS file.
    
    SrcSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If SrcSliceNum = None, a relationship-preserving copy will be made.
        
    SrcRoiName : string
        All or part of the ROIName of the ROI containing the contour(s) to be
        copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target RTS file that the contour(s) is/are to be copied 
        to. TrgRtsFpath != None only if an existing Target RTS exists and is to
        be added to. An existing RTS will be added to only for Direct copy 
        operations.
    
    TrgRoiName : string or None (optional; None by default) 
        All or part of the ROIName of the destination ROI.
        
    TrgSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the contour will be
        copied to (applies only for the case of direct copies of single 
        contours). The index is zero-indexed.
        If TrgSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the contour(s) will be copied to will
        not depend on user input.
    
    ResInterp : string
        The interpolator to be used for (non-registration-based) resampling of
        the Source labelmap image(s) to the Target grid (if applicable). 
        Acceptable values are:
        - 'NearestNeighbor'
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    PreResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to perform binary thresholding following
    #    resampling using a non-label (e.g. linear) interpolator. Example use is 
    #    on a labelmap image if a Gaussian blurring + linearly resampling 
    #    approach is taken to combat aliasing effects.
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
        
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        resampled labelmap image(s).
        
    ForceReg : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    SelxOrSitk : string (optional; 'Sitk' by default)
        Denotes which package to use for image registration and transformation.
        Acceptable values include:
            - 'Selx' for SimpleElastix
            - 'Sitk' for SimpleITK
     
    Transform : string (optional; 'affine' by default)
        The transformation to used for image registration (if applicable).  
        Acceptable values are:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxIters : string (optional; '512' by default)
        If 'default', the maximum number of iterations used for the optimiser
        during image registration (if applicable) will be the pre-set default
        in the parameter map for Transform. If != 'default' it must be a string 
        representation of an integer.
    
    InitMethod : string (optional, 'centerofgravity' by default)
        The initialisation method used for initial alignment prior to image
        registration (if applicable).
    
    SrcFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Source/Moving image.
    
    TrgFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The full filepath (or filename within the current working directory)
        of the list of fiducials for the Target/Fixed image.
    
    TxInterp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for registration-based resampling (i.e.
        transformation; if applicable).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image will be binary thresholded.
        
    #ThreshPostTx : float (optional; 0.05 by default)
    #    The threshold level used to perform binary thresholding after 
    #    registration transformation.  
    
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image will be Gaussian blurred.
        
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
    
    TxMatrix : list of float strings or None (optional; None by default)
        List of float strings representing the non-deformable transformation 
        that transforms the moving (Source) image to the fixed (Target) image;
        or None.  If not None, image registration will be skipped to save 
        computational time.
    
    GridDims : list of integers or None (optional; None by default)
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image; 
        or None.  If not None, image registration will be skipped to save 
        computational time.
    
    GridRes : list of floats or None (optional; None by default)
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image; or None.
        If not None, image registration will be skipped to save computational
        time.
    
    VectGridData : list of float strings or None (optional; None by default)
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image; or None.
        If not None, image registration will be skipped to save computational
        time.
    
    TxtToAddToTrgRoiName : string (optional, '' by default)
        String of text to add to the ROI Name of the new Target RTS.
                           
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary (optional; empty by default)
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings (optional; empty by default)
        A list containing the inputs that were called to CopyRoi().
        
    LogOfTimings : list of strings (optional; empty by default)
        A list of the time to execute certain tasks during the calling of 
        CopyRts().
          
                  
    Outputs:
    *******
        
    TrgRts : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target RTS object.
    
    Dro : Pydicom object or None
        DICOM registration object if image registration was used to create 
        TrgRts, None otherwise.
        
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    ListOfTimings : list of strings
        A list (for each message) of timings for individual operations.
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours will appear displaced w.r.t. 
    the anatomical features highlighted in the Source image.
    """
    
    import importlib
    #import RtsTools
    #importlib.reload(RtsTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    import DroTools
    importlib.reload(DroTools)
    import RoiCopyTools
    importlib.reload(RoiCopyTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    ##from RtsTools import GetPtsByCntByRoi
    #from RtsTools import GetRtsDataOfInterest, ProportionOfRoisInExtent
    #from RtsTools import CreateRts
    from GeneralTools import UniqueItems, PrintIndsByRoi#, PrintTitle
    from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    from GeneralTools import ShiftPtsByCntByRoi
    #from GeneralTools import ShiftFrame
    #from DicomTools import GetRoiNum
    #from ConversionTools import Points2ContourData
    from ConversionTools import PtsByCntByRoi2CntDataByCntByRoi
    from ImageTools import ImportImage
    from DroTools import CreateSpaDro, CreateDefDro
    from RoiCopyTools import CheckValidityOfInputs
    
    """ Start timing. """
    times = []
    times.append(time.time())
    #""" Log timing messages. """
    #TimingMsgs = []
    
    
    SrcRts = dcmread(SrcRtsFpath)
    if TrgRtsFpath:
        TrgRts = dcmread(TrgRtsFpath)
    else:
        TrgRts = None
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print(f'Running CopyRts():')
        
        SrcRoiLabels = GetRoiLabels(SrcRts)
        
        print(f'\nSrcRoiLabels = {SrcRoiLabels}')
        
        if TrgRtsFpath:
            TrgRoiLabels = GetRoiLabels(TrgRts)
        
            print(f'\nTrgRoiLabels = {TrgRoiLabels}')
    
    
    """ Check the validity of the combination of inputs. """
    CheckValidityOfInputs(SrcRts, SrcRoiName, SrcSliceNum, TrgRts, TrgRoiName,
                          TrgSliceNum, LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(SrcSliceNum, SrcRoiName, TrgSliceNum, 
    #                              SrcDcmDir, TrgDcmDir, ForceRegistration,
    #                              LogToConsole=True)
    
    #""" Manipulate UseCase if ForceRegistration = True and UseCase was 3 or 4. """
    #if ForceRegistration and not '5' in UseCase:
    #    msg = f'\nUseCase was {UseCase} but will be treated as UseCase '
    #    
    #    UseCase.replace('3', '5').replace('4', '5')
    #    
    #    msg = f'{UseCase}.'
    #    
    #    print(msg)
    
    """ Get the data of interest from the Source RTS. """
    SrcPtsByCntByRoi,\
    SrcC2SindsByRoi = GetRtsDataOfInterest(SrcRtsFpath, SrcSliceNum, 
                                           SrcRoiName, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source RTS file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Raise exception if there is no data of interest from the Source RTS."""
    if not UniqueItems(SrcC2SindsByRoi):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        #from RtsTools import GetNumOfContoursByRoi, GetCStoSliceIndsByRoi
        
        SrcRts = dcmread(SrcRtsFpath)
        
        Names = GetRoiLabels(SrcRts)
        N = len(Names)
        
        NumByRoi = GetNumOfContoursByRoi(SrcRts)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        C2SindsByRoi = GetCStoSliceIndsByRoi(SrcRts, SOPuids)
        
        msg = f"There are no contours on slice {SrcSliceNum} in any ROI in "\
              + f"the Source RTS. There are {N} ROIs in the Source RTS:"
        
        for i in range(N):
            msg += f"\nROI {i+1} with name '{Names[i]}' has {NumByRoi[i]} "\
                   + f"contours that correspond to slices {C2SindsByRoi[i]}" 
        
        raise Exception(msg)
    
    """
    Determine whether the contours that make up the ROI(s) intersect with the
    Target image extent (14/01/2021). 
    """
    FracProp = ProportionOfRoisInExtent(PtsByCntByRoi=SrcPtsByCntByRoi,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source ROI(s) intersect the '\
          + 'Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
        
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the points in the ROI lie ',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the points in the ROI lie within the physical extent '\
              + 'of the Target image.'
        raise Exception(msg)
        return None
    
    
    """ Import the 3D images. """
    SrcIm = ImportImage(SrcDcmDir)
    TrgIm = ImportImage(TrgDcmDir)
    
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if UseCase in ['1', '2a', '3a']: # 08/02
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a contour without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        
        """ The contour on SrcSliceNum is to be copied to TrgSliceNum. Modify 
        the contour-to-slice index (from SrcSliceNum) to TrgSliceNum. """
        SrcC2SindsByRoi = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcC2SindsByRoi,
                                                   IndToReplace=SrcSliceNum, 
                                                   ReplacementInd=TrgSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
                
        if TrgRtsFpath:
            """ Direct copy of a contour to an existing RTS. """
            
            """ An existing Target RTS is to be modified.  Since this is a  
            Direct copy, any existing contours in Target with ROI label  
            matching SrcRoiName are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPtsByCntByRoi,\
            TrgC2SindsByRoi = GetRtsDataOfInterest(TrgRtsFpath, TrgSliceNum,
                                                   #SrcRoiName, TrgDcmDir,
                                                   TrgRoiName, TrgDcmDir,
                                                   LogToConsole)
            """ 05/03/21: Previously SrcRoiName was used as an input in 
            GetRtsDataOfInterest for Target, hence requiring that the Target
            ROI have the same name as the Source ROI. Now TrgRoiName has been
            added to the list of inputs allowing for the Target ROI to have a
            different name. """
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target RTS file for '\
                  + 'existing contours within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            
            if LogToConsole:
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}.')
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The contour will be copied to a different slice location, so the
        z-components of the coordinates will need to be modified. """
        
        """ Shift the out-of-plane (z) components of the points in
        SrcPtsByCntByRoi. """
        SrcPtsByCntByRoi = ZshiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi,
                                               NewZind=TrgSliceNum,
                                               DicomDir=SrcDcmDir)
        
        """ Shift the points to account for any difference in the origin of
        SrcIm and TrgIm, and to apply the necessary shift from SrcSliceNum to
        TrgSliceNum. """
        SrcPtsByCntByRoi,\
        SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                             C2SindsByRoi=SrcC2SindsByRoi, 
                                             SrcImage=SrcIm, 
                                             SrcSliceNum=SrcSliceNum, 
                                             TrgImage=TrgIm, 
                                             TrgSliceNum=TrgSliceNum, 
                                             RefImage=TrgIm, ShiftInX=False,  
                                             ShiftInY=False, ShiftInZ=True, 
                                             Fractional=False, 
                                             LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
        
        
        if TrgRtsFpath:
            """ An existing Target RTS is to be modified. 
            Concatenate TrgC2SindsByRoi and SrcC2SindsByRoi, and 
            TrgPtsByCntByRoi and SrcPtsByCntByRoi. """
            
            TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + SrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
            
            TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + SrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
        
        else:
            """ A Target RTS was not provided. """
            
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
            if LogToConsole:
                print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
           
        
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi, and 
            TrgC2SindsByRoi is SrcC2SindsByRoi. """
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
        
        
        
        
        if False:
            """ Shift the points to account for any difference in the origin (i.e.
            slice number 0) of SrcIm and TrgIm, and to apply the necessary shift in 
            the z components of the points, and to the indices in SrcC2SindsByRoi. """
            SrcPtsByCntByRoi,\
            SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                                 C2SindsByRoi=SrcC2SindsByRoi, 
                                                 SrcImage=SrcIm, 
                                                 SrcSliceNum=0, 
                                                 TrgImage=TrgIm, 
                                                 TrgSliceNum=0, 
                                                 RefImage=TrgIm, ShiftInX=False,  
                                                 ShiftInY=False, ShiftInZ=True, 
                                                 Fractional=False, 
                                                 LogToConsole=LogToConsole)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        #TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1(SrcC2SindsByRoi, 
        #                                                        SrcIm, TrgIm)
        
        TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcC2SindsByRoi, 
                                                                SrcIm, TrgIm)
        
        """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi. """
        TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        
        if LogToConsole:
            print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
            print(f'TrgC2SindsByRoi (= shifted SrcC2SindsByRoi) =',
                  f'{SrcC2SindsByRoi}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        from ImageTools import GetImageInfo, ResampleLabImByRoi
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PtsByCntByRoi2PixArrByRoi
        from ConversionTools import PixArrByRoi2LabImByRoi
        from ConversionTools import PixArrByRoi2PtsByCntByRoi
        #from RtsTools import GetMaskFromContoursForSliceNum
        #from SegTools import ProportionOfSegsInExtent
        #from GeneralTools import NumOfListsAtDepthTwo
        
        
        #""" Import the 3D images. """
        #SrcIm = ImportImage(SrcDcmDir)
        #TrgIm = ImportImage(TrgDcmDir)
        
        """ Convert the RTS data of interest to a list of 3D pixel arrays for 
        each ROI in SrcPtsByCntByRoi. """
        SrcPixArrByRoi = PtsByCntByRoi2PixArrByRoi(SrcPtsByCntByRoi, SrcIm, 
                                                   LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source points to pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Determine whether the pixels that make up SrcPixArrByRoi intersect 
        with the Target image extent. 
        
        Note:
            This is redundant since this was already checked for the points.  
            Only useful as a confirmation that the conversion from points to 
            pixel arrays was correct.
        
        FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrByRoi, 
                                            F2SindsBySeg=SrcC2SindsByRoi,
                                            SrcDicomDir=SrcDcmDir,
                                            TrgDicomDir=TrgDcmDir,
                                            LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'\n{round(FracProp*100, 1)}% of the voxels in the',
                  'segmentation to be copied lie within the physical extent',
                  'of the Target image.')
        
        if FracProp == 0:
            msg = 'None of the voxels in the segmentation to be copied lie '\
                  + 'within the physical extent of the Target image.'
            raise Exception(msg)
            return None
        """
        
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabImByRoi = PixArrByRoi2LabImByRoi(PixArrByRoi=SrcPixArrByRoi, 
                                               F2SindsByRoi=SrcC2SindsByRoi, 
                                               RefIm=SrcIm,
                                               LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              + 'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the frame-to-slice indices in the order associated with
        SrcLabImByRoi.
        Comment:
            Not sure the potential ordering difference between SrcC2SindsByRoi
            and SrcF2SindsByRoi is important for further steps.. """
        SrcF2SindsByRoi = []

        for SrcLabIm in SrcLabImByRoi:
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(SrcLabIm, LogToConsole=False)
            
            SrcF2SindsByRoi.append(F2Sinds)
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        
        """ Resample the source labelmaps.
        Comment:
            Should use SrcF2SindsByRoi rather than SrcC2SindsByRoi?..."""
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = ResampleLabImByRoi(LabImByRoi=SrcLabImByRoi, 
                             #F2SindsByRoi=SrcC2SindsByRoi,
                             F2SindsByRoi=SrcF2SindsByRoi,
                             SrcIm=SrcIm, 
                             TrgIm=TrgIm,
                             Interp=ResInterp,
                             PreResVariance=PreResVariance,
                             ApplyPostResBlur=ApplyPostResBlur,
                             PostResVariance=PostResVariance,
                             LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the interpolation that was used (the actual interpolation used 
        might have been different from the interpolation set) and the threshold
        used to re-binarise the resampled labelmap image (if applicable): """
        
        #print(f'\n\n\nMetadata keys in ResSrcLabImByRoi[0]:',
        #      ResSrcLabImByRoi[0].GetMetaDataKeys())
        
        for i in range(len(ResSrcLabImByRoi)):
            Interp = ResSrcLabImByRoi[i].GetMetaData("ResInterpUsed")
            
            ListOfInputs.append(f'ResInterpUsedForRoi{i} = {Interp}')
            DictOfInputs[f'ResInterpUsedForRoi{i}'] = Interp
            
            if Interp == 'BlurThenLinear':
                Thresh = ResSrcLabImByRoi[i].GetMetaData("PostResThreshUsed")
                
                ListOfInputs.append(f'PostResThreshUsedForRoi{i} = {Thresh}')
                DictOfInputs[f'PostResThreshUsedForRoi{i}'] = Thresh
    
        
        
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from ImageTools import RegisterImagesSelx, RegisterImagesSitk
        from ImageTools import TransformLabImByRoi, CreateSitkTxFromTxMatrix
        from ImageTools import GetParamFromSelxImFilt
        
        if not Transform in ['rigid', 'affine', 'bspline']:
            msg = f'The chosen transform, {Transform}, is not one of the '\
                  + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
            
            raise Exception(msg)
        
        if (Transform in ['rigid', 'affine'] and TxMatrix == None) \
        or (Transform == 'bspline' and VectGridData == None):
            """ The transform parameters required to register SrcIm to TrgIm 
            were not found from any matching subject assessors in XNAT. Perform
            image registration: """
            
            if LogToConsole == False:
                print(f'Performing image registration using {Transform}',
                      'transform...\n')
            
            times.append(time.time())
            
            """ Register SrcIm to TrgIm. """
            if SelxOrSitk == 'Selx':
                RegIm, SelxImFiltOrSitkTx, TxParamMapDict\
                = RegisterImagesSelx(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform,
                                     MaxNumOfIters=MaxIters,
                                     InitMethod=InitMethod,
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=SrcFidsFpath,
                                     FlipK=False, LogToConsole=LogToConsole)
                
                #""" Get the registration transform parameters: """
                #TxParams = list(TxParamMapDict['TransformParameters'])
                
                """ Get the registration transform parameters for the final
                transform parameter map: (10/05/21) """
                TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                                  Param='TransformParameters')
                
                
            if SelxOrSitk == 'Sitk':
                RegIm, SelxImFiltOrSitkTx\
                = RegisterImagesSitk(FixIm=TrgIm, MovIm=SrcIm, 
                                     Transform=Transform, InitMethod=InitMethod,
                                     FinalInterp='Linear',
                                     FixFidsFpath=TrgFidsFpath, 
                                     MovFidsFpath=TrgFidsFpath,
                                     FlipK=False, LogToConsole=LogToConsole)
                
                """ Get the registration transform parameters: """
                TxParams = [str(item) for item in SelxImFiltOrSitkTx.GetParameters()]
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to register the Source image to the Target '\
                  + f'image using {Transform} transform.\n'
            ListOfTimings.append(msg)
            if LogToConsole == False:
                print(f'*{msg}')
            
            #print(f'TxParams = {TxParams}')
            
            #""" Add the row vector [0, 0, 0, 1] to TxParams to complete the 
            #Frame of Reference Transformation Matrix 
            #(http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
            #"""
            ##TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
            #TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
            #            str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
            #            str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
            #            str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
            
            
        else:
            """ The transform parameters required to register SrcIm to TrgIm 
            are obtained from the transform matrix, found from the DRO (subject 
            assessor) in XNAT. Create a SimpleITK transform from TxMatrix: 
            
            Note: 'bspline' transform is not currently catered for.
            """
            
            if Transform in ['rigid' 'affine']:
                SelxImFiltOrSitkTx = CreateSitkTxFromTxMatrix(TrgIm, TxMatrix, 
                                                              Transform)
                
                """ TxMatrix has additional bottom row [0 0 0 1] (see 
                http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
                
                Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. 
                Hence:
                TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
                """
                TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 06/05/21
                
                print(f'TxParams from TxMatrix (from DRO) = {TxParams}\n')
                
                print(f'\n')
            else:
                """ 10/05/21: Still need to create CreateSelxImFiltFromParams()
                """
                SelxImFiltOrSitkTx = CreateSelxImFiltFromParams(TrgIm, GridDims, 
                                                                GridRes, 
                                                                VectGridData,
                                                                Transform)
            
        
        
        """ Transform SrcLabImByRoi using the registration transformation. 
        
        Note:
            Even though the data is being transformed the prefix 'Res' will be
            used since further operations below will need to be applied to
            either resampled or transformed data. 
            
            Note: 'bspline' transform is not currently catered for.
            """
        
        """ 10/05/21: Need to test whether TransformLabImByRoi() works with
        SelxImFilt passed by yet-to-be-developed GetSelxImFilterFromParams() """
        
        ResSrcLabImByRoi, ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = TransformLabImByRoi(LabImByRoi=SrcLabImByRoi,
                              F2SindsByRoi=SrcF2SindsByRoi,
                              TrgIm=TrgIm,
                              SelxImFiltOrSitkTx=SelxImFiltOrSitkTx,
                              Interp=TxInterp,
                              ApplyPostTxBlur=ApplyPostTxBlur,
                              PostTxVariance=PostTxVariance,
                              ApplyPostTxBin=ApplyPostTxBin,
                              LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
        """ Get the threshold used to re-binarise the transformed labelmap 
        image: """
        for i in range(len(ResSrcLabImByRoi)):
            Thresh = ResSrcLabImByRoi[i].GetMetaData("PostTxThreshUsed")
            
            ListOfInputs.append(f'PostTxThreshUsedForRoi{i} = {Thresh}')
            DictOfInputs[f'PostTxThreshUsedForRoi{i}'] = Thresh
    else:
        TxParams = None
        
        
    
    if UseCaseToApply in ['3a', '4a', '5a']:
        """ Direct copy with resampling/registration and averaging of multi-
        framed pixel arrays.
        
        Note:
        
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled/registered labelmap image).  So if after resampling there are 
        more than one frame, the frames will be averaged or a logical OR 
        operation will be performed.
        """
        
        from GeneralTools import ShiftFramesInPixArrBySeg
        from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrByRoi) # number of ROIs
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            """ Perform averaging or OR operation? """
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                f'Prior to {ReduceUsing} operation')
            
            
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5
    
                ResSrcPixArrByRoi, ResF2SindsByRoi\
                = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                         F2SindsBySeg=ResSrcF2SindsByRoi,
                                         MakeBinary=True, 
                                         BinaryThresh=BinaryThresh,
                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                #print(f'\n\n\nBefore running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
                    
                ResSrcPixArrByRoi, ResF2SindsByRoi\
                = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                       F2SindsBySeg=ResSrcF2SindsByRoi,
                                       LogToConsole=LogToConsole)
                
                #print(f'\nAfter running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
            
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment.')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                'After {ReduceUsing} operation')
        
        
        
        """ Shift the out-of-plane elements in ResSrcPixArrByRoi to account for
        the shift from SrcSliceNumber to TrgSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be SrcSliceNum but instead
        ResSrcF2SindsByRoi[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ SrcSliceNum --> ResSrcSliceNum following resampling or 
        transformation. """
        ResSrcSliceNum = ResSrcF2SindsByRoi[0][0]
        
        if True:#LogToConsole:
            #print(f'\nSrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
            #      f'{ResSrcSliceNum}')
            #print(f'TrgSliceNum = {TrgSliceNum}')
            print(f'SrcSliceNum = {SrcSliceNum} --> ResSrcSliceNum =',
                  f'{ResSrcSliceNum} following resampling/transformation -->',
                  f'TrgSliceNum = {TrgSliceNum} for direct copy\n')
            print('ResSrcF2SindsByRoi prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsByRoi}\n')
        
        ResSrcPixArrByRoi, ResSrcF2SindsByRoi\
        = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi, 
                                   F2SindsBySeg=ResSrcF2SindsByRoi, 
                                   ##SrcImage=SrcIm, 
                                   ##SrcSliceNum=SrcSliceNum,
                                   ##TrgImage=TrgIm,
                                   #SrcImage=ResSrcLabImByRoi[0], 
                                   #SrcSliceNum=ResSrcSliceNum,
                                   #TrgImage=ResSrcLabImByRoi[0],
                                   SrcImage=TrgIm, 
                                   SrcSliceNum=ResSrcSliceNum,
                                   TrgImage=TrgIm,
                                   TrgSliceNum=TrgSliceNum, 
                                   RefImage=TrgIm,
                                   ShiftInX=False,
                                   ShiftInY=False,
                                   ShiftInZ=True,
                                   Fractional=False,
                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'After running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
            print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
            for r in range(len(ResSrcPixArrByRoi)):
                print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}\n')
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
        
        """ Convert ResSrcPixArrByRoi to a list (for each ROI) of a list (for
        each contour) of ContourData, and a list (for each ROI) of a list (for
        each contour) of a list (for each point) of coordinates. """
        
        #print(f'\nBefore running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
        #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
        #for r in range(len(ResSrcPixArrByRoi)):
        #    #print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
        #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
        
        """ Is Thresh = 0.5 below too high? """
            
        ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi, ResSrcC2SindsByRoi\
        = PixArrByRoi2PtsByCntByRoi(PixArrByRoi=ResSrcPixArrByRoi,
                                    F2SindsByRoi=ResSrcF2SindsByRoi,
                                    DicomDir=TrgDcmDir,
                                    Thresh=0.5,
                                    LogToConsole=LogToConsole)
        
        #print(f'\nAfter running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcC2SindsByRoi = {ResSrcC2SindsByRoi}')
        #print(f'len(ResSrcPtsByCntByRoi) = {len(ResSrcPtsByCntByRoi)}')
        #for r in range(len(ResSrcPtsByCntByRoi)):
        #    C = len(ResSrcPtsByCntByRoi[r])
        #    print(f'   There are {C} contours in the {r}^th ROI')
        #    for c in range(C):
        #        P = len(ResSrcPtsByCntByRoi[r][c])
        #        print(f'      There are {P} points in the {c}^th contour')
            
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source pixel arrays to points.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if LogToConsole:
            #print('\n\n', '-'*120)
            print(f'After converting ResPixArrToCopy to ContourData and',
                  f'points, ResSrcC2SindsByRoi =')
            PrintIndsByRoi(ResSrcC2SindsByRoi)
            print('')
        
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgRtsFpath:
                """ An existing Target RTS is to be modified. 
                Concatenate TrgC2SindsByRoi and ResSrcC2SindsByRoi, and
                TrgPtsByCntByRoi and ResSrcPtsByCntByRoi. """
                
                TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + ResSrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
                
                TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + ResSrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
                
            else:
                """ A Target RTS was not provided. """
                
                TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
                
                TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
                
            
            """ Convert points to contour data format. """
            TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        else:
            """ UseCase is '3b', '4b' or '5b'. """
            
            TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
            
            TrgCntDataByCntByRoi = deepcopy(ResSrcCntDataByCntByRoi)
            
            TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'Prior to running CreateRts():')
        print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
        for r in range(len(TrgPtsByCntByRoi)):
            C = len(TrgPtsByCntByRoi[r])
            print(f'   There are {C} in the {r}^th ROI')
            for c in range(len(TrgPtsByCntByRoi[r])):
                P = len(TrgPtsByCntByRoi[r][c])
                print(f'      There are {P} pts in the {c}^th contour')
        print('')
    
    
    """ Create the Target RTS. """
    
    print('*Creating new Target RTS object...\n')
    
    TrgRts = CreateRts(SrcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi,
                       TrgPtsByCntByRoi, TrgC2SindsByRoi, TrgDcmDir,
                       TxtToAddToTrgRoiName, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the RTS object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Create a DICOM Registration Object (if applicable). """
    if UseCaseToApply in ['5a', '5b']:
        # ContentDescription:
        ContentDesc = ''
        
        if Transform in ['rigid', 'affine']:
            """ Create a Spatial Registration DRO. """
            Dro = CreateSpaDro(SrcDcmDir, TrgDcmDir, TxParams, ContentDesc, 
                               LogToConsole)
        else:
            """ Get the bspline grid related registration transform parameters 
            for the final transform parameter map: (10/05/21) """
            GridOrig = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1, 
                                              'GridOrigin')
            
            GridDir = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                             'GridDirection')
            
            GridDims = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                              'GridSize')
            GridDims = [int(item) for item in GridDims]
            
            GridRes = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                             'GridSpacing')
            GridRes = [float(item) for item in GridRes]
            
            VectGridData = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -1,
                                                  'TransformParameters')
            VectGridData = [float(item) for item in VectGridData]
            
            """ Get the non-deformable registration transform parameters 
            for the second-last transform parameter map: (12/05/21) """
            TxParams = GetParamFromSelxImFilt(SelxImFiltOrSitkTx, -2,
                                              'TransformParameters')
            if TxParams != None:
                TxParams = [float(item) for item in TxParams]
            
            """ Create a Deformable Spatial Registration DRO. """
            Dro = CreateDefDro(SrcDcmDir, TrgDcmDir, GridOrig, GridDir, 
                               GridDims, GridRes, VectGridData, TxParams,
                               ContentDesc, LogToConsole)
        
        
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to create the DRO.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
    else:
        Dro = None
    
    
    if LogToConsole:
        print('\nEnd of CopyRts().')
        print('-'*120)
        
    return TrgRts, Dro, DictOfInputs, ListOfInputs, ListOfTimings
    #return TrgRts, TrgPtsByCntByRoi, TrgC2SindsByRoi, ListOfTimings




