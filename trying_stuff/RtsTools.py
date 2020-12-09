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
    Get the list of the Referenced ROI number for each ContourSequence in a
    RTS.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
    ConNumsByRoi - List of list of integers
        The Referenced ROI numbers for all Contour Sequences grouped by ROI
        number.  The integers are from 0, unlike ROI Number which begin from 1. 
        
    Note:
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
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
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
    ------
    
    ListToGroup : list
        List to be grouped.
        
    NumOfConsByRoi : list integers 
        The number of contours in each ContourSequence (counting from 0).
        
        
    Outputs:
    -------
    
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
    ------
    
    Rts : Pydicom object
        ROI object from a RTS file.
    
        
    Outputs:
    -------
    
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
    ------
    
    Rts : Pydicom object
        ROI object from a RTS file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    -------
    
    CIStoSliceInds : list of integers
        List of slice numbers that correspond to each sequence in
        ContourImageSequence.
    """
    
    # Get the Referenced SOP Instance UIDs from the Contour Image Sequence:
    RSOPuids = GetRSOPuidsInCIS(Rts)
    
    CIStoSliceInds = [] 
    
    for RefUid in RSOPuids:
        # Find the matching index of RefUid in SOPuids:
        CIStoSliceInds.append(SOPuids.index(RefUid))
        
    return CIStoSliceInds






def GetRSOPuidsInCSByRoi(Rts):
    """
    Get the list of referenced SOP Instance UIDs for each Contour Sequence
    of a RTS ROI.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
    RSOPuidsByRoi : list of list of strings
        List (for each ROI) of a list (for each contour) of Referenced SOP UIDs.
        
        
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
        
        # Initialise the list of Referenced SOP UIDs for this ROI:
        RSOPuids = []
        
        for sequence in Sequences:
            uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
            
            RSOPuids.append(uid)
        
        RSOPuidsByRoi.append(RSOPuids)
        
    return RSOPuidsByRoi






def GetCStoSliceIndsByRoi(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each Contour Sequence in all ROIs
    in a RTSTRUCT ROI.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    -------
    
    CStoSliceIndsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each sequence in ContourSequence.
    """
    
    # Get the Referenced SOP Instance UIDs from the Contour Sequence by ROI:
    RSOPuidsByRoi = GetRSOPuidsInCSByRoi(Rts)
    
    CStoSliceIndsByRoi = [] 
    
    # Loop through each list of Referenced SOP Instance UIDs:
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
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
    NewRts : Pydicom object
        Modified ROI object with ContourImageSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Rts as a template for NewRts: 
    NewRts = deepcopy(Rts)
        
    # The last item in Contour Image Sequence:
    last = deepcopy(Rts.ReferencedFrameOfReferenceSequence[0]\
                       .RTReferencedStudySequence[0]\
                       .RTReferencedSeriesSequence[0]\
                       .ContourImageSequence[-1])
    
    # Append to the Contour Image Sequence:
    NewRts.ReferencedFrameOfReferenceSequence[0]\
          .RTReferencedStudySequence[0]\
          .RTReferencedSeriesSequence[0]\
          .ContourImageSequence\
          .append(last)
             
    return NewRts









def AddToCS(Rts):
    """
    Append the last item in the Contour Sequence.
    
    Inputs:
    -------
        
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
    NewRts : Pydicom object
        Modified ROI object with ContourSequence lengthened by one.
    """
    
    from copy import deepcopy
    
    # Use Rts as a template for NewRts: 
    NewRts = deepcopy(Rts)
        
    # The last item in Contour Sequence:
    last = deepcopy(Rts.ROIContourSequence[0]\
                       .ContourSequence[-1])
    
    # Append to the Contour Sequence:
    NewRts.ROIContourSequence[0]\
          .ContourSequence\
          .append(last)
             
    return NewRts








def GetPtsByRoi(Rts, DicomDir):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in all 
    contours for all ROIs.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
                       
                            
    Outputs:
    -------
    
    PtsByRoi : list of list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each 
        dimension) of points.
        
    CStoSliceIndsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each list in PtsByRoi. 
    """
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from ConversionTools import ContourData2Points
    
    # Get the DICOM SOP UIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourImageSequence-to-slice indices:
    #CIStoSliceInds = GetCIStoSliceInds(Rts, SopUids)
    
    # Get the ContourSequence-to-slice indices by ROI:
    CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    Nrois = len(CStoSliceIndsByRoi)
    
    if False:
        print(f'\nCStoSliceIndsByRoi = {CStoSliceIndsByRoi}')
        print(f'Nrois = {Nrois}')
    
    PtsByRoi = []
    
    for r in range(Nrois):
        # Number of contours in this ROI:
        NumCS2Sinds = len(CStoSliceIndsByRoi[r])
        
        # Get a list of the Contour Sequences in this ROI:
        ContourSequences = deepcopy(Rts.ROIContourSequence[r].ContourSequence)
        
        NumCS = len(ContourSequences)
        
        if NumCS2Sinds != NumCS:
            msg = f'The number of contour-sequence-to-slice-indices '\
                  + f'({NumCS2Sinds}) must be equal to the number of '\
                  + f'sequences in ROIContourSequence[{r}].ContourSequence '\
                  + f'({NumCS}).'
            
            raise Exception(msg)
        
        # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
        #RefSopUids = [sequence.ContourImageSequence[0]\
        #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        if False:
            print(f'   r = {r}')
            print(f'   NumCS2Sinds = {NumCS2Sinds}')
            print(f'   NumCS = NumCS')
        
        
        PtsByContour = []
        
        
        # Loop through each contour sequence:
        for c in range(NumCS):
            
            # Get the indeces of all matching RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first 
            # match, and there may be more than one!
            #inds = [i for i, e in enumerate(RefSopUids) if e==SopUid]
            
            # Iterate for each index in inds:
            #for ind in inds:
            
            #ContourSequence = deepcopy(ContourSequences[c])
            
            #ContourData = [float(item) for item in ContourSequence.ContourData]
            
            ContourData = Rts.ROIContourSequence[r]\
                             .ContourSequence[c]\
                             .ContourData
            
            Pts = ContourData2Points(ContourData)
            
            PtsByContour.append(Pts)
            
            if False:
                print(f'      c = {c}')
                print(f'      len(Pts) = {len(Pts)}')
            
            
        PtsByRoi.append(PtsByContour)

    return PtsByRoi, CStoSliceIndsByRoi










def GetPtsInRoi_OLD(Rts, DicomDir, SearchString):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in a RTS
    that belong to a ROI with specified label.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    SearchString : string
        All or part of the Source ROI label containing the ROI of interest.
                       
                            
    Outputs:
    -------
    
    PointsByContour : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates. 
        
    CStoSliceInds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        list of points in PointsByContour.
    """
    
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    # The ROI Contour Sequence:
    RoiContourSequence = Rts.ROIContourSequence
    
    # The number of ROIs:
    NumOfRois = len(RoiContourSequence)
    
    RoiNum = GetRoiNum(Rts, SearchString)
    
    if RoiNum > NumOfRois - 1:
        raise Exception(f"There are only {NumOfRois} ROIs in the RTS.\n",
                        f"RoiNum = {RoiNum} > {NumOfRois}.")
    
    
    # The Contour Sequences in the ROI of interest:
    ContourSequence = Rts.ROIContourSequence[RoiNum].ContourSequence
        
    CS = len(ContourSequence)
    
    PointsByContour = []
    
    for i in range(CS):
        ContourData = ContourSequence[i].ContourData
        
        ContourData = [float(item) for item in ContourData]
        
        Points = []
                
        # Iterate for all points in threes:
        for p in range(0, len(ContourData), 3):
            point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
            
            Points.append(point)
            
        PointsByContour.append(Points)
    
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir)
    
    # The ContourSequence-to-slice indices by ROI:
    CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(Rts, SOPuids)
    
    CStoSliceInds = CStoSliceIndsByRoi[RoiNum]
        
    return PointsByContour, CStoSliceInds
















def GetPtsInRoi(Rts, DicomDir, SearchString):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in a RTS
    that belong to a ROI with specified label.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    SearchString : string
        All or part of the ROI name of the ROI of interest.
                       
                            
    Outputs:
    -------
    
    PointsByContour : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates. 
        
    CStoSliceInds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        list of points in PointsByContour.
    """
    
    from DicomTools import GetRoiNum
    
    PtsByRoi, CStoSliceIndsByRoi = GetPtsByRoi(Rts, DicomDir)
    
    RoiNum = GetRoiNum(Rts, SearchString)
    
    PointsByContour = PtsByRoi[RoiNum]
    
    CStoSliceInds = CStoSliceIndsByRoi[RoiNum]
    
    #print('\n\nResults of GetPtsInRoi:')
    #print(f'   len(PtsByRoi) = {len(PtsByRoi)}')
    #print(f'   len(CStoSliceIndsByRoi) = {len(CStoSliceIndsByRoi)}')
    #print(f'   len(PointsByContour) = {len(PointsByContour)}')
    #print(f'   len(CStoSliceInds) = {len(CStoSliceInds)}')
    #print('end of GetPtsByRoi.\n')
        
    return PointsByContour, CStoSliceInds












def GetPtsInContour(Rts, SearchString, SliceNum, DicomDir):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in a 
    contour in a ROI.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        RTS object.
        
    SearchString : string
        All or part of the ROI name of the ROI containing the contour to be
        copied.
        
    SliceNum : integer
        The slice number (counting from 0) that corresponds to the contour of
        interest.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
                            
    Outputs:
    -------
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of coordinates. 
        
    CStoSliceInd : list of a single integer
        List (of length one - the one contour defined by Points) of the slice 
        number that corresponds to Points. 
    """
    
    from copy import deepcopy
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    from ConversionTools import ContourData2Points
    
    # The ROI Contour Sequence:
    RoiContourSequence = deepcopy(Rts.ROIContourSequence)
    
    # The number of ROIs:
    NumOfRois = len(RoiContourSequence)
    
    RoiNum = GetRoiNum(Rts, SearchString)
    
    if RoiNum > NumOfRois - 1:
        raise Exception(f"There are only {NumOfRois} ROIs in the RTS.\n",
                        f"RoiNum = {RoiNum} > {NumOfRois}.")
        
    # Get the DICOM SOP UIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourSequence-to-slice indices by ROI:
    CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    CStoSliceInds = CStoSliceIndsByRoi[RoiNum]
    
    if not SliceNum in CStoSliceInds:
        raise Exception(f"There is no contour on slice {SliceNum} in the RTS.")
        
    # The contour number of interest:
    ContourNum = CStoSliceInds.index(SliceNum)
    
    # The Contour Sequences in the ROI of interest:
    #ContourSequence = deepcopy(Rts.ROIContourSequence[RoiNum].ContourSequence)
    ContourSequence = deepcopy(RoiContourSequence[RoiNum].ContourSequence)
    
    # The number of contours in this ROI:
    NumOfContours = len(ContourSequence)
    
    if ContourNum > NumOfContours - 1:
        raise Exception(f"There are only {NumOfContours} contours in the RTS.",
                        f"\nContourNum = {ContourNum} > {NumOfContours}.")
        
    ContourData = deepcopy(ContourSequence[ContourNum].ContourData)
    
    Points = ContourData2Points(ContourData)
    
    #CStoSliceInd = [SliceNum]
    
    if False:
        print(f'\n\nResults of GetPtsInContour:')
        print(f'   NumOfRois = {NumOfRois}')
        print(f'   RoiNum = {RoiNum}')
        print(f'   CStoSliceIndsByRoi = {CStoSliceIndsByRoi}')
        print(f'   CStoSliceInds = CStoSliceIndsByRoi[RoiNum={RoiNum}] =',
              f'{CStoSliceInds}')
        print(f'   ContourNum = {ContourNum}')
        print(f'   len(Points) = {len(Points)}')
        print(f'   CStoSliceInds = {CStoSliceInds}')
        
        numofpts = []
        for c in range(NumOfContours):
            contourdata = deepcopy(ContourSequence[c].ContourData)
        
            pts = ContourData2Points(contourdata)
            
            numofpts.append(len(pts))
            
            
        print(f'   Number of pts in each CS = {numofpts}')
        print('end of results.')
    
    return Points#, CStoSliceInds









def AddCopiedPts(OrigPtsByCnt, PtsToAdd, CStoSliceInds, SliceNum):
    """
    Add points in a contour to an existing list of points-by-contour.
       
    
    Inputs:
    ------
    
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
    -------
    
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
    ------
    
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
    -------
    
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
        print(f'   NewCStoSliceInds after removing duplicates = {NewCStoSliceInds}')
    
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
    Get the indices of the points (i.e. in the Image Coordinate System) in all 
    contours for all ROIs.
    
    
    Inputs:
    ------
    
    Rts : Pydicom object
        RTS object.
        
    DicomDir : sString
        Directory containing the corresponding DICOMs.
                       
                            
    Outputs:
    -------
    
    IndsInRois - List of list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each 
        dimension) of indices.
    """
    
    from ImageTools import GetImageAttributes
    from ConversionTools import ContourToIndices
    
    PtsByRoi = GetPtsByRoi(Rts, DicomDir)
    
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    
    Nrois = len(PtsByRoi)
    
    IndsByRoi = []
    
    for r in range(Nrois):
        IndsByContour = []
        
        # The list of points in the contours for this ROI:
        PtsByContour = PtsByRoi[r]
        
        for contour in PtsByContour:
            Indices = ContourToIndices(contour, IPPs[0], Dirs, Spacings)
                
            IndsByContour.append(Indices)
            
        IndsByRoi.append(IndsByContour)
        
    
    return IndsByRoi









def GetMaskFromRts(Rts, SearchString, SliceNum, DicomDir, RefImage):
    """
    Get a 2D mask (pixel array) from a RTS.
       
    
    Inputs:
    ------
    
    Rts : Pydicom object
        DICOM-RTSTRUCT object.
        
    SearchString : string
        All or part of the ROI name of the ROI containing the contour to be
        copied.
        
    SliceNum : integer
        Slice index of the DICOM stack corresponding to the contour to be
        copied (counting from 0).
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
        
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
    
        
    Outputs:
    -------
    
    Mask : Numpy data array
        A 2D mask of the contour specified by FromRoiLabel and FromSliceNum.
    """
    
    #import importlib
    #import ConversionTools
    #importlib.reload(ConversionTools)
    
    from ImageTools import GetImageAttributes
    from ConversionTools import Point2Index
    from ConversionTools import Indices2Mask
    
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    # Get the contour points (in the Patient Coordinate System):
    #Points, CStoSliceInds = GetPtsInContour(Rts, SearchString, SliceNum, 
    #                                        DicomDir)
    Points = GetPtsInContour(Rts, SearchString, SliceNum, DicomDir)
    
    #print(f'\nlen(Points) = {len(Points)}')
    #print(f'\nPoints = {Points}')
    
    # Convert the points to indices (i.e. the Image Coordinate System):
    Inds = []
            
    for point in Points:
        #print(f'\nlen(point) = {len(point)}')
        
        inds = Point2Index(Point=point, Origin=IPPs[0], 
                           Directions=Dirs, Spacings=Spacings)
        
        # Convert fractional indices to integers:
        inds = [int(round(item)) for item in inds]
        
        Inds.append(inds)
    
    
    Mask = Indices2Mask(Inds, RefImage)
    
    return Mask









#def InitialiseRts(DicomDir, RtsTemplate, LogToConsole=False):
def InitialiseRts(RtsTemplate, RoiNum, CStoSliceInds, DicomDir, 
                  LogToConsole=False):
    """
    Initialise a RTS object based on an existing RTS object (RtsTemplate).
         
    
    Inputs:
    ------
                              
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
                           
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    Rts : Pydicom object
        New RTS object.
    """
    
    from pydicom.uid import generate_uid
    import time
    from copy import deepcopy
    from DicomTools import ImportDicoms
    
    # The ROI Name of the ROI containing the contour to be copied:
    RoiName = RtsTemplate.StructureSetROISequence[RoiNum].ROIName
    
    
    # Use RtsTemplate as a template for NewRts:
    Rts = deepcopy(RtsTemplate)
    
    Dicoms = ImportDicoms(DicomDir)
    
    # Modify various tags to match the values in Dicoms:
    Rts.StudyDate = Dicoms[0].StudyDate
    Rts.StudyTime = Dicoms[0].StudyTime
    Rts.SeriesDate = Dicoms[0].SeriesDate
    Rts.PatientName = Dicoms[0].PatientName
    Rts.PatientID = Dicoms[0].PatientID
    Rts.PatientBirthDate = Dicoms[0].PatientBirthDate
    Rts.PatientSex = Dicoms[0].PatientSex
    Rts.StudyInstanceUID = Dicoms[0].StudyInstanceUID
    Rts.StudyID = Dicoms[0].StudyID
    Rts.FrameOfReferenceUID = Dicoms[0].FrameOfReferenceUID
    
    
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
    
    """ 
    Modify the number of sequences in the (only) ContourSequence.  
    The number of sequences must be equal to the length of CStoSliceInds.
    """
    
    # The number of sequences in the last ContourSequence:
    OrigCS = len(Rts.ROIContourSequence[-1].ContourSequence)
    
    ReqCS = len(CStoSliceInds)
    
    if ReqCS > OrigCS:
        for i in range(ReqCS - OrigCS):
            Rts.ROIContourSequence[-1].ContourSequence\
               .append(deepcopy(Rts.ROIContourSequence[-1].ContourSequence[-1]))          
    else:
        for i in range(OrigCS - ReqCS):
            Rts.ROIContourSequence[-1].ContourSequence.pop()

    
    
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
    
    
    
    """ 
    Modify other tags.
    """
    # Generate a new SOP Instance UID:
    Rts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    Rts.SeriesInstanceUID = generate_uid()
    
    Rts.StructureSetLabel = 'Contour(s)_copied_from_' + Rts.StructureSetLabel
    
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








def ModifyRts(Rts, CntDataByCnt, PtsByCnt, CStoSliceInds, DicomDir, 
              LogToConsole=False):
    """
    Modify various tag values of the new RTS ROI so that they are consistent
    with the corresponding tag values of the DICOMs (e.g. SOPClassUID, 
    SOPInstanceUID), and so that the ROIContourSequence values are consistent
    with those of the original RTS ROI for existing contours, and with the
    contour that has been copied from the original RTS ROI.
         
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
    ------
                              
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
    ------
        
    Rts : Pydicom object
        Modified RTS ROI object.
    """
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
    
    for i in range(len(CStoSliceInds)):
        # The DICOM slice index for this sequence:
        s = CStoSliceInds[i]
        
        if LogToConsole:
            print(f'\n   i = {i}, CStoSliceInds[{i}] = {CStoSliceInds[i]}')
            #print(f'      DICOM slice number = {s}')
        
        # Modify the Referenced SOP Instance UID in ContourImageSequence:
        Rts.ReferencedFrameOfReferenceSequence[0]\
           .RTReferencedStudySequence[0]\
           .RTReferencedSeriesSequence[0]\
           .ContourImageSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # Modify the Referenced SOP Instance UID in ROIContourSequence:
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
           .ContourImageSequence[0]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # The contour points for this contour:
        Points = PtsByCnt[i]
        
        # Modify NumberOfContourPoints:
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
           .ContourImageSequence[0]\
           .NumberOfContourPoints = f"{len(Points)}"
        
        # The contour points as a flat list of coordinates as strings:
        #ContourData = Points2ContourData(Points)
        #ContourData = CntDataByCnt[i]
        
        #print(f'\nCntDataByCnt[{i}] = {CntDataByCnt[i]}')
        
        # Modify ContourData:
        #Rts.ROIContourSequence[0]\
        #   .ContourSequence[i]\
        #   .ContourImageSequence[0]\  <--- found on 8/12/20
        #   .ContourData = CntDataByCnt[i]
           
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
           .ContourData = CntDataByCnt[i]
        
        if LogToConsole:
            print(f'      len(Points) = {len(Points)}')
    
    
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










def ExportRts(TrgRts, SrcRtsFpath, NamePrefix, ExportDir):
    """
    Export contour ROI Object to disk.  The following tags will also be 
    modified: 
        StructureSetDate
        StructureSetTime
        StructureSetLabel
        StructureSetROISequence[0].ROIName
    
    
    Inputs:
        TrgRts   - Target RT-STRUCT ROI object to be exported
        
        SrcRtsFpath - (String) Full path of the Source DICOM-RTSTRUCT file
                      (used to generate the filename of the new RTS file)
        
        NamePrefix  - (String) Prefix to be added to the assigned filename and
                      ROI Name (e.g. 'MR4_S9_s23_to_MR4_S9_s22') 
                            
        ExportDir   - (String) Directory where the RTSTRUCT is to be exported
                           
    
                            
    Returns:
        TrgRtsFpath - (String) Full path of the exported Target DICOM-RTSTRUCT 
                      file
    
    """
    
    """
    The following was moved into the main function CopyContourWithinSeries:
    # Generate a new SOP Instance UID:
    Rts.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    Rts.SeriesInstanceUID = pydicom.uid.generate_uid()
    """
    
    import os
    import time
    
    # Get the filename of the original RT-Struct file:
    SrcRtsFname = os.path.split(SrcRtsFpath)[1]
    
    # Modify the Structure Set Date and Time to the present:
    NewDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    FnamePrefix = NewDateTime + '_' + NamePrefix + '_from_' 
    #RoiNamePrefix = NamePrefix + '_from_'
    
    # Modify the Structure Set Label (this appears under Name in XNAT):
    ##TrgRts.StructureSetLabel = 'Copy_of_' + TrgRts.StructureSetLabel
    #TrgRts.StructureSetLabel = RoiNamePrefix + TrgRts.StructureSetLabel
    
    # Modify the ROI Name for the Structure Set ROI Sequence:             
    #TrgRts.StructureSetROISequence[0].ROIName = RoiNamePrefix \
    #                                            + TrgRts.StructureSetLabel
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = 'New_' + SrcRtsFname + '_' + NewDate + '_' + NewTime
    TrgRtsFname = FnamePrefix + SrcRtsFname
    
    TrgRtsFpath = os.path.join(ExportDir, TrgRtsFname)
    
    TrgRts.save_as(TrgRtsFpath)
        
    print('\nRTS ROI exported to:\n\n', TrgRtsFpath)
    
    return TrgRtsFpath







def CentroidOfContour(Contour):
    """
    Compute the centroid of a list of points (or indices).  
    
    Inputs:
        Contour  - (List of list of floats/integers) List of a list of points 
                   or indices
        
    Returns:
        Centroid - (List of floats/integers) Centroid of Contour

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
        Contours - (List of list of a list of floats/integers) 
                   List of a list of points/indices for each contour
        
    Returns:
        RoiNums  - (List of integers) Integer label assigning each contour to
                   a ROI

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