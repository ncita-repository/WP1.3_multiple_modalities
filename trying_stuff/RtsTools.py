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
#from ConversionTools import Inds2Mask


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






def GetRSOPuidsInCS(Rts):
    """
    Get the list of referenced SOP Instance UIDs for each Contour Sequence
    of a RTS ROI.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    
    Outputs:
    -------
    
    RSOPuids : list of list of strings
        List of Referenced SOP UIDs.
        
        
    Notes:
    -----
    
    If there is only one ROI Contour Sequence RSOPuids will be a list of
    length 1 of a list of strings, e.g.
        
    [['uid1', 'uid2', ...]]
    """
    
    RSOPuids = []
    
    rois = Rts.ROIContourSequence
    
    for roi in rois:
        sequences = roi.ContourSequence
        
        # Initialise the list of Referenced SOP UIDs for this ROI:
        rsopuids = []
        
        for sequence in sequences:
            uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
            
            rsopuids.append(uid)
        
        RSOPuids.append(rsopuids)
        
    return RSOPuids






def GetCStoSliceInds(Rts, SOPuids):
    """
    Get the slice numbers that correspond to each Contour Sequence in a 
    RTSTRUCT ROI.
    
    Inputs:
    ------
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    SOPuids : list of strings
        List of SOP UIDs of the DICOMs.
        
    
    Outputs:
    -------
    
    CStoSliceInds : list of integers
        List of slice numbers that correspond to each sequence in 
        ContourSequence.
    """
    
    # Get the Referenced SOP Instance UIDs from the Contour Sequence:
    RSOPuids = GetRSOPuidsInCS(Rts)
    
    CStoSliceInds = [] 
    
    # Loop through each list of Referenced SOP Instance UIDs (there may be only
    # a single list of lists of UIDs):
    for i in range(len(RSOPuids)):
        inds = []
    
        for RefUid in RSOPuids[i]:
            # Find the matching index of RefUid in SOPuids:
            inds.append(SOPuids.index(RefUid))
            
        CStoSliceInds.append(inds)
        
    return CStoSliceInds








def AddToRtsSequences_OLD(Rts):
    """
    Append the last item in the following sequences in the RTS Object:
        Contour Image Sequence
        Contour Sequence
    
    Inputs:
        Rts    - ROI Object from an RTSTRUCT file
        
    Returns:
        NewRts - Modified ROI Object with sequences lengthened by one
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

    # The last item in Contour Sequence:
    last = deepcopy(Rts.ROIContourSequence[0]\
                       .ContourSequence[-1])
    
    # Append to the Contour Sequence:
    NewRts.ROIContourSequence[0]\
          .ContourSequence\
          .append(last)
             
    return NewRts








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








def GetPtsInContour(RtsFpath, FromRoiLabel, FromSliceNum, DicomDir):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in the 
    contour on slice FromSliceNum in ROI FromRoiLabel.
    
    
    Inputs:
    ------
    
    RtsFpath : string
        Full path of the DICOM-RTSTRUCT file.
        
    FromRoiLabel : string
        All or part of the Source ROI label containing the contour to be
        copied.
        
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour to
        be copied (counting from 0).
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
                       
                            
    Outputs:
    -------
    
    Points : list of list of floats
        List (for each polygon) of a list (for each dimension) of points in a
        contour.
    """
    
    from pydicom import dcmread
    from DicomTools import GetRoiLabels
    from DicomTools import GetDicomSOPuids
    
    # Import the RTSTRUCT ROI:
    Rts = dcmread(RtsFpath)
    
    # Get the ROI labels:
    RoiLabels = GetRoiLabels(Roi=Rts)
    
    if FromRoiLabel in RoiLabels:
        RoiNum = RoiLabels.index(FromRoiLabel)
    else:
        raise Exception("There is no ROI with label" + FromRoiLabel + ".")
        
    
    # The ROI Contour Sequence containing the contour to be copied:
    #Sequence = Rts.ROIContourSequence[RoiNum]
    
    # The Referenced SOP Instance UIDs in the Contour Sequence grouped by ROI:
    RefSOPsByRoi = GetRSOPuidsInCS(Rts)
    
    # The Referenced SOP Instance UIDs in the Contour Sequence containing the 
    # contour to be copied:
    RefSOPs = RefSOPsByRoi[RoiNum]
    
    # The DICOM SOP UIDs:
    SOPs = GetDicomSOPuids(DicomDir)
    
    # The SOP corresponding to the contour to be copied:
    SOP = SOPs[FromSliceNum]
    
    if SOP in RefSOPs:
        # The Contour Sequence number for the contour to be copied:
        SeqNum = RefSOPs.index(SOP)
    else:
        raise Exception("There is no contour data for slice " + FromSliceNum + ".")
    
    ContourData = Rts.ROIContourSequence[RoiNum]\
                     .ContourSequence[SeqNum]\
                     .ContourData
    
    ContourData = [float(item) for item in ContourData]
    
    
    Points = []
    
    # Iterate for all points in threes:
    for p in range(0, len(ContourData), 3):
        point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        Points.append(point)
    

    return Points









def GetMaskFromRts(RtsFpath, FromRoiLabel, FromSliceNum, DicomDir, RefImage):
    """
    Get a 2D mask (pixel array) from a RTS.
       
    
    Inputs:
    ------
    
    RtsFpath : string
        Full path of the DICOM-RTSTRUCT file.
        
    FromRoiLabel : string
        All or part of the Source ROI label containing the contour to be
        copied.
        
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour to
        be copied (counting from 0).
        
    DicomDir : string
        Directory containing the corresponding DICOMs.
        
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
                       
                            
    Outputs:
    -------
    
    Points : list of list of floats
        List (for each polygon) of a list (for each dimension) of points in a
        contour.
        
    
        
    Outputs:
    -------
    
    Mask : Numpy data array
        A 2D mask of the contour specified by FromRoiLabel and FromSliceNum.
    """
    
    from ImageTools import GetImageAttributes
    from ConversionTools import Point2Index
    from ConversionTools import Inds2Mask
    
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    # Get the contour points (in the Patient Coordinate System):
    Points = GetPtsInContour(RtsFpath, FromRoiLabel, FromSliceNum, DicomDir)
    
    # Convert the points to indices (i.e. the Image Coordinate System):
    Inds = []
            
    for point in Points:
        inds = Point2Index(Point=point, Origin=IPPs[0], 
                           Directions=Dirs, Spacings=Spacings)
        
        # Convert fractional indices to integers:
        inds = [int(round(item)) for item in inds]
        
        Inds.append(inds)
    
    
    Mask = Inds2Mask(Inds, RefImage)
    
    return Mask








def ChangeNumOfRtsSequences(Rts, CIStoSliceInds, CStoSliceIndsByRoi, 
                            LogToConsole=False):
    """
    Adjust the number of sequences in a RTS object based so that they are 
    consistent with the number of DICOMs in the series and the number of
    contours in the contour data.
         
    
    Inputs:
    ------
                              
    Rts : Pydicom object
        RTS object.
    
    NewTrgCIStoSliceInds : List of integers
        List (for each sequence) of slice numbers that correspond to each 
        ContourImageSequence.
    
    NewTrgCStoSliceIndsByRoi : List of a list of integers
        List (for each ROI) of a list (for each sequence) of slice numbers that 
        correspond to each ContourSequence.
        
    NumOfContours : integer
        The number of contours in the RTS object's contour data.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    -------
        
    Rts : Pydicom object
        RTS object with modified number of sequences.
    """
    

    
    """ 
    Modify the number of sequences in ContourImageSequence. 
    The number of sequences must be equal to the length of CIStoSliceInds.
    """
    
    from copy import deepcopy
    
    # The original and required number of sequences in ContourImageSequence:
    OrigCIS = len(Rts.ReferencedFrameOfReferenceSequence[0]\
                     .RTReferencedStudySequence[0]\
                     .RTReferencedSeriesSequence[0]\
                     .ContourImageSequence)
    
    ReqCIS = len(CIStoSliceInds)
    
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
    The number of sequences must be equal to the length of CStoSliceIndsByRoi.
    """
    
    OrigS = len(Rts.StructureSetROISequence)
    
    ReqS = len(CStoSliceIndsByRoi) 
    
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
    The number of sequences must be equal to the length of CStoSliceIndsByRoi.
    """
    
    OrigR = len(Rts.ROIContourSequence)
    
    ReqR = len(CStoSliceIndsByRoi) 
    
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
    Modify the number of sequences in the final ContourSequence (the other
    sequences should remain unchanged as the new contour(s) are copied to a new
    contour sequence - i.e. a new ROI).  The number of sequences must be equal 
    to the length of CStoSliceIndsByRoi[-1].
    """
    
    # The number of sequences in the last ContourSequence:
    OrigCS = len(Rts.ROIContourSequence[-1].ContourSequence)
    
    ReqCS = len(CStoSliceIndsByRoi[-1])
    
    if ReqCS > OrigCS:
        for i in range(ReqCS - OrigCS):
            Rts.ROIContourSequence[-1].ContourSequence\
               .append(deepcopy(Rts.ROIContourSequence[-1].ContourSequence[-1]))          
    else:
        for i in range(OrigCS - ReqCS):
            Rts.ROIContourSequence[-1].ContourSequence.pop()

    
    
    """ 
    Modify the number of sequences in RTROIObservationsSequence. 
    The number of sequences must be equal to the length of CStoSliceIndsByRoi.
    """
    
    OrigO = len(Rts.RTROIObservationsSequence)
    
    ReqO = len(CStoSliceIndsByRoi)
    
    if ReqO > OrigO:
        for i in range(ReqO - OrigO):
            Rts.RTROIObservationsSequence\
               .append(deepcopy(Rts.RTROIObservationsSequence[-1]))          
    else:
        for i in range(OrigO - ReqO):
            Rts.RTROIObservationsSequence.pop()
    
    
    
    
    if True:#LogToConsole:
        NewCIS = len(Rts.ReferencedFrameOfReferenceSequence[0]\
                        .RTReferencedStudySequence[0]\
                        .RTReferencedSeriesSequence[0]\
                        .ContourImageSequence)
        
        NewS = len(Rts.StructureSetROISequence)
        
        NewR = len(Rts.ROIContourSequence)
        
        NewCS = len(Rts.ROIContourSequence[-1].ContourSequence)
        
        NewO = len(Rts.RTROIObservationsSequence)
        
        print(f'\nThere are {len(CIStoSliceInds)} images')
        
        print(f'\nThere are {len(CStoSliceIndsByRoi[-1])} contours in the',
              'final contour sequence')
        
        print(f'\nThere were {OrigCIS} sequences in ContourImageSequence. ',
              f'Now there are {NewCIS} sequences.')
        
        print(f'\nThere were {OrigS} sequences in StructureSetROISequence.',
              f'Now there are {NewS} sequences.')
        
        print(f'\nThere were {OrigR} sequences in ROIContourSequence.',
              f'Now there are {NewR} sequences.')
        
        print(f'\nThere were {OrigCS} sequences in the final ContourSequence.',
              f' Now there are {NewCS} sequences.')
        
        print(f'\nThere were {OrigO} sequences in RTROIObservationsSequence.',
              f'Now there are {NewO} sequences.')
        
    return Rts








def InitialiseRts(Dicoms, Image, RtsTemplate, LogToConsole=False):
    """
    Initialise a RTS object based on an existing RTS object (RtsTemplate).
         
    
    Inputs:
    ------
                              
    Dicoms : List of Pydicom objects
        List of DICOM Objects that relate to the RTS object to be created.
    
    CStoSliceInds : List of integers
        Slice numbers that correspond to each frame in PixArr.
        
    Image : SimpleITK image
        3D SimpleITK image of the DICOM stack that relates to the new RTS 
        object (required to obtain the SliceThickness and Spacing Between
        Slices in Pixel Measures Sequence)
                           
    RtsTemplate : Pydicom object
        RTS object that will be used as a template to create the new object.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    Rts : Pydicom object
        New RTS object.
    """
    
    from copy import deepcopy
    
    # Use RtsTemplate as a template for NewRts:
    Rts = deepcopy(RtsTemplate)
    
    """
    26/11: Continue modifying this
    """
    
    # Modify various tags to match the values in Dicoms:
    """
    Note:
    
    The following changes are made in ModifyRtsTagVals():
        - a new SOP Instance UID is generated
        - a new Series Instance UID is generated
        - Series Description will be modified
        - Referenced SOP Instance UIDs in Referenced Series Sequence will be
        modified
        - Segment Label for the new (last) segment will be modified
        - Image Orientation Patient (in Shared Functional Groups Sequence) will
        be modified
    """
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
    
    """ This was moved to ModifyRtsTagVals: """
    # Modify the number of seguences in NewRts:
    #Rts = ChangeNumOfRtsSequences(Rts=Rts, NumOfContours=len(CStoSliceInds), 
    #                              LogToConsole=LogToConsole)

    return Rts







def ModifyRtsTagVals(TrgRts, TrgDicoms, NewTrgPtsByRoi, NewTrgRoiNums,
                     NewTrgCIStoSliceInds, NewTrgCStoSliceIndsByRoi, 
                     FromRoiNum, FromSliceNum, SrcRts,
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
                              
    TrgRts : Pydicom object
        Original Target RTS ROI Object.
    
    TrgDicoms : List of Pydicom objects
        List of DICOM Objects that include the DICOMs that TrgRts relate to.
    
    NewTrgPtsByRoi : list of a list of a list of floats
        List (for each ROI) of a list (for each contour sequence) of a list 
        (for each dimension) of physical coordinates for Target.
    
    NewTrgRoiNums : List of integers
        ROI numbers that correspond to each ROI in NewTrgRts.
    
    NewTrgCIStoSliceInds : List of integers
        List (for each sequence) of slice numbers that correspond to each 
        ContourImageSequence.
    
    NewTrgCStoSliceIndsByRoi : List of a list of integers
        List (for each ROI) of a list (for each sequence) of slice numbers that 
        correspond to each ContourSequence.
                           
    FromRoiNum : integer
        Source ROI number of the contour that was copied (counting from 0) 
        (this is only used for generating a new Series Description).
    
    FromSliceNum : integer
        Source slice index in the DICOM stack corresponding to the contour 
        that was copied (counting from 0) (this is only used for generating a 
        new Series Description).
                           
    SrcRts : Pydicom object
        Source RTS ROI object (this is only used for generating a new Series 
        Description).
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    NewTrgRts : Pydicom object
        New Target RTS ROI object.
    """
    
    """
    26/11: Continue modifying this
    """
    
    from copy import deepcopy
    from pydicom.uid import generate_uid
    import time
    from ConversionTools import Pts2ContourData
    
    # Use TrgRts as a template for NewTrgRts:
    NewTrgRts = deepcopy(TrgRts)
    
    # Modify the number of seguences in NewRts:
    NewTrgRts = ChangeNumOfRtsSequences(Rts=NewTrgRts, 
                                        CIStoSliceInds=NewTrgCIStoSliceInds,
                                        CStoSliceIndsByRoi=NewTrgCStoSliceIndsByRoi, 
                                        LogToConsole=LogToConsole)
    
    
    # Generate a new SOP Instance UID:
    NewTrgRts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRts.SeriesInstanceUID = generate_uid()
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    NewTrgRts.StructureSetDate = NewDate
    NewTrgRts.StructureSetTime = NewTime
    
    
    # Modify the ContourImageSequence:
    for i in range(len(NewTrgCIStoSliceInds)):
        # The DICOM slice index for this sequence:
        s = NewTrgCIStoSliceInds[i]
        
        if LogToConsole:
            print(f'\nNewTrgCIStoSliceInds[{i}] = {NewTrgCIStoSliceInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Instance UID:
        NewTrgRts.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
    
    
    # Modify the StructureSetROISequence:
    for i in range(len(NewTrgRts.StructureSetROISequence)):
        NewTrgRts.StructureSetROISequence[i].RoiNumber = f"i + 1"
        
        NewTrgRts.StructureSetROISequence[i]\
                 .ReferencedFrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
        
    
    FromRoiName = deepcopy(SrcRts.StructureSetROISequence[FromRoiNum].RoiName)
    
    #ToRoiName = f'Copy_of_{FromSliceNum}_' + FromRoiName
    ToRoiName = 'New_' + FromRoiName
             
    ToRoiNum = len(NewTrgCStoSliceIndsByRoi)
    
    NewTrgRts.StructureSetROISequence[ToRoiNum].RoiName = ToRoiName
    
    
    # Modify ROIContourSequence:
    for r in range(len(NewTrgCStoSliceIndsByRoi)):
        # r = ROI number
        for c in range(len(NewTrgCStoSliceIndsByRoi[r])):
            # c = contour sequence
            
            # The DICOM slice index for this sequence:
            s = NewTrgCStoSliceIndsByRoi[r][c]
            
            # Modify the Referenced SOP Instance UID:
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[c]\
                     .ContourImageSequence[0]\
                     .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
            
            # The contour points as a list (for each point) of a list (for each
            # dimension) of coordinates:
            Points = NewTrgContourData[r][c]
            
            # Modify NumberOfContourPoints:
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[c]\
                     .ContourImageSequence[0]\
                     .NumberOfContourPoints = f"{int(len(NewTrgPtsByRoi))}"
            
            # The contour points as a flat list of coordinates as strings:
            ContourData = Pts2ContourData(Points)
            
            # Modify ContourData:
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[c]\
                     .ContourImageSequence[0]\
                     .ContourData = ContourData
    
    
    
    # Modify ROIContourSequence:
    for r in range(len(NewTrgCStoSliceIndsByRoi)):
        NewTrgRts.ObservationNumber = f"{r + 1}"
        
        NewTrgRts.ReferencedROINumber = f"{r + 1}"
    
                                          
    return NewTrgRts









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