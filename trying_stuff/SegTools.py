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



def GetRSOPuidsInRIS(Seg):
    """
    Get the list of referenced SOP Instance UIDs in the Referenced Instance
    Sequence of a SEG ROI.
    
    Inputs:
    ------
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    -------
        
    RSOPuids : list of strings
        Referenced SOP UIDs.
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
    Get the slice numbers that correspond to each Referenced Instance Sequence
    in a SEG ROI.
    
    Inputs:
    ------
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOP UIDs of the DICOMs.
        
    
    Outputs:
    -------
    
    RIStoSliceInds : list of integers
        Slice numbers that correspond to each Referenced Instance Sequence in 
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
    Get the list of referenced SOP Instance UIDs in the Per-Frame Functional
    Groups Sequence of a SEG ROI.
    
    Inputs:
    ------
        
    Seg : Pydicom object
        ROI object from a SEG file.
    
        
    Outputs:
    -------
    
    RSOPuids : list of strings
        Referenced SOP UIDs.
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
    Get the slice numbers that correspond to each Per-frame Functional Groups
    Sequence in a SEG ROI.
    
    Inputs:
    ------
        
    Seg : Pydcom object
        ROI object from a SEG file.
        
    SOPuids : list of strings
        SOP UIDs of the DICOMs.
        
    
    Outputs:
    -------
        
    PFFGStoSliceInds : list of integers
        Slice numbers that correspond to each Per-frame Functional Groups 
        Sequence in Seg.
    """
    
    # Get the list of Referenced SOP Instance UIDs:
    RSOPuids = GetRSOPuidsInPFFGS(Seg)
    
    PFFGStoSliceInds = [] 
    
    for Ruid in RSOPuids:
        # Find the matching index of Ruid in SOPuids:
        PFFGStoSliceInds.append(SOPuids.index(Ruid))
        
    return PFFGStoSliceInds





def GetDIVs(Seg):
    """
    Get the DimensionIndexValues in a SEG ROI.
    
    Inputs:
    ------
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
        
    Outputs:
    -------
    
    DIVs : list of lists of integers
        Dimension Index Values in Seg.
        
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
    Group a list of items by common segment/ROI using the Dimension Index 
    Values.
    
    Inputs:
    ------
    
    ListToGroup : list
        The list to be grouped.
        
    DIVs : list of lists of integers
        List of Dimension Index Values.
        
        
    Outputs:
    -------
    
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
    Append the last item in the Per-Frame Functional Groups Sequence of the
    SEG Object to itself (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ------
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    -------
    
    NewSeg : Pydicom object
        Modified ROI object with Per-Frame Functional Groups Sequence 
        lengthened by one.
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
    Append the last item in the Segment Sequence of the SEG Object to itself 
    (i.e. increase the length of the sequence by one). 
    
    Inputs:
    ------
    
    Seg : Pydicom object
        ROI object from a SEG file.
        
    
    Outputs:
    -------
    
    NewSeg : Pydicom object
        Modified ROI Object with Segment Sequence lengthened by one.
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
    ------
                              
    Seg : Pydicom object
        SEG object.
        
    Dicoms : List of Pydicom objects
        List of DICOM Objects that relate to the SEG object to be created.
    
    NumOfFrames : integer
        The number of frames in the SEG object's pixel array.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    -------
        
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








def InitialiseSeg(Dicoms, PFFGStoSliceInds, Image, SegTemplate,
                  LogToConsole=False):
    """
    Initialise a SEG object based on an existing SEG object (SegTemplate).
         
    
    Inputs:
    ------
                              
    Dicoms : List of Pydicom objects
        List of DICOM Objects that relate to the SEG object to be created.
    
    PFFGStoSliceInds : List of integers
        Slice numbers that correspond to each frame in PixArr.
        
    Image : SimpleITK image
        3D SimpleITK image of the DICOM stack that relates to the new SEG 
        object (required to obtain the SliceThickness and Spacing Between
        Slices in Pixel Measures Sequence)
                           
    SegTemplate : Pydicom object
        SEG object that will be used as a template to create the new object.
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    Seg : Pydicom object
        New SEG object.
    """
    
    from copy import deepcopy
    
    # Use SegTemplate as a template for NewSeg:
    Seg = deepcopy(SegTemplate)
    
    # Modify various tags to match the values in Dicoms:
    """
    Note:
    
    The following changes are made in ModifySegTagVals():
        - a new SOP Instance UID is generated
        - a new Series Instance UID is generated
        - Series Description will be modified
        - Referenced SOP Instance UIDs in Referenced Series Sequence will be
        modified
        - Segment Label for the new (last) segment will be modified
        - Image Orientation Patient (in Shared Functional Groups Sequence) will
        be modified
    """
    Seg.StudyDate = Dicoms[0].StudyDate
    Seg.SeriesDate = Dicoms[0].SeriesDate
    Seg.ContentDate = Dicoms[0].ContentDate
    Seg.StudyTime = Dicoms[0].StudyTime
    Seg.SeriesTime = Dicoms[0].SeriesTime
    Seg.ContentTime = Dicoms[0].ContentTime
    
    # Modify the number of seguences in NewSeg:
    Seg = ChangeNumOfSegSequences(Seg=Seg, Dicoms=Dicoms, 
                                  NumOfFrames=len(PFFGStoSliceInds), 
                                  LogToConsole=LogToConsole)
    
    Seg.PatientName = Dicoms[0].PatientName
    Seg.PatientID = Dicoms[0].PatientID
    Seg.PatientBirthDate = Dicoms[0].PatientBirthDate
    Seg.PatientSex = Dicoms[0].PatientSex
    Seg.PatientAge = Dicoms[0].PatientAge
    Seg.StudyInstanceUID = Dicoms[0].StudyInstanceUID
    Seg.StudyID = Dicoms[0].StudyID
    Seg.FrameOfReferenceUID = Dicoms[0].FrameOfReferenceUID
    Seg.NumberOfFrames = f"{len(PFFGStoSliceInds)}"
    Seg.Rows = Dicoms[0].Rows
    Seg.Columns = Dicoms[0].Columns
    
    
    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .SliceThickness = f"{Image.GetSpacing()[2]}"
    
    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .SpacingBetweenSlices = f"{Image.GetSpacing()[2]}"

    Seg.SharedFunctionalGroupsSequence[0]\
       .PixelMeasuresSequence[0]\
       .PixelSpacing = Dicoms[0].PixelSpacing
       
    return Seg
    
    






def ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, NewTrgSegNums,  
                     NewTrgPFFGStoSliceInds, FromSegNum, FromSliceNum, SrcSeg,
                     LogToConsole=False):
    """
    Modify various tag values of the new SEG ROI so that they are consistent
    with the corresponding tag values of the DICOMs.
         
    
    Inputs:
    ------
                              
    TrgSeg : Pydicom object
        Original Target SEG ROI Object.
    
    TrgDicoms : List of Pydicom objects
        List of DICOM Objects that include the DICOMs that TrgSeg relate to.
    
    NewTrgPixArr : Numpy array
        New Target SEG pixel array.
    
    NewTrgSegNums : List of integers
        Segment numbers that correspond to each frame in NewTrgSeg.
    
    NewTrgPFFGStoSliceInds : List of integers
        Slice numbers that correspond to each frame in NewTrgSeg.
                           
    FromSegNum : integer
        Source segment number of the segmentation that was copied (counting  
        from 0) (this is only used for generating a new Series Description).
    
    FromSliceNum : integer
        Source slice index in the DICOM stack corresponding to the segmentation 
        that was copied (counting from 0) (this is only used for generating a 
        new Series Description).
                           
    SrcSeg : Pydicom object
        Source SEG ROI object (this is only used for generating a new Series 
        Description).
                         
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Output:
    ------
        
    NewTrgSeg : Pydicom object
        Modified Target SEG ROI object.
    """
    
    from copy import deepcopy
    from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    from pydicom.uid import generate_uid
    
    # Use TrgSeg as a template for NewTrgSeg:
    NewTrgSeg = deepcopy(TrgSeg)
    
    # The original and new number of frames in pixel array:
    OrigF = TrgSeg.pixel_array.shape[0]
    NewF = NewTrgPixArr.shape[0]
    
    # Increase the length of the Per-Frame Functional Groups Sequence in 
    # NewTrgSeg:
    for i in range(NewF - OrigF):
        NewTrgSeg = AddToPFFGS(NewTrgSeg)
            
    
    # Loop through each frame/sequence in NewTrgSeg, modifying the relevant tag 
    # values:
    for i in range(NewF):
        # The slice index for the i^th sequence:
        s = NewTrgPFFGStoSliceInds[i]
        
        # The new Target segment number:
        r = NewTrgSegNums[i]
        
        if LogToConsole:
            print(f'\nNewTrgPFFGStoSliceInds[{i}] = {NewTrgPFFGStoSliceInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Instance UID:
        NewTrgSeg.ReferencedSeriesSequence[0]\
                 .ReferencedInstanceSequence[i]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .DerivationImageSequence[0]\
                 .SourceImageSequence[0]\
                 .ReferencedSOPInstanceUID = TrgDicoms[s].SOPInstanceUID
        
        # Modify the Dimension Index Values using the RoiNum and SliceNum:
        """
        Note: r and s are integers starting from 0.  Need to add 1 to them so 
        that the integers in DIVs start from 1.
        """        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .FrameContentSequence[0]\
                 .DimensionIndexValues = [r + 1, s + 1]
                     
        # Modify the ImagePositionPatient: 
        IPP = deepcopy(TrgDicoms[s].ImagePositionPatient)
        
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .PlanePositionSequence[0]\
                 .ImagePositionPatient = IPP
                 
        # Modify the Referenced Segment Number (= RoiNum + 1):
        NewTrgSeg.PerFrameFunctionalGroupsSequence[i]\
                 .SegmentIdentificationSequence[0]\
                 .ReferencedSegmentNumber = r + 1
                 
        
        print(f'\nFrame number       = {i}')
        print(f'DICOM slice number = {s}')
        print(f'Segment number     = {r}')
        
        
    
    
    # The original and new number of segments:
    """
    Since the copied segmentation will be copied to a new segment NewS is 
    always OrigS + 1.
    """
    OrigS = len(TrgSeg.SegmentSequence)
    NewS = OrigS + 1
    
    # Increase the length of the Segment Sequence in NewTrgSeg:
    for i in range(NewS - OrigS):
        NewTrgSeg = AddToSS(NewTrgSeg)
        
        
        
    #print(f'The new Target pixel array has shape {NewTrgPixArr.shape}')
    
    # Ravel (to 1D array) NewTrgPixArr and pack bits 
    # (see https://github.com/pydicom/pydicom/issues/1230):
    packed = pack_bits(NewTrgPixArr.ravel())
    
    # Convert NewTrgPixArr to bytes:
    #NewSegRoi.PixelData = NewTrgPixArr.tobytes() <-- doesn't work
    NewTrgSeg.PixelData = packed + b'\x00' if len(packed) % 2 else packed
    
    # Modify the Number of Frames:
    NewTrgSeg.NumberOfFrames = f"{len(NewTrgPFFGStoSliceInds)}"
    
    # Generate a new SOP Instance UID:
    #NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
    NewTrgSeg.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    #NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
    NewTrgSeg.SeriesInstanceUID = generate_uid()
    
    # The segmentation was copied to a new segment with number:
    ToSegNum = NewTrgSegNums[-1] 
    """ verify the above """
    
    # Modify the Series Description:
    NewTrgSegSD = 'RPCopy_of_Seg' + str(FromSegNum) + '_s' + str(FromSliceNum) \
                  + '_in_' + SrcSeg.SeriesDescription + '_to_Seg' \
                  + str(ToSegNum) + '_in_' + TrgSeg.SeriesDescription
                  
    #NewTrgSegSD = 'Case3d_RPCopy'
    
    NewTrgSeg.SeriesDescription = NewTrgSegSD
    
    # Modify the Segment Label of the last (new) segment:
    FromSegLabel = deepcopy(NewTrgSeg.SegmentSequence[FromSegNum].SegmentLabel)
    
    NewSegLabel = 'NEW_' + FromSegLabel
    
    #print(f'\n--> FromSegNum = {FromSegNum}, ToSegNum = {ToSegNum}')
    #print(f'\n--> len(NewTrgSeg.SegmentSequence) = {len(NewTrgSeg.SegmentSequence)}')
    
    NewTrgSeg.SegmentSequence[ToSegNum].SegmentLabel = NewSegLabel
                                          
    return NewTrgSeg






def ExportSegRoi(NewTrgSeg, SrcSegFpath, TrgSegFpath, NamePrefix, FromSliceNum,
                 ExportDir):
    """
    Export contour ROI Object to disc.
    
    
    Inputs:
        NewTrgSeg      - New Target SEG ROI object to be exported
        
        SrcSegFpath    - (String) Full path of the Source DICOM SEG file
        
        TrgSegFpath    - (String) Full path of the Target DICOM SEG file
        
        NamePrefix     - (String) Prefix to be added to the assigned filename
                         after the DateTime stamp (e.g. 'Case3d_RPCopy') 
                            
        ExportDir      - (String) Directory where the SEG is to be exported
                           
    
                            
    Returns:
        NewTrgSegFpath - (String) Full path of the exported Target DICOM SEG 
                         file
    
    """
    
    import os
    import time

    # Modify the Series Description to a much shorter one (the one assigned in
    # ModifySegTagVals is very long):    
    NewTrgSeg.SeriesDescription = NamePrefix
    
    # The current DateTime:
    CDT = time.strftime("%Y%m%d_%H%M%S_", time.gmtime())
    
    # The filename of the Source and Target SEGs:
    SrcSegFname = os.path.split(SrcSegFpath)[1]
    TrgSegFname = os.path.split(TrgSegFpath)[1]
    
    NewTrgSegFname = CDT + NamePrefix + '_s' + str(FromSliceNum) + '_from_' \
                     + os.path.splitext(SrcSegFname)[0] + '_to_' + TrgSegFname
    
    NewTrgSegFpath = os.path.join(ExportDir, NewTrgSegFname)
    
    
    NewTrgSeg.save_as(NewTrgSegFpath)
        
    print('\nNew Target SEG ROI exported to:\n\n', NewTrgSegFpath)
    
    return NewTrgSegFpath

