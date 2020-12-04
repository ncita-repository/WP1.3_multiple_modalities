# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:30:31 2020

@author: ctorti
"""




def ModifyRtsTagVals(TrgRts, TrgDicoms, NewTrgContourData, NewTrgRoiNums,
                     NewTrgCStoSliceInds, FromRoiNum, FromSliceNum, SrcRts,
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
    
    NewTrgContourData : list of floats
        Flat list of coordinates for Target containing the new ContourData.
    
    NewTrgRoiNums : List of integers
        ROI numbers that correspond to each ROI in NewTrgRts.
    
    NewTrgCStoSliceInds : List of integers
        Slice numbers that correspond to each contour in NewTrgRts.
                           
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
    
    # Get the sequence number that corresponds to the slice number to be 
    # copied, and the slice number to be copied to:
    FromSeqNum = OrigCIStoDcmInds.index(FromSliceNum)
    ToSeqNum = NewCIStoDcmInds.index(ToSliceNum)
    
    if LogToConsole:
        print(f'\nFromOrigSeqNum = {FromOrigSeqNum} -->',
              f'OrigCIStoDcmInds[{FromOrigSeqNum}] =',
              f'{OrigCIStoDcmInds[FromOrigSeqNum]}')
        
        print(f'ToNewSeqNum   = {ToNewSeqNum} -->',
              f'NewCIStoDcmInds[{ToNewSeqNum}] =',
              f'{NewCIStoDcmInds[ToNewSeqNum]}')
    
    # Loop through each index in NewCIStoDcmInds and modify the values:
    for i in range(len(NewCIStoDcmInds)):
        # The DICOM slice index for the i^th sequence:
        s = NewCIStoDcmInds[i]
        
        if LogToConsole:
            print(f'\nNewCIStoDcmInds[{i}] = {NewCIStoDcmInds[i]}')
            print(f'DICOM slice number = {s}')
        
        # Modify the Referenced SOP Class UID and Referenced SOP Instance UID 
        # so that they are consistent with the SOP Class UID and SOP Instance  
        # UID of their corresponding DICOMs:
        """ SOPClassUIDs will need to be modified for cases where the contour
        is being copied across different modality images. Not required for
        copying within the same series:
        NewRtsRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewRts.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        """ Ditto:
        NewRts.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPClassUID = Dicoms[s].SOPClassUID
        """
        
        NewRts.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
                 
        
        
        # Modify the Contour Number to match the value of the Source contour
        #  being copied:
        NewRts.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourNumber = f"{i + 1}"
                     
                     
        # Modify select values in ROIContourSequence:         
        if i == ToNewSeqNum: 
            if LogToConsole:
                print('\nThis is a new sequence.')
                print(f'Copying the {FromOrigSeqNum}th sequence from Source',
                      f'to the {i}th sequence in Target.')
                
            # This is the new sequence. Copy from FromSeqNum.
            # Modify the Number of Contour Points to match the value of the 
            # Source contour being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[FromOrigSeqNum]\
                                  .NumberOfContourPoints)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[FromOrigSeqNum]\
                                  .ContourData)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourData = deepcopy(val)
        
        else:
            # This is an existing sequence. Find the index of OrigCIStoDcmInds 
            # for this sequence number (of NewCIStoDcmInds):
            ind = OrigCIStoDcmInds.index(NewCIStoDcmInds[i])
            
            if LogToConsole:
                print('\nThis is an existing sequence.')
                print(f'Copying the {ind}th sequence from Source to the {i}th',
                      f'sequence in Target.')
                                          
            # Modify the Number of Contour Points to match the value of the 
            # Source contour being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[ind]\
                                  .NumberOfContourPoints)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = deepcopy(val)
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            val = deepcopy(OrigRts.ROIContourSequence[0]\
                                  .ContourSequence[ind]\
                                  .ContourData)
            
            NewRts.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourData = deepcopy(val)
                  
                  
    # Generate a new SOP Instance UID:
    NewTrgRts.SOPInstanceUID = generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRts.SeriesInstanceUID = generate_uid()
                                          
    return NewRts









def VerifyCISandCS(CIStoDcmInds, CStoDcmInds, LogToConsole=False):
    """
    CIS and CS will only be the same if there is only one ROI so this function
    is of limited use.
    
    
    
    Verify that the ContourImageSequence-to-DICOM slice indeces and the
    ContourSequence-to-DICOM slice indeces are the same.
    
    Inputs:
        CIStoDcmInds - List of integers of the DICOM slice numbers that 
                      correspond to each Contour Image Sequence
        
        CStoDcmInds  - List of integers of the DICOM slice numbers that 
                      correspond to each Contour Sequence
        
    Returns:
        (Boolean) True/False if they match/do not match  
    """
    
    if CIStoDcmInds != CStoDcmInds:
        if LogToConsole:
            print('\nThe ContourImageSequence-to-DICOM slice indices are not',
                  'the same as the ContourSequence-to-DICOM slice indices.')
            
        return False
    
    else:
        return True
    
    





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







def GetPtsInContour_OLD(RtsFpath, FromRoiLabel, FromSliceNum, DicomDir):
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
    from DicomTools import GetDicomSOPuids
    
    # Import the RTSTRUCT ROI:
    Rts = dcmread(RtsFpath)
    
    RoiNum = GetRoiNum(Rts, FromRoiLabel)
        
    # The ROI Contour Sequence containing the contour to be copied:
    #Sequence = Rts.ROIContourSequence[RoiNum]
    
    # The Referenced SOP Instance UIDs in the Contour Sequence grouped by ROI:
    RefSOPsByRoi = GetRSOPuidsInCSByRoi(Rts)
    
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







def ChangeNumOfRtsSequences_OLD(Rts, CIStoSliceInds, CStoSliceIndsByRoi, 
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
    
    
    
    """
    Output some results to the console.
    """
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









def ModifyRtsTagVals_OLD(TrgRts, TrgDicoms, NewTrgPtsByRoi, NewTrgRoiNums,
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
            
            # The contour points for this ROI and this contour:
            Points = NewTrgPtsByRoi[r][c]
            
            # Modify NumberOfContourPoints:
            NewTrgRts.ROIContourSequence[r]\
                     .ContourSequence[c]\
                     .ContourImageSequence[0]\
                     .NumberOfContourPoints = f"{len(Points)}"
            
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








