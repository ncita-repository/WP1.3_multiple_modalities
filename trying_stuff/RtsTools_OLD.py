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








def ExportRts(TrgRts, SrcRtsFpath, NamePrefix, ExportDir):
    """
    Export RTS to disk.  
    
    Inputs:
    ------
    
    TrgRts : Pydicom object
        Target RT-STRUCT ROI object to be exported.
        
    SrcRtsFpath : string
        Full path of the Source DICOM-RTSTRUCT file (used to generate the 
        filename of the new RTS file).
        
    NamePrefix : string
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
                            
    ExportDir : string
        Directory where the RTSTRUCT is to be exported.
                           
    
    Outputs:
    -------
    
    TrgRtsFpath : string
        Full path of the exported Target DICOM-RTSTRUCT file.
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
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    FnamePrefix = CurrentDateTime + '_' + NamePrefix + '_from_' 
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




def GetPtsByRoi_OLD(Rts, DicomDir):
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
    
    # Get the DICOM SOPInstanceUIDs:
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
        NumC2Sinds = len(CStoSliceIndsByRoi[r])
        
        # Get a list of the ContourSequences in this ROI:
        ContourSequences = deepcopy(Rts.ROIContourSequence[r].ContourSequence)
        
        NumCS = len(ContourSequences)
        
        if NumC2Sinds != NumCS:
            msg = f'The number of contour-sequence-to-slice-indices '\
                  + f'({NumC2Sinds}) must be equal to the number of '\
                  + f'sequences in ROIContourSequence[{r}].ContourSequence '\
                  + f'({NumCS}).'
            
            raise Exception(msg)
        
        # Get a list of all ReferencedSOPInstanceUIDs from ContourSequences:
        #RefSopUids = [sequence.ContourImageSequence[0]\
        #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        if False:
            print(f'   r = {r}')
            print(f'   NumC2Sinds = {NumC2Sinds}')
            print(f'   NumCS = NumCS')
        
        
        PtsByContour = []
        
        
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
    
    PtsByCnt : list of a list of a list of floats
        List (for each contour) of a list (for each point) of a list (for each 
        dimension) of coordinates. 
        
    C2Sinds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        list of points in PointsByContour.
    """
    
    from DicomTools import GetDicomSOPuids
    from DicomTools import GetRoiNum
    
    RoiContourSequence = Rts.ROIContourSequence
    
    NumOfRois = len(RoiContourSequence)
    
    RoiNum = GetRoiNum(Rts, SearchString)
    
    if RoiNum > NumOfRois - 1:
        raise Exception(f"There are only {NumOfRois} ROIs in the RTS.\n",
                        f"RoiNum = {RoiNum} > {NumOfRois}.")
    
    
    # The ContourSequences in the ROI of interest:
    ContourSequence = Rts.ROIContourSequence[RoiNum].ContourSequence
        
    CS = len(ContourSequence)
    
    PtsByCnt = []
    
    for i in range(CS):
        ContourData = ContourSequence[i].ContourData
        
        ContourData = [float(item) for item in ContourData]
        
        Points = []
                
        # Iterate for all points in threes:
        for p in range(0, len(ContourData), 3):
            point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
            
            Points.append(point)
            
        PtsByCnt.append(Points)
    
    
    SOPuids = GetDicomSOPuids(DicomDir)
    
    # The ContourSequence-to-slice indices by ROI:
    C2SindsByRoi = GetCStoSliceIndsByRoi(Rts, SOPuids)
    
    C2Sinds = C2SindsByRoi[RoiNum]
        
    return PtsByCnt, C2Sinds








def GetPtsInContour_OLD(Rts, SearchString, SliceNum, DicomDir):
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
    
    RoiContourSequence = deepcopy(Rts.ROIContourSequence)
    
    NumOfRois = len(RoiContourSequence)
    
    RoiNum = GetRoiNum(Rts, SearchString)
    
    if RoiNum > NumOfRois - 1:
        msg = f'There are only {NumOfRois} ROIs in the RTS.\n'\
              + f'RoiNum = {RoiNum} > {NumOfRois}.'
              
        raise Exception(msg)
    
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourSequence-to-slice indices by ROI:
    CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(Rts, SopUids)
    
    CStoSliceInds = CStoSliceIndsByRoi[RoiNum]
    
    if not SliceNum in CStoSliceInds:
        raise Exception(f'There is no contour on slice {SliceNum} in the RTS.')
        
    # The contour number of interest:
    ContourNum = CStoSliceInds.index(SliceNum)
    
    # The ContourSequences in the ROI of interest:
    #ContourSequence = deepcopy(Rts.ROIContourSequence[RoiNum].ContourSequence)
    ContourSequence = deepcopy(RoiContourSequence[RoiNum].ContourSequence)
    
    # The number of contours in this ROI:
    NumOfContours = len(ContourSequence)
    
    if ContourNum > NumOfContours - 1:
        msg = f'There are only {NumOfContours} contours in the RTS.\n'\
              + f'ContourNum = {ContourNum} > {NumOfContours}.'
                        
        raise Exception(msg)
        
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








def GetMaskFromRts_OLD(Rts, SearchString, SliceNum, DicomDir, RefImage):
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
def InitialiseRts_OLD(RtsTemplate, RoiNum, CStoSliceInds, DicomDir, 
                  NamePrefix='', LogToConsole=False):
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
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
    
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
    from DicomTools import ImportDicom
    
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






def ModifyRts_OLD(Rts, CntDataByCnt, PtsByCnt, CStoSliceInds, DicomDir, 
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
        
        # Modify the ReferencedSOPInstanceUID in ContourImageSequence:
        Rts.ReferencedFrameOfReferenceSequence[0]\
           .RTReferencedStudySequence[0]\
           .RTReferencedSeriesSequence[0]\
           .ContourImageSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # Modify the ReferencedSOPInstanceUID in ROIContourSequence:
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
           .ContourImageSequence[0]\
           .ReferencedSOPInstanceUID = Dicoms[s].SOPInstanceUID
        
        # The contour points for this contour:
        Points = PtsByCnt[i]
        
        # Modify NumberOfContourPoints:
        """
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
           .ContourImageSequence[0]\    <--- error found 14/12/2020
           .NumberOfContourPoints = f"{len(Points)}"
        """
        Rts.ROIContourSequence[0]\
           .ContourSequence[i]\
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








def GetContoursFromListOfRtss_OLD(ListOfRtss, ListOfDicomDirs, SearchString, 
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
    #from RtsTools import GetPtsInRoi
    from ImageTools import GetImageAttributes
    
    ListOfContours = []
    ListOfC2Sinds = []
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
            
            Contours, C2Sinds = GetPtsInRoi(ListOfRtss[i], 
                                            ListOfDicomDirs[i], 
                                            SearchString)
            
            if LogToConsole:      
                print(f'len(Contours) in "{SearchString}" = {len(Contours)}')
                print(f'C2Sinds = {C2Sinds}')
                [print(f'len(Contours{i}) = {len(Contours[i])}') for i in range(len(Contours))]
            
        else:
            RoiNum = None
            Contours = None
            C2Sinds = None
            
            if LogToConsole:      
                print(f'No contours for this dataset')
        
        
        ListOfRoiNums.append(RoiNum)
        
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(ListOfDicomDirs[i])
        
        ListOfContours.append(Contours)
        ListOfC2Sinds.append(C2Sinds)
        ListOfDcmFpaths.append(DcmFpaths)
        ListOfOrigins.append(IPPs[0])
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
        
        if LogToConsole:      
            #print(f'len(Contours) in "{SearchString}" = {len(Contours)}')
            print(f'len(ListOfContours) = {len(ListOfContours)}')
            #print(f'C2Sinds = {C2Sinds}')
            #[print(f'len(Contours{i}) = {len(Contours[i])}') for i in range(len(Contours))]
        
        
        #SOPuids = GetDicomSOPuids(ListOfDicomDirs[i])
        
        #CStoSliceIndsByRoi = GetCStoSliceIndsByRoi(ListOfRtss[i], SOPuids)
        
        #ListOfC2Sinds.append(CStoSliceIndsByRoi[RoiNum])
        
        #if LogToConsole:
        #    print(f'ListOfC2Sinds[{i}] = {ListOfC2Sinds[i]}')
    
            
    return ListOfContours, ListOfC2Sinds, ListOfRoiNums, ListOfDcmFpaths,\
           ListOfOrigins, ListOfDirections, ListOfSpacings






def InitialiseRts_OLD(RtsTemplate, RoiNum, CStoSliceInds, CntDataByObjByFrame,
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









def ModifyRts_OLD(Rts, CntDataByObjByFrame, PtsByObjByFrame, CStoSliceInds, 
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







def ErrorCheckRts_OLD(Rts, DicomDir, LogToConsole=False):
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
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=DicomDir, Package='sitk')
    
    
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





def GetRtsDataOfInterest_OLD(Rts, FromSliceNum, FromRoiLabel, DicomDir, 
                         LogToConsole=False):
    """
    Get data of interest from an RTS.
    
    Inputs:
    ******
    
    Rts : Pydicom object
        ROI object from an RTSTRUCT file.
        
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
        
    FromRoiLabel : string
        All or part of the Source ROI Name of the ROI containing the contour(s)
        to be copied.
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates to be copied from 
        Rts.
        
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour to be copied from Rts.
        
        
    Notes:
    *****
    
    There are 3 possible main use cases:
        1. Copy a single contour
        2. Copy an entire ROI
        3. Copy all ROIs
        
    Which main case applies depends on the inputs:
        FromSliceNum
        FromRoiLabel
    
    If FromSliceNum != None (i.e. if a slice number is defined) a single 
    contour will be copied from the ROI that matches FromRoiLabel (--> Main Use
    Case 1).  
    
    Need to decide what to do if there are multiple ROIs and 
    FromRoiLabel = None.  Possible outcomes:

        - Raise exception with message "There are multiple ROIs so the ROI
        label must be provided"
        - Copy the contour to all ROIs that have a contour on slice 
        FromSliceNum (* preferred option?)
        
    If FromSliceNum = None but FromRoiLabel != None, all contours within the 
    ROI given by FromRoiLabel will be copied (--> Main Use Case 2).  
    
    If FromSliceNum = None and FromRoiLabel = None, all contours within all
    ROIs will be copied (--> Main Use Case 3).
    
    
    If FromSliceNum = None and FromRoiLabel = None, all ROIs will be copied.  
    If not, get the ROI
    number (zero-indexed) to be copied, and reduce SrcPtsByCntByRoi and
    SrcC2SindsByRoi to the corresponding item.
    """
    
    from DicomTools import GetRoiNums, GetRoiLabels
    from copy import deepcopy
    
    AllPtsByCntByRoi, AllC2SindsByRoi = GetPtsByCntByRoi(Rts, DicomDir, 
                                                         LogToConsole)
    
    FromRoiLabels = GetRoiLabels(Roi=Rts)
    
    if FromRoiLabel:
        """ Limit data to those which belong to the chosen ROI(s). """
        
        PtsByCntByRoi = []
        C2SindsByRoi = []
        
        # Get the ROI number(s) whose name matches FromRoiLabel:
        FromRoiNums = GetRoiNums(Roi=Rts, SearchString=FromRoiLabel)
        
        FromRoiLabels = [FromRoiLabels[i] for i in FromRoiNums]
        
        for r in range(len(FromRoiNums)):
            # All points-by-contour and contour-to-slice indices for this ROI:
            AllPtsByCnt = deepcopy(AllPtsByCntByRoi[FromRoiNums[r]])
            AllC2Sinds = deepcopy(AllC2SindsByRoi[FromRoiNums[r]])
                
            if FromSliceNum:
                """ Limit data to those which belong to the chosen slice 
                number. """
                
                # Get the contour number(s) that relate to this slice (there 
                # may be more than one contour on FromSliceNum):
                FromCntNums = [i for i, e in enumerate(AllC2Sinds) if e==FromSliceNum]
                
                if FromCntNums:
                    # Keep only the indeces and points that relate to 
                    # FromCntNums:
                    PtsByCnt = [AllPtsByCnt[c] for c in FromCntNums]
                    C2Sinds = [AllC2Sinds[c] for c in FromCntNums]
                else:
                    PtsByCnt = []
                    C2Sinds = []
                
                PtsByCntByRoi.append(PtsByCnt)
                C2SindsByRoi.append(C2Sinds)
                
                if LogToConsole:
                    print(f'\nThe {len(FromCntNums)} contour(s) on slice',
                          f'{FromSliceNum} will be copied from ROI(s)',
                          f'{FromRoiLabels}.')
                    
            else:
                PtsByCntByRoi.append(AllPtsByCnt)
                C2SindsByRoi.append(AllC2Sinds)
                
                if LogToConsole:
                    print(f'\nAll {len(AllC2Sinds)} contour(s) in ROI(s)',
                          f'{FromRoiLabels} will be copied.')
    
    else:
        """ No limit will be put on the ROIs. """
        
        if FromSliceNum:
            """ Limit data to those which belong to the chosen slice number. """
            
            PtsByCntByRoi = []
            C2SindsByRoi = []
            
            for r in range(len(FromRoiLabels)):
                # All points-by-contour and contour-to-slice indices for this
                # ROI:
                AllPtsByCnt = deepcopy(AllPtsByCntByRoi[r])
                AllC2Sinds = deepcopy(AllC2SindsByRoi[r])
                
                # Get the contour number(s) that relate to this slice (there 
                # may be more than one contour on FromSliceNum) in all ROIs:
                FromCntNums = [[i for i, e in enumerate(AllC2Sinds) if e==FromSliceNum] for r in range(len(C2SindsByRoi))]
                
                PtsByCnt = []
                C2Sinds = []
                
                if FromCntNums:
                    for c in range(len(FromCntNums)):
                        i = FromCntNums[c]
                        
                        PtsByCnt.append(AllPtsByCnt[i])
                        C2Sinds.append(AllC2Sinds[i])
                        
                PtsByCntByRoi.append(PtsByCnt)
                C2SindsByRoi.append(C2Sinds)
                
                if LogToConsole:
                    print(f'\nAll contour(s) on slice {FromSliceNum} will be',
                          f'copied from all ROIs.')
                
        else:
            PtsByCntByRoi = deepcopy(AllPtsByCntByRoi)
            C2SindsByRoi = deepcopy(AllC2SindsByRoi)
            
            if LogToConsole:
                print('\nAll contour(s) on all slices in all ROIs will be',
                      'copied.')
    
    if LogToConsole:
        print(f'\nC2SindsByCntByRoi = {C2SindsByRoi}')
        
    return PtsByCntByRoi, C2SindsByRoi







def GetRtsDataFromListOfRtss_OLD(ListOfRtss, ListOfDicomDirs, SearchString='', 
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
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
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





def CreateRts_OLD(SrcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi, TrgPtsByCntByRoi, # 26/03/2021
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
    import time
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
    
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
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
    
    """ Initialise the number of ROIs with non-zero TrgC2Sinds: """
    R = 0
    
    """ Loop through each TrgC2SindsByRoi: """
    for r in range(len(TrgC2SindsByRoi)):
        """ The number of contours in this ROI: """
        C = len(TrgC2SindsByRoi[r])
        
        if C == 0:
            """ There are no contours for this ROI, so pop this sequence: """
            NewTrgRts.StructureSetROISequence.pop(r)
        else:
            R += 1
    
    
    """ Modify the ROINumbers, ReferencedFrameOfReferenceUIDs and ROINames: """
    for s in range(len(NewTrgRts.StructureSetROISequence)):
        NewTrgRts.StructureSetROISequence[s].ROINumber = f"{s+1}"
        
        NewTrgRts.StructureSetROISequence[s]\
                 .ReferencedFrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
           
        NewTrgRts.StructureSetROISequence[s]\
                 .ROIName = NewTrgRts.StructureSetROISequence[s]\
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
        
        
        #""" The number of ContourSequences: """
        #N = len(NewTrgRts.ROIContourSequence[r].ContourSequence)
        
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
    
    """ Initialise the number of ROIs with non-zero TrgC2Sinds: """
    R = 0
    
    """ Loop through each TrgC2SindsByRoi: """
    for r in range(len(TrgC2SindsByRoi)):
        """ The number of contours in this ROI: """
        C = len(TrgC2SindsByRoi[r])
        
        if C == 0:
            """ There are no contours for this ROI, so pop this sequence: """
            NewTrgRts.RTROIObservationsSequence.pop(r)
        else:
            R += 1
    
    
    """ Modify the ObservationNumbers and ReferencedROINumbers: """
    for s in range(len(NewTrgRts.RTROIObservationsSequence)):
        NewTrgRts.RTROIObservationsSequence[s].ObservationNumber = f"{s+1}"
        
        NewTrgRts.RTROIObservationsSequence[s].ReferencedROINumber = f"{s+1}"
           
    
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






