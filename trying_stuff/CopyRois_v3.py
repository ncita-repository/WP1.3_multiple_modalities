# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:30:44 2020

@author: ctorti
"""





def CopyRois_v3(SourceDicomDir, SourceRoiFpath, TargetDicomDir, TargetRoiDir,
                CopyFromSliceNum, CopyToSliceNum, TargetStructureSetLabel,
                ExportRoi, LogToConsole):
    
    # Import packages and functions:
    import pydicom, copy, os, time
    from GetDicoms import GetDicoms
    from GetImageAttributes import GetImageAttributes
    
    
    
    
    
    """ 
    **************************************************************************
    Get the indices of the slices in the Source DICOMs that have contours. 
    **************************************************************************
    """
    
    # Get all DICOMs in SourceDicomDir:
    SourceDicomFpaths, SourceDicoms = GetDicoms(DirPath=SourceDicomDir, 
                                                SortMethod='slices', 
                                                Debug=False)
    
    # If the copying is across different series/scans:
    if TargetDicomDir != SourceDicomDir:
        # Get all DICOMs in TargetDicomDir:
        TargetDicomFpaths, TargetDicoms = GetDicoms(DirPath=TargetDicomDir, 
                                                    SortMethod='slices', 
                                                    Debug=False)
        
        # Get the Image Attributes for Source and Target (just for info):
        SourceOrigin, SourceDirs,\
        SourceSpacings, SourceDims = GetImageAttributes(DicomDir=SourceDicomDir, 
                                                        Package='pydicom')
        
        TargetOrigin, TargetDirs,\
        TargetSpacings, TargetDims = GetImageAttributes(DicomDir=TargetDicomDir, 
                                                        Package='pydicom')
        
        # Get the FrameOfReferenceUIDs:
        """ Will be used later.. """
        SourceFORuid = SourceDicoms[0].FrameOfReferenceUID
        TargetFORuid = TargetDicoms[0].FrameOfReferenceUID
        
        print(f'The Source scan has image dimensions    {SourceDims}')
        print(f'The Target scan has image dimensions    {TargetDims}')
        print(f'The Source scan has pixel spacings      {SourceSpacings}')
        print(f'The Target scan has pixel spacings      {TargetSpacings}')
        #print(f'The Source scan has patient position    {SourceOrigin}')
        #print(f'The Target scan has patient position    {TargetOrigin}')
        #print(f'The Source scan has patient orientation {SourceDirs}')
        #print(f'The Target scan has patient orientation {TargetDirs}')
        if SourceFORuid == TargetFORuid:
            print('The Source and Target scans have the same Frame Of Reference.')
        else:
            print('The Source and Target scans have different Frame Of References.')
            
    else:
        TargetDicoms = copy.deepcopy(SourceDicoms)
    
    
    # Get the Source ROI:
    SourceRoi = pydicom.dcmread(SourceRoiFpath)
    
    
    # Get the list of the Source DICOMs SOP Instance UIDs:
    SourceSOPuids = [] 
    
    for dicom in SourceDicoms:
        SourceSOPuids.append(dicom.SOPInstanceUID)
        
    
    
    # Loop through each Contour Image Sequence in SourceRoi, and find the index  
    # of the matching SOP Instance UID from the DICOMs:
    
    # Initialise the Contour Image Sequence to DICOM slice indices:
    SourceImageInds = [] 
    
    # Loop through each Contour Image Sequence:
    sequences = SourceRoi.ReferencedFrameOfReferenceSequence[0]\
                         .RTReferencedStudySequence[0]\
                         .RTReferencedSeriesSequence[0]\
                         .ContourImageSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SourceSOPuids:
        SourceImageInds.append(SourceSOPuids.index(uid))
        
    
    if LogToConsole:
        print('\nThere are contours for Contour Image Sequences:    ', 
              SourceImageInds)
          
    
    # Loop through each Contour Sequence and find the index of the 
    # matching SOP Instance UID from the DICOM:
    
    # Initialise the Contour Sequence to DICOM slice indices:
    SourceContourInds = [] 
    
    # Loop through each Contour Sequence:
    sequences = SourceRoi.ROIContourSequence[0].ContourSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SourceSOPuids:
        SourceContourInds.append(SourceSOPuids.index(uid))
    
    if LogToConsole:
        print('\nThere are contours for Contour Sequence Sequences: ', 
              SourceContourInds)
    
    
    
    """ 
    **************************************************************************
    Make sure that the index of the slice to be copied (CopyFromSliceNum) has  
    a contour in SourceRoi. 
    **************************************************************************
    """
    
    
    if LogToConsole:
        print('\nChecking if a contour exists for slice CopyFromSliceNum:', 
              CopyFromSliceNum)
    
    if not CopyFromSliceNum in SourceImageInds:
        print(f'\nA contour does not exist on slice number {CopyFromSliceNum}.')
        
        return
                
                
        
    """ 
    **************************************************************************
    Create new ROI object and modify it. 
    **************************************************************************
    """
    
    # Use SourceRoi as a template for TargetRoi: 
    TargetRoi = copy.deepcopy(SourceRoi)
    
    
    
    """ Add/remove Contour Image Sequences and Contour Sequences if there are 
    fewer/more sequences in SourceRoi than the number required in TargetRoi.
    
    Note:  For the moment there will always only be one contour copied over to
           TargetRoi.
    """
    
    # The number of contour image sequences in SourceRoi (= S) and the number
    # required in TargetRoi (= T = 1 for now):
    S = len(SourceImageInds)
    #T = len(CopySliceNums)
    T = 1
    
    if LogToConsole:       
        print(f'There are {S} contour sequences in SourceRoi.')
        print(f'There are {T} contour sequences required for TargetRoi.')
        
    # Return if T = 0 or T > 1:
    if (T == 0) or (T > 1):
        print('CopySliceNums must be a list containing a single integer value.')
        
        return
    
    
    # Add/remove contour image sequences and contour sequences by 
    # appending/popping the last sequence N times, where N is the difference
    # between S and T:
    if S != T:
        if T > S:
            if LogToConsole:       
                print(f'{T - S} more contour sequences are required in TargetRoi.')
                
            for i in range(T - S):
                TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                         .RTReferencedStudySequence[0]\
                         .RTReferencedSeriesSequence[0]\
                         .ContourImageSequence\
                         .append(copy.deepcopy(TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                                                        .RTReferencedStudySequence[0]\
                                                        .RTReferencedSeriesSequence[0]\
                                                        .ContourImageSequence[-1]))

                TargetRoi.ROIContourSequence[0].ContourSequence\
                         .append(copy.deepcopy(TargetRoi.ROIContourSequence[0]\
                                                        .ContourSequence[-1]))
                
        else:
            if LogToConsole:       
                print(f'There are {S - T} excess contour sequences in SourceRoi.')
                
            for i in range(S - T):
                TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                         .RTReferencedStudySequence[0]\
                         .RTReferencedSeriesSequence[0]\
                         .ContourImageSequence.pop()
                
                TargetRoi.ROIContourSequence[0].ContourSequence.pop()
                

    T = len(TargetRoi.ROIContourSequence[0].ContourSequence)
    
    if LogToConsole:       
        print(f'There are now {T} contour sequences in TargetRoi.')
        
        
        
    
    
    
    
    """ Modify TargetRoi. """

    # Generate a new SOP Instance UID:
    TargetRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    TargetRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    """ 
    Some tags will not need modifying if copying is within the same 
    session/scan. 
    """
    # If the copying is across different scans:
    if TargetDicomDir != SourceDicomDir:
    
        # Modify the Study Date and Time:
        TargetRoi.StudyDate = TargetDicoms[0].StudyDate
        """ 
        For some data, the Series Time and Instance Creation Time in the DICOM 
        both match the Study Time in the ROI, whereas the DICOM's Study Time 
        does not.
        """
        #TargetRoi.StudyTime = TargetDicoms[0].StudyTime 
        TargetRoi.StudyTime = TargetDicoms[0].SeriesTime
        
        # Modify the Patient's Name and ID:
        TargetRoi.PatientName = TargetDicoms[0].PatientName
        
        TargetRoi.PatientID = TargetDicoms[0].PatientID
        
        """ Skipping other non-critical tags... """
        
        # Modify the Study Instance UID:
        TargetRoi.StudyInstanceUID = TargetDicoms[0].StudyInstanceUID
        
        # Modify the Frame of Reference UID:
        TargetRoi.FrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
        
        # Modify the Referenced Frame of Reference Sequence:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .FrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
        
        # Modify the Referenced SOP Instance UID in the RT Referenced Study 
        # Sequence:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .ReferencedSOPInstanceUID = TargetDicoms[0].StudyInstanceUID
        
        # Modify the Series Instance UID in the RT Referenced Series Sequence:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .SeriesInstanceUID = TargetDicoms[0].SeriesInstanceUID
        
        # Modify the Referenced Frame of Reference UID for the Structure Set  
        # ROI Sequence:
        TargetRoi.StructureSetROISequence[0]\
                 .ReferencedFrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
    
        
        # Modify the RT ROI Observations Sequence:
        """ Not sure if these need to be changed """
        #TargetRoi.RTROIObservationsSequence[0].ObservationNumber = 
        #TargetRoi.RTROIObservationsSequence[0].ReferencedROINumber = 
        #TargetRoi.RTROIObservationsSequence[0].RTROIInterpretedType = 
        #TargetRoi.RTROIObservationsSequence[0].ROIInterpreter = 
    
    
    # Modify the Structure Set Label (this appears under Name in XNAT):
    #TargetRoi.StructureSetLabel = 'Copy_of_' + SourceRoi.StructureSetLabel
    TargetRoi.StructureSetLabel = TargetStructureSetLabel
    
    # Modify the ROI Name for the Structure Set ROI Sequence:             
    TargetRoi.StructureSetROISequence[0]\
             .ROIName = TargetRoi.StructureSetLabel
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    TargetRoi.StructureSetDate = NewDate
    TargetRoi.StructureSetTime = NewTime
    
    
    # Simplify variable names for indexing.
    # The slice index to be copied:
    #s = CopyFromSliceNum
    
    # The slice index to be copied to:
    #t = CopyToSliceNum
    
    # Modify the Referenced SOP Class UID, the SOP Instance UID and the
    # Referenced SOP Class UID: 
    
    # Cycle through each sequence in TargetRoi:
    """ For now T = 1. """
    for i in range(T):
    
        # Modify the Referenced SOP Class UID to match the SOP Class UID of the
        # Target DICOM the ROI is to be copied to:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPClassUID = TargetDicoms[CopyToSliceNum].SOPClassUID
        
        # Modify the SOP Instance UID to match the SOP Instance UID of the
        # Target DICOM the ROI is to be copied to:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i]\
                 .ReferencedSOPInstanceUID = TargetDicoms[CopyToSliceNum].SOPInstanceUID
        
        # Modify the Referenced SOP Class UID to match the SOP Class UID of the
        # Target DICOM the ROI is to be copied to:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPClassUID = TargetDicoms[CopyToSliceNum].SOPClassUID
        
        # Modify the SOP Instance UID to match the SOP Instance UID of the 
        # Target DICOM the ROI is to be copied to:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPInstanceUID = TargetDicoms[CopyToSliceNum].SOPInstanceUID
        
        
        #print(f'\nTargetRoi.ROIContourSequence[0].ContourSequence[i].NumberOfContourPoints = {TargetRoi.ROIContourSequence[0].ContourSequence[i].NumberOfContourPoints}')
        
        #print(f's = {s}')
        
        #print(f'\nSourceRoi.ROIContourSequence[0].ContourSequence[s].NumberOfContourPoints = {SourceRoi.ROIContourSequence[0].ContourSequence[s].NumberOfContourPoints}')
        
        # Get the sequence number that corresponds to the slice number to be
        # copied:
        CopyFromSequenceNum = SourceImageInds.index(CopyFromSliceNum)
        
        # Modify the Number of Contour Points to match the value of the Source
        # contour being copied:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .NumberOfContourPoints = SourceRoi.ROIContourSequence[0]\
                                                   .ContourSequence[CopyFromSequenceNum]\
                                                   .NumberOfContourPoints
        
        # Modify the Contour Number to match the value of the Source contour
        # being copied:
        #TargetRoi.ROIContourSequence[0]\
        #         .ContourSequence[i]\
        #         .ContourNumber = SourceRoi.ROIContourSequence[0]\
        #                                    .ContourSequence[CopyFromSequenceNum]\
        #                                    .ContourNumber
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourNumber = f"{i + 1}"
        
        # Modify the Contour Data to match the value of the Source contour
        # being copied:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[i]\
                 .ContourData = SourceRoi.ROIContourSequence[0]\
                                         .ContourSequence[CopyFromSequenceNum]\
                                         .ContourData
        
           
    # Get the filename of the original RT-Struct file:
    SourceRoiFname = os.path.split(SourceRoiFpath)[1]
    
    # Create a new filename (this will appear under Label in XNAT):
    TargetRoiFname = 'AIM_' + NewDate + '_' + NewTime + '_copied_from_' \
                     + SourceRoiFname
    
    TargetRoiFpath = os.path.join(TargetRoiDir, TargetRoiFname)
    
    # Save the new DICOM-RTSTRUCT file:
    if ExportRoi:
        TargetRoi.save_as(TargetRoiFpath)
        
        print('\nROI exported to:\n', TargetRoiFpath)
    
    return TargetRoi




"""
    
    ----------------------------------------------------------------------------------------------
    Replace in TargetROI:                                      The value of the tag in the DICOM:
    ----------------------------------------------------------------------------------------------
    Study Instance UID (0020, 000d)                            Study Instance UID (0020, 000d)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Study Instance UID (0020, 000d)
    -> RT Referenced Study Sequence[0] (3006, 0012)
      -> Referenced SOP Instance UID (0008, 1155)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Series Instance UID (0020, 000e)
    -> RT Referenced Study Sequence[0] (3006, 0012)
      -> RT Referenced Series Sequence[0] (3006, 0014)
        -> Series Instance UID (0020, 000e)
    
    Frame of Reference UID (0020, 0052)                        Frame of Reference UID (0020, 0052)   
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Frame of Reference UID (0020, 0052)   
    -> Frame of Reference UID (0020, 0052)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Frame of Reference UID (0020, 0052)
    -> Structure Set ROI Sequence[0] (3006, 0020)
      -> Referenced Frame of Reference UID (3006, 0024)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     SOP Instance UID (0008, 0018)
    -> RT Referenced Study Sequence[0] (3006, 0012) 
      -> RT Referenced Series Sequence[0] (3006, 0014)
        -> Contour Image Sequence[i] (3006, 0016)
          -> Referenced SOP Instance UID (0008, 1155)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     SOP Instance UID (0008, 0018)
    -> ROI Contour Sequence[0] (3006, 0039)
      -> Contour Sequence[i] (3006, 0040)
        -> Contour Image Sequence[0] (3006, 0016)
          -> Referenced SOP Instance UID (0008, 1155)
"""