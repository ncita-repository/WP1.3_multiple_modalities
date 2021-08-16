# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 09:14:37 2020

@author: ctorti
"""





def CopyRois_v2(SourceDicomDir, SourceRoiFpath, TargetDicomDir, TargetRoiDir,
                CopySliceNums, LogToConsole):
    
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
    Make sure that the index of the slice to be copied (CopySliceNums) has a 
    contour. 
    
    For now CopySliceNums must be a list containing only one integer.
    **************************************************************************
    """
    
    
    if LogToConsole:
        print('\nChecking if a contour exists for slices CopySliceNums:', 
              CopySliceNums[0])
    
    if not CopySliceNums[0] in SourceImageInds:
        print(f'\nContours do not exist on slice number {CopySliceNums[0]}.')
        
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
    T = len(CopySliceNums)
    
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
    
    # Modify the Study Date and Time:
    TargetRoi.StudyDate = TargetDicoms[0].StudyDate
    """ 
    For current data, the Series Time and Instance Creation Time in the DICOM 
    both match the Study Time in the ROI, whereas the DICOM's Study Time does not.
    """
    #TargetRoi.StudyTime = TargetDicoms[0].StudyTime 
    TargetRoi.StudyTime = TargetDicoms[0].SeriesTime
    
    # Modify the Patient's Name and ID:
    TargetRoi.PatientName = TargetDicoms[0].PatientName
    TargetRoi.PatientID = TargetDicoms[0].PatientID
    
    """ Skipping other non-critical tags... """
    
    # Modify the Study Instance UID:
    TargetRoi.StudyInstanceUID = TargetDicoms[0].StudyInstanceUID
    
    # Generate a new Series Instance UID:
    TargetRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    # Modify the Frame of Reference UID:
    TargetRoi.FrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
    
    # Modify the Structure Set Label:
    TargetRoi.StructureSetLabel = 'Copy_of_' + TargetRoi.StructureSetLabel
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    TargetRoi.StructureSetDate = NewDate
    TargetRoi.StructureSetTime = NewTime
    
    # Modify the Referenced Frame of Reference Sequence:
    TargetRoi.ReferencedFrameOfReferenceSequence[0]\
             .FrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
    
    # Modify the Referenced SOP Instance UID in the RT Referenced Study Sequence:
    TargetRoi.ReferencedFrameOfReferenceSequence[0]\
             .RTReferencedStudySequence[0]\
             .ReferencedSOPInstanceUID = TargetDicoms[0].StudyInstanceUID
    
    # Modify the Series Instance UID in the RT Referenced Series Sequence:
    TargetRoi.ReferencedFrameOfReferenceSequence[0]\
             .RTReferencedStudySequence[0]\
             .RTReferencedSeriesSequence[0]\
             .SeriesInstanceUID = TargetDicoms[0].SeriesInstanceUID
    
    # Modify the Referenced Frame of Reference UID and ROI Name for the 
    # Structure Set ROI Sequence:
    TargetRoi.StructureSetROISequence[0]\
             .ReferencedFrameOfReferenceUID = TargetDicoms[0].FrameOfReferenceUID
    TargetRoi.StructureSetROISequence[0]\
             .ROIName = TargetRoi.StructureSetLabel
    
    # Modify the RT ROI Observations Sequence:
    """ Not sure if these need to be changed """
    #TargetRoi.RTROIObservationsSequence[0].ObservationNumber = 
    #TargetRoi.RTROIObservationsSequence[0].ReferencedROINumber = 
    #TargetRoi.RTROIObservationsSequence[0].RTROIInterpretedType = 
    #TargetRoi.RTROIObservationsSequence[0].ROIInterpreter = 
    
    
    # Modify the Referenced SOP Class UID, the SOP Instance UID and the
    # Referenced SOP Class UID: 
    
    # Cycle through each sequence in TargetRoi:
    for t in range(T):
        # The Source index to of the slice to be copied from is:
        s = CopySliceNums[t]
        
        if LogToConsole:
            print(f't = {t}, CopySliceNums[{i}] = {s}')
        
        # Modify the Referenced SOP Class UID:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[t]\
                 .ReferencedSOPClassUID = TargetDicoms[s].SOPClassUID
        
        # Modify the SOP Instance UID:
        TargetRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0]\
                 .RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[t]\
                 .ReferencedSOPInstanceUID = TargetDicoms[s].SOPInstanceUID
        
        # Modify the Referenced SOP Class UID:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[t]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPClassUID = TargetDicoms[s].SOPClassUID
        
        # Modify the SOP Instance UID:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[t]\
                 .ContourImageSequence[0]\
                 .ReferencedSOPInstanceUID = TargetDicoms[s].SOPInstanceUID
        
        
        if False:
            # Modify the Contour Geometric Type:
            """ Since the contours are not closed change from 'CLOSED_PLANAR' to 
            'OPEN_PLANAR' """
            #print('\nBefore:', TargetRoi.ROIContourSequence[0]\
            #                            .ContourSequence[s]\
            #                            .ContourImageSequence[0]\
            #                            .ReferencedSOPClassUID)
            TargetRoi.ROIContourSequence[0]\
                     .ContourSequence[t]\
                     .ContourGeometricType = 'OPEN_PLANAR'
            #print('\nAfter:', TargetRoi.ROIContourSequence[0]\
            #                           .ContourSequence[s]\
            #                           .ContourGeometricType)
            
            #print(TargetRoi.ROIContourSequence[0]\
            #               .ContourSequence[s].ContourGeometricType)
        
        # Modify the Number of Contour Points:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[t]\
                 .NumberOfContourPoints = SourceRoi.ROIContourSequence[0]\
                                                   .ContourSequence[s]\
                                                   .NumberOfContourPoints
        
        # Modify the Contour Number:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[t]\
                 .ContourNumber = SourceRoi.ROIContourSequence[0]\
                                            .ContourSequence[s]\
                                            .ContourNumber 
        
        # Modify the Contour Data:
        TargetRoi.ROIContourSequence[0]\
                 .ContourSequence[t]\
                 .ContourData = SourceRoi.ROIContourSequence[0]\
                                         .ContourSequence[s]\
                                         .ContourData
        
           
    # Get the filename of the original RT-Struct file:
    SourceRoiFname = os.path.split(SourceRoiFpath)[1]
    
    # Create a new filename:
    TargetRoiFname = 'AIM_' + NewDate + '_' + NewTime + '_transformed_from_' \
                     + SourceRoiFname
    
    TargetRoiFpath = os.path.join(TargetRoiDir, TargetRoiFname)
    
    # Save the new DICOM-RTSTRUCT file:
    #TargetRoi.save_as(TargetRoiFpath)
    
    #print('\nROI for Target scan exported to:\n', TargetRoiFpath)
    
    return TargetRoi