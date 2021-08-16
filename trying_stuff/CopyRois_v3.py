# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:30:44 2020

@author: ctorti
"""





def CopyRois_v3(SrcDcmDir, SrcRoiFpath, TrgDcmDir, TrgRoiDir,
                CopyFromSliceNum, CopyToSliceNum, CopyType,
                TrgSSLabel, ExportRoi, LogToConsole):
    """
    Create a new ROI Collection (exported as a DICOM-RTSTRUCT file) for a
    Target DICOM stack, having copied a contour from the Source ROI Collection 
    that corresponds to a particular slice in the Source DICOM stack.
    
    How the contour is copied depends on the following cases:
        1. The Source and Target DICOMs have the same Frame of Reference.
        
        2. The Source and Target DICOMs have different Frames of Reference and
           the contour is to be copied literally.
           
        3. The Source and Target DICOMs have different Frames of Reference and
           the contour is to be copied exactly.
    
    
    Cases 1 and 2 are similar in that the contour at a particular slice in the 
    Source DICOM stack (CopyFromSliceNum) is copied to a particular slice in
    the Target DICOM stack (CopyToSliceNum).  
    
    For Case 2, the Source and Target image stacks are registered (since they 
    have different Frames of Reference).  The points are copied in a literal 
    sense since the 3D location of the transformed points is ignored, and the 
    transformed points will generally not lie in the Target image planes.
    
    Note:
        I'm not certain how the OHIF Viewer will interpret contour points that
        do not lie in the actual planes.  This remains to be seen...
        
    For Case 3, the Source and Target image stacks are also registered, but the
    the locations where the transformed points project onto (intersect with) 
    the Target image planes is determined, so that the resulting points lie on
    the planes.
    
    Cases 1 and 2 are also similar in that the contour that is copied will not
    necessarily represent the physical structure of the ROI. Whereas in Case 3
    the physicality of the contour points is preserved.
    
    
    Inputs:
        SrcDcmDir        - (String) Directory containing the Source DICOMs.
                            
        SrcRoiFpath      - (String) Full path of the Source DICOM-RTSTRUCT file.
        
        TrgDcmDir        - (String) Directory containing the Target DICOMs.
        
        TrgRoiDir        - (String) Directory where the Target DICOM-RTSTRUCT 
                           file will be exported.
                            
                            
                             
        CopyFromSliceNum - (Integer) Slice index to be copied from the Source
                           DICOM stack (counting from 0).
                             
        CopyToSliceNum   - (Integer) Slice index in the Target DICOM stack
                           where the contour will be copied to (counting from 
                           0).
                            
                           If the Frame of Reference of the Source and Target
                           DICOMs are different this input will be ignored.
                         
        CopyType         - (String) Denotes whether the contour will be copied
                           literally or exactly.  Acceptable inputs are
                           'literal' and 'exact'
                             
                           If the Frame of Reference of the Source and Target
                           DICOMs are the same this input will be ignored.
                            
        TrgSSLabel       - (String) Defines the Structure set label to assign 
                           to the Target ROI Collection.
    
        ExportRoi        - (Boolean) Denotes whether or not the Target ROI
                           ROI Collection will be exported.
                            
        LogToConsole     - (Boolean) Denotes whether or not some results will
                           be logged to the console.
                            
        
        
    Returns:
        TrgRoi        - Target ROI object
    
    """
    
    # Import packages and functions:
    import pydicom
    import copy
    import os
    import time
    from GetDicoms import GetDicoms
    from GetImageAttributes import GetImageAttributes
    
    
    # Start timing:
    times = []
    times.append(time.time())
    
    
    """ 
    **************************************************************************
    Get the indices of the slices in the Source DICOMs that have contours. 
    **************************************************************************
    """
    
    # Get all DICOMs in SrcDcmDir:
    SrcDcmFpaths, SrcDcms = GetDicoms(DirPath=SrcDcmDir, 
                                      SortMethod='slices', 
                                      Debug=False)
    
    """ Get the Source and Target Frames of reference and other info. """
    
    #if TrgDcmDir != SrcDcmDir:
    
    # Get all DICOMs in TrgDcmDir:
    TrgDcmFpaths, TrgDcms = GetDicoms(DirPath=TrgDcmDir, 
                                      SortMethod='slices', 
                                      Debug=False)
    
    """ Following section just for info. """
    # Get the Image Attributes for Source and Target:
    SrcOrigin, SrcDirs,\
    SrcSpacings, SrcDims = GetImageAttributes(DicomDir=SrcDcmDir, 
                                              Package='pydicom')
    
    TrgOrigin, TrgDirs,\
    TrgSpacings, TrgDims = GetImageAttributes(DicomDir=TrgDcmDir, 
                                              Package='pydicom')
    
    # Get the FrameOfReferenceUIDs:
    """ Will be used later.. """
    SrcFORuid = SrcDcms[0].FrameOfReferenceUID
    TrgFORuid = TrgDcms[0].FrameOfReferenceUID
    
    print(f'The Source scan has image dimensions       {SrcDims}')
    print(f'The Target scan has image dimensions       {TrgDims}\n')
    print(f'The Source scan has pixel spacings         {SrcSpacings}')
    print(f'The Target scan has pixel spacings         {TrgSpacings}\n')
    print(f'The Source scan has patient origin         {SrcOrigin}')
    print(f'The Target scan has patient origin         {TrgOrigin}\n')
    print(f'The Source scan has patient directions     [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},')
    print(f'                                            {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},')
    print(f'                                            {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]\n')
    print(f'The Target scan has patient directions     [{TrgDirs[0]}, {TrgDirs[1]}, {TrgDirs[2]},')
    print(f'                                            {TrgDirs[3]}, {TrgDirs[4]}, {TrgDirs[5]},')
    print(f'                                            {TrgDirs[6]}, {TrgDirs[7]}, {TrgDirs[8]}]\n')
    #print(f'The Source scan has patient orientation \n{SrcDirs}')
    #print(f'The Target scan has patient orientation  \n{TrgDirs}')
    print(f'The Source scan has Frame of reference UID {SrcFORuid}')
    print(f'The Target scan has Frame of reference UID {TrgFORuid}\n')

    if SrcFORuid == TrgFORuid:
        print('(The Source and Target scans have the same Frame of reference.)')
        print(f'\nThe contour on Source slice {CopyFromSliceNum} will be',
                  f'copied to Target slice {CopyToSliceNum}.')
    else:
        print('(The Source and Target scans have different Frames of reference.)')
        if CopyType == 'literal':
            print(f'\nThe contour on Source slice {CopyFromSliceNum} will be',
                  f'copied to Target slice {CopyToSliceNum} ({CopyType} copy).')
        else: # CopyType = 'exact'
            print(f'\nThe contour points on Source slice {CopyFromSliceNum}',
                  'will be projected onto the nearest slice\n(point by point)',
                  f'in the Target stack ({CopyType} copy).')
        
    """ End of section just for info. """
            
    #else:
    #    TrgDcms = copy.deepcopy(SrcDcms)
    
    if TrgDcmDir == SrcDcmDir:
        TrgDcms = copy.deepcopy(SrcDcms)
    
    
    # Get the Source ROI:
    SrcRoi = pydicom.dcmread(SrcRoiFpath)
    
    
    # Get the list of the Source DICOMs SOP Instance UIDs:
    SrcSOPuids = [] 
    
    for dicom in SrcDcms:
        SrcSOPuids.append(dicom.SOPInstanceUID)
        
    
    
    # Loop through each Contour Image Sequence in SrcRoi, and find the index  
    # of the matching SOP Instance UID from the DICOMs:
    
    # Initialise the Source Contour Image Sequences to DICOM (slice) indices:
    SrcCntrImSeqsToDcmInds = [] 
    
    # Loop through each Contour Image Sequence:
    sequences = SrcRoi.ReferencedFrameOfReferenceSequence[0]\
                      .RTReferencedStudySequence[0]\
                      .RTReferencedSeriesSequence[0]\
                      .ContourImageSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SrcSOPuids:
        SrcCntrImSeqsToDcmInds.append(SrcSOPuids.index(uid))
        
    
    if LogToConsole:
        print('\nThere are contours for Contour Image Sequences:    ', 
              SrcCntrImSeqsToDcmInds)
          
    
    # Loop through each Contour Sequence and find the index of the 
    # matching SOP Instance UID from the DICOM:
    
    # Initialise the Source Contour Sequences to DICOM (slice) indices:
    SrcCntrSeqsToDcmInds = [] 
    
    # Loop through each Contour Sequence:
    sequences = SrcRoi.ROIContourSequence[0].ContourSequence
    
    for sequence in sequences:
        # Get the Reference SOP Instance UID:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in SrcSOPuids:
        SrcCntrSeqsToDcmInds.append(SrcSOPuids.index(uid))
    
    if LogToConsole:
        print('\nThere are contours for Contour Sequence Sequences: ', 
              SrcCntrSeqsToDcmInds)
    
    
    
    """
    Note:
        SrcCntrImSeqsToDcmInds and SrcCntrSeqsToDcmInds should be identical.
    """
    # Check that they are:
    if SrcCntrImSeqsToDcmInds != SrcCntrSeqsToDcmInds:
        print(f'\nSrcCntrImSeqsToDcmInds ({SrcCntrImSeqsToDcmInds}) !=',
              f'SrcCntrSeqsToDcmInds ({SrcCntrSeqsToDcmInds})')
    
    
    
    """ 
    **************************************************************************
    Make sure that the index of the slice to be copied (CopyFromSliceNum) has  
    a contour in SrcRoi. 
    **************************************************************************
    """
    
    
    if LogToConsole:
        print('\nChecking if a contour exists for Source slice',
              f'{CopyFromSliceNum}')
    
    if not CopyFromSliceNum in SrcCntrImSeqsToDcmInds:
        print(f'\nA contour does not exist on slice {CopyFromSliceNum}.')
        
        return
                
             
    
    
    """ 
    **************************************************************************
    If the Source and Target Frame of Reference are not the same, register the
    images, transform the contour points, determine which Target image plane is
    closest to each transformed point, and find the projection of that point 
    onto the nearest plane. 
    **************************************************************************
    """
    if SrcFORuid != TrgFORuid:
        from GetInputPoints import GetInputPoints
        from ContourInterpolatingFuncs import RegisterImageStacks
        from ContourInterpolatingFuncs import TransformFixedContours_v2
        from ContourInterpolatingFuncs import ResampleAllTransformedContours
        from ContourInterpolatingFuncs import GetIndsOfContourDataWithContourType
        #from ContourInterpolatingFuncs import CreateContourSweepData_v2
        from ContourInterpolatingFuncs import GetProjectionOfContourPtsOnPlanes as GetProjContourPts 
        from ContourInterpolatingFuncs import FlattenListOfPts
        
        
        # Get contour points from the RT-STRUCT file:
        FixPtsPCS, FixPtsBySliceAndContourPCS,\
        FixPtsICS, FixPtsBySliceAndContourICS,\
        LUT, PointData, FixContourData,\
        ContourData = GetInputPoints(FixedDicomDir=SrcDcmDir, 
                                     FixedRoiFpath=SrcRoiFpath)
        
        # Register the image stacks:
        if LogToConsole:
            print('\nRegistering image stacks...')
     
        FixIm, MovIm, RegIm,\
        ElastixImFilt = RegisterImageStacks(FixedDicomDir=SrcDcmDir, 
                                            MovingDicomDir=TrgDcmDir)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        #if LogToConsole:
        print(f'\nTook {Dtime} s to register the Moving (Target)', 
              f'image stack to the Fixed (Source) image stack.')
            
            
        
        # Transform the fixed contours:
        if LogToConsole:
            print('\nTransforming points...')
     
        ContourData = TransformFixedContours_v2(ContourData, ElastixImFilt, 
                                                MovIm)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        #if LogToConsole:
        print(f'\nTook {Dtime} s to transform the Fixed contour points.')

            
        
        
        """
        Comment 23/09:  
            For now I won't resample contours since resampling is only useful
            if interested in joining points along adjacent contours ("sweep 
            lines").  So enclosed the code below in "if False".
        """
        if False:
            # Resample the fixed contours:
            if LogToConsole:
                print('\nResampling the transformed Fixed contours...')
                
            # Chose the minimum inter-node spacing for resampling:
            dP = 0.2 # mm
            
            # Resample all transformed contours:
            ContourData = ResampleAllTransformedContours(ContourData, dP)
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            #if LogToConsole:
            # Get the indices of the slices that contain contours in the Fixed image:
            Cinds = GetIndsOfContourDataWithContourType(ContourData=ContourData, 
                                                        ContourTypeNo=1)
        
            print(f'\nTook {Dtime} s to re-sample {len(Cinds)} contours.')

            
            
            
            # Create the dictionary ContourSweepData:
            """ Might will be required later. Commenting out for now. """
            #ContourSweepData = CreateContourSweepData_v2(ContourData)
        
        
        
        
        # Decide on which set of transformed points (which key in ContourData) 
        # to either copy (if CopyType = 'literal') or project onto the Target  
        # image planes (if CopyType = 'exact'), i.e. 'TraPtsPCS', 'TraSSPtsPCS', 
        # 'TraOSPtsPCS' or any combination:
        #KeysIn = ['TraPtsPCS', 'TraSSPtsPCS', 'TraOSPtsPCS'] 
        """ Comment 23/09: 
            Commented above line since skipped resampling. 
        """
        KeysIn = ['TraPtsPCS']
            
        #print('\nkeys in ContourData:', list(ContourData.keys()))
        #print(f'\nContourData[{KeysIn[0]}] = {ContourData[KeysIn[0]]}')
            
        """
        Need to deal with a literal v exact copy differently.
        """
        if CopyType == 'exact':
            # Find the projection of the transformed points onto the Moving 
            # image planes:
            if LogToConsole:
                print('\nFinding projected/intersection points...')
                
            
            for KeyIn in KeysIn:
                """
                Note: For now the copying will work with single contours only, so
                      index the first contour in 
                      ContourData[KeyIn][CopyFromSliceNum]:
                """
                ContourPts = ContourData[KeyIn][CopyFromSliceNum][0]
                
                #print(f'\nContourPts = {ContourPts}')
                
                ProjPtsBySlice = GetProjContourPts(ContourPts=ContourPts,
                                                   FixOrigin=SrcOrigin, 
                                                   FixDirs=SrcDirs, 
                                                   FixSpacings=SrcSpacings, 
                                                   FixDims=SrcDims, 
                                                   MovOrigin=TrgOrigin, 
                                                   MovDirs=TrgDirs, 
                                                   MovSpacings=TrgSpacings, 
                                                   MovDims=TrgDims,
                                                   MovIm=MovIm, 
                                                   ElastixImFilt=ElastixImFilt, 
                                                   LogToConsole=LogToConsole)
                
                # Add the projected/intersection points to ContourData for 
                # KeyOut:
                KeyOut = KeyIn.replace('Tra', 'Mov')
                
                ContourData.update({KeyOut : ProjPtsBySlice})
                
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            #if LogToConsole:
            print(f'\nTook {Dtime} s to find the projections of {len(KeysIn)}',
                  'contours.')

                
                
                
            # Decide on which set of projected/intersection points to copy into
            # the new ROI Collection (i.e. 'MovPtsPCS', 'MovSSPtsPCS' or 
            # 'MovOSPtsPCS'):
            key = 'MovPtsPCS'
            #key = 'MovSSPtsPCS'
            #key = 'MovOSPtsPCS'
            
            TrgCntrPtsBySlice = ContourData[key]
            
            #print(f'\nTrgCntrPtsBySlice = {TrgCntrPtsBySlice}')
            
        
        else: # CopyType = 'literal'
            # Use the transformed points (need to index the 0th contour of
            # ContourData['TraPtsPCS'][CopyFromSliceNum] since the points are 
            # in a nested list):
            TrgCntrPts = ContourData['TraPtsPCS'][CopyFromSliceNum][0]
             
    
        
        
    """ 
    **************************************************************************
    Create new ROI object and modify it. 
    **************************************************************************
    """
    
    # Use SrcRoi as a template for TrgRoi: 
    TrgRoi = copy.deepcopy(SrcRoi)
    
    
    
    """ Add/remove Contour Image Sequences and Contour Sequences if there are 
    fewer/more sequences in SrcRoi than the number required in TrgRoi.
    
    Note:  For the moment there will always only be one contour copied over to
           TrgRoi.
    """
    
    # The number of contour image sequences in SrcRoi (= S) and the number
    # required in TrgRoi (= T = 1 for now):
    S = len(SrcCntrImSeqsToDcmInds)
    #T = len(CopySliceNums)
    if (SrcFORuid == TrgFORuid) or (SrcFORuid != TrgFORuid and CopyType=='literal'):
        # For now a single contour is copied requiring a single contour 
        # sequence: 
        T = 1
    else:
        # The number of contour sequences required for the Target ROI is given
        # by the number of non-empty lists in TrgCntrPtsBySlice:
        from ContourInterpolatingFuncs import GetIndsOfNonEmptyItems
        
        Tinds = GetIndsOfNonEmptyItems(List=TrgCntrPtsBySlice)
        
        T = len(Tinds)
        
        print(f'\nThe points in Source slice {CopyFromSliceNum}',
              f'project onto Target slices {Tinds}.')
        
        
    
    if LogToConsole:       
        print(f'\nThere are {S} contour sequences in SrcRoi.')
        print(f'There are {T} contour sequences required for TrgRoi.')
        
    # Return if T = 0 or T > 1:
    #if (T == 0) or (T > 1):
    if T == 0:
        #print('CopySliceNums must be a list containing a single integer value.')
        print(f'\nT = {T} - there are no contour points. Exiting.')
        
        return
    
    
    """ Adjust the number of contour sequences and contour image sequences. """
    
    # Add/remove contour image sequences and contour sequences by 
    # appending/popping the last sequence N times, where N is the difference
    # between S and T:
    if S != T:
        if T > S:
            #if LogToConsole:       
            #    print(f'{T - S} more contour sequences are required in TrgRoi.')
                
            for i in range(T - S):
                TrgRoi.ReferencedFrameOfReferenceSequence[0]\
                         .RTReferencedStudySequence[0]\
                         .RTReferencedSeriesSequence[0]\
                         .ContourImageSequence\
                         .append(copy.deepcopy(TrgRoi.ReferencedFrameOfReferenceSequence[0]\
                                                     .RTReferencedStudySequence[0]\
                                                     .RTReferencedSeriesSequence[0]\
                                                     .ContourImageSequence[-1]))

                TrgRoi.ROIContourSequence[0].ContourSequence\
                      .append(copy.deepcopy(TrgRoi.ROIContourSequence[0]\
                                                  .ContourSequence[-1]))
                
        else:
            #if LogToConsole:       
            #    print(f'There are {S - T} excess contour sequences in SrcRoi.')
                
            for i in range(S - T):
                TrgRoi.ReferencedFrameOfReferenceSequence[0]\
                      .RTReferencedStudySequence[0]\
                      .RTReferencedSeriesSequence[0]\
                      .ContourImageSequence.pop()
                
                TrgRoi.ROIContourSequence[0].ContourSequence.pop()
                

    T = len(TrgRoi.ROIContourSequence[0].ContourSequence)
    
    if LogToConsole:       
        print(f'There are now {T} contour sequences in TrgRoi.')
        
        
        
    
    
    
    
    """ Modify TrgRoi. """

    # Generate a new SOP Instance UID:
    TrgRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    TrgRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    """ 
    Some tags will not need modifying if copying is within the same 
    session/scan. 
    """
    # If the copying is across different scans:
    if TrgDcmDir != SrcDcmDir:
    
        # Modify the Study Date and Time:
        TrgRoi.StudyDate = TrgDcms[0].StudyDate
        """ 
        For some data, the Series Time and Instance Creation Time in the DICOM 
        both match the Study Time in the ROI, whereas the DICOM's Study Time 
        does not.
        """
        #TrgRoi.StudyTime = TrgDcms[0].StudyTime 
        TrgRoi.StudyTime = TrgDcms[0].SeriesTime
        
        # Modify the Patient's Name and ID:
        TrgRoi.PatientName = TrgDcms[0].PatientName
        
        TrgRoi.PatientID = TrgDcms[0].PatientID
        
        """ Skipping other non-critical tags... """
        
        # Modify the Study Instance UID:
        TrgRoi.StudyInstanceUID = TrgDcms[0].StudyInstanceUID
        
        # Modify the Frame of Reference UID:
        TrgRoi.FrameOfReferenceUID = TrgDcms[0].FrameOfReferenceUID
        
        # Modify the Referenced Frame of Reference Sequence:
        TrgRoi.ReferencedFrameOfReferenceSequence[0]\
              .FrameOfReferenceUID = TrgDcms[0].FrameOfReferenceUID
        
        # Modify the Referenced SOP Instance UID in the RT Referenced Study 
        # Sequence:
        TrgRoi.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .ReferencedSOPInstanceUID = TrgDcms[0].StudyInstanceUID
        
        # Modify the Series Instance UID in the RT Referenced Series Sequence:
        TrgRoi.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .SeriesInstanceUID = TrgDcms[0].SeriesInstanceUID
        
        # Modify the Referenced Frame of Reference UID for the Structure Set  
        # ROI Sequence:
        TrgRoi.StructureSetROISequence[0]\
              .ReferencedFrameOfReferenceUID = TrgDcms[0].FrameOfReferenceUID
    
        
        # Modify the RT ROI Observations Sequence:
        """ Not sure if these need to be changed """
        #TrgRoi.RTROIObservationsSequence[0].ObservationNumber = 
        #TrgRoi.RTROIObservationsSequence[0].ReferencedROINumber = 
        #TrgRoi.RTROIObservationsSequence[0].RTROIInterpretedType = 
        #TrgRoi.RTROIObservationsSequence[0].ROIInterpreter = 
    
    
    # Modify the Structure Set Label (this appears under Name in XNAT):
    #TrgRoi.StructureSetLabel = 'Copy_of_' + SrcRoi.StructureSetLabel
    TrgRoi.StructureSetLabel = TrgSSLabel
    
    # Modify the ROI Name for the Structure Set ROI Sequence:             
    TrgRoi.StructureSetROISequence[0]\
          .ROIName = TrgRoi.StructureSetLabel
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    TrgRoi.StructureSetDate = NewDate
    TrgRoi.StructureSetTime = NewTime
    
    
    
    # Cycle through each sequence in TrgRoi:
    """ For now T = 1. """ 
    """ 23/09 Above comment no longer valid since T can be > 1. Need to verify
    that code below can deal with that. """
    for i in range(T):
        if T == 1:
            # Either TrgDcmDir = SrcDcmDir or TrgDcmDir != SrcDcmDir and 
            # CopyType = 'literal'.
            s = CopyToSliceNum
            
        else:
            # TrgDcmDir != SrcDcmDir and CopyType = 'exact'.
            s = Tinds[i]
    
        # Modify the Referenced SOP Class UID to match the SOP Class UID of the
        # Target DICOM the ROI is to be copied to:
        TrgRoi.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPClassUID = TrgDcms[s].SOPClassUID
        
        # Modify the SOP Instance UID to match the SOP Instance UID of the
        # Target DICOM the ROI is to be copied to:
        TrgRoi.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0]\
              .RTReferencedSeriesSequence[0]\
              .ContourImageSequence[i]\
              .ReferencedSOPInstanceUID = TrgDcms[s].SOPInstanceUID
        
        # Modify the Referenced SOP Class UID to match the SOP Class UID of the
        # Target DICOM the ROI is to be copied to:
        TrgRoi.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPClassUID = TrgDcms[s].SOPClassUID
        
        # Modify the SOP Instance UID to match the SOP Instance UID of the 
        # Target DICOM the ROI is to be copied to:
        TrgRoi.ROIContourSequence[0]\
              .ContourSequence[i]\
              .ContourImageSequence[0]\
              .ReferencedSOPInstanceUID = TrgDcms[s].SOPInstanceUID
        
        
        #print(f'\nTrgRoi.ROIContourSequence[0].ContourSequence[i].NumberOfContourPoints = {TrgRoi.ROIContourSequence[0].ContourSequence[i].NumberOfContourPoints}')
        
        #print(f's = {s}')
        
        #print(f'\nSrcRoi.ROIContourSequence[0].ContourSequence[s].NumberOfContourPoints = {SrcRoi.ROIContourSequence[0].ContourSequence[s].NumberOfContourPoints}')
        
        
        """ 
        Proceed differently if the Source and Target Frame of Reference are the 
        same or not. 
        """
        if SrcFORuid == TrgFORuid:
        
            # Get the sequence number that corresponds to the slice number to
            # be copied:
            CopyFromSequenceNum = SrcCntrImSeqsToDcmInds.index(CopyFromSliceNum)
            
            # Modify the Number of Contour Points to match the value of the 
            # Source contour being copied:
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = SrcRoi.ROIContourSequence[0]\
                                                 .ContourSequence[CopyFromSequenceNum]\
                                                 .NumberOfContourPoints
            
            # Modify the Contour Number to match the value of the Source 
            # contour being copied:
            #TrgRoi.ROIContourSequence[0]\
            #      .ContourSequence[i]\
            #      .ContourNumber = SrcRoi.ROIContourSequence[0]\
            #                             .ContourSequence[CopyFromSequenceNum]\
            #                             .ContourNumber
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourNumber = f"{i + 1}"
            
            # Modify the Contour Data to match the value of the Source contour
            # being copied:
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourData = SrcRoi.ROIContourSequence[0]\
                                       .ContourSequence[CopyFromSequenceNum]\
                                       .ContourData
                                       
        else: # SrcFORuid != TrgFORuid 
            """ 
            If CopyType = 'literal', there is a single list of contour points:
            TrgCntrPts.
            
            If CopyType = 'exact', there is a list of contours: 
            TrgCntrPtsBySlice, which may contain any number of empty lists ([])
            and non-empty lists of contours (of length 1) with a sub-list of 
            points.
            """
            if T > 1:
                # Update TrgCntrPts:
                """ Note:
                    Because each list in TrgCntrPtsBySlice represents a list of
                    contours, the 0th contour must be indexed (there are at 
                    most 1 contour for now).
                """
                TrgCntrPts = TrgCntrPtsBySlice[s][0]
            
                
            # Modify the Number of Contour Points:
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .NumberOfContourPoints = f"{len(TrgCntrPts)}"
            
            # Modify the Contour Number:
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourNumber = f"{i + 1}"
            
            # Flatten the list of contour points to a list of coordinates
            # (required format of ContourData):
            TrgContourData = FlattenListOfPts(ListOfPts=TrgCntrPts)
            
            # Modify the Contour Data:
            TrgRoi.ROIContourSequence[0]\
                  .ContourSequence[i]\
                  .ContourData = TrgContourData
    
        
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    #if LogToConsole:
    print(f'\nTook {Dtime} s to create and modify a new ROI.')
        
           
    # Get the filename of the original RT-Struct file:
    SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    # Create a new filename (this will appear under Label in XNAT):
    TrgRoiFname = 'AIM_' + NewDate + '_' + NewTime + '_copied_from_' \
                  + SrcRoiFname
    
    TrgRoiFpath = os.path.join(TrgRoiDir, TrgRoiFname)
    
    # Save the new DICOM-RTSTRUCT file:
    if ExportRoi:
        TrgRoi.save_as(TrgRoiFpath)
        
        print('\nROI exported to:\n', TrgRoiFpath)
        
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    #if LogToConsole:
    print(f'\nDone.  Took {Dtime} s in total.')
    
    return TrgRoi




"""
    
    ----------------------------------------------------------------------------------------------
    Replace in TrgRoi:                                      The value of the tag in the DICOM:
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