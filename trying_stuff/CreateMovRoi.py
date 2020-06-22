# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 08:48:55 2020

@author: ctorti
"""




def CreateMovRoi(FixRoiFpath, MovDicomDir, MovPtsArr_ICS, Debug):
    # Import packages:
    import pydicom
    import copy
    import time
    import os
    from GetDicoms import GetDicoms
    
    """ 
    This will need to be tested for cases where there's more than one contour 
    in an image! 
    """
    
    
    " Create ContourData in the required format. """

    # Get the indices of the non-empty arrays in OutPtsArr_ICS:
    MovInds = []
    
    for i in range(len(MovPtsArr_ICS)):
        if MovPtsArr_ICS[i]:
            MovInds.append(i)
            
    MovN = len(MovInds)
     
    if Debug:       
        print(f'MovInds = {MovInds}')
        print(f'MovN = {MovN}')
    
    # Create a new array of contour points in the format required for the ROI object:
    ContourData = []
    Npts = []
    
    for i in range(len(MovPtsArr_ICS)):
        if MovPtsArr_ICS[i]:
            points = MovPtsArr_ICS[i]
    
            Npts.append(len(points))
    
            #print(f'Number of points = {len(points)}')
    
            PtsFlat = []
    
            for point in points:
                # Convert coordinates of point to strings:
                PtsFlat.extend([str(point[p]) for p in range(3)])
    
            ContourData.append(PtsFlat)
            
        else:
            ContourData.append([])
            
            Npts.append(0)
            
    ContourData
    #print(f'Number of arrays in ContourData = {len(ContourData)}')
    
    
    """ Get the DICOMs. """
    
    MovDicomFpaths, MovDicoms = GetDicoms(DirPath=MovDicomDir, 
                                          SortMethod='slices', 
                                          Debug=False)
    
    """ Create the ROI object. """
    
    # Use the fixed ROI as a template for the moving image:
    FixRoi = pydicom.dcmread(FixRoiFpath)
    
    MovRoi = copy.deepcopy(FixRoi)
    
    
    """ Add additional (or remove) Contour Image Sequences and Contour Sequences 
    if there are more (or fewer) arrays in ContourData than the number of 
    Contour Sequences in MovRoi.
    """
    
    # Get the number of contour image sequences in the copied ROI:
    FixN = len(MovRoi.ROIContourSequence[0].ContourSequence)
    
    if Debug:       
        print(f'Number of contour sequences in the original ROI  = {FixN}')
        print(f'Number of contour sequences required for new ROI = {MovN}')
    
    # Add/remove contour image sequences and contour sequences by 
    # appending/popping the last sequence N times, where N is the difference
    # between MovN and FixN:
    if FixN != MovN:
        if MovN > FixN:
            if Debug:       
                print(f'There are {MovN - FixN} more contour sequences required in the ROI. ')
                
            for i in range(MovN - FixN):
                MovRoi.ReferencedFrameOfReferenceSequence[0]\
                .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                .ContourImageSequence\
                .append(MovRoi.ReferencedFrameOfReferenceSequence[0]\
                .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                .ContourImageSequence[-1])

                MovRoi.ROIContourSequence[0].ContourSequence\
                .append(MovRoi.ROIContourSequence[0].ContourSequence[-1])
                
        else:
            if Debug:       
                print(f'There are {FixN - MovN} excess contour sequences in the ROI. ')
                
            for i in range(FixN - MovN):
                MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].\
                RTReferencedSeriesSequence[0].ContourImageSequence.pop()
                
                MovRoi.ROIContourSequence[0].ContourSequence.pop()
                
    #MovN = len(MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].\
    #RTReferencedSeriesSequence[0].ContourImageSequence)
    MovN = len(MovRoi.ROIContourSequence[0].ContourSequence)
    
    if Debug:       
        print(f'Number of contour sequences in new ROI           = {MovN}')
    
    
    
    """ Modify MovRoi. """

    # Generate a new SOP Instance UID:
    MovRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Modify the Study Date and Time:
    MovRoi.StudyDate = MovDicoms[0].StudyDate
    """ 
    For current data, the Series Time and Instance Creation Time in the DICOM 
    both match the Study Time in the ROI, whereas the DICOM's Study Time does not.
    """
    #MovRoi.StudyTime = MovDicoms[0].StudyTime 
    MovRoi.StudyTime = MovDicoms[0].SeriesTime
    
    # Modify the Patient's Name and ID:
    MovRoi.PatientName = MovDicoms[0].PatientName
    MovRoi.PatientID = MovDicoms[0].PatientID
    
    """ Skipping other non-critical tags... """
    
    # Modify the Study Instance UID:
    MovRoi.StudyInstanceUID = MovDicoms[0].StudyInstanceUID
    
    # Generate a new Series Instance UID:
    MovRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    # Modify the Frame of Reference UID:
    MovRoi.FrameOfReferenceUID = MovDicoms[0].FrameOfReferenceUID
    
    # Modify the Structure Set Label:
    MovRoi.StructureSetLabel = 'Copy_of_' + MovRoi.StructureSetLabel
    
    # Modify the Structure Set Date and Time to the present:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    MovRoi.StructureSetDate = NewDate
    MovRoi.StructureSetTime = NewTime
    
    # Modify the Referenced Frame of Reference Sequence:
    MovRoi.ReferencedFrameOfReferenceSequence[0]\
    .FrameOfReferenceUID = MovDicoms[0].FrameOfReferenceUID
    
    # Modify the Referenced SOP Instance UID in the RT Referenced Study Sequence:
    MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
    .ReferencedSOPInstanceUID = MovDicoms[0].StudyInstanceUID
    
    # Modify the Series Instance UID in the RT Referenced Series Sequence:
    MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
    .RTReferencedSeriesSequence[0].SeriesInstanceUID = MovDicoms[0].SeriesInstanceUID
    
    # Modify the Referenced Frame of Reference UID and ROI Name for the 
    # Structure Set ROI Sequence:
    MovRoi.StructureSetROISequence[0].ReferencedFrameOfReferenceUID = MovDicoms[0].FrameOfReferenceUID
    MovRoi.StructureSetROISequence[0].ROIName = MovRoi.StructureSetLabel
    
    # Modify the RT ROI Observations Sequence:
    """ Not sure if these need to be changed """
    #MovRoi.RTROIObservationsSequence[0].ObservationNumber = 
    #MovRoi.RTROIObservationsSequence[0].ReferencedROINumber = 
    #MovRoi.RTROIObservationsSequence[0].RTROIInterpretedType = 
    #MovRoi.RTROIObservationsSequence[0].ROIInterpreter = 
    
    
    # Modify the Contour Image Sequences in the Referenced Frame of Reference Sequence: 
    
    # Cycle through each Contour Image Sequence and modify the Referenced SOP Class UID
    # and Referenced SOP Instance UID with those of MovDicoms[0]:
    for i in range(len(MovInds)):
        if Debug:
            print(f'i = {i}, MovInds[{i}] = {MovInds[i]}')
        
        # Modify the Referenced SOP Class UID:
        MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .RTReferencedSeriesSequence[0].ContourImageSequence[i]\
        .ReferencedSOPClassUID = MovDicoms[MovInds[i]].SOPClassUID
        
        # Modify the SOP Instance UID:
        MovRoi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .RTReferencedSeriesSequence[0].ContourImageSequence[i]\
        .ReferencedSOPInstanceUID = MovDicoms[MovInds[i]].SOPInstanceUID
        
    
    # Modify the Contour Sequences in the ROI Contour Sequence: 
    
    # Cycle through each Contour Sequence and modify the Referenced SOP Class UID and
    # Referenced SOP Instance UID in the Contour Image Sequence with those of MovDicoms[0],
    # and modify the Contour Geometric Type, Number of Contour Points, Contour Number and
    # Contour Data:
    for i in range(len(MovInds)):
        if Debug:       
            print(f'i = {i}, MovInds[{i}] = {MovInds[i]}')
        
        # Modify the Referenced SOP Class UID:
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .ContourImageSequence[0].ReferencedSOPClassUID = MovDicoms[MovInds[i]].SOPClassUID
        
        # Modify the SOP Instance UID:
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .ContourImageSequence[0].ReferencedSOPInstanceUID = MovDicoms[MovInds[i]].SOPInstanceUID
        
        # Modify the Contour Geometric Type:
        """ Since the contours are not closed change from 'CLOSED_PLANAR' to 'OPEN_PLANAR' """
        #print('\nBefore:', MovRoi.ROIContourSequence[0].ContourSequence[s].ContourImageSequence[0].ReferencedSOPClassUID)
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .ContourGeometricType = 'OPEN_PLANAR'
        #print('\nAfter:', MovRoi.ROIContourSequence[0].ContourSequence[s].ContourGeometricType)
        
        #print(MovRoi.ROIContourSequence[0].ContourSequence[s].ContourGeometricType)
        
        # Modify the Number of Contour Points:
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .NumberOfContourPoints = f"{Npts[MovInds[i]]}"  
        
        
        # Modify the Contour Number:
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .ContourNumber = f"{i+1}" 
        
        # Modify the Contour Data:
        MovRoi.ROIContourSequence[0].ContourSequence[i]\
        .ContourData = ContourData[MovInds[i]]
        
           
    # Get the filename of the original RT-Struct file:
    FixRoiFname = os.path.split(FixRoiFpath)[1]
    
    # Create a new filename:
    MovRoiFname = 'AIM_' + NewDate + '_' + NewTime + '_transformed_from_' + FixRoiFname
    
    # Get the current working directory:
    CWD = os.getcwd()
    
    MovRoiFpath = os.path.join(CWD, MovRoiFname)
    
    # Save the new DICOM-RTSTRUCT file:
    #MovRoi.save_as(MovRoiFpath)
    
    #print('\nROI for Moving image exported to:\n', MovRoiFpath)
    
    return MovRoi
    
    



    
    
    
    
    
    
    