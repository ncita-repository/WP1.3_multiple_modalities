# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:07:41 2020

@author: ctorti
"""




""" Function:
    CheckRoiToDicomsCompatibility()
    
Purpose:
    Check compatibility of an ROI Collection for a DICOM series

Input:
    DataDict - dictionary containing the DICOM and RTSTRUCT filepaths
    
    CheckKey - key that determine which values in DataDict to check

    
    
Returns:
    N/A
"""




def CheckRoiToDicomsCompatibility(DataDict, CheckKey):
    
    # Import packages:
    import pydicom
    from GetDicoms import GetDicoms
    
    # Get the DICOMs:
    DicomFpaths, Dicoms = GetDicoms(DataDict[CheckKey]['DicomDir'], \
                                            SortMethod='slices')
    
    # Get the ROI Collection:
    Roi = pydicom.dcmread(DataDict[CheckKey]['RoiFpath'])
    
    
    # Count the number of tags that do not match:
    errors = 0
    
    """ 
    Check Study Instance UIDs, Series Instance UIDs and Frame of Reference UIDs
    (assuming that all Dicoms have the same values): 
    """
    
    # Check Study Instance UIDs:
    if Dicoms[0].StudyInstanceUID != Roi.StudyInstanceUID:
        print('\nThe Study Instance UIDs are not the same!')
        print('DICOM:', Dicoms[0].StudyInstanceUID)
        print('ROI  :', Roi.StudyInstanceUID)
        
        # Increment errors:
        errors = errors + 1
        
    
    # Get the Roi's ReferencedSOPInstanceUID:
    RefSopUid = Roi.ReferencedFrameOfReferenceSequence[0]\
                   .RTReferencedStudySequence[0].ReferencedSOPInstanceUID
                 
    # Check the DICOM's Study InstanceUID Instance UID against the ROI's 
    # Referenced SOP Instance UID:
    if Dicoms[0].StudyInstanceUID != RefSopUid:
        print('\nThe ROI Referenced SOP Instance UID and DICOM Study', \
              'Instance UID are not the same!')
        print('DICOM:', Dicoms[0].SeriesInstanceUID)
        print('ROI  :', RefSopUid)
        
        # Increment errors:
        errors = errors + 1
        
    
    # Get the Roi's RT Referenced Series Sequence Series Instance UID:
    RTRefSeriesSeqSeriesUid = Roi.ReferencedFrameOfReferenceSequence[0]\
                                 .RTReferencedStudySequence[0]\
                                 .RTReferencedSeriesSequence[0]\
                                 .SeriesInstanceUID
                      
    # Check the Series Instance UIDs:
    if Dicoms[0].SeriesInstanceUID != RTRefSeriesSeqSeriesUid:
        print('\nThe ROI RT Referenced Series Sequence Series Instance UID', \
              'and DICOM Series Instance UID are not the same!')
        print('DICOM:', Dicoms[0].SeriesInstanceUID)
        print('ROI  :', RTRefSeriesSeqSeriesUid)
        
        # Increment errors:
        errors = errors + 1
        
        
    # Check the Frame of Reference UIDs:
    if Dicoms[0].FrameOfReferenceUID != Roi.FrameOfReferenceUID:
        print('\nThe Frame Of Reference UIDs are not the same!')
        print('DICOM:', Dicoms[0].FrameOfReferenceUID)
        print('ROI  :', Roi.FrameOfReferenceUID)
        
        # Increment errors:
        errors = errors + 1
        
        
    # Get the ROI's Referenced Frame of Reference UID:
    RefFORUid = Roi.StructureSetROISequence[0].ReferencedFrameOfReferenceUID
    
    # Check the DICOM's Frame of Reference UIDs against the ROI's Referenced
    # Frame of Reference UID:
    if Dicoms[0].FrameOfReferenceUID != RefFORUid:
        print('\nThe ROI Referenced Frame Of Reference UID and DICOM', \
              'Frame Of Reference UID are not the same!')
        print('DICOM:', Dicoms[0].FrameOfReferenceUID)
        print('ROI  :', RefFORUid)
        
        # Increment errors:
        errors = errors + 1
        
        
    """ 
    Check Contour Image Sequence Referenced SOP Instance UIDs against Dicom's 
    SOP Instance UIDs:
    """
    
    # Get the Dicom's SOP Instance UIDs:
    SopUids = [Dicom.SOPInstanceUID for Dicom in Dicoms]
    
    # Get the Contour Image Sequences:
    CISs = Roi.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
              .ContourImageSequence
        
    # Get the Referenced SOP Instance UIDs from the CISs:
    CISRefSopUids = [CIS.ReferencedSOPInstanceUID for CIS in CISs]
    
    
    if errors > 0:
        print('')
        
        
    # Iterate through each CISRefSopUids:
    for i in range(len(CISRefSopUids)):
        # Check if CISRefSopUid is not in SopUids:
        if not CISRefSopUids[i] in SopUids:
            print('The ROI Referenced SOP UID in Contour Image', \
                  f'Sequence {i} does not match to any SOP UID in the DICOMs!')
            
            # Increment errors:
            errors = errors + 1
            
    
    """ 
    Check Contour Sequence Referenced SOP Instance UIDs against Dicom's 
    SOP Instance UIDs:
    """
    
    # Get the Contour Sequences:
    CSs = Roi.ROIContourSequence[0].ContourSequence
        
    # Get the Referenced SOP Instance UIDs from the CSs:
    CSRefSopUids = [CS.ContourImageSequence[0].ReferencedSOPInstanceUID for CS in CSs]
    
    
    if errors > 0:
        print('')
        
        
    # Iterate through each CSRefSopUids:
    for i in range(len(CSRefSopUids)):
        # Check if CSRefSopUid is not in SopUids:
        if not CSRefSopUids[i] in SopUids:
            print(f'The ROI Referenced SOP UID in Contour Sequence {i}', \
                  'does not match to any SOP UID in the DICOMs!')
            
            # Increment errors:
            errors = errors + 1
            
    
    # Final report:
    if errors > 0:
        print(f'\n--> There were {errors} tag mismatches!!')
    else:
        print('--> OK.')
            
    