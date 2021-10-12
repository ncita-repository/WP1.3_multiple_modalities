# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:48:48 2020

@author: ctorti
"""


"""
Function:
    ModifyRoi()
    
Purpose:
    Starting with an ROI object roi1 (from an DICOM-RTSTRUCT file corresponding
    to the 1st series of DICOMs) and dicomDict1_to_2 (the mapping from the 1st
    to the 2nd DICOM scan), copy roi1 and modify it accordingly so that the 
    relevant metadata corresponds to the DICOMs in the 2nd scan, and so, can
    be overlaid on the 2nd DICOM scan.
    
Input:
    roiDict1        - Dictionary of ROI data for the 1st DICOM scan
    
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
                      
    xnatDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths 
                    (see DownloadFilesForRoiCopy() for details)
                    
    regMethod      - Type of registration to perform on the DICOM scan
                     (collection) that the ROI Collection is to be copied to.  
                     Acceptable values include:
                        -- 'rigid' = perform rigid registration
                        -- 'affine' = perform affine registration
                        -- 'non-rigid' = perform non-rigid registration
                        -- any other character string (including '') = do not
                            perform any registration
                     This is needed because the ROI tags will either need to 
                     reflect the original tags of the 2nd DICOM scan or the
                     new tags of the new DICOMs of the registered scan.
    
    
Returns:
    roi2           - ROI object that contains the contours in roi1 but 
                    modified for overlaying over the 2nd DICOM scan;
                   - A new DICOM-RTSTRUCT file will be saved to the location
                   specified
                   
    xnatDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths
                    (see DownloadFilesForRoiCopy() for details)
                   - The dictionary is exported as a JSON file
    
"""

def ModifyRoi(xnatDict, dicomDict1_to_2, regMethod):
    # Import packages:
    import copy
    import time
    import json
    import pydicom
    import os
    
    # Load the DICOM-RTSTRUCT file whose ROIs correspond to the 1st DICOM scan:
    #roi1 = pydicom.dcmread(roiDict1['Roi fpath'])
    roi1 = pydicom.dcmread(xnatDict['fromRoiFpath'])
    
    # Use roi1 as a template for roi2:
    roi2 = copy.deepcopy(roi1)
    
    # Get the keys in roiDict1 and dicomDict1_to_2:
    #Rkeys = list(roiDict1.keys())
    D1to2keys = list(dicomDict1_to_2.keys())
    
    # First make changes to roi2 that only need to be made once (don't need
    # to be made for each contour).
    
    # Read in the first DICOM from the 2nd scan:
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        fpath2 = dicomDict1_to_2[D1to2keys[0]]['Registered closest Dicom fpath']
    else:
        fpath2 = dicomDict1_to_2[D1to2keys[0]]['Closest Dicom fpath']
    
    # Load the DICOM:
    dicom2 = pydicom.read_file(fpath2)
    
    # Create a dictionary to store the changes made to the ROI that will 
    # later be added to xnatDict:
    changesDict = {}
            
            
    # Make changes to the metadata in roi2 to reflect the corresponding 
    # metadata contained in the DICOM. First essential changes - some of which
    # are in the firt tier of attributes, but some are duplicated in nested 
    # attributes:
    
    # 1) Change the Study Instance UID in the first tier:                  
    if roi2.StudyInstanceUID==dicom2.StudyInstanceUID:
        print('\nThe Study Instance UID are the same:\n', roi2.StudyInstanceUID)
    else:
        print('\nChanging the Referenced SOP Instance UID from:\n', \
              roi2.StudyInstanceUID, '\nto:\n', dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyInstanceUID':roi2.StudyInstanceUID})
        
        # Change the StudyInstanceUID:
        roi2.StudyInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'toStudyInstanceUID':roi2.StudyInstanceUID})
    
     
    # 2) Change the Study Instance UID in the Referenced Frame Of Reference
    # Sequence:
    roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0].ReferencedSOPInstanceUID     
              
    if roiSIuid==dicom2.StudyInstanceUID:
        print('\nThe Study Instance UID are the same:\n', roiSIuid)
    else:
        print('\nChanging the Referenced SOP Instance UID from:\n', \
              roiSIuid, '\nto:\n', dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.RTRSS.RefSOPInstanceUID':roiSIuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .ReferencedSOPInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
        
        roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0].ReferencedSOPInstanceUID   
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.RTRSS.RefSOPInstanceUID':roiSIuid})
    
    
    
    # 3) Change the Series Instance UID in the Referenced Frame Of Reference
    # Sequence:
    roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
               .SeriesInstanceUID
               
    if roiSIuid==dicom2.SeriesInstanceUID:
        print('\nThe Series Instance UID are the same:\n', roiSIuid)
    else:
        print('\nChanging the Series Instance UID from:\n', \
              roiSIuid, '\nto:\n', dicom2.SeriesInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.RTRSS.RTRSS.SeriesInstanceUID':\
                            roiSIuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .RTReferencedSeriesSequence[0].SeriesInstanceUID \
        = copy.deepcopy(dicom2.SeriesInstanceUID)
        
        roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
               .SeriesInstanceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.RTRSS.RTRSS.SeriesInstanceUID':roiSIuid})
        
        
    # 4) Change the Frame of Reference UID:
    if roi2.FrameOfReferenceUID==dicom2.FrameOfReferenceUID:
        print('\nThe Frame Of ReferenceUID UID are the same:\n', \
              roi2.FrameOfReferenceUID)
    else:
        print('\nChanging the Frame Of ReferenceUID UID from:\n', \
              roi2.FrameOfReferenceUID, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromFrameOfReferenceUID':roi2.FrameOfReferenceUID})
        
        roi2.FrameOfReferenceUID = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'toFrameOfReferenceUID':roi2.FrameOfReferenceUID})
    
    
    # 5) Change the Frame of Reference UID in the Referenced Frame Of Reference
    # Sequence:
    roiFRuid = roi2.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
    
    if roiFRuid==dicom2.FrameOfReferenceUID:
        print('\nThe Frame Of ReferenceUID UID are the same:\n', roiFRuid)
    else:
        print('\nChanging the Frame Of ReferenceUID UID from:\n', \
              roiFRuid, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.FrameOfReferenceUID':roiFRuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID \
        = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        roiFRuid = roi2.ReferencedFrameOfReferenceSequence[0]\
        .FrameOfReferenceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.FrameOfReferenceUID':roiFRuid})
    
    
    # 6) Change the Frame of Reference UID in the Structure Set ROI Sequence:
    roiFRuid = roi2.StructureSetROISequence[0].ReferencedFrameOfReferenceUID
    
    if roiFRuid==dicom2.FrameOfReferenceUID:
        print('\nThe Referenced Frame Of ReferenceUID UID is the same as', \
              'the Frame Of Reference UID:\n', roiFRuid)
    else:
        print('\nChanging the Referenced Frame Of ReferenceUID UID from:\n', \
              roiFRuid, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromSSROIS.RefFrameOfReferenceUID':roiFRuid})
        
        roi2.StructureSetROISequence[0].ReferencedFrameOfReferenceUID \
        = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        roiFRuid = roi2.StructureSetROISequence[0]\
        .ReferencedFrameOfReferenceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toSSROIS.RefFrameOfReferenceUID':roiFRuid})
        
        
        
    """
    Now optional attributes that are contained in the first tier and which have
    the same tag labels in both the DICOM and DICOM-RTSTRUCT objects. Hence, 
    the attributes can be efficiently updated in a loop through a list of 
    common tag names.  
    """
    
    # The list of metadata to be altered with common tag names are:
    changesTags = ['StudyDate', 'StudyTime', 'PatientName', 'PatientID',\
                   'PatientSex']
    
    """ 
    Note:  Since some of these attributes are optional, not all DICOMs will 
    have the tags (e.g. Patient's Birth Date), so use try/except.
    """
    
    # Loop for each tag in changesTags:
    for tag in changesTags:
        try:
            if roi2[tag].value==dicom2[tag].value:
                print('\n' + tag + ' are the same:\n', roi2[tag].value)
            else:
                print('\nChanging '+ tag + ' from:\n', roi2[tag].value, \
                      '\nto:\n', dicom2[tag].value)
        
            # Update the dictionary changesDict:
            changesDict.update({ 'from' + tag : roi2[tag].value})
            
            roi2[tag].value = copy.deepcopy(dicom2[tag].value)
            
            # Update the dictionary changesDict:
            changesDict.update({ 'to' + tag : roi2[tag].value})
            
        except AttributeError:
            print('\n' + tag + ' does not exist in the DICOMs')
    


    """
    The following section of code is commented out since Study Time is dealt
    with in changesTags, even though it doesn't deal with decimals as was 
    previously done.  If this turns out to be a problem, Study Time can be 
    dealt with outside of the above for loop.
    """    
    # - Change the Study Time:
    """ Note: Some DICOM Study Times have a decimal with trailing zeros,
        e.g. '164104.000000', so convert to float, then integer, then back 
        to string. 
    """
    
    """
    newTime = str(int(float(copy.deepcopy(dicom2.StudyTime))))

    if roi2.StudyTime==newTime:
        print('\nThe Study Time are the same:\n', roi2.StudyTime)
    else:
        print('\nChanging the Study Time from:\n', roi2.StudyTime, '\nto:\n', \
              newTime)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyTime':roi2.StudyTime})
        
        roi2.StudyTime = newTime
        
        # Update the dictionary changesDict:
        changesDict.update({'toStudyTime':roi2.StudyTime})
    """
        
        
    """
    Now make additional changes that can't be efficiently put in a loop.
    """
    
    # 1) Modify the Structure Set Label by adding 
    # "_copy_from_Scan_A_to_B_YYYYMMDD_HHMMSS" to the Structure Set Label of 
    # the original RTSTRUCT, where A and B are the Scan Labels from the two 
    # scans (collections):
    newDate = time.strftime("%Y%m%d", time.gmtime())
    newTime = time.strftime("%H%M%S", time.gmtime())
    
    #newLabel = roi1.StructureSetLabel + '_copy_' + newDate \
    #+ '_' + newTime
    
    newLabel = roi1.StructureSetLabel + '_copy_from_Scan_' \
    + xnatDict['fromDicomScanLab'] + '_to_' + xnatDict['toDicomScanLab'] \
    + '_' + newDate + '_' + newTime
    
    print('\nChanging the Structure Set Label from:\n', \
          roi2.StructureSetLabel, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromStructureSetLabel':roi2.StructureSetLabel})
    
    roi2.StructureSetLabel = newLabel  
    
    # Update the dictionary changesDict:
    changesDict.update({'toStructureSetLabel':roi2.StructureSetLabel})
    
    
    # 2) Modify the ROI Name in the Structure Set ROI Sequence:
    """
    Assume that the characters before the first "_" is always the Patient Name.
    """
    roiValue = roi2.StructureSetROISequence[0].ROIName
    
    # Find the first "_":
    ind = roiValue.index('_')
    
    #"""
    #Modify the ROI Name using the Patient Name from the DICOM followed by the
    #characters including and following the first "_" in roiValue.
    #"""
    #newRoiName = str(copy.deepcopy(dicom2.PatientName)) + roiValue[ind::]
    
    """
    Modify the ROI Name by adding "Copy_of_" to the ROI Name of the original 
    RTSTRUCT and the current date and time.
    """
    #newDate = time.strftime("%Y%m%d", time.gmtime())
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    #newRoiName = 'Copy_of_' + roi1.StructureSetROISequence[0].ROIName \
    #+ newDate + '_' + newTime
    
    """
    Modify the ROI Name by making it the same as the new Structure Set Label.
    """
    
    print('\nChanging the ROI Name from:\n', \
          roi2.StructureSetROISequence[0].ROIName, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromSSROIS.ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
    
    roi2.StructureSetROISequence[0].ROIName = newLabel
    
    # Update the dictionary changesDict:
    changesDict.update({'toSSROIS.ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
                    
    
    
    """
    Now make essential changes to the Referenced Frame of Reference UID that 
    need to be made for each Contour Image Sequence:
    """
    
    # Get the Contour Image Sequences:
    sequences = roi2.ReferencedFrameOfReferenceSequence[0]\
        .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
        .ContourImageSequence

    #print('\nLooping through each Contour Image Sequence to:'\
    #      '\n- get the Referenced SOP Instance UID belonging to the 1st DICOM '\
    #      'scan (collection)'\
    #      '\n- find the index of the matching SOP Instance UID in '\
    #      'dicomDict1_to_2'\
    #      '\n- use the index to get the filepath of the DICOM slice in the '\
    #      '2nd scan (collection) that is closest in position to that of the '\
    #      '1st DICOM scan (collection)'\
    #      '\n- modify the (ROI\'s) Referenced SOP Instance UID to the '\
    #      'corresponding SOP Instance UID from the closest DICOM in the 2nd '\
    #      'scan (collection)')
    
    print('\nLooping through each Contour Image Sequence to modify the '\
          'Referenced SOP Instance UIDs..')
    
    
    # Store the original and modified Referenced SOP Instance UIDs from the
    # Referenced Frame Of Reference Sequence:
    fromRFORSrefsopuids = []
    toRFORSrefsopuids = []
    
    # Cycle through each Contour Image Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ReferencedSOPInstanceUID
        
        # Store this fromRFORSrefsopuid:
        fromRFORSrefsopuids.append(refsopuid)
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # The filepath of the slice from the 2nd DICOM scan that is closest in 
        # postion to the slice in the 1st DICOM scan with this SOP UID is:
        if regMethod in ['rigid', 'affine', 'non-rigid']:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
        else:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
        
        # Load the DICOM from the 2nd scan:
        dicom2 = pydicom.read_file(fpath2)
        
        # Store this toRFORSrefsopuid:
        toRFORSrefsopuids.append(dicom2.SOPInstanceUID)
    
        if refsopuid==dicom2.SOPInstanceUID:
            print('\n- The DICOM\'s SOP Instance UID is the same as the ', \
                  'ROI\'s Referenced SOP Instance UID:\n', refsopuid)
        else:
            print('\n- Changing the Referenced SOP Instance UID from:\n', \
                  refsopuid, '\nto:\n', dicom2.SOPInstanceUID)
            
            roi2.ReferencedFrameOfReferenceSequence[0]\
            .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
            .ContourImageSequence[i].ReferencedSOPInstanceUID \
            = copy.deepcopy(dicom2.SOPInstanceUID)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromRFORS.RTRSS.RTRSS.CIS.RefSOPInstanceUIDs':\
                      fromRFORSrefsopuids})
    changesDict.update({'toRFORS.RTRSS.RTRSS.CIS.RefSOPInstanceUIDs':\
                      toRFORSrefsopuids})
                    
      
              
    """
    Now make essential changes to the Referenced SOP Instance UID that need to 
    be made for each Contour Sequence:
    """
    
    # Get the Contour Sequences:
    sequences = roi2.ROIContourSequence[0].ContourSequence

    #print('\nLooping through each Contour Sequence to:'\
    #      '\n- get the Referenced SOP Instance UID belonging to the 1st DICOM '\
    #      'scan (collection)'\
    #      '\n- find the index of the matching SOP Instance UID in '\
    #      'dicomDict1_to_2'\
    #      '\n- use the index to get the filepath of the DICOM slice in the '\
    #      '2nd scan (collection) that is closest in position to that of the '\
    #      '1st DICOM scan (collection)'\
    #      '\n- modify the (ROI\'s) Referenced SOP Instance UID to the '\
    #      'corresponding SOP Instance UID from the closest DICOM in the 2nd '\
    #      'scan (collection)')
    
    print('\nLooping through each Contour Sequence to modify the '\
          'Referenced SOP Instance UIDs..')
       
    # Store the original and modified Referenced SOP Instance UIDs from the
    # ROI Contour Sequence:
    fromROICSrefsopuids = []
    toROICSrefsopuids = []
    
    # Cycle through each Contour Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ContourImageSequence[0]\
        .ReferencedSOPInstanceUID
        
        # Store this fromROICSrefsopuid:
        fromROICSrefsopuids.append(refsopuid)
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # The filepath of the slice from the 2nd DICOM scan that is closest in 
        # postion to the slice in the 1st DICOM scan with this SOP UID is:
        if regMethod in ['rigid', 'affine', 'non-rigid']:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
        else:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
        
        # Load the DICOM from the 2nd scan:
        dicom2 = pydicom.read_file(fpath2)
        
        # Store this toROICSrefsopuid:
        toROICSrefsopuids.append(dicom2.SOPInstanceUID)
    
        if refsopuid==dicom2.SOPInstanceUID:
            print('\n- The DICOM\'s SOP Instance UID is the same as the',\
                  'ROI\'s Referenced SOP Instance UID:\n', refsopuid)
        else:
            print('\n- Changing the Referenced SOP Instance UID from:\n', \
                  refsopuid, '\nto:\n', dicom2.SOPInstanceUID)
            
            roi2.ROIContourSequence[0].ContourSequence[i]\
            .ContourImageSequence[0].ReferencedSOPInstanceUID \
            = copy.deepcopy(dicom2.SOPInstanceUID)    
        
    # Update the dictionary changesDict:
    changesDict.update({'fromROICS.CS.CIS.RefSOPInstanceUIDs':\
                      fromROICSrefsopuids})
    changesDict.update({'toROICS.CS.CIS.RefSOPInstanceUIDs':\
                      toROICSrefsopuids})
            
            
    """
    Make optional changes to the Structure Set Date and Time to help identify 
    the ROI Collection from others:
    """
    # Set Date and Structure Set Time.
    # - Update the Structure Set Date:
    #roi2.StudyDate = time.strftime("%Y%m%d", time.gmtime())
    #newDate = time.strftime("%Y%m%d", time.gmtime())
    
    if roi2.StructureSetDate==newDate:
        print('\nThe Structure Set Date is today:\n', newDate)
    else:
        print('\nChanging the Structure Set Date from:\n', \
              roi2.StructureSetDate, '\nto:\n', newDate)
        
        # Update the dictionary xnatDict:
        xnatDict.update({'changesToRoi':\
                         {'fromStructureSetDate':roi2.StructureSetDate}})
    
        roi2.StructureSetDate = newDate
        
        # Update the dictionary changesDict:
        changesDict.update({'toStructureSetDate':roi2.StructureSetDate})
    
    # - Update the Structure Set Time:
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    print('\nChanging the Structure Set Time from:\n', roi2.StructureSetTime,\
          '\nto:\n', newTime)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromStructureSetTime':roi2.StructureSetTime})
    
    roi2.StructureSetTime = newTime
    
    # Update the dictionary changesDict:
    changesDict.update({'toStructureSetTime':roi2.StructureSetTime})
    
    
    
    """
    Create a new filename for the DICOM-RTSTRUCT file.  Start by getting the 
    first part of the filename of the original DICOM-RTSTRUCT file:
    """
    #origRoiDir, origRoiFname = os.path.split(roiDict1['Roi fpath'])
    fromRoiFname = xnatDict['fromRoiFname']
    
    """
    The file name format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
    "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
    or similar, find the first "_" character and use the characters that 
    precede it in the new file name.
    """
    # Find the index of the '_' character in the file name:
    #ind = origRoiFname.index('_')
    ind = fromRoiFname.index('_')
    
    """
    The new filename will be the same string of characters in origRoiFname
    including the '_', followed by the recently modified Structure Set Date
    and Structure Set Time, followed by some characters that indicate if the
    positions were "corrected" (e.g. 'zPosCorr'):
    """
    # If xnatDict['scanPosCorr'] is not equal to zero, add the character in
    # xnatDict['searchBy'] (e.g. 'z') followed by 'PosCorr':
    if xnatDict['scanPosCorr']: 
        toRoiFname = fromRoiFname[0:ind+1] + roi2.StructureSetDate + '_' \
                    + roi2.StructureSetTime + '_' + xnatDict['searchBy'] \
                    + 'PosCorr.dcm'
    else:
        toRoiFname = fromRoiFname[0:ind+1] + roi2.StructureSetDate + '_' \
                    + roi2.StructureSetTime + '.dcm'
    
    # The directory where the ROI is to be saved to:
    toRoiDir = xnatDict['toRoiDir']
    
    """
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
    else:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
    """
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    toRoiFpath = os.path.join(toRoiDir, toRoiFname)
    
   
    # The dictionary xnatDict will be updated (below) and exported to a JSON 
    # file in the path defined by xnatDict['toRoiDir'].
    # - Create a filename for the JSON file:
    """
    Use the same filename as the RTSTRUCT file but change the file extension.
    """
    jsonFname = os.path.splitext(toRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    jsonFpath = os.path.join(xnatDict['toRoiDir'], jsonFname)
    
    
    # Update the dictionary xnatDict:
    xnatDict.update({'changesToRoi':changesDict})
    xnatDict.update({'toRoiFname':toRoiFname})
    xnatDict.update({'toRoiFpath':toRoiFpath})
    xnatDict.update({'jsonFname':jsonFname})
    xnatDict.update({'jsonFpath':jsonFpath})
    
    
    # Print some useful info:
    #fromRFORSrefsopuids
    #toRFORSrefsopuids
    
    # Cycle through each Reference SOP Instance UID for each Contour Image
    # Sequence:
    for i in range(len(fromRFORSrefsopuids)):
        refsopuid = fromRFORSrefsopuids[i]
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # Get the slice no, Image Position, Slice Location and distances in
        # x, y, z and r for the 1st and 2nd DICOM scan (collection):
        sliceNo1 = dicomDict1_to_2[D1to2keys[ind]]['Dicom slice no']
        sliceNo2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom slice no']
        
        imPos1 = dicomDict1_to_2[D1to2keys[ind]]['Image pos']
        imPos2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Image pos']
        
        sliceLoc1 = dicomDict1_to_2[D1to2keys[ind]]['Slice loc']
        sliceLoc2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Slice loc']
    
        dx = dicomDict1_to_2[D1to2keys[ind]]['x2 - x1']
        dy = dicomDict1_to_2[D1to2keys[ind]]['y2 - y1']
        dz = dicomDict1_to_2[D1to2keys[ind]]['z2 - z1']
        dr = dicomDict1_to_2[D1to2keys[ind]]['r2 - r1']
        
        print(f'\n\nContour Image Sequence no {i+1}:\n', \
              f'\n   The contour(s) was copied from slice no {sliceNo1} in',\
              f'the original DICOM scan (collection)',\
              f'\n   to slice no {sliceNo2} in the new scan (collection):',\
              f'\n   - from Image Position {imPos1}',\
              f'\n                      to {imPos2}',\
              f'\n   - from Slice Location {sliceLoc1} to {sliceLoc2}',\
              f'\n   - with distances:',\
              f'\n       dx = {dx}',\
              f'\n       dy = {dy}',\
              f'\n       dz = {dz}',\
              f'\n       dr = {dr}')

    
    
    # Export the dictionary xnatDict to a JSON file:
    json.dump(xnatDict, open(jsonFpath, 'w'))
    
    print('\nJSON file created:\n', jsonFpath, '\n\n')
    
    
    # Save the new DICOM-RTSTRUCT file:
    roi2.save_as(toRoiFpath)
    
    #return roi2
    #return xnatDict
    #return roi2, xnatDict
    return roi1, roi2, xnatDict # 10/03/2020
    
