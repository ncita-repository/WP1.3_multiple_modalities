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
    from WalkThroughAttributes import WalkThroughAttributes
    
    # Load the DICOM-RTSTRUCT file whose ROIs correspond to the 1st DICOM scan:
    #roi1 = pydicom.dcmread(roiDict1['Roi fpath'])
    roi1 = pydicom.dcmread(xnatDict['fromRoiFpath'])
    
    # Use roi1 as a template for roi2:
    roi2 = copy.deepcopy(roi1)
    
    # Get the keys in roiDict1 and dicomDict1_to_2:
    #Rkeys = list(roiDict1.keys())
    D1to2keys = list(dicomDict1_to_2.keys())
    
    
    # ***********************************************************************
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
    
    # The list of metadata to be altered with are the following in pairs with
    # the DICOM tag followed by the DICOM-RTSTRUCT tag name in each pair:
    tagPairs = [['StudyInstanceUID', ['StudyInstanceUID']
                                      ], \
    
                ['StudyInstanceUID', ['ReferencedFrameOfReferenceSequence', 0, \
                                      'RTReferencedStudySequence', 0, \
                                      'ReferencedSOPInstanceUID']
                                      ], \
    
                ['SeriesInstanceUID', ['ReferencedFrameOfReferenceSequence', 0, \
                                      'RTReferencedStudySequence', 0, \
                                       'RTReferencedSeriesSequence', 0, \
                                       'SeriesInstanceUID']
                                       ], \
    
                ['FrameOfReferenceUID', ['FrameOfReferenceUID']
                                         ], \
    
                ['FrameOfReferenceUID', ['StructureSetROISequence', 0, \
                                         'ReferencedFrameOfReferenceUID']
                                        ]#, \
    
                #['StudyDate', ['StudyDate']
                #               ], \
    
                #['StudyTime', ['StudyTime']
                #               ], \
    
                #['PatientName', ['PatientName']
                #                 ], \
    
                #['PatientID', ['PatientID']
                #               ], \
    
                #['PatientSex', ['PatientSex']
                #                ]
               ]
    
    # Make changes to the metadata in roi2 to reflect the corresponding 
    # metadata contained in the DICOM. Use WalkThroughAttributes() to run 
    # through each pair of tags in tagPairs, making changes to all attributes
    # that don't have multiple items (e.g. sequences).
    
    """
    Note:
    I was going to proceed with developing WalkThroughAttributes() so that it
    could iterate through a list of multiple itemed attributes (e.g. sequences)
    but since it took longer than expected to deal with single-itemed 
    attributes I reasoned that the time it would take to deal with multiple 
    items would potentially exceed the value added in doing so. 
    
    Moreover, it might not even be possible, or even if it is, the code will 
    likely be very tedious and messy. So this is something I may consider 
    continuing with later.
    """
    
    # Iterate for each tagPair in tagPairs:
    for tagPair in tagPairs:

        # Get the top level ROI tag name by copying the second tag in tagPair 
        # and further indexing the last item in that list:
        topLevelRoiTag = copy.deepcopy(tagPair)[1][-1] 
    
        # Use WalkThroughAttributes() to make changes for this pair:
        oldVal, newVal, roi = WalkThroughAttributes(dicom2, \
                                                    roi2, \
                                                    tagPair, 
                                                    debug=debug)
    
        #print('\n', topLevelRoiTag, 'before:\n', oldVal)
        #print('\n', topLevelRoiTag, 'after:\n', newVal)
    
        # Update the dictionary changesDict if changes were made:
        if oldVal!=newVal:
            changesDict.update({'old' + topLevelRoiTag : oldVal})
            changesDict.update({'new' + topLevelRoiTag : newVal})

    # ***********************************************************************
    # ***********************************************************************
    
    
    """
    By taking the above approach rather than the previous one, which involved
    making changes to each attribute one by one has resulted in the code 
    shrinking from ~130 lines of code to ~50 if changing only the essential
    tags in tagPairs, or from ~260 to ~50 if also changing the optional tags, 
    which at the moment are commented out.
    """
    
    

    """ 
    **************************************************************************
    The following section of code is commented out since Study Time is dealt
    with in tagPairs, even though it doesn't deal with decimals as was 
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
    **************************************************************************
    """
        
        
    """
    Now make essential changes to the multi-itemed attributes that are not
    currently dealt with in WalkThroughAttributes().
    """
    
    # ***********************************************************************
    # 1) Modify the Referenced Frame of Reference UID for each Contour Image 
    # Sequence:

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
    oldRFORSrefsopuids = []
    newRFORSrefsopuids = []
    
    # Cycle through each Contour Image Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ReferencedSOPInstanceUID
        
        # Store this oldRFORSrefsopuid:
        oldRFORSrefsopuids.append(refsopuid)
        
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
        
        # Store this newRFORSrefsopuid:
        newRFORSrefsopuids.append(dicom2.SOPInstanceUID)
    
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
    changesDict.update({'oldRFORS/RTRSS/RTRSS/CIS/ReferencedSOPInstanceUIDs':\
                      oldRFORSrefsopuids})
    changesDict.update({'newRFORS/RTRSS/RTRSS/CIS/ReferencedSOPInstanceUIDs':\
                      newRFORSrefsopuids})
    # ***********************************************************************
    # ***********************************************************************                
      
              
    
    # ***********************************************************************
    # 2) Modify the Referenced SOP Instance UID for each Contour Sequence:
    
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
    oldROICSrefsopuids = []
    newROICSrefsopuids = []
    
    # Cycle through each Contour Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ContourImageSequence[0]\
        .ReferencedSOPInstanceUID
        
        # Store this oldROICSrefsopuid:
        oldROICSrefsopuids.append(refsopuid)
        
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
        
        # Store this newROICSrefsopuid:
        newROICSrefsopuids.append(dicom2.SOPInstanceUID)
    
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
    changesDict.update({'oldROICS/CS/CIS/ReferencedSOPInstanceUIDs':\
                      oldROICSrefsopuids})
    changesDict.update({'newROICS/CS/CIS/ReferencedSOPInstanceUIDs':\
                      newROICSrefsopuids})
    
    # ***********************************************************************
    # ***********************************************************************
    
    
    
    
    
    """
    Now make non-essential changes that cannot be implemented in an iterative
    fashion but are still useful:
    """
    
    # ***********************************************************************
    # 1) Modify the Structure Set Label by adding 
    # "_copy_from_Scan_A_to_B_YYYYMMDD_HHMMSS" to the Structure Set Label of 
    # the original RTSTRUCT, where A and B are the Scan Labels from the two 
    # scans (collections):
    newDate = time.strftime("%Y%m%d", time.gmtime())
    newTime = time.strftime("%H%M%S", time.gmtime())
    
    #newLabel = roi1.StructureSetLabel + '_copy_' + newDate \
    #+ '_' + newTime
    
    newLabel = roi1.StructureSetLabel + '_copy_from_Scan_' \
    + xnatDict['oldDicomScanLab'] + '_to_' + xnatDict['newDicomScanLab'] \
    + '_' + newDate + '_' + newTime
    
    if xnatDict['scanPosCorr']: 
        newLabel = newLabel + '_' + xnatDict['searchBy'] + 'PosCorr'
    
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        newLabel = newLabel + '_' + regMethod + '-registered'
    
    print('\nChanging the Structure Set Label from:\n', \
          roi2.StructureSetLabel, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'oldStructureSetLabel':roi2.StructureSetLabel})
    
    roi2.StructureSetLabel = newLabel  
    
    # Update the dictionary changesDict:
    changesDict.update({'newStructureSetLabel':roi2.StructureSetLabel})
    
    # ***********************************************************************
    # ***********************************************************************
    
    

    # ***********************************************************************
    # 2) Modify the ROI Name in the Structure Set ROI Sequence:
    """
    Modify the ROI Name by making it the same as the new Structure Set Label.
    """
    
    print('\nChanging the ROI Name from:\n', \
          roi2.StructureSetROISequence[0].ROIName, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'oldSSROIS/ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
    
    roi2.StructureSetROISequence[0].ROIName = newLabel
    
    # Update the dictionary changesDict:
    changesDict.update({'newSSROIS/ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
                    
    # ***********************************************************************
    # ***********************************************************************
    
    
            
            
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
                         {'oldStructureSetDate':roi2.StructureSetDate}})
    
        roi2.StructureSetDate = newDate
        
        # Update the dictionary changesDict:
        changesDict.update({'newStructureSetDate':roi2.StructureSetDate})
    
    # - Update the Structure Set Time:
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    print('\nChanging the Structure Set Time from:\n', roi2.StructureSetTime,\
          '\nto:\n', newTime)
    
    # Update the dictionary changesDict:
    changesDict.update({'oldStructureSetTime':roi2.StructureSetTime})
    
    roi2.StructureSetTime = newTime
    
    # Update the dictionary changesDict:
    changesDict.update({'newStructureSetTime':roi2.StructureSetTime})
    
    
    
    """
    Create a new filename for the DICOM-RTSTRUCT file.  Start by getting the 
    first part of the filename of the original DICOM-RTSTRUCT file:
    """
    #origRoiDir, origRoiFname = os.path.split(roiDict1['Roi fpath'])
    oldRoiFname = xnatDict['fromRoiFname']
    
    """
    The file name format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
    "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
    or similar, find the first "_" character and use the characters that 
    precede it in the new file name.
    """
    # Find the index of the '_' character in the file name:
    #ind = origRoiFname.index('_')
    ind = oldRoiFname.index('_')
    
    """
    The new filename will be the same string of characters in oldRoiFname
    including the '_', followed by the recently modified Structure Set Date
    and Structure Set Time, followed by some characters that indicate if the
    positions were "corrected" (e.g. 'zPosCorr'):
    """
    
    newRoiFname = oldRoiFname[0:ind+1] + roi2.StructureSetDate + '_' \
                    + roi2.StructureSetTime
                    
    # If xnatDict['scanPosCorr'] is not equal to zero, add the character in
    # xnatDict['searchBy'] (e.g. 'z') followed by 'PosCorr':
    if xnatDict['scanPosCorr']: 
        newRoiFname = newRoiFname + '_' + xnatDict['searchBy'] + 'PosCorr'
    
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        newRoiFname = newRoiFname + '_' + regMethod + '-registered'
        
    # Finally add the file extension:
    newRoiFname = newRoiFname + '.dcm'
        
    # The directory where the ROI is to be saved to:
    toRoiDir = xnatDict['toRoiDir']
    
    """
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
    else:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
    """
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    toRoiFpath = os.path.join(toRoiDir, newRoiFname)
    
   
    # The dictionary xnatDict will be updated (below) and exported to a JSON 
    # file in the path defined by xnatDict['toRoiDir'].
    # - Create a filename for the JSON file:
    """
    Use the same filename as the RTSTRUCT file but change the file extension.
    """
    jsonFname = os.path.splitext(newRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    jsonFpath = os.path.join(xnatDict['toRoiDir'], jsonFname)
    
    
    # Update the dictionary xnatDict:
    xnatDict.update({'changesToRoi':changesDict})
    xnatDict.update({'newRoiFname':newRoiFname})
    xnatDict.update({'toRoiFpath':toRoiFpath})
    xnatDict.update({'jsonFname':jsonFname})
    xnatDict.update({'jsonFpath':jsonFpath})
    
    
    #************************************************************************
    # Print some useful info:
    
    #fromRFORSrefsopuids
    #toRFORSrefsopuids
    
    # Cycle through each Reference SOP Instance UID for each Contour Image
    # Sequence:
    for i in range(len(oldRFORSrefsopuids)):
        refsopuid = oldRFORSrefsopuids[i]
        
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
        
    #************************************************************************
    #************************************************************************

    
    
    # Export the dictionary xnatDict to a JSON file:
    json.dump(xnatDict, open(jsonFpath, 'w'))
    
    print('\nJSON file created:\n', jsonFpath, '\n\n')
    
    
    # Save the new DICOM-RTSTRUCT file:
    roi2.save_as(toRoiFpath)
    
    #return roi2
    #return xnatDict
    #return roi2, xnatDict
    return roi1, roi2, xnatDict # 10/03/2020
    
