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
                      
    DataDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths 
                    (see DownloadFilesForRoiCopy() for details)
                    
    RegMethod      - Type of registration to perform on the DICOM scan
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
    ToRoi           - ROI object that contains the contours in roi1 but 
                    modified for overlaying over the 2nd DICOM scan;
                   - A new DICOM-RTSTRUCT file will be saved to the location
                   specified
                   
    DataDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths
                    (see DownloadFilesForRoiCopy() for details)
                   - The dictionary is exported as a JSON file
    
"""

#def ModifyRoi(FromRoiFpath, ToRoiDir, FromDicoms, ToDicoms, Inds, DataDict):
def ModifyRoi(FromDicoms, ToDicoms, DataDict, FromKey, ToKey):
    # Import packages:
    import copy
    import time
    import json
    import pydicom
    import os
    from WalkThroughAttributes import WalkThroughAttributes
    from CreateDir import CreateDir
    #import ConvertNumpyDtype
    import importlib
    import ModifyRoiDirs
    importlib.reload(ModifyRoiDirs)
    from ModifyRoiDirs import ModifyRoiDirs
    
    # Get necessary into from DataDict:
    FromRoiFpath = DataDict[FromKey]['RoiFpath']
    FromSeriesNo = DataDict[FromKey]['SeriesNo']
    FromSliceNos = DataDict[FromKey]['SliceNos']
    
    #ToRoiDir = DataDict[ToKey]['RoiDir'] # this gets updated below!
    ToRoiDir = DataDict[ToKey]['RoiDir']
    ToSeriesNo = DataDict[ToKey]['SeriesNo']
    ToSliceNos = DataDict[ToKey]['SliceNos']
    #Inds = DataDict[ToKey]['Inds']
    SearchBy = DataDict['SearchBy']
    Shift = DataDict['Shift']
    RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug']
    
    # This is used to get slice numbers later on:
    #if 'Map' in ToKey:
    #    Inds = DataDict[ToKey]['Inds']
    #else:
    #    #Inds = DataDict[ToKey]['SliceNos']
    #    Inds = DataDict[FromKey]['SliceNos']
    #    #Inds = range(len(FromDicoms))
    if 'Mapping' in DataDict.keys():
        MappingInds = DataDict['Mapping']['MappingInds']
    
    #print('\nToKey =', ToKey)
    
    #print('\nDataDict[' + ToKey + ']:\n', DataDict[ToKey])
    
    if False:
        # First modify RoiDir in DataDict:
        DataDict = ModifyRoiDirs(DataDict)
        
        # Get the updated ToRoiDir:
        ToRoiDir = DataDict[ToKey]['RoiDir']
    
    #print('\nToKey =', ToKey)
    
    #print('\nDataDict[' + ToKey + ']:\n', DataDict[ToKey])
    
    # Load the DICOM-RTSTRUCT file whose ROIs correspond to the 1st DICOM scan:
    #roi1 = pydicom.dcmread(roiDict1['Roi fpath'])
    FromRoi = pydicom.dcmread(FromRoiFpath)
    
    # Use FromRoi as a template for ToRoi:
    ToRoi = copy.deepcopy(FromRoi)
    
    # Get the SOP Instance UIDs from the FromDicoms and ToDicoms.  If the key
    # 'Mapping' is in DataDict.keys(), replace FromDicoms with MappedDicoms:
    #if 'Mapping' in DataDict.keys():
    #    print('\nMapping in DataDict.keys()')
    #    MappedDicoms = [FromDicoms[ind] for ind in MappingInds]
    #    FromSopUids = [item.SOPInstanceUID for item in MappedDicoms]
    #else:
    #    FromSopUids = [item.SOPInstanceUID for item in FromDicoms]
    FromSopUids = [item.SOPInstanceUID for item in FromDicoms]
    ToSopUids = [item.SOPInstanceUID for item in ToDicoms]
    
    if True:#False:
        print(f'\nlen(FromDicoms)  = {len(FromDicoms)}')
        print(f'len(FromSopUids) = {len(FromSopUids)}')
        print(f'len(ToDicoms)    = {len(ToDicoms)}')
        print(f'len(ToSopUids)   = {len(ToSopUids)}')
        print(f'len(MappingInds) = {len(MappingInds)}')
        
    
    print('\n--- Modifying a copy of the ROI Collection from', 
          DataDict[FromKey]['SeriesDesc'], 'to', \
          DataDict[ToKey]['SeriesDesc'], '...')
    
    
    # ***********************************************************************
    # First make changes to ToRoi that only need to be made once (don't need
    # to be made for each contour).
    
    
    # Create a dictionary to store the changes made to the ROI that will 
    # later be added to xnatDict:
    ChangesDict = {}
    
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
                
    
    """ 
    April 3:  
        
    "Missing dependencies" error when trying to upload ROIs to XNAT so added 
    additional tags that I thought were non-essential are are commented out  
    above.  But error remained.
    """
    
    
    """
    Make changes to the metadata in ToRoi to reflect the corresponding 
    metadata contained in the DICOM. Use WalkThroughAttributes() to run 
    through each pair of tags in tagPairs, making changes to all attributes
    that don't have multiple items (e.g. sequences).
    """
    
    """
    March 23 Note:
        
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
    
        # Use WalkThroughAttributes() to make changes for this pair using the
        # first DICOM in ToDicoms:
        oldVal, newVal, roi = WalkThroughAttributes(ToDicoms[0], \
                                                    ToRoi, \
                                                    tagPair, 
                                                    Debug)
    
        #print('\n', topLevelRoiTag, 'before:\n', oldVal)
        #print('\n', topLevelRoiTag, 'after:\n', newVal)
    
        # Update the dictionary ChangesDict if changes were made:
        if oldVal!=newVal:
            ChangesDict.update({'old' + topLevelRoiTag : oldVal})
            ChangesDict.update({'new' + topLevelRoiTag : newVal})

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

    if ToRoi.StudyTime==newTime:
        print('\nThe Study Time are the same:\n', ToRoi.StudyTime)
    else:
        print('\nChanging the Study Time from:\n', ToRoi.StudyTime, '\nto:\n', \
              newTime)
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'fromStudyTime':ToRoi.StudyTime})
        
        ToRoi.StudyTime = newTime
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'toStudyTime':ToRoi.StudyTime})
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
    sequences = ToRoi.ReferencedFrameOfReferenceSequence[0]\
                .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                .ContourImageSequence
    
    if True:#Debug:
        print('\n--- Looping through each Contour Image Sequence to modify'\
              'the Referenced SOP Instance UIDs..')
    
    
    # Store the original and modified Referenced SOP Instance UIDs from the
    # Referenced Frame Of Reference Sequence:
    oldRFORSrefsopuids = []
    newRFORSrefsopuids = []
    
    # Store the indeces that link each Ref SOP Instance UID from the Contour
    # Image Sequence to the corresponding SOP Instance UID from the DICOMs:
    SopUids_ind_of_RefSopUids_in_CIS = []
    
    # Cycle through each Contour Image Sequence:
    for i in range(len(sequences)):
        print(f'\nContour Image Sequence #{i}:\n')
        print('\nOriginal Ref SOP Instance UID =\n', \
              ToRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i].ReferencedSOPInstanceUID)   
        
        # Get the Referenced SOP Instance UID for this sequence:
        RefSopUid = sequences[i].ReferencedSOPInstanceUID
        
        # Store this oldRFORSrefsopuid:
        oldRFORSrefsopuids.append(RefSopUid)
        
        # Find the index of the DICOM in FromDicoms whose SOP Instance UID
        # matches the ROI's Referenced SOP Instance UID:
        ind = FromSopUids.index(RefSopUid)
        
        # Store this ind:
        SopUids_ind_of_RefSopUids_in_CIS.append(ind)
        
        # Use ind to get the SOP Instance UID of ToDicoms that is to replace
        # the old Ref SOP Instance UID:
        SopUid = FromSopUids[ind] 
        
        print('\nTo be replaced with SOP Instance UID =\n', SopUid)     
        
        
        # Change the ROI's Referenced SOP Instance UID to that of the DICOM's 
        # SOP Instance UID:
    
        ToRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i].ReferencedSOPInstanceUID \
                 = copy.deepcopy(SopUid)
        
        print('\nNew Ref SOP Instance UID =\n', \
              ToRoi.ReferencedFrameOfReferenceSequence[0]\
                 .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                 .ContourImageSequence[i].ReferencedSOPInstanceUID) 
        
        # Store this newRFORSrefsopuid:
        newRFORSrefsopuids.append(SopUid)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'oldRFORS/RTRSS/RTRSS/CIS/ReferencedSOPInstanceUIDs':\
                      oldRFORSrefsopuids})
    ChangesDict.update({'newRFORS/RTRSS/RTRSS/CIS/ReferencedSOPInstanceUIDs':\
                      newRFORSrefsopuids})
    # ***********************************************************************
    # ***********************************************************************                
      
              
    
    # ***********************************************************************
    # 2) Modify the Referenced SOP Instance UID for each Contour Sequence:
    
    # Get the Contour Sequences:
    sequences = ToRoi.ROIContourSequence[0].ContourSequence
    
    if Debug:
        print('\n--- Looping through each Contour Sequence to modify the '\
          'Referenced SOP Instance UIDs..')
       
    # Store the original and modified Referenced SOP Instance UIDs from the
    # ROI Contour Sequence:
    oldROICSrefsopuids = []
    newROICSrefsopuids = []
    
    # Store the indeces that link each Ref SOP Instance UID from the Contour
    # Sequence to the corresponding SOP Instance UID from the DICOMs:
    SopUids_ind_of_RefSopUids_in_CS = []
    
    # Cycle through each Contour Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        RefSopUid = sequences[i].ContourImageSequence[0]\
                    .ReferencedSOPInstanceUID
        
        # Store this oldROICSrefsopuid:
        oldROICSrefsopuids.append(RefSopUid)
        
        # Find the index of the DICOM in FromDicoms whose SOP Instance UID
        # matches the ROI's Referenced SOP Instance UID:
        ind = FromSopUids.index(RefSopUid)
        
        # Store this ind:
        SopUids_ind_of_RefSopUids_in_CS.append(ind)
        
        # Use ind to get the SOP Instance UID:
        SopUid = FromSopUids[ind] 

        
        # Change the ROI's Referenced SOP Instance UID to that of the DICOM's 
        # SOP Instance UID:
        ToRoi.ROIContourSequence[0].ContourSequence[i]\
                 .ContourImageSequence[0].ReferencedSOPInstanceUID \
                 = copy.deepcopy(SopUid) 
        
        # Store this newROICSrefsopuid:
        newROICSrefsopuids.append(SopUid)
        
    # Update the dictionary ChangesDict:
    ChangesDict.update({'oldROICS/CS/CIS/ReferencedSOPInstanceUIDs':\
                      oldROICSrefsopuids})
    ChangesDict.update({'newROICS/CS/CIS/ReferencedSOPInstanceUIDs':\
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
    
    newLabel = FromRoi.StructureSetLabel + '_copy_from_SeriesNo_' \
               + FromSeriesNo + '_to_' + ToSeriesNo + '_' \
               + newDate + '_' + newTime
    
    if Shift: 
        #newLabel = newLabel + '_' + SearchBy + 'PosCorr'
        newLabel = newLabel + '_Shifted'
    
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        newLabel = newLabel + '_' + RegMethod + '-Reg'
    
    if True:#Debug:
        print('\n--- Changing the Structure Set Label from:\n', \
          ToRoi.StructureSetLabel, '\nto:\n', newLabel)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'oldStructureSetLabel':ToRoi.StructureSetLabel})
    
    # Modify the Structure Set Label:
    ToRoi.StructureSetLabel = newLabel  
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'newStructureSetLabel':ToRoi.StructureSetLabel})
    
    # ***********************************************************************
    # ***********************************************************************
    
    

    # ***********************************************************************
    # 2) Modify the ROI Name in the Structure Set ROI Sequence:
    """
    Modify the ROI Name by making it the same as the new Structure Set Label.
    """
    
    if True:#Debug:
        print('\n--- Changing the ROI Name from:\n', \
              ToRoi.StructureSetROISequence[0].ROIName, '\nto:\n', newLabel)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'oldSSROIS/ROIName':\
                        ToRoi.StructureSetROISequence[0].ROIName})
    
    # Modify the ROI Name:
    ToRoi.StructureSetROISequence[0].ROIName = newLabel
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'newSSROIS/ROIName':\
                        ToRoi.StructureSetROISequence[0].ROIName})
                    
    # ***********************************************************************
    # ***********************************************************************
    
    
            
            
    """
    Make optional changes to the Structure Set Date and Time to help identify 
    the ROI Collection from others:
    """
    # Set Date and Structure Set Time.
    # - Update the Structure Set Date:
    #ToRoi.StudyDate = time.strftime("%Y%m%d", time.gmtime())
    #newDate = time.strftime("%Y%m%d", time.gmtime())
    
    if ToRoi.StructureSetDate==newDate:
        if True:#Debug:
            print('\n--- The Structure Set Date is today:\n', newDate)
    else:
        if True:#Debug:
            print('\n--- Changing the Structure Set Date from:\n', \
                  ToRoi.StructureSetDate, '\nto:\n', newDate)
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'oldStructureSetDate':ToRoi.StructureSetDate})
    
        # Modify the Structure Set Date:
        ToRoi.StructureSetDate = newDate
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'newStructureSetDate':ToRoi.StructureSetDate})
    
    # - Update the Structure Set Time:
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    if True:#Debug:
        print('\n--- Changing the Structure Set Time from:\n', \
              ToRoi.StructureSetTime, '\nto:\n', newTime)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'oldStructureSetTime':ToRoi.StructureSetTime})
    
    # Modify the Structure Set Time:
    ToRoi.StructureSetTime = newTime
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'newStructureSetTime':ToRoi.StructureSetTime})
    
    
    
    """
    Create a new filename for the DICOM-RTSTRUCT file.  
    Start by getting the filename of the CopyFrom DICOM-RTSTRUCT file:
    """
    # Get the CopyFrom filename:
    FromRoiFname = os.path.basename(FromRoiFpath)
    
    # Get filename (with extension):
    FromRoiFname = os.path.splitext(FromRoiFname)
    
    # Get filename without extension:
    FromRoiFname = FromRoiFname[0]
    
    """
    The file name format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
    "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
    or similar, find the first "_" character and use the characters that 
    precede it in the new file name.
    """
    # Find the index of the '_' character in the file name:
    #ind = origRoiFname.index('_')
    #ind = FromRoiFname[0].index('_')
    ind = FromRoiFname.index('_')
    
    """
    ToRoiFname will be the same string of characters in FromRoiFname
    including the '_', followed by the recently modified Structure Set Date
    and Structure Set Time, followed by some characters that indicate if the
    positions were "corrected" (e.g. 'zPosCorr'):
    """
    
    #ToRoiFname = FromRoiFname[0][0:ind+1] + ToRoi.StructureSetDate \
    #                 + '_' + ToRoi.StructureSetTime
    #ToRoiFname = FromRoiFname[0][0:ind+1]
    ToRoiFname = FromRoiFname[0:ind+1]
    
    """
    Note: The first index [0] gets the filename without extension. The second
    index [0:ind+1] gets all characters to and including '_'. 
    """
    
    # If Shift is not equal to zero, add the character in SearchBy 
    # (e.g. 'z') followed by 'PosCorr':
    if Shift: 
        ToRoiFname = ToRoiFname + '_' + SearchBy + 'PosCorr'
    
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        ToRoiFname = ToRoiFname + '_' + RegMethod + '-Reg'
        
    
    print('\nFromRoiFname =', FromRoiFname)
    print('ToRoi.StructureSetDate =', ToRoi.StructureSetDate)
    print('ToRoi.StructureSetTime =', ToRoi.StructureSetTime)
    
    # Add the ROI Structure Set Date and Time:
    ToRoiFname = FromRoiFname + '_' + ToRoi.StructureSetDate \
                     + '_' + ToRoi.StructureSetTime
        
    # Finally add the file extension:
    ToRoiFname = ToRoiFname + '.dcm'
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    ToRoiFpath = os.path.join(ToRoiDir, ToRoiFname)
    
   
    """
    The dictionary DataDict will be updated (below) and exported to a JSON 
    file in the path defined by DataDict[ToKey]['RoiDir'].
    """
    # - Create a filename for the JSON file, using the same base filename as 
    # the RTSTRUCT file:
    JsonFname = os.path.splitext(ToRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    JsonFpath = os.path.join(ToRoiDir, JsonFname)
    
    # Create the directory ToRoiDir:
    CreateDir(ToRoiDir, Debug)
    
    #print('\nChangesDict:\n', ChangesDict)
    
    #print('\nToKey =', ToKey)
    
    #print('\nDataDict[' + ToKey + ']:\n', DataDict[ToKey])
    
    # Update the dictionary DataDict:
    DataDict[ToKey].update({'ChangesToRoi':ChangesDict})
    #DataDict[ToKey].update({'newRoiFname':newRoiFname})
    DataDict[ToKey].update({'RoiFpath':ToRoiFpath})
    #DataDict[ToKey].update({'jsonFname':jsonFname})
    DataDict[ToKey].update({'JsonFpath':JsonFpath})
    
    
    #************************************************************************
    # Print some useful info:
    
    """
    March 28 Note:
        The code below will need to be checked and corrected.
    """
    
    #fromRFORSrefsopuids
    #toRFORSrefsopuids
    
    
    # Cycle through each Reference SOP Instance UID for each Contour Image
    # Sequence:
    for i in range(len(oldRFORSrefsopuids)):
        RefSopUid = oldRFORSrefsopuids[i]
        
        # Get the index of the this RefSopUid:
        #ToInd = ToRefSopUids.index(RefSopUid)
        
        # Find the index of the DICOM in FromDicoms whose SOP Instance UID
        # matches this refsopuid:
        FromInd = FromSopUids.index(RefSopUid)
        
        # Get the slice numbers:
        FromSliceNo = FromSliceNos[FromInd]
        ToSliceNo = ToSliceNos[MappingInds[i]]
        ##ToSliceNo = ToSliceNos[ToInd]
        
        if False:
            print('\nFromSliceNo =', FromSliceNo)
            print(f'\MappingInds (N = {len(MappingInds)}) =', MappingInds)
        
        # Get the slice no, Image Position, Slice Location and distances in
        # x, y, z and r for the 1st and 2nd DICOM scan (collection):
        #sliceNo1 = dicomDict1_to_2[D1to2keys[ind]]['Dicom slice no']
        #sliceNo2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom slice no']
        #ToSliceNo = Inds[FromSliceNo] (March 31)
        
        #print(f'\nFromSliceNo = {FromSliceNo}, ToSliceNo = {ToSliceNo}')
    
        
        #FromImPos = FromDicoms[FromSliceNo].ImagePositionPatient[dim]
        #ToImPos = ToDicoms[ToSliceNo].ImagePositionPatient[dim]
        #FromImPos = FromDicoms[FromSliceNo].ImagePositionPatient
        #ToImPos = ToDicoms[ToSliceNo].ImagePositionPatient
        FromImPos = FromDicoms[i].ImagePositionPatient #(March 31)
        ToImPos = ToDicoms[i].ImagePositionPatient #(March 31)
        #FromImPos = FromDicoms[FromInd].ImagePositionPatient
        #ToImPos = ToDicoms[ToInd].ImagePositionPatient 
        
        if SearchBy=='x':
            dim = 0
        elif SearchBy=='y':
            dim = 1
        elif SearchBy=='z':
            dim = 2
        """
        Above code doesn't deal with case where SearchBy='r'!
        """
        
        FromScanPos = FromImPos[dim]
        ToScanPos = ToImPos[dim]
        
        #dPos = DataDict[ToKey]['AbsMinDpos']
        dx = round(ToImPos[0] - FromImPos[0], 2)
        dy = round(ToImPos[1] - FromImPos[1], 2)
        dz = round(ToImPos[2] - FromImPos[2], 2)
        dr = round((dx**2 + dy**2 + dz**2)**.5, 2)
        
        if True:#Debug:
            print(f'\n\nContour Image Sequence #{i+1}:\n', \
                  f'\n   The contour(s) was copied from slice #{FromSliceNo}' \
                  f' in DICOM series #{FromSeriesNo} to slice #{ToSliceNo}' \
                  f' in DICOM series #{ToSeriesNo}:\n', \
                  f'\n   - from Image Position', FromImPos, \
                  f'\n                      to', ToImPos, '\n', \
                  f'\n   - from Scan Position', FromScanPos, 'mm to', \
                  ToScanPos, 'mm\n', \
                  f'\n   - with positional differences:',\
                  f'\n       dx = {dx} mm',\
                  f'\n       dy = {dy} mm',\
                  f'\n       dz = {dz} mm',\
                  f'\n       dr = {dr} mm')
        
    #************************************************************************
    #************************************************************************
    
    
    
    
    
    
    # Export the dictionary DataDict to a JSON file:
    json.dump(DataDict, open(JsonFpath, 'w'))
    
    print('\n--- JSON file created:\n', JsonFpath, '\n\n')
    
    
    # Save the new DICOM-RTSTRUCT file:
    ToRoi.save_as(ToRoiFpath)
    
    return FromRoi, ToRoi, DataDict # 10/03/2020
    
