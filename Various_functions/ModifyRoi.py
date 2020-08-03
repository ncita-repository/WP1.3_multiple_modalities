# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:48:48 2020

@author: ctorti
"""


"""
Function:
    ModifyRoi()
    
Purpose:
    Copy ROI Collection from source and modify so they may be overlaid on
    target DICOMs.
    
Input:
    SourceDicoms - array of DICOM objects whose ROIs are to be copied
    
    TargetDicoms - array of DICOM objects that the copied ROIs are to be 
                   overlaid on
                      
    DataDict     - Dictionary containing several variables 
                  (see DownloadFilesForRoiCopy() for details)
                    
    SourceKey    - key in DataDict of source data
                     
    TargetKey    - key in DataDict of target data
    
    
Returns:
    SourceRoi    - ROI object belonging to SourceDicoms
                  
    TargetRoi    - ROI object belonging to TargetDicoms (copied from SourceRoi)
                   
    DataDict     - updated dictionary
    
"""

#def ModifyRoi(SourceRoiFpath, TargetRoiDir, SourceDicoms, TargetDicoms, Inds, DataDict):
def ModifyRoi(SourceDicoms, TargetDicoms, DataDict, SourceKey, TargetKey):
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
    
    print('\nlen(SourceDicoms) =', len(SourceDicoms))
    print('\nlen(TargetDicoms) =', len(TargetDicoms))
    
    # Get necessary into from DataDict:
    SourceRoiFpath = DataDict[SourceKey]['RoiFpath']
    SourceSopUids = DataDict[SourceKey]['SopUids']
    SourceSeriesNo = DataDict[SourceKey]['SeriesNo']
    SourceSliceNos = DataDict[SourceKey]['SliceNos']
    
    #TargetRoiDir = DataDict[TargetKey]['RoiDir'] # this gets updated below!
    TargetRoiDir = DataDict[TargetKey]['RoiDir']
    TargetSopUids = DataDict[TargetKey]['SopUids']
    TargetSeriesNo = DataDict[TargetKey]['SeriesNo']
    TargetSliceNos = DataDict[TargetKey]['SliceNos']
    #Inds = DataDict[TargetKey]['Inds']
    SearchBy = DataDict['SearchBy']
    Shift = DataDict['Shift']
    RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug']
    
    # This is used to get slice numbers later on:
    #if 'Map' in TargetKey:
    #    Inds = DataDict[TargetKey]['Inds']
    #else:
    #    #Inds = DataDict[TargetKey]['SliceNos']
    #    Inds = DataDict[SourceKey]['SliceNos']
    #    #Inds = range(len(SourceDicoms))
    #print('\nKeys in DataDict:\n', DataDict.keys())
    if 'ToRemapKey' in DataDict.keys():
        #MappingInds = DataDict['Mapping']['MappingInds']
        SourceToTargetMapping = DataDict['SourceToTargetMapping']['MappingInds']
        TargetToSourceMapping = DataDict['TargetToSourceMapping']['MappingInds']
        
        # The key of the series that was remapped:
        RemappedKey = DataDict['RemappedKey']
        
    else:
        print('\n\nNo remapped data found.')
        
        #SourceToTargetMapping = DataDict['SourceToTargetMapping']['MappingInds']
        #TargetToSourceMapping = DataDict['TargetToSourceMapping']['MappingInds']
        
    
    #print('\nTargetKey =', TargetKey)
    
    #print('\nDataDict[' + TargetKey + ']:\n', DataDict[TargetKey])
    
    if False:
        # First modify RoiDir in DataDict:
        DataDict = ModifyRoiDirs(DataDict)
        
        # Get the updated TargetRoiDir:
        TargetRoiDir = DataDict[TargetKey]['RoiDir']
    
    #print('\nTargetKey =', TargetKey)
    
    #print('\nDataDict[' + TargetKey + ']:\n', DataDict[TargetKey])
    
    """ Load the Source ROIs: """
    SourceRoi = pydicom.dcmread(SourceRoiFpath)
    
    """ Use SourceRoi as a template for TargetRoi: """
    TargetRoi = copy.deepcopy(SourceRoi)
    
    if False:
        # Get the SOP Instance UIDs from the SourceDicoms and TargetDicoms.  If the key
        # 'Mapping' is in DataDict.keys(), replace SourceDicoms with MappedDicoms:
        #if 'Mapping' in DataDict.keys():
        #    print('\nMapping in DataDict.keys()')
        #    MappedDicoms = [SourceDicoms[ind] for ind in MappingInds]
        #    SourceSopUids = [item.SOPInstanceUID for item in MappedDicoms]
        #else:
        #    SourceSopUids = [item.SOPInstanceUID for item in SourceDicoms]
        SourceSopUids = [item.SOPInstanceUID for item in SourceDicoms]
        TargetSopUids = [item.SOPInstanceUID for item in TargetDicoms]
    
    if False:
        print(f'\nlen(SourceDicoms)  = {len(SourceDicoms)}')
        print(f'len(SourceSopUids) = {len(SourceSopUids)}')
        print(f'len(TargetDicoms)    = {len(TargetDicoms)}')
        print(f'len(TargetSopUids)   = {len(TargetSopUids)}')
        #print(f'len(MappingInds) = {len(MappingInds)}')
        print(f'len(SourceToTargetMapping) = {len(SourceToTargetMapping)}')
        print(f'len(TargetToSourceMapping) = {len(TargetToSourceMapping)}')
        
    
    print('\n--- Modifying a copy of the ROI Collection from', 
          DataDict[SourceKey]['SeriesDesc'], 'to', \
          DataDict[TargetKey]['SeriesDesc'], '...')
    
    
    """
    ***********************************************************************
    STEP 1:
        Make changes to TargetRoi that only need to be made once (don't need
        to be made for each contour).
    ***********************************************************************
    """
    
    
    # Create a dictionary to store the changes made to the ROI that will 
    # later be added to xnatDict:
    ChangesDict = {}
    
    # The list of metadata to be altered with are the following in pairs with
    # the DICOM tag followed by the DICOM-RTSTRUCT tag name in each pair:
    TagPairs = [['StudyInstanceUID', ['StudyInstanceUID']
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
    Make changes to the metadata in TargetRoi to reflect the corresponding 
    metadata contained in the DICOM. Use WalkThroughAttributes() to run 
    through each pair of tags in TagPairs, making changes to all attributes
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
    
    # Iterate for each TagPair in TagPairs:
    for TagPair in TagPairs:

        # Get the top level ROI tag name by copying the second tag in TagPair 
        # and further indexing the last item in that list:
        TopLevelRoiTag = copy.deepcopy(TagPair)[1][-1] 
    
        # Use WalkThroughAttributes() to make changes for this pair using the
        # first DICOM in TargetDicoms:
        OldVal, NewVal, roi = WalkThroughAttributes(TargetDicoms[0], \
                                                    TargetRoi, \
                                                    TagPair, 
                                                    Debug)
    
        #print('\n', TopLevelRoiTag, 'before:\n', OldVal)
        #print('\n', TopLevelRoiTag, 'after:\n', NewVal)
    
        # Update the dictionary ChangesDict if changes were made:
        if OldVal!=NewVal:
            ChangesDict.update({'Old' + TopLevelRoiTag : OldVal})
            ChangesDict.update({'New' + TopLevelRoiTag : NewVal})

    """
    ***********************************************************************
    END OF STEP #1.
    ***********************************************************************
    """
    
    
    """
    By taking the above approach rather than the previous one, which involved
    making changes to each attribute one by one has resulted in the code 
    shrinking from ~130 lines of code to ~50 if changing only the essential
    tags in TagPairs, or from ~260 to ~50 if also changing the optional tags, 
    which at the moment are commented out.
    """
    
    

    """ 
    The following section of code is commented out since Study Time is dealt
    with in TagPairs, even though it doesn't deal with decimals as was 
    previously done.  If this turns out to be a problem, Study Time can be 
    dealt with outside of the above for loop.
    """  
    if False:
        # - Change the Study Time:
        """ Note: Some DICOM Study Times have a decimal with trailing zeros,
            e.g. '164104.000000', so convert to float, then integer, then back 
            to string. 
        """
        
        NewTime = str(int(float(copy.deepcopy(TargetDicoms[0].StudyTime))))
    
        if TargetRoi.StudyTime==NewTime:
            print('\nThe Study Time are the same:\n', TargetRoi.StudyTime)
        else:
            print('\nChanging the Study Time from:\n', TargetRoi.StudyTime, '\nto:\n', \
                  NewTime)
            
            # Update the dictionary ChangesDict:
            ChangesDict.update({'OldStudyTime':TargetRoi.StudyTime})
            
            TargetRoi.StudyTime = NewTime
            
            # Update the dictionary ChangesDict:
            ChangesDict.update({'NewStudyTime':TargetRoi.StudyTime})
   
        
        
    """
    ***********************************************************************
    STEP 2:
        Make essential changes to the multi-itemed attributes that are not
        currently dealt with in WalkThroughAttributes().
    ***********************************************************************
    
    
    ***********************************************************************
    STEP 2a:
        Modify the Referenced SOP Instance UID in each Contour Image Sequence. 
    ***********************************************************************
    
    For each Referenced SOP Instance UID in each Contour Image Sequence in 
    SourceRoi get the index of the SourceDicom with matching SOP Instance UID. 
    
    """        
    
    """ April 14: SOP Instance UIDs now stored in DataDict (see above): """
    # Get the SOP Instance UIDs from the SourceDicoms and TargetDicoms:
    #SourceSopUids = [item.SOPInstanceUID for item in SourceDicoms]
    #TargetSopUids = [item.SOPInstanceUID for item in TargetDicoms]
    
    #print('\nSourceSopUids =\n', SourceSopUids)
    #print('\nTargetSopUids =\n', TargetSopUids)
                
    # Store the Referenced SOP Instance UIDs:
    SourceCISRefSopUids = []
    
    # Store the indeces that map each Source Contour Image Sequence's 
    # Referenced SOP Instance UID to the matching SOP Instance UID of the 
    # Source DICOMs:
    #SourceCISRefSopUidToSourceSopUidMapping = []
    
    # Get the Source Contour Image Sequences:
    Sequences = SourceRoi.ReferencedFrameOfReferenceSequence[0]\
                .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
                .ContourImageSequence


    """ For debugging: """
    if False:
        print('\nSource SOP UIDs:\n')
        [print(item) for item in SourceSopUids]
        print('\n\nSource CIS Ref SOP UIDs:\n')
        [print(item) for item in [Sequences[i].ReferencedSOPInstanceUID for i in range(len(Sequences))]]
    
    
    
    # Get the Source Contour Image Sequences:
    SourceCISsequences = SourceRoi.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence
                               
    # Initialise the Source Contour Image Sequence to Source Dicom Mapping,
    # Target Contour Image Sequences, and Target Referenced SOP Instance UIDs:
    SourceCISToSourceDicomMapping = []
    TargetCISsequences = []
    TargetCISRefSopUids = []
                
    # Cycle through each Source Contour Image Sequence:
    for i in range(len(SourceCISsequences)):
        #print(f'\nContour Image Sequence #{i}:\n')
        
        # This Source Contour Image Sequence:
        SourceCISsequence = SourceCISsequences[i]
              
        # Get the Referenced SOP Instance UID for this SourceSequence:
        RefSopUid = SourceCISsequence.ReferencedSOPInstanceUID
        
        # Store this RefSopUid:
        SourceCISRefSopUids.append(RefSopUid)
    
        # Find the index of the matching SourceSopUid:
        ind = SourceSopUids.index(RefSopUid)
        
        """ April 6:
        I'll probably need to add some logic above to deal with the case 
        where there isn't a match (for example, if the original SourceDicoms
        had greater length than the original TargetDicoms) """
        
        # Store this ind:
        SourceCISToSourceDicomMapping.append(ind)
        
        """ 
        Note:
            Each element in SourceCISToSourceDicomMapping relates a Source
            Contour Image Sequence to the Source DICOM it belongs to.
            
            This DICOM may not be included in SourceToTargetMapping, if, for
            example, there were many more DICOMs in the original SourceDicoms
            than the number of DICOMs in the original TargetDicoms.
        """
        # If this ind is SourceToTargetMapping, get the SOP UID of this 
        # TargetDicom, and copy this SourceSequence to TargetSequence:
        if ind in SourceToTargetMapping:
            # Copy this sequence:
            TargetCISsequence = copy.deepcopy(SourceCISsequence)
            
            # Modify the Referenced SOP Instance UID for this sequence, using
            # the TargetSopUid (the SOP Instance that best matches this 
            # SourceDicom):
            TargetCISsequence.ReferencedSOPInstanceUID = TargetSopUids[ind]
            
            # Store this TargetCISRefSopUid:
            TargetCISRefSopUids.append(TargetSopUids[ind])
            
            # Append this TargetSequence:
            TargetCISsequences.append(TargetCISsequence)
            
    # Modify the Contour Image Sequences of TargetRoi:
    TargetRoi.ReferencedFrameOfReferenceSequence[0]\
             .RTReferencedStudySequence[0]\
             .RTReferencedSeriesSequence[0]\
             .ContourImageSequence = copy.deepcopy(TargetCISsequences)
    
    """ ***************************************************************** """
    
    
    
    """
    ***********************************************************************
    END OF STEP 2a.
    ***********************************************************************                
    """
    
              
    
    """
    ***********************************************************************
    STEP 2b:
        Modify the Referenced SOP Instance UID in each Contour Sequence.
    ***********************************************************************
    
    For each Referenced SOP Instance UID in each Contour Sequence in 
    SourceRoi get the index of the SourceDicom with matching SOP Instance UID.
    
    """
       
    # Get the Source Contour Sequences:
    SourceCSsequences = SourceRoi.ROIContourSequence[0].ContourSequence
                               
    # Initialise the Source Contour Sequence to Source Dicom Mapping,
    # Target Contour Sequences, and Target Referenced SOP Instance UIDs:
    SourceCSToSourceDicomMapping = []
    TargetCSsequences = []
    TargetCSRefSopUids = []
                
    # Cycle through each Source Contour Sequence:
    for i in range(len(SourceCSsequences)):
        #print(f'\nContour Sequence #{i}:\n')
        
        # This Source Contour Sequence:
        SourceCSsequence = SourceCSsequences[i]
              
        # Get the Referenced SOP Instance UID for this SourceSequence:
        RefSopUid = SourceCSsequence.ContourImageSequence[0]\
                                    .ReferencedSOPInstanceUID
        
        # Store this RefSopUid:
        SourceCISRefSopUids.append(RefSopUid)
    
        # Find the index of the matching SourceSopUid:
        ind = SourceSopUids.index(RefSopUid)
        
        """ April 6:
        I'll probably need to add some logic above to deal with the case 
        where there isn't a match (for example, if the original SourceDicoms
        had greater length than the original TargetDicoms) """
        
        # Store this ind:
        SourceCSToSourceDicomMapping.append(ind)
        
        """ 
        Note:
            Each element in SourceCSToSourceDicomMapping relates a Source
            Contour Sequence to the Source DICOM it belongs to.
            
            This DICOM may not be included in SourceToTargetMapping, if, for
            example, there were many more DICOMs in the original SourceDicoms
            than the number of DICOMs in the original TargetDicoms.
        """
        # If this ind is SourceToTargetMapping, get the SOP UID of this 
        # TargetDicom, and copy this SourceSequence to TargetSequence:
        if ind in SourceToTargetMapping:
            # Copy this sequence:
            TargetCSsequence = copy.deepcopy(SourceCSsequence)
            
            # Modify the Referenced SOP Instance UID for this sequence, using
            # the TargetSopUid (the SOP Instance that best matches this 
            # SourceDicom):
            TargetCSsequence.ContourImageSequence[0]\
                            .ReferencedSOPInstanceUID = TargetSopUids[ind]
            
            # Store this TargetCSRefSopUid:
            TargetCSRefSopUids.append(TargetSopUids[ind])
            
            # Append this TargetSequence:
            TargetCSsequences.append(TargetCSsequence)
            
    # Modify the Contour Sequences of TargetRoi:
    TargetRoi.ROIContourSequence[0]\
             .ContourSequence = copy.deepcopy(TargetCSsequences)
    
    """ ***************************************************************** """
    
    
    
    
    
    
    
    
    """
    ***********************************************************************
    END OF STEP 2b.
    ***********************************************************************                
    """
    
    
    
    
    """
    ***********************************************************************
    STEP 3:
        Make non-essential changes that cannot be implemented in an iterative
        fashion but are still useful.
    ***********************************************************************
    
    
    ***********************************************************************
    STEP 3a:
        Modify the Structure Set Label.
    ***********************************************************************
    """
        
    # Add "_copy_from_Scan_A_to_B_YYYYMMDD_HHMMSS" to the original Structure 
    # Set Label, where A and B are the Scan Labels of the Source and Target:
    NewDate = time.strftime("%Y%m%d", time.gmtime())
    NewTime = time.strftime("%H%M%S", time.gmtime())
    
    #NewLabel = roi1.StructureSetLabel + '_copy_' + NewDate \
    #+ '_' + NewTime
    
    if False:
        NewLabel = SourceRoi.StructureSetLabel + '_copy_from_SeriesNo_' \
                   + SourceSeriesNo + '_to_' + TargetSeriesNo + '_' \
                   + NewDate + '_' + NewTime
        
        if Shift: 
            #NewLabel = NewLabel + '_' + SearchBy + 'PosCorr'
            NewLabel = NewLabel + '_Shifted'
        
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
            NewLabel = NewLabel + '_' + RegMethod + '-Reg'
    
    
    NewLabel = SourceRoi.StructureSetLabel + '(Series_no_' + SourceSeriesNo \
                + ')_copied_for_' + TargetKey + 'Dicoms_(Series_no_' \
                + TargetSeriesNo + ')_' + NewDate + '_' + NewTime
               
    if Debug:
        print('\n--- Changing the Structure Set Label from:\n', \
          TargetRoi.StructureSetLabel, '\nto:\n', NewLabel)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'OldStructureSetLabel':TargetRoi.StructureSetLabel})
    
    # Modify the Structure Set Label:
    TargetRoi.StructureSetLabel = NewLabel  
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'NewStructureSetLabel':TargetRoi.StructureSetLabel})
    
    """
    ***********************************************************************
    END OF STEP 3a.
    ***********************************************************************
    """
    

    """
    ***********************************************************************
    STEP 3b:
        Modify the ROI Name in the Structure Set ROI Sequence
    ***********************************************************************
    """
    
    # Modify the ROI Name by making it the same as the new Structure Set Label:
    
    if Debug:
        print('\n--- Changing the ROI Name from:\n', \
              TargetRoi.StructureSetROISequence[0].ROIName, '\nto:\n', NewLabel)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'OldROIName':\
                        TargetRoi.StructureSetROISequence[0].ROIName})
    
    # Modify the ROI Name:
    TargetRoi.StructureSetROISequence[0].ROIName = NewLabel
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'NewROIName':\
                        TargetRoi.StructureSetROISequence[0].ROIName})
                    
    """
    ***********************************************************************
    END OF STEP 3b.
    ***********************************************************************
    """
    
    
    """
    ***********************************************************************
    STEP 3c:
        Modify the Structure Set Date and Time to help identify the new ROI
        Collection from others.
    ***********************************************************************
    """
    
    # Only make changes if the Structure Set Date != NewDate:
    if TargetRoi.StructureSetDate==NewDate:
        if Debug:
            print('\n--- The Structure Set Date is today:\n', NewDate)
    else:
        if Debug:
            print('\n--- Changing the Structure Set Date from:\n', \
                  TargetRoi.StructureSetDate, '\nto:\n', NewDate)
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'OldStructureSetDate':TargetRoi.StructureSetDate})
    
        # Modify the Structure Set Date:
        TargetRoi.StructureSetDate = NewDate
        
        # Update the dictionary ChangesDict:
        ChangesDict.update({'OewStructureSetDate':TargetRoi.StructureSetDate})
    
    # - Update the Structure Set Time:
    #NewTime = time.strftime("%H%M%S", time.gmtime())
    
    if Debug:
        print('\n--- Changing the Structure Set Time from:\n', \
              TargetRoi.StructureSetTime, '\nto:\n', NewTime)
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'OldStructureSetTime':TargetRoi.StructureSetTime})
    
    # Modify the Structure Set Time:
    TargetRoi.StructureSetTime = NewTime
    
    # Update the dictionary ChangesDict:
    ChangesDict.update({'NewStructureSetTime':TargetRoi.StructureSetTime})
    
    """
    ***********************************************************************
    END OF STEP 3c.
    ***********************************************************************
    """
    
    
    """
    ***********************************************************************
    STEP 3d:
        Create a new filename for the DICOM-RTSTRUCT file.
    ***********************************************************************
    """
      
    # Start by getting the filename of the CopyFrom DICOM-RTSTRUCT file:
    SourceRoiFname = os.path.basename(SourceRoiFpath)
    
    # Get filename (with extension):
    SourceRoiFname = os.path.splitext(SourceRoiFname)
    
    # Get filename without extension:
    SourceRoiFname = SourceRoiFname[0]
    
    if False:
        """
        The file name format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
        "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
        or similar, find the first "_" character and use the characters that 
        precede it in the new file name.
        """
        # Find the index of the '_' character in the file name:
        #ind = origRoiFname.index('_')
        #ind = SourceRoiFname[0].index('_')
        ind = SourceRoiFname.index('_')
        
        """
        TargetRoiFname will be the same string of characters in SourceRoiFname
        including the '_', followed by the recently modified Structure Set Date
        and Structure Set Time, followed by some characters that indicate if the
        positions were "corrected" (e.g. 'zPosCorr'):
        """
        
        #TargetRoiFname = SourceRoiFname[0][0:ind+1] + TargetRoi.StructureSetDate \
        #                 + '_' + TargetRoi.StructureSetTime
        #TargetRoiFname = SourceRoiFname[0][0:ind+1]
        TargetRoiFname = SourceRoiFname[0:ind+1]
        
        """
        Note: The first index [0] gets the filename without extension. The second
        index [0:ind+1] gets all characters to and including '_'. 
        """
        
        # If Shift is not equal to zero, add the character in SearchBy 
        # (e.g. 'z') followed by 'PosCorr':
        if Shift: 
            TargetRoiFname = TargetRoiFname + '_' + SearchBy + 'PosCorr'
        
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
            TargetRoiFname = TargetRoiFname + '_' + RegMethod + '-Reg'
        
        # Add the ROI Structure Set Date and Time:
        TargetRoiFname = TargetRoiFname + '_' + TargetRoi.StructureSetDate \
                         + '_' + TargetRoi.StructureSetTime
    
    
    # Define the target ROI filename:             
    TargetRoiFname = SourceRoiFname + '_copied_for_' + TargetKey \
               + 'Dicoms_(Series_no_' + TargetSeriesNo + ')_' \
               + NewDate + '_' + NewTime + ')'
    
        
    # Finally add the file extension:
    TargetRoiFname = TargetRoiFname + '.dcm'
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    TargetRoiFpath = os.path.join(TargetRoiDir, TargetRoiFname)
    
    if Debug:
        print('\nSourceRoiFname =', SourceRoiFname)
        print('TargetRoi.StructureSetDate =', TargetRoi.StructureSetDate)
        print('TargetRoi.StructureSetTime =', TargetRoi.StructureSetTime)
        
        
    
    """
    ***********************************************************************
    END OF STEP 3d.
    ***********************************************************************
    """
    
    
    """
    ***********************************************************************
    STEP 4:
        Create the TargetRoiDir, the JSON filepath and update DataDict.
    ***********************************************************************
    
    
    Note:
        The dictionary DataDict will be updated (below) and exported to a JSON 
        file in the path defined by DataDict[TargetKey]['RoiDir'].
    """
    # Create a filename for the JSON file, using the same base filename as 
    # the RTSTRUCT file:
    JsonFname = os.path.splitext(TargetRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    JsonFpath = os.path.join(TargetRoiDir, JsonFname)
    
    # Create the directory TargetRoiDir:
    CreateDir(TargetRoiDir, Debug)
    
    #print('\nChangesDict:\n', ChangesDict)
    
    #print('\nTargetKey =', TargetKey)
    
    #print('\nDataDict[' + TargetKey + ']:\n', DataDict[TargetKey])
    
    # Update the dictionary DataDict:
    DataDict[TargetKey].update({'ChangesToTargetRoi':ChangesDict})
    #DataDict[TargetKey].update({'newRoiFname':newRoiFname})
    DataDict[TargetKey].update({'RoiFpath':TargetRoiFpath})
    #DataDict[TargetKey].update({'jsonFname':jsonFname})
    DataDict[TargetKey].update({'JsonFpath':JsonFpath})
    
    """
    ***********************************************************************
    END OF STEP 4.
    ***********************************************************************
    """
    
    
    
    """
    ***********************************************************************
    STEP 5:
        Print some useful info.
    ***********************************************************************
    
    
    March 28 Note:
        The code below will need to be checked and corrected.
        
    April 6:
        The numbers don't seem to agree so ignoring this section for now.
    """ 
   
    if False:
        # Cycle through each Reference SOP Instance UID for each Contour Image
        # Sequence:
        for i in range(len(OldTargetCISRefSopUids)):
            #OldTargetRefSopUid = OldTargetCISRefSopUids[i]
            
            # Get the index of the this RefSopUid:
            #ToInd = ToRefSopUids.index(RefSopUid)
            
            # Find the index of the DICOM in SourceDicoms whose SOP Instance UID
            # matches this refsopuid:
            #FromInd = SourceSopUids.index(RefSopUid)
            
            # Use SourceCISRefSopUidToSourceSopUidMapping to get the Source DICOM
            # index:
            SourceInd = SourceCISRefSopUidtoSourceSopUidMapping[i]
            
            # Use SourceToTargetMapping to get the Target DICOM index:
            TargetInd = SourceToTargetMapping[SourceInd]
            
            # Use SourceInd and TargetInd to get the slice numbers:
            SourceSliceNo = SourceSliceNos[SourceInd]
            #TargetSliceNo = TargetSliceNos[MappingInds[i]]
            ##TargetSliceNo = TargetSliceNos[ToInd]
            TargetSliceNo = TargetSliceNos[TargetInd]
            
            if False:
                print('\nSourceSliceNo =', SourceSliceNo)
                #print(f'\MappingInds (N = {len(MappingInds)}) =', MappingInds)
                print('TargetSliceNo =', TargetSliceNo)
            
            # Get the Image Positions:
            SourceImPos = SourceDicoms[SourceInd].ImagePositionPatient 
            TargetImPos = TargetDicoms[TargetInd].ImagePositionPatient
            
            # Get the scanning dimension:
            if SearchBy=='x':
                dim = 0
            elif SearchBy=='y':
                dim = 1
            elif SearchBy=='z':
                dim = 2
            """
            Above code doesn't deal with case where SearchBy='r'!
            """
            
            # Get the position along the scanning dimension:
            SourceScanPos = SourceImPos[dim]
            TargetScanPos = TargetImPos[dim]
            
            # Get the positional differences along x, y, z and r: 
            #dPos = DataDict[TargetKey]['AbsMinDpos']
            dx = round(TargetImPos[0] - SourceImPos[0], 2)
            dy = round(TargetImPos[1] - SourceImPos[1], 2)
            dz = round(TargetImPos[2] - SourceImPos[2], 2)
            dr = round((dx**2 + dy**2 + dz**2)**.5, 2)
            
            # Get the positional scan difference:
            if SearchBy=='x':
                dpos = dx
            elif SearchBy=='y':
                dpos = dy
            elif SearchBy=='z':
                dpos = dz
            """
            Above code doesn't deal with case where SearchBy='r'!
            """
            
            if True:#Debug:
                print(f'\n\nContour Image Sequence #{i} of {len(OldTargetCISRefSopUids)}:\n', \
                      f'\n   Copied from slice #{SourceSliceNo} in series' \
                      f' #{SourceSeriesNo} to slice #{TargetSliceNo} in series' \
                      f' #{TargetSeriesNo}:\n', \
                      f'\n   - from Image Position', SourceImPos, \
                      f'\n                      to', TargetImPos, '\n', \
                      #f'\n   - from Scan Position', SourceScanPos, 'mm to', \
                      #TargetScanPos, 'mm\n', \
                      f'\n   - with positional difference {dpos} mm.')
        
    """
    ***********************************************************************
    END OF STEP 5.
    ***********************************************************************
    """
    
    
    """
    ***********************************************************************
    STEP 6:
        Export the new RTSTRUCT and JSON files.
    ***********************************************************************
    """ 
    
    # Export the dictionary DataDict to a JSON file:
    json.dump(DataDict, open(JsonFpath, 'w'))
    
    print('\n--- JSON file created:\n', JsonFpath, '\n\n')
    
    
    # Save the new DICOM-RTSTRUCT file:
    TargetRoi.save_as(TargetRoiFpath)
    
    """
    ***********************************************************************
    END OF STEP 6.
    ***********************************************************************
    """
    
    return SourceRoi, TargetRoi, DataDict 
    
