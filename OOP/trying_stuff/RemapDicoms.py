# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:44:12 2020

@author: ctorti
"""





def RemapDicoms(PreRemappedDicoms, DataDict, KeyToRemap):
    
    # Import packages and functions:
    
    import os
    import copy
    import pydicom
    #import time
    #import numpy as np
    #import importlib
    from CreateDir import CreateDir
    #from GetScanPositions import GetScanPositions
    
    
    # Get some necessary info from DataDict:
    #SearchBy = DataDict['SearchBy']
    Debug = DataDict['Debug']
    #KeyToRemap = DataDict['KeyToRemap']
    
    # The filepaths and slice numbers of the DICOMs to be remapped: 
    PreRemappedFpaths = DataDict[KeyToRemap]['DicomFpaths']
    PreRemappedSliceNos = DataDict[KeyToRemap]['SliceNos']
    PreRemappedSeriesDesc = DataDict[KeyToRemap]['SeriesDesc']
    PreRemappedScanPos = DataDict[KeyToRemap]['ScanPos']

    
    """ Note April 27:
       Previously RemapDicoms() was only used to remap the longer series, i.e.
       if there were more slices in SourceDicoms than in TargetDicoms,
       KeyToRemap = 'Source' and the indeces in 
       DataDict['SourceToTarget']['MappingInds'] was used.
       
       Problems occured with the registration however when the two slices from
       Source and Target were vastly different.  So the key 'SourceToTarget'
       was replaced with 'SourceToTargetNearest' (which was equivalent to the
       previous approach) and 'SourceToTargetNearOnly'. 
       
       The latter is a subset of the former where all indeces with scan
       positional differences greater than a threshold distance (e.g. 1 mm)
       are omitted.
       
       Hence, rather than only remapping the longer of SourceDicoms and 
       TargetDicoms, now the shorter series also needs to be remapped to 
       account for the indeces that were omitted.
       
       GetClosestSlices() not only outputs 'SourceToTargetNearOnly' and
       'TargetToSourceNearOnly', but also 'SourceNearOnly' and 'TargetNearOnly'.
       The latter two contain the indeces that the shorter series need to be
       remapped to.
       
       Defining the following:
           
       S2Tinds = DataDict['SourceToTargetNearOnly']['MappingInds']
       T2Sinds = DataDict['TargetToSourceNearOnly']['MappingInds']
       
       Sinds = DataDict['SourceNearOnly']['MappingInds']
       Tinds = DataDict['TargetNearOnly']['MappingInds']
        
       S2T = len(S2Tinds)
       T2S = len(T2Sinds)
        
       If KeyToRemap = 'Source' and S2T > T2S, then 
       SourceDicoms needs to be remapped using T2Sinds, and 
       TargetDicoms needs to be remapped using Tinds.
       
       If KeyToRemap = 'Source' and S2T < T2S, then
       SourceDicoms needs to be remapped using Sinds, and
       TargetDicoms needs to be remapped using S2Tinds.
       
       If KeyToRemap = 'Target' and S2T > T2S, then
       SourceDicoms needs to be remapped using T2Sinds, and 
       TargetDicoms needs to be remapped using Tinds.
       
       If KeyToRemap = 'Target' and S2T < T2S, then
       SourceDicoms needs to be remapped using Sinds, and
       TargetDicoms needs to be remapped using S2Tinds.
       
       Each operation shown under the If statement will be done by separate
       calls to this function RemapDicoms.  So the following are the list of
       possible inputs and outputs:
       
       If KeyToRemap = 'Source' and S2T > T2S, 
       Then SourceDicoms needs to be remapped using T2Sinds.
       
       If KeyToRemap = 'Source' and S2T < T2S,
       Then SourceDicoms needs to be remapped using Sinds.
       
       If KeyToRemap = 'Target' and S2T > T2S,
       Then TargetDicoms needs to be remapped using Tinds.
       
       If KeyToRemap = 'Target' and S2T < T2S, 
       Then TargetDicoms needs to be remapped using S2Tinds.
    
    
    """
    
    # Get the mapping indeces (use the "NearOnly" indeces that correspond to
    # only slices that are within the maximum allowed scan position difference
    # defined in GetClosestSlices()):
    # S2Tinds = SourceToTargetInds and T2Sinds = TargetToSourceInds:
    S2Tinds = DataDict['SourceToTargetNearOnly']['MappingInds']
    T2Sinds = DataDict['TargetToSourceNearOnly']['MappingInds']
    
    Sinds = DataDict['SourceNearOnly']['MappingInds']
    Tinds = DataDict['TargetNearOnly']['MappingInds']
    
    S2T = len(S2Tinds)
    T2S = len(T2Sinds)
    
    # The scan positional differences are:
    S2Tdiffs = DataDict['SourceToTargetNearOnly']['ScanPosDiffs']
    T2Sdiffs = DataDict['TargetToSourceNearOnly']['ScanPosDiffs']
    
    Sdiffs = DataDict['SourceNearOnly']['ScanPosDiffs']
    Tdiffs = DataDict['TargetNearOnly']['ScanPosDiffs']
    
    
    # If SourceDicoms are to be remapped...
    if 'Source' in KeyToRemap:
        if S2T > T2S:
            # ... and S2T > T2S:
            MappingInds = copy.deepcopy(T2Sinds)
            RemappedScanPosDiffs = copy.deepcopy(T2Sdiffs)
        
        else:
            # ... and S2T < T2S:
            MappingInds = copy.deepcopy(Sinds)
            RemappedScanPosDiffs = copy.deepcopy(Sdiffs)
        
        
    # If TargetDicoms are to be remapped...
    if 'Target' in KeyToRemap:
        if S2T > T2S:
            # ... and S2T < T2S:
            MappingInds = copy.deepcopy(Tinds)
            RemappedScanPosDiffs = copy.deepcopy(Tdiffs)
            
        else:
            # ... and S2T < T2S:
            MappingInds = copy.deepcopy(S2Tinds)
            RemappedScanPosDiffs = copy.deepcopy(S2Tdiffs)
        
    
    
    # Use MappingInds to get a list of remapped filepaths in PreRemappedFpaths:
    RemappedPreRemappedFpaths = [PreRemappedFpaths[ind] for ind in MappingInds]
    
    # Use MappingInds to get a list of remapped DICOMS and SliceNos:
    RemappedDicoms = [PreRemappedDicoms[ind] for ind in MappingInds]
    RemappedSliceNos = [PreRemappedSliceNos[ind] for ind in MappingInds]
    RemappedScanPos = [PreRemappedScanPos[ind] for ind in MappingInds]
    
    # Get the scan positions of the remapped DICOMs:
    #RemappedScanPos = GetScanPositions(RemappedDicoms, SearchBy)
    
    
    if Debug:
        print('\nPreRemappedSliceNos =', PreRemappedSliceNos)
        print('\nlen(PreRemappedSliceNos) =', len(PreRemappedSliceNos))
        print('\nRemappedSliceNos =', RemappedSliceNos)
        print('\nlen(RemappedSliceNos) =', len(RemappedSliceNos))
        print('\nRemappedScanPos =', RemappedScanPos)
        print('\nlen(RemappedScanPos) =', len(RemappedScanPos))
        print('\nRemappedScanPosDiffs =', RemappedScanPosDiffs)
        print('\nlen(RemappedScanPosDiffs) =', len(RemappedScanPosDiffs))
    
    
    # The key of the remapped series:
    #RemappedKey = DataDict['RemappedKey']
    RemappedKey = 'Remapped' + KeyToRemap
    
    print('\nKeyToRemap =', KeyToRemap)
    print('RemappedKey =', RemappedKey)
    
    # And the Series No for the remapped series, and the directory where the 
    # remapped DICOMs are to be saved:
    RemappedDicomDir = DataDict[RemappedKey]['DicomDir']
    RemappedSeriesNo = DataDict[RemappedKey]['SeriesNo']
    RemappedSeriesDesc = DataDict[RemappedKey]['SeriesDesc']
    
    
    
    """
    The Study Instance UID and Frame of Reference UID of RegDicoms will be kept 
    unchanged from those of MovingDicoms, but a new Series Instance UID will 
    need to be generated.  
    
    Unique SOP Instance UIDs will be generated for each RegDicom.  
    
    I'm not sure if the Frame of Reference UID of RegDicoms should be changed
    from that of MovingDicoms.  For now I'll assume not.
    """
    
    # Start by using the first DICOM:
    Dicom = PreRemappedDicoms[0]
    
    # Get the Series Description:
    #OrigSeriesDesc = Dicom.SeriesDescription
    
    # Create a new Series Description by appending 'Remapped' to the original
    # Series Description:
    #NewSeriesDesc = OrigSeriesDesc + ' - Remapped'
    RemappedSeriesDesc = 'Remapped ' + PreRemappedSeriesDesc
    
    """ missing stuff here?... """
    
    # Create the directory:
    CreateDir(RemappedDicomDir, Debug)
    
    # Create a new Series Instance UID:
    RemappedSeriesUid = pydicom.uid.generate_uid()
    
    # Initialise a list to store the filepaths and SOP Instance UIDs of the 
    # exported remapped DICOMs:
    RemappedFpaths = []
    RemappedSopUids = []
    
    F = len(RemappedDicoms)
    
    # Iterate for each DICOM in RemappedDicoms:
    for f in range(F):
        # Get the remapped DICOM:
        Dicom = RemappedDicoms[f]
        
        # Get the remapped original filename with extension:
        OrigFname = os.path.basename(RemappedPreRemappedFpaths[f])
        
        # Get filename with extension:
        OrigFname = os.path.splitext(OrigFname)
        
        # Get filename without extension:
        OrigFnameNoExt = OrigFname[0]
        
        # Create a new filename, adding ' - Remapped f+1' to it, where f+1 is
        # an incremented integer:
        RemappedFname = OrigFnameNoExt + f'-Remapped-{f+1}.dcm'
        
        if Debug:
            print('\n--- Changing the file name from:\n', \
              OrigFname, '\nto:\n', RemappedFname)
        
        # Create new filepath:
        RemappedFpath = os.path.join(RemappedDicomDir, RemappedFname)
        
        # Append RemappedFpath to RemappedFpaths:
        RemappedFpaths.append(RemappedFpath)
        
        
        """ Amend the DICOM tags: """
        
        # Get the original SOP Instance UID:
        PreRemappedSopUid = Dicom.SOPInstanceUID

        # Generate a new SOP Instance UID:
        RemappedSopUid = pydicom.uid.generate_uid() 
        
        # Store the new SOP Instance UID:
        RemappedSopUids.append(RemappedSopUid)
        
        if Debug:
            print('\n--- Changing the SOP Instance UID from:\n', \
              PreRemappedSopUid, '\nto:\n', RemappedSopUid)
        

    
        # Make changes to the tags:
        Dicom.SOPInstanceUID = RemappedSopUid
        Dicom.SeriesInstanceUID = RemappedSeriesUid
        #Dicom.SeriesDescription = NewSeriesDesc
        Dicom.SeriesDescription = RemappedSeriesDesc
        Dicom.SeriesNumber = RemappedSeriesNo
        #RegDicom.StudyInstanceUID = postRegStudyUID
        #RegDicom.SeriesInstanceUID = postRegSeriesUID
        #RegDicom.FrameOfReferenceUID = postRegFrameOfRefUID
        
        
        if Debug:
            print('\n--- Exporting the registered DICOM to:\n', \
                  RemappedFpath)
        
        # Export the DICOM file:
        Dicom.save_as(RemappedFpath)
        
        
        
    # Update the dictionary DataDict:
    DataDict[RemappedKey].update({'DicomFpaths':RemappedFpaths})
    DataDict[RemappedKey].update({'SopUids':RemappedSopUids})
    DataDict[RemappedKey].update({'SliceNos':RemappedSliceNos})
    DataDict[RemappedKey].update({'ScanPos':RemappedScanPos})
    DataDict[RemappedKey].update({'ScanPosDiffs':RemappedScanPosDiffs})
    
    
    
    return RemappedDicoms, DataDict 
