# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:43:08 2020

@author: ctorti
"""



"""
April 2:  
    Consider setting to zero all pixels in the pixel arrays of the DICOMs 
    that have large positional differences (say > 0.5 mm or > 1 mm).
"""



def RemapDicoms(PreRemappedDicoms, DataDict):
    
    # Import packages and functions:
    
    import os
    #import copy
    import pydicom
    #import time
    #import numpy as np
    #import importlib
    from CreateDir import CreateDir
    
    
    # Get some necessary info from DataDict:
    Debug = DataDict['Debug']
    #MappingInds = DataDict['Mapping']['MappingInds']
    
    # The key of the series to be remapped:
    ToRemapKey = DataDict['ToRemapKey']
    
    # The filepaths and slice numbers of the DICOMs to be remapped: 
    PreRemappedFpaths = DataDict[ToRemapKey]['DicomFpaths']
    PreRemappedSliceNos = DataDict[ToRemapKey]['SliceNos']
    PreRemappedSeriesDesc = DataDict[ToRemapKey]['SeriesDesc']
    
    # Get the mapping indeces:
    if 'Source' in ToRemapKey:
        MappingInds = DataDict['TargetToSourceMapping']['MappingInds']
    else:
        MappingInds = DataDict['SourceToTargetMapping']['MappingInds']
    
    # Use MappingInds to get a list of remapped PreRemappedFpaths:
    RemappedPreRemappedFpaths = [PreRemappedFpaths[ind] for ind in MappingInds]
    
    # Use MappingInds to get a list of remapped DICOMS and SliceNos:
    RemappedDicoms = [PreRemappedDicoms[ind] for ind in MappingInds]
    RemappedSliceNos = [PreRemappedSliceNos[ind] for ind in MappingInds]
    
    print('\nPreRemappedSliceNos =', PreRemappedSliceNos)
    print('\nlen(PreRemappedSliceNos) =', len(PreRemappedSliceNos))
    print('\nRemappedSliceNos =', RemappedSliceNos)
    print('\nlen(RemappedSliceNos) =', len(RemappedSliceNos))
    
    
    # The key of the remapped series:
    RemappedKey = DataDict['RemappedKey']
    
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
    
    
    # Create a dictionary to store the changes made to the DICOM of the  
    # registered slices:
    #changesDict = {}
    
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
    
    
    return RemappedDicoms, DataDict 
