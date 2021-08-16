# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 13:44:00 2020

@author: ctorti
"""





def ShiftDicoms(PreShifteddDicoms, DataDict):
    
    # Import packages:
    import copy
    import pydicom
    import os
    from CreateDir import CreateDir
    
    
    # Get some necessary info from DataDict:
    SearchBy = DataDict['SearchBy']
    Shift = DataDict['Shift']
    Debug = DataDict['Debug']
    #OldFpaths = DataDict['Target']['DicomFpaths'] # assume shift will always be
    # applied to PreShifteddDicoms in 'To' key
    #ShiftedDicomDir = DataDict['Shifted']['DicomDir']
    #ShiftedDicomDir = DataDict['To']['DicomDir'] # this will be modified later
    #""" actually ShiftedDicomDir is not modified later """
    #ShiftedSeriesNo = DataDict['Shifted']['SeriesNo']
    
    # The key of the series that was shifted:
    ToShiftKey = DataDict['ToShiftKey']
    
    # The directory of the series to be shifted:
    #SourceDicomDir = DataDict[ShiftKey]['DicomDir']
    
    # The filepaths of the DICOMs to be shifted:
    PreShiftedFpaths = DataDict[ToShiftKey]['DicomFpaths']
    #PreShiftedSeriesNo = DataDict[ToShiftKey]['SeriesNo']
    PreShiftedSeriesDesc = DataDict[ToShiftKey]['SeriesDesc']
    
    # The key for the shifted series:
    ShiftedKey = DataDict['ShiftedKey']
    
    # The DICOM directory and Series no for the shifted series:
    ShiftedDicomDir = DataDict[ShiftedKey]['DicomDir']
    ShiftedSeriesNo = DataDict[ShiftedKey]['SeriesNo']
    
    
    # Get the Image Position (Patient) for all PreShifteddDicoms:
    PreShiftedPositions = [PreShifteddDicoms[i].ImagePositionPatient for i in range(len(PreShifteddDicoms))]
    
    # Copy the old positions:
    #ShiftedPositions = copy.deepcopy(PreShiftedPositions)
        
    # Apply shift to ShiftedPositions:
    """
    Note this doesn't currently handle the case SearchBy = 'r'...
    """
    if SearchBy=='x':
        ShiftedPositions = [[PreShiftedPositions[i][0] + Shift, \
                             PreShiftedPositions[i][1], \
                             PreShiftedPositions[i][2]
                             ] for i in range(len(PreShiftedPositions))]
    
    if SearchBy=='y':
        ShiftedPositions = [[PreShiftedPositions[i][0], \
                             PreShiftedPositions[i][1] + Shift, \
                             PreShiftedPositions[i][2]
                             ] for i in range(len(PreShiftedPositions))]

    if SearchBy=='z':
        ShiftedPositions = [[PreShiftedPositions[i][0], \
                             PreShiftedPositions[i][1], \
                             PreShiftedPositions[i][2] + Shift
                             ] for i in range(len(PreShiftedPositions))]
    
    
    if Debug:
        print('\nPosition changed from:\n', PreShiftedPositions, '\nto:\n',\
              ShiftedPositions)
            
        
    # Create new list of DICOM objects by copying the original list:
    ShiftedDicoms = copy.deepcopy(PreShifteddDicoms)
        
    """
    The Study Instance UID and Frame of Reference UID of the scan-position-
    shifted DICOM will be kept unchanged, but a new Series Instance UID 
    will need to be generated.  
    
    Unique SOP Instance UIDs will be generated for each individual DICOM.  
    I'm not sure if the Frame of Reference UID should also be changed.
    For now I'll assume not.
    """
    
    # Create a new Series Instance UID:
    ShiftedSeriesUid = pydicom.uid.generate_uid()
    
    # Get the Series Description of the first DICOM:
    #OldSeriesDesc = PreShifteddDicoms[0].SeriesDescription
    
    # Create a new Series Description by appending info about the shift that
    # was applied to the original Series Description:
    #NewSeriesDesc = OldSeriesDesc + ' - ' + SearchBy + '-Shifted'
    ShiftedSeriesDesc = 'Shifted ' + PreShiftedSeriesDesc
    
    # Modify ShiftedDicomDir to include the Series Description in 
    # parentheses:
    #ShiftedDicomDir = ShiftedDicomDir + ' (' + NewSeriesDesc + ')'
    
    # Create the directory:
    CreateDir(ShiftedDicomDir, Debug)
    
    
    # Initialise list of ShiftedFpaths and ShiftedSopUids:
    ShiftedFpaths = []
    ShiftedSopUids = []
    
    # Loop for each DICOM in ShiftedDicoms:
    for i in range(len(ShiftedDicoms)):
        # Modify the original positions with the new positions:
        ShiftedDicoms[i].ImagePositionPatient = ShiftedPositions[i]
        
        # Generate a new SOP Instance UID:
        ShiftedSopUid = pydicom.uid.generate_uid()
        
        # Append to ShiftedSopUids:
        ShiftedSopUids.append(ShiftedSopUid)
        
        # Modify the old SOP UID with the new one:
        ShiftedDicoms[i].SOPInstanceUID = ShiftedSopUid
        
        # Modify the old Series UID with the new one:
        ShiftedDicoms[i].SeriesInstanceUID = ShiftedSeriesUid
        
        # Modify the Series Number:
        ShiftedDicoms[i].SeriesNumber = ShiftedSeriesNo
        
        # Modify the Series Description:
        #ShiftedDicoms[i].SeriesDescription = NewSeriesDesc
        ShiftedDicoms[i].SeriesDescription = ShiftedSeriesDesc

        # Get the original file name with extension from toFpaths:
        PreShiftedFname = os.path.basename(PreShiftedFpaths[i]) 
        
        # Get the file name without extension:
        PreShiftedFnameNoExt = os.path.splitext(PreShiftedFname)[0]
        
        # Create a new filename, adding some info about the scan position 
        # shift:
        ShiftedFname = PreShiftedFnameNoExt + '-' + SearchBy + '-Shifted.dcm'
        
        if Debug:
            print('\n   Changing the file name from:\n', \
                  PreShiftedFname, '\nto:\n', ShiftedFname)
        
        # Create new filepath:
        ShiftedFpath = os.path.join(ShiftedDicomDir, ShiftedFname)
        
        # Append to ShiftedFpaths:
        ShiftedFpaths.append(ShiftedFpath)
        
        
        if Debug:
            print('\n    Exporting the position corrected DICOM to:\n', \
                  ShiftedFpath)
        
        # Export the DICOM file:
        ShiftedDicoms[i].save_as(ShiftedFpath)
        
        
        
    # Update DataDict:
    DataDict[ShiftedKey]['DicomDir'] = ShiftedDicomDir
    DataDict[ShiftedKey]['DicomFpaths'] = ShiftedFpaths
    DataDict[ShiftedKey]['SopUids'] = ShiftedSopUids
    #DataDict[ShiftedKey]['SeriesDesc'] = NewSeriesDesc
    DataDict[ShiftedKey]['SeriesDesc'] = ShiftedSeriesDesc
    DataDict[ShiftedKey]['SliceNos'] = DataDict[ToShiftKey]['SliceNos'] # make ='ToShiftKey'
    
    
    return ShiftedDicoms, DataDict