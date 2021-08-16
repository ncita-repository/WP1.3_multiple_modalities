# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 08:41:43 2020

@author: ctorti
"""


"""
Function:
    GetClosestSlices()
    
Purpose:
    Starting with two DICOM objects FromDicoms and ToDicoms, find the indeces
    that map each DICOM in FromDicoms to the DICOM in ToDicoms with nearest
    position along the x-/y-/z-direction or (x,y,z) position. 
    
Input:
    FromDicoms   - list of DICOM objects to be mapped to ToDicoms
    
    ToDicoms     - list of DICOM objects 
    
    DataDict     - dictionary containing several variables
                
    
Returns:
    ShiftedDicoms - list of all DICOM objects from ToDicoms with modified 
                    positions given by Shift (if Shift = 0, an
                    empty array will be returned for ShiftedDicoms)
                   
    MappedDicoms  - list of DICOM objects from ToDicoms (if Shift = 0) or 
                    ShiftedDicoms (if Shift != 0) whose scan positions 
                    most closely match that of the DICOMs in FromDicoms
                      
    DataDict      - dictionary containing several variables
"""


def GetClosestSlices(FromDicoms, ToDicoms, DataDict):
    
    # Import packages:
    import copy
    import pydicom
    import os
    from IndecesCount import IndecesCount
    from CreateDir import CreateDir
    
    # Get some necessary info from DataDict:
    SearchBy = DataDict['SearchBy']
    Shift = DataDict['Shift']
    Debug = DataDict['Debug']
    ToDicomFpaths = DataDict['To']['DicomFpaths']
    PosCorrDicomDir = DataDict['PosCorr']['DicomDir']
    #PosCorrScanLabel = DataDict['PosCorr']['ScanLabel']
    PosCorrSeriesNo = DataDict['PosCorr']['SeriesNo']
    if Shift:
        #MapPosCorrDicomDir = DataDict['MapPosCorr']['DicomDir']
        MappedDicomDir = DataDict['MapPosCorr']['DicomDir']
        #MapPosCorrScanLabel = DataDict['MapPosCorr']['ScanLabel']
        #MapPosCorrSeriesNo = DataDict['MapPosCorr']['SeriesNo']
        MappedSeriesNo = DataDict['MapPosCorr']['SeriesNo']
    else:
        MappedDicomDir = DataDict['Mapped']['DicomDir']
        #MappedScanLabel = DataDict['Mapped']['ScanLabel']
        MappedSeriesNo = DataDict['Mapped']['SeriesNo']
    
    
    
    
    # Get the Image Position (Patient) of each DICOM object in FromDicoms and
    # ToDicoms:
    #fromPositions = [FromDicoms[i].ImagePositionPatient for i in range(len(FromDicoms))]
    ToPositions = [ToDicoms[i].ImagePositionPatient for i in range(len(ToDicoms))]
    
    # Apply correction to scan direction if Shift is non-zero:
    if not Shift: # i.e. if Shift = 0
        ShiftedDicoms = []
        
    else: # i.e. if Shift != 0
        # Apply correction to scan position based on the chosen SearchBy:
        if SearchBy=='x':
            CorrPositions = [[ToPositions[i][0] + Shift, \
                              ToPositions[i][1], \
                              ToPositions[i][2]
                              ] for i in range(len(ToPositions))]
        
        if SearchBy=='y':
            CorrPositions = [[ToPositions[i][0], \
                              ToPositions[i][1] + Shift, \
                              ToPositions[i][2]
                              ] for i in range(len(ToPositions))]
    
        if SearchBy=='z':
            CorrPositions = [[ToPositions[i][0], \
                              ToPositions[i][1], \
                              ToPositions[i][2] + Shift
                              ] for i in range(len(ToPositions))]
        
        if Debug:
            print('\nPositions changed from:\n', ToPositions, '\nto:\n',\
                  CorrPositions)
                
            
        # Create new DICOM objects that are copies of ToDicoms but will have
        # modified Image Positions:
        ShiftedDicoms = copy.deepcopy(ToDicoms)
        
        """
        The Study Instance UID and Frame of Reference UID of the position-
        corrected DICOMs will be kept unchanged from those of ToDicoms, but a 
        new Series Instance UID will need to be generated.  
        Unique SOP Instance UIDs will be generated for each individual DICOM.  
        I'm not sure if the Frame of Reference UID of ShiftedDicoms should be 
        made equal to the FoR UID of the FromDicoms.  For now I'll assume not.
        """
        
        # Create a new Series Instance UID:
        NewSeriesUid = pydicom.uid.generate_uid()
        
        # Get the Series Description of the first DICOM in ToDicoms:
        ToSeriesDesc = ToDicoms[0].SeriesDescription
        
        # Create a new Series Description by appending ' - PosCorr' to the 
        # Series Description of ToDicoms:
        PosCorrSeriesDesc = ToSeriesDesc + ' - ' + SearchBy + 'PosCorr'
        
        # Modify PosCorrDicomDir to include the Series Description in 
        # parentheses:
        PosCorrDicomDir = PosCorrDicomDir + ' (' + PosCorrSeriesDesc + ')'
        
        # Modify PosCorrDicomDir in DataDict:
        DataDict['PosCorr']['DicomDir'] = PosCorrDicomDir
        
        # Create the directory:
        CreateDir(PosCorrDicomDir, Debug)
        
        # Initialise list of posCorrFpaths:
        PosCorrDicomFpaths = []
        
        # Loop for each DICOM in ShiftedDicoms:
        for i in range(len(ShiftedDicoms)):
            # Modify the original positions with the new positions:
            ShiftedDicoms[i].ImagePositionPatient = CorrPositions[i]
            
            # Generate a new SOP Instance UID:
            NewSopUid = pydicom.uid.generate_uid()
            
            # Modify the old SOP UIDs with the new ones:
            ShiftedDicoms[i].SOPInstanceUID = NewSopUid
            
            # Modify the old Series UIDs with the new one:
            ShiftedDicoms[i].SeriesInstanceUID = NewSeriesUid
            
            # Modify the Series Number:
            ShiftedDicoms[i].SeriesNumber = PosCorrSeriesNo
            
            # Modify the Series Description:
            ShiftedDicoms[i].SeriesDescription = PosCorrSeriesDesc
    
            # Get the original file name with extension from toFpaths:
            ToDicomFname = os.path.basename(ToDicomFpaths[i]) 
            
            # Get the file name without extension:
            ToDicomFnameNoExt = os.path.splitext(ToDicomFname)[0]
            
            # Create a new filename, adding some info about the scan position 
            # correction:
            PosCorrDicomFname = ToDicomFnameNoExt + '-' + SearchBy + 'PosCorr.dcm'
            
            if Debug:
                print('\n   Changing the file name from:\n', \
                      ToDicomFname, '\nto:\n', PosCorrDicomFname)
            
            # Create new filepath:
            PosCorrDicomFpath = os.path.join(PosCorrDicomDir, PosCorrDicomFname)
            
            # Append to posCorrFpath:
            PosCorrDicomFpaths.append(PosCorrDicomFpath)
            
            
            if Debug:
                print('\n    Exporting the position corrected DICOM to:\n', \
                      PosCorrDicomFpath)
            
            # Export the DICOM file:
            ShiftedDicoms[i].save_as(PosCorrDicomFpath)
            
            
        # Add the position-corrected DICOM fpaths to DataDict:
        DataDict['PosCorr'].update({'DicomFpaths':PosCorrDicomFpaths})
        DataDict['PosCorr'].update({'SeriesDesc':PosCorrSeriesDesc})
    
    
    if Debug:
        print(f'\nThere are {len(FromDicoms)} slices in FromDicoms.')
        print(f'There are {len(ToDicoms)} slices in ToDicoms.')
        if Shift:
            print(f'There are {len(ShiftedDicoms)} slices in ShiftedDicoms.')
            print(f'There are {len(PosCorrDicomFpaths)} fpaths in PosCorrDicomFpaths.')
            
    
    # Find the index of the nearest slice in the list of DICOMs ToDicoms or
    # MappedDicoms (depending on whether Shift is zero) to each slice in
    # FromDicoms.
    
    # Initialise the list of indeces that map each DICOM in FromDicoms to 
    # the DICOM in ToDicoms (or MappedDicoms) with nearest scan position:
    Inds = []
    
    # Initialise a list of the scan position differences between each Pos1 and 
    # Pos2, based on the absolute minimum positional differences between each 
    # Pos1 and and all Pos2:
    ScanPosDiffs = []
    
    # Loop for every DICOM object in FromDicoms:
    for FromDicom in FromDicoms:
        # Get the Image Position Patient in FromDicom:
        Pos1 = FromDicom.ImagePositionPatient
        
        ## Initialise the list of indeces that map each DICOM in FromDicoms to 
        ## the DICOM in ToDicoms (or MappedDicoms) with nearest scan position:
        #Inds = []
        
        # Initialise a list of the positional differences between Pos1 and 
        # every Pos2:
        Dpos = []

            
        # Loop through each DICOM in ToDicoms or MappedDicoms:
        """ Note: ToDicoms and MappedDicoms will have equal length. """
        for i in range(len(ToDicoms)):
            if Shift:
                Pos2 = ShiftedDicoms[i].ImagePositionPatient
            else:
                Pos2 = ToDicoms[i].ImagePositionPatient
                
            ## Initialise a list of the positional differences between Pos1 and 
            ## every Pos2:
            #Dpos = []
        
            # Calculate the distance between Pos1 and Pos2 along x, y, z and r:
            dx = Pos2[0] - Pos1[0]
            dy = Pos2[1] - Pos1[1]
            dz = Pos2[2] - Pos1[2]
            dr = (dx**2 + dy**2 + dz**2)**.5
            
            # Append the positional difference along x, y, z or r, depending
            # on SearchBy, to Dpos:
            if SearchBy=='x':
                Dpos.append(dx)
            if SearchBy=='y':
                Dpos.append(dy)
            if SearchBy=='z':
                Dpos.append(dz)
            if SearchBy=='r':
                Dpos.append(dr)
            
            
            
        # Find the index of the absolution minimum positional difference along 
        # the chosen direction:
        Ind = [abs(dpos) for dpos in Dpos].index(min([abs(dpos) for dpos in Dpos]))
        
        # The scan positional difference:
        ScanPosDiffs.append(Dpos[Ind])
        
        # Store the index:
        Inds.append(Ind)


    # The key that will be used to update DataDict:
    if Shift:
        UpdateKey = 'MapPosCorr'
    else:
        UpdateKey = 'Mapped'
        
    # Store the indeces and the absolute minimum positional differences in 
    # DataDict:
    #DataDict[UpdateKey].update({'Inds':Inds})
    #DataDict[UpdateKey].update({'AbsMinDpos':AbsMinDpos})
    
    if Debug:
        print(f'\nInds (N = {len(Inds)}) =', Inds)
        print(f'\nScanPosDiffs (N = {len(ScanPosDiffs)}) =', [round(item,2) for item in ScanPosDiffs])
        
    # Create a new list of DICOM objects that are the same length as FromDicoms
    # with the closest DICOM in ToDicoms or ShiftedDicoms as determined from
    # the list of indeces Inds, and corresponding filepaths:
    MappedDicoms = []
    MappedDicomFpaths = []
    
    """
    The Study Instance UID and Frame of Reference UID of the mapped DICOMs will 
    be kept unchanged from those of the position-corrected DICOMs, but a new
    Series Instance UID will need to be generated.  
    Unique SOP Instance UIDs will be generated for each individual DICOM.  
    I'm not sure if the Frame of Reference UID of MappedDicoms should be 
    made equal to the FoR UID of the FromDicoms.  For now I'll assume not.
    """
        
    # Create a new Series Instance UID:
    NewSeriesUid = pydicom.uid.generate_uid()
    
    # Create a new Series Description by appending ' - PosCorr' to the 
    # Series Description of ToDicoms:
    if Shift:
        MappedSeriesDesc = ToSeriesDesc + '- ' + SearchBy + 'PosCorr - Mapped'
    else:
        MappedSeriesDesc = ToSeriesDesc + ' - Mapped'
        
    
    # Modify PosCorrDicomDir to include the Series Description in parentheses:
    MappedDicomDir = MappedDicomDir + ' (' + MappedSeriesDesc + ')'
    
    # Modify PosCorrDicomDir in DataDict:
    if Shift:
        DataDict['MapPosCorr']['DicomDir'] = MappedDicomDir
    else:
        DataDict['Mapped']['DicomDir'] = MappedDicomDir   
        
    
    # Create the directory:
    CreateDir(MappedDicomDir, Debug)
        
    """
    Some slices in MappedDicoms may be repeated because of repeated indeces in
    Inds. So use IndecesCount() to create a list of unique indeces in string
    format.
    """
    # Use IndecesCount() to create a list of unique indeces in string format
    # with '-1', '-2', '-3', etc. appended to the index for each occurrance:
    IndsStr = IndecesCount(Inds)
    
    # Loop for every index in Inds:
    for i in range(len(Inds)):
        if Shift:
            MappedDicom = copy.deepcopy(ShiftedDicoms[Inds[i]])
            
            #print(f'\nInds[{i}] = {Inds[i]}')
            #print(f'\nPosCorrDicomFpaths[Inds[{i}]] = {PosCorrDicomFpaths[Inds[i]]}')
            
            # Get the original file name with extension:
            OrigFname = os.path.basename(PosCorrDicomFpaths[Inds[i]])
        
        else:
            MappedDicom = copy.deepcopy(ToDicoms[Inds[i]])
            
            # Get the original file name with extension:
            OrigFname = os.path.basename(ToDicomFpaths[Inds[i]])
            
        # Get the file name without extension:
        OrigFnameNoExt = os.path.splitext(OrigFname)[0]
            
        # Generate a new SOP Instance UID:
        NewSopUid = pydicom.uid.generate_uid()
        
        # Modify the old SOP UID with the new one:
        MappedDicom.SOPInstanceUID = NewSopUid
        
        # Modify the old Series UIDs with the new one:
        MappedDicom.SeriesInstanceUID = NewSeriesUid
        
        # Modify the Series Number:
        MappedDicom.SeriesNumber = MappedSeriesNo
        
        # Modify the Series Description:
        MappedDicom.SeriesDescription = MappedSeriesDesc
        
        # Append MappedDicom to MappedDicoms:
        MappedDicoms.append(MappedDicom)
        
        # Create a new filename, adding information about the mapping index:
        #MappedDicomFname = OrigFnameNoExt + '-mapInd-' + str(Inds[i]) + '.dcm'
        MappedDicomFname = OrigFnameNoExt + '-mapInd-' + IndsStr[i] + '.dcm'
        
        if Debug:
            print('\n   Changing the file name from:\n', \
                  OrigFname, '\nto:\n', MappedDicomFname)
        
        # Create new filepath:
        #if Shift:
        #    MappedDicomFpath = os.path.join(MapPosCorrDicomDir, MappedDicomFname)
        #else:
        #    MappedDicomFpath = os.path.join(MappedDicomDir, MappedDicomFname)
        MappedDicomFpath = os.path.join(MappedDicomDir, MappedDicomFname)
        
        # Append to mappedFpaths:
        MappedDicomFpaths.append(MappedDicomFpath)
        
        
        if Debug:
            print('\n    Exporting the mapped position corrected DICOM to:\n', \
                  MappedDicomFpath)
        
        # Export the DICOM file:
        MappedDicom.save_as(MappedDicomFpath)
        
        
    # The key that will be used to update DataDict:
    if Shift:
        UpdateKey = 'MapPosCorr'
    else:
        UpdateKey = 'Mapped'
    
    
    # Slice numbers from ToDicoms:
    ToSliceNos = DataDict['To']['SliceNos']
    
    # The slice numbers for the MappedDicoms:
    MappedSliceNos = [ToSliceNos[i] for i in Inds]
    
    # Update DataDict:
    DataDict[UpdateKey].update({'Inds':Inds})
    DataDict[UpdateKey].update({'IndsStr':IndsStr})
    DataDict[UpdateKey].update({'ScanPosDiffs':ScanPosDiffs})
    DataDict[UpdateKey].update({'DicomFpaths':MappedDicomFpaths})
    DataDict[UpdateKey].update({'SeriesDesc':MappedSeriesDesc})
    DataDict[UpdateKey].update({'SliceNos':MappedSliceNos})
        
    
    
    
    return ShiftedDicoms, MappedDicoms, DataDict

