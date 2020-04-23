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


def GetClosestSlices(SourceDicoms, TargetDicoms, DataDict):
    
    # Import packages:
    #import copy
    #import pydicom
    #import os
    #from IndecesCount import IndecesCount
    #from CreateDir import CreateDir
    
    # Get some necessary info from DataDict:
    SearchBy = DataDict['SearchBy']
    #Shift = DataDict['Shift']
    Debug = DataDict['Debug']
    #ToDicomFpaths = DataDict['To']['DicomFpaths']
    #ShiftedDicomDir = DataDict['Shifted']['DicomDir']
    ##ShiftedScanLabel = DataDict['Shifted']['ScanLabel']
    #ShiftedSeriesNo = DataDict['Shifted']['SeriesNo']
    #if Shift:
    #    #MapShiftedDicomDir = DataDict['MapShifted']['DicomDir']
    #    MappedDicomDir = DataDict['MapShifted']['DicomDir']
    #    #MapShiftedScanLabel = DataDict['MapShifted']['ScanLabel']
    #    #MapShiftedSeriesNo = DataDict['MapShifted']['SeriesNo']
    #    MappedSeriesNo = DataDict['MapShifted']['SeriesNo']
    #else:
    #    MappedDicomDir = DataDict['Mapped']['DicomDir']
    #    #MappedScanLabel = DataDict['Mapped']['ScanLabel']
    #    MappedSeriesNo = DataDict['Mapped']['SeriesNo']
    
    

    
    
    if Debug:
        print(f'\nThere are {len(SourceDicoms)} slices in SourceDicoms.')
        print(f'There are {len(TargetDicoms)} slices in TargetDicoms.')
        #if Shift:
        #    print(f'There are {len(ShiftedDicoms)} slices in ShiftedDicoms.')
        #    print(f'There are {len(ShiftedDicomFpaths)} fpaths in ShiftedDicomFpaths.')
            
    
    """ FIRST GET MAPPING FROM TARGET TO SOURCE : """
    
    # For every DICOM in TargetDicoms, find the index of the nearest slice in  
    # the list of DICOMs SourceDicoms:
    
    # Initialise the list of indeces that map the DICOMs in TargetDicoms to the 
    # nearest DICOM in SourceDicoms:
    TargetToSourceMapping = []
    
    # Initialise a list of the scan position differences between each DICOM in
    # TargetDicoms and the nearest DICOM in SourceDicoms, based on the absolute 
    # minimum positional differences between each set of positions:
    TargetToSourceScanPosDiffs = []
    
    # Loop for every DICOM object in TargetDicoms:
    for TargetDicom in TargetDicoms:
        # Get the Image Position Patient in TargetDicom:
        TargetPosition = TargetDicom.ImagePositionPatient
        
        # Initialise a list of the positional differences between 
        # TargetPosition and every position in SourceDicoms:
        Dpos = []
            
        # Loop through each DICOM in SourceDicoms:
        for i in range(len(SourceDicoms)):
            SourcePosition = SourceDicoms[i].ImagePositionPatient
        
            # Calculate the distance from TargetPosition to SourcePosition 
            # along x, y, z and r:
            dx = SourcePosition[0] - TargetPosition[0]
            dy = SourcePosition[1] - TargetPosition[1]
            dz = SourcePosition[2] - TargetPosition[2]
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
        
        # Store the mapping index:
        TargetToSourceMapping.append(Ind)
        
        # and the scan positional difference:
        TargetToSourceScanPosDiffs.append(Dpos[Ind])
        
        
        
        
    """ NEXT GET MAPPING FROM SOURCE TO TARGET : """
    
    # For every DICOM in SourceDicoms, find the index of the nearest slice in  
    # the list of DICOMs TargetDicoms:
    
    # Initialise the list of indeces that map the DICOMs in SourceDicoms to the 
    # nearest DICOM in TargetDicoms:
    SourceToTargetMapping = []
    
    # Initialise a list of the scan position differences between each DICOM in
    # SourceDicoms and the nearest DICOM in TargetDicoms, based on the absolute 
    # minimum positional differences between each set of positions:
    SourceToTargetScanPosDiffs = []
    
    # Loop for every DICOM object in SourceDicoms:
    for SourceDicom in SourceDicoms:
        # Get the Image Position Patient in SourceDicom:
        SourcePosition = SourceDicom.ImagePositionPatient
        
        # Initialise a list of the positional differences between 
        # SourcePosition and every position in TargetDicoms:
        Dpos = []
            
        # Loop through each DICOM in TargetDicoms:
        for i in range(len(TargetDicoms)):
            TargetPosition = TargetDicoms[i].ImagePositionPatient
        
            # Calculate the distance from SourcePosition to TargetPosition 
            # along x, y, z and r:
            dx = TargetPosition[0] - SourcePosition[0]
            dy = TargetPosition[1] - SourcePosition[1]
            dz = TargetPosition[2] - SourcePosition[2]
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
        
        # Store the mapping index:
        SourceToTargetMapping.append(Ind)
        
        # and the scan positional difference:
        SourceToTargetScanPosDiffs.append(Dpos[Ind])
        
        
        

    if True:#Debug:
        print(f'\nThere are {len(SourceDicoms)} slices in SourceDicoms,', \
              f'and {len(TargetDicoms)} in TargetDicoms.')
        print(f'\nSourceToTargetMapping (N = {len(SourceToTargetMapping)}):\n\n',\
              SourceToTargetMapping)
        print(f'\nSourceToTargetScanPosDiffs (N = {len(SourceToTargetScanPosDiffs)}):\n\n',\
              [round(item,2) for item in SourceToTargetScanPosDiffs])
        print(f'\nTargetToSourceMapping (N = {len(TargetToSourceMapping)}):\n\n',\
              TargetToSourceMapping)
        print(f'\nTargetToSourceScanPosDiffs (N = {len(TargetToSourceScanPosDiffs)}):\n\n',\
              [round(item,2) for item in TargetToSourceScanPosDiffs])
      
        
    # Update DataDict:
    DataDict.update({'SourceToTargetMapping':{'MappingInds':SourceToTargetMapping, \
                                              'ScanPosDiffs':SourceToTargetScanPosDiffs
                                              }, \
                    'TargetToSourceMapping':{'MappingInds':TargetToSourceMapping, \
                                              'ScanPosDiffs':TargetToSourceScanPosDiffs
                                             }                     
                    })
   
    
    return DataDict

