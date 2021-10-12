# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 08:41:43 2020

@author: ctorti
"""


"""
Function:
    
    GetClosestSlices()
    
Purpose:
    
    Starting with two DICOM objects SourceDicoms and TargetDicoms, find the 
    indeces that map each DICOM in SourceDicoms to the DICOM in TargetDicoms 
    with the nearest position along the direction specified by 
    DataDict['SearchBy']. Likewise, find the indeces that map each DICOM in 
    TargetDicoms to the nearest DICOM in SourceDicoms.
    
    Two sets of mappings will be created - one which lists all indeces as
    described above.  A second will omit any indeces that correspond to scan
    positional differences greater than a set threshold distance.
    
Input:
    
    SourceDicoms  - list of DICOM objects to be mapped to TargetDicoms
    
    TargetDicoms  - list of DICOM objects 
    
    DataDict      - dictionary containing several variables
                
    
Returns:
                      
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
    
    
    # Define the maximum scan positional difference allowed for the reduced
    # mapping indeces:
    MaxPosDiff = 1 # mm

    
    
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
    # nearest DICOM in SourceDicoms for both mapping sets:
    TargetToSourceNearestInds = [] # nearest indeces
    TargetToSourceNearOnlyInds = [] # "near" (within MaxPosDiff) only indeces 
    
    # Initialise a list of the scan position differences between each DICOM in
    # TargetDicoms and the nearest DICOM in SourceDicoms, based on the absolute 
    # minimum positional differences between each set of positions:
    TargetToSourceNearestDiffs = []
    TargetToSourceNearOnlyDiffs = []
    
    # Initialise a list of indeces for TargetDicoms for which there is a slice
    # within SourceDicoms whose scan position is within MaxPosDiff, and the
    # corresponding positional differences:
    TargetNearOnlyInds = []
    TargetNearOnlyDiffs = []
    
    # Get the number of DICOMs in SourceDicoms and TargetDicoms:
    S = len(SourceDicoms)
    T = len(TargetDicoms)
    
    # Loop for every DICOM object in TargetDicoms:
    for t in range(T):
        # Get the Image Position Patient in TargetDicoms[t]:
        TargetPosition = TargetDicoms[t].ImagePositionPatient
        
        # Initialise a list of the positional differences between 
        # TargetPosition and every position in SourceDicoms:
        Dpos = []
            
        # Loop through each DICOM in SourceDicoms:
        for s in range(S):
            SourcePosition = SourceDicoms[s].ImagePositionPatient
        
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
        
        # Append Ind and Dpos[Ind] to TargetToSourceNearestInds and 
        # TargetToSourceNearestDiffs:
        TargetToSourceNearestInds.append(Ind)
        
        TargetToSourceNearestDiffs.append(Dpos[Ind])
        
        # Append Ind TargetToSourceNearOnlyInds, Dpos[Ind] to 
        # TargetToSourceNearOnlyDiffs, and t to TargetNearOnlyInds only if 
        # abs(Dpos[Ind]) <= MaxPosDiff:
        if abs(Dpos[Ind]) <= MaxPosDiff:
            TargetToSourceNearOnlyInds.append(Ind)
            
            TargetToSourceNearOnlyDiffs.append(Dpos[Ind])
            
            TargetNearOnlyInds.append(t)
            
            TargetNearOnlyDiffs.append(Dpos[Ind])
        
        
        
        
    """ NEXT GET MAPPING FROM SOURCE TO TARGET : """
    
    # For every DICOM in SourceDicoms, find the index of the nearest slice in  
    # the list of DICOMs TargetDicoms:
    
    # Initialise the list of indeces that map the DICOMs in SourceDicoms to the 
    # nearest DICOM in TargetDicoms:
    SourceToTargetNearestInds = []
    SourceToTargetNearOnlyInds = []
    
    # Initialise a list of the scan position differences between each DICOM in
    # SourceDicoms and the nearest DICOM in TargetDicoms, based on the absolute 
    # minimum positional differences between each set of positions:
    SourceToTargetNearestDiffs = []
    SourceToTargetNearOnlyDiffs = []
    
    # Initialise a list of indeces for SourceDicoms for which there is a slice
    # within TargetDicoms whose scan position is within MaxPosDiff, and the
    # corresponding positional differences:
    SourceNearOnlyInds = []
    SourceNearOnlyDiffs = []
    
    # Loop for every DICOM object in SourceDicoms:
    for s in range(S):
        # Get the Image Position Patient in SourceDicoms[s]:
        SourcePosition = SourceDicoms[s].ImagePositionPatient
        
        # Initialise a list of the positional differences between 
        # SourcePosition and every position in TargetDicoms:
        Dpos = []
            
        # Loop through each DICOM in TargetDicoms:
        for t in range(T):
            TargetPosition = TargetDicoms[t].ImagePositionPatient
        
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
        
        # Append Ind and Dpos[Ind] to SourceToTargetNearestInds and 
        # SourceToTargetNearestDiffs:
        SourceToTargetNearestInds.append(Ind)
        
        SourceToTargetNearestDiffs.append(Dpos[Ind])
        
        # Append Ind to SourceToTargetNearOnlyInds, Dpos[Ind] to 
        # SourceToTargetNearOnlyDiffs, and s to SourceNearOnlyInds only if 
        # abs(Dpos[Ind]) <= MaxPosDiff:
        if abs(Dpos[Ind]) <= MaxPosDiff:
            SourceToTargetNearOnlyInds.append(Ind)
            
            SourceToTargetNearOnlyDiffs.append(Dpos[Ind])
            
            SourceNearOnlyInds.append(s)
            
            SourceNearOnlyDiffs.append(Dpos[Ind])
            
        

    if True:#Debug:
        print(f'\nThere are {len(SourceDicoms)} slices in SourceDicoms,', \
              f'and {len(TargetDicoms)} in TargetDicoms.')
        
        # Print the "Nearest" indeces and scan positional differences:
        
        print(f'\nSourceToTargetNearestInds (N = {len(SourceToTargetNearestInds)}):\n\n',\
              SourceToTargetNearestInds)
        
        print(f'\nSourceToTargetNearestDiffs (N = {len(SourceToTargetNearestDiffs)}):\n\n',\
              [round(item,2) for item in SourceToTargetNearestDiffs])
        
        print(f'\nTargetToSourceNearestInds (N = {len(TargetToSourceNearestInds)}):\n\n',\
              TargetToSourceNearestInds)
        
        print(f'\nTargetToSourceNearestDiffs (N = {len(TargetToSourceNearestDiffs)}):\n\n',\
              [round(item,2) for item in TargetToSourceNearestDiffs])
        
        # Print the "NearOnly" indeces and scan positional differences:
        
        print(f'\nSourceToTargetNearOnlyInds (N = {len(SourceToTargetNearOnlyInds)}):\n\n',\
              SourceToTargetNearOnlyInds)
        
        print(f'\nSourceToTargetNearOnlyDiffs (N = {len(SourceToTargetNearOnlyDiffs)}):\n\n',\
              [round(item,2) for item in SourceToTargetNearOnlyDiffs])
        
        print(f'\nTargetToSourceNearOnlyInds (N = {len(TargetToSourceNearOnlyInds)}):\n\n',\
              TargetToSourceNearOnlyInds)
        
        print(f'\nTargetToSourceNearOnlyDiffs (N = {len(TargetToSourceNearOnlyDiffs)}):\n\n',\
              [round(item,2) for item in TargetToSourceNearOnlyDiffs])
        
        print(f'\nSourceNearOnlyInds (N = {len(SourceNearOnlyInds)}):\n\n',\
              SourceNearOnlyInds)
        
        print(f'\nSourceNearOnlyDiffs (N = {len(SourceNearOnlyDiffs)}):\n\n',\
              [round(item,2) for item in SourceNearOnlyDiffs])
        
        print(f'\nTargetNearOnlyInds (N = {len(TargetNearOnlyInds)}):\n\n',\
              TargetNearOnlyInds)
        
        print(f'\nTargetNearOnlyDiffs (N = {len(TargetNearOnlyDiffs)}):\n\n',\
              [round(item,2) for item in TargetNearOnlyDiffs])
        
      
        
    # Store SourceToTargetNearestInds, SourceToTargetNearestDiffs 
    # TargetToSourceNearestInds and TargetToSourceNearestDiffs in DataDict 
    # with the keys 'SourceToTargetNearest' and 'TargetToSourceNearest', and
    # store SourceNearOnlyInds and TargetNearOnlyInds with keys of the same 
    # name:
    DataDict.update({'SourceToTargetNearest':{'MappingInds':SourceToTargetNearestInds, \
                                              'ScanPosDiffs':SourceToTargetNearestDiffs
                                              }, \
                    'TargetToSourceNearest':{'MappingInds':TargetToSourceNearestInds, \
                                              'ScanPosDiffs':TargetToSourceNearestDiffs
                                             }, \
                    'SourceToTargetNearOnly':{'MappingInds':SourceToTargetNearOnlyInds, \
                                              'ScanPosDiffs':SourceToTargetNearOnlyDiffs
                                              }, \
                    'TargetToSourceNearOnly':{'MappingInds':TargetToSourceNearOnlyInds, \
                                              'ScanPosDiffs':TargetToSourceNearOnlyDiffs
                                              }, \
                    'SourceNearOnly':{'MappingInds':SourceToTargetNearOnlyInds, \
                                      'ScanPosDiffs':SourceToTargetNearOnlyDiffs
                                      }, \
                    'TargetNearOnly':{'MappingInds':TargetNearOnlyInds, \
                                      'ScanPosDiffs':TargetNearOnlyDiffs
                                      }
                    })
   
    
    return DataDict

