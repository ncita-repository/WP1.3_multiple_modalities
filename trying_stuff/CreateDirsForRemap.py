# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:52:34 2020

@author: ctorti
"""



def CreateDirsForRemap(DataDict):
    # Import packages:
    import os
    
    
    # Get some necessary info from dataDict:
    RootDir = DataDict['RootDir']
    ToRemapKey = DataDict['ToRemapKey']
    ForRemapKey = DataDict['ForRemapKey']
    #Shift = DataDict['Shift']
    #searchBy = dataDict['SearchBy']
    #RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug']
    
    #ToRemapKey = DataDict['ToRemapKey']
    
    # Use ToRemapKey to index the series that will be remapped:
    PreRemappedSeriesNo = DataDict[ToRemapKey]['SeriesNo']
    PreRemappedSeriesDesc = DataDict[ToRemapKey]['SeriesDesc']
    PreRemappedSliceNos = DataDict[ToRemapKey]['SliceNos']
    
    # Use ForRemapKey to get details of the series that the remapping is to
    # match in slice number:
    RemapForSeriesNo = DataDict[ForRemapKey]['SeriesNo']
    RemapForSeriesDesc = DataDict[ForRemapKey]['SeriesDesc']
    
  
    """
    For efficiency, create a list of all directories that need to be 
    created, AllDirs, and add the directories defined below:
    """
    AllDirs = []
    
    
    
    """
    DEFINE THE SCAN LABELS AND DIRECTORY NAMES:
    """

    """
    Note: 
        Format of Series Numbers will be the PreRemappedSeriesNo + '0' + '1' to 
        indicate that a remapping was applied:
    """
    
    # DICOMs whose scan positions are shifted:
    #RemappedSeriesNo = PreRemappedSeriesNo + '0' + '1'
    #RemappedSeriesNo = PreRemappedSeriesNo + '1'
    #RemappedSeriesNo = PreRemappedSeriesNo + '010' + RemapForSeriesNo
    RemappedSeriesNo = PreRemappedSeriesNo + '0' + RemapForSeriesNo
    #RemappedSeriesDesc = PreRemappedSeriesDesc + ' Remapped' 
    #RemappedSeriesDesc = 'Remapped ' + PreRemappedSeriesDesc
    RemappedSeriesDesc = PreRemappedSeriesDesc + ' Remapped for ' + RemapForSeriesDesc
    RemappedDirName = RemappedSeriesNo + ' ' + RemappedSeriesDesc
    RemappedDicomDir = os.path.join(RootDir, RemappedDirName + ' DICOMs')
    RemappedRoiDir = os.path.join(RootDir, RemappedDirName + ' ROIs')
        
    # Add to AllDirs:
    AllDirs.append(RemappedRoiDir)
    AllDirs.append(RemappedDicomDir)
    
    
    # Create a new key for DataDict that indicates which key was shifted:
    RemappedKey = 'Remapped' + ToRemapKey
    
    # Update DataDict:     
    DataDict.update({'ToRemapKey':ToRemapKey, \
                     'RemappedKey':RemappedKey,})      

    DataDict.update({RemappedKey:{#'ScanLabel':MapPosCorrScanLabel,\
                                  'RoiDir':RemappedRoiDir,\
                                  'DicomDir':RemappedDicomDir,\
                                  'SeriesNo':RemappedSeriesNo,\
                                  'SeriesDesc':RemappedSeriesDesc,\
                                  'SliceNos':PreRemappedSliceNos # make ='PreRemapped' 
                                    }
                     })
            
    
    
    # CREATE THE DIRECTORIES IF THEY DON'T EXIST:
    for dirPath in AllDirs:
        if Debug:
            print('\nChecking if directory', dirPath, 'exists')
        try:
            os.mkdir(dirPath)
            if Debug:
                print('\nDirectory', dirPath, 'created')
        except FileExistsError:
            if Debug:
                print('\nDirectory', dirPath, 'already exists')

    
    
    return DataDict