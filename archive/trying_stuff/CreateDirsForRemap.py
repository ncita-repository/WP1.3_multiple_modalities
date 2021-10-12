# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:32:20 2020

@author: ctorti
"""



def CreateDirsForRemap(DataDict, KeyToRemap):
    # Import packages:
    import os
    
    
    # Get some necessary info from dataDict:
    RootDir = DataDict['RootDir']
    Debug = DataDict['Debug']
    #KeyToRemap = DataDict['KeyToRemap']
    
    # Use KeyToRemap to index the series that will be remapped:
    PreRemappedSeriesNo = DataDict[KeyToRemap]['SeriesNo']
    PreRemappedSeriesDesc = DataDict[KeyToRemap]['SeriesDesc']
    
  
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
        Format of Series Numbers will be the PreRemappedSeriesNo + '1' to 
        indicate that a remapping was applied:
    """
    
    # Define the Series number and Series description, Directory name, DICOM
    # directory and ROI directory for the remapped DICOMs:
    RemappedSeriesNo = PreRemappedSeriesNo + '1'
    RemappedSeriesDesc = PreRemappedSeriesDesc + ' Remapped'
    RemappedDirName = RemappedSeriesNo + ' ' + RemappedSeriesDesc
    RemappedDicomDir = os.path.join(RootDir, RemappedDirName + ' DICOMs')
    RemappedRoiDir = os.path.join(RootDir, RemappedDirName + ' ROIs')
        
    # Add to AllDirs:
    AllDirs.append(RemappedRoiDir)
    AllDirs.append(RemappedDicomDir)
    
    
    # Create a new key for DataDict that indicates which key was remapped:
    RemappedKey = 'Remapped' + KeyToRemap
    
    # Store RemappedKey in DataDict:     
    #DataDict.update({'RemappedKey':RemappedKey})    
    
    # Store RoiDir, DicomDir, SeriesNo, SeriesDesc and SliceNos in DataDict for
    # the key given by RemappedKey:
    DataDict.update({RemappedKey:{#'ScanLabel':MapPosCorrScanLabel,\
                                  'RoiDir':RemappedRoiDir,\
                                  'DicomDir':RemappedDicomDir,\
                                  'SeriesNo':RemappedSeriesNo,\
                                  'SeriesDesc':RemappedSeriesDesc
                                    }
                     })
            
    
    
    # CREATE THE DIRECTORIES IF THEY DON'T EXIST:
    for DirPath in AllDirs:
        if Debug:
            print('\nChecking if directory', DirPath, 'exists')
        try:
            os.mkdir(DirPath)
            if Debug:
                print('\nDirectory', DirPath, 'created')
        except FileExistsError:
            if Debug:
                print('\nDirectory', DirPath, 'already exists')

    
    
    return DataDict