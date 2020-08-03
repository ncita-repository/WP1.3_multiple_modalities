# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:08:26 2020

@author: ctorti
"""





def CreateDirsForDownloads(DataDict):
    # Import packages:
    import os
    
    
    # Get some necessary info from dataDict:
    RootDir = DataDict['RootDir']
    #FromScanLabel = DataDict['From']['ScanLabel']
    #ToScanLabel = DataDict['To']['ScanLabel']
    SourceSeriesNo = DataDict['Source']['SeriesNo']
    TargetSeriesNo = DataDict['Target']['SeriesNo']
    #SourceSeriesDesc = 'Source'
    #TargetSeriesDesc = 'Target'
    SourceSeriesDesc = DataDict['Source']['SeriesDesc']
    TargetSeriesDesc = DataDict['Target']['SeriesDesc']
    #Shift = DataDict['Shift']
    #searchBy = dataDict['SearchBy']
    #RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug']
    
    
    """
    DEFINE THE SCAN LABELS AND DIRECTORY NAMES:
    """
    
    """
    For the DICOMs and DICOM-RTSTRUCT of the scan whose ROIs are to be copied:
    """
    
    # Where the DICOM-RTSTRUCT to be copied is to be saved to:
    #FromRoiDir = os.path.join(RootDir, FromScanLabel + ' ROIs')
    #SourceRoiDir = os.path.join(RootDir, SourceSeriesNo + ' ROIs')
    SourceRoiDir = os.path.join(RootDir, SourceSeriesNo + ' ' + \
                                SourceSeriesDesc + ' ROIs (source)')
    
    # Where the DICOM scan that relate to the ROIs to be copied are to be saved 
    # to:
    #FromDicomDir = os.path.join(RootDir, FromScanLabel + ' DICOMs')
    #SourceDicomDir = os.path.join(RootDir, SourceSeriesNo + ' DICOMs')
    SourceDicomDir = os.path.join(RootDir, SourceSeriesNo + ' ' + \
                                  SourceSeriesDesc + ' DICOMs (source)')
    
    # Add the above directories to dataDict:
    DataDict['Source'].update({'RoiDir':SourceRoiDir,\
                             'DicomDir':SourceDicomDir,\
                             #'ScanLabel':FromScanLabel,\
                             'SeriesNo':SourceSeriesNo,\
                             'SeriesDesc':SourceSeriesDesc,\
                             'SliceNos':[]
                             })
    
    """
    Note:  The entry 'SeriesDesc' is a place-holder that will be added to all
    of the main keys in DataDict.  
    It will later be populated with the DICOMs' Series Description. 
    """
    
    
    """ 
    For the DICOMs and DICOM-RTSTRUCT of the scan that the ROIs are to be 
    copied to:
    """
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    #ToDicomDir = os.path.join(RootDir, ToScanLabel + ' DICOMs')
    #TargetDicomDir = os.path.join(RootDir, TargetSeriesNo + ' DICOMs')
    TargetDicomDir = os.path.join(RootDir, TargetSeriesNo + ' ' + \
                                  TargetSeriesDesc + ' DICOMs (target)')
    
    # The new (modified copy) ROIs are to be stored at:
    #ToRoiDir = os.path.join(RootDir, ToScanLabel + ' ROIs')
    #TargetRoiDir = os.path.join(RootDir, TargetSeriesNo + ' ROIs')
    TargetRoiDir = os.path.join(RootDir, TargetSeriesNo + ' ' + \
                                TargetSeriesDesc + ' ROIs (target)')
    
    
    # Add the above directories to dataDict:
    DataDict['Target'].update({#'ScanLabel':ToScanLabel,\
                           'RoiDir':TargetRoiDir,\
                           'DicomDir':TargetDicomDir,\
                           'SeriesNo':TargetSeriesNo,\
                           'SeriesDesc':TargetSeriesDesc,\
                           'SliceNos':[]
                           })
    

    
    """
    For efficiency, create a list of all directories that need to be created,
    AllDirs, and add the directories defined above.
    """
    AllDirs = [RootDir, SourceRoiDir, SourceDicomDir, TargetRoiDir, TargetDicomDir]
    
    
        
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