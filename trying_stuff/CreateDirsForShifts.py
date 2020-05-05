# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 22:05:11 2020

@author: ctorti
"""



def CreateDirsForShifts(DataDict):
    # Import packages:
    import os
    
    
    # Get some necessary info from dataDict:
    RootDir = DataDict['RootDir']
    Shift = DataDict['Shift'] # e.g. -158.199231349734
    ShiftFeature = DataDict['ShiftFeature'] # e.g. 'using lens of eye'
    ToShiftKey = DataDict['ToShiftKey'] # e.g. 'Target'
    #searchBy = dataDict['SearchBy']
    #RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug'] # e.g. False
    
    #FromSeriesNo = DataDict['From']['SeriesNo']
    #ToSeriesNo = DataDict['To']['SeriesNo']
    #ToSeriesDesc = DataDict['To']['SeriesDesc']
    #ToSliceNos = DataDict['To']['SliceNos']
    
    # Use ToShiftKey to index the series that will be shifted:
    PreShiftedSeriesNo = DataDict[ToShiftKey]['SeriesNo']
    PreShiftedSeriesDesc = DataDict[ToShiftKey]['SeriesDesc']
    PreShiftedSliceNos = DataDict[ToShiftKey]['SliceNos']
    
    if Debug:
        print('\n PreShiftedSeriesDesc =', PreShiftedSeriesDesc)
    
    # Proceed only if Shift is non-zero:
    if Shift:
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
            Format of Series Numbers will be the PreShiftedSeriesNo + '0' + a number
            of binary digits which indicate that scan-position-shifting was 
            applied:
        """
        
        # DICOMs whose scan positions are shifted:
        #ShiftedSeriesNo = PreShiftedSeriesNo + '0' + '1'
        ShiftedSeriesNo = PreShiftedSeriesNo + '0'
        #ShiftedSeriesDesc = 'Shifted ' + PreShiftedSeriesDesc
        ShiftedSeriesDesc = PreShiftedSeriesDesc + ' Shifted using ' + ShiftFeature
        ShiftedDirName = ShiftedSeriesNo + ' ' + ShiftedSeriesDesc
        ShiftedDicomDir = os.path.join(RootDir, ShiftedDirName + ' DICOMs')
        ShiftedRoiDir = os.path.join(RootDir, ShiftedDirName + ' ROIs')
        
        if Debug:
            print('\n ShiftedSeriesDesc =', ShiftedSeriesDesc)
            
        # Add to AllDirs:
        AllDirs.append(ShiftedRoiDir)
        AllDirs.append(ShiftedDicomDir)
        
        
        # Create a new key for DataDict that indicates which key was shifted:
        ShiftedKey = 'Shifted' + ToShiftKey # e.g. 'Shifted' + 'Target' = 'ShiftedTarget'
        
        # Update DataDict:     
        DataDict.update({'ToShiftKey':ToShiftKey, \
                         'ShiftedKey':ShiftedKey,})      
    
        DataDict.update({ShiftedKey:{#'ScanLabel':MapPosCorrScanLabel,\
                                     'RoiDir':ShiftedRoiDir,\
                                     'DicomDir':ShiftedDicomDir,\
                                     'SeriesNo':ShiftedSeriesNo,\
                                     'SeriesDesc':ShiftedSeriesDesc,\
                                     'SliceNos':PreShiftedSliceNos # make ='PreShifted'
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