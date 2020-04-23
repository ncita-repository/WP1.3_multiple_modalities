# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:46:00 2020

@author: ctorti
"""



def CreateDirsForReg(DataDict):
    # Import packages:
    import os
    
    
    # Get some necessary info from DataDict:
    RootDir = DataDict['RootDir']
    Shift = DataDict['Shift']
    #searchBy = dataDict['SearchBy']
    RegMethod = DataDict['RegMethod']
    FixedKey = DataDict['FixedKey']
    MovingKey = DataDict['MovingKey']
    RegisteredKey = DataDict['RegisteredKey']
    Debug = DataDict['Debug']
    
    if False:
        #FromSeriesNo = DataDict['From']['SeriesNo']
        ToSeriesNo = DataDict['To']['SeriesNo']
        ToSeriesDesc = DataDict['To']['SeriesDesc']
        ToSliceNos = DataDict['To']['SliceNos']
    
    FixedSeriesNo = DataDict[FixedKey]['SeriesNo']
    FixedSeriesDesc = DataDict[FixedKey]['SeriesDesc']
    FixedSliceNos = DataDict[FixedKey]['SliceNos']
    
    MovingSeriesNo = DataDict[MovingKey]['SeriesNo']
    MovingSeriesDesc = DataDict[MovingKey]['SeriesDesc']
    MovingSliceNos = DataDict[MovingKey]['SliceNos']
    
    
    # Proceed only if RegfMethod one of the expected registration methods:
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        
        """
        For efficiency, create a list of all directories that need to be created,
        AllDirs, and add the directories defined below:
        """
        AllDirs = []
        
        
        
        """
        DEFINE THE SCAN LABELS AND DIRECTORY NAMES:
        """
    
        """ Below section in if statement simplified below.. """
        if False:
            """
            Note: 
                Format of Series Numbers will be the ToSeriesNo + '0' + a number of
                binaries which indicate whether position-correction was done,
                mapping, registration, etc.
            """
                
            # If Shift is non-zero:
            if Shift:
                # Get some necessary info from DataDict:
                ShiftedSeriesNo = DataDict['Shifted']['SeriesNo']
                ShiftedSeriesDesc = DataDict['Shifted']['SeriesDesc']
                ShiftedSliceNos = DataDict['Shifted']['SliceNos']
                
                # Add 1 after ShiftedSeriesNo then 1/2/3 for RegMethod:
                if RegMethod=='rigid':
                    RegisteredSeriesNo = ShiftedSeriesNo + '11'
                elif RegMethod=='affine':
                    RegisteredSeriesNo = ShiftedSeriesNo + '12'
                elif RegMethod=='non-rigid':
                    RegisteredSeriesNo = ShiftedSeriesNo + '13'
                
                RegisteredSeriesDesc = ShiftedSeriesDesc + ' ' + RegMethod + '-reg'
                
            else:
                # Add 0 (to represent that Shift=0) after ToSeriesNo, then 1/2/3 
                # for RegMethod:
                if RegMethod=='rigid':
                    RegisteredSeriesNo = ToSeriesNo + '01'
                elif RegMethod=='affine':
                    RegisteredSeriesNo = ToSeriesNo + '02'
                elif RegMethod=='non-rigid':
                    RegisteredSeriesNo = ToSeriesNo + '03'
                
                RegisteredSeriesDesc = ToSeriesDesc + ' - ' + RegMethod + '-reg'
            
            
            
        """ April 6:  Simplify above: """
        # Let RegisteredSeriesNo be = MovingSeriesNo + '0' + N + '0' + 
        # FixedSeriesNo, 
        # where N = 1/2/3 for RegMethod = 'fixed'/'affine'/'non-rigid':
        if RegMethod=='rigid':
            N = '010'
            #N = '333'
        elif RegMethod=='affine':
            N = '020'
            #N = '444'
        elif RegMethod=='non-rigid':
            N = '030'
            #N = '555'
            
        RegisteredSeriesNo = MovingSeriesNo + N + FixedSeriesNo
        
        RegisteredSeriesDesc = MovingSeriesDesc + ' ' + RegMethod \
                               + ' registered to ' + FixedSeriesDesc
        
            
        # Create the registration directories:
        RegisteredDirName = RegisteredSeriesNo + ' ' + RegisteredSeriesDesc
        RegisteredDicomDir = os.path.join(RootDir, RegisteredDirName + ' DICOMs')
        RegisteredRoiDir = os.path.join(RootDir, RegisteredDirName + ' ROIs')
        
        # Add to AllDirs:
        AllDirs.append(RegisteredRoiDir)
        AllDirs.append(RegisteredDicomDir)
        
        # Update DataDict:      
        DataDict.update({'FixedKey':FixedKey, \
                         'MovingKey':MovingKey, \
                         'RegisteredKey':RegisteredKey, \
                         'Fixed':[], \
                         'Moving':[]
                         })    
               
        DataDict.update({RegisteredKey:{#'ScanLabel':MapPosCorrScanLabel,\
                                        'RoiDir':RegisteredRoiDir,\
                                        'DicomDir':RegisteredDicomDir,\
                                        'SeriesNo':RegisteredSeriesNo,\
                                        'SeriesDesc':RegisteredSeriesDesc,\
                                        'SliceNos':MovingSliceNos # make ='Moving'
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