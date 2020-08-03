# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 11:53:54 2020

@author: ctorti
"""



def ModifyRoiDirs(DataDict):
    
    # Import packages:
    import os
    
    # Get necessary into from DataDict:
    RootDir = DataDict['RootDir']
    Debug = DataDict['Debug'] 
    
    # Get keys in DataDict:
    keys = list(DataDict.keys())
    
    if True:#Debug:
        print('\nkeys before running ModifyRoiDirs = ', keys)
    
    # Iterate through each key in DataDict:
    for key in keys:
        #print('\nkey = ', key)
    
        #print('\nDataDict[' + key + '] =', DataDict[key])
        
        # If DataDict[key] is a dictionary, try to get the key 'DicomDir':
        if type(DataDict[key])==dict:
            try:
                DicomDirName = os.path.basename(DataDict[key]['DicomDir'])
                OldRoiDirName = os.path.basename(DataDict[key]['RoiDir'])
                
                if True:#Debug:
                    print('\nDicomDirName =', DicomDirName)
                    print('OldRoiDirName =', OldRoiDirName)
    
                # The string to search for in DicomDirName:
                string = 'DICOMs'
    
                # The index of DicomDirName after the string is:
                ind = DicomDirName.index(string) + len(string)
    
                #print(ind)
                #print(len(DicomDirName))
    
                """
                If DicomDirName has more characters than the number of 
                characters up to and including the 'DICOMs', then DicomDirName 
                has additional characters that need to be added to RoiDirName.
                """
                if len(DicomDirName) > ind:
                    NewRoiDirName = OldRoiDirName + DicomDirName[ind::]
                    
                    if True:#Debug:
                        print('NewRoiDirName =', NewRoiDirName)
    
                    # The new Roi Directory is:
                    NewRoiDir = os.path.join(RootDir, NewRoiDirName)
                    
                    # Update DataDict:
                    DataDict[key].update({'RoiDir':NewRoiDir})
    
            except KeyError:
                #print('\n   This key doesn\'t have \'DicomDir\'')
                continue
    
    if True:#Debug:
        # Get keys in DataDict:
        keys = list(DataDict.keys())
        
        print('\nkeys after running ModifyRoiDirs = ', keys)
    
    return DataDict