# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:47:48 2020

@author: ctorti
"""


"""
Function:
    GetClosestSlices()
    
Purpose:
    Starting with two DICOM dictionaries dicomDict1 and dicomDict2, find the 
    SOP Instance UID, filepath, Image Position, Image Orientation, Slice
    Location, etc. of the slices in the 2nd series that are closest to those
    in the 1st series. Also calculate the positional differences along x, y,
    z and r.
    
Input:
    dicomDict1 - dictionary of DICOM data of the 1st series
    
    dicomDict2 - dictionary of DICOM data of the 2nd series
    
    dataDict   - dictionary containing several variables
                
    
Returns:
    dicomDict1_to_2 - dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
                      
    dicomDict2posCorr - dicomDict2 with correction applied to the positional
                        dimension given by searchBy (if scanPosCorr=0, it will
                        simply be a copy of dicomDict2)
"""

#def GetClosestSlices(dicomDict1, dicomDict2, searchBy, scanPosCorr, debug):
def GetClosestSlices(dicomDict1, dicomDict2, dataDict):
    
    # Import packages:
    import copy
    import pydicom
    import os
    
    # Get some necessary info from dataDict:
    searchBy = dataDict['SearchBy']
    scanPosCorr = dataDict['ScanPosCorr']
    toPosCorrDicomDir = dataDict['To']['PosCorrDicomDir']
    debug = dataDict['Debug']
    
    
    # Initialise dicomDict1_to_2 as a copy of dicomDict1:
    dicomDict1_to_2 = copy.deepcopy(dicomDict1)
    
    # Get the list of keys of the two DICOM dictionaries:
    uids1 = list(dicomDict1.keys())
    oldUids2 = list(dicomDict2.keys())
    

    # First make corrections to the Image Positions of a copy of dicomDict2
    dicomDict2posCorr = copy.deepcopy(dicomDict2)
    
    #  If scanPosCorr is non-zero, loop through each key in 
    # dicomDict2posCorr, amending the Image Position:
        
    if scanPosCorr:
        for oldUid2 in oldUids2:
            # For debugging purposes, store pos2 before and after the change:
            pos2_debug = []
            
            # Get the Image Position:
            pos2 = dicomDict2posCorr[oldUid2]['Image pos']
            
            #print('\npos2 was', pos2)
            pos2_debug.append(pos2)
            
            # Apply correction to scan position based on the chosen searchBy:
            if searchBy=='x':
                pos2[0] = pos2[0] + scanPosCorr
            if searchBy=='y':
                pos2[1] = pos2[1] + scanPosCorr
            if searchBy=='z':
                pos2[2] = pos2[2] + scanPosCorr
                
            #print('\npos2 is now', pos2)
            pos2_debug.append(pos2)
            
            # Apply the corrected position to the dictionary:
            dicomDict2posCorr[oldUid2]['Image pos'] = pos2
            
            if debug:
                print('\npos2 was', pos2_debug[0], ', and changed to',\
                      pos2_debug[1])
                
                
            """ 
            Added on 21/03:
            Create new DICOM objects containing the new Image Positions,
            create new SOP Instance UIDs and save to new files.
            """
            
            # Generate a new SOP Instance UID:
            newUid2 = pydicom.uid.generate_uid()
            
            # Replace the old UID with the new one:
            dicomDict2posCorr[newUid2] = dicomDict2posCorr.pop(oldUid2)
            
            # Get the filepath of the original DICOM:
            oldFpath2 = dicomDict2posCorr[newUid2]['Dicom fpath']
            
            # Read in the original DICOM:
            dicom2 = pydicom.read_file(oldFpath2)
            
            # Convert the elements from float to strings:
            pos2str = [str(item) for item in pos2]
            
            # Amend the Image Position in the DICOM object:
            dicom2.ImagePositionPatient = pos2str
            
            # Amend the SOP Instance UID in the DICOM object:
            dicom2.SOPInstanceUID = newUid2
            
            # Get the old file name with extension:
            oldFname2 = os.path.basename(oldFpath2)
                
            # Get filename (without extension):
            #oldFname2 = os.path.splitext(oldFname2)
            
            # Create a new filename, adding some info about the position 
            # correction:
            newFname2 = oldFname2[0] + '-' + searchBy + '-PosCorr.dcm'
            
            if debug:
                print('\n   Changing the file name from:\n', \
                      oldFname2, '\nto:\n', newFname2)
            
            # Create new filepath:
            newFpath2 = os.path.join(toPosCorrDicomDir, newFname2)
            
            # Replace the old fpath in the dictionary with the new one:
            dicomDict2posCorr[newUid2]['Dicom fpath'] = newFpath2
            
            
            if True: #debug:
                print('\n    Exporting the position corrected DICOM to:\n', \
                      newFpath2)
            
            # Export the DICOM file:
            dicom2.save_as(newFpath2)
    
    
    # For every key (SOP UID) in dicomDict1 find the *new* SOP UID of the 
    # nearest slice in dicomDict2posCorr.
    
    # Get the new SOP Instance UIDs of dicomDict2posCorr:
    newUids2 = list(dicomDict2posCorr.keys())

    # Store the indices of the keys of dicomDict2posCorr whose positions are 
    # closest to the positions of the slices in dicomDict1:
    indsX = []
    indsY = []
    indsZ = []
    indsR = []
    
    for uid1 in uids1:
        pos1 = dicomDict1[uid1]['Image pos']
        
        #print('\npos1 =', pos1, '\n')
        
        # Create a list of the positional differences between pos2 and every 
        # slice in dicomDict1 along x, y, z and the vector:
        dX = []
        dY = []
        dZ = []
        dR = []
        
        for newUid2 in newUids2:
            #pos2 = dicomDict2[key2]['Image pos'] # this changes dicomDict2!
            pos2 = dicomDict2posCorr[newUid2]['Image pos']
            
            #print('pos2 =', pos2)
            
            # Calculate the distance between pos1 and pos2 along x, y and z:
            dx = pos2[0] - pos1[0]
            dy = pos2[1] - pos1[1]
            dz = pos2[2] - pos1[2]
            
            #print('\ndz =', dz)
            
            dr = (dx**2 + dy**2 + dz**2)**.5
            
            dX.append(dx)
            dY.append(dy)
            dZ.append(dz)
            dR.append(dr)
            
            #print('dx, dy, dz, dr =', dx, dy, dz, dr)
            
        # Find the index of the absolution minimum positional difference along 
        # different directions:
        indX = [abs(dx) for dx in dX].index(min([abs(dx) for dx in dX]))
        indY = [abs(dy) for dy in dY].index(min([abs(dy) for dy in dY]))
        indZ = [abs(dz) for dz in dZ].index(min([abs(dz) for dz in dZ]))
        indR = [abs(dr) for dr in dR].index(min([abs(dr) for dr in dR]))
        
        # Store these indices:
        indsX.append(indX)
        indsY.append(indY)
        indsZ.append(indZ)
        indsR.append(indR)
        
        # Update dicomDict1_to_2 with the details of dicomDict2posCorr that 
        # most closesly matchs the position along the dimension indicated by
        # searchBy:
        if searchBy=='x':
            ind = copy.deepcopy(indX)
        if searchBy=='y':
            ind = copy.deepcopy(indY)
        if searchBy=='z':
            ind = copy.deepcopy(indZ)
        if searchBy=='r':
            ind = copy.deepcopy(indR)
            
        if False: #debug:    
            print('\nind =', ind)
            
        
        # Simply indexing of the keys for below:
        key2 = newUids2[ind]
        
        # Update the dictionary:
        dicomDict1_to_2[uid1].update({'Closest Dicom SOP UID':\
                                      key2, \
                                      'Closest Dicom slice no':\
                                      dicomDict2posCorr[key2]['Dicom slice no'],\
                                      'Closest Dicom fpath':\
                                      dicomDict2posCorr[key2]['Dicom fpath'],\
                                      'Closest Image pos':\
                                      dicomDict2posCorr[key2]['Image pos'],\
                                      'Closest Image orient':\
                                      dicomDict2posCorr[key2]['Image orient'],\
                                      'Closest Slice loc':\
                                      dicomDict2posCorr[key2]['Slice loc'],\
                                      'Closest Pix spacing':\
                                      dicomDict2posCorr[key2]['Pix spacing'],\
                                      'x2 - x1':dX[ind],
                                      'y2 - y1':dY[ind],
                                      'z2 - z1':dZ[ind],
                                      'r2 - r1':dR[ind], \
                                      'Closest Pix array':\
                                      dicomDict2posCorr[key2]['Pix array'],\
                                      'Closest Dicom':\
                                      dicomDict2posCorr[key2]['Dicom'],\
                                     })
            
    #return dicomDict1_to_2
    return dicomDict1_to_2, dicomDict2posCorr

