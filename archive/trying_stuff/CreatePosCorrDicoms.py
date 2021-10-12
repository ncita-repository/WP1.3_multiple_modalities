# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:05:36 2020

@author: ctorti
"""


"""
It probably makes more sense to do this all within GetClosestSlices()
rather than as a separate function...
"""


def CreatePosCorrDicoms(dicomDictPosCorr, xnatDict, debug):
    # Import packages:
    import os
    
    # Get some necessary info from xnatDict:
    searchBy = xnatDict['searchBy']
    #scanPosCorr = xnatDict['scanPosCorr']
    toPosCorrDicomDir = xnatDict['toPosCorrDicomDir']
    
    
    # Get keys in dicomDictPosCorr:
    oldUids = list(dicomDictPosCorr.keys())
    
    for oldUid in oldUids:
        
        # Get the filepath of the original DICOM:
        oldFpath = dicomDictPosCorr[oldUid]['Dicom fpath']
        
        # Read in the original DICOM:
        dicom = pydicom.read_file(oldFpath)
            
        # Get the amended Image Position from the dictionary:
        pos2 = dicomDictPosCorr[oldUid]['Image pos']
        
        # Convert the elements from float to strings:
        pos2 = [str(item) for item in pos2]
        
        # Generate a new SOP Instance UID:
        newUid = pydicom.uid.generate_uid()
        
        # Replace the oldKey with newKey:
        dicomDictPosCorr[newUid] = dicomDictPosCorr.pop(oldUid)
        
        # Amend the Image Position in the DICOM object:
        dicom.ImagePositionPatient = pos2
        
        # Amend the SOP Instance UID in the DICOM object:
        dicom.SOPInstanceUID = newUid
        
        # Get the old file name with extension:
        oldFname = os.path.basename(oldFpath)
            
        # Get filename (without extension):
        #oldFname = os.path.splitext(oldFname)
        
        # Create a new filename, adding some info about the position correction:
        newFname = oldFname[0] + '-' + searchBy + '-PosCorr.dcm'
        
        if debug:
            print('\n   Changing the file name from:\n', \
                  oldFname, '\nto:\n', newFname)
        
        # Create new filepath:
        newFpath = os.path.join(toPosCorrDicomDir, newFname)
        
        # Replace the oldFpath with newFpath:
        dicomDictPosCorr[newUid]['Dicom fpath'] = newFpath
        
        
        if debug:
            print('\n    Exporting the position corrected DICOM to:\n', \
                  newFpath)
        
        # Export the DICOM file:
        dicom.save_as(newFpath)
        
        return dicomDictPosCorr
    
    
