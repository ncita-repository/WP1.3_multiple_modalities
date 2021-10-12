# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:46:56 2020

@author: ctorti
"""


"""
Function:
    CombineDicomAndRoiData()
    
Purpose:
    Combine the dictionaries of DICOM and ROI data.

Input:
    dicomDict - dictionary of DICOM data
    
    roiDict   - dictionary DICOM-RTSTRUCT (ROI Collection) data
    
Returns:
    dicomAndRoiDict - dictionary combining DICOM and DICOM-RTSRUCT data with 
                      SOP Instance UID (from dicomDict) or Referenced SOP 
                      Instance UID (from roiDict) as keys
"""

def CombineDicomAndRoiData(dicomDict, roiDict):
    # Import packages:
    import copy
    
    # Get a list of keys from dicomDict (Dkeys) and roiDict (Rkeys):
    Dkeys = list(dicomDict.keys())
    Rkeys = list(roiDict.keys())
    
    # dicomAndRoiDict will start off as a copy of dicomDict:
    dicomAndRoiDict = copy.deepcopy(dicomDict)
    
    # Cycle through each key in dicomDict:
    for Dkey in Dkeys:
        # Check if Dkey is in Rkeys:
        if Dkey in Rkeys:
            # There is contour data for this Dkey with index:
            #ind = Rkeys.index(Dkey) # index of matching key in roiDict
            
            # Update dicomAndRoiDict[Dkey] with roiDict[Rkey]:
            dicomAndRoiDict[Dkey].update(roiDict[Dkey]) # Note: Dkey=Ckey
    
        else:
            # There is no contour data for this Dkey, so fill in the 
            # missing values:
            dicomAndRoiDict[Dkey].update({'Roi fpath':'', \
                                          'No of contours':0, \
                                          'Contour nos':[], \
                                          'No of contour pts':0, \
                                          'Contour data':[]
                                          })
        
    return dicomAndRoiDict

