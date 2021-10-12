# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 15:49:45 2020

@author: ctorti
"""



"""
Function:
    GetSopUids()
    
Purpose:
    Get SOP Instance UID from a list of DICOM objects.

Input:
    Dicoms  - array of DICOM objects
    
    
Returns:
    SopUids - array of SOP Instance UIDs
                
"""


def GetSopUids(Dicoms):
    
    # Import packages:
    import pydicom
    
    # Initialise a list to store the SOP Instance UIDs:
    SopUids = []
    
    # Loop though each DICOM:
    for Dicom in Dicoms:
        # Get the SOP Instance UID:
        SopUid = Dicom.SOPInstanceUID
    
        # Store the SopUid:
        SopUids.append(SopUid)
        
    return SopUids