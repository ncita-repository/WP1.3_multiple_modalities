# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 08:37:20 2020

@author: ctorti
"""

"""
Function:
    GetScanPositions()
    
Purpose:
    Get scan positions from a list of DICOM objects.

Input:
    Dicoms        - array of DICOM objects
    
    ScanAxis      - the scan axis with acceptable values 'x', 'y' or 'z'
    
    
Returns:
    ScanPositions - array of scan positions as float data type
                
"""


def GetScanPositions(Dicoms, ScanAxis):
    
    # Import packages:
    import pydicom
    
    #print('\nScanAxis =', ScanAxis)
    
    # Initialise a list to store the SOP Instance UIDs:
    ScanPositions = []
    
    # Loop though each DICOM:
    for Dicom in Dicoms:
        # Get the Image Position Patient:
        Position = Dicom.ImagePositionPatient
        
        # Store the scan position from Position to ScanPositions according to
        # the dimension defined by ScanAxis:
        if ScanAxis == 'x':
            ScanPositions.append(float(Position[0]))
            
        elif ScanAxis == 'y':
            ScanPositions.append(float(Position[1]))
            
        elif ScanAxis == 'z':
            ScanPositions.append(float(Position[2]))
            
        else:
            print('\nThe input argument "ScanAxis" must be "x", "y" or "z"')
    
        
    return ScanPositions