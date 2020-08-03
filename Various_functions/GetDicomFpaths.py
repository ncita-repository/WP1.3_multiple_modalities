# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:44:41 2020

@author: ctorti
"""



"""
Function:
    GetDicomFpaths()
    
Purpose:
    To get list of full file paths for all DICOM files in a given directory
    sorted using natsort.
    
    11/03/20 Note: Even though the DICOMs are sorted in a natural way (e.g.
    0000.dcm, 0001.dcm, 0002.dcm) it doesn't guarantee that the files are
    sorted in the sequence of the slices!  So rather than relying on the 
    file names I will need to use the Image Position (Patient) to sort. In the
    modifications to the function below sorting along the slice direction will 
    be an option. 
    

Input:
    dirPath    - path to directory containing DICOM files
    
    sortMethod - method to use for sorting of files with acceptable values:
                -- 'none' (or 'None' or any other strings) = no sorting
                -- 'natural' (or 'Natural') = using natsort
                -- 'slices' (or 'Slices') = sorted along the slice direction 
    
Returns:
    filePaths - full file paths of all DICOM files in dirPath either:
                -- without sorting
                -- sorted in a natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ...
                rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
                -- sorted along the scan direction (The slice direction
                will be inferred by considering the min/max positions along the 
                x, y and z directions.)
"""


def GetDicomFpaths(dirPath, sortMethod):
    
    # Import packages:
    import os
    
    # Create an empty list of all DICOM file paths:
    filePaths = []

    # Use os to get list of file paths of DICOM-only files:
    for dirName, subdirList, fileList in os.walk(dirPath):
        for fileName in fileList:
            if '.dcm' in fileName.lower():  # check for DICOM files only
                filePaths.append(os.path.join(dirName,fileName))
    
    if sortMethod in ['natural', 'Natural']:
        # Import natsort:
        import natsort
        
        # Sort files in natural way:
        filePaths = natsort.natsorted(filePaths)
        
    if sortMethod in ['slices', 'Slices']:
        # Import packages:
        import pydicom
        import numpy as np
        
        # Create an empty array to store the positions:
        positions = []
        
        for filePath in filePaths:
            #print('\nfilePath =', filePath)
            
            # Read in DICOM:
            dicom = pydicom.read_file(filePath)
            
            # Parse the Image Position (Patient):
            position = [float(item) for item in dicom.ImagePositionPatient]
            
            # Add position for this DICOM to positions:
            positions.append(position)
            
        # Convert to numpy array:
        positions = np.array(positions)

        # Get the maximum change of positions along x, y and z:
        dPositions = np.max(positions, axis=0) - np.min(positions, axis=0)

        # Get the index of the axis with the greatest change of position:
        ind = np.where(dPositions==np.max(dPositions))[0][0]
        
        """
        Solution on how to sort along a specific column:
        
        https://stackoverflow.com/questions/22698687/how-to-sort-2d-array-numpy-ndarray-based-to-the-second-column-in-python
        
        The following are equivalent:
        
        positions[positions[:, 2].argsort()]
        positions[np.argsort(positions[:, 1])]
        
        """
        
        sortInds = positions[:, ind].argsort()
        
        #positionsSorted = positions[sortInds]
        
        # Sort filepaths with sortInds:
        #filePaths = filePaths[sortInds] # <-- TypeError: list indices must be 
        # integers or slices, not list
        filePaths = [filePaths[i] for i in sortInds]
        
    
    return filePaths