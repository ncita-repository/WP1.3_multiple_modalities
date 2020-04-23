# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:44:41 2020

@author: ctorti
"""



"""
Function:
    GetDicomFpaths()
    
Purpose:
    To get list of full file paths for all Dicom files in a given directory
    sorted using natsort.
    
    11/03/20 Note: Even though the Dicoms are sorted in a natural way (e.g.
    0000.dcm, 0001.dcm, 0002.dcm) it doesn't guarantee that the files are
    sorted in the sequence of the slices!  So rather than relying on the 
    file names I will need to use the Image Position (Patient) to sort. In the
    modifications to the function below sorting along the slice direction will 
    be an option. 
    

Input:
    DirPath    - path to directory containing Dicom files
    
    SortMethod - method to use for sorting of files with acceptable values:
                -- 'none' (or 'None' or any other strings) = no sorting
                -- 'natural' (or 'Natural') = using natsort
                -- 'slices' (or 'Slices') = sorted along the slice direction 
    
Returns:
    FilePaths - full file paths of all Dicom files in DirPath either:
                -- without sorting
                -- sorted in a natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ...
                rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
                -- sorted along the scan direction (The slice direction
                will be inferred by considering the min/max Positions along the 
                x, y and z directions.)
"""


def GetDicomFpaths(DirPath, SortMethod):
    
    # Import packages:
    import os
    
    #print('\nSortMethod (input for GetDicomFpaths()) =', SortMethod)
    
    # Create an empty list of all Dicom file paths:
    FilePaths = []

    # Use os to get list of file paths of Dicom-only files:
    for DirName, SubDirList, FileList in os.walk(DirPath):
        for FileName in FileList:
            if '.dcm' in FileName.lower():  # check for Dicom files only
                FilePaths.append(os.path.join(DirName,FileName))
    
    if SortMethod in ['natural', 'Natural']:
        # Import natsort:
        import natsort
        
        # Sort files in natural way:
        FilePaths = natsort.natsorted(FilePaths)
        
    if SortMethod in ['slices', 'Slices']:
        # Import packages:
        import pydicom
        import numpy as np
        
        # Create an empty array to store the Positions:
        Positions = []
        
        for FilePath in FilePaths:
            #print('\nFilePath =', FilePath)
            
            # Read in Dicom:
            Dicom = pydicom.read_file(FilePath)
            
            # Parse the Image Position (Patient):
            Position = [float(item) for item in Dicom.ImagePositionPatient]
            
            # Add Position for this Dicom to Positions:
            Positions.append(Position)
            
        # Convert to numpy array:
        Positions = np.array(Positions)

        # Get the maximum change of Positions along x, y and z:
        dPositions = np.max(Positions, axis=0) - np.min(Positions, axis=0)

        # Get the index of the axis with the greatest change of Position:
        ind = np.where(dPositions==np.max(dPositions))[0][0]
        
        # Convert ind to 'x', 'y' or 'z':
        if ind == 0:
            SortAxis = 'x'
        if ind == 1:
            SortAxis = 'y'
        if ind == 2:
            SortAxis = 'z'
        
        """
        Solution on how to sort along a specific column:
        
        https://stackoverflow.com/questions/22698687/how-to-sort-2d-array-numpy-ndarray-based-to-the-second-column-in-python
        
        The following are equivalent:
        
        Positions[Positions[:, 2].argsort()]
        Positions[np.argsort(Positions[:, 1])]
        
        """
        
        SortInds = Positions[:, ind].argsort()
        
        #PositionsSorted = Positions[SortInds]
        
        print('\n   The Dicom FilePaths will be sorted along the', SortAxis, \
              'positions.')
        
        # Sort FilePaths with SortInds:
        #FilePaths = FilePaths[SortInds] # <-- TypeError: list indices must be 
        # integers or slices, not list
        FilePaths = [FilePaths[i] for i in SortInds]
        
    
    return FilePaths