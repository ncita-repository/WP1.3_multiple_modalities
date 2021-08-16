# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:25:06 2020

@author: ctorti
"""



""" Function:
    CheckDicomPixValAndDtypeCompatibility()
    
Purpose:
    Check compatibility of an ROI Collection for a DICOM series

Input:
    DataDict - dictionary containing the DICOM and RTSTRUCT filepaths
    
    CheckKey - key that determine which values in DataDict to check

    
    
Returns:
    N/A
"""




def CheckDicomPixValAndDtypeCompatibility(FixedDir, MovingDir, RegDir):
    
    # Import packages:
    import pydicom
    import SimpleITK as sitk
    import numpy as np
    from GetDicomFpaths import GetDicomFpaths
    
    # Get the fixed, moving and registered DICOM filepaths:
    FixedFpaths = GetDicomFpaths(FixedDir, 'slices')
    MovingFpaths = GetDicomFpaths(MovingDir, 'slices')
    RegFpaths = GetDicomFpaths(RegDir, 'slices')
        
        
    # Read in the first DICOM in FixedDir, MovingDir and RegDir using both
    # Pydicom and SimpleITK:
    FixedPyd = pydicom.read_file(FixedFpaths[0])
    MovingPyd = pydicom.read_file(MovingFpaths[0])
    RegPyd = pydicom.read_file(RegFpaths[0])
    
    FixedItk = sitk.ReadImage(FixedFpaths[0])
    MovingItk = sitk.ReadImage(MovingFpaths[0])
    RegItk = sitk.ReadImage(RegFpaths[0])
    
        
    # Convert from SimpleITK to Numpy data arrays:
    FixedNda = sitk.GetArrayFromImage(FixedItk)
    MovingNda = sitk.GetArrayFromImage(MovingItk)
    RegNda = sitk.GetArrayFromImage(RegItk)
    


    # Combine all Pydicom, SimpleITK and Numpy data array DICOMs:
    AllPyds = [FixedPyd, MovingPyd, RegPyd]
    #AllItks = [FixedItk, MovingItk, RegItk]
    AllNdas = [FixedNda, MovingNda, RegNda]
    
    # Create an array of strings for the three different DICOMs:
    AllNames = ['Fixed', 'Moving', 'Registered']
    
    # Create an array of strings of the data types for all three DICOMs:
    AllDataTypes = []

    
    # Loop for each DICOM:
    for i in range(len(AllNames)):
        print(AllNames[i] + ':\n')
        
        # Get the Pydicom, SimpleITK and Nda:
        DicomNda = AllNdas[i]
        DicomPyd = AllPyds[i]
        #DicomItk = AllItks[i]
        
        # Count the number mismatches:
        errors = 0
        
        # Get datatype, min and max pixel values:
        DataType = str(DicomNda.dtype)
        MinValue = np.min(DicomNda)
        MaxValue = np.max(DicomNda)
        
        # Store the data type:
        AllDataTypes.append(DataType)
        
        # Define the minimum and maximum allowable values depending on the datatype:
        if DataType == 'uint8':
            MinAllowedValue = 0
            MaxAllowedValue = 255
        
        if DataType == 'uint16':
            MinAllowedValue = 0
            MaxAllowedValue = 65535
        
        if DataType == 'int8':
            MinAllowedValue = -128
            MaxAllowedValue = 127
        
        if DataType == 'int16':
            MinAllowedValue = -32768
            MaxAllowedValue = 32767
        
        print('Bits Allocated       =', DicomPyd.BitsAllocated)
        print('Bits Stored          =', DicomPyd.BitsStored)
        print('High Bit             =', DicomPyd.HighBit)
        print('Pixel Representation =', DicomPyd.PixelRepresentation, '(' + DataType + ')')
        print('Minimum pixel value  =', MinValue)
        print('Maximum pixel value  =', MaxValue)
        if 'int' in DataType:
            print('\nDtype Info           =', np.iinfo(DicomNda.dtype))
        if 'float' in DataType:
            print('\nDtype Info           =', np.finfo(DicomNda.dtype))
        
        # Check if MinValue < MinAllowedValue:
        if MinValue < MinAllowedValue:
            print(f'\nMinimum pixel value of {MinValue} does not match datatype', DataType)
            
            # Increment errors:
            errors = errors + 1
            
        # Check if MaxValue > MaxAllowedValue:
        if MaxValue > MaxAllowedValue:
            print(f'\nMaximum pixel value of {MaxValue} does not match datatype', DataType)
            
            # Increment errors:
            errors = errors + 1
            
        # Check if DataType is unsigned but MinValue < 0:
        if 'u' in DataType and MinValue < 0:
            print(f'\nMinimum pixel value of {MinValue} does not match datatype', DataType)
            
            # Increment errors:
            errors = errors + 1
            
        # Check if DataType is signed but MinValue not < 0:
        if not 'u' in DataType and not MinValue < 0:
            print(f'\nMinimum pixel value of {MinValue} does not match datatype', DataType)
            
            # Increment errors:
            errors = errors + 1

    
        # Print some of the pixel array:
        print('\nPixel Array:\n\n', DicomNda)
        
        # Report for this DICOM:
        if errors > 0:
            print(f'\n--> There were data type mismatches!!\n\n\n')
        else:
            print('\n--> OK.\n\n\n')
            
    