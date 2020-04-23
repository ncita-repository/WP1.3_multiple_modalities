# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:52:44 2020

@author: ctorti
"""


"""
Function:
    GetDicoms()
    
Purpose:
    Get DICOM objects and the filepaths of the DICOM files.

Input:
    DirPath    - path to directory containing DICOM files
    
    SortMethod - method to use for sorting of files using GetDicomFpaths() with 
                acceptable values:
                -- 'none' (or 'None' or any other strings) = no sorting
                -- 'natural' (or 'Natural') = using natsort
                -- 'slices' (or 'Slices') = sorted along the slice direction 
    
Returns:
    fpaths    - full filepaths of DICOM files
                
    dicoms    - array of DICOM objects
"""


def GetDicoms(DirPath, SortMethod):
    
    # Import packages:
    import pydicom
    import importlib
    #from DicomHelperFuncs import GetDicomFpaths
    import GetDicomFpaths
    importlib.reload(GetDicomFpaths)
    from GetDicomFpaths import GetDicomFpaths
    
    #print('\nSortMethod (input for GetDicoms()) =', SortMethod)
    
    # Get the filepaths of the DICOM files:
    fpaths = GetDicomFpaths(DirPath, SortMethod)
    
    # Initialise a dictionary to store everything:
    #dicomDict = {}
    
    # Initialise a list to store the DICOM objects:
    dicoms = []
    
    # Loop though each DICOM filepath:
    for f in range(len(fpaths)):
        fpath = fpaths[f]
        
        #print('\nfpath =', fpath)
        
        # Read in the DICOM:
        dicom = pydicom.read_file(fpath)
        
        # Append this DICOM object:
        dicoms.append(dicom)
        
        # Get the SOP Instance UID:
        #uid = dicom.SOPInstanceUID
        
        # Get the Image Position (Patient):
        #position = [float(item) for item in dicom.ImagePositionPatient]
        
        # Get the Image Orientation Patient:
        #orientation = [float(item) for item in dicom.ImageOrientationPatient]
        
        # Get the Slice Location:
        #sliceLocation = float(dicom.SliceLocation)
        
        # Get the Pixel Spacing:
        #pixSpacing = [float(item) for item in dicom.PixelSpacing]
        
        # Get the Pixel Array:
        #pixelArray = dicom.pixel_array
    
        # Store the values in the dictionary using the SOP Instance UID as keys:
        #dicomDict[uid] = {'Dicom slice no':f+1, \
        #                  'Dicom fpath':fpath, \
        #                  'Image pos':position, \
        #                  'Image orient':orientation, \
        #                  'Slice loc':sliceLocation, \
        #                  'Pix spacing':pixSpacing, \
        #                  'Pix array':pixelArray, \
        #                  'Dicom':dicom
        #                 }
    
    #return dicomDict
    return fpaths, dicoms