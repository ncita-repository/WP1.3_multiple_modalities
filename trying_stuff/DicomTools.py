# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:44:06 2020

@author: ctorti
"""



# Import packages:
import os
from pydicom import read_file
import numpy as np
from natsort import natsorted

from ImageTools import GetImageAttributes



"""
******************************************************************************
******************************************************************************
GENERAL DICOM FUNCTIONS
******************************************************************************
******************************************************************************
"""


def GetDicomFpaths(DicomDir, SortMethod='slices', LogToConsole=None):
    """
    Get a sorted list of full filepaths for all DICOM files in a directory
    
    Input:
        DicomDir   - (String) Path to directory containing DICOMs
        
        SortMethod - (String) Method to use for sorting of files;
                    Acceptable values are:
                    -- 'none' (or 'None' or any other strings) = no sorting
                    -- 'natural' (or 'Natural') = using natsort
                    -- 'slices' (or 'Slices') = sorted along the slice 
                    direction in ascending order
        
    Returns:
        FilePaths - full filepaths of all Dicom files in DicomDir either:
                    -- without sorting
                    -- sorted in a natural way, 
                    e.g. 3-1, 3-2, ..., 3-9, 3-10, 3-11, ...
                    rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...
                    -- sorted along the scan direction in ascending order
                    
    """
    
    #print('\nSortMethod (input for GetDicomFpaths()) =', SortMethod)
    
    # Create an empty list of all Dicom file paths:
    FilePaths = []

    # Use os to get list of file paths of Dicom-only files:
    for DirName, SubDirList, FileList in os.walk(DicomDir):
        for FileName in FileList:
            if '.dcm' in FileName.lower():  # check for Dicom files only
                FilePaths.append(os.path.join(DirName,FileName))
                
    #print('FilePaths =', FilePaths)
    
    if FilePaths==[]:
        print('No DICOMs found.')
    else:
        if SortMethod in ['natural', 'Natural']:
            # Sort files in natural way:
            FilePaths = natsorted(FilePaths)
            
        if SortMethod in ['slices', 'Slices']:
            # Create an empty array to store the Positions:
            Positions = []
            
            for FilePath in FilePaths:
                #print('\nFilePath =', FilePath)
                
                # Read in Dicom:
                Dicom = read_file(FilePath)
                
                # Parse the Image Position (Patient):
                Position = [float(item) for item in Dicom.ImagePositionPatient]
                
                # Add Position for this Dicom to Positions:
                Positions.append(Position)
                
            # Convert to numpy array:
            Positions = np.array(Positions)
            
            #print('Positions =', Positions)
    
            # Get the maximum change of Positions along x, y and z:
            dPositions = np.max(Positions, axis=0) - np.min(Positions, axis=0)
            
            if LogToConsole:
                print('\n   The maximum change of Positions along x, y and z',
                      f'= {dPositions}')
    
            # Get the index of the axis with the greatest change of Position:
            ind = np.where(dPositions==np.max(dPositions))[0][0]
            
            if LogToConsole:
                print('\n   The index of the axis with the greatest change of',
                      f'Position = {ind}')
            
            # Convert ind to 'x', 'y' or 'z':
            if ind == 0:
                SortAxis = 'x'
            if ind == 1:
                SortAxis = 'y'
            if ind == 2:
                SortAxis = 'z'
                
            if LogToConsole:
                print('\n   The DICOM filepaths will be sorted along the ' \
                      + SortAxis + '-axis.')
            
                print('\n   The pre-sorted ' + SortAxis + '-positions are:\n')
                print([round(Positions[i][ind], 2) for i in range(len(Positions))])
            
            """
            Solution on how to sort along a specific column:
            
            https://stackoverflow.com/questions/22698687/how-to-sort-2d-array-numpy-ndarray-based-to-the-second-column-in-python
            
            The following are equivalent:
            
            Positions[Positions[:, 2].argsort()]
            Positions[np.argsort(Positions[:, 1])]
            
            """
            
            SortInds = Positions[:, ind].argsort()
            
            #PositionsSorted = Positions[SortInds]
            
            if LogToConsole:
                print('\n   The sorted ' + SortAxis + '-positions are:\n')
                print([round(Positions[i][ind], 2) for i in SortInds])
            
            # Sort FilePaths with SortInds:
            #FilePaths = FilePaths[SortInds] # <-- TypeError: list indices must be 
            # integers or slices, not list
            FilePaths = [FilePaths[i] for i in SortInds]
        
    
    return FilePaths





def ImportDicoms(DicomDir, SortMethod='slices', LogToConsole=False):
    """
    Import DICOM objects from a directory containing DICOM files.
    
    Input:
        DicomDir   - (String) Path to directory containing DICOM files
        
        SortMethod - (String) Method to use for sorting of files using 
                    GetDicomFpaths() with acceptable values:
                    -- 'none' (or 'None' or any other strings) = no sorting
                    -- 'natural' (or 'Natural') = using natsort to sort the 
                    filenames
                    -- 'slices' (or 'Slices') = sorted along the slice direction 
        
    Returns:
        Dicoms    - List of DICOM objects
    """
    
    
    #print('\nSortMethod (input for GetDicoms()) =', SortMethod)
    
    # Get the filepaths of the DICOM files:
    fpaths = GetDicomFpaths(DicomDir, SortMethod, LogToConsole)
    
    # Initialise a dictionary to store everything:
    #dicomDict = {}
    
    # Initialise a list to store the DICOM objects:
    Dicoms = []
    
    # Loop though each DICOM filepath:
    for f in range(len(fpaths)):
        fpath = fpaths[f]
        
        #print('\nfpath =', fpath)
        
        # Read in the DICOM:
        dicom = read_file(fpath)
        
        # Append this DICOM object:
        Dicoms.append(dicom)
        
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
    #return fpaths, dicoms
    return Dicoms






def GetDicomSOPuids(DicomDir):
    """
    Get a list of SOP Instance UIDs for a DICOM series.
    
    Inputs:
        DicomDir - (String) Directory containing the DICOMs
        
    Returns:
        SOPuids  - List of strings of the SOP UIDs of the DICOMs 
    """
    
    Dicoms = ImportDicoms(DicomDir=DicomDir, SortMethod='slices', 
                          LogToConsole=False)
    
    SOPuids = [] 
    
    for dicom in Dicoms:
        SOPuids.append(dicom.SOPInstanceUID)
        
    return SOPuids






def InspectDicomStudyDir(StudyDir):
    # First get a list of DICOM scan directories
    ScanDirs = []
    ScanDirNames = []

    # Use os to get list of file paths of Dicom-only files:
    for DirName, SubDirList, FileList in os.walk(StudyDir):
        #print('DirName =', DirName)
        #print('SubDirList =', SubDirList)
        #print('FileList =', FileList)
        #print('\n\n')
        
        DicomsFound = 0
        
        for file in FileList:
            if '.dcm' in file:
                #print(FileList, 'contains at least 1 DICOM')
                DicomsFound += 1
                
        
        if DicomsFound:
            ScanDirs.append(DirName)
            
            # DirName is a directory containing DICOMs
            RootDir, ScanDir = os.path.split(DirName)
    
            ScanDirNames.append(ScanDir)
    
    
    # Sort the list of directory names using natsort:
    ScanDirNames = natsorted(ScanDirNames)
    
    """ 
    The DICOM scan directory names have the form:
        'N-Label'
        
    e.g.:
        '2-T1 SE SAG IPATX2-05296'
        
    where scan number = 2
    and scan label = 'T1 SE SAG IPATX2-05296'
    
    Get the scan numbers and labels separately.
    """
    ScanNums = []
    ScanLabels = []
    
    for ScanDirName in ScanDirNames:
        # Index of first '-':
        ind = ScanDirName.index('-')
        
        ScanNums.append(ScanDirName[0:ind])
        
        ScanLabels.append(ScanDirName[ind+1::])
    
    
    # Initialise the Origins, Directions, Spacings, Dims, and FORuids:
    Origins = []
    Dirs = []
    Spacings = []
    Dims = []
    FORuids = []
    
    # Loop through each Scan directory:
    for ScanDir in ScanDirs:
    
        # Get the Image Attributes:
        origin, dirs,\
        spacings, dims = GetImageAttributes(DicomDir=ScanDir, 
                                            Package='pydicom')
        
        # Get the DICOMs:
        DicomFpaths, Dicoms = ImportDicoms(DirPath=ScanDir, 
                                           SortMethod='slices', 
                                           Debug=False)
        
        FORuid = Dicoms[0].FrameOfReferenceUID
        
        Origins.append(origin)
        Dirs.append(dirs)
        Spacings.append(spacings)
        Dims.append(dims)
        FORuids.append(FORuid)
        
    
    # Create a dictionary:
    StudyDict = {'Scan num':ScanNums,
                 'Scan label':ScanLabels,
                 #'FOR UID':FORuids,
                 'Origin':Origins,
                 'Directions':Dirs,
                 'Dimensions':Dims,
                 'Spacings':Spacings
                 }
    
    #return StudyDict, ScanNums, ScanLabels, FORuids, Origins, Dirs, Dims, Spacings
    return StudyDict
    
    



def InspectDicomSubjectDir(SubjectDir):
    """
    This is not complete.
    """
    
    # Get a list of study directories:
    StudyDirs = os.list(SubjectDir)
    
    # Sort using natsort:
    StudyDirs = natsorted(StudyDirs)
    
    
    """
    The directory names of the study directories have the form:
        'MM-DD-YYYY-MOD StudyLabel'
        
    where MOD is the modality (e.g. 'MRI', 'PET', 'CT')
    
    Get the date, modality and study label separately.
    """
    
    StudyDates = []
    Modalities = []
    StudyLabels = []
    
    for StudyDir in StudyDirs:
        # Get the study directory names:
        StudyDirRoot, StudyDirName = os.path.split(StudyDir)

        #StudyDates.append()
    

    return
