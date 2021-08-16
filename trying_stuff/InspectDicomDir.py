# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 12:55:32 2020

@author: ctorti
"""



# Import packages and functions:
#import pydicom
import os
from GetDicoms import GetDicoms
from GetImageAttributes import GetImageAttributes
import natsort
    
    


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
    ScanDirNames = natsort.natsorted(ScanDirNames)
    
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
        DicomFpaths, Dicoms = GetDicoms(DirPath=ScanDir, 
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
    # Get a list of study directories:
    StudyDirs = os.list(SubjectDir)
    
    # Sort using natsort:
    StudyDirs = natsort.natsorted(StudyDirs)
    
    
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