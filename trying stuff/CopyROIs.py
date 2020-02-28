# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:08:42 2020

@author: ctorti
"""


""" Function:
    
    CopyROIs()
    
Purpose:
    Copy ROIs created for one series of DICOMs to a different series, so they 
    can be overlaid onto that different series.
    
Details:
    This function requires the REST variables of:
        
    i) The ROI to be copied and modified, i.e. project, subject, session 
    (experiment) and RTSTRUCT file name); and  
    
    ii) The DICOMs that the ROI is to be applied to, i.e. project, subject,
    session (experiment) and scan label.
    
    are required.
    
    The above files will be downloaded to a temporary folder, the data
    read in, ROI copied, modified and saved, then the temporary files will 
    be deleted.
    
    *** NOTES: ***
    1) It is assumed that there are equal scan numbers (equal file numbers) in 
    the two DICOM series.
    
Input:
    xnatDict - Dictionary containing XNAT login details, REST path variables 
               and directory paths, i.e.:
                   
                xnatDict = {# Define the XNAT login details:
                            'xnatAddress':'http://10.1.1.17',\
                            'username':'admin',\
                            'password':'password',\
                            # Define the REST variables that point to the ROIs 
                            # to be copied:
                            'RoiProjLab':'BrainTumorProg',\
                            'RoiSubjLab':'PGBM-002',\
                            'RoiExpLab':'PGBM-002_MR_1',\
                            'RoiLab':'RTSTRUCT_20200123_145731',\
                            # Define the REST variables that point to the DICOM 
                            # files that the ROI applies to:
                            'fromDicomProjLab':'BrainTumorProg',\
                            'fromDicomSubjLab':'PGBM-002',\
                            'fromDicomExpLab':'PGBM-002_MR_1',\
                            'fromDicomScanLab':'T1post',\
                            # Define the REST variables that point to the DICOM 
                            # files that the ROI is to be applied to:
                            'toDicomProjLab':'BrainTumorProg',\
                            'toDicomSubjLab':'PGM-002',\
                            'toDicomExpLab':'PGM-002_MR_1',\
                            'toDicomScanLab':'T1post',\
                            # Define the directories to download files to:
                            'fromDicomDir':r'C:\Temp\From DICOM',\
                            'fromRoiDir':r'C:\Temp\From RTSTRUCT',\
                            'toDicomDir':r'C:\Temp\To DICOM',\
                            'toRoiDir':r'C:\Temp\To RTSTRUCT'
                           }
    
    debug    - Print results if True
    
Returns:
    xnatDict - with the added full filepath of a new DICOM-RTSTRUCT file in the 
               directory specified by xnatDict['toRoiDir'], with a new file 
               name generated with the current date and time
"""


import xnat
#import pydicom
#import DicomHelperFuncs
#from DicomHelperFuncs import GetDicomFpaths
from ModifyROIs import ModifyROIs
#import copy
import os, shutil
import json


def CopyROIs(xnatDict, debug):
    if debug:
        print('The input variables (xnatDict):\n\n')
        print(xnatDict)
        print('\n\n')
    
    # Create XNAT session:
    # Note:  xnatDict['xnatAddress'] = 'http://10.1.1.17'
    session = xnat.connect(xnatDict['xnatAddress'], \
                           user=xnatDict['username'], \
                           password=xnatDict['password'])
    
    # Get the REST variables that point to the DICOM-RTSTRUCT file to be copied:
    RoiProjLab = xnatDict['RoiProjLab'] # e.g. 'BrainTumorProg'
    RoiSubjLab = xnatDict['RoiSubjLab'] # e.g. 'PGBM-002'
    RoiExpLab = xnatDict['RoiExpLab']   # e.g. 'PGBM-002_MR_1'
    RoiLab = xnatDict['RoiLab']         # e.g. 'RTSTRUCT_20200123_145731'
    
    # Define the REST variables that point to the DICOM files that the ROI applies to:  
    fromDicomProjLab = xnatDict['fromDicomProjLab'] # e.g. 'BrainTumorProg'
    fromDicomSubjLab = xnatDict['fromDicomSubjLab'] # e.g. 'PGBM-002'
    fromDicomExpLab = xnatDict['fromDicomExpLab']   # e.g. 'PGBM-002_MR_1'
    fromDicomScanLab = xnatDict['fromDicomScanLab'] # e.g. 'T1post'
    
    # Define the REST variables that point to the DICOM files that the ROI is to be applied to:
    toDicomProjLab = xnatDict['toDicomProjLab'] # e.g. 'BrainTumorProg'
    toDicomSubjLab = xnatDict['toDicomSubjLab'] # e.g. 'PGM-002'
    toDicomExpLab = xnatDict['toDicomExpLab']   # e.g. 'PGM-002_MR_1'
    toDicomScanLab = xnatDict['toDicomScanLab'] # e.g. 'T1post'
    
    # Get the directories to download files to:
    fromDicomDir = xnatDict['fromDicomDir'] # e.g. r'C:\Temp\From DICOM'
    fromRoiDir = xnatDict['fromRoiDir']     # e.g. r'C:\Temp\From RTSTRUCT' 
    toDicomDir = xnatDict['toDicomDir']     # e.g. r'C:\Temp\To DICOM'
    toRoiDir = xnatDict['toRoiDir']         # e.g. r'C:\Temp\To RTSTRUCT'  
    
    # Combine all directory paths:
    allDirs = [fromRoiDir, fromDicomDir, toRoiDir, toDicomDir]
    
    # Create directories if they don't exist:
    for dirPath in allDirs:
        try:
            os.mkdir(dirPath)
            if debug:
                print('Directory', dirPath, 'created')
        except FileExistsError:
            if debug:
                print('Directory', dirPath, 'already exists')
    
    """ Note:
    Using .download_dir(directory) downloads to a complicated XNAT-like folder structure,
    so instead will use .download(full_file_path) """
    
    # Get the DICOM-RTSTRUCT object:
    roi = session.projects[RoiProjLab].subjects[RoiSubjLab].experiments[RoiExpLab].assessors[RoiLab].resources['RTSTRUCT'].files[0]
    
    # Get the filename:
    roiId = roi.id
    
    # The ROI filepath:
    fromRoiFpath = os.path.join(fromRoiDir, roiId)
    
    # Download the DICOM-RTSTRUCT file:
    print('\nDownloading the DICOM-RTSTRUCT file to be copied:', roiId)
    roi.download(fromRoiFpath)
    
    
    # Get the DICOM objects of the series the ROIs relate to:
    #dicoms = session.projects[DicomProjLab].subjects[DicomSubjLab].experiments[DicomExpLab].scans[DicomScanLab].files
    dicoms = session.projects[fromDicomProjLab].subjects[fromDicomSubjLab].experiments[fromDicomExpLab].scans[fromDicomScanLab]
    
    # Iterate for each file in dicomObjects:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file relate to..')
    for f in range(len(dicoms.files)):
        # Get the filename:
        dicomId = dicoms.files[f].id
        
        # Download the DICOM:
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        dicoms.files[f].download(os.path.join(fromDicomDir, dicomId))
        
        
    # Get the DICOM objects of the series the ROIs are to be overlaid to:
    #dicoms = session.projects[DicomProjLab].subjects[DicomSubjLab].experiments[DicomExpLab].scans[DicomScanLab].files
    dicoms = session.projects[toDicomProjLab].subjects[toDicomSubjLab].experiments[toDicomExpLab].scans[toDicomScanLab]
    
    # Download each DICOM by iterating for each file in dicomObjects:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file is to be overlaid on..')
    for f in range(len(dicoms.files)):
        # Get the filename:
        dicomId = dicoms.files[f].id
        
        # Download the DICOM:
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        dicoms.files[f].download(os.path.join(toDicomDir, dicomId))
        
    
    # Use ModifyROIs to copy the contour data, modify and save to toRoiDir:
    toRoiFpath = ModifyROIs(fromDicomDir, fromRoiFpath, toDicomDir, toRoiDir)
    
    print('New DICOM-RTSTRUCT file saved to:\n\n', toRoiFpath, '\n\n')
    
    # Get the filename from toRoiFpath:
    toRoiFname = os.path.split(toRoiFpath)[1]
    
    
    # The dictionary will be updated (below) and exported to a JSON file in 
    # the path defined by toRoiDir.
    # Create a filename. Use the same filename as the RTSTRUCT file but change
    # the file extension:
    jsonFname = os.path.splitext(toRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    jsonFpath = os.path.join(toRoiDir, jsonFname)
    
    
    # Update the dictionary:
    #xnatDict['toRoiFname'] = toRoiFname
    #xnatDict['toRoiFpath'] = toRoiFpath
    xnatDict.update({'toRoiFname':toRoiFname})
    xnatDict.update({'toRoiFpath':toRoiFpath})
    xnatDict.update({'jsonFname':jsonFname})
    xnatDict.update({'jsonFpath':jsonFpath})
    
    
    # Export the updated dictionary to a JSON file:
    json.dump(xnatDict, open(jsonFpath, 'w'))
    
    print('JSON file created:\n', jsonFpath, '\n\n')
    
    
    # Do some housekeeping - delete the DICOM files in fromDicomDir and 
    # toDicomDir:
    print('Deleting temporary DICOM files..')
    for directory in [fromDicomDir, toDicomDir]:
        for root, dirs, files in os.walk(directory):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                shutil.rmtree(os.path.join(root, d))
    
    
    return xnatDict
    