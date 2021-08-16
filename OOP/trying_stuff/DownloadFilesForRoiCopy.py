# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:45:38 2020

@author: ctorti
"""


""" Function:
    
    DownloadFilesForRoiCopy()
    
Purpose:
    Download from XNAT:
        1) The ROIs from a DICOM-RTSTRUCT file (to be copied), 
        2) All DICOMs from the same scan as the ROIs
        3) All DICOMs from a different scan for which the ROIs are to be
        overlaid onto 

    
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
    DataDict - Dictionary containing XNAT login details, REST path variables 
               and directory paths, etc.
    
    
Returns:
    DataDict - Dictionary as described above but with added entires
"""


def DownloadFilesForRoiCopy(DataDict):
    
    # Import packages and functions:
    import xnat
    import os
    #import json
    #import importlib
    
    
    # Get some necessary info from DataDict:
    Debug = DataDict['Debug']
    
    # Get the directories to download files to:
    #rootDir = DataDict['rootDir']           # e.g. r'C:\Temp\2020-02-28'
    SourceDicomDir = DataDict['Source']['DicomDir'] # e.g. r'C:\Temp\2020-02-28\From DICOM'
    SourceRoiDir = DataDict['Source']['RoiDir']     # e.g. r'C:\Temp\2020-02-28\From RTSTRUCT' 
    TargetDicomDir = DataDict['Target']['DicomDir']     # e.g. r'C:\Temp\2020-02-28\To DICOM'
    #toRoiDir = DataDict['To']['RoiDir']         # e.g. r'C:\Temp\2020-02-28\To RTSTRUCT' 

    if Debug:
        print('The input variables (xnatDict):\n\n')
        print(DataDict)
        print('\n\n')
    
    # Create XNAT session:
    """ Note:  XNAT Address is http://10.1.1.17 """
    session = xnat.connect(DataDict['XnatAddress'], \
                           user=DataDict['XnatUsername'], \
                           password=DataDict['XnatPassword'])
    
     
    
    
    """ 
    Note:
    Using .download_dir(directory) downloads to a complicated XNAT-like folder 
    structure, so instead will use .download(full_file_path) 
    """
    
    # Get the DICOM-RTSTRUCT object:
    RoiProj = session.projects[DataDict['Source']['RoiProjLabel']]
    RoiSubj = RoiProj.subjects[DataDict['Source']['RoiSubjLabel']]
    RoiExp = RoiSubj.experiments[DataDict['Source']['RoiExpLabel']]

    """
    It seems that the OHIF-Viewer exports ROI Collections with labels of 
    the form "AIM_YYYYMMDD_HHMMSS", whilst copies that are made and 
    exported using pydicom have the form "RTSTRUCT_YYYYMMDD_HHMMSS".
    
    Rather than finding the reason why and forcing the new RTSTRUCT file
    to have the same label, I've changed the key in xnatDict from
    'SourceRoiLabel', which is a specific label, to
    'SourceRoiLabelDateTime', which is only the DateTime part of the label.
    Below I'll loop through all ROI IDs and chose the one that contains
    the DateTime specified in 'SourceRoiLabelDateTime'.
    """
    
    # Get all assessors for this RoiExp:
    assessors = RoiExp.assessors
    
    # loop through each assessor and look for the one whose filename 
    # containes the DateTime contained in xnatDict['SourceRoiLabDateTime']:
    for i in range(len(assessors)):
        # This ROI object:
        thisRoi = assessors[i].resources['RTSTRUCT'].files[0]
        thisRoiID = thisRoi.id
        
        if DataDict['Source']['RoiLabelDateTime'] in thisRoiID:
            # The matching ROI object:
            roi = assessors[i].resources['RTSTRUCT'].files[0]
    
    # Get the filename:
    SourceRoiFname = roi.id
    
    # The ROI filepath:
    SourceRoiFpath = os.path.join(SourceRoiDir, SourceRoiFname)
    
    # Download the DICOM-RTSTRUCT file if it hasn't already been:
    if os.path.isfile(SourceRoiFpath):
        print('\nThe DICOM-RTSTRUCT file has already been downloaded.')
    else:
        print('\nDownloading the DICOM-RTSTRUCT file to be copied:', SourceRoiFname)
    
        roi.download(SourceRoiFpath)
    
    
    
    # Get the scan data of the series the ROIs relate to:
    #FromScanData = RoiExp.scans[DataDict['From']['ScanLabel']]
    SourceScanData = RoiExp.scans[DataDict['Source']['SeriesNo']]
    
    #return FromScanData
    
    # Store the DICOM filenames and filepaths in SourceScanData in lists:
    #FromDicomFnames = []
    #FromDicomFpaths = []

    # Iterate for each file in SourceScanData:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file relate to..')
    for f in range(len(SourceScanData.files)):
        # Get the filename:
        SourceDicomFname = SourceScanData.files[f].id
        
        # Store the filename in FromDicomFnames:
        #FromDicomFnames.append(FromDicomFname)
        
        # Create the filepath:
        SourceDicomFpath = os.path.join(SourceDicomDir, SourceDicomFname)
        
        # Store the filepath in FromDicomFpaths:
        #FromDicomFpaths.append(FromDicomFpath)
        
        if Debug:
            print('\n     Downloading DICOM file', SourceDicomFname)
        
        # Download the DICOM (if it hasn't already been):
        if not os.path.isfile(SourceDicomFpath):
            SourceScanData.files[f].download(SourceDicomFpath)
        
        
    # Get the scan data of the series the ROIs are to be overlaid to:
    TargetDicomsProj = session.projects[DataDict['Target']['DicomProjLabel']]
    TargetDicomsSubj = TargetDicomsProj.subjects[DataDict['Target']['DicomSubjLabel']] 
    TargetDicomsExp = TargetDicomsSubj.experiments[DataDict['Target']['DicomExpLabel']]
    #ToScanData = ToDicomsExp.scans[DataDict['To']['ScanLabel']]
    TargetScanData = TargetDicomsExp.scans[DataDict['Target']['SeriesNo']]
    
    # Store the DICOM filenames and filepaths in ToScanData in lists:
    #ToDicomFnames = []
    #ToDicomFpaths = []
    
    # Download each DICOM by iterating for each file in ToScanData:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file is to be overlaid on..')
    for f in range(len(TargetScanData.files)):
        # Get the filename:
        TargetDicomFname = TargetScanData.files[f].id
        
        # Store the filename in ToDicomFnames:
        #ToDicomFnames.append(ToDicomFname)
        
        # Create the filepath:
        TargetDicomFpath = os.path.join(TargetDicomDir, TargetDicomFname)
        
        # Store the filepath in ToDicomFpaths:
        #ToDicomFpaths.append(ToDicomFpath)
        
        if Debug:
            print('\n     Downloading DICOM file', TargetDicomFname)
        
        # Download the DICOM (if it hasn't already been):
        if not os.path.isfile(TargetDicomFpath):
            TargetScanData.files[f].download(TargetDicomFpath)
        
    
    # Update the dictionary:
    """
    Filenames commented out since redundant (can get filenames from filepaths).
    And DICOM filepaths commented out since function GetDicomFpaths() is more
    useful as it can produce a sorted list.
    """
    #DataDict['From'].update({'RoiFname':FromRoiFname})
    DataDict['Source'].update({'RoiFpath':SourceRoiFpath})
    #DataDict['From'].update({'DicomFnames':FromDicomFnames})
    #DataDict['From'].update({'DicomFpaths':FromDicomFpaths})
    #DataDict['To'].update({'DicomFnames':ToDicomFnames})
    #DataDict['To'].update({'DicomFpaths':ToDicomFpaths})
    """ April 14: The two lines below commented out since they account for all
    files, including non DICOM ones: """
    #DataDict['Source'].update({'SliceNos':list(range(len(SourceScanData.files)))})
    #DataDict['Target'].update({'SliceNos':list(range(len(TargetScanData.files)))})
    
    return DataDict
    