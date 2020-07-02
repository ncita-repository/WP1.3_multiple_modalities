# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:40:08 2020

@author: ctorti
"""



def DownloadRoiFromXnat(ProjectLabel, SubjectLabel, ExperimentLabel,
                        DateTimeInLabel, DownloadDir):
    
    # Import packages and functions:
    import xnat
    import os
    import pydicom
    #import json
    #import importlib
    
    
    # Define the XNAT address, username and password:
    XnatAddress = 'http://10.1.1.17'
    XnatUsername = 'admin'
    XnatPassword = 'admin'
    
    
    if DownloadDir == 'cwd':
        DownloadDir = os.getcwd()

    
    # Create XNAT session:
    """ Note:  XNAT Address is http://10.1.1.17 """
    session = xnat.connect(XnatAddress, \
                           user=XnatUsername, \
                           password=XnatPassword)
    
    """ 
    Note:
    Using .download_dir(directory) downloads to a complicated XNAT-like folder 
    structure, so instead will use .download(full_file_path) 
    """
    
    # Get the DICOM-RTSTRUCT object:
    Project = session.projects[ProjectLabel]
    Subject = Project.subjects[SubjectLabel]
    Experiment = Subject.experiments[ExperimentLabel]

    """
    It seems that the OHIF-Viewer exports ROI Collections with labels of 
    the form "AIM_YYYYMMDD_HHMMSS", whilst copies that are made and 
    exported using pydicom have the form "RTSTRUCT_YYYYMMDD_HHMMSS".
    
    Rather than finding the reason why and forcing the new RTSTRUCT file
    to have the same label, I'll loop through all ROI IDs and chose the one 
    that contains the DateTime specified in RoiLabelDateTime.
    """
    
    # Get all assessors for this RoiExp:
    Assessors = Experiment.assessors
    
    N = len(Assessors)
    
    print('\nSearching for ROI Collection with', DateTimeInLabel, 
          'in the ROI ID..')
    match = False
    
    # loop through each assessor and look for the one whose filename 
    # containes the DateTime contained in xnatDict['SourceRoiLabDateTime']:
    for i in range(N):
        # This ROI object:
        ThisRoi = Assessors[i].resources['RTSTRUCT'].files[0]
        ThisRoiID = ThisRoi.id
        
        print('\nFound ROI Collection with ID:', ThisRoiID)
        
        if DateTimeInLabel in ThisRoiID:
            # The matching ROI file:
            RoiFile = Assessors[i].resources['RTSTRUCT'].files[0]
            
            match = True
    
    if match:
        # Get the filename:
        Fname = RoiFile.id
        
        # The ROI filepath:
        Fpath = os.path.join(DownloadDir, Fname)
        
        # Download the DICOM-RTSTRUCT file if it hasn't already been:
        if os.path.isfile(Fpath):
            print('\nThe DICOM-RTSTRUCT file has already been downloaded to:\n',
                  Fpath)
        else:
            print('\nDownloading the DICOM-RTSTRUCT file to be copied:\n', Fpath)
        
            RoiFile.download(Fpath)
            
            
        RoiObject = pydicom.dcmread(Fpath)
            
        return RoiObject
    
    else:
        print('\n\nNo matching ROI Collection with the DateTime provided:',
              DateTimeInLabel)
        
        return []
    