# -*- coding: utf-8 -*-
"""
Created on Mon May 11 07:36:34 2020

@author: ctorti
"""

"""
Function:
    GetXnatData()
    
Purpose:
    Download DICOMs and ROI Collection from XNAT prior to running CopyRoi().

Inputs:
    SourceRoiProjLabel     - Project Label of the source ROI Collection
    
    SourceRoiSubjLabel     - Subject Label of the source ROI Collection
    
    SourceRoiExpLabel      - Experiment Label of the source ROI Collection
    
    SourceRoiLabelDateTime - DateTime within the ROI Label of the source ROI 
                           Collection
    
    SourceSeriesNo         - Series number of the DICOM scan (collection) with
                           source ROI Collection
    
    SourceSeriesDesc       - Series description of the DICOM scan (collection) 
                           with source ROI Collection (Note: This is only used
                           for directory names of downloaded and new created
                           DICOM and ROI Collections)
    
    TargetDicomProjLabel   - Project Label of the target DICOM scan 
                             (collection) the ROI Collection is to be copied to
    
    TargetDicomSubjLabel   - Subject Label of the target DICOM scan 
                            (collection) the ROI Collection is to be copied to
    
    TargetDicomExpLabel    - Experiment Label of the target DICOM scan 
                            (collection) the ROI Collection is to be copied to
    
    TargetSeriesNo         - Series number of the target DICOM scan (collection) 
                            the ROI Collection is to be copied to
                            
    TargetSeriesDesc       - Series description of the target DICOM scan 
                           (collection) the ROI Collection is to be copied to
                           (Note: This is only used for directory names of 
                           downloaded and new created DICOM and ROI Collections)
    
    Debug                  - Print results if True

    
Returns:               
    DataDict               - Dictionary containing XNAT login details, REST path 
                            variables, directory paths, file paths, and other
                            details organised in keys (e.g. 'Source', 'Target',
                            'Shifted', 'Remapped', 'Registered')
                           - The dictionary is exported as a JSON file 
"""


    
def GetXnatData(XnatAddress, XnatUsername, XnatPassword,  
                RootDir, CopyFrom, CopyTo, 
                SourceRoiProjLabel, SourceRoiSubjLabel, SourceRoiExpLabel,
                SourceRoiLabelDateTime, SourceSeriesNo, SourceSeriesDesc,
                TargetDicomProjLabel, TargetDicomSubjLabel, TargetDicomExpLabel,
                TargetSeriesNo, TargetSeriesDesc,
                Debug):
    
    # IMPORT PACKAGES AND FUNCTIONS:
    import time, os#, shutil
    import importlib
    #import copy
    
    import CreateDir
    importlib.reload(CreateDir)
    from CreateDir import CreateDir
    
    import CreateDirsForDownloads
    importlib.reload(CreateDirsForDownloads)
    from CreateDirsForDownloads import CreateDirsForDownloads
    
    import DownloadFilesForRoiCopy
    importlib.reload(DownloadFilesForRoiCopy)
    from DownloadFilesForRoiCopy import DownloadFilesForRoiCopy
    
    import GetDicomFpaths
    importlib.reload(GetDicomFpaths)
    from GetDicomFpaths import GetDicomFpaths
    
    
    
    # Keep track of the time it takes to execute certain functions:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    # Get the current date:
    todaysDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    """
    Define the root directory (RootDir) to be a sub-directory of C:\Temp with
    todaysDate + SourceRoiProjLabel + SourceRoiSubjLabel + SourceRoiExpLabel 
    + ' copied to ' + TargetDicomExpLabel.
    """
    # Assign the root directory:
    #RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel \
    #                       + ' - ' + SourceRoiExpLabel + ' - ' + TargetDicomExpLabel)
    #RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel)
    RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel \
                           + ' ' + CopyFrom + ' to ' + CopyTo)
    
    
    # Create a dictionary to store various variables:
    #"""
    #Only run the following if DataDict doesn't already exist:
    #"""
    #if not 'DataDict' in locals(): 
    ##if not 'DataDict' in globals(): 
    
    DataDict = {# Define the XNAT login details:
                'XnatAddress':XnatAddress,
                'XnatUsername':XnatUsername,
                'XnatPassword':XnatPassword,
                
                # Add the other variables:
                'RootDir':RootDir,
                'Debug':Debug,
                
                
                # Define the REST variables that point to the ROIs to be 
                # copied and the DICOM files that the ROI applies to:
                'Source':{'RoiProjLabel':SourceRoiProjLabel,
                        'RoiSubjLabel':SourceRoiSubjLabel,
                        'RoiExpLabel':SourceRoiExpLabel,
                        'RoiLabelDateTime':SourceRoiLabelDateTime,
                        #'ScanLabel':SourceScanLabel
                        'SeriesNo':SourceSeriesNo,
                        'SeriesDesc':SourceSeriesDesc
                        },
        
                # Define the REST variables that point to the DICOM files  
                # that the ROI is to be applied to:
                'Target':{'DicomProjLabel':TargetDicomProjLabel,
                      'DicomSubjLabel':TargetDicomSubjLabel,
                      'DicomExpLabel':TargetDicomExpLabel,
                      #'ScanLabel':TargetScanLabel
                      'SeriesNo':TargetSeriesNo,
                      'SeriesDesc':TargetSeriesDesc
                      }
                }


    
  
    """ CREATE DIRECTORIES WHERE FILES WILL BE DOWNLOADED FROM XNAT: """   
    #DataDict = AddToDataDict(DataDict)
    DataDict = CreateDirsForDownloads(DataDict)
    
    
    
    """ DOWNLOAD THE DICOM-RTSTRUCT AND DICOM FILES FROM XNAT: """
    DataDict = DownloadFilesForRoiCopy(DataDict)

    
    # times[1] = time after downloading of files from XNAT:
    times.append(time.time())
    
    # The time taken to download the data from XNAT:
    Dtime = round(times[-1] - times[-2], 1)
    
    print(f'\nTook {Dtime} s to download data from XNAT.')
    
    return DataDict