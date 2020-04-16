# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



"""
Function:
    GetDicomFpaths()
    
Purpose:
    To get list of full file paths for all DICOM files in a given directory
    sorted using natsort.

Input:
    dirPath   - path to directory containing DICOM files
    
Returns:
    filePaths - full file paths of all DICOM files in dirPath sorted in a 
                natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ... rather than
                3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
"""


def GetDicomFpaths(dirPath):
    
    # Import packages:
    import os
    import natsort
    
    filePaths = []  # create an empty list of all DICOM file paths

    # Use os to get list of file paths of DICOM-only files:
    for dirName, subdirList, fileList in os.walk(dirPath):
        for fileName in fileList:
            if '.dcm' in fileName.lower():  # check for DICOM files only
                filePaths.append(os.path.join(dirName,fileName))
                 
    # Sort files in natural way:
    filePaths = natsort.natsorted(filePaths)
    
    return filePaths



""" 
NEW VERSION:
    
    There's no need to define the Project Label, Subject Label, and Experiment
    Label for both the original ROI and the DICOMs that relate to it.  It's
    enough to use the REST variables for the ROI to point to the DICOMs.
    
    
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
    xnatDict - Dictionary containing XNAT login details, REST path variables 
               and directory paths, i.e.:
                   
                xnatDict = {# Define the XNAT login details:
                            'xnatAddress':'http://10.1.1.17',\
                            'username':'admin',\
                            'password':'password',\
                            
                            # Define the REST variables that point to the ROIs 
                            # to be copied:
                            'roiProjLab':'BrainTumorProg',\
                            'roiSubjLab':'PGBM-002',\
                            'roiExpLab':'PGBM-002_MR_1',\
                            'roiLab':'RTSTRUCT_20200123_145731',\
                            
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
    
    debug      - Print results if True
    
Returns:
    xnatDict - with the added full filepath of a new DICOM-RTSTRUCT file in the 
               directory specified by xnatDict['toRoiDir'], with a new file 
               name generated with the current date and time
"""


def DownloadFilesForRoiCopy(xnatDict, debug):
    
    # Import packages and functions:
    import xnat
    import os
    #import json
    #import importlib
    

    if debug:
        print('The input variables (xnatDict):\n\n')
        print(xnatDict)
        print('\n\n')
    
    # Create XNAT session:
    # Note:  xnatDict['xnatAddress'] = 'http://10.1.1.17'
    session = xnat.connect(xnatDict['xnatAddress'], \
                           user=xnatDict['username'], \
                           password=xnatDict['password'])
    
    # Get the directories to download files to:
    rootDir = xnatDict['rootDir']           # e.g. r'C:\Temp\2020-02-28'
    fromDicomDir = xnatDict['fromDicomDir'] # e.g. r'C:\Temp\2020-02-28\From DICOM'
    fromRoiDir = xnatDict['fromRoiDir']     # e.g. r'C:\Temp\2020-02-28\From RTSTRUCT' 
    toDicomDir = xnatDict['toDicomDir']     # e.g. r'C:\Temp\2020-02-28\To DICOM'
    toRoiDir = xnatDict['toRoiDir']         # e.g. r'C:\Temp\2020-02-28\To RTSTRUCT'  
    
    # Combine all directory paths:
    allDirs = [rootDir, fromRoiDir, fromDicomDir, toRoiDir, toDicomDir]
    
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
    #roi = session.projects[RoiProjLab].subjects[RoiSubjLab].experiments[RoiExpLab].assessors[RoiLab].resources['RTSTRUCT'].files[0]
    roiProj = session.projects[xnatDict['fromRoiProjLab']]
    roiSubj = roiProj.subjects[xnatDict['fromRoiSubjLab']]
    roiExp = roiSubj.experiments[xnatDict['fromRoiExpLab']]

    """
    It seems that the OHIF-Viewer exports ROI Collections with labels of 
    the form "AIM_YYYYMMDD_HHMMSS", whilst copies that are made and 
    exported using pydicom have the form "RTSTRUCT_YYYYMMDD_HHMMSS".
    
    Rather than finding the reason why and forcing the new RTSTRUCT file
    to have the same label, I've changed the key in xnatDict from
    'fromRoiLab', which is a specific label, to
    'fromRoiLabDateTime', which is only the DateTime part of the label.
    Below I'll loop through all ROI IDs and chose the one that contains
    the DateTime specified in 'fromRoiLabDateTime'.
    """
    
    # Get all assessors for this roiExp:
    assessors = roiExp.assessors
    
    # loop through each assessor and look for the one whose filename 
    # containes the DateTime contained in xnatDict['fromRoiLabDateTime']:
    for i in range(len(assessors)):
        # This ROI object:
        thisRoi = assessors[i].resources['RTSTRUCT'].files[0]
        thisRoiID = thisRoi.id
        
        if xnatDict['fromRoiLabDateTime'] in thisRoiID:
            # The matching ROI object:
            roi = assessors[i].resources['RTSTRUCT'].files[0]
    
    # Get the filename:
    fromRoiFname = roi.id
    
    # The ROI filepath:
    fromRoiFpath = os.path.join(fromRoiDir, fromRoiFname)
    
    # Download the DICOM-RTSTRUCT file:
    print('\nDownloading the DICOM-RTSTRUCT file to be copied:', fromRoiFname)
    
    roi.download(fromRoiFpath)
    
    
    # Get the DICOM objects of the series the ROIs relate to:
    dicoms = roiExp.scans[xnatDict['fromDicomScanLab']]
    
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
    #dicoms = session.projects[toDicomProjLab].subjects[toDicomSubjLab].experiments[toDicomExpLab].scans[toDicomScanLab]
    dicomsProj = session.projects[xnatDict['toDicomProjLab']]
    dicomsSubj = dicomsProj.subjects[xnatDict['toDicomSubjLab']] 
    dicomsExp = dicomsSubj.experiments[xnatDict['toDicomExpLab']]
    dicoms = dicomsExp.scans[xnatDict['toDicomScanLab']]
    
    # Download each DICOM by iterating for each file in dicomObjects:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file is to be overlaid on..')
    for f in range(len(dicoms.files)):
        # Get the filename:
        dicomId = dicoms.files[f].id
        
        # Download the DICOM:
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        dicoms.files[f].download(os.path.join(toDicomDir, dicomId))
        
    
    # Update the dictionary:
    xnatDict.update({'fromRoiFname':fromRoiFname})
    xnatDict.update({'fromRoiFpath':fromRoiFpath})
    
    return xnatDict
    



"""
Function:
    GetDicomData()
    
Purpose:
    Get Image Position (Patient), Image Orientation (Patient), Slice Location,
    Pixel Spacing and Pixel Array from a directory of DICOMs, and store values
    in a dictionary using the SOP Instance UID as keys.

Input:
    dirPath  - path to directory containing DICOM files
    
Returns:
    dicomDict - dictionary containing Image Position (Patient), Image 
                Orientation (Patient), Slice Location, Pixel Spacing and Pixel 
                Array with SOP Instance UID as keys
"""

def GetDicomData(dirPath):
    
    # Import packages:
    import pydicom
    from DicomHelperFuncs import GetDicomFpaths
    
    # Get the filepaths of the DICOM files:
    dicomsFpaths = GetDicomFpaths(dirPath)
    
    # Initialise dictionary to store everything:
    dicomDict = {}
    
    # Loop though each DICOM filepath:
    for f in range(len(dicomsFpaths)):
        fpath = dicomsFpaths[f]
        
        # Read in the DICOM:
        dicom = pydicom.read_file(fpath)
        
        # Get the SOP Instance UID:
        uid = dicom.SOPInstanceUID
        
        # Get the Image Position (Patient):
        position = [float(item) for item in dicom.ImagePositionPatient]
        
        # Get the Image Orientation Patient:
        orientation = [float(item) for item in dicom.ImageOrientationPatient]
        
        # Get the Slice Location:
        sliceLocation = float(dicom.SliceLocation)
        
        # Get the Pixel Spacing:
        pixSpacing = [float(item) for item in dicom.PixelSpacing]
        
        # Get the Pixel Array:
        pixelArray = dicom.pixel_array
    
        # Store the values in the dictionary using the SOP Instance UID as keys:
        dicomDict[uid] = {'Dicom frame no':f+1, \
                          'Dicom fpath':fpath, \
                          'Image pos':position, \
                          'Image orient':orientation, \
                          'Slice loc':sliceLocation, \
                          'Pix spacing':pixSpacing, \
                          'Pix array':pixelArray
                         }
    
    return dicomDict



"""
Function:
    GetRoiData()
    
Purpose:
    Get No of contours (e.g. 3), Contour nos (e.g. [1, 2, 3]), No of contour
    pts (e.g. [44, 45, 116]), and Contour Data from a ROI filepath, and store 
    values in a dictionary using the Reference SOP Instance UID as keys.

Input:
    filePath  - filepath to a DICOM-RTSTRUCT (ROI Collection) file
    
Returns:
    roiDict - dictionary containing No of contours, Contour nos, No of contour
              pts, and Contour Data with Referenced SOP Instance UID as keys
"""

def GetRoiData(filePath):
    
    # Import packages:
    import pydicom
    
    # Read in the DICOM-RTSTRUCT file:
    roi = pydicom.read_file(filePath)
    
    # Get the Contour Sequences:
    sequences = roi.ROIContourSequence[0].ContourSequence

    # Initialise dictionary to store everything:
    roiDict = {}
    
    # Loop for each contour sequence:
    for sequence in sequences:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
    
        # Get the keys (UIDs) in contoursDict so far:
        keys = list(roiDict.keys())
    
        # Check if this uid is in keys. If it isn't this is the first time
        # for this uid.  Initialise the values:
        if uid not in keys:
            Ncontours = 0
            contourNo = [] 
            NcontourPts = []
            contourData = []
    
        # Append the following values for this sequence:
        
        # - the number of contours in this sequence:
        Ncontours = Ncontours + 1 # increment
        
        # -the Contour Number:
        contourNo.append(int(sequence.ContourNumber))
        
        # - the Number Of Contour Points:
        NcontourPts.append(int(sequence.NumberOfContourPoints))
        
        # - the Contour Data:
        #contourData.append(np.array(sequence.ContourData).astype(np.float))
        contourData.append([float(item) for item in sequence.ContourData])
    
        # Store the values for this uid. Values will be over-written if there's more than
        # one contour for this uid, but the lists have been appended already:
        roiDict[uid] = {'Roi fpath':filePath, \
                        'No of contours':Ncontours, \
                        'Contour nos':contourNo, \
                        'No of contour pts':NcontourPts, \
                        'Contour data':contourData
                       }
    
    return roiDict


"""
Function:
    CombineDicomAndRoiData()
    
Purpose:
    Combine the dictionaries of DICOM and ROI data.

Input:
    dicomDict - dictionary of DICOM data
    
    roiDict   - dictionary DICOM-RTSTRUCT (ROI Collection) data
    
Returns:
    dicomAndRoiDict - dictionary combining DICOM and DICOM-RTSRUCT data with 
                      SOP Instance UID (from dicomDict) or Referenced SOP 
                      Instance UID (from roiDict) as keys
"""

def CombineDicomAndRoiData(dicomDict, roiDict):
    # Import packages:
    import copy
    
    # Get a list of keys from dicomDict (Dkeys) and roiDict (Rkeys):
    Dkeys = list(dicomDict.keys())
    Rkeys = list(roiDict.keys())
    
    # dicomAndRoiDict will start off as a copy of dicomDict:
    dicomAndRoiDict = copy.deepcopy(dicomDict)
    
    # Cycle through each key in dicomDict:
    for Dkey in Dkeys:
        # Check if Dkey is in Rkeys:
        if Dkey in Rkeys:
            # There is contour data for this Dkey with index:
            #ind = Rkeys.index(Dkey) # index of matching key in roiDict
            
            # Update dicomAndRoiDict[Dkey] with roiDict[Rkey]:
            dicomAndRoiDict[Dkey].update(roiDict[Dkey]) # Note: Dkey=Ckey
    
        else:
            # There is no contour data for this Dkey, so fill in the 
            # missing values:
            dicomAndRoiDict[Dkey].update({'Roi fpath':'', \
                                          'No of contours':0, \
                                          'Contour nos':[], \
                                          'No of contour pts':0, \
                                          'Contour data':[]
                                          })
        
    return dicomAndRoiDict





"""
Function:
    GetClosestDicomPositions()
    
Purpose:
    Starting with two DICOM dictionaries dicomDict1 and dicomDict2, find the 
    SOP Instance UID, filepath, Image Position, Image Orientation, Slice
    Location, etc. of the slices in the 2nd series that are closest to those
    in the 1st series. Also calculate the positional differences along x, y,
    z and r.
    
Input:
    dicomDict1 - dictionary of DICOM data of the 1st series
    
    dicomDict2 - dictionary of DICOM data of the 2nd series
    
    searchBy   - how to define the closest slice with possible values:
                -- 'x' 
                -- 'y'
                -- 'z'
                -- 'r'
                to find closest slice along x/y/z direction or r, where
                r = sqrt(dx^2 + dy^2 + dz^2)
    
Returns:
    dicomDict1_to_2 - dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
"""

def GetClosestDicomPositions(dicomDict1, dicomDict2, searchBy):
    
    # Import packages:
    import copy
    
    # Initialise dicomDict1_to_2 as a copy of dicomDict1:
    dicomDict1_to_2 = copy.deepcopy(dicomDict1)
    
    # For every key (SOP UID) in dicomDict1 find the SOP UID of the nearest 
    # slice in dicomDict2:
    keys1 = list(dicomDict1.keys())
    keys2 = list(dicomDict2.keys())
    
    # Store the indices of the keys of dicomDict2 whose positions are closest 
    # to the positions of the slices in dicomDict1 for both the z direction 
    # and overall:
    indsX = []
    indsY = []
    indsZ = []
    indsR = []
    
    for key1 in keys1:
        pos1 = dicomDict1[key1]['Image pos']
        
        #print('\npos1 =', pos1, '\n')
        
        # Create a list of the positional differences between pos2 and every 
        # slice in dicomDict1 along x, y, z and the vector:
        dX = []
        dY = []
        dZ = []
        dR = []
        
        for key2 in keys2:
            pos2 = dicomDict2[key2]['Image pos']
            
            #print('pos2 =', pos2)
            
            # Calculate the distance between pos1 and pos2 along x, y and z:
            dx = pos2[0] - pos1[0]
            dy = pos2[1] - pos1[1]
            dz = pos2[2] - pos1[2]
            
            dr = (dx**2 + dy**2 + dz**2)**.5
            
            dX.append(dx)
            dY.append(dy)
            dZ.append(dz)
            dR.append(dr)
            
            #print('dx, dy, dz, dr =', dx, dy, dz, dr)
            
        # Find the index of the absolution minimum positional difference along 
        # different directions:
        indX = [abs(dx) for dx in dX].index(min([abs(dx) for dx in dX]))
        indY = [abs(dy) for dy in dY].index(min([abs(dy) for dy in dY]))
        indZ = [abs(dz) for dz in dZ].index(min([abs(dz) for dz in dZ]))
        indR = [abs(dr) for dr in dR].index(min([abs(dr) for dr in dR]))
        
        # Store these indices:
        indsX.append(indX)
        indsY.append(indY)
        indsZ.append(indZ)
        indsR.append(indR)
        
        # Update dicomDict1_to_2 with the details of dicomDict2 that most 
        # closesly match in position based on the chosen searchBy:
        if searchBy=='x':
            ind = copy.deepcopy(indX)
        if searchBy=='y':
            ind = copy.deepcopy(indY)
        if searchBy=='z':
            ind = copy.deepcopy(indZ)
        if searchBy=='r':
            ind = copy.deepcopy(indR)
            
        dicomDict1_to_2[key1].update({'Closest Dicom SOP UID':\
                                      keys2[ind], \
                                      'Closest Dicom frame no':\
                                      dicomDict2[keys2[ind]]['Dicom frame no'],\
                                      'Closest Dicom fpath':\
                                      dicomDict2[keys2[ind]]['Dicom fpath'],\
                                      'Closest Image pos':\
                                      dicomDict2[keys2[ind]]['Image pos'],\
                                      'Closest Image orient':\
                                      dicomDict2[keys2[ind]]['Image orient'],\
                                      'Closest Slice loc':\
                                      dicomDict2[keys2[ind]]['Slice loc'],\
                                      'Closest Pix spacing':\
                                      dicomDict2[keys2[ind]]['Pix spacing'],\
                                      'x2 - x1':dX[ind],
                                      'y2 - y1':dY[ind],
                                      'z2 - z1':dZ[ind],
                                      'r2 - r1':dR[ind],
                                     })
            
    return dicomDict1_to_2




"""
Function:
    ModifyRoi()
    
Purpose:
    Starting with an ROI object roi1 (from an DICOM-RTSTRUCT file corresponding
    to the 1st series of DICOMs) and dicomDict1_to_2 (the mapping from the 1st
    to the 2nd DICOM scan), copy roi1 and modify it accordingly so that the 
    relevant metadata corresponds to the DICOMs in the 2nd scan, and so, can
    be overlaid on the 2nd DICOM scan.
    
Input:
    roiDict1        - Dictionary of ROI data for the 1st DICOM scan
    
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
                      
    xnatDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths 
                    (see DownloadFilesForRoiCopy() for details)
    
    
Returns:
    roi2           - ROI object that contains the contours in roi1 but 
                    modified for overlaying over the 2nd DICOM scan;
                   - A new DICOM-RTSTRUCT file will be saved to the location
                   specified
                   
    xnatDict       - The updated dictionary
                   - The dictionary is exported as a JSON file
    
"""

def ModifyRoi(xnatDict, dicomDict1_to_2):
    # Import packages:
    import copy
    import time
    import json
    import pydicom
    import os
    
    # Load the DICOM-RTSTRUCT file whose ROIs correspond to the 1st DICOM scan:
    #roi1 = pydicom.dcmread(roiDict1['Roi fpath'])
    roi1 = pydicom.dcmread(xnatDict['fromRoiFpath'])
    
    # Use roi1 as a template for roi2:
    roi2 = copy.deepcopy(roi1)
    
    # Get the keys in roiDict1 and dicomDict1_to_2:
    #Rkeys = list(roiDict1.keys())
    D1to2keys = list(dicomDict1_to_2.keys())
    
    # First make changes to roi2 that only need to be made once (don't need
    # to be made for each contour).
    
    # Read in the first DICOM from the 2nd scan:
    fpath2 = dicomDict1_to_2[D1to2keys[0]]['Closest Dicom fpath']
    
    # Load the DICOM:
    dicom2 = pydicom.read_file(fpath2)
    
    # Create a dictionary to store the changes made to the ROI that will 
    # later be added to xnatDict:
    changesDict = {}
    
    # Make changes to the metadata in roi2 to reflect the corresponding 
    # metadata contained in the DICOM:
    
    # - Change the Study Date:

    if roi2.StudyDate==dicom2.StudyDate:
        print('\nThe Study Date are the same:\n', roi2.StudyDate)
    else:
        print('\nChanging the Study Date from:\n', roi2.StudyDate, '\nto:\n', \
              dicom2.StudyDate)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyDate':roi2.StudyDate})
        
        roi2.StudyDate = copy.deepcopy(dicom2.StudyDate)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyDate':dicom2.StudyDate})

        
    # - Change the Study Time:
    """ Note: Some DICOM Study Times have a decimal with trailing zeros,
        e.g. '164104.000000', so convert to float, then integer, then back 
        to string. """
    
    newTime = str(int(float(copy.deepcopy(dicom2.StudyTime))))

    if roi2.StudyTime==newTime:
        print('\nThe Study Time are the same:\n', roi2.StudyTime)
    else:
        print('\nChanging the Study Time from:\n', roi2.StudyTime, '\nto:\n', \
              newTime)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyTime':roi2.StudyTime})
        
        roi2.StudyTime = newTime
        
        # Update the dictionary changesDict:
        changesDict.update({'toStudyTime':roi2.StudyTime})
            
    # - Change the Patient Name:
    if roi2.PatientName==dicom2.PatientName:
        print('\nThe Patient Name are the same:\n', roi2.PatientName)
    else:
        print('\nChanging the Patient Name from:\n', roi2.PatientName, \
              '\nto:\n', dicom2.PatientName)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromPatientName':roi2.PatientName})
        
        roi2.PatientName = copy.deepcopy(dicom2.PatientName)
        
        # Update the dictionary changesDict:
        changesDict.update({'toPatientName':roi2.PatientName})
        
    
    # - Change the Patient ID:
    if roi2.PatientID==dicom2.PatientID:
        print('\nThe Patient ID are the same:\n', roi2.PatientID)
    else:
        print('\nChanging the Patient ID from:\n', roi2.PatientID, '\nto:\n', \
              dicom2.PatientID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromPatientID':roi2.PatientID})
    
        roi2.PatientID = copy.deepcopy(dicom2.PatientID)
        
        # Update the dictionary changesDict:
        changesDict.update({'toPatientID':roi2.PatientID})
    
    # - Change the Patient Birth Date:
    """ Not all DICOMs have Patient's Birth Date, so.."""
    try:
        if roi2.PatientBirthDate==dicom2.PatientBirthDate:
            print('\nThe Patient Birth Date are the same:\n', \
                  roi2.PatientBirthDate)
        else:
            print('\nChanging the Patient Birth Date from:\n', \
                  roi2.PatientBirthDate, '\nto:', dicom2.PatientBirthDate)
            
            # Update the dictionary changesDict:
            changesDict.update({'fromPatientBirthDate':roi2.PatientBirthDate})
            
            roi2.PatientBirthDate = copy.deepcopy(dicom2.PatientBirthDate)
            
            # Update the dictionary changesDict:
            changesDict.update({'toPatientBirthDate':roi2.PatientBirthDate})
            
    except AttributeError:
        print('\nPatient Birth Date does not exist in the DICOMs')
     
    # - Change the Patient Sex:
    if roi2.PatientSex==dicom2.PatientSex:
        print('\nThe Patient Sex are the same:\n', roi2.PatientSex)
    else:
        print('\nChanging the Patient Sex from:\n', roi2.PatientSex, \
              '\nto:\n', dicom2.PatientSex)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromPatientSex':roi2.PatientSex})
        
        roi2.PatientSex = copy.deepcopy(dicom2.PatientSex)
        
        # Update the dictionary changesDict:
        changesDict.update({'toPatientSex':roi2.PatientSex})
    
    # - Change the Study Instance UID:
    if roi2.StudyInstanceUID==dicom2.StudyInstanceUID:
        print('\nThe Study Instance UID are the same:\n', \
              roi2.StudyInstanceUID)
    else:
        print('\nChanging the Study Instance UID from:\n', \
              roi2.StudyInstanceUID, '\nto:\n', dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromStudyInstanceUID':roi2.StudyInstanceUID})
        
        roi2.StudyInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'toStudyInstanceUID':roi2.StudyInstanceUID})
        
    # - Change the Study Instance UID in the Referenced Frame Of Reference
    # Sequence:
    roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0].ReferencedSOPInstanceUID     
              
    if roiSIuid==dicom2.StudyInstanceUID:
        print('\nThe Study Instance UID are the same:\n', roiSIuid)
    else:
        print('\nChanging the Referenced SOP Instance UID from:\n', \
              roiSIuid, '\nto:\n', dicom2.StudyInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.RTRSS.RefSOPInstanceUID':roiSIuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .ReferencedSOPInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
        
        roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
              .RTReferencedStudySequence[0].ReferencedSOPInstanceUID   
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.RTRSS.RefSOPInstanceUID':roiSIuid})
    
    # - Change the Series Instance UID:
    roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
               .SeriesInstanceUID
               
    if roiSIuid==dicom2.SeriesInstanceUID:
        print('\nThe Series Instance UID are the same:\n', roiSIuid)
    else:
        print('\nChanging the Series Instance UID from:\n', \
              roiSIuid, '\nto:\n', dicom2.SeriesInstanceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.RTRSS.RTRSS.SeriesInstanceUID':\
                            roiSIuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0]\
        .RTReferencedSeriesSequence[0].SeriesInstanceUID \
        = copy.deepcopy(dicom2.SeriesInstanceUID)
        
        roiSIuid = roi2.ReferencedFrameOfReferenceSequence[0]\
               .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
               .SeriesInstanceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.RTRSS.RTRSS.SeriesInstanceUID':roiSIuid})
    
    # - Change the Frame of Reference UID:
    if roi2.FrameOfReferenceUID==dicom2.FrameOfReferenceUID:
        print('\nThe Frame Of ReferenceUID UID are the same:\n', \
              roi2.FrameOfReferenceUID)
    else:
        print('\nChanging the Frame Of ReferenceUID UID from:\n', \
              roi2.FrameOfReferenceUID, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'changesToRoi':\
                         {'fromFrameOfReferenceUID':roi2.FrameOfReferenceUID}})
        
        roi2.FrameOfReferenceUID = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'toFrameOfReferenceUID':roi2.FrameOfReferenceUID})
    
    # - Change the Frame of Reference UID in the Referenced Frame Of Reference
    # Sequence:
    roiFRuid = roi2.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
    
    if roiFRuid==dicom2.FrameOfReferenceUID:
        print('\nThe Frame Of ReferenceUID UID are the same:\n', roiFRuid)
    else:
        print('\nChanging the Frame Of ReferenceUID UID from:\n', \
              roiFRuid, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromRFORS.FrameOfReferenceUID':roiFRuid})
        
        roi2.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID \
        = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        roiFRuid = roi2.ReferencedFrameOfReferenceSequence[0]\
        .FrameOfReferenceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toRFORS.FrameOfReferenceUID':roiFRuid})
    
    # - Change the Frame of Reference UID in the Structure Set ROI Sequence:
    roiFRuid = roi2.StructureSetROISequence[0].ReferencedFrameOfReferenceUID
    
    if roiFRuid==dicom2.FrameOfReferenceUID:
        print('\nThe Referenced Frame Of ReferenceUID UID is the same as', \
              'the Frame Of Reference UID:\n', roiFRuid)
    else:
        print('\nChanging the Referenced Frame Of ReferenceUID UID from:\n', \
              roiFRuid, '\nto:\n', dicom2.FrameOfReferenceUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'fromSSROIS.RefFrameOfReferenceUID':roiFRuid})
        
        roi2.StructureSetROISequence[0].ReferencedFrameOfReferenceUID \
        = copy.deepcopy(dicom2.FrameOfReferenceUID)
        
        roiFRuid = roi2.StructureSetROISequence[0]\
        .ReferencedFrameOfReferenceUID
        
        # Update the dictionary changesDict:
        changesDict.update({'toSSROIS.RefFrameOfReferenceUID':roiFRuid})
    
    # - Change the Structure Set label:
    
    #""" 
    #Assume that the characters before the first "_" is always the Patient Name. 
    #"""
    #roiValue = roi2.StructureSetLabel
    
    # Find the first "_":
    #ind = roiValue.index('_')
    
    #"""
    ## Modify the Structure Set Label using the Patient Name from the DICOM 
    ## followed by the characters including and following the first "_" in 
    ## roiValue:
    #newLabel = str(copy.deepcopy(dicom2.PatientName)) + roiValue[ind::]
    #"""
    
    # - Modify the Structure Set Label by adding 
    # "_copy_from_Scan_A_to_B_YYYYMMDD_HHMMSS" to the Structure Set Label of 
    # the original RTSTRUCT:
    """ Note: A and B are the Scan Labels from the two scans (collections)."""
    newDate = time.strftime("%Y%m%d", time.gmtime())
    newTime = time.strftime("%H%M%S", time.gmtime())
    
    #newLabel = roi1.StructureSetLabel + '_copy_' + newDate \
    #+ '_' + newTime
    
    newLabel = roi1.StructureSetLabel + '_copy_from_Scan_' \
    + xnatDict['fromDicomScanLab'] + '_to_' + xnatDict['toDicomScanLab'] \
    + '_' + newDate + '_' + newTime
    
    print('\nChanging the Structure Set Label from:\n', \
          roi2.StructureSetLabel, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromStructureSetLabel':roi2.StructureSetLabel})
    
    roi2.StructureSetLabel = newLabel  
    
    # Update the dictionary changesDict:
    changesDict.update({'toStructureSetLabel':roi2.StructureSetLabel})
    
    
    # - Change the ROI Name in the Structure Set ROI Sequence:
    """
    Assume that the characters before the first "_" is always the Patient Name.
    """
    roiValue = roi2.StructureSetROISequence[0].ROIName
    
    # Find the first "_":
    ind = roiValue.index('_')
    
    #"""
    #Modify the ROI Name using the Patient Name from the DICOM followed by the
    #characters including and following the first "_" in roiValue.
    #"""
    #newRoiName = str(copy.deepcopy(dicom2.PatientName)) + roiValue[ind::]
    
    """
    Modify the ROI Name by adding "Copy_of_" to the ROI Name of the original 
    RTSTRUCT and the current date and time.
    """
    #newDate = time.strftime("%Y%m%d", time.gmtime())
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    #newRoiName = 'Copy_of_' + roi1.StructureSetROISequence[0].ROIName \
    #+ newDate + '_' + newTime
    
    """
    Modify the ROI Name by making it the same as the new Structure Set Label.
    """
    
    print('\nChanging the ROI Name from:\n', \
          roi2.StructureSetROISequence[0].ROIName, '\nto:\n', newLabel)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromSSROIS.ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
    
    roi2.StructureSetROISequence[0].ROIName = newLabel
    
    # Update the dictionary changesDict:
    changesDict.update({'toSSROIS.ROIName':\
                     roi2.StructureSetROISequence[0].ROIName})
                    
    
    # Now make changes to roi2 that need to be made for each Contour Image
    # Sequence:
    # Get the Contour Image Sequences:
    sequences = roi2.ReferencedFrameOfReferenceSequence[0]\
        .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
        .ContourImageSequence

    #print('\nLooping through each Contour Image Sequence to:'\
    #      '\n- get the Referenced SOP Instance UID belonging to the 1st DICOM '\
    #      'scan (collection)'\
    #      '\n- find the index of the matching SOP Instance UID in '\
    #      'dicomDict1_to_2'\
    #      '\n- use the index to get the filepath of the DICOM slice in the '\
    #      '2nd scan (collection) that is closest in position to that of the '\
    #      '1st DICOM scan (collection)'\
    #      '\n- modify the (ROI\'s) Referenced SOP Instance UID to the '\
    #      'corresponding SOP Instance UID from the closest DICOM in the 2nd '\
    #      'scan (collection)')
    
    print('\nLooping through each Contour Image Sequence to modify the '\
          'Referenced SOP Instance UIDs..')
    
    
    # Store the original and modified Referenced SOP Instance UIDs from the
    # Referenced Frame Of Reference Sequence:
    fromRFORSrefsopuids = []
    toRFORSrefsopuids = []
    
    # Cycle through each Contour Image Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ReferencedSOPInstanceUID
        
        # Store this fromRFORSrefsopuid:
        fromRFORSrefsopuids.append(refsopuid)
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # The filepath of the slice from the 2nd DICOM scan that is closest in 
        # postion to the slice in the 1st DICOM scan with this SOP UID is:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
        
        # Load the DICOM from the 2nd scan:
        dicom2 = pydicom.read_file(fpath2)
        
        # Store this toRFORSrefsopuid:
        toRFORSrefsopuids.append(dicom2.SOPInstanceUID)
    
        if refsopuid==dicom2.SOPInstanceUID:
            print('\n- The DICOM\'s SOP Instance UID is the same as the ', \
                  'ROI\'s Referenced SOP Instance UID:\n', refsopuid)
        else:
            print('\n- Changing the Referenced SOP Instance UID from:\n', \
                  refsopuid, '\nto:\n', dicom2.SOPInstanceUID)
            
            roi2.ReferencedFrameOfReferenceSequence[0]\
            .RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\
            .ContourImageSequence[i].ReferencedSOPInstanceUID \
            = copy.deepcopy(dicom2.SOPInstanceUID)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromRFORS.RTRSS.RTRSS.CIS.RefSOPInstanceUIDs':\
                      fromRFORSrefsopuids})
    changesDict.update({'toRFORS.RTRSS.RTRSS.CIS.RefSOPInstanceUIDs':\
                      toRFORSrefsopuids})
                    
                    
    # Now make changes to roi2 that need to be made for each Contour Sequence:
    # Get the Contour Sequences:
    sequences = roi2.ROIContourSequence[0].ContourSequence

    #print('\nLooping through each Contour Sequence to:'\
    #      '\n- get the Referenced SOP Instance UID belonging to the 1st DICOM '\
    #      'scan (collection)'\
    #      '\n- find the index of the matching SOP Instance UID in '\
    #      'dicomDict1_to_2'\
    #      '\n- use the index to get the filepath of the DICOM slice in the '\
    #      '2nd scan (collection) that is closest in position to that of the '\
    #      '1st DICOM scan (collection)'\
    #      '\n- modify the (ROI\'s) Referenced SOP Instance UID to the '\
    #      'corresponding SOP Instance UID from the closest DICOM in the 2nd '\
    #      'scan (collection)')
    
    print('\nLooping through each Contour Sequence to modify the '\
          'Referenced SOP Instance UIDs..')
       
    # Store the original and modified Referenced SOP Instance UIDs from the
    # ROI Contour Sequence:
    fromROICSrefsopuids = []
    toROICSrefsopuids = []
    
    # Cycle through each Contour Sequence:
    for i in range(len(sequences)):
        # Get the Referenced SOP Instance UID for this sequence:
        refsopuid = sequences[i].ContourImageSequence[0]\
        .ReferencedSOPInstanceUID
        
        # Store this fromROICSrefsopuid:
        fromROICSrefsopuids.append(refsopuid)
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # The filepath of the slice from the 2nd DICOM scan that is closest in 
        # postion to the slice in the 1st DICOM scan with this SOP UID is:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
        
        # Load the DICOM from the 2nd scan:
        dicom2 = pydicom.read_file(fpath2)
        
        # Store this toROICSrefsopuid:
        toROICSrefsopuids.append(dicom2.SOPInstanceUID)
    
        if refsopuid==dicom2.SOPInstanceUID:
            print('\n- The DICOM\'s SOP Instance UID is the same as the',\
                  'ROI\'s Referenced SOP Instance UID:\n', refsopuid)
        else:
            print('\n- Changing the Referenced SOP Instance UID from:\n', \
                  refsopuid, '\nto:\n', dicom2.SOPInstanceUID)
            
            roi2.ROIContourSequence[0].ContourSequence[i]\
            .ContourImageSequence[0].ReferencedSOPInstanceUID \
            = copy.deepcopy(dicom2.SOPInstanceUID)    
        
    # Update the dictionary changesDict:
    changesDict.update({'fromROICS.CS.CIS.RefSOPInstanceUIDs':\
                      fromROICSrefsopuids})
    changesDict.update({'toROICS.CS.CIS.RefSOPInstanceUIDs':\
                      toROICSrefsopuids})
            
            
    # To help identify the ROI Collection from others update the Structure
    # Set Date and Structure Set Time.
    # - Update the Structure Set Date:
    #roi2.StudyDate = time.strftime("%Y%m%d", time.gmtime())
    #newDate = time.strftime("%Y%m%d", time.gmtime())
    
    if roi2.StructureSetDate==newDate:
        print('\nThe Structure Set Date is today:\n', newDate)
    else:
        print('\nChanging the Structure Set Date from:\n', \
              roi2.StructureSetDate, '\nto:\n', newDate)
        
        # Update the dictionary xnatDict:
        xnatDict.update({'changesToRoi':\
                         {'fromStructureSetDate':roi2.StructureSetDate}})
    
        roi2.StructureSetDate = newDate
        
        # Update the dictionary changesDict:
        changesDict.update({'toStructureSetDate':roi2.StructureSetDate})
    
    # - Update the Structure Set Time:
    #newTime = time.strftime("%H%M%S", time.gmtime())
    
    print('\nChanging the Structure Set Time from:\n', roi2.StructureSetTime, \
          '\nto:\n', newTime)
    
    # Update the dictionary changesDict:
    changesDict.update({'fromStructureSetTime':roi2.StructureSetTime})
    
    roi2.StructureSetTime = newTime
    
    # Update the dictionary changesDict:
    changesDict.update({'toStructureSetTime':roi2.StructureSetTime})
    
    # Create a new filename for the DICOM-RTSTRUCT file:
    #newRoiFname = 'AIM_' + roi2.StructureSetDate + '_' \
    #+ roi2.StructureSetTime + '.dcm'
    
    # Create a new filename for the DICOM-RTSTRUCT file:
    # - Start by getting the first part of the filename of the original 
    # DICOM-RTSTRUCT file:
    #origRoiDir, origRoiFname = os.path.split(roiDict1['Roi fpath'])
    fromRoiFname = xnatDict['fromRoiFname']
    
    
    # - Since the format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
    # "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
    # or similar, find the first "_" character and use the characters that 
    # precede it in the new file name:
    #ind = origRoiFname.index('_')
    ind = fromRoiFname.index('_')
    
    # - The new filename will be the same string of characters in origRoiFname
    # including the '_', followed by the recently modified Structure Set Date
    # and Structure Set Time:
    toRoiFname = fromRoiFname[0:ind+1] + roi2.StructureSetDate + '_' \
    + roi2.StructureSetTime + '.dcm'
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    toRoiFpath = os.path.join(xnatDict['toRoiDir'], toRoiFname)
    
   
    # The dictionary xnatDict will be updated (below) and exported to a JSON 
    # file in the path defined by xnatDict['toRoiDir'].
    # - Create a filename for the JSON file:
    """
    Use the same filename as the RTSTRUCT file but change the file extension.
    """
    jsonFname = os.path.splitext(toRoiFname)[0] + '.json'
    
    # Create a full filepath for the JSON file:
    jsonFpath = os.path.join(xnatDict['toRoiDir'], jsonFname)
    
    
    # Update the dictionary xnatDict:
    xnatDict.update({'changesToRoi':changesDict})
    xnatDict.update({'toRoiFname':toRoiFname})
    xnatDict.update({'toRoiFpath':toRoiFpath})
    xnatDict.update({'jsonFname':jsonFname})
    xnatDict.update({'jsonFpath':jsonFpath})
    
    
    # Print some useful info:
    #fromRFORSrefsopuids
    #toRFORSrefsopuids
    
    # Cycle through each Reference SOP Instance UID for each Contour Image
    # Sequence:
    for i in range(len(fromRFORSrefsopuids)):
        refsopuid = fromRFORSrefsopuids[i]
        
        # Find the index of the key (SOP Instance UID) in dicomDict1_to_2 that 
        # matches the Reference SOP Instance UID in roiDict1:
        ind = D1to2keys.index(refsopuid)
        
        # Get the slice no, Image Position, Slice Location and distances in
        # x, y, z and r for the 1st and 2nd DICOM scan (collection):
        sliceNo1 = dicomDict1_to_2[D1to2keys[ind]]['Dicom frame no']
        sliceNo2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom frame no']
        
        imPos1 = dicomDict1_to_2[D1to2keys[ind]]['Image pos']
        imPos2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Image pos']
        
        sliceLoc1 = dicomDict1_to_2[D1to2keys[ind]]['Slice loc']
        sliceLoc2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Slice loc']
    
        dx = dicomDict1_to_2[D1to2keys[ind]]['x2 - x1']
        dy = dicomDict1_to_2[D1to2keys[ind]]['y2 - y1']
        dz = dicomDict1_to_2[D1to2keys[ind]]['z2 - z1']
        dr = dicomDict1_to_2[D1to2keys[ind]]['r2 - r1']
        
        print(f'\n\nContour Image Sequence no {i+1}:\n', \
              f'\n   The contour(s) was copied from slice no {sliceNo1} in',\
              f'the original DICOM scan (collection)',\
              f'\n   to slice no {sliceNo2} in the new scan (collection):',\
              f'\n   - from Image Position {imPos1}',\
              f'\n                      to {imPos2}',\
              f'\n   - from Slice Location {sliceLoc1} to {sliceLoc2}',\
              f'\n   - with distances:',\
              f'\n       dx = {dx}',\
              f'\n       dy = {dy}',\
              f'\n       dz = {dz}',\
              f'\n       dr = {dr}')

    
    
    # Export the dictionary xnatDict to a JSON file:
    json.dump(xnatDict, open(jsonFpath, 'w'))
    
    print('\nJSON file created:\n', jsonFpath, '\n\n')
    
    
    # Save the new DICOM-RTSTRUCT file:
    roi2.save_as(toRoiFpath)
    
    #return roi2
    #return xnatDict
    return roi2, xnatDict
    



"""
Function:
    CopyRoi()
    
Purpose:
    Copy a ROI from one DICOM scan (collection) to another.

Input:
    fromRoiProjLab     - Project Label of the ROI Collection
    
    fromRoiSubjLab     - Subject Label of the ROI Collection
    
    fromRoiExpLab      - Experiment Label of the ROI Collection
    
    fromRoiLabDateTime - DateTime within the ROI Label of the ROI Collection
    
    fromDicomScanLab   - Scan Label of the DICOM scan (collection) with
                        ROI Collection
    
    toDicomProjLab     - Project Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomSubjLab     - Subject Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomExpLab      - Experiment Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomScanLab     - Scan Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
                        
    searchBy           - How to define the closest slice with possible values:
                        -- 'x' 
                        -- 'y'
                        -- 'z'
                        -- 'r'
                        to find closest slice along x/y/z direction or r, where
                        r = sqrt(dx^2 + dy^2 + dz^2)
    
    del_dicoms         - Delete downloaded DICOMs if True
    
    debug              - Print results if True

    
Returns:
    roi2               - ROI object that contains the contours in roi1 but 
                        modified for overlaying over the 2nd DICOM scan;
                       - A new DICOM-RTSTRUCT file will be saved to the location
                       specified
                   
    xnatDict           - The updated dictionary
                       - The dictionary is exported as a JSON file 
"""



def CopyRoi(fromRoiProjLab, fromRoiSubjLab, fromRoiExpLab, fromRoiLabDateTime,\
            fromDicomScanLab, toDicomProjLab, toDicomSubjLab, toDicomExpLab, \
            toDicomScanLab, searchBy, del_dicoms, debug):
    
    # Import packages and functions:
    
    #import importlib
    #import DicomHelperFuncs
    #importlib.reload(DicomHelperFuncs)
    from DicomHelperFuncs import DownloadFilesForRoiCopy
    from DicomHelperFuncs import GetDicomData
    #from DicomHelperFuncs import GetRoiData
    from DicomHelperFuncs import GetClosestDicomPositions
    from DicomHelperFuncs import ModifyRoi
    import time
    import os
    import shutil
    
    # Keep track of the time it takes to execute certain functions:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    # Get the current date:
    todaysDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    # Assign the root directory:
    rootDir = os.path.join(r'C:\Temp', todaysDate)
    
    # Create a dictionary to store various variables:
    xnatDict = {# Define the XNAT login details:
                'xnatAddress':'http://10.1.1.17',\
                'username':'admin',\
                'password':'admin',\
        
                # Define the REST variables that point to the ROIs to be copied:
                'fromRoiProjLab':fromRoiProjLab,\
                'fromRoiSubjLab':fromRoiSubjLab,\
                'fromRoiExpLab':fromRoiExpLab,\
                #'fromRoiLab':'AIM_20200304_110009',\
                #'fromRoiLab':'AIM_20200306_114104',\
                'fromRoiLabDateTime':fromRoiLabDateTime,\
        
                # Define the REST variables that point to the DICOM files that the ROI applies to:
                'fromDicomScanLab':fromDicomScanLab,\
        
                # Define the REST variables that point to the DICOM files that the ROI is to be applied to:
                'toDicomProjLab':toDicomProjLab,\
                'toDicomSubjLab':toDicomSubjLab,\
                'toDicomExpLab':toDicomExpLab,\
                'toDicomScanLab':toDicomScanLab,\
        
                # Define the label to assign to Structure Set Label and to ROI Name of the new RTSTUCT file:
                #'toRoiStructSetLab' = 'ROIs_t2_blade_tra_slice_5'
        
        
                # Define the directories to download files to:
                'rootDir':rootDir, \
                'fromDicomDir':os.path.join(rootDir, 'From DICOM'),\
                'fromRoiDir':os.path.join(rootDir, 'From RTSTRUCT'),\
                'toDicomDir':os.path.join(rootDir, 'To DICOM'),\
                'toRoiDir':os.path.join(rootDir, 'To RTSTRUCT')
               }
    
    
    # Download the DICOM-RTSTRUCT and DICOM files from XNAT:
    xnatDict = DownloadFilesForRoiCopy(xnatDict, debug=debug)
    
    # times[1] = time after downloading of files from XNAT:
    times.append(time.time())
    
    print(f'\nTook {round(times[1] - times[0], 1)} s to download data from ',\
                    'XNAT')
    
    # Parse the ROI Collection to be copied:
    """ This is not needed following changes to the helper functions """
    #roiDict1 = GetRoiData(xnatDict['fromRoiFpath'])


    # Parse the DICOM data that corresponds to the ROI Collection to be copied:
    dicomDict1 = GetDicomData(xnatDict['fromDicomDir'])
    
    
    # Parse the DICOM data that the ROI Collection is to be copied to:
    dicomDict2 = GetDicomData(xnatDict['toDicomDir'])

    
    # Get the mapping from the 1st to 2nd DICOM scans using the z-coordinate 
    # to find the nearest slice:
    dicomDict1_to_2 = GetClosestDicomPositions(dicomDict1, dicomDict2, \
                                               searchBy=searchBy)
    
    
    
    # Print some results for diagnostics purposes:
    if debug:
        from DicomHelperFuncs import IndexDictInTier2
        
        key1 = IndexDictInTier2(dicomDict1, tier2key='Dicom frame no', tier2val=5)


        print('Contour copied from Image Position =', \
              dicomDict1[key1]['Image pos'], '\n')
        
        # loop through each key in dicomDict2:
        keys2 = list(dicomDict2.keys())
        
        for i in range(len(keys2)):
            print(f'Slice {i+1} Image Position =', \
                  dicomDict2[keys2[i]]['Image pos'])
            
            
        
    # Copy and modify the ROI Collection (from the original DICOM scan) so
    # they may be overlaid on the new DICOM scan:
    roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2)
    
    # times[2] = time after parsing files, finding closest slice locations 
    # and modifying the ROIs:
    times.append(time.time())
    
    print(f'\nTook {round(times[2] - times[1], 1)} s to copy and modify the ',\
                    'ROIs')
    
    # Do some housekeeping - delete the DICOM files in fromDicomDir and 
    # toDicomDir:
    if del_dicoms:
        print('\nDeleting temporary DICOM files..')
        for directory in [xnatDict['fromDicomDir'], xnatDict['toDicomDir']]:
            for root, dirs, files in os.walk(directory):
                for f in files:
                    os.unlink(os.path.join(root, f))
                for d in dirs:
                    shutil.rmtree(os.path.join(root, d))
                    
    
    print(f'\nTook {round(times[2] - times[0], 1)} s to do everything')
    
    return roi2, xnatDict



"""
Function:
    IndexDictInTier2()
    
Purpose:
    Find the key in a dictionary with a particular tier 2 value in the list
    of tier 2 keys.

Input:
    dictionary - a dictionary of the form {tier1keys:{tier2keys:tier2values}}
    
    tier2key   - tier 2 key that the value is to be found in
    
    tier2val   - tier 2 value to be searched for
    
    
Returns:
    tier1key   - tier 1 key that has a tier 2 key with desired tier 2 value 
"""


def IndexDictInTier2(dictionary, tier2key, tier2val):
    for tier1key in dictionary.keys():
        if dictionary[tier1key][tier2key]==tier2val:
            return tier1key
    return -1
    

"""
Function:
    CompareROIandDICOMtags()
    
Purpose:
    Print tags from DICOM-RTSTRUCT and from DICOM files for comparison.

Input:
    dicomFpaths - Full file paths of DICOM files
    roi         - ROI object from DICOM-RTSTRUCT file
    imageInds   - indeces that link each Contour Image Sequence with its
                  corresponding DICOM file number
    contourInds - indeces that link each Contour Sequence with its
                  corresponding DICOM file number
    
Returns:
    None
"""


def CompareROIandDICOMtags(dicomFpaths, roi, imageInds, contourInds):
    
    # Import packages:
    import pydicom
    
    
    print('The SOP Instance UID from the ROI and from DICOMs indexed by', \
          'imageInds:\n')

    print('roi.SOPInstanceUID   =', roi.SOPInstanceUID)
    print('')
    for i in range(len(imageInds)):
        fpath = dicomFpaths[imageInds[i]]
        # Read in the DICOMs:
        dicom = pydicom.read_file(fpath)
        if i==0:
            print('dicom.SOPInstanceUID =', dicom.SOPInstanceUID)
        else:
            print('                      =', dicom.SOPInstanceUID)
    
    print('')
    print('')
    print('**************************************************************************************************************')
    print('')
    print('')
    
    print('Study Instance UID and Ref SOP Instance (unexpected match) from the ROI, and DICOM Study Instance UID:\n')
    
    print('roi.StudyInstanceUID   =', roi.StudyInstanceUID, '\n')
    print('roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID \n', \
          '                       =', \
          roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID, '\n')
    print('dicom.StudyInstanceUID =', dicom.StudyInstanceUID)
    
    print('')
    print('')
    print('**************************************************************************************************************')
    print('')
    print('')
    
    print('Series Instance UID from ROI and DICOM:\n')
    
    print('roi.SeriesInstanceUID   =', roi.SeriesInstanceUID, '\n')
    print('roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID \n', \
          '                        =', \
          roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID, '\n')
    print('dicom.SeriesInstanceUID =', dicom.SeriesInstanceUID)
    
    print('')
    print('')
    print('**************************************************************************************************************')
    print('')
    print('')
    
    print('Frame of Reference UID from ROI and DICOM:\n')
    
    print('roi.FrameOfReferenceUID                                       \n=', roi.FrameOfReferenceUID, '\n')
    print('roi.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID \n=', roi.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID, '\n')
    print('roi.StructureSetROISequence[0].ReferencedFrameOfReferenceUID  \n=', roi.StructureSetROISequence[0].ReferencedFrameOfReferenceUID, '\n')
    print('dicom.FrameOfReferenceUID                                     \n=', dicom.FrameOfReferenceUID)
    print('')
    print('')
    
    print('')
    print('')
    print('**************************************************************************************************************')
    print('')
    print('')
    
    print('Referenced SOP Instance UID from Contour Image Sequence and SOP Instance UID from DICOMs \n', 
          '(both indexed by imageInds in turn):\n')
    
    print('roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0]\n', \
          '.ContourImageSequence[i].ReferencedSOPInstanceUID')
    print('dicom.SOPInstanceUID \n')
    for i in range(len(imageInds)):
        print(roi.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence[i].ReferencedSOPInstanceUID)
        
        # Load the imageInds[i]^th DICOM file in dicomFpaths:
        fpath = dicomFpaths[imageInds[i]]
        # Read in the DICOM:
        dicom = pydicom.read_file(fpath)
        print(dicom.SOPInstanceUID)
        print('')
    
    print('')
    print('')
    print('**************************************************************************************************************')
    print('')
    print('')
            
    print('Referenced SOP Instance UID from Contour Sequence and SOP Instance UID from DICOMs \n', 
          '(both indexed by contourInds in turn):\n')
    
    print('roi.ROIContourSequence[0].ContourSequence[i].ContourImageSequence[0].ReferencedSOPInstanceUID')
    print('dicom.SOPInstanceUID \n')
    for i in range(len(contourInds)):
        print(roi.ROIContourSequence[0].ContourSequence[i].ContourImageSequence[0].ReferencedSOPInstanceUID)
        
        # Load the imageInds[i]^th DICOM file in dicomFpaths:
        fpath = dicomFpaths[contourInds[i]]
        # Read in the DICOM:
        dicom = pydicom.read_file(fpath)
        print(dicom.SOPInstanceUID)
        print('')
        
    return






#############################################################################
#############################################################################
    """ OLD VERSIONS OF FUNCTIONS ABOVE """
#############################################################################
#############################################################################
    

""" PREVIOUS VERSION """

""" Function:
    
    DownloadFilesForRoiCopy_OLD()
    
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
    xnatDict - Dictionary containing XNAT login details, REST path variables 
               and directory paths, i.e.:
                   
                xnatDict = {# Define the XNAT login details:
                            'xnatAddress':'http://10.1.1.17',\
                            'username':'admin',\
                            'password':'password',\
                            
                            # Define the REST variables that point to the ROIs 
                            # to be copied:
                            'roiProjLab':'BrainTumorProg',\
                            'roiSubjLab':'PGBM-002',\
                            'roiExpLab':'PGBM-002_MR_1',\
                            'roiLab':'RTSTRUCT_20200123_145731',\
                            
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
    
    del_dicoms - Delete downloaded DICOMs if True
    
    debug      - Print results if True
    
Returns:
    xnatDict - with the added full filepath of a new DICOM-RTSTRUCT file in the 
               directory specified by xnatDict['toRoiDir'], with a new file 
               name generated with the current date and time
"""


def DownloadFilesForRoiCopy_OLD(xnatDict, del_dicoms, debug):
    
    # Import functions and packages:
    import xnat
    import os, shutil
    #import json
    #import importlib
    #import ModifyROIs
    #importlib.reload(ModifyROIs)
    from ModifyROIs import ModifyROIs

    if debug:
        print('The input variables (xnatDict):\n\n')
        print(xnatDict)
        print('\n\n')
    
    # Create XNAT session:
    # Note:  xnatDict['xnatAddress'] = 'http://10.1.1.17'
    session = xnat.connect(xnatDict['xnatAddress'], \
                           user=xnatDict['username'], \
                           password=xnatDict['password'])
    
    # Get the directories to download files to:
    rootDir = xnatDict['rootDir']           # e.g. r'C:\Temp\2020-02-28'
    fromDicomDir = xnatDict['fromDicomDir'] # e.g. r'C:\Temp\2020-02-28\From DICOM'
    fromRoiDir = xnatDict['fromRoiDir']     # e.g. r'C:\Temp\2020-02-28\From RTSTRUCT' 
    toDicomDir = xnatDict['toDicomDir']     # e.g. r'C:\Temp\2020-02-28\To DICOM'
    toRoiDir = xnatDict['toRoiDir']         # e.g. r'C:\Temp\2020-02-28\To RTSTRUCT'  
    
    # Combine all directory paths:
    allDirs = [rootDir, fromRoiDir, fromDicomDir, toRoiDir, toDicomDir]
    
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
    #roi = session.projects[RoiProjLab].subjects[RoiSubjLab].experiments[RoiExpLab].assessors[RoiLab].resources['RTSTRUCT'].files[0]
    roiProj = session.projects[xnatDict['fromRoiProjLab']]
    roiSubj = roiProj.subjects[xnatDict['fromRoiSubjLab']]
    roiExp = roiSubj.experiments[xnatDict['fromRoiExpLab']]
    
    # The ROI object:
    roi = roiExp.assessors[xnatDict['fromRoiLab']].resources['RTSTRUCT'].files[0]
    
    # Get the filename:
    fromRoiFname = roi.id
    
    # The ROI filepath:
    fromRoiFpath = os.path.join(fromRoiDir, fromRoiFname)
    
    # Download the DICOM-RTSTRUCT file:
    print('\nDownloading the DICOM-RTSTRUCT file to be copied:', fromRoiFname)
    roi.download(fromRoiFpath)
    
    
    # Get the DICOM objects of the series the ROIs relate to:
    #dicoms = session.projects[DicomProjLab].subjects[DicomSubjLab].experiments[DicomExpLab].scans[DicomScanLab].files
    #dicoms = session.projects[fromDicomProjLab].subjects[fromDicomSubjLab].experiments[fromDicomExpLab].scans[fromDicomScanLab]
    dicomsProj = session.projects[xnatDict['fromDicomProjLab']]
    dicomsSubj = dicomsProj.subjects[xnatDict['fromDicomSubjLab']] 
    dicomsExp = dicomsSubj.experiments[xnatDict['fromDicomExpLab']]
    dicoms = dicomsExp.scans[xnatDict['fromDicomScanLab']]
    
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
    #dicoms = session.projects[toDicomProjLab].subjects[toDicomSubjLab].experiments[toDicomExpLab].scans[toDicomScanLab]
    dicomsProj = session.projects[xnatDict['toDicomProjLab']]
    dicomsSubj = dicomsProj.subjects[xnatDict['toDicomSubjLab']] 
    dicomsExp = dicomsSubj.experiments[xnatDict['toDicomExpLab']]
    dicoms = dicomsExp.scans[xnatDict['toDicomScanLab']]
    
    # Download each DICOM by iterating for each file in dicomObjects:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file is to be overlaid on..')
    for f in range(len(dicoms.files)):
        # Get the filename:
        dicomId = dicoms.files[f].id
        
        # Download the DICOM:
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        dicoms.files[f].download(os.path.join(toDicomDir, dicomId))
        
    
    # Update the dictionary:
    xnatDict.update({'fromRoiFname':fromRoiFname})
    xnatDict.update({'fromRoiFpath':fromRoiFpath})

    
    
    """
    04/03:
    The call of ModifyROIs should be done outside of this function. 
    This function should strictly download the files from XNAT.
    The filepath to the JSON file can be stored in the xnatDict so that 
    the additions to the JSON file can be done within ModifyROIs.
    
    
    
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
    """
    
    # Do some housekeeping - delete the DICOM files in fromDicomDir and 
    # toDicomDir:
    if del_dicoms:
        print('Deleting temporary DICOM files..')
        for directory in [fromDicomDir, toDicomDir]:
            for root, dirs, files in os.walk(directory):
                for f in files:
                    os.unlink(os.path.join(root, f))
                for d in dirs:
                    shutil.rmtree(os.path.join(root, d))
    
    
    return xnatDict
  