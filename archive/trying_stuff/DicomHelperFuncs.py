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
    
    11/03/20 Note: Even though the DICOMs are sorted in a natural way (e.g.
    0000.dcm, 0001.dcm, 0002.dcm) it doesn't guarantee that the files are
    sorted in the sequence of the slices!  So rather than relying on the 
    file names I will need to use the Image Position (Patient) to sort. In the
    modifications to the function below sorting along the slice direction will 
    be an option. 
    

Input:
    dirPath    - path to directory containing DICOM files
    
    sortMethod - method to use for sorting of files with acceptable values:
                -- 'none' (or 'None' or any other strings) = no sorting
                -- 'natural' (or 'Natural') = using natsort
                -- 'slices' (or 'Slices') = sorted along the slice direction 
    
Returns:
    filePaths - full file paths of all DICOM files in dirPath either:
                -- without sorting
                -- sorted in a natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ...
                rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
                -- sorted along the scan direction (The slice direction
                will be inferred by considering the min/max positions along the 
                x, y and z directions.)
"""


def GetDicomFpaths(dirPath, sortMethod):
    
    # Import packages:
    import os
    
    # Create an empty list of all DICOM file paths:
    filePaths = []

    # Use os to get list of file paths of DICOM-only files:
    for dirName, subdirList, fileList in os.walk(dirPath):
        for fileName in fileList:
            if '.dcm' in fileName.lower():  # check for DICOM files only
                filePaths.append(os.path.join(dirName,fileName))
    
    if sortMethod=='natural' or sortMethod=='Natural':
        # Import natsort:
        import natsort
        
        # Sort files in natural way:
        filePaths = natsort.natsorted(filePaths)
        
    if sortMethod=='slices' or sortMethod=='Slices':
        # Import packages:
        import pydicom
        import numpy as np
        
        # Create an empty array to store the positions:
        positions = []
        
        for filePath in filePaths:
            #print('\nfilePath =', filePath)
            
            # Read in DICOM:
            dicom = pydicom.read_file(filePath)
            
            # Parse the Image Position (Patient):
            position = [float(item) for item in dicom.ImagePositionPatient]
            
            # Add position for this DICOM to positions:
            positions.append(position)
            
        # Convert to numpy array:
        positions = np.array(positions)

        # Get the maximum change of positions along x, y and z:
        dPositions = np.max(positions, axis=0) - np.min(positions, axis=0)

        # Get the index of the axis with the greatest change of position:
        ind = np.where(dPositions==np.max(dPositions))[0][0]
        
        """
        Solution on how to sort along a specific column:
        
        https://stackoverflow.com/questions/22698687/how-to-sort-2d-array-numpy-ndarray-based-to-the-second-column-in-python
        
        The following are equivalent:
        
        positions[positions[:, 2].argsort()]
        positions[np.argsort(positions[:, 1])]
        
        """
        
        sortInds = positions[:, ind].argsort()
        
        #positionsSorted = positions[sortInds]
        
        # Sort filepaths with sortInds:
        #filePaths = filePaths[sortInds] # <-- TypeError: list indices must be 
        # integers or slices, not list
        filePaths = [filePaths[i] for i in sortInds]
        
    
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
    
    debug      - Prints some useful results if True
    
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
    #toPosCorrDicomDir = xnatDict['toPosCorrDicomDir']
    #toRegDicomDir = xnatDict['toRegDicomDir']
    toRoiDir = xnatDict['toRoiDir']         # e.g. r'C:\Temp\2020-02-28\To RTSTRUCT'  
    
    # Combine all directory paths:
    #allDirs = [rootDir, fromRoiDir, fromDicomDir, toRoiDir, toDicomDir]
    allDirs = [rootDir, fromDicomDir, toDicomDir, fromRoiDir, toRoiDir]
    #allDirs = [rootDir, fromDicomDir, toDicomDir, toPosCorrDicomDir, \
    #           toRegDicomDir, fromRoiDir, toRoiDir]
    
    # Create directories if they don't exist:
    for dirPath in allDirs:
        if debug:
            print('\nChecking if directory', dirPath, 'exists')
        try:
            os.mkdir(dirPath)
            if debug:
                print('\nDirectory', dirPath, 'created')
        except FileExistsError:
            if debug:
                print('\nDirectory', dirPath, 'already exists')
    
    """ 
    Note:
    Using .download_dir(directory) downloads to a complicated XNAT-like folder 
    structure, so instead will use .download(full_file_path) 
    """
    
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
    
    # Download the DICOM-RTSTRUCT file if it hasn't already been:
    if os.path.isfile(fromRoiFpath):
        print('\nThe DICOM-RTSTRUCT file has already been downloaded.')
    else:
        print('\nDownloading the DICOM-RTSTRUCT file to be copied:', fromRoiFname)
    
        roi.download(fromRoiFpath)
    
    
    # Get the DICOM objects of the series the ROIs relate to:
    dicoms = roiExp.scans[xnatDict['fromDicomScanLab']]
    
    # Iterate for each file in dicomObjects:
    print('\nDownloading the DICOM files that the DICOM-RTSTRUCT file relate to..')
    for f in range(len(dicoms.files)):
        # Get the filename:
        dicomId = dicoms.files[f].id
        
        # Create the filepath:
        dicomFpath = os.path.join(fromDicomDir, dicomId)
        
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        
        # Download the DICOM (if it hasn't already been):
        if not os.path.isfile(dicomFpath):
            dicoms.files[f].download(dicomFpath)
        
        
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
        
        # Create the filepath:
        dicomFpath = os.path.join(toDicomDir, dicomId)
        
        if debug:
            print('\n     Downloading DICOM file', dicomId)
        
        # Download the DICOM (if it hasn't already been):
        if not os.path.isfile(dicomFpath):
            dicoms.files[f].download(dicomFpath)
        
    
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
    dirPath    - path to directory containing DICOM files
    
    sortMethod - method to use for sorting of files using GetDicomFpaths() with 
                acceptable values:
                -- 'none' (or 'None' or any other strings) = no sorting
                -- 'natural' (or 'Natural') = using natsort
                -- 'slices' (or 'Slices') = sorted along the slice direction 
    
Returns:
    dicomDict  - dictionary containing Image Position (Patient), Image 
                Orientation (Patient), Slice Location, Pixel Spacing and Pixel 
                Array with SOP Instance UID as keys
"""

def GetDicomData(dirPath, sortMethod):
    
    # Import packages:
    import pydicom
    from DicomHelperFuncs import GetDicomFpaths
    
    # Get the filepaths of the DICOM files:
    dicomsFpaths = GetDicomFpaths(dirPath, sortMethod)
    
    # Initialise dictionary to store everything:
    dicomDict = {}
    
    # Loop though each DICOM filepath:
    for f in range(len(dicomsFpaths)):
        fpath = dicomsFpaths[f]
        
        #print('\nfpath =', fpath)
        
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
        dicomDict[uid] = {'Dicom slice no':f+1, \
                          'Dicom fpath':fpath, \
                          'Image pos':position, \
                          'Image orient':orientation, \
                          'Slice loc':sliceLocation, \
                          'Pix spacing':pixSpacing, \
                          'Pix array':pixelArray, \
                          'Dicom':dicom
                         }
    
    return dicomDict



"""
Function:
    GetRoiData()
    
Purpose:
    Get No of contours (e.g. 3), Contour nos (e.g. [1, 2, 3]), No of contour
    pts (e.g. [44, 45, 116]), and Contour Data from a ROI filepath, and store 
    values in a dictionary using the Reference SOP Instance UID as keys.
    
Changes:
    10/03/20:
        Allow input to be either a filpath or a directory. Originally the 
        input was a filepath to avoid ambiguity for cases where more than
        DICOM-RTSTRUCT file may be present.  But it's more convenient to 
        point to a directory, as in most cases, there will only be one file.
        
        Warning: If there's more than one file the first file will be used.

Input:
    filePath_or_Dir  - filepath to a DICOM-RTSTRUCT (ROI Collection) file
                        or
                       directory to a DICOM-RTSTRUCT file
    
Returns:
    roiDict - dictionary containing No of contours, Contour nos, No of contour
              pts, and Contour Data with Referenced SOP Instance UID as keys
"""

def GetRoiData(filePath_or_Dir):
    
    # Import packages:
    import pydicom
    import os
    
    # Check if input is a filepath:
    if os.path.isfile(filePath_or_Dir):
        filePath = filePath_or_Dir
        
        # Read in the DICOM-RTSTRUCT file:
        roi = pydicom.read_file(filePath)
        
    # Check if input is a directory:
    if os.path.isdir(filePath_or_Dir):
        # Get list of file names:
        fileNames = os.listdir(filePath_or_Dir)
        
        if len(fileNames)==0:
            print('No files in', filePath_or_Dir)
        
        dcm_fileNames = []
        
        # Check for DICOM files:
        for fileName in fileNames:
            if fileName.endswith('.dcm'):
                dcm_fileNames.append(fileName)
                
        if len(dcm_fileNames)==0:
            print('No .dcm files in', filePath_or_Dir)
          
        # Deal with case of multiple .dcm files:
        elif len(dcm_fileNames)>1:
            print('Multiple .dcm files in', filePath_or_Dir, '.', \
                  'Reading in the first file.')
            
            # Complete the filepath:
            filePath = os.path.join(filePath_or_Dir, dcm_fileNames[0])
            
            # Read in the DICOM-RTSTRUCT file:
            roi = pydicom.read_file(filePath)
            
        # Else there is only one .dcm file.
        else: 
            # Complete the filepath:
            #roiFpath = os.path.join(filePath_or_Dir, dcm_filenames)
            filePath = os.path.join(filePath_or_Dir, dcm_fileNames[0])
            
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
        try:
            contourNo.append(int(sequence.ContourNumber))
        except AttributeError:
            #print('\nDICOM-RTSTRUCT file does not contain Contour Number tag.')
            contourNo = []
        
        # - the Number Of Contour Points:
        NcontourPts.append(int(sequence.NumberOfContourPoints))
        
        # - the Contour Data:
        #contourData.append(np.array(sequence.ContourData).astype(np.float))
        contourData.append([float(item) for item in sequence.ContourData])
    
        # Store the values for this uid. Values will be over-written if there's more than
        # one contour for this uid, but the lists have been appended already:
        roiDict[uid] = {'Roi fpath':filePath, \
                        'No of contours':Ncontours, \
                        'Contour no':contourNo, \
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
    GetClosestSlices()
    
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
    
    scanPosCorr - Scan position correction that needs to be applied to the DICOM 
                scan (collection) to account for different positional reference 
                points (as is the case for different imaging modalities).  The
                correction will be applied to the same dimension as defined by 
                searchBy.  Accepatable values are any float. Note that if
                searchBy='r' the correction will not be applied.
                
    debug       - Print results if True
    
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
                      
    dicomDict2posCorr - dicomDict2 with correction applied to the positional
                        dimension given by searchBy (if scanPosCorr=0, it will
                        simply be a copy of dicomDict2)
"""

def GetClosestSlices(dicomDict1, dicomDict2, searchBy, scanPosCorr, debug):
    
    # Import packages:
    import copy
    
    # Initialise dicomDict1_to_2 as a copy of dicomDict1:
    dicomDict1_to_2 = copy.deepcopy(dicomDict1)
    
    # Get the list of keys of the two DICOM dictionaries:
    keys1 = list(dicomDict1.keys())
    keys2 = list(dicomDict2.keys())
    

    # First make corrections to the Image Positions of a copy of dicomDict2
    dicomDict2posCorr = copy.deepcopy(dicomDict2)
    
    #  If scanPosCorr is non-zero, loop through each key in 
    # dicomDict2posCorr, amending the Image Position:
        
    if scanPosCorr:
        for key2 in keys2:
            # For debugging purposes, store pos2 before and after the change:
            pos2_debug = []
            
            # Get the Image Position:
            pos2 = dicomDict2posCorr[key2]['Image pos']
            
            #print('\npos2 was', pos2)
            pos2_debug.append(pos2)
            
            # Apply correction to scan position based on the chosen searchBy:
            if searchBy=='x':
                pos2[0] = pos2[0] + scanPosCorr
            if searchBy=='y':
                pos2[1] = pos2[1] + scanPosCorr
            if searchBy=='z':
                pos2[2] = pos2[2] + scanPosCorr
                
            #print('\npos2 is now', pos2)
            pos2_debug.append(pos2)
            
            # Apply the corrected position to the dictionary:
            dicomDict2posCorr[key2]['Image pos'] = pos2
            
            if debug:
                print('\npos2 was', pos2_debug[0], ', and changed to',\
                      pos2_debug[1])
    
    
    # For every key (SOP UID) in dicomDict1 find the SOP UID of the nearest 
    # slice in dicomDict2posCorr.

    # Store the indices of the keys of dicomDict2posCorr whose positions are 
    # closest to the positions of the slices in dicomDict1:
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
            #pos2 = dicomDict2[key2]['Image pos'] # this changes dicomDict2!
            pos2 = dicomDict2posCorr[key2]['Image pos']
            
            #print('pos2 =', pos2)
            
            # Calculate the distance between pos1 and pos2 along x, y and z:
            dx = pos2[0] - pos1[0]
            dy = pos2[1] - pos1[1]
            dz = pos2[2] - pos1[2]
            
            #print('\ndz =', dz)
            
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
        
        # Update dicomDict1_to_2 with the details of dicomDict2posCorr that 
        # most closesly matchs the position along the dimension indicated by
        # searchBy:
        if searchBy=='x':
            ind = copy.deepcopy(indX)
        if searchBy=='y':
            ind = copy.deepcopy(indY)
        if searchBy=='z':
            ind = copy.deepcopy(indZ)
        if searchBy=='r':
            ind = copy.deepcopy(indR)
            
        print('\nind =', ind)
            
        dicomDict1_to_2[key1].update({'Closest Dicom SOP UID':\
                                      keys2[ind], \
                                      'Closest Dicom slice no':\
                                      #dicomDict2[keys2[ind]]['Dicom slice no'],\
                                      dicomDict2posCorr[keys2[ind]]['Dicom slice no'],\
                                      'Closest Dicom fpath':\
                                      #dicomDict2[keys2[ind]]['Dicom fpath'],\
                                      dicomDict2posCorr[keys2[ind]]['Dicom fpath'],\
                                      'Closest Image pos':\
                                      #dicomDict2[keys2[ind]]['Image pos'],\
                                      dicomDict2posCorr[keys2[ind]]['Image pos'],\
                                      'Closest Image orient':\
                                      #dicomDict2[keys2[ind]]['Image orient'],\
                                      dicomDict2posCorr[keys2[ind]]['Image orient'],\
                                      'Closest Slice loc':\
                                      #dicomDict2[keys2[ind]]['Slice loc'],\
                                      dicomDict2posCorr[keys2[ind]]['Slice loc'],\
                                      'Closest Pix spacing':\
                                      #dicomDict2[keys2[ind]]['Pix spacing'],\
                                      dicomDict2posCorr[keys2[ind]]['Pix spacing'],\
                                      'x2 - x1':dX[ind],
                                      'y2 - y1':dY[ind],
                                      'z2 - z1':dZ[ind],
                                      'r2 - r1':dR[ind], \
                                      'Closest Pix array':\
                                      #dicomDict2[keys2[ind]]['Pix array'],\
                                      dicomDict2posCorr[keys2[ind]]['Pix array'],\
                                      'Closest Dicom':\
                                      #dicomDict2[keys2[ind]]['Dicom'],\
                                      dicomDict2posCorr[keys2[ind]]['Dicom'],\
                                     })
            
    #return dicomDict1_to_2
    return dicomDict1_to_2, dicomDict2posCorr




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
                    
    regMethod      - Type of registration to perform on the DICOM scan
                     (collection) that the ROI Collection is to be copied to.  
                     Acceptable values include:
                        -- 'rigid' = perform rigid registration
                        -- 'affine' = perform affine registration
                        -- 'non-rigid' = perform non-rigid registration
                        -- any other character string (including '') = do not
                            perform any registration
                     This is needed because the ROI tags will either need to 
                     reflect the original tags of the 2nd DICOM scan or the
                     new tags of the new DICOMs of the registered scan.
    
    
Returns:
    roi2           - ROI object that contains the contours in roi1 but 
                    modified for overlaying over the 2nd DICOM scan;
                   - A new DICOM-RTSTRUCT file will be saved to the location
                   specified
                   
    xnatDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths
                    (see DownloadFilesForRoiCopy() for details)
                   - The dictionary is exported as a JSON file
    
"""

def ModifyRoi(xnatDict, dicomDict1_to_2, regMethod):
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
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        fpath2 = dicomDict1_to_2[D1to2keys[0]]['Registered closest Dicom fpath']
    else:
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
        changesDict.update({'fromFrameOfReferenceUID':roi2.FrameOfReferenceUID})
        
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
        if regMethod in ['rigid', 'affine', 'non-rigid']:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
        else:
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
        if regMethod in ['rigid', 'affine', 'non-rigid']:
            fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
        else:
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
    
    print('\nChanging the Structure Set Time from:\n', roi2.StructureSetTime,\
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
    
    # The directory where the ROI is to be saved to:
    toRoiDir = xnatDict['toRoiDir']
    
    """
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Registered closest Dicom fpath']
    else:
        fpath2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom fpath']
    """
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    toRoiFpath = os.path.join(toRoiDir, toRoiFname)
    
   
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
        sliceNo1 = dicomDict1_to_2[D1to2keys[ind]]['Dicom slice no']
        sliceNo2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom slice no']
        
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
    #return roi2, xnatDict
    return roi1, roi2, xnatDict # 10/03/2020
    



"""
Function:
    ConvertNumpyDtype()
    
Purpose:
    Convert numpy array from one data type to another, preserving the
    dynamic range.

Input:
    arrayIn  - Numpy array whose data type is to be converted; acceptable data
               types include:
                  -- uint8
                  -- uint16
                  -- int8
                  -- int16
                  -- float
                        
    convertTo - Data type to convert arrayIn to
    
    debug     - If True, some results will be printed

    
Returns:
    arrayOut  - Numpy array with data type determined by convertTo   
    

Note:

At the moment this function only handles a limited number of possible cases of
"from-to" data types: uint8, uint16, int8, int16 and float.                 
                  
"""

def ConvertNumpyDtype(arrayIn, convertTo, debug):
    # Import packages:
    import numpy as np
    
    # Get info on input array, first assuming that it is an int dtype:
    try:
        info = np.iinfo(arrayIn.dtype)
        
    # If ValueError arises, assume it is a float:
    except ValueError:
        info = np.finfo(arrayIn.dtype)
    
    if debug:    
        print('\nBefore dtype:', arrayIn.dtype)
        print(f'\nBefore: {info}')
        print('\nconvertTo = ', convertTo)
      
    
    # Determine value to scale array up to:
    if 'uint8' in convertTo:
        maxValue = 255
    
    if 'uint16' in convertTo:
        maxValue = 65535
        
    if 'int8' in convertTo:
        maxValue = 127
    
    if 'int16' in convertTo:
        maxValue = 32767
        
    if 'float' in convertTo:
        maxValue = 3.4028235e+38
        
    if debug:    
        print('\nMax value possible for ' + convertTo + f' = {maxValue}')
    
    # Get min/max values:
    arrayMin = np.min(arrayIn)
    arrayMax = np.max(arrayIn)
    
    if debug:
        print(f'\nMin/Max of input array = [{arrayMin}, {arrayMax}]')
    
    #if arrayMin < 0:
    
    # Deal with negative values if convertTo is an unsigned integer
    # (by subtracting the min value):
    if 'u' in convertTo:
        array = arrayIn - arrayMin

        # Get min/max values:
        arrayMin = np.min(array)
        arrayMax = np.max(array)
        
        if debug:
            print(f'\nMin/Max of array after subtracting Min = [{arrayMin}, {arrayMax}]')
        
    else:
        array = arrayIn - 0
    
    # Normalise array to [0, 1] or [-1, 1] (depends on whether
    # the input array is signed or unsigned):
    #array = arrayIn.astype(np.float64) / info.max # is conversion to float useful?
    #array = arrayIn / info.max 
    array = array / maxValue
    
    # Get min/max values:
    arrayMin = np.min(array)
    arrayMax = np.max(array)
    
    if debug:
        print(f'\nMin/Max of array after normalising = [{arrayMin}, {arrayMax}]')
        
    
    # Scale array by maxValue:
    array = maxValue * array
    
    # Get min/max values:
    arrayMin = np.min(array)
    arrayMax = np.max(array)
    
    if debug:
        print(f'\nMin/Max of array after scaling by {maxValue} = [{arrayMin}, {arrayMax}]')
    
    # Convert array to desired dtype:
    if 'uint8' in convertTo:
        if debug:
            print('\nConverting to uint8')
        arrayOut = array.astype(np.uint8)
        
    if 'uint16' in convertTo:
        if debug:
            print('\nConverting to uint16')
        arrayOut = array.astype(np.uint16)
        
    if 'int8' in convertTo:
        if debug:
            print('\nConverting to int8')
        arrayOut = array.astype(np.int8)
        
    if 'int16' in convertTo:
        if debug:
            print('\nConverting to int16')
        arrayOut = array.astype(np.int16)
        
    if 'float' in convertTo:
        if debug:
            print('\nConverting to float')
        arrayOut = array.astype(np.float)
        
        
    # For some reason the attempts to convert to 
        
    
    # Get min/max values:
    arrayMin = np.min(arrayOut)
    arrayMax = np.max(arrayOut)
    
    if debug:
            print('\nMin/Max of array after converting to', convertTo, \
                  f'= [{arrayMin}, {arrayMax}]')
    
    # Re-check info:
    try:
        info = np.iinfo(arrayOut.dtype)
        
    # If ValueError arises, assume it is a float:
    except ValueError:
        info = np.finfo(arrayOut.dtype)
        
    if debug:
            print('\nAfter dtype:', arrayOut.dtype)
    if debug:
            print(f'\nAfter: {info}')
    
    return arrayOut



"""
Function:
    RegisterSlices()
    
Purpose:
    Apply registration to slices from a dictionary resulting from calling the
    function GetClosestSlices().

Input:
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
                        
    regMethod       - Type of registration to perform on the DICOM scan
                      (collection) that the ROI Collection is to be copied
                      to.  Acceptable values include:
                        -- 'rigid' = perform rigid registration
                        -- 'affine' = perform affine registration
                        -- 'non-rigid' = perform non-rigid registration
                        -- any other character string (including '') = do not
                            perform any registration
                        
    xnatDict        - Dictionary containing XNAT login details, REST path 
                      variables and directory paths
                      (see DownloadFilesForRoiCopy() for details)
                      
    debug           - Prints some useful results if True

    
Returns:
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- Registration method used
                      -- Registered closest Dicom SOP Instance UID
                      -- Registered closest Dicom filepath
                      -- Registered closest pixel array 
                      -- Registered closest Dicom 
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series 
    
    xnatDict        - Dictionary containing XNAT login details, REST path 
                      variables and directory paths
                      (see DownloadFilesForRoiCopy() for details)
                  
"""



def RegisterSlices(dicomDict1_to_2, regMethod, xnatDict, debug):
    
    # Import packages and functions:
    
    #import importlib
    #import DicomHelperFuncs
    #importlib.reload(DicomHelperFuncs)
    import os
    import SimpleITK as sitk
    import copy
    import pydicom
    import time
    import numpy as np
    from DicomHelperFuncs import ConvertNumpyDtype
    
    #print(regMethod)
    
    # If regMethod is not 'rigid', 'affine' or 'non-rigid', simply return an
    # empty array for imReg:
    if not regMethod in ['rigid', 'affine', 'non-rigid']:
        print('Not performing registration')
        return [], xnatDict
    
    else:
        """
        # Get the root directory, searchBy and scanPosCorr variables:
        rootDir = xnatDict['rootDir']
        searchBy = xnatDict['searchBy']
        scanPosCorr = xnatDict['scanPosCorr']
        
        # Store the registered DICOMs to a new directory.
        # Create a new variable for the directory name:
        dirName = ''
        
        # Start by adding the Scan Label:
        dirName = dirName + xnatDict['toDicomScanLab']
        
        # If scanPosCorr is non-zero:
        if scanPosCorr:
            dirName = dirName + ' ' + searchBy + 'PosCorr'
          
        # Add the registration method to be applied:
        dirName = dirName + ' ' + regMethod + 'Reg'
        
        # The registered DICOMs that relate to the new ROIs are to be stored at:
        toDicomRegDir = os.path.join(rootDir, dirName + ' DICOMs')
        
        # Update the dictionary xnatDict:
        xnatDict.update({'toDicomRegDir':toDicomRegDir})
        """
        
        # Get the directory to export the registered DICOMs:
        toRegDicomDir = xnatDict['toRegDicomDir']
        
        # Create the directory to store the registered DICOMs if they don't 
        # exist:
        try:
            os.mkdir(toRegDicomDir)
            if debug:
                print('Directory', toRegDicomDir, 'created')
        except FileExistsError:
            if debug:
                print('Directory', toRegDicomDir, 'already exists')
                
        
        # Create a dictionary to store the changes made to the DICOM of the  
        # registered slices:
        changesDict = {}
        
        #print('\ndicomDict1_to_2:\n', dicomDict1_to_2, '\n\n')
        
        # Get the keys (SOP UIDs) in dicomDict1_to_2:
        keys = list(dicomDict1_to_2.keys())
        K = len(keys)
        
        
        """
        The new DICOMS will need some new (unique) tags. New SOP Instance UIDs
        will be generated for each key in dicomDict1_to_2, in the for loop 
        below.  Other tags only need generating once, so will use the values 
        from the first key (keys[0]).
        """
        # Load the first DICOM in the set to be registered:
        dicom2 = dicomDict1_to_2[keys[0]]['Closest Dicom']
        
        # Add some details of the type of registration used in the Series
        # Description:
        preRegSeriesDesc = dicom2.SeriesDescription
        postRegSeriesDesc = preRegSeriesDesc + ' ' + regMethod + ' registered'
        
        # Add some details of the type of registration used in the Series
        # Number:
        preRegSeriesNo = dicom2.SeriesNumber
        postRegSeriesNo = str(preRegSeriesNo) + '-' + regMethod + 'Reg'
        """
        Maybe Series Number cannot have strings, so just change the number.
        """
        postRegSeriesNo = str(preRegSeriesNo) + '-2'
        """
        That was a problem too.
        """
        postRegSeriesNo = preRegSeriesNo + 10
        
        # Create a different series number for different registration methods
        # (so that multiple different registrations can be uploaded to XNAT):
        if 'affine' in regMethod:
            postRegSeriesNo = preRegSeriesNo + 100
        elif 'affine' in regMethod:
            postRegSeriesNo = preRegSeriesNo + 200
        else:
            postRegSeriesNo = preRegSeriesNo + 300
        
        # Create a new Study Instance UID:
        preRegStudyUID = dicom2.StudyInstanceUID
        postRegStudyUID = pydicom.uid.generate_uid()   
        
        # Create a new Series Instance UID:
        preRegSeriesUID = dicom2.SeriesInstanceUID
        postRegSeriesUID = pydicom.uid.generate_uid()  
        
        # Create a new Frame of Reference UID:
        preRegFrameOfRefUID = dicom2.FrameOfReferenceUID
        postRegFrameOfRefUID = pydicom.uid.generate_uid() 
        
        # Display changes:    
        print('\nChanging the Series Description from:\n', \
          preRegSeriesDesc, '\nto:\n', postRegSeriesDesc)
        print('\nChanging the Series Number from:\n', \
          preRegSeriesNo, '\nto:\n', postRegSeriesNo)
        print('\nChanging the Study Instance UID from:\n', \
          preRegStudyUID, '\nto:\n', postRegStudyUID)
        print('\nChanging the Series Instance UID from:\n', \
          preRegSeriesUID, '\nto:\n', postRegSeriesUID)
        print('\nChanging the Frame of Reference UID from:\n', \
          preRegFrameOfRefUID, '\nto:\n', postRegFrameOfRefUID)
            
        
        # Get the details of the pixel array of the first DICOM from the set 
        # to be registered (not necessary but useful for debugging):
        preRegPixelSpacing = dicom2.PixelSpacing
        #preRegRows = dicom2.Rows
        #preRegColumns = dicom2.Columns
        preRegPixelArray = dicom2.pixel_array
        
        print(f'\nPre-registered pixel array dims = {np.shape(preRegPixelArray)}')
        #print(f'Pre-registered no of rows = {preRegRows}')
        #print(f'Pre-registered no of columns = {preRegColumns}')
        print(f'\nPre-registered pixel spacing = {preRegPixelSpacing}')
        
        
        # Store the original and modified SOP Instance UIDs of the registered
        # sliced:
        preRegSOPUIDs = []
        postRegSOPUIDs = []
        
        # Keep track of the time it takes to execute each loop:
        times = []
        
        # times[0] = start time:
        times.append(time.time())
        
        # Iterate for each key:
        for k in range(K):
            key = keys[k]
            
            print('\nkey =', key)
            
            print(f'\nPerforming registration on slice {k+1}/{K}..')
            
            """
            Although the pixel arrays are already contained within 
            dicomDict1_to_2, it seems that SimpleITK needs to read and image 
            directly from a DICOM file.
            """
            # So get the filepaths:
            fpath1 = dicomDict1_to_2[key]['Dicom fpath']
            fpath2 = dicomDict1_to_2[key]['Closest Dicom fpath']
            
            print('\nfpath2 =', fpath2)
            
            """
            Before doing anything serious (like registrations) check to see
            if the files already exist (there may be reasons to run this 
            function but without over-writting the already registered images).
            """
            
            # Get current filename:
            preRegFname = os.path.basename(fpath2)
            
            # Get filename (with extension):
            preRegFname = os.path.splitext(preRegFname)
            
            # Create a new filename, adding some info about the registration
            # method used:
            postRegFname = preRegFname[0] + '-' + regMethod + '-registered.dcm'
            
            print('\n   Changing the file name from:\n', \
              preRegFname[0], '\nto:\n', postRegFname)
            
            # Create new filepath:
            postRegFpath = os.path.join(toRegDicomDir, postRegFname)
            
            
            """
            # Decide whether to over-write the registered files (if they 
            # already exist) or not:
            overWrite = False
            
            # Check if the registered DICOM file already exists:
            if os.path.isfile(postRegFpath):
                print('\nThe registered DICOM file already exists.')
                
                if not overWrite:
                    print('\nThe registered DICOM file will not be overwritten.')
                else:
                    print('\nThe registered DICOM file will be overwritten.')
            """
        
            
            # Read in the 3D images using SimpleITK:
            im3D1 = sitk.ReadImage(fpath1)
            im3D2 = sitk.ReadImage(fpath2)
            
            """
            SimpleITK treats the images as 3D so need to extract the 2D images.
            """
            # Extract the 2D images from the 3D images:
            """
            Note: When loading DICOMs with 512x512 pixel arrays the resulting
            2D images end up with 384x384 pixels for some unknown reason.
            If I replace (im3D1.GetWidth(), im3D1.GetHeight(), 0) with
            (512, 512, 0) I get the following error:
                RuntimeError: Exception thrown in SimpleITK Extract: 
                    C:\SimpleElastix\build_20200108\ITK\Modules\Core\Common\src\itkDataObject.cxx:393:
                    Requested region is (at least partially) outside the 
                    largest possible region.
            So for now will have to accept that the image resolutions will be
            different.
            """
            im1 = sitk.Extract(im3D1, (im3D1.GetWidth(), im3D1.GetHeight(), 0), (0, 0, 0))
            im2 = sitk.Extract(im3D2, (im3D2.GetWidth(), im3D2.GetHeight(), 0), (0, 0, 0))
            
            
            # Perform registration depending on regMethod:
            if regMethod in ['rigid', 'affine']:
                # Perform rigid or affine registration:
                imFilter = sitk.ElastixImageFilter()
                imFilter.SetFixedImage(im1)
                imFilter.SetMovingImage(im2)
                imFilter.SetParameterMap(sitk.GetDefaultParameterMap(regMethod))
                imReg = imFilter.Execute()
            
            else:
                # Perform non-rigid registration:
                imReg = sitk.Elastix(im1, im2)
                
                
            # Convert the SimpleITK image imReg to a numpy data array:
            imRegNda = sitk.GetArrayFromImage(imReg)
            
            # Since registered SimpleITK objects seem to return floats after
            # converting to numpy, convert to the same data type of the pre-
            # registered DICOM.
            
            # Convert the pre-registered image to numpty array:
            im2Nda = sitk.GetArrayFromImage(im2)
            
            # Get the data type:
            im2dtype = str(im2Nda.dtype)
            
            if debug:
                print('\nData type of pre-registered image was ', im2dtype)
            
            # Convert the numpy data array from float to im2dtype:
            imRegNdaConverted = ConvertNumpyDtype(imRegNda,\
                                                  convertTo=im2dtype,\
                                                  debug=debug)
                
            
            # Create a new DICOM to modify with the registered pixel array:
            dicomReg = copy.deepcopy(dicomDict1_to_2[key]['Closest Dicom'])
            
            # Change the pixel array to imRegNdaConverted:
            #dicomReg.pixel_array = imRegNda # <- AttributeError: can't set attribute
            #dicomReg['PixelArray'] = imRegNda # <-TypeError: Dataset contents must be DataElement instances
            #dicomReg.pixel_array = preRegPixelArray # <- AttributeError: can't set attribute
            #dicomReg.pixel_array = dicomReg.pixel_array # <- same error
            #dicomReg.PixelData = imRegNda.tostring()
            dicomReg.PixelData = imRegNdaConverted.tostring()
            
            """
            Since the array dimensions have changed, update the relevant 
            metadata
            """
            #dicomReg.Rows, dicomReg.Columns = imRegNda.shape
            dicomReg.Rows, dicomReg.Columns = imRegNdaConverted.shape
            
            # Get the original (pre-registered) SOP Instance UID:
            #preRegSOPUID = dicomDict1_to_2[key]['Closest Dicom SOP UID']
            preRegSOPUID = dicomReg.SOPInstanceUID
            
            # Append to postRegSOPUIDs:
            preRegSOPUIDs.append(preRegSOPUID)
            
            
            # Generate a new SOP Instance UID:
            """
            Note that there may be slices that are repeated in the second set 
            due to differences in slice thickness, e.g. slice no 1 in the first 
            set may be closest to slices 1 and 2 in the second set.  Hence the
            same slice in the second set will have a unique SOP Instance UID
            if it is repeated.  In any case, each slice will need a unique SOP
            Instance UID to avoid throwing up errors when it comes to uploading
            DICOMs to XNAT.
            """
            postRegSOPUID = pydicom.uid.generate_uid() 
            
            # Append to postRegSOPUIDs:
            postRegSOPUIDs.append(postRegSOPUID)
            
            print('\n   Changing the SOP Instance UID from:\n', \
              preRegSOPUID, '\nto:\n', postRegSOPUID)
            

        
            # Make changes to the tags:
            dicomReg.SOPInstanceUID = postRegSOPUID
            dicomReg.SeriesDescription = postRegSeriesDesc
            dicomReg.SeriesNumber = postRegSeriesNo
            dicomReg.StudyInstanceUID = postRegStudyUID
            dicomReg.SeriesInstanceUID = postRegSeriesUID
            dicomReg.FrameOfReferenceUID = postRegFrameOfRefUID
            
            
            """
            # Get current filename:
            preRegFname = os.path.basename(fpath2)
            
            # Get filename (with extension):
            preRegFname = os.path.splitext(preRegFname)
            
            # Create a new filename, adding some info about the registration
            # method used:
            postRegFname = preRegFname[0] + '-' + regMethod + '-registered.dcm'
            
            print('\n   Changing the file name from:\n', \
              preRegFname[0], '\nto:\n', postRegFname)
            
            # Create new filepath:
            postRegFpath = os.path.join(toDicomRegDir, postRegFname)
            """
            
            
            # Export the DICOM file:
            print('\nExporting the registered DICOM to:\n', postRegFpath)
            dicomReg.save_as(postRegFpath)

        
            # Update the dictionary dicomDict1_to_2:
            dicomDict1_to_2[key].update(
                    {
                    'Registration method':regMethod,\
                    'Registered closest Dicom SOP UID':postRegSOPUID,\
                    'Registered closest Dicom fpath':postRegFpath,\
                    #'Registered closest Pix array':imRegNda,\
                    'Registered closest Pix array':imRegNdaConverted,\
                    'Registered closest Dicom':dicomReg
                    }
                    )
            
            # times[k] = this loop time:
            times.append(time.time())
            
            # Time taken to complete this loop:
            dt = times[-1] - times[-2]
            
            # Total time taken so far:
            Dt = times[-1] - times[0]
            
            # Average time taken per loop so far:
            period = Dt/(k+1)
            
            # Number of iterations remaining:
            rem = K - (k+1)
            
            # Anticipated time remaining to completion:
            t_rem = period*rem
            
            print(f'\n   Took {round(dt, 1)} s for this registration.')
            if Dt < 60:
                print(f'\n   Total time so far {round(Dt, 1)} s.')
            else:
                print(f'\n   Total time so far {round(Dt/60, 1)} min.')
            if t_rem < 60:
                print(f'\n   Expected time to completion {round(t_rem, 1)} s.')
            else:
                print(f'\n   Expected time to completion {round(t_rem/60, 1)} min.')
            print('')
        
        
            
        # Display changes:    
        print('\nChanging the Series Description from:\n', \
          preRegSeriesDesc, '\nto:\n', postRegSeriesDesc)
        print('\nChanging the Series Number from:\n', \
          preRegSeriesNo, '\nto:\n', postRegSeriesNo)
        print('\nChanging the Study Instance UID from:\n', \
          preRegStudyUID, '\nto:\n', postRegStudyUID)
        print('\nChanging the Series Instance UID from:\n', \
          preRegSeriesUID, '\nto:\n', postRegSeriesUID)
        print('\nChanging the Frame of Reference UID from:\n', \
          preRegFrameOfRefUID, '\nto:\n', postRegFrameOfRefUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'preRegSeriesDescription':preRegSeriesDesc})
        changesDict.update({'postRegSeriesDescription':postRegSeriesDesc})
        changesDict.update({'preRegSeriesNumber':preRegSeriesNo})
        changesDict.update({'postRegSeriesNumber':postRegSeriesNo})
        changesDict.update({'preRegStudyInstanceUID':preRegStudyUID})
        changesDict.update({'postRegStudyInstanceUID':postRegStudyUID})
        changesDict.update({'preRegSeriesInstanceUID':preRegSeriesUID})
        changesDict.update({'postRegSeriesInstanceUID':postRegSeriesUID})
        changesDict.update({'preRegFrameOfReferenceUID':preRegFrameOfRefUID})
        changesDict.update({'postRegFrameOfReferenceUID':postRegFrameOfRefUID})
        changesDict.update({'preRegSOPInstanceUIDs':preRegSOPUIDs})
        changesDict.update({'postRegSOPInstanceUIDs':postRegSOPUIDs})
        
        return dicomDict1_to_2, xnatDict 




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
                        
    scanPosCorr        - Scan position correction that needs to be applied to
                        the DICOM scan (collection) to account for different
                        positional reference points (as is the case for 
                        different imaging modalities). The correction will be
                        applied to the same dimension as defined by searchBy.
                        Accepatable values are any float.  Note that if
                        searchBy='r' the correction will not be applied.
                        
    regMethod          - Type of registration to perform on the DICOM scan
                        (collection) that the ROI Collection is to be copied
                        to.  Acceptable values include:
                            -- 'rigid'
                            -- 'affine'
                            -- 'non-rigid'
                            -- 'none' (or 'None' or '')
    
    del_dicoms         - Delete downloaded DICOMs if True
    
    debug              - Print results if True

    
Returns:
    roi2               - ROI object that contains the contours in roi1 but 
                        modified for overlaying over the 2nd DICOM scan;
                       - A new DICOM-RTSTRUCT file will be saved to the location
                       specified
                   
    xnatDict           - Dictionary containing XNAT login details, REST path 
                        variables and directory paths
                        (see DownloadFilesForRoiCopy() for details)
                       - The dictionary is exported as a JSON file 
"""



def CopyRoi(fromRoiProjLab, fromRoiSubjLab, fromRoiExpLab, fromRoiLabDateTime,\
            fromDicomScanLab,\
            toDicomProjLab, toDicomSubjLab, toDicomExpLab, toDicomScanLab,\
            searchBy, scanPosCorr, regMethod, del_dicoms, debug):
    
    # IMPORT PACKAGES AND FUNCTIONS:
    
    #import importlib
    #import DicomHelperFuncs
    #importlib.reload(DicomHelperFuncs)
    from DicomHelperFuncs import DownloadFilesForRoiCopy
    from DicomHelperFuncs import GetDicomData
    #from DicomHelperFuncs import GetRoiData
    from DicomHelperFuncs import GetClosestSlices
    from DicomHelperFuncs import ModifyRoi
    from DicomHelperFuncs import RegisterSlices
    import time, os, shutil
    #import SimpleITK as sitk
    
    # Keep track of the time it takes to execute certain functions:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    # Get the current date:
    todaysDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    # Assign the root directory:
    rootDir = os.path.join(r'C:\Temp', todaysDate)
    
    """
    If registration is to be applied to the second DICOM scan, the ROIs from
    the first DICOM scan will be copied and modified so that they will overlay
    onto the registered series. To help make things clear, the exported 
    RTSTRUCT file will be stored in a folder with the e.g. directory structure:
        10 zPosCorr affineReg ROIs
        
    where 10 is the Series Number, 'zPosCorr' indicates that the Image 
    Positions of the second scan were corrected, and affineReg indicates that 
    an affine registration method was applied.
    
    If no position correction needs to be made, or no registration applied, the
    structure will be simpler, e.g.:
        10 ROIs
    """
    
    
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
                'toDicomScanLab':toDicomScanLab
                }
    
    # First define the folder structures of the DICOMs and DICOM-RTSTRUCT of
    # the scan whose ROIs are to be copied:
    
    # Where the DICOM-RTSTRUCT to be copied is to be saved to:
    fromRoiDir = os.path.join(rootDir, fromDicomScanLab + ' ROIs')
    
    # Where the DICOM scan that relate to the ROIs to be copies are to be saved 
    # to:
    fromDicomDir = os.path.join(rootDir, fromDicomScanLab + ' DICOMs')
    
    
    # Now deal with the folder structures for the DICOMs and DICOM-RTSTRUCT of
    # the scan that the ROIs are to be copied to:
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    #toDicomDir = os.path.join(rootDir, toDicomScanLab + ' DICOMs')
    
    # For the directory names to export the position-corrected DICOMs, the 
    # position-corrected-and-registered DICOMs, and the DICOM-RTSTRUCT for the
    # latter DICOMs, start by creating a empty character string that will be 
    # appended to:
    dirName = ''
    
    # Add the Scan Label:
    dirName = dirName + toDicomScanLab
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    toDicomDir = os.path.join(rootDir, toDicomScanLab + ' DICOMs')
    
    # If scanPosCorr is non-zero:
    if scanPosCorr:
        dirName = dirName + ' ' + searchBy + 'PosCorr'
      
    # The DICOMs of the series that the ROIs are to be copied to AFTER a change
    # to the Image Position but PRIOR to image registration:
    # If scanPosCorr is non-zero:
    if scanPosCorr:
        toPosCorrDicomDir = os.path.join(rootDir, dirName + ' DICOMs')
    else:
        toPosCorrDicomDir = ''
        
    # The DICOMs that relate to the new ROIs (which may be resulting from 
    # a registration) are to be stored at:
    # If registration is to be applied:
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        dirName = dirName + ' ' + regMethod + 'Reg'
        
        toRegDicomDir = os.path.join(rootDir, dirName + ' DICOMs')
    else:
        toRegDicomDir = ''
    
    # The new (modified copy) ROIs are to be stored at:
    toRoiDir = os.path.join(rootDir, dirName + ' ROIs')
    
    
    

    
    
    
    # Add the above variables + searchBy and scanPosCorr to xnatDict:
    xnatDict.update({
            # Define the directories to download files to:
            'rootDir':rootDir, \
            'fromDicomDir':fromDicomDir,\
            'fromRoiDir':fromRoiDir,\
            'toDicomDir':toDicomDir,\
            'toPosCorrDicomDir':toPosCorrDicomDir,\
            'toRegDicomDir':toRegDicomDir,\
            'toRoiDir':toRoiDir,\
                
            # Store some of the input variables:
            'searchBy':searchBy,
            'scanPosCorr':scanPosCorr
            })
        
    ## If registration is to be done, toDicomRegDir needs to be added too:
    #if regMethod in ['rigid', 'affine', 'non-rigid']:
    #    xnatDict.update({
    #            'toDicomRegDir':toDicomRegDir
    #            })
    #    
    ## Finally add the remaining variables:
    #xnatDict.update({
    #        'toRoiDir':toRoiDir,\
    #            
    #        # Store some of the input variables:
    #        'searchBy':searchBy,
    #        'scanPosCorr':scanPosCorr
    #        })
    
    #print('\nxnatDict:\n\n', xnatDict)
    
    # DOWNLOAD THE DICOM-RTSTRUCT AND DICOM FILES FROM XNAT:
    xnatDict = DownloadFilesForRoiCopy(xnatDict, debug=debug)
    
    # times[1] = time after downloading of files from XNAT:
    times.append(time.time())
    
    print(f'\nTook {round(times[-1] - times[-2], 1)} s to download data from ',\
                    'XNAT')
    
    # Parse the ROI Collection to be copied:
    """ This is not needed following changes to the helper functions """
    #roiDict1 = GetRoiData(xnatDict['fromRoiFpath'])


    # PARSE THE DICOM DATA THAT CORRESPONDS TO ROI COLLECTION TO BE COPIED:
    dicomDict1 = GetDicomData(xnatDict['fromDicomDir'], sortMethod='slices')
    
    
    # PARSE THE DICOM DATA THAT THE ROI COLLECTION IS TO BE COPIED TO:
    dicomDict2 = GetDicomData(xnatDict['toDicomDir'], sortMethod='slices')

    
    # GET THE MAPPING FROM THE 1st TO 2nd DICOM SCANS 
    # (using the z-coordinate to find the nearest slice): 
    dicomDict1_to_2, dicomDict2posCorr = GetClosestSlices(dicomDict1,\
                                                          dicomDict2,\
                                                          searchBy=searchBy,\
                                                          scanPosCorr=scanPosCorr,\
                                                          debug=debug)
    
    
    # Print some results for diagnostics purposes:
    """
    There are problems with this:
        KeyError: 'Dicom frame no'

    So for now commenting it out.
    """
    if False: #debug:
        from DicomHelperFuncs import IndexDictInTier2
        
        key1 = IndexDictInTier2(dicomDict1, tier2key='Dicom frame no',\
                                tier2val=5)


        print('Contour copied from Image Position =', \
              dicomDict1[key1]['Image pos'], '\n')
        
        # loop through each key in dicomDict2:
        #keys2 = list(dicomDict2.keys())
        keys2 = list(dicomDict2posCorr.keys())
        
        for i in range(len(keys2)):
            print(f'Slice {i+1} Image Position =', \
                  #dicomDict2[keys2[i]]['Image pos'])
                  dicomDict2posCorr[keys2[i]]['Image pos'])
            
    
    
    # APPLY REGISTRATION IF regMethod INDICATES TO:
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        dicomDict1_to_2, xnatDict = RegisterSlices(dicomDict1_to_2,\
                                                   regMethod=regMethod,\
                                                   xnatDict=xnatDict,\
                                                   debug=debug)
        
        # times[2] = time after performing registration:
        times.append(time.time())
        
        print(f'\nTook {round(times[-1] - times[-2], 1)} s to register images')
        
        
    # COPY AND MODIFY THE ROI COLLECTION 
    # ((from the original DICOM scan) so they may be overlaid on the new DICOM 
    # scan):
    #roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2)
    #roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2) # 10/03/2020
    roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2, regMethod) # 13/03/2020
    
    # times[3] = time after parsing files, finding closest slice locations 
    # and modifying the ROIs:
    times.append(time.time())
    
    print(f'\nTook {round(times[-1] - times[-2], 1)} s to copy and modify the ',\
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
                    
    
    print(f'\nTook {round(times[-1] - times[0], 1)} s to do everything')
    
    #return roi2, xnatDict
    #return dicomDict1, roi1, dicomDict2, roi2, dicomDict1_to_2, xnatDict # 10/03/2020
    return dicomDict1, roi1, dicomDict2posCorr, roi2, dicomDict1_to_2, xnatDict # 15/03/2020



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
    CombineDicomAndRoiDicts()
    
Purpose:
    Combine DICOM and ROI dictionaries.
    
Input:
    dicomDict  - dictionary of DICOM data
    
    roiDict    - dictionary of ROI data 
    
Returns:
    dicomAndRoiDict - dictionary containing the DICOM pixel arrays and contour 
                      data from the DICOM-RTSTRUCT file; the keys of the 
                      dictionary are the SOP Instance UIDs of the DICOMs (and 
                      the Reference SOP Instance UIDs of the ROIs)
"""

def CombineDicomAndRoiDicts(dicomDict, roiDict):
    
    # Import packages:
    #import numpy as np
    import copy
    #import os
    
    # Get the keys of the dictionaries dicomDict and roiDict:
    Ckeys = list(roiDict.keys())
    Dkeys = list(dicomDict.keys())
    
    # Create a new dictionary called dicomAndRoiDict that combines both 
    # dictionaries. 
    # Start by copying dicomDict:
    dicomAndRoiDict = copy.deepcopy(dicomDict)
    
    for Dkey in Dkeys: # cycle through each key in dicomDict
        # Check if Dkey is in Ckeys:
        if Dkey in Ckeys:
            # There is contour data for this Dkey, so update 
            # dicomAndRoiDict[Dkey] with roiDict[Ckey]:
            #ind = Ckeys.index(Dkey) # index of matching key in roiDict
            dicomAndRoiDict[Dkey].update(roiDict[Dkey]) # Note: Dkey=Ckey
        
        else:
            # There is no contour data for this Dkey, so fill in the 
            # missing values:
            dicomAndRoiDict[Dkey].update({'No of contours':0, \
                                          'Contour no':[], \
                                          'No of contour pts':0, \
                                          'Contour data':[]
                                          })
        
            
    # Create a new array "contourPts" containing the (x,y) image 
    # coordinates of each contour point and save it in dicomAndRoiDict:
    
    # Note:
    # Conversion of contour points from patient coordinates to image coordinates
    # adapted from James Petts' (ICR) Java code: 
    
    for uid, values in dicomAndRoiDict.items():
        # Get the values for this uid:
        # Note that for any given uid, NcontourPts and contourData could
        # be lists of values, if multiple contours existed for that uid.  
        # Even if there was only one contour, the single value (or single 
        # list) will be enclosed in [].
        position = values['Image pos'] # e.g. [x, y, z]
        orientation = values['Image orient'] # e.g. [row_x, row_y, row_z, col_x, col_y, col_z]
        pixSpacing = values['Pix spacing'] # e.g. [dx, dy]
        Ncontours = values['No of contours'] # e.g. 1 or 2, etc.
        contourData = values['Contour data'] # e.g. [[x1,y1,z1,...]] or [[x1,y1,z1,...], [x1,y1,z1,...]], etc.
    
        # The patient coordinates:
        x0 = position[0]
        y0 = position[1]
        z0 = position[2]
        
        # The direction cosines along rows and columns:
        #row_x, row_y, row_z, col_x, col_y, col_z = orientation
        theta_r = orientation[0:3]
        theta_c = orientation[3:6]
        
        # The indeces of the largest direction cosines along rows and columns:
        #print('theta_r =', theta_r)
        #print('theta_c =', theta_c)
        #ind_r = np.where(theta_r==max(theta_r))[0][0]
        #ind_c = np.where(theta_c==max(theta_c))[0][0]
        #ind_r = np.where(theta_r==max(theta_r))
        #ind_c = np.where(theta_c==max(theta_c))
        ind_r = theta_r.index(max(theta_r))
        ind_c = theta_c.index(max(theta_c))
        #print('ind_r =', ind_r)
        #print('ind_c =', ind_c)
        
        # The pixel spacings:
        dx = pixSpacing[0]
        dy = pixSpacing[1]
        
        # Create an array of (x,y) points from the flattened contour data:
        contourPts = []
        
        # Iterate for each contour for this uid:
        for c in range(Ncontours):
            # Store the contour points for each contour within any given uid:
            contourPts_c = []
        
            # Iterate for all countour points in threes:
            for p in range(0, len(contourData[c]), 3):
                # Define vector v as the the patient position subtracted from 
                # the contour coordinates:
                v = contourData[c][p] - x0, \
                contourData[c][p+1] - y0, \
                contourData[c][p+2] - z0
                
                # Convert from patient coord system to image coord system:          
                i = (v[ind_r] - theta_c[ind_r]*v[ind_c]/theta_c[ind_c]) / \
                (theta_r[ind_r]*dx * (1 - (theta_c[ind_r]*theta_r[ind_c])/(theta_r[ind_r]*theta_c[ind_c])))
                
                j = (v[ind_c] - theta_r[ind_c]*i*dx) / (theta_c[ind_c]*dy)
    
                contourPts_c.append([i,j])
                
            # Append contourPts with contourPts_c:
            contourPts.append(contourPts_c)
    
        # Add contourPts to dicomAndRoiDict for this uid (key):
        dicomAndRoiDict[uid].update({'Contour pts':contourPts})
    
    
    return dicomAndRoiDict
    
    

"""
Function:
    PlotDicomsAndRoisFromDicts()
    
Purpose:
    Plot DICOMS and ROIs from dictionaries.
    
Input:
    dicomAndRoiDict - dictionary of DICOM and ROI data combined
    
    restrictBy      - option to restrict plot to limited range, with example
                    acceptable values:
                        -- [] = plot all slices
                        -- ['slices', 4, 9] = plot slices 4 to 9 inclusive
                        -- ['z-position', -87.5, 21.5] = plot slices whose 
                        z-position are within the range [-87,21]
                        -- similarly replace 'z-position' with 'x-position' or 
                        'y-position'
    
    exportPlot - If =True export the plot
    
    exportDir  - path to directory where plot is to be exported to
    
Returns:
    dicomAndRoiDict - dictionary containing the DICOM pixel arrays and contour 
                      data from the DICOM-RTSTRUCT file; the keys of the 
                      dictionary are the SOP Instance UIDs of the DICOMs (and 
                      the Reference SOP Instance UIDs of the ROIs)
"""

def PlotDicomsAndRoisFromDicts(dicomAndRoiDict, \
                               restrictBy, \
                               exportPlot, \
                               exportDir):
    
    # Import packages:
    import numpy as np
    import matplotlib.pyplot as plt
    #import copy
    #import os

    
    # Get the keys of the dictionary dicomAndRoiDict:
    keys = list(dicomAndRoiDict.keys())
    
    
    # Need the determine the indices (inds) that fulfill the restrictBy 
    # requirement.  Start by getting the min, max values of the range to 
    # restrict by if restrictBy is not empty:
    if restrictBy!=[]:
        # Get the limits:
        limits = restrictBy[1::]
        
    # Get the indices (inds) either by restricting by slice no or position:
    if 'slices' in restrictBy[0]:
        # Create a list of slice numbers:
        allSlices = [dicomAndRoiDict[key]['Dicom slice no'] for key in keys]
        
        # Use the limits to get the desired slice numbers:
        someSlices = [i for i in allSlices if i>=limits[0] and i<=limits[1]]
        
        #print('allSlices =', allSlices)
        #print('someSlices =', someSlices)
        
        # Get the indeces where someSlices are contained within allSlices:
        inds = [allSlices.index(someSlices[i]) for i in range(len(someSlices))]
        
        #print('inds =', inds)
        print(f'\nThere are {len(inds)} slices within the range {limits}.')
        
        # Use inds to get the restricted keys:
        keys = keys[inds]
        
    elif 'position' in restrictBy[0]:
        #print('restrictBy[0] =', restrictBy[0][0])
        
        # Get the position index (0, 1 or 2 for x, y or z) depending on
        # what the first character in restrictBy[0] is:
        if restrictBy[0][0]=='x':
            posInd = 0
        elif restrictBy[0][0]=='y':
            posInd = 1
        elif restrictBy[0][0]=='z':
            posInd = 2
        
        # Create a list of positions:
        allPos = [dicomAndRoiDict[key]['Image pos'][posInd] for key in keys]
        
        # Use the limits to get the desired positions:
        somePos = [i for i in allPos if i>=limits[0] and i<=limits[1]]
        
        #print('\nallPos =', allPos)
        #print('\nsomePos =', somePos)
        
        # Get the indeces where somePos are contained within allPos:
        inds = [allPos.index(somePos[i]) for i in range(len(somePos))]
        
        #print('\ninds =', inds)
        print(f'\nThere are {len(inds)} positions within the range', \
              restrictBy[0][0], f'= {limits}.')
        
        # Use inds to get the restricted keys:
        #keys = keys[inds] # AttributeError: 'NoneType' object has no attribute 
        # 'keys'
        keys = [keys[i] for i in inds]


    
    
    # Plot the DICOM scans with contours overlayed:
    # Configure plot:

    # Set the title font:
    fontSize=11
    
    Nuids = len(keys)
    
    # Set the number of subplot rows and columns:
    rows = np.int8(np.ceil(np.sqrt(Nuids)))
    cols = np.int8(np.ceil(np.sqrt(Nuids)))
    
    # Limit the number of rows to 4 otherwise it's difficult to make sense of
    # the images:
    if cols>4:
        cols = 4
        rows = np.int8(np.ceil(Nuids/cols))
    
    # Decide whether to flip the images up-down so that the orientation is
    # the same as shown in the OHIF-Viewer:
    #flip = False
    flip = True
    
    # Plot results:
    
    #plt.figure(figsize=(15,15), dpi=300);
    plt.figure(figsize=(15,20), dpi=300);
    #plt.figure(dpi=300);
    
    i = 0 # for subplot pos
     
    # Look for each key in keys:
    for key in keys:
        values = dicomAndRoiDict[key]
        
        sliceNo = values['Dicom slice no']
        pixArray = values['Pix array']
        Ncontours = values['No of contours']
        contourPts = values['Contour pts']
        zPos = values['Image pos'][2]
        
        # Plot the pixel array:
        i = i + 1    
        plt.subplot(rows,cols,i, aspect='equal')
        if flip:
            plt.pcolormesh(np.flipud(pixArray));
        else:
            plt.pcolormesh(pixArray);
        
        plt.title(f'Slice {sliceNo} at z = {zPos}', size=fontSize);
        plt.axis('off');
        
        # Plot the contours:
        if Ncontours > 0:   
            for c in range(Ncontours):
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
                for x, y in contourPts[c]:
                    X.append(x)
                    Y.append(y)
    
                # Create 2D arrays for the X and Y values using meshgrid:
                #X, Y = np.meshgrid(xList, yList)
                
                # Flip Y pixel coordinates if flip=True:
                if flip:
                    N_r, N_c = np.shape(pixArray)
                    
                    Y = N_c - np.array(Y)
    
                plt.plot(X, Y, linewidth=0.5, c='red');
    
    
    
    # 10./03/20: This needs to be changes since dicom doesn't exist...
    """
    if exportPlot:
        # Set up details for exporting of plot:
        #seriesDescription = dicom.SeriesDescription
        
        
        print('Study description =', dicom.StudyDescription)
        print('Patient ID =', dicom.PatientID)
        print('Series description =', dicom.SeriesDescription)
        
        
        exportFname = dicom.StudyDescription + ' ' + dicom.PatientID + ' ' \
        + dicom.SeriesDescription + '.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        # Export the plot:
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print('Plot exported to:\n\n', exportFpath)
    """
    
    return 
    
    
    
    
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
    
    
"""
Function:
    GetDicomFpaths_OLD()
    
Purpose:
    To get list of full file paths for all DICOM files in a given directory
    sorted using natsort.
    
    11/03/20 Note: Even though the DICOMs are sorted in a natural way (e.g.
    0000.dcm, 0001.dcm, 0002.dcm) it doesn't guarantee that the files are
    sorted in the sequence of the slices!  So rather than relying on the 
    file names I will need to use the Image Position (Patient) to sort. This
    will be done in a separate function.
    

Input:
    dirPath   - path to directory containing DICOM files
    
Returns:
    filePaths - full file paths of all DICOM files in dirPath sorted in a 
                natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ... rather than
                3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
"""


def GetDicomFpaths_OLD(dirPath):
    
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
  
    

"""
Function:
    PlotDicomsAndRoisFromDicts_OLD()
    
Purpose:
    Plot DICOMS and ROIs from dictionaries.
    
Input:
    dicomDict  - dictionary of DICOM data
    
    roiDict    - dictionary of ROI data 
    
    exportPlot - If =True export the plot
    
    exportDir  - path to directory where plot is to be exported to
    
Returns:
    dicomAndRoiDict - dictionary containing the DICOM pixel arrays and contour 
                      data from the DICOM-RTSTRUCT file; the keys of the 
                      dictionary are the SOP Instance UIDs of the DICOMs (and 
                      the Reference SOP Instance UIDs of the ROIs)
"""

def PlotDicomsAndRoisFromDicts_OLD(dicomDict, roiDict, exportPlot, exportDir):
    
    # Import packages:
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    import os
    
    # Get the keys of the dictionaries dicomDict and roiDict:
    Ckeys = list(roiDict.keys())
    Dkeys = list(dicomDict.keys())
    
    # Create a new dictionary called dicomAndRoiDict that combines both 
    # dictionaries. 
    # Start by copying dicomDict:
    dicomAndRoiDict = copy.deepcopy(dicomDict)
    
    for Dkey in Dkeys: # cycle through each key in dicomDict
        # Check if Dkey is in Ckeys:
        if Dkey in Ckeys:
            # There is contour data for this Dkey, so update 
            # dicomAndRoiDict[Dkey] with roiDict[Ckey]:
            ind = Ckeys.index(Dkey) # index of matching key in roiDict
            dicomAndRoiDict[Dkey].update(roiDict[Dkey]) # Note: Dkey=Ckey
        
        else:
            # There is no contour data for this Dkey, so fill in the 
            # missing values:
            dicomAndRoiDict[Dkey].update({'No of contours':0, \
                                          'Contour no':[], \
                                          'No of contour pts':0, \
                                          'Contour data':[]
                                          })
        
            
    # Create a new array "contourPts" containing the (x,y) image 
    # coordinates of each contour point and save it in dicomAndRoiDict:
    
    # Note:
    # Conversion of contour points from patient coordinates to image coordinates
    # adapted from James Petts' (ICR) Java code: 
    
    for uid, values in dicomAndRoiDict.items():
        # Get the values for this uid:
        # Note that for any given uid, NcontourPts and contourData could
        # be lists of values, if multiple contours existed for that uid.  
        # Even if there was only one contour, the single value (or single 
        # list) will be enclosed in [].
        position = values['Image pos'] # e.g. [x, y, z]
        orientation = values['Image orient'] # e.g. [row_x, row_y, row_z, col_x, col_y, col_z]
        pixSpacing = values['Pix spacing'] # e.g. [dx, dy]
        Ncontours = values['No of contours'] # e.g. 1 or 2, etc.
        contourData = values['Contour data'] # e.g. [[x1,y1,z1,...]] or [[x1,y1,z1,...], [x1,y1,z1,...]], etc.
    
        # The patient coordinates:
        x0 = position[0]
        y0 = position[1]
        z0 = position[2]
        
        # The direction cosines along rows and columns:
        #row_x, row_y, row_z, col_x, col_y, col_z = orientation
        theta_r = orientation[0:3]
        theta_c = orientation[3:6]
        
        # The indeces of the largest direction cosines along rows and columns:
        #print('theta_r =', theta_r)
        #print('theta_c =', theta_c)
        #ind_r = np.where(theta_r==max(theta_r))[0][0]
        #ind_c = np.where(theta_c==max(theta_c))[0][0]
        #ind_r = np.where(theta_r==max(theta_r))
        #ind_c = np.where(theta_c==max(theta_c))
        ind_r = theta_r.index(max(theta_r))
        ind_c = theta_c.index(max(theta_c))
        #print('ind_r =', ind_r)
        #print('ind_c =', ind_c)
        
        # The pixel spacings:
        dx = pixSpacing[0]
        dy = pixSpacing[1]
        
        # Create an array of (x,y) points from the flattened contour data:
        contourPts = []
        
        # Iterate for each contour for this uid:
        for c in range(Ncontours):
            # Store the contour points for each contour within any given uid:
            contourPts_c = []
        
            # Iterate for all countour points in threes:
            for p in range(0, len(contourData[c]), 3):
                # Define vector v as the the patient position subtracted from 
                # the contour coordinates:
                v = contourData[c][p] - x0, \
                contourData[c][p+1] - y0, \
                contourData[c][p+2] - z0
                
                # Convert from patient coord system to image coord system:          
                i = (v[ind_r] - theta_c[ind_r]*v[ind_c]/theta_c[ind_c]) / \
                (theta_r[ind_r]*dx * (1 - (theta_c[ind_r]*theta_r[ind_c])/(theta_r[ind_r]*theta_c[ind_c])))
                
                j = (v[ind_c] - theta_r[ind_c]*i*dx) / (theta_c[ind_c]*dy)
    
                contourPts_c.append([i,j])
                
            # Append contourPts with contourPts_c:
            contourPts.append(contourPts_c)
    
        # Add contourPts to dicomAndRoiDict for this uid (key):
        dicomAndRoiDict[uid].update({'Contour pts':contourPts})
    
    
    #return dicomAndRoiDict

#"""
    
    # Plot the DICOM scans with contours overlayed:
    # Configure plot:

    # Set the title font:
    fontSize=4
    
    Nuids = len(dicomAndRoiDict.keys())
    
    # Set the number of subplot rows and columns:
    rows = np.int8(np.ceil(np.sqrt(Nuids)))
    cols = np.int8(np.ceil(np.sqrt(Nuids)))
    
    # Limit the number of rows to 5 otherwise it's difficult to make sense of
    # the images:
    if cols>5:
        cols = 5
        rows = np.int8(np.ceil(Nuids/cols))
    
    # Decide whether to flip the images up-down so that the orientation is
    # the same as shown in the OHIF-Viewer:
    #flip = False
    flip = True
    
    # Plot results:
    
    #plt.figure(figsize=(15,15), dpi=300);
    plt.figure(figsize=(15,30), dpi=300);
    #plt.figure(dpi=300);
    
    i = 0 # for subplot pos
    
    for uid, values in dicomAndRoiDict.items(): 
        frameNo = values['Dicom frame no']
        pixArray = values['Pix array']
        Ncontours = values['No of contours']
        contourPts = values['Contour pts']
        zPos = values['Image pos'][2]
        
        # Plot the pixel array:
        i = i + 1    
        plt.subplot(rows,cols,i, aspect='equal')
        if flip:
            plt.pcolormesh(np.flipud(pixArray));
        else:
            plt.pcolormesh(pixArray);
        
        plt.title(f'Slice {frameNo} at z = {zPos}', size=fontSize);
        plt.axis('off');
        
        # Plot the contours:
        if Ncontours > 0:   
            for c in range(Ncontours):
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
                for x, y in contourPts[c]:
                    X.append(x)
                    Y.append(y)
    
                # Create 2D arrays for the X and Y values using meshgrid:
                #X, Y = np.meshgrid(xList, yList)
                
                # Flip Y pixel coordinates if flip=True:
                if flip:
                    N_r, N_c = np.shape(pixArray)
                    
                    Y = N_c - np.array(Y)
    
                plt.plot(X, Y, linewidth=0.5, c='red');
    
    
    
    # 10./03/20: This needs to be changes since dicom doesn't exist...
    if exportPlot:
        # Set up details for exporting of plot:
        #seriesDescription = dicom.SeriesDescription
        
        
        print('Study description =', dicom.StudyDescription)
        print('Patient ID =', dicom.PatientID)
        print('Series description =', dicom.SeriesDescription)
        
        
        exportFname = dicom.StudyDescription + ' ' + dicom.PatientID + ' ' \
        + dicom.SeriesDescription + '.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        # Export the plot:
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print('Plot exported to:\n\n', exportFpath)
    
    return dicomAndRoiDict
    
#"""    
    

"""
Function:
    CopyRoi_OLD()
    
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



def CopyRoi_OLD(fromRoiProjLab, fromRoiSubjLab, fromRoiExpLab, fromRoiLabDateTime,\
            fromDicomScanLab,\
            toDicomProjLab, toDicomSubjLab, toDicomExpLab, toDicomScanLab,\
            searchBy, del_dicoms, debug):
    
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
    dicomDict1 = GetDicomData(xnatDict['fromDicomDir'], sortMethod='slices')
    
    
    # Parse the DICOM data that the ROI Collection is to be copied to:
    dicomDict2 = GetDicomData(xnatDict['toDicomDir'], sortMethod='slices')

    
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
    #roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2)
    roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2) # 10/03/2020
    
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
    
    #return roi2, xnatDict
    return dicomDict1, roi1, dicomDict2, roi2, dicomDict1_to_2, xnatDict # 10/03/2020



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
                                      'Closest Dicom slice no':\
                                      dicomDict2[keys2[ind]]['Dicom slice no'],\
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
                                      'r2 - r1':dR[ind], \
                                      'Closest Pix array':\
                                      dicomDict2[keys2[ind]]['Pix array'],\
                                      'Closest Dicom':\
                                      dicomDict2[keys2[ind]]['Dicom'],\
                                     })
            
    return dicomDict1_to_2



"""
Function:
    ModifyRoi_OLD()
    
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
                   
    xnatDict       - Dictionary containing XNAT login details, REST path 
                    variables and directory paths
                    (see DownloadFilesForRoiCopy() for details)
                   - The dictionary is exported as a JSON file
    
"""

def ModifyRoi_OLD(xnatDict, dicomDict1_to_2):
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
        changesDict.update({'fromFrameOfReferenceUID':roi2.FrameOfReferenceUID})
        
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
    
    print('\nChanging the Structure Set Time from:\n', roi2.StructureSetTime,\
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
        sliceNo1 = dicomDict1_to_2[D1to2keys[ind]]['Dicom slice no']
        sliceNo2 = dicomDict1_to_2[D1to2keys[ind]]['Closest Dicom slice no']
        
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
    #return roi2, xnatDict
    return roi1, roi2, xnatDict # 10/03/2020
    

