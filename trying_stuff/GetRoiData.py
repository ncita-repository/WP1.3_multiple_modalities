# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:46:39 2020

@author: ctorti
"""


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
