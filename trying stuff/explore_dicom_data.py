# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:10:50 2019

@author: Cristiano Torti


Basic functions to load DICOM data and output a 
selection of elements of interest (defined below).

The DICOM elements that are of interest (for now) are as follows:
- Series Instance UID
- Instance Creation Date
- Instance Create Time
- Modality
- Study ID
- Series Number
- Instance Number
- Image Position (Patient)
- Image Orientation (Patient)
- Frame of Reference UID


Revisions:

05/12/2019:
Revisions section added to comments section while in new branch.

"""


def load_dicoms_natsort(rootDir):
    # This function loads DICOM data from rootDir.
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Output: data = list of DICOM data for each DICOM file found in rootDir; 
    #               each row in data for a different .dcm file.
    
    # Import packages:
    import pydicom 
    import os
    import natsort
    
    # Generate list of DICOM file names (fnames), base directories (fdirs),
    # and full file paths (fpaths):
    # (Note: Only fpaths is used but others may be useful later..)
    
    Nfiles = 0
    fpaths = [] # (full) file paths
    fdirs = [] # file directories (base directories) for each path in fpaths
    fnames = [] # file names for each path in fpaths
    
    for root, dirs, files in os.walk(rootDir):
        for file in files:
    
            fpath = os.path.join(root, file)
    
            filepath_no_ext, file_ext = os.path.splitext(fpath)
    
            if 'dcm' in file_ext:
                fpaths.append(fpath) # append file path to fpaths
    
                dir_basename = os.path.basename(root) # get the dir basename
                fdirs.append(dir_basename) # append the directory name to fdirs
    
                fname = os.path.split(fpath)[1] # get the file name only
                fnames.append(fname) # append file name to fnames
        
    
    # Sort files in natural way:
    fpaths = natsort.natsorted(fpaths)
    
     # Read in all DICOM data using fpaths:
    data = []
    
    for fpath in fpaths:
        #data.append(pydicom.read_file(fpath))  
        data.append(pydicom.dcmread(fpath))  
        
    return data



def load_dicom(rootDir):
    # This function loads DICOM data from rootDir.
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Output: data = list of DICOM data for each DICOM file found in rootDir; 
    #               each row in data for a different .dcm file.
    
    # Import packages:
    import pydicom 
    import os
    
    # Generate list of DICOM file names (fnames), base directories (fdirs),
    # and full file paths (fpaths):
    # (Note: Only fpaths is used but others may be useful later..)
    
    Nfiles = 0
    fpaths = [] # (full) file paths
    fdirs = [] # file directories (base directories) for each path in fpaths
    fnames = [] # file names for each path in fpaths
    
    for root, dirs, files in os.walk(rootDir):
        for file in files:
    
            fpath = os.path.join(root, file)
    
            filepath_no_ext, file_ext = os.path.splitext(fpath)
    
            if 'dcm' in file_ext:
                fpaths.append(fpath) # append file path to fpaths
    
                dir_basename = os.path.basename(root) # get the dir basename
                fdirs.append(dir_basename) # append the directory name to fdirs
    
                fname = os.path.split(fpath)[1] # get the file name only
                fnames.append(fname) # append file name to fnames
        
    
     # Read in all DICOM data using fpaths:
    data = []
    
    for fpath in fpaths:
        data.append(pydicom.read_file(fpath))  
        
    return data



def list_dicom_EOI(rootDir):
    # This function creates a list of lists of some select elements in a list 
    # of DICOM data.  The DICOM elements that are of interest (for now) are as 
    # follows:
    # - Series Instance UID
    # - Instance Creation Date
    # - Instance Create Time
    # - Modality
    # - Study ID
    # - Series Number
    # - Instance Number
    # - Image Position (Patient)
    # - Image Orientation (Patient)
    # - Frame of Reference UID
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Ouput: dataList = list of select elements from data
    
    # Import packages:
    import pydicom 
    
    # Load the DICOM data:
    data = load_dicom(rootDir)
    
    dataList = [] # initialise list

    # Loop for each data set:
    for f in range(len(data)):
        thisSeries = [] # initialise series for this f
        thisSeries.append(data[f].SeriesInstanceUID)
        thisSeries.append(data[f].InstanceCreationDate)
        thisSeries.append(data[f].InstanceCreationTime)
        thisSeries.append(data[f].Modality)
        thisSeries.append(data[f].StudyID)
        thisSeries.append(data[f].SeriesNumber)
        thisSeries.append(data[f].InstanceNumber)
        thisSeries.append(data[f].ImagePositionPatient)
        thisSeries.append(data[f].ImageOrientationPatient)
        thisSeries.append(data[f].FrameOfReferenceUID)
        
        dataList.append(thisSeries) # append thisSeries to dataList

    return dataList



def list_SIUID_from_dicom_data(data):
    # This function creates a list of SeriesInstanceUIDs from a list of DICOM 
    # data.
    
    # Input: data = list of DICOM data
    
    # Output: SIUIDs = list of SeriesInstanceUIDs from data
    
    # Import packages:
    import pydicom
    
    SIUIDs = [] # initialise list

    for f in range(len(data)):
        SIUIDs.append(data[f].SeriesInstanceUID)
        
    return SIUIDs


def list_dicom_EOI_groupby_SIUID(rootDir):
    # This function creates a list of lists of some select elements in a list 
    # of DICOM data, grouped by SeriesInstanceUID. The DICOM elements that are 
    # of interest (for now) are as follows:
    # - Series Instance UID
    # - Instance Creation Date
    # - Instance Create Time
    # - Modality
    # - Study ID
    # - Series Number
    # - Instance Number
    # - Image Position (Patient)
    # - Image Orientation (Patient)
    # - Frame of Reference UID
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Ouput: dataListGroupBySIUID = list of select elements from data grouped 
    #           by common SeriesInstanceUID
    
    # Import packages:
    import pydicom 
    
    # Load the DICOM data:
    data = load_dicom(rootDir)
    
    # Get list of SIUIDs:
    SIUIDs = list_SIUID_from_dicom_data(data)
    
    # Get list of unique SeriesInstanceUIDs:
    USIUIDs = list(set(SIUIDs))
        
    dataListGroupBySIUID = [] # initialise list

    # Loop for each unique SeriesInstanceUID:
    for u in range(len(USIUIDs)):
        USIUID = USIUIDs[u] # this unique SeriesInstanceUID
        
        ICDs = [] # initialise list of InstanceCreationDate
        ICTs = [] # initialise list of InstanceCreationTime
        Ms = [] # initialise list of Modality
        SIDs = [] # initialise list of StudyID
        SNs = [] # initialise list of SeriesNumber
        INs = [] # initialise list of InstanceNumber
        IPPs = [] # initialise list of ImagePositionPatient
        IOPs = [] # initialise list of ImageOrientationPatient
        FORUIDs = [] # initialise list of FrameOfReferenceUID
        
        # Loop for each data set:
        for f in range(len(data)):
            # If the SeriesInstanceUID for this data[f] = 
            # this unique SeriesInstanceUID (USIUID):
            if data[f].SeriesInstanceUID == USIUID:
                ICDs.append(data[f].InstanceCreationDate)
                ICTs.append(data[f].InstanceCreationTime)
                Ms.append(data[f].Modality)
                SIDs.append(data[f].StudyID)
                ICTs.append(data[f].SeriesNumber)
                INs.append(data[f].InstanceNumber)
                IPPs.append(data[f].ImagePositionPatient)
                IOPs.append(data[f].ImageOrientationPatient)
                FORUIDs.append(data[f].FrameOfReferenceUID)
                
        # Append items to dataList for this unique SeriesInstanceUID: 
        dataThisUID = [USIUID, ICDs, ICTs, Ms, SIDs, ICTs, INs, IPPs, \
                       IOPs, FORUIDs]
        
        dataListGroupBySIUID.append(dataThisUID)

    return dataListGroupBySIUID



def display_dataframe_dicom_EOI(rootDir):
    # This function will display a Pandas dataframe of select elements from
    # a list of DICOM data. The DICOM elements that are of interest (for now) 
    # are as follows:
    # - Series Instance UID
    # - Instance Creation Date
    # - Instance Create Time
    # - Modality
    # - Study ID
    # - Series Number
    # - Instance Number
    # - Image Position (Patient)
    # - Image Orientation (Patient)
    # - Frame of Reference UID
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Ouput: df = Pandas dataframe of select elements from data
    
    # Import packages:
    import pandas
       
    # Get the dataList:
    dataList = list_dicom_EOI(rootDir)
    
    # Create a dataframe for easy visualisation:
    headers = ['Series Instance UID', 'Instance Creation Date', \
               'Instance Create Time', 'Modality', 'Study ID', \
               'Series Number', 'Instance Number', \
               'Image Position (Patient)', 'Image Orientation (Patient)', \
               'Frame of Reference UID']
    df = pandas.DataFrame(dataList, columns=headers)
    
    return df
    


def display_dataframe_dicom_EOI_groupby_SIUID(rootDir):
    # This function will display a Pandas dataframe of select elements from
    # a list of DICOM data grouped by SeriesInstanceUID. The DICOM elements 
    # that are of interest (for now) are as follows:
    # - Series Instance UID
    # - Instance Creation Date
    # - Instance Create Time
    # - Modality
    # - Study ID
    # - Series Number
    # - Instance Number
    # - Image Position (Patient)
    # - Image Orientation (Patient)
    # - Frame of Reference UID
    
    # Input: rootDir = filepath of root folder where DICOM data is located
    
    # Ouput: df = Pandas dataframe of select elements from data grouped by
    #           SeriesInstanceUID
    
    # Import packages:
    import pandas
       
    # Get the dataList:
    dataList = list_dicom_EOI_groupby_SIUID(rootDir)
    
    # Create a dataframe for easy visualisation:
    headers = ['Series Instance UID', 'Instance Creation Date', \
               'Instance Create Time', 'Modality', 'Study ID', \
               'Series Number', 'Instance Number', \
               'Image Position (Patient)', 'Image Orientation (Patient)', \
               'Frame of Reference UID']
    df = pandas.DataFrame(dataList, columns=headers)
    
    return df
    