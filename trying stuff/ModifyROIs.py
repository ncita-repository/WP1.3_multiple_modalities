# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:02:30 2020

@author: ctorti
"""


""" Function:
    
    ModifyROIs()
    
Purpose:
    Modify ROIs created for one series of DICOMs so they can be overlaid on
    a different series of DICOMs.  
    
Details:
    The function requires that the DICOMs and DICOM-RTSTRUCT files are already 
    stored on disc.  
    
    The function "CopyROIs" calls this function so that only the REST variables 
    of:
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
    origDicomsDir - directory containing original series of DICOM files
    
    origRoiPath   - file path of original DICOM-RTSTRUCT file that corresponds 
                    to the original DICOM series
    
    newDicomsDir  - directory containing a new series of DICOM files
    
    newRoiDir     - directory where the new DICOM-RTSTRUCT file is to be
                    exported to
    
    
Returns:
    newRoiFpath - the full filepath of a new DICOM-RTSTRUCT file in the 
                  directory specified by newRoiDir and a new file name 
                  generated with the current date and time
"""



#import numpy as np
import pydicom
#import DicomHelperFuncs
from DicomHelperFuncs import GetDicomFpaths
import copy
import os, time

def ModifyROIs(origDicomsDir, origRoiFpath, newDicomsDir, newRoiDir):
    # The filepaths of the original DICOMs:
    fpaths1 = GetDicomFpaths(origDicomsDir)
    
    # The ROIs that belong to the DICOMs in fpaths1 is:
    roi1 = pydicom.dcmread(origRoiFpath)
    
    # The filepaths of the new DICOMs:
    fpaths2 = GetDicomFpaths(newDicomsDir)
    
    
    # Load each DICOM in fpaths1, reading in the SOP Instance UIDs and store
    # them:
    sopuid1 = [] # initialise array
    
    for f in range(len(fpaths1)):
        fpath1 = fpaths1[f]
        
        # Read in the DICOMs:
        dicom1 = pydicom.read_file(fpath1)
        
        # Get the SOP Instance UIDs and store them in sopuid1:
        sopuid1.append(dicom1.SOPInstanceUID)
    
    
    # Since the number (n) of images (DICOMs) that have one or more contours 
    # is not necessarily the same as the total number (m) of contours, store 
    # the indeces that relate each "Contour Image Sequence"
    # (i.e. each image that has one or more contours) to the file 
    # number of the DICOM they relate to; and likewise, the indeces that relate 
    # each "Contour Sequence" (i.e. each individual contour) to the 
    # file number of the DICOM it relates to.
    
    # First iterate through each Contour Image Sequence and find the index of 
    # the matching SOP Instance UID from the DICOM:
    
    imageInds = [] # Contour Image Sequence to DICOM file number
    
    # Loop through each Contour Image Sequence:
    for sequence in roi1.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence:
        # Get the Reference SOP Instance UID:
        uid = sequence.ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in sopuid1:
        imageInds.append(sopuid1.index(uid))
            
    # Iterate through each Contour Sequence and find the index of the 
    # matching SOP Instance UID from the DICOM:
    
    contourInds = [] # Contour Sequence to DICOM file number
    
    # Loop through each Contour Sequence:
    for sequence in roi1.ROIContourSequence[0].ContourSequence:
        # Get the Reference SOP Instance UID:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
        
        # Find the matching index of this UID in sopuid1:
        contourInds.append(sopuid1.index(uid))
    
    
    """ 
    Start by making a copy of roi1, roi2, then:
    
    ----------------------------------------------------------------------------------------------
    Replace in roi2:                                           The value of the tag in the DICOM:
    ----------------------------------------------------------------------------------------------
    Study Instance UID (0020, 000d)                            Study Instance UID (0020, 000d)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Study Instance UID (0020, 000d)
    -> RT Referenced Study Sequence[0] (3006, 0012)
      -> Referenced SOP Instance UID (0008, 1155)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Series Instance UID (0020, 000e)
    -> RT Referenced Study Sequence[0] (3006, 0012)
      -> RT Referenced Series Sequence[0] (3006, 0014)
        -> Series Instance UID (0020, 000e)
    
    Frame of Reference UID (0020, 0052)                        Frame of Reference UID (0020, 0052)   
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Frame of Reference UID (0020, 0052)   
    -> Frame of Reference UID (0020, 0052)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     Frame of Reference UID (0020, 0052)
    -> Structure Set ROI Sequence[0] (3006, 0020)
      -> Referenced Frame of Reference UID (3006, 0024)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     SOP Instance UID (0008, 0018)
    -> RT Referenced Study Sequence[0] (3006, 0012) 
      -> RT Referenced Series Sequence[0] (3006, 0014)
        -> Contour Image Sequence[i] (3006, 0016)
          -> Referenced SOP Instance UID (0008, 1155)
    
    Referenced Frame of Reference Sequence[0] (3006, 0010)     SOP Instance UID (0008, 0018)
    -> ROI Contour Sequence[0] (3006, 0039)
      -> Contour Sequence[i] (3006, 0040)
        -> Contour Image Sequence[0] (3006, 0016)
          -> Referenced SOP Instance UID (0008, 1155)
    
    where the DICOM is the corresponding DICOM file number in the second series indexed by
    the indeces that link the DICOMs in the first series to its RTSTRUCT file.
    
    Note:  Even though the Study Instance and Frame of Reference UIDs in the ROI from the first
    DICOM series match the second DICOM series I'll include them below to keep this code as
    general as possible.
    """
    
    # Use roi1 as a template for roi2:
    roi2 = copy.deepcopy(roi1)
    
    # First modify the tags that are the same for all Contour or Contour Image 
    # Sequences:
    
    # Load the first DICOM from the new series:
    dicom2 = pydicom.read_file(fpaths2[0])
    
    # Change the Study Date and Time:
    roi2.StudyDate = copy.deepcopy(dicom2.StudyDate)
    #roi2.StudyTime = copy.deepcopy(dicom2.StudyTime) # <-- some DICOM Study
    # Times seem to have a decimal with trailing zeros,e.g. '164104.000000',
    # so convert to float, then integer, then back to string:
    roi2.StudyTime = str(int(float(copy.deepcopy(dicom2.StudyTime))))
    
    # Change the Patient's Name, ID, Birth Date and Sex:
    roi2.PatientName = copy.deepcopy(dicom2.PatientName)
    roi2.PatientID = copy.deepcopy(dicom2.PatientID)
    roi2.PatientBirthDate = copy.deepcopy(dicom2.PatientBirthDate)
    roi2.PatientSex = copy.deepcopy(dicom2.PatientSex)
    
    # Change the Study Instance UIDs:
    roi2.StudyInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
    roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID = copy.deepcopy(dicom2.StudyInstanceUID)
    
    # Change the Series Instance UID:
    roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID = copy.deepcopy(dicom2.SeriesInstanceUID)
    
    # Change the Frame of Reference UIDs:
    roi2.FrameOfReferenceUID = copy.deepcopy(dicom2.FrameOfReferenceUID)
    roi2.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID = copy.deepcopy(dicom2.FrameOfReferenceUID)
    roi2.StructureSetROISequence[0].ReferencedFrameOfReferenceUID = copy.deepcopy(dicom2.FrameOfReferenceUID)
    
    # Change the Structure Set label.  Assume that the characters before the 
    # first "_" is always the Patient Name:
    roiValue = roi2.StructureSetLabel
    # Find the first "_":
    ind = roiValue.index('_')
    # Modify the Structure Set Label using the Patient Name from the DICOM 
    # followed by the characters including and following the first "_" in 
    # roiValue:
    roi2.StructureSetLabel = str(copy.deepcopy(dicom2.PatientName)) \
    + roiValue[ind::]  
    
    # Change the ROI Name in the Structure Set ROI Sequence.  Assume that the 
    # characters before the first "_" is always the Patient Name:
    roiValue = roi2.StructureSetROISequence[0].ROIName
    # Find the first "_":
    ind = roiValue.index('_')
    # Modify the ROI Name using the Patient Name from the DICOM followed by
    # the characters including and following the first "_" in roiValue:
    roi2.StructureSetROISequence[0].ROIName = str(copy.deepcopy(dicom2.PatientName)) \
    + roiValue[ind::]
    
    # Loop through each Contour Image Sequence and each corresponding DICOM:
    for i in range(len(imageInds)):
        # Load the imageInds[i]^th DICOM file in fpaths2:
        fpath2 = fpaths2[imageInds[i]]
        
        # Read in the DICOM:
        dicom2 = pydicom.read_file(fpath2)
        
        # Replace the DICOM tags in roi2:
        roi2.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence[i].ReferencedSOPInstanceUID = copy.deepcopy(dicom2.SOPInstanceUID)
        
    # Loop through each Contour Sequence and each corresponding DICOM:
    for i in range(len(contourInds)):
        # Load the contourInds[i]^th DICOM file in fpaths2:
        fpath2 = fpaths2[contourInds[i]]
        
        # Read in the DICOM:
        dicom2 = pydicom.read_file(fpath2)
        
        # Replace the DICOM tags in roi2:
        roi2.ROIContourSequence[0].ContourSequence[i].ContourImageSequence[0].ReferencedSOPInstanceUID = copy.deepcopy(dicom2.SOPInstanceUID)
        
    
    # To help identify the ROI Collection from others update dates and times:
    #roi2.StudyDate = time.strftime("%Y%m%d", time.gmtime())
    roi2.StructureSetDate = time.strftime("%Y%m%d", time.gmtime())
    roi2.StructureSetTime = time.strftime("%H%M%S", time.gmtime())
    
    # Create a new filename for the DICOM-RTSTRUCT file:
    #newRoiFname = 'AIM_' + roi2.StructureSetDate + '_' \
    #+ roi2.StructureSetTime + '.dcm'
    
    # Create a new filename for the DICOM-RTSTRUCT file:
    # Start by getting the first part of the filename of the original RTSTRUCT 
    # file:
    origRoiDir, origRoiFname = os.path.split(origRoiFpath)
    
    # Since the format is usually something like "AIM_YYYYMMDD_HHMMSS.dcm" or
    # "RTSTRUCT_YYYYMMDD_HHMMSS.dcm".  Assuming the format will be one of these
    # or similar, find the first "_" character and use the characters that 
    # precede it in the new file name:
    ind = origRoiFname.index('_')
    
    # The new filename:
    newRoiFname = origRoiFname[0:ind+1] + roi2.StructureSetDate + '_' \
    + roi2.StructureSetTime + '.dcm'
    
    # The full filepath of the new DICOM-RTSTRUCT file:
    newRoiFpath = os.path.join(newRoiDir, newRoiFname)
    
    # Save the new file:
    roi2.save_as(newRoiFpath)
    
    return newRoiFpath