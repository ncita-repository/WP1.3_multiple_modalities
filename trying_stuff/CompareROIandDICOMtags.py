# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

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

import pydicom


def CompareROIandDICOMtags(dicomFpaths, roi, imageInds, contourInds):
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
