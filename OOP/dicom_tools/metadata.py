# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:07:41 2021

@author: ctorti
"""

import SimpleITK as sitk


def get_dcm_uids(dicomDir):
    """
    Get the Study, Series, Frame of Reference and SOP Instance UIDs for a DICOM 
    series.
    
    Parameters
    ----------
    dicomDir : str
        Directory containing the DICOMs.
        
    Returns
    -------
    studyUID : str
        Study Instance UID of the DICOM series / XNAT scan.
    seriesUID : str
        Series Instance UID of the DICOM series / XNAT scan.
    FORUID : str
        Frame of reference UID of the DICOM series / XNAT scan.
    SOPUIDs : list of strs
        List of the SOP UIDs of the DICOM series / XNAT scan.
    """
    
    seriesReader = sitk.ImageSeriesReader()
    
    fpaths = seriesReader.GetGDCMSeriesFileNames(dicomDir)
    
    fileReader = sitk.ImageFileReader()
    
    SOPUIDs = [] 
    
    for fpath in fpaths:
        fileReader.SetFileName(fpath)

        fileReader.ReadImageInformation()
    
        SOPUIDs.append(fileReader.GetMetaData('0008|0018'))
        
    studyUID = fileReader.GetMetaData('0020|000d')
    seriesUID = fileReader.GetMetaData('0020|000e')
    FORUID = fileReader.GetMetaData('0020|0052')
        
    return studyUID, seriesUID, FORUID, SOPUIDs

def get_roicol_labels(roicol):
    """
    Get the list of contour/segment labels in an ROI Collection (RTS/SEG).
    
    Parameters
    ----------
    roicol : Pydicom Object
        Object from a RTS/SEG file.
        
    Returns
    -------
    labels : list of strs
        List (for each ROI/segment) of ROIName/SegmentLabel tag values.
    """
    
    if roicol == None:
        return None
    
    if str(type(roicol)) != "<class 'pydicom.dataset.FileDataset'>":
        msg = "The input argument is not of the required Pydicom FileDataset."
        raise Exception(msg)
    
    try:
        mod = roicol.Modality
        
        labels = [] 
        
        if 'RTSTRUCT' in mod:
            sequences = roicol.StructureSetROISequence
            
            for sequence in sequences:
                label = sequence.ROIName
                
                labels.append(label)
                
        elif 'SEG' in mod:
        
            sequences = roicol.SegmentSequence
            
            for sequence in sequences:
                label = sequence.SegmentLabel
                
                labels.append(label) 
                
        else:
            msg = "The ROI collection is not recognised as being RTS or SEG."
            raise Exception(msg)
        
        return labels
            
    except AttributeError:
        # The Pydicom object does not have a Modality attribute.
        return None

def get_roicol_nums(roicol, searchStr):
    """
    Get the ROI/segment number(s) that matches searchString.
    
    Parameters
    ----------
    roicol : Pydicom Object
        RTS or SEG object.
    searchStr : str
        All or part of the ROIName/SegmentLabel of the ROI/segment of interest.
        If searchString is an empty string return the first ROI/segment number. 
                       
    Returns
    -------
    roiNums : list of ints
        A list of indices corresponding to the ROI(s)/segment(s) that matches 
        searchString.
    """
    
    roiLabels = get_roicol_labels(roicol)
    
    roiNums = []
    
    if not roiLabels:
        msg = 'There are no ROIs in the RTS/SEG.'
        raise Exception(msg)
    
    if not searchStr in ["", None]: 
        roiNums = [i for i, roiLabel in enumerate(roiLabels) if searchStr in roiLabel]
        
        #print(f'\nroiNum = {roiNum}')
        
        #if roiNum:
        #    roiNum = roiNum[0] # 13/01/2021
        #    
        #else:
        #    msg = "There is no ROI/segment in the RTS/SEG containing "\
        #          + f"names/labels {roiLabels} \nthat matches 'searchStr' "\
        #          + f"= '{searchStr}'."
        #    raise Exception(msg)
            
        if not roiNums:
            msg = "There is no ROI/segment in the RTS/SEG containing "\
                  + f"names/labels {roiLabels} \nthat matches 'searchStr' "\
                  + f"= '{searchStr}'."
            
            raise Exception(msg)
            
    return roiNums