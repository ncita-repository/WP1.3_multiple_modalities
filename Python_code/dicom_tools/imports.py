# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:28:28 2021

@author: ctorti
"""


import SimpleITK as sitk
from pydicom import read_file



def get_dcm_fpaths(dicomDir):
    """
    Get a list of full filepaths for all DICOM files in a directory,
    
    Based on SimpleITK's ImageSeriesReader.
    
    Parameters
    ----------
    dicomDir : str
        Path to directory containing DICOMs.
    
    Returns
    -------
    filepaths : list of strs
        List of the full filepaths of all DICOM files in dicomDir.
    """
    
    reader = sitk.ImageSeriesReader()
    
    filepaths = reader.GetGDCMSeriesFileNames(dicomDir)
    
    return filepaths

def import_dcms(dicomDir, sortMethod='slices', p2c=False):
    """
    Import DICOM objects from a directory.
    
    Parameters
    ----------
    dicomDir : str
        Path to directory containing DICOM files.
    sortMethod : str, optional ('slices' by default) 
        Denotes the method used for sorting of files fetched using 
        GetDicomFpaths() with acceptable values:
            -- 'none' (or 'None' or any other strings) = no sorting
            -- 'natural' (or 'Natural') = using natsort to sort the filenames
            -- 'slices' (or 'Slices') = sorted along the slice direction 
    p2c : bool, optional (False by default)
        If True results will be printed to the console.
    
    Returns
    -------
    Dicoms : list of Pydicom objects
        List of DICOM objects.
    """
    
    #from pydicom import read_file
    
    #print('\nsortMethod (input for GetDicoms()) =', sortMethod)
    
    #fpaths = get_dcm_fpaths(dicomDir, sortMethod, p2c)
    fpaths = get_dcm_fpaths(dicomDir)
    
    dicoms = []
    
    for f in range(len(fpaths)):
        fpath = fpaths[f]
        
        dicom = read_file(fpath)
        
        dicoms.append(dicom)
        
    return dicoms

def import_dcm(dicomDir, sortMethod='slices', inStackNum=0, p2c=False):
    """
    Import a single DICOM object from a directory.
    
    Parameters
    ----------
    dicomDir : str
        Path to directory containing DICOM files.
    sortMethod : str, optional ('slices' by default) 
        Denotes the method used for sorting of files fetched using 
        GetDicomFpaths() with acceptable values:
            -- 'none' (or 'None' or any other strings) = no sorting
            -- 'natural' (or 'Natural') = using natsort to sort the filenames
            -- 'slices' (or 'Slices') = sorted along the slice direction
    inStackNum : int, optional (0 by default)
        The in-stack index of the file within the sorted list of DICOM 
        filepaths to import.  The first file is imported by default.
    p2c : bool, optional (False by default)
        If True results will be printed to the console.
    
    Returns
    -------
    dicom : Pydicom Object
    """
    
    #from pydicom import read_file
    
    #fpaths = get_dcm_fpaths(dicomDir, sortMethod, p2c)
    fpaths = get_dcm_fpaths(dicomDir)
    
    if inStackNum < len(fpaths):
        dicom = read_file(fpaths[inStackNum])
        
        return dicom
    
    else:
        msg = 'It is not possible to import a DICOM at inStackNum = '\
              + f'{inStackNum} since there are only {len(fpaths)} DICOMs in '\
              + 'the specified directory.'
        raise Exception(msg)
        
        return None