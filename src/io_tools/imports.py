# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:12:04 2021

@author: ctorti
"""

"""
from importlib import reload
import dicom_tools.dcm_metadata
reload(dicom_tools.dcm_metadata)
"""

import json
import SimpleITK as sitk
from pydicom import read_file
from dicom_tools.dcm_metadata import get_dcm_fpaths


def import_dcms(dicomDir, sortMethod='slices', p2c=False):
    """
    Import DICOM objects from a directory as a list of Pydicom Objects.
    
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

def import_dicoms_as_im(dicomDir, dtype='float32'):
    """
    Import a DICOM series as a SimpleITK image.
    
    Parameters
    ----------
    dicomDir : str
        Directory containing the DICOMs.
    dtype : string, optional ('float32' by default)
        The pixel ID type as defined at the link in the Notes section. 
        Acceptable values are:
            - 'SameAsRefImage' (default)
            - 'uint8'
            - 'int8'
            - 'uint16'
            - 'int16'
            - 'uint32'
            - 'int32'
            - 'uint64'
            - 'int64'
            - 'float32'
            - 'float64'
    
    Returns
    -------
    im : SimpleITK Object
        SimpleITK 3D image representation of the DICOM stack.
    pixarr : Numpy Data array
        Numpy 3D array representation of the DICOM stack.
    
    Notes
    -----
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    if dtype == 'uint8':
        pixIDtype = sitk.sitkUInt8
    elif dtype == 'int8':
        pixIDtype = sitk.sitkInt8
    elif dtype == 'uint16':
        pixIDtype = sitk.sitkUInt16
    elif dtype == 'int16':
        pixIDtype = sitk.sitkInt16
    elif dtype == 'uint32':
        pixIDtype = sitk.sitkUInt32
    elif dtype == 'int32':
        pixIDtype = sitk.sitkInt32
    elif dtype == 'uint64':
        pixIDtype = sitk.sitkUInt64
    elif dtype == 'int64':
        pixIDtype = sitk.sitkInt64
    elif dtype == 'float32':
        pixIDtype = sitk.sitkFloat32
    elif dtype == 'float64':
        pixIDtype = sitk.sitkFloat64
    else:
        msg = f"dtype '{dtype}' is not a valid input."
        raise Exception(msg)
    
    reader = sitk.ImageSeriesReader()
    fnames = reader.GetGDCMSeriesFileNames(dicomDir)
    reader.SetFileNames(fnames)
    #reader.ReadImageInformation() # doesn't exist for ImageSeriesReader; works
    # for ImageFileReader
    im = reader.Execute()
    #im = sitk.Cast(im, sitk.sitkFloat32)
    im = sitk.Cast(im, pixIDtype)
    
    return im, sitk.GetArrayViewFromImage(im)

def import_nifti(filepath):
    """
    Import a NIFTI as a SimpleITK image.
    
    Parameters
    ---------- 
    Filepath : str
        Filepath of NIFTI file.
        
    Returns
    -------
    im : SimpleITK object
    """
    
    reader = sitk.ImageFileReader()
    reader.SetImageIO('NiftiImageIO')
    reader.SetFileName(filepath)
    im = reader.Execute()
    
    return im

def import_im(filepath, fileFormat):
    """ 
    Export a SimpleITK Image in a desired format. 
    
    Parameters
    ----------
    filepath : str
        The filepath of the file to import.
    fileFormat : str
        The file format can be any format specified here:
            https://simpleitk.readthedocs.io/en/master/IO.html#transformations
        e.g. 'HDF5ImageIO', 'NiftiImageIO', 'NrrdImageIO'. 
    
    Returns
    -------
    im : SimpleITK Image 
        The imported image.
    
    Note
    ----
    The file extension will be defined based on fileFormat. Note that this part
    of this function is incomplete.
    """
    
    reader = sitk.ImageFileReader()
    reader.SetImageIO(fileFormat)
    reader.SetFileName(filepath)
    im = reader.Execute()
    
    return im

def import_dict_from_json(filepath):
    """
    Import a dictionary from a JSON file from the current working directory.
    
    Parameters
    ----------
    filepath : str
        The file path of the JSON file (which must be in the current working 
        directory).  filepath does not need to include the .json extension.
        
    Returns
    -------
    out : dict
    """
    
    if not '.json' in filepath:
        filepath += '.json'
    
    with open(filepath, 'r') as file:
        out = json.load(file)
    
    return out

def import_list_from_txt(filepath):
    """
    Import a list from a text file from the current working directory.
    
    Parameters
    ----------   
    filepath : str
        The file path of the text file (which must be in the current working 
        directory).  filepath does not need to include the .txt extension.
    
    Returns
    -------  
    out : list
    """
    
    if not '.txt' in filepath:
        filepath += '.txt'
    
    out = []
    
    with open(filepath, 'r') as file:
        for line in file:
            out.append(line.replace('\n', '').split(' '))
    
    return out
