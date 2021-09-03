# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:12:04 2021

@author: ctorti
"""

import json
import SimpleITK as sitk


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
    im : SimpleITK object
        SimpleITK 3D image of the DICOM stack.
    
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
    
    return im

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
