# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:39:09 2020

@author: ctorti
"""




def ExportSitkImageToNumpyArray(SitkIm, FileName):
    """
    Export a SimpleITK image to a Numpy array to the current working directory.
    
    Inputs:
        SitkIm   - A SimpleITK image
        
        FileName - The file name to assign to the exported file.  If FileName
                   does not include the .npy extension it will be added
                   automatically.
        
    Returns:
        Nothing returned; 
        The file will be saved in the current working directory
    """
    
    import numpy as np
    import SimpleITK as sitk
    
    # Convert Sitk image to numpy array
    Npa = sitk.GetArrayFromImage(SitkIm)
    
    np.save(FileName, Npa)
    
    return