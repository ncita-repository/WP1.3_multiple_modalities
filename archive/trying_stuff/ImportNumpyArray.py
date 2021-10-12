# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:54:46 2020

@author: ctorti
"""




def ImportNumpyArray(FileName):
    """
    Import a Numpy array from a .npy file in the current working directory.
    
    Inputs:
        FileName - The file name of the .npy file (which must be in the current  
                   working directory).  If FileName does not include the .npy
                   extension it will be added automatically.       
        
    Returns:
        Npa      - The Numpy array
    """
    
    import numpy as np
    
    if not '.npy' in FileName:
        FileName += '.npy'
    
    Npa = np.load(FileName)
    
    return Npa