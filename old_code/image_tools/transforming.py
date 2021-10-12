# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 13:07:07 2021

@author: ctorti
"""


""" Transforming functions for SimpleITK images. 
    
These are effectively resampling functions but separated from those contained
in image_tools.resampling.py to keep package sizes more manageable.

MAYBE NOT NEEDED?
"""

import SimpleITK as sitk


def transform_im(MovIm, FixIm, SitkTransform, Interp='Linear', 
                 LogToConsole=False):
    """
    Transform (i.e. resample) a 3D SimpleITK image using SimpleITK.
    
    Parameters
    ----------
    MovIm : SimpleITK Image 
        The 3D image to be transformed.
    FixIm : SimpleITK Image
        The 3D image whose grid MovIm is to be transformed to.
    SitkTransform : SimpleITK Transform
        The transform (e.g. image registration transform) that maps Im to TrgIm.
    Interp : str, optional ('Linear' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Linear'
        - 'NearestNeighbor'
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    TxIm : SimpleITK Image
        The 3D transformed image.
    """
            
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of TransformImageSitk():')
        print('\n\n', '-'*120)
        
    
    if Interp == 'NearestNeighbor':
        #sitkInterp = sitk.sitkLinear # 13/05/21
        sitkInterp = sitk.sitkNearestNeighbor # 13/05/21
    else:
        #sitkInterp = sitk.sitkNearestNeighbor # 13/05/21
        sitkInterp = sitk.sitkLinear # 13/05/21
    
    
    TxIm = sitk.Resample(MovIm, FixIm, SitkTransform, sitkInterp)
    
    if LogToConsole:
        print('-'*120)
    
    return TxIm