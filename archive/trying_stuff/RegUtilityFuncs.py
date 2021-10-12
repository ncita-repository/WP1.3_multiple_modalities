# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:18:10 2020

@author: ctorti
"""



""" Define utility functions for plotting of registration results """


""" Source:
https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb
"""


# Import packages:
import matplotlib.pyplot as plt
import SimpleITK as sitk
#import time
import numpy as np
#import os
#import re
from plotly.subplots import make_subplots 
import plotly.express as px
import plotly.graph_objects as go
    



def ShowDefFieldInfo(DeformationField):
    # Get some info on DeformationField:
    DefFieldDtype = DeformationField.GetPixelIDTypeAsString()
    DefFieldSize = DeformationField.GetSize()
    DefFieldSpacing = DeformationField.GetSpacing()
    DefFieldOrigin = DeformationField.GetOrigin()
    """ Note: 
        Executing this kills the kernel in Anaconda Jupyter. 
    """
    DefFieldNda = sitk.GetArrayFromImage(DeformationField)
    
    print('\nDefField Dtype    =', DefFieldDtype)
    print('DefField Size     =', DefFieldSize)
    print('DefField Spacing  =', DefFieldSpacing)
    print('DefField Origin   =', DefFieldOrigin)
    print(f'DefField Min, Max = {np.min(DefFieldNda)}, {np.max(DefFieldNda)}')

    return


def ShowImagesInfo(FixIm, MovIm):
    #FixDtype = str(FixIm.dtype)
    FixDtype = FixIm.GetPixelIDTypeAsString()
    #MovDtype = str(MovIm.dtype)
    MovDtype = MovIm.GetPixelIDTypeAsString()
    
    print('\nFixIm Dtype =', FixDtype)
    print('MovIm Dtype =', MovDtype)
    
    FixSize = FixIm.GetSize()
    MovSize = MovIm.GetSize()
    
    print('\nFixIm Size =', FixSize)
    print('MovIm Size =', MovSize)
    
    FixSpacing = FixIm.GetSpacing()
    MovSpacing = MovIm.GetSpacing()
    
    print('\nFixIm Spacing =', FixSpacing)
    print('MovIm Spacing =', MovSpacing)
    
    FixOrigin = FixIm.GetOrigin()
    MovOrigin = MovIm.GetOrigin()
    
    print('\nFixIm Origin =', FixOrigin)
    print('MovIm Origin =', MovOrigin)
    
    FixNda = sitk.GetArrayFromImage(FixIm)
    MovNda = sitk.GetArrayFromImage(MovIm)
    
    print(f'\nFixIm Min, Max value = {np.min(FixNda)}, {np.max(FixNda)}')
    print(f'MovIm Min, Max value = {np.min(MovNda)}, {np.max(MovNda)}')
    
    print(f'\nFixIm Min, Max x-pos = {FixOrigin[0]}, {FixOrigin[0] + FixSize[0]*FixSpacing[0]}')
    print(f'FixIm Min, Max y-pos = {FixOrigin[1]}, {FixOrigin[1] + FixSize[1]*FixSpacing[1]}')
    if len(FixSize) == 3:
        print(f'FixIm Min, Max z-pos = {FixOrigin[2]}, {FixOrigin[2] + FixSize[2]*FixSpacing[2]}')
    
    print(f'\nMovIm Min, Max x-pos = {MovOrigin[0]}, {MovOrigin[0] + MovSize[0]*MovSpacing[0]}')
    print(f'MovIm Min, Max y-pos = {MovOrigin[1]}, {MovOrigin[1] + MovSize[1]*MovSpacing[1]}')
    if len(MovSize) == 3:
        print(f'MovIm Min, Max z-pos = {MovOrigin[2]}, {MovOrigin[2] + MovSize[2]*MovSpacing[2]}')
    
    return
        
