# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:13:25 2020

@author: ctorti
"""

def test():

    import SimpleITK as sitk
    
    sitk.PrintParameterMap(sitk.GetDefaultParameterMap('affine'))
    
    return