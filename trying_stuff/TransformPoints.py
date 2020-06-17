# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:01:12 2020

@author: ctorti
"""

def TransformPoints(MovingImage, TransformFilter):
    # Import packages:
    import SimpleITK as sitk
    
    # Initiate TransformixImageFilter:
    TransformixImFilt = sitk.TransformixImageFilter()
    #TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TransformixImFilt.LogToConsoleOff() # <-- no output in Jupyter
    TransformixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    # Get the transform parameter map. Start by getting the transfer parameter map
    # used for the registration (= ElastixParamMap):
    #TransformixParamMap = ImageFilter.GetParameterMap() 
    
    # Set the parameter map:
    TransformixImFilt.SetTransformParameterMap(TransformFilter.GetTransformParameterMap())
    #TransformixImFilt.SetTransformParameterMap(Transform)
    
    # Set the input points (assumed to be in the current working directory):
    TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
    
    # Need to explicitely tell elastix that the image is 3D since by default it 
    # will only transform to 2D:
    TransformixImFilt.SetMovingImage(MovingImage)
    
    # Set up for computing deformation field:
    #TransformixImFilt.ComputeDeformationFieldOn()
    
    # Transform the contour points:
    TransformixImFilt.Execute()
    
    #CWD = os.getcwd()
    
    print('\nThe points in inputpoints.txt have been transformed to outputpoints.txt.')
    
    return