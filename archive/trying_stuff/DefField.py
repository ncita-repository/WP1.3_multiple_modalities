# -*- coding: utf-8 -*-
"""
Created on Fri May 29 12:23:10 2020

@author: ctorti
"""


# Import packages:
import SimpleITK as sitk
#print(sitk.Version())
import numpy as np
import time
import os
#import copy
#import importlib
import json
#import pydicom
#import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('matplotlib', 'notebook')
#from ipywidgets import interact, fixed
#from IPython.display import clear_output

#import gui
#import registration_gui as rgui

#from GetDicomFpaths import GetDicomFpaths
#from GetDicoms import GetDicoms
#from GetContourPoints3D import GetContourPoints3D
#from GetAllContourPoints3D import GetAllContourPoints3D
#from CreateInputFileForElastix import CreateInputFileForElastix
#from ParseTransformixOutput import ParseTransformixOutput
#import RegUtilityFuncs as ruf



# Get the current working directory:
CurrentWorkingDir = os.getcwd()

# Import DataDict:
ImportDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12'
JsonFpath = os.path.join(ImportDir, 'DataDict.json')

file = open(JsonFpath)
DataDict = json.load(file)
file.close()



# Define FixedDir and MovingDir based on the directories contained in 
# DataDict:
FixedDir = DataDict['Source']['DicomDir']
MovingDir = DataDict['Target']['DicomDir']

SwapFixedMoving = False
#SwapFixedMoving = True

#if SwapFixedMoving:
#    print('FixedDir and MovingDir swapped.')
#else:
#    print('FixedDir and MovingDir not swapped.')

# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
#FixedReader = sitk.ImageSeriesReader(sitk.sitkFloat32)
if SwapFixedMoving:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(MovingDir)
else:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()
FixedIm = sitk.Cast(FixedIm, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
#MovingReader = sitk.ImageSeriesReader(sitk.sitkFloat32)
if SwapFixedMoving:
    MovingNames = MovingReader.GetGDCMSeriesFileNames(FixedDir)
else:
    MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()
#MovingIm = MovingReader.Execute(sitk.sitkFloat32)
MovingIm = sitk.Cast(MovingIm, sitk.sitkFloat32)


# Register MovingIm to FixedIm using SimpleElastix:

# Start timing:
times = []
times.append(time.time())

# Initiate ElastixImageFilter:
ElastixImFilt = sitk.ElastixImageFilter()
#ElastixImFilt.LogToConsoleOn() # <-- no output in Jupyter
ElastixImFilt.LogToConsoleOff()
ElastixImFilt.LogToFileOn() # <-- outputs ONLY ONCE per kernel in Jupyter

# Define the fixed and moving images:
ElastixImFilt.SetFixedImage(FixedIm)
ElastixImFilt.SetMovingImage(MovingIm)

# Get the default parameter map template:
#ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ParamMap = sitk.GetDefaultParameterMap('affine')

# Set registration parameters:
ParamMap['AutomaticTransformInitialization'] = ['true']
ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
ParamMap['WriteIterationInfo'] = ['true']
ParamMap['MaximumNumberOfIterations'] = ['512']
ParamMap['UseDirectionCosines'] = ['true']

""" Note:
    Steve mentioned increasing the resolution of the BSpline but did not 
    elaborate.
    
    Was he referring to:
        
        "FixedImageBSplineInterpolationOrder" (default value is 1)
        "FixedKernelBSplineOrder" (default value is 0)
        "FinalBSplineInterpolationOrder" (default value is 3) <-- this is in transform parameter file
        
    or a different parameter within the parameter map?
    
    Seems to me that FinalBSplineInterpolationOrder is the most likely parameter
    of relevance, so I'll begin by trying to increase the order of that.
"""
# The default value for FinalBSplineInterpolationOrder is 3.  Increase it:
""" Note that Spline order must be between 0 and 5. """
#ParamMap['FinalBSplineInterpolationOrder'] = ['5'] 
""" The above doesn't seem to have helped. """

# Print the parameters:
#for keys,values in ParamMap.items():
#    print(keys, '=', values)

# Set the parameter map:
ElastixImFilt.SetParameterMap(ParamMap)

# Register the 3D images:
ElastixImFilt.Execute()

# Get the registered image:
RegIm = ElastixImFilt.GetResultImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')



# Compute the deformation field:

""" 
For backwards compatibility, only one of ComputeDeformationFieldOn() 
or SetFixedPointSetFileName() can be active at any one time. 
"""



# Get the transform parameter map:
#TransformParameterMap = ElastixImFilt.GetParameterMap()

# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
#TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter 
TransformixImFilt.LogToConsoleOff()
TransformixImFilt.LogToFileOn() 
TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output

# Set up for computing deformation field:
TransformixImFilt.ComputeDeformationFieldOn() # <-- no output
# Set up other computations to see if they produce output:
TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'

# Set the transform parameter map:
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())

# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)

# Transform the contour points:
TransformixImFilt.Execute()
#TransformixIm = TransformixImFilt.Execute()

# Get the deformation field:
DefField = TransformixImFilt.GetDeformationField()


# Get some info on DefField:
DefFieldDtype = DefField.GetPixelIDTypeAsString()
DefFieldSize = DefField.GetSize()
DefFieldSpacing = DefField.GetSpacing()
DefFieldOrigin = DefField.GetOrigin()
""" Note: 
    Executing this kills the kernel in Anaconda Jupyter. 
"""
DefFieldNda = sitk.GetArrayFromImage(DefField)

print('\nDefField Dtype    =', DefFieldDtype)
print('DefField Size     =', DefFieldSize)
print('DefField Spacing  =', DefFieldSpacing)
print('DefField Origin   =', DefFieldOrigin)
print(f'DefField Min, Max = {np.min(DefFieldNda)}, {np.max(DefFieldNda)}')


# Export the deformation field image to .nii format:
""" Note: 
    Executing this kills the kernel in Anaconda Jupyter. 
"""
DefFieldNiiFpath = os.path.join(CurrentWorkingDir, 'defField.nii')
print('\nExporting deformation field as NII to', DefFieldNiiFpath)
sitk.WriteImage(DefField, DefFieldNiiFpath)

# Export the numpy data array as a .csv:

# Take the product of the dimensions of DefField:
dim2 = int(DefFieldSize[0])*int(DefFieldSize[1])*int(DefFieldSize[2])

# Flatten the 3D array to a 2D array:
DefFieldNda2D = DefFieldNda.reshape(1, dim2)

DefFieldCsvFpath = os.path.join(CurrentWorkingDir, 'defField.csv')
print('\nExporting deformation field as CSV to', DefFieldCsvFpath)
np.savetxt(DefFieldCsvFpath, DefFieldNda2D, delimiter=',')


