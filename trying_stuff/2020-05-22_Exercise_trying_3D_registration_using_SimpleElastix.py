#!/usr/bin/env python
# coding: utf-8

# # 3D Registration
# 
# Work continued from previous notebook dated May 7th but removed older stuff.  Retained only the bits that worked properly.

# ### Import packages and functions

# In[1]:


import SimpleITK as sitk
print(sitk.Version())
import numpy as np
import time
import os
import copy
import importlib
import json
import pydicom
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')
from ipywidgets import interact, fixed
from IPython.display import clear_output

import gui
import registration_gui as rgui

from GetDicomFpaths import GetDicomFpaths
from GetDicoms import GetDicoms
#from GetContourPoints3D import GetContourPoints3D
from GetAllContourPoints3D import GetAllContourPoints3D
from CreateInputFileForElastix import CreateInputFileForElastix
from ParseTransformixOutput import ParseTransformixOutput
import RegUtilityFuncs as ruf


# ### Download XNAT data (DICOMs, ROI Collection):
# 
# --> Marked down since already downloaded.

# import importlib
# 
# import GetXnatData
# importlib.reload(GetXnatData)
# from GetXnatData import GetXnatData
# 
# import CopyRoi
# importlib.reload(CopyRoi)
# from CopyRoi import CopyRoi
# 
# 
# """ 
# Make selections below by un-commenting the desired option.
# """
# 
# # Define the XNAT address, username and password:
# XnatAddress = 'http://10.1.1.17'
# XnatUsername = 'admin'
# XnatPassword = 'admin'
# 
# # Decide on the subject:
# Subject = 'ACRIN-FMISO-Brain-011'
# #Subject = ''
# 
# # Decide whether to copy to select PET or CT scans:
# #CopyFrom = 'MR'
# #CopyTo = 'PET'
# #CopyTo = 'CT'
# CopyFrom = 'MR_4'# 1960-04-10 69626
# CopyTo = 'MR_12' # 1961-06-11 79433
# 
# # Decide on what feature to use to manually determine the required shift 
# # along the scan direction:
# #ShiftFeature = 'peak of nose'
# #ShiftFeature = 'infinity-like structure in brain'
# #ShiftFeature = 'lens of eye'
# #ShiftFeature = 'butterfly-like structure in brain'
# ShiftFeature = '' # i.e. no shifting to be done
# 
# # Decide whether or not to print some intermediate results for debugging purposes:
# """ Some of the more important/useful intermediate results are printed by default. """
# #Debug = True
# Debug = False
# 
# 
# 
# # Define the REST variables that correspond to the Source ROI Collection,
# # the corresponding Source DICOMs, and the Target DICOMs that the ROI
# # Collection is to be copied to:
# if Subject == 'ACRIN-FMISO-Brain-011' and CopyFrom == 'MR':
#     SourceRoiProjLabel = 'ACRIN'
#     SourceRoiSubjLabel = 'ACRIN-FMISO-Brain-011'
#     SourceRoiExpLabel = 'ACRIN-FMISO-Brain-011_MR_2'
#     SourceRoiLabelDateTime = '20200330'
#     SourceSeriesNo = '5'
#     SourceSeriesDesc = 'T2 3D Slab'
# 
# if Subject == 'ACRIN-FMISO-Brain-011' and CopyTo == 'PET':
#     TargetDicomProjLabel = 'ACRIN'
#     TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
#     TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_CT_PET_1'
#     TargetSeriesNo = '3'
#     TargetSeriesDesc = '3D_Brain miso '
#     
# if Subject == 'ACRIN-FMISO-Brain-011' and CopyTo == 'CT':
#     TargetDicomProjLabel = 'ACRIN'
#     TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
#     TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_CT_19610312'
#     TargetSeriesNo = '2'
#     TargetSeriesDesc = 'HEAD 2.4 H40 SOFT TISSUE'
#     
# if Subject == 'ACRIN-FMISO-Brain-011' and CopyFrom == 'MR_4' and CopyTo == 'MR_12':
#     SourceRoiProjLabel = 'ACRIN'
#     SourceRoiSubjLabel = 'ACRIN-FMISO-Brain-011'
#     SourceRoiExpLabel = 'ACRIN-FMISO-Brain-011_MR_4'
#     SourceRoiLabelDateTime = '20200511'
#     SourceSeriesNo = '8'
#     SourceSeriesDesc = 'T1 SE AXIAL POST FS FC'
#     
#     TargetDicomProjLabel = 'ACRIN'
#     TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
#     TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_MR_12'
#     TargetSeriesNo = '8'
#     TargetSeriesDesc = 'T1 SE AXIAL POST FS FC'
#     
# 
# # Define the root directory where DICOMs will be downloaded from XNAT to:
# """ Note: Files will be downloaded to RootDir\TodaysDate SourceRoiProjLabel"""
# RootDir = r'C:\Temp'
# 
# # Get the data from XNAT:
# DataDict = GetXnatData(XnatAddress=XnatAddress,
#                        XnatUsername=XnatUsername,
#                        XnatPassword=XnatPassword,
#                        RootDir=RootDir,
#                        CopyFrom=CopyFrom,
#                        CopyTo=CopyTo,
#                        SourceRoiProjLabel=SourceRoiProjLabel,
#                        SourceRoiSubjLabel=SourceRoiSubjLabel,
#                        SourceRoiExpLabel=SourceRoiExpLabel,
#                        SourceRoiLabelDateTime=SourceRoiLabelDateTime,
#                        SourceSeriesNo=SourceSeriesNo,
#                        SourceSeriesDesc=SourceSeriesDesc,
#                        TargetDicomProjLabel=TargetDicomProjLabel,
#                        TargetDicomSubjLabel=TargetDicomSubjLabel,
#                        TargetDicomExpLabel=TargetDicomExpLabel,
#                        TargetSeriesNo=TargetSeriesNo,
#                        TargetSeriesDesc=TargetSeriesDesc,
#                        Debug=Debug)

# DataDict

# ### Export DataDict:
# 
# --> Marked down since already exported.

# import json
# 
# # Create filepath for json file:
# JsonFpath = os.path.join(DataDict['RootDir'], 'DataDict.json')
# 
# # Export DataDict:
# data = json.dumps(DataDict)
# file = open(JsonFpath, 'w')
# file.write(data)
# file.close()

# ### Import DataDict:

# In[2]:


ImportDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12'
JsonFpath = os.path.join(ImportDir, 'DataDict.json')

file = open(JsonFpath)
DataDict = json.load(file)
file.close()

DataDict


# ### Define FixedDir and MovingDir based on the directories contained in DataDict:

# In[3]:


FixedDir = DataDict['Source']['DicomDir']
MovingDir = DataDict['Target']['DicomDir']

SwapFixedMoving = False
#SwapFixedMoving = True

if SwapFixedMoving:
    print('FixedDir and MovingDir swapped.')
else:
    print('FixedDir and MovingDir not swapped.')

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

# Get some info on FixedIm and MovingIm:

#FixedDtype = str(FixedIm.dtype)
FixedDtype = FixedIm.GetPixelIDTypeAsString()
#MovingDtype = str(MovingIm.dtype)
MovingDtype = MovingIm.GetPixelIDTypeAsString()
print('\nFixedIm Dtype  =', FixedDtype)
print('MovingIm Dtype =', MovingDtype)

FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
print('\nFixedIm Spacing  =', FixedSpacing)
print('MovingIm Spacing =', MovingSpacing)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
print('\nFixedIm Origin  =', FixedOrigin)
print('MovingIm Origin =', MovingOrigin)

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
print(f'\nFixedIm Min, Max  = {np.min(FixedNda)}, {np.max(FixedNda)}')
print(f'MovingIm Min, Max = {np.min(MovingNda)}, {np.max(MovingNda)}')


# ### Display images using interactive plot:

# gui.MultiImageDisplay(image_list = [FixedIm, MovingIm],                   
#                       title_list = ['Fixed', 'Moving'], figure_size=(8,4));

# interact(display_images_v2,
#          fixed_image_z=(0,FixedIm.GetSize()[2]-1),
#          moving_image_z=(0,MovingIm.GetSize()[2]-1),
#          fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
#          moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
#          fixed_title = 'Fixed image',
#          moving_title = 'Moving image');

# ### Get the contour points:

# In[4]:


FixedContourPts = GetAllContourPoints3D(DicomDir=DataDict['Source']['DicomDir'],
                                        RoiFpath=DataDict['Source']['RoiFpath'])

FixedContourPts


# ### Create inputpoints.txt for Elastix, containing the contour points to be deformed:

# In[5]:


# Create inputpoints.txt for Elastix, containing the contour points to be
# translated:
CreateInputFileForElastix(DicomDir=DataDict['Source']['DicomDir'],
                          RoiFpath=DataDict['Source']['RoiFpath'])


# In[6]:


### View the created file inputpoints.pcs


# In[7]:


# Open the text file:
#TextFile = open('inputpoints.pcs', 'r')
TextFile = open('inputpoints.txt', 'r')

for line in TextFile:
    print(line)


# ### Register MovingIm to FixedIm using SimpleElastix:

# In[17]:


# Start timing:
times = []
times.append(time.time())

# Initiate ElastixImageFilter:
ElastixImFilt = sitk.ElastixImageFilter()
ElastixImFilt.LogToConsoleOn() # <-- no output in Jupyter
ElastixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter

# Define the fixed and moving images:
ElastixImFilt.SetFixedImage(FixedIm)
ElastixImFilt.SetMovingImage(MovingIm)

# Get the default parameter map template for affine transformation:
#ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ElastixParamMap = sitk.GetDefaultParameterMap('affine')

# Re-assign some parameters:
ElastixParamMap['AutomaticTransformInitialization'] = ['true']
ElastixParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
ElastixParamMap['WriteIterationInfo'] = ['true']
ElastixParamMap['MaximumNumberOfIterations'] = ['512']
ElastixParamMap['UseDirectionCosines'] = ['true']
""" 29/05: Trying this instead of trying to change it for Transformix """
ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0'] 

# Print the parameters:
for keys,values in ElastixParamMap.items():
    print(keys, '=', values)
    
# Set the parameter map:
ElastixImFilt.SetParameterMap(ElastixParamMap)

# Register the 3D images:
ElastixImFilt.Execute()
# Get the registered image:
RegIm = ElastixImFilt.GetResultImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Get some info on FixedIm, MovingIm and RegIm:

# In[5]:


# Get some info on FixedIm, MovingIm and RegIm:

FixedDtype = FixedIm.GetPixelIDTypeAsString()
MovingDtype = MovingIm.GetPixelIDTypeAsString()
RegDtype = RegIm.GetPixelIDTypeAsString()
print('\nFixedIm Dtype  =', FixedDtype)
print('MovingIm Dtype =', MovingDtype)
print('RegIm Dtype    =', RegDtype)

FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
RegSize = RegIm.GetSize()
print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)
print('RegIm Size    =', RegSize)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
RegSpacing = RegIm.GetSpacing()
print('\nFixedIm Spacing  =', FixedSpacing)
print('MovingIm Spacing =', MovingSpacing)
print('RegIm Spacing    =', RegSpacing)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
RegOrigin = RegIm.GetOrigin()
print('\nFixedIm Origin  =', FixedOrigin)
print('MovingIm Origin =', MovingOrigin)
print('RegIm Origin    =', RegOrigin)

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
RegNda = sitk.GetArrayFromImage(RegIm)
print(f'\nFixedIm Min, Max  = {np.min(FixedNda)}, {np.max(FixedNda)}')
print(f'MovingIm Min, Max = {np.min(MovingNda)}, {np.max(MovingNda)}')
print(f'RegIm Min, Max    = {np.min(RegNda)}, {np.max(RegNda)}')


# In[10]:


print(FixedIm)


# In[11]:


print(MovingIm)


# In[12]:


print(RegIm)


# ### Plot the results:

# # Interactive feature doesn't seem to work:
# gui.MultiImageDisplay(image_list = [MovingIm, RegIm],                   
#                       title_list = ['Moving', 'Registered'], figure_size=(8,4));

# # Interactive feature doesn't seem to work:
# gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
#                                     known_transformation=ElastixImFilt);

# In[13]:


interact(ruf.display_images_v3,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         registered_image_z=(0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         registered_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         registered_title = 'Registered image');


# In[14]:


interact(ruf.display_images_and_contours,
         fixed_ind = (0,FixedIm.GetSize()[2]-1),
         moving_ind = (0,MovingIm.GetSize()[2]-1),
         reg_ind = (0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         reg_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_im = fixed(FixedIm),
         moving_im = fixed(MovingIm),
         reg_im = fixed(RegIm),
         fixed_contour_pts = fixed(FixedContourPts),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         reg_title = 'Registered image');


# In[15]:


ImFilter.PrintParameterMap()


# ### Transform the contour points:
# 
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[25]:


# Get the transform parameter map. Start by getting the transfer parameter map
# used for the registration (= ElastixParamMap):
TransformixParamMap = ElastixImFilt.GetParameterMap()

""" 29/05:
    Having read section 4.3 "The transform parameter file" in the Elastix manual 
    I found this interesting:
    
    "An important parameter in the transform parameter files is the FinalBSplineInterpolationOrder.
    Usually it is set to 3, because that produces the best quality result image after registration, see Sec 5.3.4.
    However, if you use transformix to deform a segmentation of the moving image (so, a binary image), you
    need to manually change the FinalBSplineInterpolationOrder to 0. This will make sure that the deformed
    segmentation is still a binary label image. If third order interpolation is used, the deformed segmentation
    image will contain garbage. This is related to the “overshoot-property” of higher-order B-spline interpolation."
    
    So perhaps because I'm looking to transform points, and points are sort of binary, maybe I need to set the
    FinalBSplineInterpolationOrder to 0.
"""
# Set the FinalBSplineInterpolationOrder to 0:
TransformixParamMap['FinalBSplineInterpolationOrder'] = ['0']
""" For some reason the above line produces the error:

TypeError: 'tuple' object does not support item assignment

so I re-assigned the FinalBSplineInterpolationOrder parameter for ElastixParamMap instead (see previous cells).
"""

# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
TransformixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter

# Set the parameter map:
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
#TransformixImFilt.SetTransformParameterMap(TransformixParamMap)
""" For some reason the above line produces the error:
Description: itk::ERROR: Self(00000159EFE661E0): No entry Spacing found in transformParameterMap
"""

# Set the input points:
TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')

# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)

# Set up for computing deformation field:
#TransformixImFilt.ComputeDeformationFieldOn()

# Transform the contour points:
TransformixImFilt.Execute()

# Get the deformation field:
#DefField = TransformixImFilt.GetDeformationField()

# Print the deformation field parameters:
#for keys,values in DeformationField.items():
#    print(keys, '=', values)
# """ AttributeError: 'Image' object has no attribute 'items' """


# ### 29/05:  I don't know why I was unable to assign the FinalBSplineInterpolationOrder for TransformixParamMap, since I was able to re-assign other parameters in ElastixParamMap...
# 
# ### So instead start from scratch, using a parameter map template as was done for ElastixParamMap:

# In[16]:


# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
TransformixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter

# Get the default parameter map template for affine transformation:
TransformixParamMap = sitk.GetDefaultParameterMap('affine')

# Re-assign some parameters:
TransformixParamMap['AutomaticTransformInitialization'] = ['true']
TransformixParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
TransformixParamMap['WriteIterationInfo'] = ['true']
TransformixParamMap['MaximumNumberOfIterations'] = ['512']
TransformixParamMap['UseDirectionCosines'] = ['true']
""" 29/05:
    Having read section 4.3 "The transform parameter file" in the Elastix manual 
    I found this interesting:
    
    "An important parameter in the transform parameter files is the FinalBSplineInterpolationOrder.
    Usually it is set to 3, because that produces the best quality result image after registration, see Sec 5.3.4.
    However, if you use transformix to deform a segmentation of the moving image (so, a binary image), you
    need to manually change the FinalBSplineInterpolationOrder to 0. This will make sure that the deformed
    segmentation is still a binary label image. If third order interpolation is used, the deformed segmentation
    image will contain garbage. This is related to the “overshoot-property” of higher-order B-spline interpolation."
    
    So perhaps because I'm looking to transform points, and points are sort of binary, maybe I need to set the
    FinalBSplineInterpolationOrder to 0.
"""
# Set the FinalBSplineInterpolationOrder to 0:
TransformixParamMap['FinalBSplineInterpolationOrder'] = ['0']

# Print the parameters:
for keys,values in TransformixParamMap.items():
    print(keys, '=', values)
    
# Set the parameter map:
TransformixImFilt.SetParameterMap(TransformixParamMap)
""" AttributeError: 'TransformixImageFilter' object has no attribute 'SetParameterMap' """

# Set the input points:
TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')

# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)

# Set up for computing deformation field:
#TransformixImFilt.ComputeDeformationFieldOn()

# Transform the contour points:
TransformixImFilt.Execute()


# In[32]:


# Perform the registration first (copied here since many cells up):

# Start timing:
times = []
times.append(time.time())

# Initiate ElastixImageFilter:
ElastixImFilt = sitk.ElastixImageFilter()
ElastixImFilt.LogToConsoleOn() # <-- no output in Jupyter
ElastixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter

# Define the fixed and moving images:
ElastixImFilt.SetFixedImage(FixedIm)
ElastixImFilt.SetMovingImage(MovingIm)

# Get the default parameter map template for affine transformation:
#ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ElastixParamMap = sitk.GetDefaultParameterMap('affine')

# Re-assign some parameters:
ElastixParamMap['AutomaticTransformInitialization'] = ['true']
ElastixParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
ElastixParamMap['WriteIterationInfo'] = ['true']
ElastixParamMap['MaximumNumberOfIterations'] = ['512']
ElastixParamMap['UseDirectionCosines'] = ['true']
""" 29/05: Trying this instead of trying to change it for Transformix """
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0'] 
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['1'] 
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['2'] 
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['3'] 
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['4'] 
#ElastixParamMap['FinalBSplineInterpolationOrder'] = ['5'] 

# Print the parameters:
for keys,values in ElastixParamMap.items():
    print(keys, '=', values)
    
# Set the parameter map:
ElastixImFilt.SetParameterMap(ElastixParamMap)

# Register the 3D images:
ElastixImFilt.Execute()
# Get the registered image:
RegIm = ElastixImFilt.GetResultImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# Now use Transformix to transform the points:

# Get the transform parameter map. Start by getting the transfer parameter map
# used for the registration (= ElastixParamMap):
TransformixParamMap = ElastixImFilt.GetParameterMap()

""" 29/05:
    Having read section 4.3 "The transform parameter file" in the Elastix manual 
    I found this interesting:
    
    "An important parameter in the transform parameter files is the FinalBSplineInterpolationOrder.
    Usually it is set to 3, because that produces the best quality result image after registration, see Sec 5.3.4.
    However, if you use transformix to deform a segmentation of the moving image (so, a binary image), you
    need to manually change the FinalBSplineInterpolationOrder to 0. This will make sure that the deformed
    segmentation is still a binary label image. If third order interpolation is used, the deformed segmentation
    image will contain garbage. This is related to the “overshoot-property” of higher-order B-spline interpolation."
    
    So perhaps because I'm looking to transform points, and points are sort of binary, maybe I need to set the
    FinalBSplineInterpolationOrder to 0.
"""
# Set the FinalBSplineInterpolationOrder to 0:
#TransformixParamMap['FinalBSplineInterpolationOrder'] = ['0']
""" For some reason the above line produces the error:

TypeError: 'tuple' object does not support item assignment

so I re-assigned the FinalBSplineInterpolationOrder parameter for ElastixParamMap instead (see previous cells).
"""

# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
TransformixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter

# Set the parameter map:
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
#TransformixImFilt.SetTransformParameterMap(TransformixParamMap)
""" For some reason the above line produces the error:
Description: itk::ERROR: Self(00000159EFE661E0): No entry Spacing found in transformParameterMap
"""

# Set the input points:
TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')

# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)

# Set up for computing deformation field:
#TransformixImFilt.ComputeDeformationFieldOn()

# Transform the contour points:
TransformixImFilt.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to transform the points.')

# Get the deformation field:
#DefField = TransformixImFilt.GetDeformationField()

# Print the deformation field parameters:
#for keys,values in DeformationField.items():
#    print(keys, '=', values)
# """ AttributeError: 'Image' object has no attribute 'items' """


# ### 29/05:  So after succeeding in re-assigning the FinalBSplineInterpolationOrder to 0, I can see that the outputpoints.txt indices still exceed the image dimentions.
# 
# ### So I'll try all allowed values for FinalBSplineInterpolationOrder:  0 to 5 and see if any produce sensible indices in outputpoints.txt. They all seemed to have the same issue.

# ### Compute the deformation field:
# 
# """ For backwards compatibility, only one of ComputeDeformationFieldOn() or SetFixedPointSetFileName() can be active at any one time. """

# In[16]:


# Get the current working directory:
CurrentWorkingDir = os.getcwd()

# Get the transform parameter map:
#TransformParameterMap = ElastixImFilt.GetParameterMap()

# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.LogToConsoleOn() # <-- no output
TransformixImFilt.LogToFileOn() 
TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output
# Set up for computing deformation field:
TransformixImFilt.ComputeDeformationFieldOn() # <-- no output
# Set up other computations to see if they produce output:
TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'
# Set the transform parameter map:
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())

#TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)

# Transform the contour points:
TransformixImFilt.Execute()
#TransformixIm = TransformixImFilt.Execute()

# Get the deformation field:
DeformationField = TransformixImFilt.GetDeformationField()

DeformationField


# # Rather than using the same Parameter map that was used for ElastixImFilt,
# # create a new one so that this cell can be the first execution of .Execute()
# # in the hope that elastix.log will be revealing:
# 
# # Initiate ElastixImageFilter:
# ElastixImFilt = sitk.ElastixImageFilter()
# ElastixImFilt.LogToConsoleOn() # <-- no output
# ElastixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel
# # Define the fixed and moving images:
# ElastixImFilt.SetFixedImage(FixedIm)
# ElastixImFilt.SetMovingImage(MovingIm)
# # Get the default parameter map template:
# #ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
# ParamMap = sitk.GetDefaultParameterMap('affine')
# # Set registration parameters:
# ParamMap['AutomaticTransformInitialization'] = ['true']
# ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
# ParamMap['WriteIterationInfo'] = ['true']
# ParamMap['MaximumNumberOfIterations'] = ['512']
# ParamMap['UseDirectionCosines'] = ['true']
# # Print the parameters:
# for keys,values in ParamMap.items():
#     print(keys, '=', values)
# # Set the parameter map:
# ElastixImFilt.SetParameterMap(ParamMap)
# 
# # Initiate TransformixImageFilter:
# TransformixImFilt = sitk.TransformixImageFilter()
# TransformixImFilt.LogToConsoleOn() # <-- no output
# TransformixImFilt.LogToFileOn() 
# TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output
# # Set up for computing deformation field:
# TransformixImFilt.ComputeDeformationFieldOn() # <-- no output
# # Set up other computations to see if they produce output:
# TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
# TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'
# # Set the transform parameter map:
# TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
# 
# #TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
# # Need to explicitely tell elastix that the image is 3D since by default it will
# # only transform to 2D:
# TransformixImFilt.SetMovingImage(MovingIm)
# 
# # Transform the contour points:
# TransformixImFilt.Execute()
# #TransformixIm = TransformixImFilt.Execute()
# 
# # Get the deformation field:
# DeformationField = TransformixImFilt.GetDeformationField()
# 
# DeformationField

# # 27/05: More requests for help posted:
# 
# https://github.com/SuperElastix/SimpleElastix/issues/181
# 
# https://github.com/SuperElastix/SimpleElastix/issues/182
# 
# https://github.com/SuperElastix/SimpleElastix/issues/259

# # 27/05: This looks potentially useful:
# 
# https://github.com/SuperElastix/SimpleElastix/issues/259

# In[ ]:


# Executing this kills the kernel:
# Convert the deformation field to a numpy array:
#DeformationFieldNda = sitk.GetArrayFromImage(DeformationField) 

# Executing this also kills the kernel:
# So write the image instead:
#sitk.WriteImage(DeformationField, 'defField.nii')


# In[17]:


get_ipython().run_line_magic('pinfo', 'DeformationField')


# ### View the transformed points:

# In[13]:


# Open the text file:
TextFile = open('outputpoints.txt', 'r')

for line in TextFile:
    print(line)


# ### Parse outputpoints.txt:

# In[17]:


# Run ParseTransformixOutput() to parse the output from Transformix:
PtNos, InInds, InPts, FixOutInds, OutPts, Defs, MovOutInds = ParseTransformixOutput()

"""
Note: The relationship between the above variable names and the items in outputpoints.txt are as follows:

PtNos      = Point
InInds     = InputIndex
InPts      = InputPoint
FixOutInds = OutputIndexFixed
OutPts     = OutputPoint
Defs       = Deformation
MovOutInds = OutputIndexMoving
"""

#[print(f'\nPoint {PtNos[i]}:\nInputIndex = {InInds[i]}\nFixedOutputIndex = {FixOutInds[i]}') for i in range(len(PtNos))]

for i in range(len(PtNos)):
    print(f'\nPoint {PtNos[i]}:')
    print(f'InputIndex        = {InInds[i]}')
    print(f'InputPoint        = {InPts[i]}')
    print(f'FixedOutputIndex  = {FixOutInds[i]}')
    print(f'OutputPoint       = {OutPts[i]}')
    print(f'Deformation       = {Defs[i]}')
    print(f'MovingOutputIndex = {MovOutInds[i]}')


# In[26]:


FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
RegSize = RegIm.GetSize()

print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)
print('RegIm Size    =', RegSize)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
RegSpacing = RegIm.GetSpacing()
print('\nFixedIm Spacing  =', FixedSpacing)
print('MovingIm Spacing =', MovingSpacing)
print('RegIm Spacing    =', RegSpacing)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
RegOrigin = RegIm.GetOrigin()
print('\nFixedIm Origin  =', FixedOrigin)
print('MovingIm Origin =', MovingOrigin)
print('RegIm Origin    =', RegOrigin)

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
RegNda = sitk.GetArrayFromImage(RegIm)
print(f'\nFixedIm (min, max) x-pos  = {FixedOrigin[0]}, {FixedOrigin[0] + FixedSize[0]*FixedSpacing[0]}')
print(f'FixedIm (min, max) y-pos  = {FixedOrigin[1]}, {FixedOrigin[1] + FixedSize[1]*FixedSpacing[1]}')
print(f'FixedIm (min, max) z-pos  = {FixedOrigin[2]}, {FixedOrigin[2] + FixedSize[2]*FixedSpacing[2]}')

print(f'\nMovingIm (min, max) x-pos  = {MovingOrigin[0]}, {MovingOrigin[0] + MovingSize[2]*MovingSpacing[0]}')
print(f'MovingIm (min, max) y-pos  = {MovingOrigin[1]}, {MovingOrigin[1] + MovingSize[1]*MovingSpacing[1]}')
print(f'MovingIm (min, max) z-pos  = {MovingOrigin[2]}, {MovingOrigin[2] + MovingSize[2]*MovingSpacing[2]}')

print(f'\nRegIm (min, max) x-pos  = {RegOrigin[0]}, {RegOrigin[0] + RegSize[0]*RegSpacing[0]}')
print(f'RegIm (min, max) y-pos  = {RegOrigin[1]}, {RegOrigin[1] + RegSize[1]*RegSpacing[1]}')
print(f'RegIm (min, max) z-pos  = {RegOrigin[2]}, {RegOrigin[2] + RegSize[2]*RegSpacing[2]}')


# # Interactive feature doesn't seem to work:
# gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
#                                     known_transformation=TransformixImFilt);

# ### If I'm interpretting Kasper's suggestion correcly, I can use the z index from MovingOutputIndex to select the appropriate slice to plot the (x,y) points in OutputPoint.
# 
# 

# In[113]:


# The form of FixedContourPts:
FixedContourPts


# In[19]:


# Create MovingContourPts using the z index from MovingOutputIndex and the OutputPoint:

# Initialise MovingContourPts by creating an array of empty arrays with length equal to
# the number of slices in MovingIm:
MovingContourPts = []
[MovingContourPts.append([]) for i in range(MovingIm.GetSize()[2])]

for i in range(len(PtNos)):
    # Get the point:
    point = OutPts[i]
    
    #print(f'\npoint = {point}')
    
    # Get the slice number for this point:
    sliceno = MovOutInds[i][2]
    
    #print(f'sliceno = {sliceno}')
    
    # Add point to MovingContourPts[sliceno]:
    MovingContourPts[sliceno].append(point)
    
MovingContourPts


# In[21]:


interact(ruf.display_images_and_contours_v2,
         fixed_Zind = (0,FixedIm.GetSize()[2]-1),
         moving_Zind = (0,MovingIm.GetSize()[2]-1),
         reg_Zind = (0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         reg_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_im = fixed(FixedIm),
         moving_im = fixed(MovingIm),
         reg_im = fixed(RegIm),
         fixed_contour_pts = fixed(FixedContourPts),
         moving_contour_pts = fixed(MovingContourPts),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         reg_title = 'Registered image');


# ### 22/05:  Still haven't resolved issue with outputpoints.txt having indices that exceed the image dimensions.
# 
# Hoping that someone will reply to my posts here:
# 
# https://github.com/SuperElastix/SimpleElastix/issues/377
# 
# and here:
# 
# https://discourse.itk.org/t/indices-of-transformixs-outputpoints-txt-exceeds-image-dimensions/3105

# ### 26/05:  Trying advice of Ibrahim here:
# 
# #### https://groups.google.com/forum/#!searchin/elastix-imageregistration/transformix%7Csort:date/elastix-imageregistration/zMs_jA3YZ9w/Gc0ZbJ_RBQAJ
#     
# #### So will need to:
# #### 1) Generate a deformation field from Elastix
# #### 2) Convert contour points to fiducials using Slicer
# #### 3) Load the deformation field in Slicer and apply to fiducials
# #### 4) Convert fiducials back to contour points (?)

# In[18]:


Zind = 14

plt.subplots(1,4,figsize=(17,8))

# draw the fixed image in the first subplot
plt.subplot(1,4,1)
plt.imshow(sitk.GetArrayFromImage(FixedIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('fixed image')
plt.axis('off')

# draw the moving image in the second subplot
plt.subplot(1,4,2)
plt.imshow(sitk.GetArrayFromImage(MovingIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('moving image')
plt.axis('off')

# draw the registered image in the third subplot
plt.subplot(1,4,3)
plt.imshow(sitk.GetArrayFromImage(RegIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('registered image')
plt.axis('off')

# draw the transformix image in the fourth subplot
plt.subplot(1,4,4)
plt.imshow(sitk.GetArrayFromImage(TransformixIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('transformix image')
plt.axis('off')


# In[11]:


help(sitk.TransformixImageFilter())


# In[24]:


print('GetLogFileName =', TransformixImFilt.GetLogFileName())
print('GetLogToConsole =', TransformixImFilt.GetLogToConsole())
print('GetLogToFile =', TransformixImFilt.GetLogToFile())
print('GetOutputDirectory =', TransformixImFilt.GetOutputDirectory())


# In[25]:


print('GetDeformationField =', TransformixImFilt.GetDeformationField())


# In[26]:


print('GetTransformParameter =', TransformixImFilt.GetTransformParameter())


# In[27]:


print('GetTransformParameterMap =', TransformixImFilt.GetTransformParameterMap())


# In[28]:


print('PrintParameterMap =', TransformixImFilt.PrintParameterMap())


# In[30]:


ElastixImFilt.PrintParameterMap()


# In[33]:


sitk.PrintParameterMap(ElastixImFilt.GetParameterMap())


# In[41]:


import sys
sys.stdout.flush()

sitk.PrintParameterMap(ElastixImFilt.GetParameterMap())


# ### 28/05: Wrapped some code into a python file to get around issue with things not working in Jupyter.

# In[ ]:


from GetDeformationField import GetDeformationField

DefField, DefFieldNda = GetDeformationField()


# ### the kernel still dies.  So will have to find another way.
# 
# ### I used the console in Anaconda Spyder and although the kernel didn't die, after executing:
# 
# from GetDeformationField import GetDeformationField
# 
# DefField, DefFieldNda = GetDeformationField()
# 
# Took 8.1 s to register the 3D image stacks.
# 
# ### I executed DefField and DefFieldNda, which resulted in:
# 
# NameError: name 'DefField' is not defined
# 
# NameError: name 'DefFieldNda' is not defined

# ### Then I tried PyCharm.  First I tried importing GetDeformationField in the the console:
# 
# from GetDeformationField import DeformationField
# 
# Traceback (most recent call last):
# 
#   File "C:\Users\ctorti\Anaconda37\lib\site-packages\IPython\core\interactiveshell.py", line 3326, in run_code
#     exec(code_obj, self.user_global_ns, self.user_ns)
#     
#   File "<ipython-input-2-3387b7b22da0>", line 1, in <module>
#     from GetDeformationField import DeformationField
#     
# ImportError: cannot import name 'DeformationField' from 'GetDeformationField' (C:\Code\WP1.3_multiple_modalities\trying_stuff\GetDeformationField.py)

# ### Then in PyCharm's terminal console:

# (base) C:\Code\WP1.3_multiple_modalities\trying_stuff>python
# Python 3.7.4 (default, Aug  9 2019, 18:34:13) [MSC v.1915 64 bit (AMD64)] :: Anaconda, Inc. on win32
# Type "help", "copyright", "credits" or "license" for more information.
# >>> from GetDeformationField import GetDeformationField
# >>> nii, nda = GetDeformationField()
# Installing all components.
# InstallingComponents was successful.
# 
# ELASTIX version: 4.900
# Command line options from ElastixBase:
# -fMask    unspecified, so no fixed mask used
# -mMask    unspecified, so no moving mask used
# -out      ./
# -priority unspecified, so NORMAL process priority
# -threads  unspecified, so all available threads are used
# Command line options from TransformBase:
# -t0       unspecified, so no initial transform used
# 
# Reading images...
# Reading images took 0 ms.
# 
# See C:\Code\WP1.3_multiple_modalities\trying_stuff\2020-05-28_Outputs\2020-05-28_Pycharm_terminal_commands.txt
# for remainder of output.

# ### If I run the same script in Anaconda Powershell and Anaconda Prompt.  It appears to run as before, but does not update defField.nii or create the new file defField.csv, despite no errors.
# 
# ### So it appears that for some reason it has to be run in Pycharm.  Although I've only succeeded in getting the deformation field exported once.  Re-running the function GetDeformationField(), or more recently a script version called DefField.py, has failed to overwrite the original .nii file.
# 
# ### I shut down PyCharm, then repeated script DefField.py.  Again it didn't produce a .nii or .csv file.
# 
# ### Then I ran:
# 
# #### python
# #### from GetDeformationField import GetDeformationField
# #### nii, nda = GetDeformationField()
# 
# ### As before it generated output but no update to .nii or new .csv file.
# 
# ### Then I repeated the above commands this time in PyCharm's Python console.  Same result.
# 
# ### I don't know what else to try...So instead I'll see if I can load the .nii file and examine the contents.

# In[8]:


# Get the current working directory:
CurrentWorkingDir = os.getcwd()

# Read in the deformation field:
DefFieldFpath = os.path.join(CurrentWorkingDir, 'defField.nii')

# Read in the deformation field:
DefField = sitk.ReadImage(DefFieldFpath)

# Convert to numpy data array:
DefFieldNda = sitk.GetArrayFromImage(DefField)


# In[15]:


print(DefField)


# In[11]:


np.shape(DefFieldNda)


# In[16]:


# Get some info on FixedIm, MovingIm, RegIm and DefField:

FixedDtype = FixedIm.GetPixelIDTypeAsString()
MovingDtype = MovingIm.GetPixelIDTypeAsString()
RegDtype = RegIm.GetPixelIDTypeAsString()
DefFieldDtype = DefField.GetPixelIDTypeAsString()
print('\nFixedIm Dtype  =', FixedDtype)
print('MovingIm Dtype =', MovingDtype)
print('RegIm Dtype    =', RegDtype)
print('DefField Dtype =', DefFieldDtype)

FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
RegSize = RegIm.GetSize()
DefFieldSize = DefField.GetSize()
print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)
print('RegIm Size    =', RegSize)
print('DefField Size =', DefFieldSize)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
RegSpacing = RegIm.GetSpacing()
DefFieldSpacing = DefField.GetSpacing()
print('\nFixedIm Spacing  =', FixedSpacing)
print('MovingIm Spacing =', MovingSpacing)
print('RegIm Spacing    =', RegSpacing)
print('DefField Spacing =', DefFieldSpacing)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
RegOrigin = RegIm.GetOrigin()
DefFieldOrigin = DefField.GetOrigin()
print('\nFixedIm Origin  =', FixedOrigin)
print('MovingIm Origin =', MovingOrigin)
print('RegIm Origin    =', RegOrigin)
print('DefField Origin =', DefFieldOrigin)

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
RegNda = sitk.GetArrayFromImage(RegIm)
print(f'\nFixedIm Min, Max  = {np.min(FixedNda)}, {np.max(FixedNda)}')
print(f'MovingIm Min, Max = {np.min(MovingNda)}, {np.max(MovingNda)}')
print(f'RegIm Min, Max    = {np.min(RegNda)}, {np.max(RegNda)}')
print(f'DefField Min, Max = {np.min(DefFieldNda)}, {np.max(DefFieldNda)}')


# In[17]:


Zind = 14

plt.subplots(1,4,figsize=(17,8))

# draw the fixed image in the first subplot
plt.subplot(1,4,1)
plt.imshow(sitk.GetArrayFromImage(FixedIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('fixed image')
plt.axis('off')

# draw the moving image in the second subplot
plt.subplot(1,4,2)
plt.imshow(sitk.GetArrayFromImage(MovingIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('moving image')
plt.axis('off')

# draw the registered image in the third subplot
plt.subplot(1,4,3)
plt.imshow(sitk.GetArrayFromImage(RegIm)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('registered image')
plt.axis('off')

# draw the deformation field in the fourth subplot
plt.subplot(1,4,4)
plt.imshow(sitk.GetArrayFromImage(DefField)[Zind,:,:], cmap=plt.cm.Greys_r);
plt.title('deformation field')
plt.axis('off')


# ### As I had suspected from the output in the console, the deformation field values are all zero!

# ### 29/05: Found a Stackoverflow thread where someone referred to needing to increase the resolution of the BSpline resolution, but lacked details:
# 
# https://stackoverflow.com/questions/60854606/simple-itk-problem-when-saving-deformation-field-the-deformation-field-is-empt/62080612#62080612
# 
# Some egotistic users didn't like my little rant and one even decided to delete my post, so have had to post a new question:
# 
# Thanks to stupid egomaniacs I wasted another half hour preparing a new post:
# 
# https://stackoverflow.com/questions/62083178/simpleelastix-deformation-field-is-null

# ### I increased the value of FinalBSplineInterpolationOrder from 3 (default) to 5 (maximum allowed value).  Try to import the .nii exported and see what it looks like:

# In[8]:


# Get the current working directory:
CurrentWorkingDir = os.getcwd()

# Read in the deformation field:
DefFieldFpath = os.path.join(CurrentWorkingDir, 'defField.nii')

# Read in the deformation field:
DefField = sitk.ReadImage(DefFieldFpath)

# Convert to numpy data array:
DefFieldNda = sitk.GetArrayFromImage(DefField)


# In[11]:


# Get some info on DefField:
DefFieldDtype = DefField.GetPixelIDTypeAsString()
DefFieldSize = DefField.GetSize()
DefFieldSpacing = DefField.GetSpacing()
DefFieldOrigin = DefField.GetOrigin()

print('\nDefField Dtype    =', DefFieldDtype)
print('DefField Size     =', DefFieldSize)
print('DefField Spacing  =', DefFieldSpacing)
print('DefField Origin   =', DefFieldOrigin)
print(f'DefField Min, Max = {np.min(DefFieldNda)}, {np.max(DefFieldNda)}')


# In[ ]:




