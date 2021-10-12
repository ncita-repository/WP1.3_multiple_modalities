#!/usr/bin/env python
# coding: utf-8

# # 3D Registration
# 
# Work continued from previous notebook dated May 7th but removed older stuff.  Retained only the bits that worked properly.

# ### Import packages and functions

# In[1]:


import SimpleITK as sitk
#print(sitk.Version())
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

# In[5]:


FixedContourPts = GetAllContourPoints3D(DicomDir=DataDict['Source']['DicomDir'],
                                        RoiFpath=DataDict['Source']['RoiFpath'])

FixedContourPts


# ### Create inputpoints.txt for Elastix, containing the contour points to be deformed:

# In[6]:


# Create inputpoints.txt for Elastix, containing the contour points to be
# translated:
CreateInputFileForElastix(DicomDir=DataDict['Source']['DicomDir'],
                          RoiFpath=DataDict['Source']['RoiFpath'])


# In[7]:


### View the created file inputpoints.pcs


# In[8]:


# Open the text file:
TextFile = open('inputpoints.pcs', 'r')

for line in TextFile:
    print(line)


# ### Register MovingIm to FixedIm using SimpleElastix:

# In[9]:


# Start timing:
times = []
times.append(time.time())

# Initiate ElastixImageFilter:
ElastixImFilt = sitk.ElastixImageFilter()
ElastixImFilt.SetFixedImage(FixedIm)
ElastixImFilt.SetMovingImage(MovingIm)
ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ElastixImFilt.LogToConsoleOn()
ElastixImFilt.LogToFileOn()
# Register the 3D images:
ElastixImFilt.Execute()
# Get the registered image:
RegIm = ElastixImFilt.GetResultImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Get some info on FixedIm, MovingIm and RegIm:

# In[10]:


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

# Get the Transform Parameter Map:
#TransformParameterMap = sitk.GetTransformParameterMap() # AttributeError: module 'SimpleITK' 
# has no attribute 'GetTransformParameterMap'
TransformParameterMap = ElastixImFilt.GetTransformParameterMap()
TransformParameterMap


# In[10]:


print(FixedIm)


# In[11]:


print(MovingIm)


# In[14]:


print(RegIm)


# ### Plot the results:

# # Interactive feature doesn't seem to work:
# gui.MultiImageDisplay(image_list = [MovingIm, RegIm],                   
#                       title_list = ['Moving', 'Registered'], figure_size=(8,4));

# # Interactive feature doesn't seem to work:
# gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
#                                     known_transformation=ElastixImFilt);

# In[12]:


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


# In[11]:


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


# In[74]:


ImFilter.PrintParameterMap()


# ### Transform the contour points:
# 
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[12]:


# Get the transform parameter map:
TransformParameterMap = ElastixImFilt.GetParameterMap()

# Initiate TransformixImageFilter:
TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)
# Transform the contour points:
TransformixImFilt.Execute()


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

#[print(f'\nPoint {PtNos[i]}:\nInputIndex = {InInds[i]}\nFixedOutputIndex = {FixOutInds[i]}') for i in range(len(PtNos))]

for i in range(len(PtNos)):
    print(f'\nPoint {PtNos[i]}:')
    print(f'InputIndex        = {InInds[i]}')
    print(f'InputPoint        = {InPts[i]}')
    print(f'FixedOutputIndex  = {FixOutInds[i]}')
    print(f'OutputPoint       = {OutPts[i]}')
    print(f'Deformation       = {Defs[i]}')
    print(f'MovingOutputIndex = {MovOutInds[i]}')


# In[18]:


print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)
print('RegIm Size    =', RegSize)


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

# In[ ]:




