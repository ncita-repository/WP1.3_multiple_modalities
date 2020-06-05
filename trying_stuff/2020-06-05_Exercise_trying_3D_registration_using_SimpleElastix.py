#!/usr/bin/env python
# coding: utf-8

# # Guidance
# 
# 
# ## Data
# 
# This is a stripped down version of previous version "2020-05-22_Exercise_trying_3D_registration_using_SimpleElastix.ipynb" to make it easier to follow.  For details of what was tried, why, what worked/failed, etc. refer to the May 22 version.
# 
# The data used is from Subject 011 from the ACRIN-FMISO-Brain data set obtained here:
# 
# https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948305#a6b63241cb4843a291c45ef8a62d2073
# 
# You can selectively download data for that subject only.  Specifically two series from two sessions were used:
# 
# 8-T1 SE AXIAL POST FS FC-59362 (Series 8) within 04-10-1960-MRI Brain wwo Contrast-69626
# 
# and
# 
# 8-T1 SE AXIAL POST FC-81428 (Series 8) within 06-11-1961-MRI Brain wwo Contrast-79433
# 
# An ROI Collection was manually created within the OHIF-Viewer for the series from 1960.  The DICOM-RTSTRUCT file can be obtained here:
# 
# https://drive.google.com/drive/folders/1CDrwTYvlvBL6rpuyz29zxO6AU4Ancmd6?usp=sharing
# 
# Originally the 1960 series (with ROI Collection) was defined as the "fixed" image, and the 1961 series as "moving".  This is fine but requires the inverse transformation used to register moving to fixed to deform the fixed image's contour points to the moving image.
# 
# I'm sure this can be done but I don't yet know how, so in the meantime fixed will be defined as the 1961 series, and moving the 1960 one.
# 
# ## SimpleElastix
# 
# I recommend trying to install SimpleElastix using a pre-compiled Python module obtained here:
# 
# https://www.francescosantini.com/wp/2020/03/06/simpleelastix-precompiled-python-modules/
# 
# If it doesn't work you will have to compile it from source:
# 
# https://simpleelastix.readthedocs.io/GettingStarted.html
# 
# There's a "more official" wheel provided here:
# 
# https://sourceforge.net/projects/simpleelastix/
# 
# but it might be for an older version of SimpleElastix (and I think I had trouble with this wheel and ended up having to compile it - I hadn't yet discovered Francesco Santini's page).
# 
# ## Other packages
# 
# Most of the other packages can be pip installed straightforwardly.  
# 
# The other imports are functions that can be obtained here:
# 
# https://github.com/ncita-repository/WP1.3_multiple_modalities/tree/master/trying_stuff

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

from GetDicomFpaths import GetDicomFpaths
from GetDicoms import GetDicoms
from Get3DContourPointsFromDicoms import Get3DContourPointsFromDicoms
from Get2DContourPointsFromDicoms import Get2DContourPointsFromDicoms
import Create3DInputFileForElastix
importlib.reload(Create3DInputFileForElastix)
from Create3DInputFileForElastix import Create3DInputFileForElastix
from Create2DInputFileForElastix import Create2DInputFileForElastix
from ParseTransformixOutput import ParseTransformixOutput
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf


# In[2]:


# Define FixedDicomDir and MovingDicomDir:
FixedDicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011\06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428'
MovingDicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011\04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362'

# Define the filepath to the ROI Collection for the moving image:
MovingRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
MovingRoiFname = r'AIM_20200511_073405.dcm'
MovingRoiFpath = os.path.join(MovingRoiDir, MovingRoiFname)


# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDicomDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()
FixedIm = sitk.Cast(FixedIm, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDicomDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()
MovingIm = sitk.Cast(MovingIm, sitk.sitkFloat32)

# Get some info on FixedIm and MovingIm:
ruf.ShowImagesInfo(FixedIm, MovingIm)

# Get the contour points:
MovingContourPts = Get3DContourPointsFromDicoms(DicomDir=MovingDicomDir,
                                                RoiFpath=MovingRoiFpath)

print('\nMovingContourPts:\n')
MovingContourPts


# ### Display images using interactive plot:

# gui.MultiImageDisplay(image_list = [FixedIm, MovingIm],                   
#                       title_list = ['Fixed', 'Moving'], figure_size=(8,4));

# In[3]:


#interact(ruf.display_images,
#         fix_ind=(0,FixedIm.GetSize()[2]-1),
#         mov_ind=(0,MovingIm.GetSize()[2]-1),
#         fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
#         mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
#         fix_title = 'Fixed image',
#         mov_title = 'Moving image');

interact(ruf.display_images_with_contours,
         fix_ind=(0,FixedIm.GetSize()[2]-1),
         mov_ind=(0,MovingIm.GetSize()[2]-1),
         fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         fix_im = fixed(FixedIm),
         mov_im = fixed(MovingIm),
         mov_contour_pts = fixed(MovingContourPts),
         fix_title = 'Fixed image',
         mov_title = 'Moving image');


# ### Create inputpoints.txt for Elastix, containing the contour points to be transformed:

# In[21]:


# Create inputpoints.txt for Elastix, containing the 3D contour points to be
# transformed:
Create3DInputFileForElastix(DicomDir=MovingDicomDir,
                            RoiFpath=MovingRoiFpath)


# ### View the created file inputpoints.pcs

# In[23]:


# Open the text file:
#TextFile = open('inputpoints.pcs', 'r')
TextFile = open('inputpoints.txt', 'r')

for line in TextFile:
    print(line)


# ### Register MovingIm to FixedIm using SimpleElastix:

# In[55]:


RegIm, ElastixImFilt = ruf.RunSimpleElastixReg(FixedIm, MovingIm)


# ### Get some info on FixedIm, MovingIm and RegIm:

# In[35]:


# Get some info on FixedIm, MovingIm and RegIm:

ruf.ShowRegResults(FixedIm, MovingIm, RegIm)


# In[10]:


print(FixedIm)


# In[11]:


print(MovingIm)


# In[12]:


print(RegIm)


# ### Plot the results:

# In[50]:


#interact(ruf.display_sitk_images_and_reg_result,
#         fix_ind=(0,FixedIm.GetSize()[2]-1),
#         mov_ind=(0,MovingIm.GetSize()[2]-1),
#         reg_ind=(0,RegIm.GetSize()[2]-1),
#         fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
#         mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
#         reg_npa = fixed(sitk.GetArrayFromImage(RegIm)),
#         fix_title = 'Fixed image',
#         mov_title = 'Moving image',
#         reg_title = 'Registered image');

interact(ruf.display_sitk_images_and_reg_result_with_contours,
         fix_ind = (0,FixedIm.GetSize()[2]-1),
         mov_ind = (0,MovingIm.GetSize()[2]-1),
         reg_ind = (0,RegIm.GetSize()[2]-1),
         fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         reg_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fix_im = fixed(FixedIm),
         mov_im = fixed(MovingIm),
         reg_im = fixed(RegIm),
         mov_contour_pts = fixed(MovingContourPts),
         fix_title = 'Fixed image',
         mov_title = 'Moving image',
         reg_title = 'Registered image');


# ImFilter.PrintParameterMap()

# ### Transform the contour points:
# 
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[61]:


# Transform MovingContourPts:
ruf.TransformPoints(MovingImage=MovingIm, TransformFilter=ElastixImFilt)


# In[4]:


# Parse outputpoints.txt:
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

print('\nParsed OutputPoints.txt:\n')
for i in range(len(PtNos)):
    print(f'\nPoint {PtNos[i]}:')
    print(f'InputIndex        = {InInds[i]}')
    print(f'InputPoint        = {InPts[i]}')
    print(f'FixedOutputIndex  = {FixOutInds[i]}')
    print(f'OutputPoint       = {OutPts[i]}')
    print(f'Deformation       = {Defs[i]}')
    print(f'MovingOutputIndex = {MovOutInds[i]}')


# In[ ]:


# Compute the deformation field:
""" Running this in Anaconda Jupyter will probably kill the kernel.  If so, run in PyCharm. """
DefField = ruf.ComputeDeformationField(MovingImage=MovingIm, TransformFilter=ElastixImFilt)

