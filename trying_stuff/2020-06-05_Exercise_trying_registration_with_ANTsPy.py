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
# ## ANTsPy
# 
# Official guidance on installing ANTsPy:
# 
# https://github.com/ANTsX/ANTsPy/blob/master/tutorials/InstallingANTsPy.md
# 
# I was unsuccessful with any of the methods listed.  But was able to install a pre-compiled wheel using:
# 
# pip install https://github.com/SGotla/ANTsPy/releases/download/0.1.7Win64/antspy-0.1.7-cp37-cp37m-win_amd64.whl
# 
# See for details:
# 
# https://github.com/ANTsX/ANTsPy/issues/28
# 
# 
# ## Other packages
# 
# Most of the other packages can be pip installed straightforwardly.  
# 
# The other imports are functions that can be obtained here:
# 
# https://github.com/ncita-repository/WP1.3_multiple_modalities/tree/master/trying_stuff

# In[1]:


import ants
print('ANTsPy version', ants.__version__)
import SimpleITK as sitk
print(sitk.Version())
import numpy as np
import time
import os
import copy
import importlib
import json
import pydicom
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')

from ipywidgets import interact, fixed
from IPython.display import clear_output

from GetDicomFpaths import GetDicomFpaths
from GetDicoms import GetDicoms
#from GetContourPoints3D import GetContourPoints3D
from Get3DContourPointsFromDicoms import Get3DContourPointsFromDicoms
from Create3DInputFileForElastix import Create3DInputFileForElastix
from Create2DInputFileForElastix import Create2DInputFileForElastix
from ParseTransformixOutput import ParseTransformixOutput
from CreateDictionaryOfPoints import CreateDictionaryOfPoints

import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf


# In[4]:


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


# In[5]:


interact(ruf.display_images_with_contours,
         fix_ind = (0,FixedIm.GetSize()[2]-1),
         mov_ind = (0,MovingIm.GetSize()[2]-1),
         fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         fix_im = fixed(FixedIm),
         mov_im = fixed(MovingIm),
         mov_contour_pts = fixed(MovingContourPts),
         fix_title = 'Fixed image',
         mov_title = 'Moving image');


# In[7]:


# Read in the DICOMs as a 3D ants image:

FixedIm = ants.dicom_read(FixedDicomDir)
MovingIm = ants.dicom_read(MovingDicomDir)


# In[8]:


print(FixedIm)


# In[9]:


print(MovingIm)


# In[10]:


# Start timing:
times = []
times.append(time.time())


# Register in 3D this time:
reg = ants.registration(FixedIm, MovingIm, 'Affine')

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

RegIm = reg['warpedmovout']


# ### Took only 6-9 s to register the 3D stacks.  By comparison, SimpleElastix took 23-33 s!!

# In[11]:


print(reg)


# In[12]:


FixedIm.shape
#FixedIm.numpy()


# In[20]:


FixedIm[14].shape


# In[21]:


FixedIm[14,:,:].shape


# In[22]:


FixedIm[:,:,14].shape


# In[11]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

interact(ruf.display_ants_images_and_reg_result,
         fix_ind = (0,FixedIm.shape[2]-1),
         mov_ind = (0,MovingIm.shape[2]-1),
         reg_ind = (0,RegIm.shape[2]-1),
         fix_npa = fixed(FixedIm.numpy()),
         mov_npa = fixed(MovingIm.numpy()),
         reg_npa = fixed(RegIm.numpy()),
         fix_im = fixed(FixedIm),
         mov_im = fixed(MovingIm),
         reg_im = fixed(RegIm),
         fix_title = 'Fixed image',
         mov_title = 'Moving image',
         reg_title = 'Registered image');


# ### Something is very wrong with these images.  When scanning through FixedIm the first slice is repeated at slice 11, again at 22,...
# 
# ### It looks like the DICOM reader of ANTsPy don't order the slices along the z axis!
# 
# ### I can't see how to load the images in the correct order, or how to use .read_image in a sequential fashion, using a list of ordered filepaths, so instead, use sitk to read the DICOMs, convert to numpy, then to ants.

# In[36]:


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


# Convert from sitk image to numpy array to ants image:
FixedImAnts = ants.from_numpy(data=sitk.GetArrayFromImage(FixedIm), 
                              origin=FixedIm.GetOrigin(),
                              spacing=FixedIm.GetSpacing(),
                              direction=FixedIm.GetDirection(),
                              has_components=False,
                              is_rgb=False)

FixedIAnts


# In[34]:


FixedIm.GetOrigin()


# ### 03/06: Can't figure out why I'm getting the error above, since FixedIm.GetOrigin() does have 3 dimensions...

# In[28]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

interact(ruf.display_ants_images_and_reg_result_with_contours,
         fixed_ind = (0,FixedIm.shape[2]-1),
         moving_ind = (0,MovingIm.shape[2]-1),
         reg_ind = (0,RegIm.shape[2]-1),
         fixed_npa = fixed(FixedIm.numpy()),
         moving_npa = fixed(MovingIm.numpy()),
         reg_npa = fixed(RegIm.numpy()),
         fixed_im = fixed(FixedIm),
         moving_im = fixed(MovingIm),
         reg_im = fixed(RegIm),
         fixed_contour_pts = fixed(FixedContourPts),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         reg_title = 'Registered image');


# In[45]:


# Create a dictionary of contour points:
FixedContourPtsDict = CreateDictionaryOfPoints(DicomDir=DataDict['Source']['DicomDir'],
                                               RoiFpath=DataDict['Source']['RoiFpath'])

FixedContourPtsDict


# In[47]:


# Convert to a dataframe:
FixedPts = pd.DataFrame(data=FixedContourPtsDict)

FixedPts


# In[56]:


# Apply the inverse transformation to the fixed contour points:
MovingPts = ants.apply_transforms_to_points(3, FixedPts, reg['invtransforms'], verbose=True)


# In[57]:


MovingPts


# ### The x and y coordinates look reasonable in that they don't seem to exceed the physical dims, but the z coordinates seem wrong:  z=-54.5 mm has been deformed to z=95.1 mm!!!

# In[58]:


# Apply the forward transformation to the moving contour points to see if the result is the same as FixedPts:
FixedPts_recovered = ants.apply_transforms_to_points(3, MovingPts, reg['fwdtransforms'], verbose=True)


# In[59]:


FixedPts_recovered


# In[51]:


FixedIm.origin


# In[52]:


MovingIm.origin


# In[ ]:




