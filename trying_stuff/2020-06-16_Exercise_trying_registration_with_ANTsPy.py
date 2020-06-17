#!/usr/bin/env python
# coding: utf-8

# # Guidance
# 
# 
# ## Data
# 
# This is a stripped down version of previous version "2020-06-05_Exercise_trying_3D_registration_using_SimpleElastix.ipynb" to make it easier to follow.  For details of what was tried, why, what worked/failed, etc. refer to the May 22 version.
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
import dicom2nifti
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')

from ipywidgets import interact, fixed
from IPython.display import clear_output

from GetDicomFpaths import GetDicomFpaths
from GetDicoms import GetDicoms
from GetContourPtsForDicoms import GetContourPtsForDicoms
from ParseTransformixOutput import ParseTransformixOutput
from CreateDictionaryOfPoints import CreateDictionaryOfPoints
from CreateInputFileForElastix import CreateInputFileForElastix
from ParseTransformixOutput import ParseTransformixOutput
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf


# In[7]:


# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm'
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Get contour points into necessary arrays and dictionaries:
FixPts_PCS, FixPts_PCS_Arr, FixPts_PCS_Dict,FixPts_ICS, FixPts_ICS_Arr, FixPts_ICS_Dict = GetContourPtsForDicoms(DicomDir=FixDicomDir,
                                                                     RoiFpath=FixRoiFpath)

# Convert the dictionaries to dataframes:
FixPts_PCS_Df = pd.DataFrame(data=FixPts_PCS_Dict)
FixPts_ICS_Df = pd.DataFrame(data=FixPts_ICS_Dict)


# Read in the 3D stack of Fixed DICOMs:
FixReader = sitk.ImageSeriesReader()
FixNames = FixReader.GetGDCMSeriesFileNames(FixDicomDir)
FixReader.SetFileNames(FixNames)
FixIm_sitk = FixReader.Execute()
FixIm_sitk = sitk.Cast(FixIm_sitk, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovReader = sitk.ImageSeriesReader()
MovNames = MovReader.GetGDCMSeriesFileNames(MovDicomDir)
MovReader.SetFileNames(MovNames)
MovIm_sitk = MovReader.Execute()
MovIm_sitk = sitk.Cast(MovIm_sitk, sitk.sitkFloat32)

# Get some info on FixIm and MovIm:
#ruf.ShowImagesInfo(FixIm_sitk, MovIm_sitk)

# Convert from sitk to nifti:
FixNiftiFname = 'FixIm_itk2nii.nii'
MovNiftiFname = 'MovIm_itk2nii.nii'

writer = sitk.ImageFileWriter()
writer.SetFileName(FixNiftiFname)
writer.Execute(FixIm_sitk)

writer = sitk.ImageFileWriter()
writer.SetFileName(MovNiftiFname)
writer.Execute(MovIm_sitk)

# Read sitk2nii images into ants:
FixIm_ants = ants.image_read(FixNiftiFname)
MovIm_ants = ants.image_read(MovNiftiFname)


# Start timing:
times = []
times.append(time.time())

# Register MovIm_ants to FixIm_ants:
reg = ants.registration(FixIm_ants, MovIm_ants, 'Affine')

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

RegIm_ants = reg['warpedmovout']


verbose = True
verbose = False

# Apply the forward transformation to the FixPtsDFrame_PCS and FixPtsDFrame_ICS:
FwdTxPts_PCS_Df = ants.apply_transforms_to_points(dim=3, 
                                                 points=FixPtsDFrame_PCS, 
                                                 transformlist=reg['fwdtransforms'], 
                                                 verbose=verbose)

FwdTxPts_ICS_Df = ants.apply_transforms_to_points(dim=3, 
                                                 points=FixPtsDFrame_ICS, 
                                                 transformlist=reg['fwdtransforms'], 
                                                 verbose=verbose)

# Apply the inverse transformation to the same inputs:
RevTxPts_PCS_Df = ants.apply_transforms_to_points(dim=3, 
                                                 points=FixPtsDFrame_PCS, 
                                                 transformlist=reg['invtransforms'], 
                                                 verbose=verbose)

RevTxPts_ICS_Df = ants.apply_transforms_to_points(dim=3, 
                                                 points=FixPtsDFrame_ICS, 
                                                 transformlist=reg['invtransforms'], 
                                                 verbose=verbose)

# Convert from dataframes to lists:
FwdTxPts_PCS = FwdTxPts_PCS_Df.values.tolist()

FwdTxPts_ICS = FwdTxPts_ICS_Df.values.tolist()

RevTxPts_PCS = RevTxPts_PCS_Df.values.tolist()

RevTxPts_ICS = RevTxPts_ICS_Df.values.tolist()


# In[11]:


FwdTxPts_PCS


# In[12]:


FwdTxPts_ICS


# In[13]:


# Create a dictionary to store the output:
OutPtsDict = {'FixPts_PCS':FixPts_PCS,
              'FixPts_ICS':FixPts_ICS,
              'FwdTxPts_PCS':FwdTxPts_PCS,
              'FwdTxPts_ICS':FwdTxPts_ICS,
              'RevTxPts_PCS':RevTxPts_PCS,
              'RevTxPts_ICS':RevTxPts_ICS
             }

#OutPtsDict

# Convert the dictionary to a dataframe:
OutPtsDf = pd.DataFrame(data=OutPtsDict)

OutPtsDf[0:20]


# # 17/06:
# 
# ## The issue that I see with ANTsPy is that unlike SimpleITK, it only seems to provide the transformed points - not the indices, so working out where each point belongs on the image stack will be quite tedious.

# interact(ruf.display_images_with_contours,
#          fix_ind = (0,FixedIm_sitk.GetSize()[2]-1),
#          mov_ind = (0,MovingIm_sitk.GetSize()[2]-1),
#          fix_npa = fixed(sitk.GetArrayFromImage(FixedIm_sitk)),
#          mov_npa = fixed(sitk.GetArrayFromImage(MovingIm_sitk)),
#          fix_im = fixed(FixedIm_sitk),
#          mov_im = fixed(MovingIm_sitk),
#          mov_contour_pts = fixed(MovingContourPts),
#          fix_title = 'Fixed image',
#          mov_title = 'Moving image');

# print(FixedIm_ants.origin)
# print('\n', FixedIm_ants.range)

# import RegUtilityFuncs
# importlib.reload(RegUtilityFuncs)
# import RegUtilityFuncs as ruf
# 
# interact(ruf.display_ants_images_and_reg_result,
#          fix_ind = (0,FixedIm_ants.shape[2]-1),
#          mov_ind = (0,MovingIm_ants.shape[2]-1),
#          reg_ind = (0,RegIm_ants.shape[2]-1),
#          fix_im = fixed(FixedIm_ants),
#          mov_im = fixed(MovingIm_ants),
#          reg_im = fixed(RegIm_ants));

# import RegUtilityFuncs
# importlib.reload(RegUtilityFuncs)
# import RegUtilityFuncs as ruf
# 
# interact(ruf.display_ants_images_and_reg_result_with_contours,
#          fix_ind = (0,FixedIm_ants.shape[2]-1),
#          mov_ind = (0,MovingIm_ants.shape[2]-1),
#          reg_ind = (0,RegIm_ants.shape[2]-1),
#          fix_im = fixed(FixedIm_ants),
#          mov_im = fixed(MovingIm_ants),
#          reg_im = fixed(RegIm_ants),
#          mov_contour_pts=fixed(MovPtsArr_ICS));

# In[17]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

export_fig = False
#export_fig = True

ruf.display_all_sitk_images(fix_im=FixIm_sitk, mov_im=MovIm_sitk, export_fig=export_fig)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# # "Importantly, point mapping goes the opposite direction of image mapping, for both reasons of convention and engineering."
# 
# ### Link:
# 
# #### https://antspyx.readthedocs.io/en/latest/ants.registration.html?highlight=points#ants.registration.apply_transforms.apply_transforms_to_points
# 
# 

# In[9]:


FwdTxPts_PCS


# In[10]:


FwdTxPts_ICS


# In[11]:


RevTxPts_PCS


# In[12]:


RevTxPts_ICS


# ### The points obtained using the forward and reverse transforms look the same to me!  
# 
# ### 16/06:  Check if this is still true: Also the results using points in the PCS/ICS look similar to the results obtained using SimpleElastix.

# In[54]:


FixedIm_ants.origin


# In[55]:


MovingIm_ants.origin


# In[57]:


# Read in the DICOMs:
FixedDicomFpaths, FixedDicoms = GetDicoms(DirPath=FixedDicomDir, SortMethod='slices', Debug=False)
MovingDicomFpaths, MovingDicoms = GetDicoms(DirPath=MovingDicomDir, SortMethod='slices', Debug=False)


# In[70]:


FixedIPP = []

for Dicom in FixedDicoms:
    IPP = Dicom.ImagePositionPatient
    
    IPP = [float(value) for value in IPP]
    
    #print(IPP)
    FixedIPP.append(IPP)
    
MovingIPP = []

for Dicom in MovingDicoms:
    MovingIPP.append([float(value) for value in Dicom.ImagePositionPatient])
    
FixedIPP


# In[71]:


MovingIPP


# In[72]:


MovingIPP[13]


# In[84]:


MovPtsArr_PCS[13]


# In[86]:


MovPtsArr_ICS[13]


# In[89]:


MovingIPP


# In[88]:


MovingIm_ants.origin


# In[92]:


MovingIm_sitk.GetOrigin()


# In[94]:


MovingIm_sitk_from_nii = sitk.ReadImage(MovingNiftiFpath)
MovingIm_sitk_from_nii.GetOrigin()


# ### Same
# 
# ### Try writing to nifti using sitk:

# In[110]:


writer = sitk.ImageFileWriter()
writer.SetFileName('MovingIm_itk2nii.nii')
writer.Execute(MovingIm_sitk)


# In[111]:


# Read sitk2nii image:
MovingIm_ants2 = ants.image_read('MovingIm_itk2nii.nii')

MovingIm_ants2.origin


# In[112]:


MovingIm_sitk.GetOrigin()


# In[113]:


MovingDicoms[0].ImagePositionPatient


# ### The difference between the IPP from Pydicom and Origin from sitk:

# In[118]:


print([MovingIPP[0][i] - MovingIm_sitk.GetOrigin()[i] for i in range(3)])


# # Next repeat the work done above (loading of ants images, registration and  transforming the contour points) to see if they make more sense now.

# In[125]:


# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixedDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')
MovingDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')

# Define the filepath to the ROI Collection for the moving image:
MovingRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
MovingRoiFname = r'AIM_20200511_073405.dcm'
MovingRoiFpath = os.path.join(MovingRoiDir, MovingRoiFname)

# Get contour points into necessary arrays and dictionaries:
MovPts_PCS, MovPtsArr_PCS, MovPtsDict_PCS,MovPts_ICS, MovPtsArr_ICS, MovPtsDict_ICS = GetContourPtsForDicoms(DicomDir=MovingDicomDir,
                                                                   RoiFpath=MovingRoiFpath)

# Convert the dictionaries to dataframes:
MovPtsDFrame_PCS = pd.DataFrame(data=MovPtsDict_PCS)
MovPtsDFrame_ICS = pd.DataFrame(data=MovPtsDict_ICS)


# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDicomDir)
FixedReader.SetFileNames(FixedNames)
FixedIm_sitk = FixedReader.Execute()
FixedIm_sitk = sitk.Cast(FixedIm_sitk, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDicomDir)
MovingReader.SetFileNames(MovingNames)
MovingIm_sitk = MovingReader.Execute()
MovingIm_sitk = sitk.Cast(MovingIm_sitk, sitk.sitkFloat32)

# Convert from sitk to nifti:
FixedNiftiFname = 'FixedIm_itk2nii.nii'
MovingNiftiFname = 'MovingIm_itk2nii.nii'

writer = sitk.ImageFileWriter()
writer.SetFileName(FixedNiftiFname)
writer.Execute(FixedIm_sitk)

writer = sitk.ImageFileWriter()
writer.SetFileName(MovingNiftiFname)
writer.Execute(MovingIm_sitk)

# Read sitk2nii images into ants:
FixedIm_ants = ants.image_read(FixedNiftiFname)
MovingIm_ants = ants.image_read(MovingNiftiFname)

# Register MovingIm_ants to FixedIm_ants:
# Start timing:
times = []
times.append(time.time())


# Register in 3D this time:
reg = ants.registration(FixedIm_ants, MovingIm_ants, 'Affine')

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

RegIm_ants = reg['warpedmovout']


# Apply the forward transformation to the MovingContourPtsDFrame_PCS and MovingContourPtsDFrame_ICS:
FwdTxPts_PCS = ants.apply_transforms_to_points(dim=3, 
                                               points=MovPtsDFrame_PCS, 
                                               transformlist=reg['fwdtransforms'], 
                                               verbose=True)

FwdTxPts_ICS = ants.apply_transforms_to_points(dim=3, 
                                               points=MovPtsDFrame_ICS, 
                                               transformlist=reg['fwdtransforms'], 
                                               verbose=True)

# Apply the inverse transformation to the same inputs:
RevTxPts_PCS = ants.apply_transforms_to_points(dim=3, 
                                               points=MovPtsDFrame_PCS, 
                                               transformlist=reg['invtransforms'], 
                                               verbose=True)

RevTxPts_ICS = ants.apply_transforms_to_points(dim=3, 
                                               points=MovPtsDFrame_ICS, 
                                               transformlist=reg['invtransforms'], 
                                               verbose=True)


# In[126]:


FwdTxPts_PCS


# In[127]:


FwdTxPts_ICS


# In[128]:


RevTxPts_PCS


# In[129]:


RevTxPts_ICS


# ### Frustratingly they don't seem to be very different! 

# In[136]:


MovingIPP


# In[139]:


# The first slice in MovingIm that has a contour is the 13th slice, with IPP:
MovingIPP[13]


# In[140]:


MovPtsArr_PCS[13]


# In[141]:


MovPtsArr_ICS[13]


# In[148]:


# The position of the first contour point in the ICS w.r.t. IPP is:
print([MovPtsArr_ICS[13][0][i] - MovingIPP[13][i] for i in range(3)])


# In[149]:


# The number of pixel positions are:
print(MovingIm_ants.spacing)
print('\n', [(MovPtsArr_ICS[13][0][i] - MovingIPP[13][i])/MovingIm_ants.spacing[i] for i in range(2)])


# In[150]:


MovingIm_ants.shape


# In[132]:


FixedIPP


# In[135]:


FixPtsDFrame_PCS


# In[ ]:


### The transformed points seem like they're not in 


# In[ ]:




