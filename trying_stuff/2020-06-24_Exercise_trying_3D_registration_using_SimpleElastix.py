#!/usr/bin/env python
# coding: utf-8

# # Guidance
# 
# 
# ## Data
# 
# This is a stripped down version of previous version "2020-06-16_Exercise_trying_3D_registration_using_SimpleElastix.ipynb" to make it easier to follow.  For details of what was tried, why, what worked/failed, etc. refer to the May 22 version.
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

# 

# In[2]:


# Import packages and functions
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
import GetContourPtsForDicoms
importlib.reload(GetContourPtsForDicoms)
from GetContourPtsForDicoms import GetContourPtsForDicoms
from GetInputPoints import GetInputPoints
#import CreateInputFileForElastix
#importlib.reload(CreateInputFileForElastix)
from CreateInputFileForElastix import CreateInputFileForElastix
from ParseTransformixOutput import ParseTransformixOutput
import GetImageAttributes
importlib.reload(GetImageAttributes)
from GetImageAttributes import GetImageAttributes
import PCStoICS
importlib.reload(PCStoICS)
from PCStoICS import PCStoICS
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
from GetOutputPoints import GetOutputPoints
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm'
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose which package to use to get the Image Plane Attributes:
package = 'sitk'
package = 'pydicom'


# Get the Image Attributes for the images:
FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                              Package=package)

MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                              Package=package)

# Get contour points into necessary arrays and dictionaries:
#FixPts_PCS, FixPtsArr_PCS, FixPtsDict_PCS = GetContourPtsForDicoms(DicomDir=FixDicomDir,
#                                                                   RoiFpath=FixRoiFpath)

InPts_PCS, InPtsArr_PCS, InPtsDict_PCS,InPts_ICS, InPtsArr_ICS, InPtsDict_ICS = GetInputPoints(DicomDir=FixDicomDir, 
                                                        RoiFpath=FixRoiFpath,
                                                        Origin=FixOrigin,
                                                        Directions=FixDirs,
                                                        Spacings=FixSpacings)

# Read in the 3D stack of Fixed DICOMs:
FixReader = sitk.ImageSeriesReader()
FixNames = FixReader.GetGDCMSeriesFileNames(FixDicomDir)
FixReader.SetFileNames(FixNames)
FixIm = FixReader.Execute()
FixIm = sitk.Cast(FixIm, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovReader = sitk.ImageSeriesReader()
MovNames = MovReader.GetGDCMSeriesFileNames(MovDicomDir)
MovReader.SetFileNames(MovNames)
MovIm = MovReader.Execute()
MovIm = sitk.Cast(MovIm, sitk.sitkFloat32)

# Get some info on FixIm and MovIm:
#ruf.ShowImagesInfo(FixIm, MovIm)

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
CreateInputFileForElastix(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath, 
                          Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings,
                          CoordSys=CoordSys)

# Register:
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)

# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
#PtNos, InInds, InPts_PCS, FixOutInds,\
#OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

PtNos, InInds, InPts_PCS, InPts_ICS,FixOutInds, OutPts_PCS, OutPts_ICS,Defs_PCS, Defs_ICS, MovOutInds,OutPtsArr_PCS, OutPtsArr_ICS = GetOutputPoints(Origin=MovOrigin, Directions=MovDirs,
                                               Spacings=MovSpacings, Dimensions=MovDims)

print('\nDone.')


# In[3]:


# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]

# Choose whether to export figure:
export_fig = True
export_fig = False


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_June17_Attempt2_eqs.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                              fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                              export_fig=export_fig,
                                                              export_fname=fig_fname,
                                                              plot_slices=plot_slices)


# ## Success!
# 
# ## Next steps:
# 
# ## 1) Store the transformed points into an RTSTRUCT file for importing into the OHIF-Viewer to see how they appear there.
# 
# ## 2) Tidy up the contours:
# 
# ### i.   Some/most/all contours require "closing"
# ### ii.  Some require removal of un-closed segments (e.g. 18th slice)
# ### iii. Some require avoiding segments crossing through the outer boundary (e.g. 15-17th slice)
# ### iv. Some require some awkward closing (e.g. 13th slice)

# ## Step 1:  Create ROI for Moving Image
# 
# ### This is a bit trickier than simply modifying a copy of the fixed ROI, since there are more contour sequences within the transformed points than in the original ROI.

# In[ ]:


# Use the fixed ROI as a template for the moving image:
FixRoi = pydicom.dcmread(FixRoiFpath)

MovRoi = copy.deepcopy(FixRoi)

#MovRoi
#MovRoi.PatientName


# In[181]:


print(len(OutPtsArr_ICS))
print(MovIm.GetSize()[2])
print(len(MovDicoms))
print()


# In[ ]:


""" MovDicoms needs to be indexed by MovZinds after re-ordering... """


# In[13]:


# Read in the DICOMs:
FixDicomFpaths, FixDicoms = GetDicoms(DirPath=FixDicomDir, SortMethod='slices', Debug=False)

MovDicomFpaths, MovDicoms = GetDicoms(DirPath=MovDicomDir, SortMethod='slices', Debug=False)

FixDicoms[0]


# In[39]:


MovDicoms[0]


# In[122]:


print('FixRoi:')
FixRoi


# In[232]:


print('MovRoi:')
MovRoi


# In[42]:


import CreateMovRoi
importlib.reload(CreateMovRoi)
from CreateMovRoi import CreateMovRoi

MovRoi = CreateMovRoi(FixRoiFpath=FixRoiFpath, MovDicomDir=MovDicomDir, 
                      MovPtsArr_ICS=OutPtsArr_ICS, Debug=False)


# In[43]:


MovRoi


# # Sagittal views of Fixed and Moving Images

# FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, Package='pydicom')
# MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, Package='pydicom')
# 
# print('FixDims:', FixDims)
# print('\nMovDims:', MovDims)
# 
# print('\nFixDirs:\n', FixDirs[0:3], '\n', FixDirs[3:6], '\n', FixDirs[6:])
# print('\nMovDirs:\n', MovDirs[0:3], '\n', MovDirs[3:6], '\n', MovDirs[6:])
# 
# # Check the aspect ratios in different views:
# FixAxAR = FixSpacings[1]/FixSpacings[0]
# MovAxAR = MovSpacings[1]/MovSpacings[0]
# 
# FixSagAR = FixSpacings[1]/FixSpacings[2]
# MovSagAR = MovSpacings[1]/MovSpacings[2]
# #FixSagAR = FixSpacings[2]/FixSpacings[1]
# #MovSagAR = MovSpacings[2]/MovSpacings[1]
# 
# FixCorAR = FixSpacings[2]/FixSpacings[0]
# MovCorAR = MovSpacings[2]/MovSpacings[0]
# 
# print('\nFixAxAR =', FixAxAR)
# print('FixSagAR =', FixSagAR)
# print('FixCorAR =', FixCorAR)
# print('\nMovAxAR =', MovAxAR)
# print('MovSagAR =', MovSagAR)
# print('MovCorAR =', MovCorAR)
# 
# # Create 3D array of DICOM images:
# Fix3D = np.zeros(FixDims)
# Mov3D = np.zeros(MovDims)
# 
# for i in range(len(FixDicoms)):
#     Fix3D[:, :, i] = FixDicoms[i].pixel_array
# 
# for i in range(len(MovDicoms)):
#     Mov3D[:, :, i] = MovDicoms[i].pixel_array
#     
# 

# # plot 3 orthogonal slices
# plt.subplots(1, 3, figsize=(14, 5))
# 
# a1 = plt.subplot(1, 3, 1)
# plt.imshow(Fix3D[:, :, FixDims[2]//2])
# a1.set_aspect(FixAxAR)
# 
# a2 = plt.subplot(1, 3, 2)
# plt.imshow(Fix3D[:, FixDims[1]//2, :])
# a2.set_aspect(FixSagAR)
# 
# a3 = plt.subplot(1, 3, 3)
# plt.imshow(Fix3D[FixDims[0]//2, :, :].T)
# a3.set_aspect(FixCorAR)
# 
# plt.show()

# # Sagittal plots:
# 
# Nrows = FixDims[2] # same as the number of slices
# Nrows = FixDims[2]//2
# Ncols = 1
# 
# # Define the slice thickness along sagital cuts:
# #dSag = int(np.floor(FixDims[0]/FixDims[2])) # same as below
# dSag = FixDims[0]//FixDims[2] # same as above
# 
# # Create list of indices for the sagittal slices:
# #yInds = np.arange(0, FixDims[0]+1, dSag) # better to use linspace..
# yInds = np.linspace(0, FixDims[1]-1, num=Nrows, endpoint=True, retstep=False, dtype=int)
# 
# 
# Rotate = False
# Rotate = True
# 
# plt.subplots(Nrows, Ncols, figsize=(5, 5*Nrows))
# #plt.subplots(Nrows, Ncols, figsize=(5, 50))
# 
# i = 0
# 
# for s in range(Nrows):
#     i = i + 1 # increment sub-plot number
#     
#     #print(f'\nyInd = {yInds[s]}')
#     
#     if Rotate:
#         arr = np.rot90(Fix3D[:, yInds[s], :])
#     else:
#         arr = Fix3D[:, yInds[s], :]
#     
#     ax1 = plt.subplot(Nrows, Ncols, i)
#     plt.imshow(arr, cmap=plt.cm.Greys_r);
#     plt.title(f'Fixed slice {s+1}/{Nrows}')
#     plt.axis('off')
#     if Rotate:
#         ax1.set_aspect(1/FixSagAR)
#     else:
#         ax1.set_aspect(FixSagAR)
# 

# In[219]:


FixNda = sitk.GetArrayFromImage(FixIm)

print('FixDims         =', FixDims)
print('FixNda.shape    =', FixNda.shape)
print('FixIm.GetSize() =', FixIm.GetSize())
print('')
print('FixDirs              =\n', FixDirs)
print('\nFixIm.GetDirection() =\n', FixIm.GetDirection())
print('')
print(FixDirs[0:2])
print([FixDirs[0], FixDirs[2]])
print(FixDirs[1:3])

print('\n', type(FixIm.GetDirection()))
#print('\n', type(list(FixIm.GetDirection())))
print('\n', type(FixDirs))

#FixDirs = FixIm.GetDirection()
#FixDirs = list(FixIm.GetDirection())
FixDirs = [FixIm.GetDirection()[i] for i in range(len(FixIm.GetDirection()))]

print('\nFixDirs =\n', FixDirs)

a = [1, 2, 3, 5, 6, 7]
#a = list([1, 2, 3, 5, 6, 7])

for i in a:
    FixDirs[i] = - FixDirs[i]

print('\nFixDirs =\n', FixDirs)


# In[43]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
#plot_slices = [12, 13, 14, 15]


# Choose the perspective:
perspective = 'axial'
#perspective = 'sagittal'
#perspective = 'coronal'

# Choose whether to export figure:
export_fig = True
export_fig = False

# Choose whether to plot contour points as lines or dots:
contours_as = 'lines'
contours_as = 'dots'


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_' + perspective + '.png'

#ruf.display_all_sitk_images_and_reg_results_with_all_contours_v2(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
#                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
#                                                                 export_fig=export_fig,
#                                                                 export_fname=fig_fname,
#                                                                 plot_slices=plot_slices,
#                                                                 perspective=perspective)

ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                                 export_fig=export_fig,
                                                                 export_fname=fig_fname,
                                                                 plot_slices=plot_slices,
                                                                 perspective=perspective,
                                                                 contours_as=contours_as)


# ### June 24:  Done producing sagittal (and coronal).
# 
# ### Next work on tidying up the contours.

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




