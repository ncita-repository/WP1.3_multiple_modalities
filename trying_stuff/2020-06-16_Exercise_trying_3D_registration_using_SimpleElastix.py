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

# In[20]:


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
import PCStoICS
importlib.reload(PCStoICS)
from PCStoICS import PCStoICS
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf


# In[21]:


# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm'
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Get contour points into necessary arrays and dictionaries:
#FixPts_PCS, FixPtsArr_PCS, FixPtsDict_PCS = GetContourPtsForDicoms(DicomDir=FixDicomDir,
#                                                                   RoiFpath=FixRoiFpath)

InPts_PCS, InPtsArr_PCS, InPtsDict_PCS,InPts_ICS, InPtsArr_ICS, InPtsDict_ICS = GetInputPoints(DicomDir=FixDicomDir, 
                                                        RoiFpath=FixRoiFpath)

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
ruf.ShowImagesInfo(FixIm, MovIm)


# In[45]:


#MovPtsArr_PCS


# In[22]:


# Convert InPts_PCS from PCS to ICS:
InPts_ICS = PCStoICS(Pts_PCS=InPts_PCS, SitkIm=MovIm)
InPts_ICS


# In[5]:


InPtsArr_PCS


# In[6]:


InPtsArr_ICS


# In[7]:


InPtsArr_ICS


# ### Display images using interactive plot:

# gui.MultiImageDisplay(image_list = [FixedIm, MovingIm],                   
#                       title_list = ['Fixed', 'Moving'], figure_size=(8,4));

# if False:
#     interact(ruf.display_images,
#     fix_ind=(0,FixedIm.GetSize()[2]-1),
#     mov_ind=(0,MovingIm.GetSize()[2]-1),
#     fix_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
#     mov_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
#     fix_title = 'Fixed image',
#     mov_title = 'Moving image');
# 
# interact(ruf.display_images_with_contours,
#          fix_ind=(0,FixIm.GetSize()[2]-1),
#          mov_ind=(0,MovIm.GetSize()[2]-1),
#          fix_npa = fixed(sitk.GetArrayFromImage(FixIm)),
#          mov_npa = fixed(sitk.GetArrayFromImage(MovIm)),
#          fix_im = fixed(FixIm),
#          mov_im = fixed(MovIm),
#          mov_contour_pts = fixed(MovPtsArr_ICS),
#          fix_title = 'Fixed image',
#          mov_title = 'Moving image');

# ### Create inputpoints.txt for Elastix, containing the contour points to be transformed:

# In[23]:


import CreateInputFileForElastix
importlib.reload(CreateInputFileForElastix)
from CreateInputFileForElastix import CreateInputFileForElastix

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
CreateInputFileForElastix(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath, CoordSys=CoordSys)


# ### View the created file inputpoints.pcs

# In[11]:


# Open the text file:
#TextFile = open('inputpoints.txt', 'r')

#for line in TextFile:
#    print(line)


# ### Register MovingIm to FixedIm using SimpleElastix:

# In[4]:


RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)


# ### Get some info on FixIm, MovIm and RegIm:

# In[7]:


# Get some info on FixIm, MovIm and RegIm:

#ruf.ShowRegResults(FixIm, MovIm, RegIm)


# In[8]:


#print(FixedIm)


# In[9]:


#print(MovingIm)


# In[10]:


#print(RegIm)


# ### Plot the results:

# import RegUtilityFuncs
# importlib.reload(RegUtilityFuncs)
# import RegUtilityFuncs as ruf
# 
# 
# if False:
#     interact(ruf.display_sitk_images_and_reg_result,
#              fix_ind=(0,FixedIm.GetSize()[2]-1),
#              mov_ind=(0,MovingIm.GetSize()[2]-1),
#              reg_ind=(0,RegIm.GetSize()[2]-1),
#              fix_im = fixed(FixedIm),
#              mov_im = fixed(MovingIm),
#              reg_im = fixed(RegIm));
# 
# interact(ruf.display_sitk_images_and_reg_result_with_contours,
#          fix_ind = (0,FixIm.GetSize()[2]-1),
#          mov_ind = (0,MovIm.GetSize()[2]-1),
#          reg_ind = (0,RegIm.GetSize()[2]-1),
#          fix_im = fixed(FixIm),
#          mov_im = fixed(MovIm),
#          reg_im = fixed(RegIm),
#          mov_contour_pts = fixed(MovPtsArr_ICS));

# ImFilter.PrintParameterMap()

# ### Transform the contour points:
# 
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[24]:


# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)


# In[16]:


# Parse outputpoints.txt:
PtNos, InInds, InPts_PCS, FixOutInds,OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

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
    print(f'InputPoint_PCS    = {InPts_PCS[i]}')
    print(f'FixedOutputIndex  = {FixOutInds[i]}')
    print(f'OutputPoint_PCS   = {OutPts_PCS[i]}')
    print(f'Deformation_PCS   = {Defs_PCS[i]}')
    print(f'MovingOutputIndex = {MovOutInds[i]}')


# ### Convert output points from PCS to ICS

# In[22]:


OutPts_PCS[0]


# In[23]:


[OutPts_PCS[0]]


# In[25]:


from PCStoICS import PCStoICS

PCStoICS(Pts_PCS=[OutPts_PCS[0]], SitkIm=MovIm)


# In[23]:


PCStoICS(Pts_PCS=[OutPts_PCS[0]], SitkIm=MovIm)[0]


# In[25]:


import GetOutputPoints
importlib.reload(GetOutputPoints)
from GetOutputPoints import GetOutputPoints

# Convert InPts_PCS, OutPts_PCS and Defs_PCS to ICS:
#InPts_ICS = ruf.PCStoICS(InPts_PCS, MovIm)
#OutPts_ICS = ruf.PCStoICS(OutPts_PCS, MovIm)
#
#P = len(PtNos)
#Defs_ICS = [[OutPts_ICS[p][i] - InPts_ICS[p][i] for i in range(3)] for p in range(P)]

PtNos, InInds, InPts_PCS, InPts_ICS,FixOutInds, OutPts_PCS, OutPts_ICS,Defs_PCS, Defs_ICS, MovOutInds,OutPtsArr_PCS, OutPtsArr_ICS = GetOutputPoints(MovingIm=MovIm)

# Create a dictionary to store the output:
OutPtsDict = {'PtNo':PtNos,
              'InInd':InInds,
              'InPt_PCS':InPts_PCS,
              'InPt_ICS':InPts_ICS,
              'FixOutInd':FixOutInds,
              'OutPt_PCS':OutPts_PCS,
              'OutPt_ICS':OutPts_ICS,
              'Def_PCS':Defs_PCS,
              'Def_ICS':Defs_ICS,
              'MovOutInd':MovOutInds
             }

#OutPtsDict

# Convert the dictionary to a dataframe:
OutPtsDf = pd.DataFrame(data=OutPtsDict)

OutPtsDf[0:20]


# In[62]:


print(151.8/13.1)
print(130.3/11.7)
print(0.62/1.91)


# In[18]:


OutPtsArr_ICS


# In[17]:


### 15/06: Display result of:
# 1. Created inputpoints.txt in PCS (directly from ContourData)
# 2. Transformed points
# 3. Converted points given by OutputPoint from PCS to ICS (using 3D conversion)
# 4. Grouped converted points by z-index given by OutputIndexMoving

""" 
Although the points are not overlaying correctly, at least the contours were spread 
across a number of slices (compare to results obtained by deforming points in the ICS,
which resulted in most of the z-indices in OutputIndexMoving = 13, with some in 12 and 
some in 14).
This suggests that the points need to be transformed in the PCS.
But further work required to get the points overlayed properly.
"""

# Choose whether to export figure:
export_fig = True
export_fig = False


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_sitk.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                              fix_pts=FixPtsArr_ICS, mov_pts=TxPtsArr_ICS,
                                                              export_fig=export_fig,
                                                              export_fname=fig_fname)


# In[53]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

InputPoint_ICS = ruf.PCStoICS(InPts_PCS, MovIm)
OutputPoint_ICS = ruf.PCStoICS(OutPts_PCS, MovIm)

print(OutputPoint_ICS[0])
print(InputPoint_ICS[0])

[OutputPoint_ICS[0][i] - InputPoint_ICS[0][i] for i in range(3)]

[[OutputPoint_ICS[p][i] - InputPoint_ICS[p][i] for i in range(3)] for p in range(len(OutputPoint_ICS))]


# In[18]:


OutPtsArr_PCS


# In[26]:


OutPtsArr_ICS


# In[28]:


### 17/06: Display result of:
# 1. Created inputpoints.txt in PCS (directly from ContourData)
# 2. Transformed points
# 3. Converted points given by OutputPoint from PCS to ICS (using 3D conversion)
# 4. Grouped converted points by z-index given by OutputIndexMoving

### using latest equations derived on June 17 - Attempt #1.

""" 
Although the points are not overlaying correctly, at least the contours were spread 
across a number of slices (compare to results obtained by deforming points in the ICS,
which resulted in most of the z-indices in OutputIndexMoving = 13, with some in 12 and 
some in 14).
This suggests that the points need to be transformed in the PCS.
But further work required to get the points overlayed properly.
"""

import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Choose whether to export figure:
export_fig = True
#export_fig = False


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_sitk.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                              fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                              export_fig=export_fig,
                                                              export_fname=fig_fname)


# In[26]:


### 17/06: Display result of:
# 1. Created inputpoints.txt in PCS (directly from ContourData)
# 2. Transformed points
# 3. Converted points given by OutputPoint from PCS to ICS (using 3D conversion)
# 4. Grouped converted points by z-index given by OutputIndexMoving

### using latest equations derived on June 17 - Attempt #2.

""" 
Although the points are not overlaying correctly, at least the contours were spread 
across a number of slices (compare to results obtained by deforming points in the ICS,
which resulted in most of the z-indices in OutputIndexMoving = 13, with some in 12 and 
some in 14).
This suggests that the points need to be transformed in the PCS.
But further work required to get the points overlayed properly.
"""

import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Choose whether to export figure:
export_fig = True
#export_fig = False


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_sitk.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                              fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                              export_fig=export_fig,
                                                              export_fname=fig_fname)


# In[ ]:


print('FixIm origin  =', FixIm.GetOrigin())
print('')
print('MovIm origin =', MovIm.GetOrigin())

fix_dz = FixIm.GetSpacing()[2]
mov_dz = MovIm.GetSpacing()[2]
reg_dz = RegIm.GetSpacing()[2]

fix_z0 = FixIm.GetOrigin()[2]
mov_z0 = MovIm.GetOrigin()[2]
reg_z0 = RegIm.GetOrigin()[2]

fix_npa = sitk.GetArrayFromImage(FixIm)
mov_npa = sitk.GetArrayFromImage(MovIm)
reg_npa = sitk.GetArrayFromImage(RegIm)

fix_N = len(fix_npa)
mov_N = len(mov_npa)
reg_N = len(reg_npa)

# Compute the z position:
fix_z = [fix_z0 + fix_dz*(s) for s in range(fix_N)]
mov_z = [mov_z0 + mov_dz*(s) for s in range(mov_N)]
reg_z = [reg_z0 + reg_dz*(s) for s in range(reg_N)]

print('\nfix_z =\n', fix_z)
print('\nmov_z =\n', mov_z)
print('\nreg_z =\n', reg_z)


# In[42]:


FixSize = FixIm.GetSize()
MovSize = MovIm.GetSize()

print(f'z at FixIm[0,0,0]     = {FixIm.TransformIndexToPhysicalPoint((0, 0, 0))[2]} mm')
print(f'z at FixIm[96,128,0]  = {FixIm.TransformIndexToPhysicalPoint((int(FixSize[0]/2), int(FixSize[1]/2), 0))[2]} mm')
print(f'z at FixIm[192,256,0] = {FixIm.TransformIndexToPhysicalPoint((FixSize[0], FixSize[1], 0))[2]} mm')
print('')
print(f'z at MovIm[0,0,0]     = {MovIm.TransformIndexToPhysicalPoint((0, 0, 0))[2]} mm')
print(f'z at MovIm[96,128,0]  = {MovIm.TransformIndexToPhysicalPoint((int(MovSize[0]/2), int(MovSize[1]/2), 0))[2]} mm')
print(f'z at MovIm[192,256,0] = {MovIm.TransformIndexToPhysicalPoint((MovSize[0], MovSize[1], 0))[2]} mm')
print('')
print('FixIm z-Direction =', FixIm.GetDirection()[6:])
print('')
print('MovIm z-Direction =', MovIm.GetDirection()[6:])


# In[43]:


MovSize = MovIm.GetSize()

for z_ind in range(MovSize[2]):
    print(f'z at MovIm[96,128,{z_ind}] = {MovIm.TransformIndexToPhysicalPoint((int(MovSize[0]/2), int(MovSize[2]/2), z_ind))[2]} mm')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## Older stuff below...

# In[ ]:


# Compute the deformation field:
""" Running this in Anaconda Jupyter will probably kill the kernel.  If so, run in PyCharm. """
#DefField = ruf.ComputeDeformationField(MovingImage=MovingIm, TransformFilter=ElastixImFilt)

