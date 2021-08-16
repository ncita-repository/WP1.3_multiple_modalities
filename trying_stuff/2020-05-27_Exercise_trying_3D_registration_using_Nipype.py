#!/usr/bin/env python
# coding: utf-8

# # 3D Registration
# 
# Copied from previous notebook dated May 22nd but will try Nipype instead of SimpleElastix.

# ### Import packages and functions

# In[1]:


#import SimpleITK as sitk
#print(sitk.Version())
import nipype
print('nipype version', nipype.__version__)
from nipype.interfaces.elastix import ApplyWarp
from nipype.interfaces.elastix import PointsWarp
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


# ### 29/05:  Follow the PointsWarp example here:
# 
# https://nipype.readthedocs.io/en/1.1.5/interfaces/generated/interfaces.elastix/registration.html

# In[9]:


# Get current working directory:
CurrentWorkingDir = os.getcwd()

# Define the filepath for the output points:
OutputPtsFpath = os.path.join(CurrentWorkingDir, 'outputpoints_from_nipype.txt')

reg = PointsWarp()
reg.inputs.points_file = 'inputpoints.txt'
reg.inputs.transform_file = 'TransformParameters.0.txt'

reg.inputs.output_path = OutputPtsFpath
""" Results in error:
TraitError: The 'output_path' trait of a PointsWarpInputSpec instance must be a pathlike object 
or string representing an existing directory, but a value of 
'C:\\Code\\WP1.3_multiple_modalities\\trying_stuff\\outputpoints_from_nipype.txt' <class 'str'> 
was specified.
"""
#reg.inputs.output_path = 'outputpoints_from_nipype.txt'
""" Results in error:
TraitError: The 'output_path' trait of a PointsWarpInputSpec instance must be a pathlike object 
or string representing an existing directory, but a value of 'outputpoints_from_nipype.txt' 
<class 'str'> was specified.
"""

reg.cmdline


# ### There was no output.  I ran the same in PyCharm's Python console but also no output.
# 
# ### I've tried defining the output path but that has failed due to TraitErrors

# # End.
# 
# # Cells below were from the previous notebook file that was used as a template.

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

# In[19]:


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


# In[ ]:




