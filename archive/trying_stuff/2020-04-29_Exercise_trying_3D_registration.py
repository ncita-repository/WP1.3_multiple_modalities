#!/usr/bin/env python
# coding: utf-8

# # 3D Registration

# In[1]:


import SimpleITK as sitk
import time
import os
import copy
from GetDicomFpaths import GetDicomFpaths


# ### Working with TCIA ACRIN-FMISO Brain data:
# 
# ### Get data here:https://drive.google.com/drive/folders/1WhXPJ2LB9yuGqdPNFZkadTFEDXYgGZkZ?usp=sharing
# 
# ### The DICOMs are organised in separate folders for MR and CT imaging.  The MR dataset also contains an ROICollection (created in OHIF-Viewer) saved as a DICOM-RTSTRUCT file.
# 
# ### The intention was that MR would be the "fixed" set, and CT the "moving" set.  And once the CT images are registered to the MR, the contours from MR would be mapped to the registered CT volume.
# 
# ### Assuming that FixedIm = MR image (3D) and MovingIm = CT image (3D):
# 
# #### FixedIm Size  = (512, 512, 72)
# #### MovingIm Size = (256, 256, 160)
# 
# #### FixedIm Origin  = (-108.2880859375, -282.97186279296875, 115.51118469238281)
# #### MovingIm Origin = (-134.53753662109375, -167.07070922851562, -79.45004272460938)
# 
# #### FixedIm Spacing  = (0.423828125, 0.423828125, 2.4130020141601562)
# #### MovingIm Spacing = (1.0, 1.0, 1.0)

# ### Define the directories of Fixed and Moving:

# In[9]:


FixedDir = r'C:\Temp\ACRIN_MR_DICOMs'
MovingDir = r'C:\Temp\ACRIN_CT_DICOMs'


# ### First tried using sitk.JoinSeries(FixedVectorOfImages) and sitk.JoinSeries(MovingVectorOfImages) but not sure JoinSeries is appropriate to use here...

# In[11]:


# Use GetDicomFpaths to get list of Fixed and Moving filepaths.
FixedFpaths = GetDicomFpaths(FixedDir, 'slices', Debug=False)
MovingFpaths = GetDicomFpaths(MovingDir, 'slices', Debug=False)


# Start timing:
times = []
times.append(time.time())

# Concatenate the 2D images into one 3D image:
FixedVectorOfImages = sitk.VectorOfImage()
MovingVectorOfImages = sitk.VectorOfImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to concatenate the 2D images.')

for fpath in FixedFpaths:
    FixedVectorOfImages.push_back(sitk.ReadImage(fpath))
    
times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to read the Fixed 3D image stack.')
    
for fpath in MovingFpaths:
    MovingVectorOfImages.push_back(sitk.ReadImage(fpath))

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to read the Moving 3D image stack.')

FixedIm = sitk.JoinSeries(FixedVectorOfImages)
MovingIm = sitk.JoinSeries(MovingVectorOfImages)

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
ImFilter.SetFixedImage(FixedIm)
ImFilter.SetMovingImage(MovingIm)
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('groupwise'))

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Try registering a stack to itself to see if anything can be made to work:

# In[39]:



# Start timing:
times = []
times.append(time.time())

FixedIm = sitk.JoinSeries(FixedVectorOfImages)
#MovingIm = sitk.JoinSeries(MovingVectorOfImages)

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
ImFilter.SetFixedImage(FixedIm)
#ImFilter.SetMovingImage(MovingIm)
ImFilter.SetMovingImage(FixedIm)
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('groupwise'))

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### The error "Inputs do not occupy the same physical space" suggests to me that JoinSeries() is meant to be used to align a series of images of similar physical content - not for a stack of images that correspond to completely different images as is the case in a 3D scan.
# 
# ### This time use ImageSeriesReader():

# In[14]:


FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()

FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
print('\nFixedIm Size  =', FixedSize)
print('MovingIm Size =', MovingSize)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
print('\nFixedIm Origin  =', FixedOrigin)
print('MovingIm Origin =', MovingOrigin)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
print('\nFixedIm Spacing  =', FixedSpacing)
print('MovingIm Spacing =', MovingSpacing)


# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.SetFixedImage(FixedIm)
ImFilter.SetMovingImage(MovingIm)
#ImFilter.SetMovingImage(FixedIm)
#ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('groupwise'))
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
RegIm3D = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Try setting some additional distance metrics following the thread:
# 
# ### https://github.com/SuperElastix/SimpleElastix/issues/70

# In[41]:



# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()

# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
""" 
Setting some additional parameters following the thread:
https://github.com/SuperElastix/SimpleElastix/issues/70
"""
ImFilter.SetParameter( "Metric", "CorrespondingPointsEuclideanDistanceMetric" )
ImFilter.SetParameter( "Metric", ("NormalizedMutualInformation", "CorrespondingPointsEuclideanDistanceMetric",))
ImFilter.SetParameter("Metric0Weight", "0.0")
ImFilter.SetFixedImage(FixedIm)
ImFilter.SetMovingImage(MovingIm)
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Same as above but setting the fixed and moving images to be the same:

# In[16]:



# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
""" 
Setting some additional parameters following the thread:
https://github.com/SuperElastix/SimpleElastix/issues/70
"""
ImFilter.SetParameter( "Metric", "CorrespondingPointsEuclideanDistanceMetric" )
ImFilter.SetParameter( "Metric", ("NormalizedMutualInformation", "CorrespondingPointsEuclideanDistanceMetric",))
ImFilter.SetParameter("Metric0Weight", "0.0")
ImFilter.SetFixedImage(FixedIm)
#ImFilter.SetMovingImage(MovingIm)
ImFilter.SetMovingImage(FixedIm) # <-- SetMovingImage = FixedIm
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### So performing a self registration worked.
# 
# ### Next try setting AutomaticTransformInitializationMethod as in section 5.5.2 in the manual:

# In[42]:



# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()

# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
""" 
Setting some additional parameters from section 5.5.2 in the manual:
"""
ImFilter.SetParameter( "AutomaticTransformInitializationMethod", "GeometricalCenter" )
ImFilter.SetFixedImage(FixedIm)
ImFilter.SetMovingImage(MovingIm)
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Repeat above but perform registration of FixedIm3D to itself:

# In[21]:



# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
""" 
Setting some additional parameters from section 5.5.2 in the manual:
"""
ImFilter.SetParameter( "AutomaticTransformInitializationMethod", "GeometricalCenter" )
ImFilter.SetFixedImage(FixedIm)
#ImFilter.SetMovingImage(MovingIm)
ImFilter.SetMovingImage(FixedIm) # <-- SetMovingImage = FixedIm
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

RegSize = RegIm.GetSize()
print('\nRegIm Size  =', RegSize)

RegOrigin = RegIm.GetOrigin()
print('\nRegIm Origin  =', RegOrigin)

RegSpacing = RegIm.GetSpacing()
print('\nRegIm Spacing  =', RegSpacing)


# ### So performing a self-registration worked.

# ### Next try reversing the Fixed and Moving stack as in the following:
# 
# ### https://groups.google.com/forum/#!category-topic/elastix-imageregistration/simpleelastix/HPVGTR5V6cQ

# In[45]:


SwapFixedMoving = False
#SwapFixedMoving = True

if SwapFixedMoving:
    print('FixedDir and MovingDir swapped.')
else:
    print('FixedDir and MovingDir as originally defined.')

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

# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ImFilter = sitk.ElastixImageFilter()
ImFilter.LogToConsoleOn()
ImFilter.LogToFileOn()
""" 
Setting some additional parameters from section 5.5.2 in the manual:
"""
ImFilter.SetParameter( "AutomaticTransformInitializationMethod", "GeometricalCenter" )
ImFilter.SetFixedImage(FixedIm)
ImFilter.SetMovingImage(MovingIm)
ImFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

RegSize = RegIm.GetSize()
print('\nRegIm Size  =', RegSize)

RegSpacing = RegIm.GetSpacing()
print('\nRegIm Spacing  =', RegSpacing)

RegOrigin = RegIm.GetOrigin()
print('\nRegIm Origin  =', RegOrigin)


# ### That didn't help, even after casting to float32.
# 
# ### April 30:
# 
# ### Try following the example shown in this tutorial:
# 
# ### https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb

# In[24]:


# Define utility functions:

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from ipywidgets import interact, fixed
from IPython.display import clear_output

# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
def display_images(fixed_image_z, moving_image_z, fixed_npa, moving_npa):
    # create a figure with two subplots and the specified size
    plt.subplots(1,2,figsize=(10,8))
    
    # draw the fixed image in the first subplot
    plt.subplot(1,2,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title('fixed image')
    plt.axis('off')
    
    # draw the moving image in the second subplot
    plt.subplot(1,2,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title('moving image')
    plt.axis('off')

# callback invoked by the ipython interact method for scrolling and modifying the alpha blending
# of an image stack of two images that occupy the same physical space. 
def display_images_with_alpha(image_z, alpha, fixed, moving):
    img = (1.0 - alpha)*fixed[:,:,image_z] + alpha*moving[:,:,image_z] 
    plt.imshow(sitk.GetArrayFromImage(img),cmap=plt.cm.Greys_r);
    plt.axis('off')
    
    
# callback invoked when the StartEvent happens, sets up our new data
def start_plot():
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []

# callback invoked when the EndEvent happens, do cleanup of data and figure
def end_plot():
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    #close figure, we don't want to get a duplicate of the plot latter on
    plt.close()

# callback invoked when the IterationEvent happens, update our data and display new figure    
def plot_values(registration_method):
    global metric_values, multires_iterations
    
    metric_values.append(registration_method.GetMetricValue())                                       
    #clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    #plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
# callback invoked when the sitkMultiResolutionIterationEvent happens, update the index into the 
# metric_values list. 
def update_multires_iterations():
    global metric_values, multires_iterations
    multires_iterations.append(len(metric_values))


# In[25]:


# Read images:

#fixed_image =  sitk.ReadImage("Data/training_001_ct.mha",
#                              sitk.sitkFloat32)
#moving_image = sitk.ReadImage("Data/training_001_mr_T1.mha",
#                              sitk.sitkFloat32) 

# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingNames)
MovingIm = MovingReader.Execute()

interact(display_images,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)));


# In[26]:


# Initial alignment:

CenteringMode = sitk.CenteredTransformInitializerFilter.GEOMETRY

InitialTransform =     sitk.CenteredTransformInitializer(FixedIm, 
                                      MovingIm, 
                                      sitk.Euler3DTransform(), 
                                      CenteringMode)

MovingImResampled = sitk.Resample(MovingIm,
                                  FixedIm,
                                  InitialTransform,
                                  sitk.sitkLinear,
                                  0.0,
                                  MovingIm.GetPixelIDValue())

interact(display_images_with_alpha,
         image_z=(0,FixedIm.GetSize()[2]),
         alpha=(0.0,1.0,0.05),
         fixed = fixed(FixedIm),
         moving=fixed(MovingImResampled));


# ### Modify the code a bit based on zivy's suggested resampling code on 25 Sept 2018 under "def GlobalAffineRegistration(fixed,moving)":
# 
# ### https://github.com/SimpleITK/SimpleITK/issues/576

# In[47]:


# Initial alignment:

#CenteringMode = sitk.CenteredTransformInitializerFilter.GEOMETRY

#InitialTransform = \
#    sitk.CenteredTransformInitializer(FixedIm, 
#                                      MovingIm, 
#                                      sitk.Euler3DTransform(), 
#                                      CenteringMode)
#
#MovingImResampled = sitk.Resample(MovingIm,
#                                  FixedIm,
#                                  InitialTransform,
#                                  sitk.sitkLinear,
#                                  0.0,
#                                  MovingIm.GetPixelIDValue())

InitialTransform =     sitk.CenteredTransformInitializer(FixedIm, 
                                      MovingIm, 
                                      sitk.AffineTransform(FixedIm.GetDimension()))


Resampler = sitk.ImageRegistrationMethod()

Resampler.SetMetricAsJointHistogramMutualInformation(32)

Resampler.SetOptimizerAsGradientDescent(learningRate=1.0,
                                        numberOfIterations=100,
                                        estimateLearningRate = Resampler.EachIteration)

Resampler.SetOptimizerScalesFromPhysicalShift()
Resampler.SetInitialTransform(InitialTransform, inPlace=False)
Resampler.SetInterpolator(sitk.sitkLinear)

ResamplerTransform = Resampler.Execute(sitk.Cast(FixedIm, sitk.sitkFloat32), 
                                       sitk.Cast(MovingIm, sitk.sitkFloat32))

MovingImResampled = sitk.Resample(MovingIm, FixedIm, ResamplerTransform, 
                                  sitk.sitkLinear, 0.0, MovingIm.GetPixelID())
    
interact(display_images_with_alpha,
         image_z=(0,FixedIm.GetSize()[2]),
         alpha=(0.0,1.0,0.05),
         fixed = fixed(FixedIm),
         moving=fixed(MovingImResampled));


# ### No errors this time but not sure where to go from here.  No images displayed..
# 
# ### May 1:  Use as a template the following tutorial:
# 
# ### https://github.com/SimpleITK/SPIE2018_COURSE/blob/master/04_basic_registration.ipynb

# In[48]:


import SimpleITK as sitk
#from downloaddata import fetch_data as fdata
get_ipython().run_line_magic('matplotlib', 'notebook')
import gui
import registration_gui as rgui

import numpy as np
import os
import time
import copy
from CreateDir import CreateDir

# Define the directory to export results to:
#TodaysDate = time.strftime("%Y-%m-%d", time.gmtime())
#RootDir = 'C:\Temp'
#ExportDir = os.path.join(r'C:\Temp', TodaysDate)
#CreateDir(ExportDir, Debug=False)

print('FixedDir =', FixedDir)
print('MovingDir =', MovingDir)

# Read in the 3D stack of Fixed DICOMs:
FixedReader = sitk.ImageSeriesReader()
FixedFilenames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedFilenames)
FixedIm = FixedReader.Execute()

# Read in the 3D stack of Moving DICOMs:
MovingReader = sitk.ImageSeriesReader()
MovingFilenames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
MovingReader.SetFileNames(MovingFilenames)
MovingIm = MovingReader.Execute()

# Print original pixel data type:
FixedDtype = FixedIm.GetPixelIDTypeAsString()
MovingDtype = MovingIm.GetPixelIDTypeAsString()
print('\nOriginal FixedIm Dtype  =', FixedDtype)
print('Original MovingIm Dtype =', MovingDtype)

# Convert from int16/uint16 to float32:
FixedIm = sitk.Cast(FixedIm, sitk.sitkFloat32)
MovingIm = sitk.Cast(MovingIm, sitk.sitkFloat32)

# Print original pixel data type after casting to float32:
FixedDtype = FixedIm.GetPixelIDTypeAsString()
MovingDtype = MovingIm.GetPixelIDTypeAsString()
print('\nFixedIm3D Dtype after Cast()  =', FixedDtype)
print('MovingIm3D Dtype after Cast() =', MovingDtype)

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

gui.MultiImageDisplay(image_list = [FixedIm, MovingIm],                   
                      title_list = ['Fixed', 'Moving'], figure_size=(8,4));


# ### Perform "classical" registration:

# In[33]:


# Initialize registration by aligning the centers of the two volumes. 
# To qualitatively evaluate the result we use a linked cursor approach, 
# click on one image and the corresponding point is added to the other image.

InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# ### Poor result.
# 
# ### Try again but with reversal of Fixed/Moving:

# In[35]:


# Initialize registration by aligning the centers of the two volumes. 
# To qualitatively evaluate the result we use a linked cursor approach, 
# click on one image and the corresponding point is added to the other image.

InitialTransform = sitk.CenteredTransformInitializer(MovingIm, 
                                                     FixedIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)

gui.RegistrationPointDataAquisition(MovingIm, FixedIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# ### Swapping FixedIm/MovingIm made no difference.
# 
# ### Next try setting the similarity metric and optimiser settings:

# In[37]:


# First re-define InitialTransform as before swapping Fixed with Moving:
# Initialize registration by aligning the centers of the two volumes. 
# To qualitatively evaluate the result we use a linked cursor approach, 
# click on one image and the corresponding point is added to the other image.

InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)

# Set the similarity metric and optimiser settings:
RegMethod = sitk.ImageRegistrationMethod()

# Similarity metric settings:
RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
RegMethod.SetMetricSamplingPercentage(0.01)

# Set the interpolator:
RegMethod.SetInterpolator(sitk.sitkLinear)

# Optimiser settings:
RegMethod.SetOptimizerAsGradientDescent(learningRate=1.0, 
                                        numberOfIterations=100, 
                                        convergenceMinimumValue=1e-6, 
                                        convergenceWindowSize=10)

RegMethod.SetOptimizerScalesFromPhysicalShift()

# Setup for the multi-resolution framework:       
RegMethod.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

# Don't optimize in-place, we would possibly like to run this cell multiple times:
RegMethod.SetInitialTransform(InitialTransform, inPlace=False) 
""" Note: InitialTransform = sitk.CenteredTransformInitializerFilter.GEOMETRY """

# Connect all of the observers so that we can perform plotting during registration:
RegMethod.AddCommand(sitk.sitkStartEvent, rgui.start_plot)
RegMethod.AddCommand(sitk.sitkEndEvent, rgui.end_plot)
RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, rgui.update_multires_iterations) 
RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: rgui.plot_values(RegMethod))

FinalTransformClassical = RegMethod.Execute(FixedIm, MovingIm)

# Check the reason optimisation terminated:
print(f'Final metric value: {RegMethod.GetMetricValue()}')
print(f'Optimizer\'s stopping condition: {RegMethod.GetOptimizerStopConditionDescription()}')


# In[49]:


# Qualitatively evaluate the result using a linked cursor approach (visual evaluation):
gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=FinalTransformClassical);


# ### The points seem to correspond reasonably well, but the moving image doesn't look registered.
# ### 05/05 Correction:  After re-evaluating cells the results now look as poor as the ones above.

# In[ ]:


# Save the results to file:
#MovingResampled = sitk.Resample(MovingIm, FixedIm, FinalTransform, 
#                                sitk.sitkLinear, 0.0, MovingIm.GetPixelID())

#sitk.WriteImage(MovingResampled, os.path.join(ExportDir, 'T2_3D_Slab_resampled.mha'))
#sitk.WriteTransform(FinalTransform, os.path.join(ExportDir, 'T2_3D_Slab_reg_to_HEAD_2-4_H40_SOFT_TISSUE.tfm'))


# ### Use ITK v4 registration framework:

# In[50]:


# Define the registration method:
RegMethod = sitk.ImageRegistrationMethod()

# Similarity metric settings:
RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
RegMethod.SetMetricSamplingPercentage(0.01)

# Set the interpolator:
RegMethod.SetInterpolator(sitk.sitkLinear)

# Optimizer settings:
RegMethod.SetOptimizerAsGradientDescent(learningRate=1.0, 
                                        numberOfIterations=100, 
                                        convergenceMinimumValue=1e-6, 
                                        convergenceWindowSize=10)

RegMethod.SetOptimizerScalesFromPhysicalShift()

# Setup for the multi-resolution framework:         
RegMethod.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

""" Note the above is the same as the classical registration.  They depart from here: """

# Set the initial moving and optimised transforms:
OptimizedTransform = sitk.Euler3DTransform()    
RegMethod.SetMovingInitialTransform(InitialTransform) 
""" Note: InitialTransform = sitk.CenteredTransformInitializerFilter.GEOMETRY """
RegMethod.SetInitialTransform(OptimizedTransform, inPlace=False)

# Connect all of the observers so that we can perform plotting during registration:
RegMethod.AddCommand(sitk.sitkStartEvent, rgui.start_plot)
RegMethod.AddCommand(sitk.sitkEndEvent, rgui.end_plot)
RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, rgui.update_multires_iterations) 
RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: rgui.plot_values(RegMethod))

# Need to compose the transformations after registration.
FinalTransform_v4 = RegMethod.Execute(FixedIm, MovingIm)
#FinalTransform_v4.AddTransform(InitialTransform)

# Check the reason optimisation terminated:
print(f'Final metric value: {RegMethod.GetMetricValue()}')
print(f'Optimizer\'s stopping condition: {RegMethod.GetOptimizerStopConditionDescription()}')


# In[51]:


# Qualitatively evaluate the result using a linked cursor approach (visual evaluation):
gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=FinalTransform_v4);


# ### For some reason the gui doesn't seem to be showing the corresponding moving slice based on the cursor as it had for the classical registration...

# ### This time try using the identify transform as the initialiser:

# In[52]:


# Use the identity transform as the initial transform:
InitialTransform = sitk.Transform()

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# In[53]:


# Define the registration method:
RegMethod = sitk.ImageRegistrationMethod()

# Similarity metric settings:
RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
RegMethod.SetMetricSamplingPercentage(0.01)

# Set the interpolator:
RegMethod.SetInterpolator(sitk.sitkLinear)

# Optimizer settings:
RegMethod.SetOptimizerAsGradientDescent(learningRate=1.0, 
                                        numberOfIterations=100, 
                                        convergenceMinimumValue=1e-6, 
                                        convergenceWindowSize=10)

RegMethod.SetOptimizerScalesFromPhysicalShift()

# Setup for the multi-resolution framework:         
RegMethod.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

""" Note the above is the same as the classical registration.  They depart from here: """

# Set the initial moving and optimised transforms:
OptimizedTransform = sitk.Euler3DTransform()    
RegMethod.SetMovingInitialTransform(InitialTransform) 
""" Note: InitialTransform = sitk.Transform() i.e. identity transform """
RegMethod.SetInitialTransform(OptimizedTransform, inPlace=False)

# Connect all of the observers so that we can perform plotting during registration:
RegMethod.AddCommand(sitk.sitkStartEvent, rgui.start_plot)
RegMethod.AddCommand(sitk.sitkEndEvent, rgui.end_plot)
RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, rgui.update_multires_iterations) 
RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: rgui.plot_values(RegMethod))

# Need to compose the transformations after registration.
FinalTransform_v4 = RegMethod.Execute(FixedIm, MovingIm)
#FinalTransform_v4.AddTransform(InitialTransform)

# Check the reason optimisation terminated:
print(f'Final metric value: {RegMethod.GetMetricValue()}')
print(f'Optimizer\'s stopping condition: {RegMethod.GetOptimizerStopConditionDescription()}')


# In[126]:


# Get the patient positions:
img = sitk.ReadImage(FixedFilenames[0])
print('Fixed Patient position:  ' + img.GetMetaData('0018|5100'))
img = sitk.ReadImage(MovingFilenames[0])
print('Moving Patient position: ' + img.GetMetaData('0018|5100'))


# ### Try using CenteredTransformInitializer initialization:
# 
# Compare GEOMETRY and MOMENTS based approaches.

# In[55]:


# First GEOMETRY approach:
InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# In[56]:


# Next MOMENTS approach:
InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.MOMENTS)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# ### Both look highly inaccurate to me!

# ### Try using the Exhaustive optimizer initialization using CenteredTransformInitializer with GEOMETRY approach:
# 
# Compare GEOMETRY and MOMENTS approaches.

# In[57]:


InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)
RegMethod = sitk.ImageRegistrationMethod()
RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
RegMethod.SetMetricSamplingPercentage(0.01)
RegMethod.SetInterpolator(sitk.sitkLinear)

# The order of parameters for the Euler3DTransform is [angle_x, angle_y, angle_z, t_x, t_y, t_z]. The parameter 
# sampling grid is centered on the InitialTransform parameter values, that are all zero for the rotations. Given
# the number of steps, their length and optimizer scales we have:
# angle_x = 0
# angle_y = -pi, 0, pi
# angle_z = -pi, 0, pi
RegMethod.SetOptimizerAsExhaustive(numberOfSteps=[0,1,1,0,0,0], stepLength = np.pi)
RegMethod.SetOptimizerScales([1,1,1,1,1,1])

# Perform the registration in-place so that the InitialTransform is modified.
RegMethod.SetInitialTransform(InitialTransform, inPlace=True)
""" Note: InitialTransform = sitk.CenteredTransformInitializerFilter.GEOMETRY """
RegMethod.Execute(FixedIm, MovingIm)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4));


# ### Try MOMENTS this time:

# In[58]:


InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.MOMENTS)
RegMethod = sitk.ImageRegistrationMethod()
RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
RegMethod.SetMetricSamplingPercentage(0.01)
RegMethod.SetInterpolator(sitk.sitkLinear)

# The order of parameters for the Euler3DTransform is [angle_x, angle_y, angle_z, t_x, t_y, t_z]. The parameter 
# sampling grid is centered on the InitialTransform parameter values, that are all zero for the rotations. Given
# the number of steps, their length and optimizer scales we have:
# angle_x = 0
# angle_y = -pi, 0, pi
# angle_z = -pi, 0, pi
RegMethod.SetOptimizerAsExhaustive(numberOfSteps=[0,1,1,0,0,0], stepLength = np.pi)
RegMethod.SetOptimizerScales([1,1,1,1,1,1])

# Perform the registration in-place so that the InitialTransform is modified.
RegMethod.SetInitialTransform(InitialTransform, inPlace=True)
""" Note: InitialTransform = sitk.CenteredTransformInitializerFilter.MOMENTS """
RegMethod.Execute(FixedIm, MovingIm)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4));


# ### For some reason the interactive cursor feature is not working, so difficult to assess whether the results of the exhaustive optimiser are any better.
# 
# ### This is the end of the example covering rigid transformations.

# In[ ]:





# In[ ]:




