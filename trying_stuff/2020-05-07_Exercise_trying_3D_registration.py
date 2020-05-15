#!/usr/bin/env python
# coding: utf-8

# # 3D Registration

# In[4]:


import SimpleITK as sitk
#print(sitk.Version())
import numpy as np
import time
import os
import copy
import importlib

from GetDicomFpaths import GetDicomFpaths

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from ipywidgets import interact, fixed
from IPython.display import clear_output

get_ipython().run_line_magic('matplotlib', 'notebook')
import gui
import registration_gui as rgui


# ### Define utility functions:

# In[220]:


# Define utility functions:
""" Source:
https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb
"""

# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Original: """
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


# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Slight modification to display_images: Add figure titles as arguments: """
def display_images_v2(fixed_image_z, moving_image_z, fixed_npa, moving_npa, fixed_title, moving_title):
    # Get min, max:
    #fixed_min = np.min(fixed_npa)
    #fixed_max = np.max(fixed_npa)
    #moving_min = np.min(moving_npa)
    #moving_max = np.max(moving_npa)

    # create a figure with two subplots and the specified size
    plt.subplots(1,2,figsize=(10,8))

    # draw the fixed image in the first subplot
    plt.subplot(1,2,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, vmin=fixed_min, vmax=fixed_max);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(fixed_npa[fixed_image_z,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(fixed_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,2,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, vmin=moving_min, vmax=moving_max);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(moving_npa[moving_image_z,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(moving_title)
    plt.axis('off')

# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Slight modification to display_images_v2: Plot registered image as well: """
def display_images_v3(fixed_image_z, moving_image_z, registered_image_z, 
                      fixed_npa, moving_npa, registered_npa, 
                      fixed_title, moving_title, registered_title):

    # create a figure with two subplots and the specified size
    plt.subplots(1,3,figsize=(9,6))

    # draw the fixed image in the first subplot
    plt.subplot(1,3,1)
    plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(fixed_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,3,2)
    plt.imshow(moving_npa[moving_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(moving_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1,3,3)
    plt.imshow(registered_npa[registered_image_z,:,:],cmap=plt.cm.Greys_r);
    plt.title(registered_title)
    plt.axis('off')
    
# callback invoked by the interact ipython method for scrolling through the image stacks of
# the two images (moving and fixed)
""" Add contours to display_images_v3 """
def display_images_and_contours(fixed_ind, moving_ind, reg_ind, 
                                fixed_npa, moving_npa, reg_npa, 
                                fixed_im, moving_im, reg_im,
                                fixed_contour_pts,
                                fixed_title, moving_title, reg_title):
    
    # Get the number of slices in the images:
    fixed_Nslices = fixed_im.GetSize()[2]
    moving_Nslices = moving_im.GetSize()[2]
    reg_Nslices = reg_im.GetSize()[2]

    # Get the z positions of the origins:
    fixed_origin_z = fixed_im.GetOrigin()[2]
    moving_origin_z = moving_im.GetOrigin()[2]
    reg_origin_z = reg_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fixed_dz = fixed_im.GetSpacing()[2]
    moving_dz = moving_im.GetSpacing()[2]
    reg_dz = reg_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fixed_zpos = fixed_origin_z + fixed_ind*fixed_dz
    moving_zpos = moving_origin_z + moving_ind*moving_dz
    reg_zpos = reg_origin_z + reg_ind*reg_dz

    # Get the number of contours for this index:
    fixed_Ncontours = len(fixed_contour_pts[fixed_ind])

    print(f'\nFixed slice  {fixed_ind}/{fixed_Nslices}',
          f'at z = {round(fixed_zpos, 2)} mm')
    print(f'Moving slice {moving_ind}/{moving_Nslices}',
          f'at z = {round(moving_zpos, 2)} mm')
    print(f'Reg slice    {reg_ind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fixed_zpos, 2)} mm')
    #print(f'Moving depth position {round(moving_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    print(f'\nFixed no. of contours  = {fixed_Ncontours}')
    #print(f'Moving no. of contours = {fixed_Ncontours}')
    #print(f'Reg no. of contours    = {fixed_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,5))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(fixed_npa[fixed_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fixed_npa[fixed_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fixed_title)
    plt.axis('off')
    
    # Plot contours for fixed image if Npts is not empty:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[fixed_ind][0]:
        for x, y, z in fixed_contour_pts[fixed_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(fixed_npa[fixed_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(moving_npa[moving_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(moving_npa[moving_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(moving_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(reg_npa[reg_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    # Using reg_inds and fixed_contour_pts, plot contours belonging to the 
    # fixed image for the registered image if Npts is not empty:
    if fixed_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fixed_contour_pts[reg_ind][0]:
        for x, y, z in fixed_contour_pts[reg_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(reg_npa[reg_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

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
    

# Plot registration results as was done for DICOM objects (i.e. PlotRegResults()) but for plotting
# of SimpleITK objects:
#def PlotRegResults(FixedImage, MovingImage, RegisteredImage):
#    # Import packages:
#    import SimpleITK as sitk
#    
#    # Convert from ITK to Numpy data arrays:
#    FixedImage = sitk.GetArrayFromImage(FixedImage)
#    MovingImage = sitk.GetArrayFromImage(MovingImage)
#    RegisteredImage = sitk.GetArrayFromImage(RegisteredImage)
    
    


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

# ### Define the directories of Fixed and Moving and load the images:

# In[5]:


FixedDir = r'C:\Temp\ACRIN_MR_DICOMs'
MovingDir = r'C:\Temp\ACRIN_CT_DICOMs'

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

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
print(f'\nFixedIm Min, Max  = {np.min(FixedNda)}, {np.max(FixedNda)}')
print(f'MovingIm Min, Max = {np.min(MovingNda)}, {np.max(MovingNda)}')


# ### First tried using sitk.JoinSeries(FixedVectorOfImages) and sitk.JoinSeries(MovingVectorOfImages) but not sure JoinSeries is appropriate to use here...

# In[3]:


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

# In[5]:


SwapFixedMoving = False
#SwapFixedMoving = True

FixedReader = sitk.ImageSeriesReader()
if SwapFixedMoving:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(MovingDir)
else:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

MovingReader = sitk.ImageSeriesReader()
if SwapFixedMoving:
    MovingNames = MovingReader.GetGDCMSeriesFileNames(FixedDir)
else:
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
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### See elastix_20200506_0843.log for more details:
# 
# itk::ExceptionObject (00000099665E7978)
# 
# Location: "ElastixTemplate - Run()" 
# 
# File: C:\SimpleElastix\build_20200108\Elastix\Common\CostFunctions\itkAdvancedImageToImageMetric.hxx
# 
# Line: 1005
# 
# Description: itk::ERROR: AdvancedMattesMutualInformationMetric(000002670EFF89F0): Too many samples map outside moving image buffer: 0 / 4332
# 
# Error occurred during actual registration.

# ### See comment by salan668 here:
#     
# https://github.com/SuperElastix/SimpleElastix/issues/132
# 
# "Here the fixed image was a 3D high-resolution image and moving image was a 3D low-resolution image. I found that if the fixed image and the moving image has similar resolution, SimpleElastix works. But if resolutions are different from each other, this error appears."

# ### Try reversing FixedIm and MovingIm:

# In[3]:


#SwapFixedMoving = False
SwapFixedMoving = True

FixedReader = sitk.ImageSeriesReader()
if SwapFixedMoving:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(MovingDir)
else:
    FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
FixedReader.SetFileNames(FixedNames)
FixedIm = FixedReader.Execute()

MovingReader = sitk.ImageSeriesReader()
if SwapFixedMoving:
    MovingNames = MovingReader.GetGDCMSeriesFileNames(FixedDir)
else:
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


# ### Same error as before - see elastix_20200506_0955.log for more details:
# 
# itk::ExceptionObject (0000000F06BE75E8)
# 
# Location: "ElastixTemplate - Run()" 
# 
# File: C:\SimpleElastix\build_20200108\Elastix\Common\CostFunctions\itkAdvancedImageToImageMetric.hxx
# 
# Line: 1005
# 
# Description: itk::ERROR: AdvancedMattesMutualInformationMetric(0000018241DBEE20): Too many samples map outside moving image buffer: 0 / 4805
# 
# Error occurred during actual registration.

# ### Try setting some additional distance metrics following the thread:
# 
# ### https://github.com/SuperElastix/SimpleElastix/issues/70
# 
# ### 06/05/2020:  Although part of the error code is the same ("Location unknown") I think this error might be unrelated to the error I got.

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


# ### So performing a self registration worked.  This is consistent with salan668's comment.
# 
# ### Next try setting AutomaticTransformInitializationMethod as in section 5.5.2 in the manual:
# 
# ### Also see comment by kaspermasrstal on 01/22/2018 here:
# ### https://github.com/SuperElastix/SimpleElastix/issues/152
# 
# #### "No, this is because the images do not overlap in the world coordinate system. You can try to enable automatric transform initialization by setting "AutomaticTransformInitialization" to "true". "

# In[3]:



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


# ### See error log elastix_20200507_0956.log for details:
# 
# #### itk::ExceptionObject (00000084C4BE7398)
# #### Location: "ElastixTemplate - Run()" 
# #### File: C:\SimpleElastix\build_20200108\Elastix\Common\CostFunctions\itkAdvancedImageToImageMetric.hxx
# #### Line: 1005
# #### Description: itk::ERROR: AdvancedMattesMutualInformationMetric(000001CD3C8E1BD0): Too many samples map outside moving image buffer: 0 / 4332
# 
# #### Error occurred during actual registration.

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


# ### That didn't help, even after casting to float32.  It seems that without knowing how to resolve the issue that SimpleElastix has with registrations of unequal resolutions it won't be possible to progress.  So will move onto SimpleITK framework instead.
# 
# ### April 30:
# 
# ### Try following the example shown in this tutorial:
# 
# ### https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb

# In[ ]:


# Define utility functions:

# These were moved to the top.
"""
import matplotlib.pyplot as plt
%matplotlib inline

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
"""


# In[ ]:


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


# In[13]:


fixed(sitk.GetArrayFromImage(FixedIm))


# In[16]:


# fixed_image_z=(0,FixedIm.GetSize()[2]-1)
# = (0, 159)

# fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm))

# plt.imshow(fixed_npa[fixed_image_z,:,:],cmap=plt.cm.Greys_r);



# fixed_npa[fixed_image_z,:,:]
# =
# fixed(sitk.GetArrayFromImage(FixedIm))[(0,FixedIm.GetSize()[2]-1),:,:]
sitk.GetArrayFromImage(FixedIm)[(0,FixedIm.GetSize()[2]-1),:,:]

# = sitk.GetArrayFromImage(FixedIm)[(0,FixedIm.GetSize()[2]-1),:,:]


# In[24]:


# Initial alignment:

InitialTransform =     sitk.CenteredTransformInitializer(FixedIm, 
                                      MovingIm, 
                                      sitk.Euler3DTransform(), 
                                      sitk.CenteredTransformInitializerFilter.GEOMETRY)

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


# In[ ]:


print(InitialTransform)


# ### Modify the code a bit based on zivy's suggested resampling code on 25 Sept 2018 under "def GlobalAffineRegistration(fixed,moving)":
# 
# ### https://github.com/SimpleITK/SimpleITK/issues/576

# In[21]:


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
# ### 06/05/20:  Re-executed today and get the same error as from the cell above.
# 
# ### May 1:  Use as a template the following tutorial:
# 
# ### https://github.com/SimpleITK/SPIE2018_COURSE/blob/master/04_basic_registration.ipynb

# In[ ]:


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
print('\nFixedIm Dtype after Cast()  =', FixedDtype)
print('MovingIm Dtype after Cast() =', MovingDtype)

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

# In[ ]:


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

# In[ ]:


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


# In[ ]:


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


# In[ ]:


# Qualitatively evaluate the result using a linked cursor approach (visual evaluation):
gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=FinalTransform_v4);


# ### For some reason the gui doesn't seem to be showing the corresponding moving slice based on the cursor as it had for the classical registration...

# ### This time try using the identify transform as the initialiser:

# In[ ]:


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

# In[ ]:


# First GEOMETRY approach:
InitialTransform = sitk.CenteredTransformInitializer(FixedIm, 
                                                     MovingIm, 
                                                     sitk.Euler3DTransform(), 
                                                     sitk.CenteredTransformInitializerFilter.GEOMETRY)

gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=InitialTransform);


# In[ ]:


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

# In[ ]:


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

# In[ ]:


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

# ### 06/05/2020:  Try another example:
# 
# ### https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethod1_docs.html

# In[32]:


import sys

sys.argv


# In[15]:


def command_iteration(method) :
    print("{0:3} = {1:10.5f} : {2}".format(method.GetOptimizerIteration(),
                                   method.GetMetricValue(),
                                   method.GetOptimizerPosition()))

#if len ( sys.argv ) < 4:
#    print( "Usage: {0} <fixedImageFilter> <movingImageFile> <outputTransformFile>".format(sys.argv[0]))
#    sys.exit ( 1 )


#fixed = sitk.ReadImage(sys.argv[1], sitk.sitkFloat32)

#moving = sitk.ReadImage(sys.argv[2], sitk.sitkFloat32)

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

# Registration:

R = sitk.ImageRegistrationMethod()
R.SetMetricAsMeanSquares()
R.SetOptimizerAsRegularStepGradientDescent(4.0, .01, 200 )
R.SetInitialTransform(sitk.TranslationTransform(FixedIm.GetDimension()))
R.SetInterpolator(sitk.sitkLinear)

R.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(R) )

outTx = R.Execute(FixedIm, MovingIm)

print("\n-------")
print(outTx)
print("Optimizer stop condition: {0}".format(R.GetOptimizerStopConditionDescription()))
print(" Iteration: {0}".format(R.GetOptimizerIteration()))
print(" Metric value: {0}".format(R.GetMetricValue()))

#sitk.WriteTransform(outTx,  sys.argv[3])

if ( not "SITK_NOSHOW" in os.environ ):

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(FixedIm);
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(100)
    resampler.SetTransform(outTx)

    MovingResampledIm = resampler.Execute(MovingIm)
    simg1 = sitk.Cast(sitk.RescaleIntensity(FixedIm), sitk.sitkUInt8)
    simg2 = sitk.Cast(sitk.RescaleIntensity(MovingResampledIm), sitk.sitkUInt8)
    cimg = sitk.Compose(simg1, simg2, simg1//2.+simg2//2.)
    if False:
        sitk.Show( cimg, "ImageRegistration1 Composition" )
        """ sitk::ERROR: No appropriate executable found. 
        So rather than using procedural interface use the object oriented interface as
        shown here:
        https://simpleitk.readthedocs.io/en/master/link_ImageViewing_docs.html
        --> Now it works after having downloaded Fiji but the kernel gets stuck.
        """
    # Use the default image viewer:
    if False:
        ImViewer = sitk.ImageViewer()
        ImViewer.SetTitle('ImageRegistration1 Composition')
        ImViewer.Execute(cimg)
        """ Same as above.
        """
    # Change viewer, and display again.
    #ImViewer.SetApplication('/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP')
    #ImViewer.Execute(grid_image)

    # Change the viewer command, (use ITK-SNAP's -z option to open the image in zoomed mode)
    #ImViewer.SetCommand('/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP -z 2')
    #ImViewer.Execute(grid_image)


# Get info about FixedIm, MovingIm and MovingResampledIm:

#FixedDtype = str(FixedIm.dtype)
FixedDtype = FixedIm.GetPixelIDTypeAsString()
#MovingDtype = str(MovingIm.dtype)
MovingDtype = MovingIm.GetPixelIDTypeAsString()
MovingResampledDtype = MovingResampledIm.GetPixelIDTypeAsString()
print('\nFixedIm Dtype           =', FixedDtype)
print('MovingIm Dtype          =', MovingDtype)
print('MovingResampledIm Dtype =', MovingResampledDtype)

FixedSize = FixedIm.GetSize()
MovingSize = MovingIm.GetSize()
MovingResampledSize = MovingResampledIm.GetSize()
print('\nFixedIm Size           =', FixedSize)
print('MovingIm Size          =', MovingSize)
print('MovingResampledIm Size =', MovingResampledSize)

FixedSpacing = FixedIm.GetSpacing()
MovingSpacing = MovingIm.GetSpacing()
MovingResampledSpacing = MovingResampledIm.GetSpacing()
print('\nFixedIm Spacing           =', FixedSpacing)
print('MovingIm Spacing          =', MovingSpacing)
print('MovingResampledIm Spacing =', MovingResampledSpacing)

FixedOrigin = FixedIm.GetOrigin()
MovingOrigin = MovingIm.GetOrigin()
MovingResampledOrigin = MovingResampledIm.GetOrigin()
print('\nFixedIm Origin           =', FixedOrigin)
print('MovingIm Origin          =', MovingOrigin)
print('MovingResampledIm Origin =', MovingResampledOrigin)

FixedNda = sitk.GetArrayFromImage(FixedIm)
MovingNda = sitk.GetArrayFromImage(MovingIm)
MovingResampledNda = sitk.GetArrayFromImage(MovingResampledIm)
print(f'\nFixedIm Min, Max           = {np.min(FixedNda)}, {np.max(FixedNda)}')
print(f'MovingIm Min, Max          = {np.min(MovingNda)}, {np.max(MovingNda)}')
print(f'MovingResampledIm Min, Max = {np.min(MovingResampledNda)}, {np.max(MovingResampledNda)}')

if np.min(MovingResampledNda) == np.max(MovingResampledNda):
    print('\nMin/Max values of MovingResampledIm are the same!!!')


# 
# ### After failing to execute the Image viewing example here:
# 
# ### https://simpleitk.readthedocs.io/en/master/link_ImageViewing_docs.html
# 
# 
# ### I left a question here:
# 
# ### https://github.com/SimpleITK/SimpleITK/issues/723

# ### Check results of resampling MovingIm:

# In[7]:


# First plot FixedIm and MovingIm using display_images:

interact(display_images,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         fixed_title = 'FixedIm',
         moving_title = 'MovingIm');
    


# In[10]:


# First plot FixedIm and MovingIm using display_images_v2:

interact(display_images_v2,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         fixed_title = 'FixedIm',
         moving_title = 'MovingIm');
    


# In[13]:


# Next plot FixedIm and MovingResampledIm using display_images:

interact(display_images,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingResampledIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingResampledIm)),
         fixed_title = 'FixedIm',
         moving_title = 'MovingResampledIm');


# In[14]:


# Next plot FixedIm and MovingResampledIm using display_images_v2:

interact(display_images_v2,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingResampledIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingResampledIm)),
         fixed_title = 'FixedIm',
         moving_title = 'MovingResampledIm');


# ### I've finally realised that the min/max values of MovingResampledIm are equal - and that's why the plots of MovingResampledIm are all black - not due to issues with normalisation in the colourmap!

# In[16]:


MovingResampledNda


# ### So now need to figure out why the values are the same.
# 
# ### On page 13 here:
# 
# ### https://simpleitk.org/SPIE2018_COURSE/images_and_resampling.pdf
# 
# ### under "Common errors" there's talk of an empty (all black) image after resampling due to wrong settings for the resampling grid, or from using the inverse of the transformation.
# 
# ### I'm not sure this applies to the above because I'm not getting zero values.  But maybe it may still be the same root cause.  
# 
# ### I don't know what resampling grid was used above.  Nothing stands out in the code as having anything to do with a resampling grid.  The pdf linked above also doesn't explain where the sampling grid is within or how it's defined.  Page 12 lists 3 methods of specifying the resampling grid that are all incomprehensible.
# 

# # 07.05/2020:  Try using data that might be easier to work with.  MR and CT/PET data have completely different physical locations and different resolutions and pixel spacings.
# 
# ### Next try working with CT and PET (different resolutions and pixel spacings but same physical coordinates). 
# 
# ### If that still fails, then switch to two MR series (same resolution, same pixel spacing, same physical coordinates).

# ### Define the directories of Fixed and Moving and load the images:

# In[6]:


#FixedDir = r'C:\Temp\ACRIN_MR_DICOMs'
FixedDir = r'C:\Temp\ACRIN_PET_DICOMs'
MovingDir = r'C:\Temp\ACRIN_CT_DICOMs'

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


# ### Looks like the PET/CT data do not have the same physical coordinates.  
# 
# # Work with MR only data sets.

# ### Define the directories of Fixed and Moving and load the images:

# In[10]:



RootDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixedDir = os.path.join(RootDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovingDir = os.path.join(RootDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

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


# In[ ]:


gui.MultiImageDisplay(image_list = [FixedIm, MovingIm],                   
                      title_list = ['Fixed', 'Moving'], figure_size=(8,4));


# In[35]:


interact(display_images_v2,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image');


# ### Trying using SimpleElastix framework again:

# In[12]:


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
RegIm = ImFilter.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# In[21]:


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
print('RegIm Spacing bb =', RegSpacing)

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
TransformParameterMap = ImFilter.GetTransformParameterMap()
TransformParameterMap


# In[23]:


# Try to get the transformation (transformix) as here:
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

TransformixImFilter = sitk.TransformixImageFilter()
TransformixImFilter


# In[ ]:


gui.MultiImageDisplay(image_list = [MovingIm, RegIm],                   
                      title_list = ['Moving', 'Registered'], figure_size=(8,4));


# In[ ]:


# Qualitatively evaluate the result using a linked cursor approach (visual evaluation):
gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=ImFilter);


# In[ ]:


# Qualitatively evaluate the result using a linked cursor approach (visual evaluation):
gui.RegistrationPointDataAquisition(FixedIm, MovingIm, figure_size=(8,4), 
                                    known_transformation=TransformixImFilter);


# In[34]:


interact(display_images_v3,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         registered_image_z=(0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         registered_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         registered_title = 'Registered image');


# # 11/05/2020:  Repeat above but this time using data downloaded from XNAT with view to start working with ROI Collection:

# ### Download XNAT data (DICOMs, ROI Collection):

# In[6]:


import importlib

import GetXnatData
importlib.reload(GetXnatData)
from GetXnatData import GetXnatData

import CopyRoi
importlib.reload(CopyRoi)
from CopyRoi import CopyRoi


""" 
Make selections below by un-commenting the desired option.
"""

# Define the XNAT address, username and password:
XnatAddress = 'http://10.1.1.17'
XnatUsername = 'admin'
XnatPassword = 'admin'

# Decide on the subject:
Subject = 'ACRIN-FMISO-Brain-011'
#Subject = ''

# Decide whether to copy to select PET or CT scans:
#CopyFrom = 'MR'
#CopyTo = 'PET'
#CopyTo = 'CT'
CopyFrom = 'MR_4'# 1960-04-10 69626
CopyTo = 'MR_12' # 1961-06-11 79433

# Decide on what feature to use to manually determine the required shift 
# along the scan direction:
#ShiftFeature = 'peak of nose'
#ShiftFeature = 'infinity-like structure in brain'
#ShiftFeature = 'lens of eye'
#ShiftFeature = 'butterfly-like structure in brain'
ShiftFeature = '' # i.e. no shifting to be done

# Decide whether or not to print some intermediate results for debugging purposes:
""" Some of the more important/useful intermediate results are printed by default. """
#Debug = True
Debug = False



# Define the REST variables that correspond to the Source ROI Collection,
# the corresponding Source DICOMs, and the Target DICOMs that the ROI
# Collection is to be copied to:
if Subject == 'ACRIN-FMISO-Brain-011' and CopyFrom == 'MR':
    SourceRoiProjLabel = 'ACRIN'
    SourceRoiSubjLabel = 'ACRIN-FMISO-Brain-011'
    SourceRoiExpLabel = 'ACRIN-FMISO-Brain-011_MR_2'
    SourceRoiLabelDateTime = '20200330'
    SourceSeriesNo = '5'
    SourceSeriesDesc = 'T2 3D Slab'

if Subject == 'ACRIN-FMISO-Brain-011' and CopyTo == 'PET':
    TargetDicomProjLabel = 'ACRIN'
    TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
    TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_CT_PET_1'
    TargetSeriesNo = '3'
    TargetSeriesDesc = '3D_Brain miso '
    
if Subject == 'ACRIN-FMISO-Brain-011' and CopyTo == 'CT':
    TargetDicomProjLabel = 'ACRIN'
    TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
    TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_CT_19610312'
    TargetSeriesNo = '2'
    TargetSeriesDesc = 'HEAD 2.4 H40 SOFT TISSUE'
    
if Subject == 'ACRIN-FMISO-Brain-011' and CopyFrom == 'MR_4' and CopyTo == 'MR_12':
    SourceRoiProjLabel = 'ACRIN'
    SourceRoiSubjLabel = 'ACRIN-FMISO-Brain-011'
    SourceRoiExpLabel = 'ACRIN-FMISO-Brain-011_MR_4'
    SourceRoiLabelDateTime = '20200511'
    SourceSeriesNo = '8'
    SourceSeriesDesc = 'T1 SE AXIAL POST FS FC'
    
    TargetDicomProjLabel = 'ACRIN'
    TargetDicomSubjLabel = 'ACRIN-FMISO-Brain-011'
    TargetDicomExpLabel = 'ACRIN-FMISO-Brain-011_MR_12'
    TargetSeriesNo = '8'
    TargetSeriesDesc = 'T1 SE AXIAL POST FS FC'
    

# Define the root directory where DICOMs will be downloaded from XNAT to:
""" Note: Files will be downloaded to RootDir\TodaysDate SourceRoiProjLabel"""
RootDir = r'C:\Temp'

# Get the data from XNAT:
DataDict = GetXnatData(XnatAddress=XnatAddress,
                       XnatUsername=XnatUsername,
                       XnatPassword=XnatPassword,
                       RootDir=RootDir,
                       CopyFrom=CopyFrom,
                       CopyTo=CopyTo,
                       SourceRoiProjLabel=SourceRoiProjLabel,
                       SourceRoiSubjLabel=SourceRoiSubjLabel,
                       SourceRoiExpLabel=SourceRoiExpLabel,
                       SourceRoiLabelDateTime=SourceRoiLabelDateTime,
                       SourceSeriesNo=SourceSeriesNo,
                       SourceSeriesDesc=SourceSeriesDesc,
                       TargetDicomProjLabel=TargetDicomProjLabel,
                       TargetDicomSubjLabel=TargetDicomSubjLabel,
                       TargetDicomExpLabel=TargetDicomExpLabel,
                       TargetSeriesNo=TargetSeriesNo,
                       TargetSeriesDesc=TargetSeriesDesc,
                       Debug=Debug)



""" 11/05/2020:
This might need modifying for 3D registrations: """
if False:
    # Run CopyRoi():

    # Define the scan direction:
    ScanDirection = 'z'

    # Select the amount to shift the ScanDirection positions of the Target DICOMs by:
    """ 
    e.g. If ScanDirection = 'z', and Shift = 3, the z component of the Image Position
    Patient for the Target DICOMs will be shifted by 3 mm, and the resulting DICOMs
    stored in a new list called ShiftedDicoms. If Shift = 0, no shift will be applied. 
    """
    # For Sarcoma dataset, use the difference between the Distance Source to Detector 
    # (949.075) and Table Height (170) of the CT scan as the z-position correction:
    #Shift = 949.075 - 170

    # For ACRIN dataset, the Distance Source to Detector and Table Height tags are 
    # empty, so need to use anatomical features instead: 
    if ShiftFeature and Subject == 'ACRIN-FMISO-Brain-011':
        if CopyFrom == 'MR' and CopyTo == 'PET':
            if 'nose' in ShiftFeature:
                # Using the z positions where visually the "peak" of the patient's nose could 
                # be seen maximally:
                ##Shift = -117.82 - (- 52.45) # = z_PET - z_MR (slice 43/47 and 133/160) <-- this one is incorrect
                Shift = -52.45 - (- 117.82) # = z_MR - z_PET (slice 133/160 and 43/47) <-- this is correct

            if 'infinity' in ShiftFeature:
                # or using the infinity-like structure in the brain:
                #Shift = -101.47 - (- 36.45) # = z_PET - z_MR (slice 38/47 and 117/160) <-- this one is incorrect
                Shift = -36.45 - (- 101.47) # = z_MR - z_PET (slice 117/160 and 48/47) <-- this is correct

        if CopyFrom == 'MR' and CopyTo == 'CT':
            if 'lens' in ShiftFeature:
                # Using the slices that maximise the size of the eye's lens:
                Shift = -21.450044136494 - 136.74918721324 # z_MR - z_CT (slice 102/160 and 10/72) = -158.199231349734

            if 'butterfly' in ShiftFeature:
                # Using the slices that show the butterfly-like structure with some dark fine structures within them:   
                Shift = 14.549955863506 - 204.33718721324 # z_MR - z_CT (slice 66/160 and 39/72) = -189.787231349734

    # Don't apply any shift:
    if not ShiftFeature:
        Shift = 0

    # Select the registration method(s) to run:
    """ 
    Note: One or more registrations will be applied in a for loop depending on the 
    number of strings in RegMethods.  RegMethods must be an array of one or more
    strings.  If it contains an empty string, i.e. '', no registration will be applied. 
    """
    #RegMethods = ['rigid', 'affine', 'non-rigid']
    #RegMethods = ['rigid']
    RegMethods = ['affine']
    #RegMethods = ['non-rigid']

    # Don't apply any registration:
    #RegMethods = [''] 

    # Decide whether or not to delete the Source and Target DICOMs downloaded from
    # XNAT when done:
    #DelDicoms = True
    DelDicoms = False

    # Decide whether or not to re-run registrations if the registered DICOMs already 
    # exist (i.e. have already been exported locally):
    """ Since registrations are time consuming this can be helpful. """
    #RerunReg = True
    RerunReg = False

    # Decide whether or not to print some intermediate results for debugging purposes:
    """ Some of the more important/useful intermediate results are printed by default. """
    #Debug = True
    Debug = False

    # Run CopyRoi(), looping through each RegMethod in RegMethods:
    for RegMethod in RegMethods:
        if Shift and RegMethod:
            print('\n\n\nCopying ROI Collection from', SourceSeriesDesc, 'to',                   TargetSeriesDesc, ',\n applying shift to', ScanDirection + '-direction and',                   RegMethod, 'registration...\n\n\n')
        elif Shift and not RegMethod:
            print('\n\n\nCopying ROI Collection from', SourceSeriesDesc, 'to',                   TargetSeriesDesc, ',\n applying shift to', ScanDirection + '-direction ...\n\n\n')
        elif not Shift and RegMethod:
            print('\n\n\nCopying ROI Collection from', SourceSeriesDesc, 'to',                   TargetSeriesDesc, ',\n applying', RegMethod, 'registration...\n\n\n')
        elif not Shift and not RegMethod:
            print('\n\n\nCopying ROI Collection from', SourceSeriesDesc, 'to',                   TargetSeriesDesc, '...\n\n\n')

        SourceDicoms, SourceRois,        TargetDicoms, TargetRois,        ShiftedDicoms, ShiftedRois,        RemappedSourceDicoms, RemappedSourceRois,        RemappedTargetDicoms, RemappedTargetRois,        FixedDicoms, FixedRois,        MovingDicoms, MovingRois,        RegisteredDicoms, RegisteredRois,        DataDict = CopyRoi(DataDict=DataDict,
                           SearchBy=ScanDirection,
                           Shift=Shift,
                           ShiftFeature=ShiftFeature,
                           RegMethod=RegMethod,
                           DelDicoms=DelDicoms,
                           RerunReg=RerunReg)

        """ Note:  RemappedTargetDicoms will actually be RemappedShiftedTargetDicoms if Shift != 0. """


# In[40]:


DataDict


# ### Export DataDict:

# In[10]:


import json

# Create filepath for json file:
JsonFpath = os.path.join(DataDict['RootDir'], 'DataDict.json')

# Export DataDict:
data = json.dumps(DataDict)
file = open(JsonFpath, 'w')
file.write(data)
file.close()


# ### Import DataDict:

# In[ ]:


ImportDir = r'C:\\Temp\\2020-05-11 ACRIN MR_4 to MR_12'
JsonFpath = os.path.join(ImportDir, 'DataDict.json')

file = open(JsonFpath)
DataDict = json.load(file)
file.close()

DataDict


# In[34]:


interact(display_images_v3,
         fixed_image_z=(0,FixedIm.GetSize()[2]-1),
         moving_image_z=(0,MovingIm.GetSize()[2]-1),
         registered_image_z=(0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(MovingIm)),
         registered_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_title = 'Fixed image',
         moving_title = 'Moving image',
         registered_title = 'Registered image');


# ### Define FixedDir and MovingDir based on the directories contained in DataDict:

# In[12]:


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


# ### Use SimpleElastix to register MovingIm to FixedIm:

# In[100]:


# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ElastixImFilt = sitk.ElastixImageFilter()
ElastixImFilt.SetFixedImage(FixedIm)
ElastixImFilt.SetMovingImage(MovingIm)
#ElastixImFilt.SetMovingImage(FixedIm)
#ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('groupwise'))
ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ElastixImFilt.LogToConsoleOn()
ElastixImFilt.LogToFileOn()
RegIm = ElastixImFilt.Execute()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')


# ### Seems fine. So continue where I left off last Thursday.  The registration looked fine.  Now need to incorporate the contours in the 3D registration work.
# 
# ### Use GetDicoms() and GetContourPoints3D() to get the contour points:

# In[84]:


# Use GetDicoms() and GetContourPoints3D() to get the contour points:

import pydicom
import GetDicoms
importlib.reload(GetDicoms)
from GetDicoms import GetDicoms
import GetContourPoints3D
importlib.reload(GetContourPoints3D)
from GetContourPoints3D import GetContourPoints3D

# Get the Fixed (Source) DICOMs:
FixedDicomFpaths, FixedDicoms = GetDicoms(DataDict['Source']['DicomDir'], 'slices', Debug=False)

# Read in the ROI Collection for FixedDicoms:
FixedRois = pydicom.dcmread(DataDict['Source']['RoiFpath'])

# Create an array of contour points for each DICOM in FixedDicoms:
FixedContourPts = []

for i in range(len(FixedDicoms)):
    FixedContourPts.append(GetContourPoints3D(FixedDicoms[i], FixedRois))

FixedContourPts


# In[192]:


len(FixedContourPts)


# ### Write the points to a text file with .pts extension as used for elastix:
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[103]:


""" FixedContourPts are in the form:
[
 [],
 [],
 ...
 [
  [point1_x_contourN, point1_y_contourN, point1_z_contourN],
  [point2_x_contourN, point2_y_contourN, point2_z_contourN],
  ...
  [pointM_x_contourN, pointM_y_contourN, pointM_z_contourN]
  ], 
 [],
 ]
 
where each array corresponds to a slice in FixedIm. If any given array is empty, 
i.e. [], there is no contour for that slice.  If there is a contour, the points 

  [point1_x_contourN, point1_y_contourN, point1_z_contourN],
  [point2_x_contourN, point2_y_contourN, point2_z_contourN],
  ...
  [pointM_x_contourN, pointM_y_contourN, pointM_z_contourN]
  
 are the M points that make up a contour for the nth slice.
 
NOTE:  I'M NOT SURE WHETHER THE FUNCTION GetContourPoints3D() WILL CORRECTLY
       DEAL WITH CASES WHERE THERE ARE MORE THAN ONE CONTOUR FOR ANY GIVEN SLICE.
       LIKEWISE I'M NOT SURE THAT THE DEFINITION OF FixedContourPts ABOVE WILL
       BE CORRECT.  AND THE CODE BELOW MAY NEED MODIFICATION TOO...
 
Elastix requires the text file of points to be of the form:

point (or) index
number_of_points (or) number_of_indeces
point1_x point1_y point1_z
point2_x point2_y point2_z
...

Hence the points in FixedContourPts will need to be flattened.
"""

# Flatten FixedContourPts:
#FixedContourPtsFlat = []

# Count the number of contours and 3D points, and the mapping 
# of each point to the slice it belongs to:
FixedNoCts = 0
FixedNoPts = 0
FixedPtsToSliceNo = []

for SliceNo in range(len(FixedContourPts)):
    # The array for this slice:
    SliceArray = FixedContourPts[SliceNo]
    
    # If SliceArray is not empty there are contour points:
    if SliceArray:
        #FixedContourPtsFlat.extend(SliceArray)
        
        # Increment the number of contours:
        FixedNoCts = FixedNoCts + 1
        
        # Increment the number of points by looping through
        # each array of 3D points within each SliceArray:
        for PointsArray in SliceArray:
            FixedNoPts = FixedNoPts + 1
            
            # Append SliceNo to FixedPtsToSliceNo for this SliceNo:
            FixedPtsToSliceNo.append(SliceNo)
        
#FixedContourPtsFlat
print(FixedNCts)
print(FixedNPts)
print(FixedPtsToSliceNo)


# In[104]:


# Define the filename for the FixedContoutPts:
#FixedContourPtsFname = 'FixedContourPts.pts'
FixedContourPtsFname = 'FixedContourPts.txt'

# Open a text file:
TextFile = open(FixedContourPtsFname, 'w')

# Write 'point' to the first line since they are points, not indeces:
TextFile.write('point')

# Write FixedNPts to the second line:
TextFile.write(f'\n{FixedNPts}')

for ContourArray in FixedContourPts:
    if ContourArray:
        # Loop through each array of 3D points 
        # within each ContourArray and add the
        # points in PointsArray on a new line
        # with space between points:
        for PointsArray in ContourArray:
            TextFile.write(f'\n{PointsArray[0]} {PointsArray[1]} {PointsArray[2]}')
            
TextFile.close()


# In[54]:


FixedDicoms[0].ImageOrientationPatient


# In[55]:


FixedDicoms[0].ImagePositionPatient


# In[56]:


FixedDicoms[0].PixelSpacing


# In[32]:


FixedIm.GetOrigin()


# In[ ]:





# In[86]:


interact(display_images_and_contours,
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


# In[23]:


len(FixedContourPts)


# In[62]:


FixedIm.GetSize()


# In[63]:


FixedIm.GetSize()[2]-1


# In[65]:


len(FixedDicoms)


# In[66]:


len(FixedContourPts)


# In[74]:


ImFilter.PrintParameterMap()


# In[108]:


print(FixedIm.GetOrigin())
print(MovingIm.GetOrigin())
print(RegIm.GetOrigin())


# ### Transform the contour points:
# 
# https://simpleelastix.readthedocs.io/PointBasedRegistration.html

# In[106]:


TransformixImFilt = sitk.TransformixImageFilter()
TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
TransformixImFilt.SetFixedPointSetFileName(FixedContourPtsFname)
# Need to explicitely tell elastix that the image is 3D since by default it will
# only transform to 2D:
TransformixImFilt.SetMovingImage(MovingIm)
TransformixImFilt.Execute()


# ### Parse the transformed point set:

# In[112]:


import csv

# The filename for the deformed points is:
TransformedPtsFname = 'outputpoints.txt'

# Open the text file:
TextFile = open(TransformedPtsFname, 'r')

Reader = csv.reader(TextFile)

# Iterate through each line:
for i, line in enumerate(Reader):
    print(line[0].split('OutputPoint'))
    #print(line[0])


TextFile.close()


# In[117]:


# Open the text file:
TextFile = open(TransformedPtsFname, 'r')

for line in TextFile:
    print(line)


# In[142]:


list(TextFile)


# In[185]:


# Open the text file:
TextFile = open(TransformedPtsFname, 'r')

# Search strings:
PointSearchStr = 'OutputPoint = [ '
DeformationSearchStr = 'Deformation = [ '

# Store the transformed points and deformation values:
TransformedPts = []
DeformationVals = []

# Store all lines:
Lines = []

# Iterate through each line in TextFile:
for line in TextFile:
    #print(line)
    #print(line.find(SearchStr))
    
    Lines.append(line)
    
    # Get starting index for the point after PointSearchStr:
    StartInd = line.find(PointSearchStr) + len(PointSearchStr)
    
    # Text that remains after StartInd:
    RemText = line[StartInd::]
    
    # Get index of closing bracket for point:
    EndInd = StartInd + RemText.find(']') - 1
    
    PointsWithSpaces = line[StartInd:EndInd]
    
    #print('PointsWithSpaces =' , PointsWithSpaces, '\n')
    
    #PointsWithoutSpaces = PointsWithSpaces.replace(' ', ', ')
    
    #print('PointsWithoutSpaces =' , PointsWithoutSpaces, '\n')
    
    #TransformedPts.append([float(string) for string in PointsWithoutSpaces])
    

    # Get the points:
    points = []

    # Get the points in PointsWithSpaces until there are
    # no characters remaining:
    while PointsWithSpaces and len(points) < 3:
        #print('\nPointsWithSpaces         =' , PointsWithSpaces)
        
        # Strip away any whitespaces on the left:
        PointsWithSpaces = PointsWithSpaces.lstrip()
        
        #print('after lstrip()           =' , PointsWithSpaces)

        # Find the next whitespace:
        ind = PointsWithSpaces.find(' ')
        
        #print('ind of next whitespace   =', ind)
        
        # If no whitespaces remain, simply append all characters:
        if ind==-1:
            # Append all characters:
            points.append(PointsWithSpaces)
            
            #print('point appended to points =', PointsWithSpaces)
            
            # Effectively pop all characters in PointsWithSpaces:
            PointsWithSpaces = ''
        else:
            # Append the characters from 0 to ind:
            points.append(PointsWithSpaces[0:ind])
        
            #print('point appended to points =', PointsWithSpaces[0:ind])

            # Effectively pop PointsWithSpaces from 0 to ind+1:
            PointsWithSpaces = PointsWithSpaces[ind+1::]
        
        #print('after pop                =', PointsWithSpaces)


    # Append points to TransformedPts:
    #TransformedPts.append(points)
    TransformedPts.append([float(point) for point in points])
    
    #print('\npoints added to TransformedPts =', points)
    #print('')
        
#print('TransformedPts = ', TransformedPts)

TextFile.close()


# In[186]:


TransformedPts


# for item in enumerate(PointsWithSpaces):
#     print(item)

# list(PointsWithSpaces)

# In[162]:


# This didn't work:

data = np.genfromtxt(TransformedPtsFname, delimiter='\n')

data


# In[118]:


# This didn't work:

import re

# Open the text file:
TextFile = open(TransformedPtsFname, 'r')

for line in TextFile:
    #print(line)
    print(line.rstrip('OutputPoint = ['))
    print(line.find('OutputPoint = ['))


# ### Use mapping FixedPtsToSliceNo to convert TransformedPts to the analog of FixedContourPts.
# 
# NB Each row in FixedContourPts has an array that is either empty, i.e. [], or has an array of points.
# 
# FixedPtsToSliceNo maps the point number to the slice number. i.e. 
# 
# [13,
#  13,
#  13,
#  14,
#  14,
#  ...
#  ]
#  
#  => Points 1-3 belong to slice 13, Points 4, 5... belong to slice 14.

# In[188]:


print(FixedPtsToSliceNo)


# In[191]:


FixedPtsToSliceNo.index(13)


# In[ ]:


# Initialise MovingContourPts:
MovingContourPts = []

for SliceNo in range(len(FixedContourPts)):
    # Initialise Points for this SliceNo:
    Points = []
    
    
    if SliceNo in FixedPtsToSliceNos:
        while FixedPtsToSliceNos:
            # Get the index of the next matching SliceNo:
            ind = FixedPtsToSliceNo.index(SliceNo)

            if ind == -1:
                FixedPtsToSlices = []
            else:
                # Pop the p given by ind:
                item = FixedPtsToSlices.pop(ind)
                
                # Append the 
                Points.append(FixedPtsToSlices[ind])
                FixedPtsToSlices.pop(ind)
    else:
        MovingContourPts.append(Points)


# ### There might be an easier way...

# In[201]:


# The number of slices in FixedContourPts:
NSlices = len(FixedContourPts)

# Initialise MovingContourPts:
MovingContourPts = []

# Initialise an interger
i = 0

# Loop for each row (slice) in FixedContourPts:
for SliceNo in range(NSlices):
    # Get the points in FixedContourPts for this SliceNo:
    points = FixedContourPts[SliceNo]
    
    # If points is not []:
    if points:
        # Get the number of points that belong to this SliceNo:
        N = len(points)
        
        # Append the next N points from TransformedPts to MovingContourPts:
        MovingContourPts.append(TransformedPts[i:i+N])
        
        # Update i:
        i = i + N + 1
    
    else:
        # Append an empty array to MovingContourPts:
        MovingContourPts.append([])


if not i == Nsliprint()

MovingContourPts


# In[202]:


# Check:

print('len(FixedDicoms)      =', len(FixedDicoms))
print('len(FixedContourPts)  =', len(FixedContourPts))
print('len(MovingContourPts) =', len(MovingContourPts))

for r in range(len(FixedContourPts)):
    print('\nr =', r)
    print('len(FixedContourPts[r])  =', len(FixedContourPts[r]))
    print('len(MovingContourPts[r]) =', len(MovingContourPts[r]))


# ### 15/05: It's better than it was but r = 17 has 99 v 95 points in Fixed v Moving, so something is still not right.
# 
# ### However more importantly the approach is flawed because the contour points in any given row in FixedContourPts (that correspond to a particular slice in FixedDicoms) do not relate to the same row in MovingContourPts - and moreover, there are different numbers of slices!

# In[204]:


print('FixedIm.GetSize() =', FixedIm.GetSize())
print('MovingIm.GetSize() =', MovingIm.GetSize())


# ### So instead the transformed contours (contour points) need to be re-sliced (or re-sampled) using the voxels of MovingIm.

# In[209]:


FixedContourPts[13][0]


# In[208]:


TransformedPts[0]


# ### I doubt that the transformed points are what they should be.  There's a huge displacement in the z coordinate.

# In[212]:


print('FixedIm.GetOrigin()  =', FixedIm.GetOrigin())
print('MovingIm.GetOrigin() =', MovingIm.GetOrigin())


# In[ ]:





# ### 15/05:  Use SimpleElastix to register MovingIm to FixedIm but this time resample MovingIm:
# 
# Following the code by nwschurink on Jan 25 2019 as a guide (although he applied the resampling to the fixed mask rather than to the moving image:
# 
# https://github.com/SuperElastix/SimpleElastix/issues/261

# In[232]:


# Start timing:
times = []
times.append(time.time())

# Register the 3D images:
ElastixImFilt = sitk.ElastixImageFilter()
#ElastixImFilt = sitk.SimpleElastix() # AttributeError: 
""" module 'SimpleITK' has no attribute 'SimpleElastix' 
See https://github.com/SuperElastix/SimpleElastix/issues/86
"""
ElastixImFilt.SetFixedImage(FixedIm)
# Resample MovingIm:
Resampler = sitk.ResampleImageFilter()
Resampler.SetReferenceImage(FixedIm)
ResampledMovingIm = Resampler.Execute(MovingIm)
#ElastixImFilt.SetMovingImage(MovingIm)
ElastixImFilt.SetMovingImage(ResampledMovingIm)
#ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('groupwise'))
ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
ElastixImFilt.LogToConsoleOn()
ElastixImFilt.LogToFileOn()
#ElastixImFilt.Execute()

#ElastixImFilt.SetParameter("Registration", "AdvancedImageToImageMetric"); # <-- get error
#ElastixImFilt.SetParameter("Interpolator", "BSplineInterpolator"); # <-- ok
#ElastixImFilt.SetParameter("Metric", ("NormalizedMutualInformation", "CorrespondingPointsEuclideanDistanceMetric")); # <-- get error
#ElastixImFilt.SetParameter("Metric0Weight", "0.0");

ElastixImFilt.Execute()

RegIm = ElastixImFilt.GetResultImage()

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
print(f'\nTook {Dtime} s to register the 3D image stacks.')

# Trying this in response to the error:
"""
RuntimeError: Exception thrown in SimpleITK TransformixImageFilter_Execute: C:\SimpleElastix\source\SimpleElastix\Code\Elastix\src\sitkTransformixImageFilterImpl.cxx:116:
sitk::ERROR: 
itk::ExceptionObject (000000085FBE7D50)
Location: "unknown" 
File: C:\SimpleElastix\build_20200108\Elastix\Core\Main\elxTransformixFilter.hxx
Line: 272
Description: itk::ERROR: Self(000001C0EF06E950): Empty parameter map in parameter object.
"""
#ParameterMap = ElastixImFilt.GetParameterMap()

times.append(time.time())

if False:
    # Get the deformation field:
    TransformixImFilt = sitk.TransformixImageFilter()
    #TransformixImFilt = sitk.SimpleTransformix()
    TransformixImFilt.SetMovingImage(ResampledMovingIm)
    TransformixImFilt.ComputeDeformationFieldOn()
    #TransformixImFilt.LogToConsole() # get AttributeError:
    """ 'TransformixImageFilter' object has no attribute 'LogToConsole' """
    TransformixImFilt.LogToFileOn()
    # Trying this in response to the error:
    """
    RuntimeError: Exception thrown in SimpleITK TransformixImageFilter_Execute: C:\SimpleElastix\source\SimpleElastix\Code\Elastix\src\sitkTransformixImageFilterImpl.cxx:116:
    sitk::ERROR: 
    itk::ExceptionObject (000000085FBE7D50)
    Location: "unknown" 
    File: C:\SimpleElastix\build_20200108\Elastix\Core\Main\elxTransformixFilter.hxx
    Line: 272
    Description: itk::ERROR: Self(000001C0EF06E950): Empty parameter map in parameter object.
    """
    #TransformixImFilt.SetParameterMap(ParameterMap)
    TransformixImFilt.Execute()

    Deform = TransformixImFilt.GetDeformationField()

    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\nTook {Dtime} s to get the deformation field.')


# In[228]:


print('FixedIm.GetSize()           =', FixedIm.GetSize())
print('MovingIm.GetSize()          =', MovingIm.GetSize())
print('ResampledMovingIm.GetSize() =', ResampledMovingIm.GetSize())
print('RegIm.GetSize()             =', RegIm.GetSize())


# In[229]:


interact(display_images_and_contours,
         fixed_ind = (0,FixedIm.GetSize()[2]-1),
         moving_ind = (0,ResampledMovingIm.GetSize()[2]-1),
         reg_ind = (0,RegIm.GetSize()[2]-1),
         fixed_npa = fixed(sitk.GetArrayFromImage(FixedIm)),
         moving_npa = fixed(sitk.GetArrayFromImage(ResampledMovingIm)),
         reg_npa = fixed(sitk.GetArrayFromImage(RegIm)),
         fixed_im = fixed(FixedIm),
         moving_im = fixed(ResampledMovingIm),
         reg_im = fixed(RegIm),
         fixed_contour_pts = fixed(FixedContourPts),
         fixed_title = 'Fixed image',
         moving_title = 'Resampled Moving image',
         reg_title = 'Registered image');


# In[ ]:




