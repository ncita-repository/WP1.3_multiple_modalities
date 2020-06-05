# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:18:10 2020

@author: ctorti
"""



""" Define utility functions for plotting of registration results """


""" Source:
https://github.com/InsightSoftwareConsortium/SimpleITKWorkshopImageJ2015/blob/master/6_Similarity_Metrics_and_Optimization.ipynb
"""


# Import packages:
import matplotlib.pyplot as plt
import SimpleITK as sitk
import time
import numpy as np
import os
import re





def RunSimpleElastixReg(FixedIm, MovingIm):
    # Start timing:
    times = []
    times.append(time.time())
    
    # Initiate ElastixImageFilter:
    ElastixImFilt = sitk.ElastixImageFilter()
    #ElastixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    ElastixImFilt.LogToConsoleOff() # <-- no output in Jupyter
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
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    # Set the parameter map:
    ElastixImFilt.SetParameterMap(ElastixParamMap)
    
    print('Performing registration...')
    
    # Register the 3D images:
    ElastixImFilt.Execute()
    # Get the registered image:
    RegIm = ElastixImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\nTook {Dtime} s to register the 3D image stacks.')
    
    return RegIm, ElastixImFilt
    


def TransformPoints(MovingImage, TransformFilter):
    # Initiate TransformixImageFilter:
    TransformixImFilt = sitk.TransformixImageFilter()
    #TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TransformixImFilt.LogToConsoleOff() # <-- no output in Jupyter
    TransformixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    # Get the transform parameter map. Start by getting the transfer parameter map
    # used for the registration (= ElastixParamMap):
    #TransformixParamMap = ImageFilter.GetParameterMap() 
    
    # Set the parameter map:
    TransformixImFilt.SetTransformParameterMap(TransformFilter.GetTransformParameterMap())
    #TransformixImFilt.SetTransformParameterMap(Transform)
    
    # Set the input points (assumed to be in the current working directory):
    TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
    
    # Need to explicitely tell elastix that the image is 3D since by default it 
    # will only transform to 2D:
    TransformixImFilt.SetMovingImage(MovingImage)
    
    # Set up for computing deformation field:
    #TransformixImFilt.ComputeDeformationFieldOn()
    
    # Transform the contour points:
    TransformixImFilt.Execute()
    
    #CWD = os.getcwd()
    
    print('\nThe points in inputpoints.txt have been transformed to outputpoints.txt.')
    
    return
    

def ParseOutputPoints():
    """ THIS NEEDS TO BE CORRECTED! """
    
    # Define the patterns to search:
    """ 2D example:
        OutputPtsSearch = re.compile("OutputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s*\]")
    """
    
    Point_pattern = re.compile("Point\s*([\d]+)")
    #InputIndex_pattern = re.compile("InputIndex\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    InputIndex_pattern = re.compile("InputIndex\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s*\]")
    
    InputPoint_pattern = re.compile("InputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    OutputIndexFixed_pattern = re.compile("OutputIndexFixed\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    OutputPoint_pattern = re.compile("OutputPoint\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    Deformation_pattern = re.compile("Deformation\s*=\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]")
    OutputIndexMoving_pattern = re.compile("OutputIndexMoving\s*=\s*\[\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*\]")
    
    Point = []
    InputIndex = []
    InputPoint = []
    OutputIndexFixed = []
    OutputPoint = []
    Deformation = []
    OutputIndexMoving = []
    
    with open("outputpoints.txt", "r") as file:
        for line in file:
            Parsed = Point_pattern.search(line)
            PointNo = int(Parsed.group(1))
            Point.append(PointNo)
            
            Parsed = InputIndex_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            InputIndex.append(Indices)
        
            Parsed = InputPoint_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            InputPoint.append(Points)
            
            Parsed = OutputIndexFixed_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            OutputIndexFixed.append(Indices)
            
            Parsed = OutputPoint_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            OutputPoint.append(Points)
            
            Parsed = Deformation_pattern.search(line)
            Points = (float(Parsed.group(1)), float(Parsed.group(2)), float(Parsed.group(3)))
            Deformation.append(Points)
            
            Parsed = OutputIndexMoving_pattern.search(line)
            Indices = (int(Parsed.group(1)), int(Parsed.group(2)), int(Parsed.group(3)))
            OutputIndexMoving.append(Indices)
            
    return {'Point':Point, 'InputIndex':InputIndex, 'InputPoint':InputPoint,
            'OutputIndexFixed':OutputIndexFixed, 'OutputPoint':OutputPoint,
            'Deformation':Deformation, 'OutputIndexMoving':OutputIndexMoving}
    


def ComputeDeformationField(MovingImage, TransformFilter):
    # Get the current working directory:
    #CurrentWorkingDir = os.getcwd()

    # Initiate TransformixImageFilter:
    TransformixImFilt = sitk.TransformixImageFilter()
    #TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TransformixImFilt.LogToConsoleOff()
    TransformixImFilt.LogToFileOn() 
    #TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output
    
    # Set up for computing deformation field:
    TransformixImFilt.ComputeDeformationFieldOn() # <-- no output in Jupyter
    
    # Set up other computations:
    TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
    TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'
    
    # Set the transform parameter map:
    TransformixImFilt.SetTransformParameterMap(TransformFilter.GetTransformParameterMap())
    
    # Need to explicitely tell elastix that the image is 3D since by default it
    # will only transform to 2D:
    TransformixImFilt.SetMovingImage(MovingImage)
    
    # Transform the contour points:
    TransformixImFilt.Execute()
    #TransformixIm = TransformixImFilt.Execute()
    
    # Get the deformation field:
    DefField = TransformixImFilt.GetDeformationField()
    
    CWD = os.getcwd()
    
    # Check if deformation file was exported:
    if not os.path.isfile(os.path.join(CWD, 'deformationField.nii')):
        print('\nSeems that the deformation field was not exported.',
              '\nManually exporting file...')
        
        # Convert the deformation field to a numpy array:
        """ Note: 
            Executing this kills the kernel in Anaconda Jupyter. 
        """
        DefFieldNda = sitk.GetArrayFromImage(DefField) 
        
        # Write the deformation field to a .nii file:
        """ Note: 
            Executing this also kills the kernel in Anaconda Jupyter. 
        """
        DefFieldNiiFpath = os.path.join(CWD, 'deformationField.nii')
        
        sitk.WriteImage(DefField, DefFieldNiiFpath)
        
        print('\nDeformation field was exported to', DefFieldNiiFpath)
    
    return DefField



def ShowDefFieldInfo(DeformationField):
    # Get some info on DeformationField:
    DefFieldDtype = DeformationField.GetPixelIDTypeAsString()
    DefFieldSize = DeformationField.GetSize()
    DefFieldSpacing = DeformationField.GetSpacing()
    DefFieldOrigin = DeformationField.GetOrigin()
    """ Note: 
        Executing this kills the kernel in Anaconda Jupyter. 
    """
    DefFieldNda = sitk.GetArrayFromImage(DeformationField)
    
    print('\nDefField Dtype    =', DefFieldDtype)
    print('DefField Size     =', DefFieldSize)
    print('DefField Spacing  =', DefFieldSpacing)
    print('DefField Origin   =', DefFieldOrigin)
    print(f'DefField Min, Max = {np.min(DefFieldNda)}, {np.max(DefFieldNda)}')

    return


def ShowImagesInfo(FixedIm, MovingIm):
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
    
    print(f'\nFixedIm Min, Max value  = {np.min(FixedNda)}, {np.max(FixedNda)}')
    print(f'MovingIm Min, Max value = {np.min(MovingNda)}, {np.max(MovingNda)}')
    
    print(f'\nFixedIm Min, Max x-pos  = {FixedOrigin[0]}, {FixedOrigin[0] + FixedSize[0]*FixedSpacing[0]}')
    print(f'FixedIm Min, Max y-pos  = {FixedOrigin[1]}, {FixedOrigin[1] + FixedSize[1]*FixedSpacing[1]}')
    if FixedSize[2] == 3:
        print(f'FixedIm (min, max) z-pos  = {FixedOrigin[2]}, {FixedOrigin[2] + FixedSize[2]*FixedSpacing[2]}')
    
    print(f'\nMovingIm Min, Max x-pos  = {MovingOrigin[0]}, {MovingOrigin[0] + MovingSize[0]*MovingSpacing[0]}')
    print(f'MovingIm Min, Max y-pos  = {MovingOrigin[1]}, {MovingOrigin[1] + MovingSize[1]*MovingSpacing[1]}')
    if MovingSize[2] == 3:
        print(f'MovingIm Min, Max z-pos  = {MovingOrigin[2]}, {MovingOrigin[2] + MovingSize[2]*MovingSpacing[2]}')
    
    return
        


def ShowRegResults(FixedIm, MovingIm, RegIm):
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
    
    print(f'\nFixedIm Min, Max value  = {np.min(FixedNda)}, {np.max(FixedNda)}')
    print(f'MovingIm Min, Max value = {np.min(MovingNda)}, {np.max(MovingNda)}')
    print(f'RegIm Min, Max value    = {np.min(RegNda)}, {np.max(RegNda)}')
    
    print(f'\nFixedIm Min, Max x-pos  = {FixedOrigin[0]}, {FixedOrigin[0] + FixedSize[0]*FixedSpacing[0]}')
    print(f'FixedIm Min, Max y-pos  = {FixedOrigin[1]}, {FixedOrigin[1] + FixedSize[1]*FixedSpacing[1]}')
    if FixedSize[2] == 3:
        print(f'FixedIm Min, Max z-pos  = {FixedOrigin[2]}, {FixedOrigin[2] + FixedSize[2]*FixedSpacing[2]}')
    
    print(f'\nMovingIm Min, Max x-pos  = {MovingOrigin[0]}, {MovingOrigin[0] + MovingSize[0]*MovingSpacing[0]}')
    print(f'MovingIm Min, Max y-pos  = {MovingOrigin[1]}, {MovingOrigin[1] + MovingSize[1]*MovingSpacing[1]}')
    if MovingSize[2] == 3:
        print(f'MovingIm Min, Max z-pos  = {MovingOrigin[2]}, {MovingOrigin[2] + MovingSize[2]*MovingSpacing[2]}')
    
    print(f'\nRegIm Min, Max x-pos  = {RegOrigin[0]}, {RegOrigin[0] + RegSize[0]*RegSpacing[0]}')
    print(f'RegIm Min, Max y-pos  = {RegOrigin[1]}, {RegOrigin[1] + RegSize[1]*RegSpacing[1]}')
    if RegSize[2] == 3:
        print(f'RegIm Min, Max z-pos  = {RegOrigin[2]}, {RegOrigin[2] + RegSize[2]*RegSpacing[2]}')
    
    return
    


""" Various displaying functions modified from "display_images" here:
    # Source: https://github.com/SimpleITK/SPIE2019_COURSE
"""

def display_images(fix_ind, mov_ind, fix_npa, mov_npa, fix_title, mov_title):
    # Get min, max:
    #fix_min = np.min(fix_npa)
    #fix_max = np.max(fix_npa)
    #mov_min = np.min(mov_npa)
    #mov_max = np.max(mov_npa)

    # create a figure with two subplots and the specified size
    plt.subplots(1,2,figsize=(10,8))

    # draw the fixed image in the first subplot
    plt.subplot(1,2,1)
    plt.imshow(fix_npa[fix_ind,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r, vmin=fix_min, vmax=fix_max);
    #plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(fix_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,2,2)
    plt.imshow(mov_npa[mov_ind,:,:],cmap=plt.cm.Greys_r);
    #plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r, vmin=mov_min, vmax=mov_max);
    #plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r, norm=True);
    #plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r, vmin=0, vmax=255);
    plt.title(mov_title)
    plt.axis('off')
    
    

def display_images_with_contours(fix_ind, mov_ind,
                                 fix_npa, mov_npa,
                                 fix_im, mov_im,
                                 mov_contour_pts,
                                 fix_title, mov_title):
    
    # Get the number of slices in the images:
    fix_Nslices = fix_im.GetSize()[2]
    mov_Nslices = mov_im.GetSize()[2]

    # Get the z positions of the origins:
    fix_origin_z = fix_im.GetOrigin()[2]
    mov_origin_z = mov_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fix_dz = fix_im.GetSpacing()[2]
    mov_dz = mov_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fix_zpos = fix_origin_z + fix_ind*fix_dz
    mov_zpos = mov_origin_z + mov_ind*mov_dz

    # Get the number of contours for this index:
    mov_Ncontours = len(mov_contour_pts[mov_ind])

    print(f'\nFixed slice  {fix_ind}/{fix_Nslices}',
          f'at z = {round(fix_zpos, 2)} mm')
    print(f'Moving slice {mov_ind}/{mov_Nslices}',
          f'at z = {round(mov_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fix_zpos, 2)} mm')
    #print(f'Moving depth position {round(mov_zpos, 2)} mm')
    
    #print(f'\nFixed no. of contours  = {fix_Ncontours}')
    print(f'\nMoving no. of contours = {mov_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 2, figsize=(10,8))

    # draw the fixed image in the first subplot
    plt.subplot(1, 2, 1, aspect='equal')
    plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fix_npa[fix_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')
    


    # draw the moving image in the second subplot
    plt.subplot(1, 2, 2, aspect='equal')
    plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(mov_npa[mov_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # Plot contours for fixed image if Npts is not empty:
    if mov_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[fix_ind][0]:
        for x, y, z in mov_contour_pts[mov_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(mov_npa[mov_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
        
        
        
#def display_images_and_reg_result(fix_ind, mov_ind, reg_ind, 
#                                  fix_npa, mov_npa, reg_npa, 
#                                  fix_title, mov_title, reg_title):
#
#    # create a figure with two subplots and the specified size
#    plt.subplots(1,3,figsize=(9,6))
#
#    # draw the fixed image in the first subplot
#    plt.subplot(1,3,1)
#    plt.imshow(fix_npa[fix_ind,:,:],cmap=plt.cm.Greys_r);
#    plt.title(fix_title)
#    plt.axis('off')
#
#    # draw the moving image in the second subplot
#    plt.subplot(1,3,2)
#    plt.imshow(mov_npa[mov_ind,:,:],cmap=plt.cm.Greys_r);
#    plt.title(mov_title)
#    plt.axis('off')
#    
#    # draw the registered image in the third subplot
#    plt.subplot(1,3,3)
#    plt.imshow(reg_npa[reg_image_z,:,:],cmap=plt.cm.Greys_r);
#    plt.title(reg_title)
#    plt.axis('off')
#    
#    return


def display_sitk_images_and_reg_result(fix_ind, mov_ind, reg_ind, 
                                       fix_npa, mov_npa, reg_npa, 
                                       fix_title, mov_title, reg_title):

    # create a figure with two subplots and the specified size
    plt.subplots(1,3,figsize=(14,8))

    # draw the fixed image in the first subplot
    plt.subplot(1,3,1)
    plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,3,2)
    plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1,3,3)
    plt.imshow(reg_npa[reg_ind,:,:], cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    return


def display_ants_images_and_reg_result(fix_ind, mov_ind, reg_ind, 
                                       fix_npa, mov_npa, reg_npa, 
                                       fix_title, mov_title, reg_title):
    
    # create a figure with two subplots and the specified size
    plt.subplots(1,3,figsize=(15,6))

    # draw the fixed image in the first subplot
    plt.subplot(1,3,1)
    plt.imshow(np.fliplr(np.rot90(fix_npa[:,:,fix_ind], 3)), cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1,3,2)
    plt.imshow(np.fliplr(np.rot90(mov_npa[:,:,mov_ind], 3)), cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1,3,3)
    plt.imshow(np.fliplr(np.rot90(reg_npa[:,:,reg_ind], 3)), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    return


def display_sitk_images_and_reg_result_with_contours(fix_ind, mov_ind, reg_ind, 
                                                     fix_npa, mov_npa, reg_npa, 
                                                     fix_im, mov_im, reg_im,
                                                     mov_contour_pts,
                                                     fix_title, mov_title, reg_title):
    
    # Get the number of slices in the images:
    fix_Nslices = fix_im.GetSize()[2]
    mov_Nslices = mov_im.GetSize()[2]
    reg_Nslices = reg_im.GetSize()[2]

    # Get the z positions of the origins:
    fix_origin_z = fix_im.GetOrigin()[2]
    mov_origin_z = mov_im.GetOrigin()[2]
    reg_origin_z = reg_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fix_dz = fix_im.GetSpacing()[2]
    mov_dz = mov_im.GetSpacing()[2]
    reg_dz = reg_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fix_zpos = fix_origin_z + fix_ind*fix_dz
    mov_zpos = mov_origin_z + mov_ind*mov_dz
    reg_zpos = reg_origin_z + reg_ind*reg_dz

    # Get the number of contours for this index:
    mov_Ncontours = len(mov_contour_pts[mov_ind])

    print(f'\nFixed slice  {fix_ind}/{fix_Nslices}',
          f'at z = {round(fix_zpos, 2)} mm')
    print(f'Moving slice {mov_ind}/{mov_Nslices}',
          f'at z = {round(mov_zpos, 2)} mm')
    print(f'Reg slice    {reg_ind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fix_zpos, 2)} mm')
    #print(f'Moving depth position {round(mov_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed no. of contours  = {fix_Ncontours}')
    print(f'\nMoving no. of contours = {mov_Ncontours}')
    #print(f'Reg no. of contours    = {fix_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,8))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fix_npa[fix_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(mov_npa[mov_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # Plot contours for fixed image if Npts is not empty:
    if mov_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[fix_ind][0]:
        for x, y, z in mov_contour_pts[mov_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(mov_npa[mov_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(reg_npa[reg_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    if False:
        # Using reg_inds and fix_contour_pts, plot contours belonging to the 
        # fixed image for the registered image if Npts is not empty:
        if mov_Ncontours:
            # Unpack tuple and store each x,y tuple in arrays X and Y:
            X = []
            Y = []
            #for x, y, z in fix_contour_pts[reg_ind][0]:
            for x, y, z in mov_contour_pts[reg_ind]:
                X.append(x)
                Y.append(y)
    
            # Flip Y pixel coordinates?:
            if False:
                N_r, N_c = np.shape(reg_npa[reg_ind,:,:])
    
                Y = N_c - np.array(Y)
    
            plt.plot(X, Y, linewidth=0.5, c='red');
        

def display_ants_images_and_reg_result_with_contours(fix_ind, mov_ind, reg_ind, 
                                                     fix_npa, mov_npa, reg_npa, 
                                                     fix_im, mov_im, reg_im,
                                                     fix_contour_pts,
                                                     fix_title, mov_title, reg_title):
    
    # Get the number of slices in the images:
    fix_Nslices = fix_im.shape[2]
    mov_Nslices = mov_im.shape[2]
    reg_Nslices = reg_im.shape[2]

    # Get the z positions of the origins:
    fix_origin_z = fix_im.origin[2]
    mov_origin_z = mov_im.origin[2]
    reg_origin_z = reg_im.origin[2]

    # Get the pixel spacings along the z axis:
    fix_dz = fix_im.spacing[2]
    mov_dz = mov_im.spacing[2]
    reg_dz = reg_im.spacing[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fix_zpos = fix_origin_z + fix_ind*fix_dz
    mov_zpos = mov_origin_z + mov_ind*mov_dz
    reg_zpos = reg_origin_z + reg_ind*reg_dz

    # Get the number of contours for this index:
    fix_Ncontours = len(fix_contour_pts[fix_ind])

    print(f'\nFixed slice  {fix_ind}/{fix_Nslices}',
          f'at z = {round(fix_zpos, 2)} mm')
    print(f'Moving slice {mov_ind}/{mov_Nslices}',
          f'at z = {round(mov_zpos, 2)} mm')
    print(f'Reg slice    {reg_ind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fix_zpos, 2)} mm')
    #print(f'Moving depth position {round(mov_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    print(f'\nFixed no. of contours  = {fix_Ncontours}')
    #print(f'Moving no. of contours = {fix_Ncontours}')
    #print(f'Reg no. of contours    = {fix_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,6))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(np.fliplr(np.rot90(fix_npa[:,:,fix_ind], 3)), cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fix_npa[fix_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')
    
    # Plot contours for fixed image if Npts is not empty:
    if fix_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[fix_ind][0]:
        for x, y, z in fix_contour_pts[fix_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(fix_npa[:,:,fix_ind])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(np.fliplr(np.rot90(mov_npa[:,:,mov_ind], 3)), cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(mov_npa[mov_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(np.fliplr(np.rot90(reg_npa[:,:,reg_ind], 3)), cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    # Using reg_inds and fix_contour_pts, plot contours belonging to the 
    # fixed image for the registered image if Npts is not empty:
    if fix_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[reg_ind][0]:
        for x, y, z in fix_contour_pts[reg_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(reg_npa[:,:,reg_ind])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
        
    return
        
        
""" Similar to display_images_and_contours() but include MovingContourPts: """
def display_images_and_reg_result_with_all_contours(fix_ind, mov_ind, reg_ind, 
                                                    fix_npa, mov_npa, reg_npa, 
                                                    fix_im, mov_im, reg_im,
                                                    fix_contour_pts, mov_contour_pts,
                                                    fix_title, mov_title, reg_title):
    
    # Get the number of slices in the images:
    fix_Nslices = fix_im.GetSize()[2]
    mov_Nslices = mov_im.GetSize()[2]
    reg_Nslices = reg_im.GetSize()[2]

    # Get the z positions of the origins:
    fix_origin_z = fix_im.GetOrigin()[2]
    mov_origin_z = mov_im.GetOrigin()[2]
    reg_origin_z = reg_im.GetOrigin()[2]

    # Get the pixel spacings along the z axis:
    fix_dz = fix_im.GetSpacing()[2]
    mov_dz = mov_im.GetSpacing()[2]
    reg_dz = reg_im.GetSpacing()[2]

    # Use the z positions of the origins, the pixel spacings along the z axis,
    # and the indeces to get the z positions of the slices in the stacks:
    fix_zpos = fix_origin_z + fix_ind*fix_dz
    mov_zpos = mov_origin_z + mov_ind*mov_dz
    reg_zpos = reg_origin_z + reg_ind*reg_dz

    # Get the number of contours for this index:
    fix_Ncontours = len(fix_contour_pts[fix_ind])
    mov_Ncontours = len(mov_contour_pts[mov_ind])

    print(f'\nFixed slice  {fix_ind}/{fix_Nslices}',
          f'at z = {round(fix_zpos, 2)} mm')
    print(f'Moving slice {mov_ind}/{mov_Nslices}',
          f'at z = {round(mov_zpos, 2)} mm')
    print(f'Reg slice    {reg_ind}/{reg_Nslices}',
          f'at z = {round(reg_zpos, 2)} mm')
    
    #print(f'\nFixed depth position  {round(fix_zpos, 2)} mm')
    #print(f'Moving depth position {round(mov_zpos, 2)} mm')
    #print(f'Reg depth position    {round(reg_zpos, 2)} mm')
    
    print(f'\nFixed no. of contours  = {fix_Ncontours}')
    print(f'Moving no. of contours = {mov_Ncontours}')
    #print(f'Reg no. of contours    = {fix_Ncontours}')


    # create a figure with two subplots and the specified size
    plt.subplots(1, 3, figsize=(14,5))

    # draw the fixed image in the first subplot
    plt.subplot(1, 3, 1, aspect='equal')
    plt.imshow(fix_npa[fix_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(fix_npa[fix_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(fix_title)
    plt.axis('off')
    
    # Plot contours for fixed image if fix_Ncontours is not 0:
    if fix_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[fix_ind][0]:
        for x, y, z in fix_contour_pts[fix_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(fix_npa[fix_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');

    # draw the moving image in the second subplot
    plt.subplot(1, 3, 2, aspect='equal')
    plt.imshow(mov_npa[mov_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(mov_npa[mov_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(mov_title)
    plt.axis('off')
    
    # Plot contours for moving image if mov_Ncontours is not 0:
    if mov_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        for x, y, z in mov_contour_pts[mov_ind]:
            X.append(x)
            Y.append(y)

        # Flip Y pixel coordinates?:
        if False:
            N_r, N_c = np.shape(mov_npa[mov_ind,:,:])

            Y = N_c - np.array(Y)

        plt.plot(X, Y, linewidth=0.5, c='red');
    
    # draw the registered image in the third subplot
    plt.subplot(1, 3, 3, aspect='equal')
    plt.imshow(reg_npa[reg_ind,:,:], cmap=plt.cm.Greys_r);
    #plt.pcolormesh(np.flipud(reg_npa[reg_ind,:,:]), cmap=plt.cm.Greys_r);
    plt.title(reg_title)
    plt.axis('off')
    
    # Using reg_inds and fix_contour_pts, plot contours belonging to the 
    # fixed image for the registered image if fix_Ncontours is not 0:
    if fix_Ncontours:
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        #for x, y, z in fix_contour_pts[reg_ind][0]:
        for x, y, z in fix_contour_pts[reg_ind]:
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
    
    
    
def display_all_sitk_images(fix_im, mov_im, export_fig):
    
    fix_dz = fix_im.GetSpacing()[2]
    mov_dz = mov_im.GetSpacing()[2]
    
    fix_z0 = fix_im.GetOrigin()[2]
    mov_z0 = mov_im.GetOrigin()[2]

    fix_npa = sitk.GetArrayFromImage(fix_im)
    mov_npa = sitk.GetArrayFromImage(mov_im)
    
    fix_N = len(fix_npa)
    mov_N = len(mov_npa)
    
    Nrows = max(fix_N, mov_N)
    Ncols = 2

    # create a figure with two subplots and the specified size
    plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 5*Nrows))

    i = 0 # sub-plot number
    
    for s in range(Nrows):
        #print(f's = {s}')
        
        # Compute the z position:
        fix_z = fix_z0 + fix_dz*(s+1)
        mov_z = mov_z0 + mov_dz*(s+1)
        
        i = i + 1 # increment sub-plot number
        
        # Draw the fixed image:
        if s <= fix_N-1:
            plt.subplot(Nrows, Ncols, i)
            plt.imshow(fix_npa[s,:,:], cmap=plt.cm.Greys_r);
            #plt.title(f'Fixed slice {s+1}/{fix_N} (z = {round(fix_z,1)} mm)')
            plt.title(f'Fixed slice {s+1}/{fix_N}')
            plt.axis('off')
    
        i = i + 1 # increment sub-plot number
        
        # Draw the moving image:
        if s <= mov_N-1:
            plt.subplot(Nrows, Ncols, i)
            plt.imshow(mov_npa[s,:,:], cmap=plt.cm.Greys_r);
            #plt.title(f'Moving slice {s+1}/{mov_N} (z = {round(mov_z,1)} mm)')
            plt.title(f'Moving slice {s+1}/{mov_N}')
            plt.axis('off')
            
    if export_fig:
        plt.savefig('Images_sitk.png', bbox_inches='tight')



def display_all_ants_images(fix_im, mov_im, export_fig):
    
    fix_dz = fix_im.spacing[2]
    mov_dz = mov_im.spacing[2]
    
    fix_z0 = fix_im.origin[2]
    mov_z0 = mov_im.origin[2]

    fix_npa = fix_im.numpy()
    mov_npa = mov_im.numpy()
    
    #fix_N = len(fix_npa) # = 192
    #mov_N = len(mov_npa) # = 212
    fix_N = fix_im.shape[2]
    mov_N = mov_im.shape[2]
    
    Nrows = max(fix_N, mov_N)
    Ncols = 2

    #print(f'fix_N = {fix_N}')
    #print(f'mov_N = {mov_N}')
    
    # create a figure with two subplots and the specified size
    plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 5*Nrows))

    i = 0 # sub-plot number
    
    for s in range(Nrows):
        #print(f's = {s}')
        
        # Compute the z position:
        fix_z = fix_z0 + fix_dz*(s+1)
        mov_z = mov_z0 + mov_dz*(s+1)
        
        i = i + 1 # increment sub-plot number
        
        # Draw the fixed image:
        if s <= fix_N-1:
            plt.subplot(Nrows, Ncols, i)
            plt.imshow(np.flipud(np.rot90(fix_npa[:,:,s])), cmap=plt.cm.Greys_r);
            #plt.title(f'Fixed slice {s+1}/{fix_N} (z = {round(fix_z,1)} mm)')
            plt.title(f'Fixed slice {s+1}/{fix_N}')
            plt.axis('off')
    
        i = i + 1 # increment sub-plot number
        
        # Draw the moving image:
        if s <= mov_N-1:
            plt.subplot(Nrows, Ncols, i)
            plt.imshow(np.flipud(np.rot90(mov_npa[:,:,s])), cmap=plt.cm.Greys_r);
            #plt.title(f'Moving slice {s+1}/{mov_N} (z = {round(mov_z,1)} mm)')
            plt.title(f'Moving slice {s+1}/{mov_N}')
            plt.axis('off')

    if export_fig:
        plt.savefig('Images_ants.png', bbox_inches='tight')