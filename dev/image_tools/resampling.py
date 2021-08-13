# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 20:24:22 2021

@author: ctorti
"""


""" Resampling functions for SimpleITK images. """

import SimpleITK as sitk

from image_tools.attrs_info import get_im_info
from image_tools.operations import change_im_dtype, find_thresh, binarise_im
from image_tools.operations import gaussian_blur_im#, recursive_gaussian_blur_im
from conversion_tools.pixarrs_ims import im_to_pixarr
from general_tools.fiducials import get_landmark_tx
from plotting_tools.plotting import plot_two_ims


def resample_im(Im, RefIm, sitkTx=sitk.Transform(), Interp='Linear', 
                LogToConsole=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------   
    Im : SimpleITK Image
        The 3D image to be resampled.
    RefIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
    sitkTx : SimpleITK Transform, optional (sitk.Transform() by default)
        The SimpleITK Transform to be used. The identity transform is default.
    Interp : str, optional ('Linear' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    ResIm : SimpleITK image
        The resampled 3D image.
    
    Note
    ----
    While a linear (or BSpline) interpolator is appropriate for intensity 
    images, only a NearestNeighbor interpolator is appropriate for binary 
    images (e.g. segmentations) so that no new labels are introduced.
    """
    
    # Define which interpolator to use:
    if Interp == 'NearestNeighbor':    
        sitkInterp = sitk.sitkNearestNeighbor
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'LabelGaussian':
        sitkInterp = sitk.sitkLabelGaussian
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Linear':
        sitkInterp = sitk.sitkLinear
        #sitkPixType = sitk.sitkFloat32
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Bspline':
        sitkInterp = sitk.sitkBSpline
        sitkPixType = sitk.sitkFloat32
        
    else:
        msg = '"Interp" must be "Linear", "BSpline" or "LabelGaussian".'
        
        raise Exception(msg)
    
    #print('\nUsing', Interpolation, 'interp\n')
    
    Resampler = sitk.ResampleImageFilter()
    Resampler.SetTransform(sitkTx)
    Resampler.SetInterpolator(sitkInterp)
    
    Resampler.SetReferenceImage(RefIm)
    Resampler.SetOutputPixelType(sitkPixType)
    Resampler.SetDefaultPixelValue(0)
    
    ResIm = Resampler.Execute(Im)
    
    #if LogToConsole:
    #    ResTx = Resampler.GetTransform()
    #    print('\nResampling transform:\n', ResTx)

    return ResIm

def resample_labim_OLD(LabIm, F2Sinds, SrcIm, TrgIm, Interp, 
                   PreResVariance=(1,1,1), ApplyPostResBlur=True, 
                   PostResVariance=(1,1,1), LogToConsole=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------  
    LabIm : SimpleITK image
        The 3D label image to be resampled.
    F2Sinds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        LabIm.
    SrcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of LabIm. 
    TrgIm : SimpleITK image
        The Target 3D image whose gridspace the LabIm will be resampled to.
    Interp : str
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    PreResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    ResLabIm : SimpleITK images
        The resampled 3D label image.
    ResPixArr : Numpy arrays
        The resampled pixel array (converted from ResLabIm).
    ResF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in ResPixArr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, F2Sinds = []) even if 
    there were contours (i.e. if SrcC2Sinds != []). If this is the case try
    suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.
    
    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
    
    Link: https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print('\n\n', '-'*120)
        
    if not Interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f'The chosen interpolation, {Interp}, is not one of the '\
              + 'accepted inputs: \'NearestNeighbor\', \'LabelGaussian\', or '\
              + '\'BlurThenLinear\'.'
        
        raise Exception(msg)
    
    # Store the interpolation set as metadata:
    LabIm.SetMetaData("ResInterpSet", Interp)
    
    #SrcSpacings = SrcIm.GetSpacing()
    #TrgSpacings = TrgIm.GetSpacing()
    
    #SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
    #                                                      TrgSpacings))
    
    #SrcSpacings = np.array(SrcIm.GetSpacing())
    #TrgSpacings = np.array(TrgIm.GetSpacing())
        
    #SpacingsRatio = np.divide(TrgSpacings, SrcSpacings)
    
    #VolumeRatio = np.prod(SpacingsRatio)
    
    #""" 
    #The FWHM of the Gaussian is:
    #    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
    #    
    #Due to the Nyquist-Shannon sampling theorem, the FWHM should be at
    #least 2x the voxel spacings of TrgIm. For symmetry it makes sense for
    #one Source voxel width to be blurred to three Target voxels widths:
    #    sigma = (3*TrgSpacings)/2.355
    #"""
    
    #if LogToConsole:
    #    print(f'\nSrcSpacings = {SrcSpacings}')
    #    print(f'TrgSpacings = {TrgSpacings}')
    #    print(f'\nSpacingsRatio = {SpacingsRatio}')
    
    if Interp in ['NearestNeighbor', 'LabelGaussian']:
        if LogToConsole:
            print('Attempting to resample LabIm using the',
                  f'{Interp} interpolator...\n')
        
        """ Attempt to resample LabIm to the Target image's grid using the
        chosen labelmap interpolator. """
        ResLabIm = resample_im(Im=LabIm, RefIm=TrgIm, Interp=Interp)
        
        if LogToConsole:
            print('Image info for ResLabIm after resampling using',
                  f'{Interp} interpolator:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
        
        if ApplyPostResBlur:
            # Gaussian blur the image.
            ResLabIm = gaussian_blur_im(Im=ResLabIm, 
                                        Variance=PostResVariance)
            
            if LogToConsole:
                print('Image info for ResLabIm after Gaussian blurring:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
            if LogToConsole:
                print('')
            
            #""" The original binary pixel array: """
            #PixArr_B, F2Sinds = im_to_pixarr(LabIm)
            #
            #""" The non-binary pixel array following blur: """
            #PixArr_NB, F2Sinds = im_to_pixarr(ResLabIm)
            #
            #Ratio = GetVoxelVolumeRatio(LabIm, ResLabIm)
            #
            #if LogToConsole:
            #    print(f'PixArr_B.shape = {PixArr_B.shape}')
            #    print(f'PixArr_NB.shape = {PixArr_NB.shape}\n')
            
            # Binarise the resampled labelmap if required.
            if len(UniqueVals) != 2 or sum(UniqueVals) != 1:
                """Find suitable threshold value that approximately preserves
                the number of pre-blurred truth values scaled by VolumeRatio: 
                """
                Thresh = find_thresh(BinaryIm=LabIm, 
                                     NonBinaryIm=ResLabIm,
                                     LogToConsole=LogToConsole)
                
                # Binary threshold the image.
                ResLabIm = binarise_im(Im=ResLabIm, 
                                       #Thresh=PostResThresh)
                                       Thresh=Thresh) 
                
                if LogToConsole:
                    print('\nImage info for ResLabIm after binary',
                          #f'thresholding at {PostResThresh}:')
                          f'thresholding at {Thresh}:')
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = get_im_info(ResLabIm, LogToConsole)
                if LogToConsole:
                    print('')
        
        #print(f'\n   ResF2Sinds = {ResF2Sinds}')
        
        """ Is ResF2Sinds empty and not expected to be? If so, try the 
        "BlurThenLinear" approach.
        
        Note:
            If F2Sinds isn't empty, ResF2Sinds shouldn't be empty either.
        
            F2Sinds will be empty if there were no segmentations/contours 
            of interest for the r^th ROI. In this case an empty ResF2Sinds 
            is acceptable. """
        
        if ResF2Sinds == []:
            print(f'There are {len(F2Sinds)} non-empty masks in the input',
                  f'labelmap image but {len(ResF2Sinds)} non-empty frames in',
                  f'the resampled labelmap image using {Interp}.',
                  'Try using the \'BlurThenLinear\' approach.\n')
            
            Interp = 'BlurThenLinear'

    if Interp == 'BlurThenLinear':
        if LogToConsole:
            print('\nResampling LabIm by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        #""" The original binary pixel array: """
        #PixArr_B, F2Sinds = im_to_pixarr(LabIm)
        
        """ Convert LabIm from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        LabIm = change_im_dtype(Image=LabIm, NewPixelType='Float32')
        
        if LogToConsole:
            print('\nImage info for LabIm after converting to float:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = get_im_info(LabIm, LogToConsole)
              
        # Gaussian blur the image.
        LabIm = gaussian_blur_im(Im=LabIm, Variance=PreResVariance)
        
        # Use the RecursiveGaussian image filter:
        #ResLabIm = recursive_gaussian_blur_im(ResLabIm, 
        #                                      Sigma=3,
        #                                      Direction=2)
        
        if LogToConsole:
            print('\nImage info for LabIm after applying Gaussian blur:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = get_im_info(LabIm, LogToConsole)
        
        if LogToConsole:
            print('\nLabIm prior to resampling:')
            print(f'   LabIm.GetSize() = {LabIm.GetSize()}')
            print(f'   LabIm.GetSpacing() = {LabIm.GetSpacing()}')
            print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
            print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
        
        # Linearly resample LabIm to the Target image's grid.
        ResLabIm = resample_im(Im=LabIm, RefIm=TrgIm, Interp='Linear')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after resampling using a',
                  'linear interpolator:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
        
        if ApplyPostResBlur:
            # Gaussian blur the image.
            ResLabIm = gaussian_blur_im(Im=ResLabIm, 
                                        Variance=PostResVariance)
        
        #""" The non-binary pixel array following blur + linear resampling: """
        #PixArr_NB, F2Sinds = im_to_pixarr(ResLabIm)
        #
        #if LogToConsole:
        #        print(f'PixArr_B.shape = {PixArr_B.shape}')
        #        print(f'PixArr_NB.shape = {PixArr_NB.shape}\n')
            
        """Find suitable threshold value that approximately preserves the
        number of pre-blurred + linear resampled truth values scaled by
        VolumeRatio: """
        Thresh = find_thresh(BinaryIm=LabIm, NonBinaryIm=ResLabIm,
                             LogToConsole=LogToConsole)
        
        # Binary threshold the image. 
        ResLabIm = binarise_im(Im=ResLabIm, 
                               #Thresh=PostResThresh)
                               Thresh=Thresh) 
        
        if LogToConsole:
            print('\nImage info for ResLabIm after binary thresholding',
                  #f'at {PostResThresh}:')
                  f'at {Thresh}:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
    
    # Ensure that ResLabIm is a 32-bit unsigned integer (PixID = 5).
    if PixID != 5: 
        if LogToConsole:
            print(f'\nResLabIm has PixelID = {PixID} ({PixIDTypeAsStr})).')
        
        # Convert ResLabIm from float to 32-bit unsigned integer. 
        ResLabIm = change_im_dtype(Image=ResLabIm, NewPixelType='UInt32')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after converting to 32-bit',
                  'unsigned int:')
        #print(f'\nThe metadata keys are:', ResLabIm.GetMetaDataKeys())
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = get_im_info(ResLabIm, LogToConsole)
    
    # Convert ResLabIm to a pixel array. 
    ResPixArr, ResF2Sinds = im_to_pixarr(ResLabIm)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    ResLabIm.SetMetaData("ResInterpUsed", Interp)
    
    if LogToConsole:
        print(f'\nThe key "ResInterpUsed" with value "{Interp}" was',
              'added to the metadata of ResLabIm. \nThe metadata keys are:',
              ResLabIm.GetMetaDataKeys())
    
    if Interp == 'BlurThenLinear':
        """ Store the threshold used as metadata: """
        ResLabIm.SetMetaData("PostResThreshUsed", f"{Thresh}")
        
        if LogToConsole:
            print(f'\nThe key "PostResThreshUsed" with value "{Thresh}" was',
                  'added to the metadata of ResLabIm. \nThe metadata keys:',
                  ResLabIm.GetMetaDataKeys())
            
    if LogToConsole:
        print('\nAfter converting ResLabIm to a pixel array:')
        print(f'ResPixArr.shape = {ResPixArr.shape}')
        print(f'ResF2Sinds = {ResF2Sinds}')
        print('-'*120)
        
    return ResLabIm, ResPixArr, ResF2Sinds

def resample_labim(LabIm, F2Sinds, SrcIm, TrgIm,  
                   sitkTx=sitk.Transform(), Interp='NearestNeighbor', 
                   ApplyPreResBlur=False, PreResVariance=(1,1,1), 
                   ApplyPostResBlur=True, PostResVariance=(1,1,1), 
                   LogToConsole=False):
    """
    Resample a 3D label image.
    
    Parameters
    ----------  
    LabIm : SimpleITK image
        The 3D label image to be resampled.
    F2Sinds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        LabIm.
    SrcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of LabIm. 
    TrgIm : SimpleITK image
        The Target 3D image whose gridspace the LabIm will be resampled to.
    sitkTx : SimpleITK Transform, optional (sitk.Transform() by default)
        The SimpleITK Transform to be used. The identity transform is default.
    Interp : str, optional ('NearestNeighbor' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    ApplyPreResBlur : bool, optional (False by default)
        If True, the pre-resampled LabIm will be Gaussian blurred. This 
        argument is only applicable if Interp is 'NearestNeighbor' or
        'LabelGaussian', since LabIm will be blurred prior to resampling if
        Interp is 'BlurThenLinear' regardless of the argument's value.
    PreResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    ResLabIm : SimpleITK Images
        The resampled 3D label image.
    ResPixArr : Numpy arrays
        The pixel array representation of ResLabIm.
    ResF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in ResPixArr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, F2Sinds = []) even if 
    there were contours (i.e. if SrcC2Sinds != []). If this is the case try
    suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.
    
    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
    
    Link: https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print('\n\n', '-'*120)
        
    if not Interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f"The chosen interpolation, {Interp}, is not one of the "\
              + "accepted arguments: 'NearestNeighbor', 'LabelGaussian', or "\
              + "'BlurThenLinear'."
        raise Exception(msg)
    
    # Store the interpolation set as metadata:
    LabIm.SetMetaData("ResInterpSet", Interp)
    
    if Interp in ['NearestNeighbor', 'LabelGaussian']:
        if LogToConsole:
            print(f'Attempting to resample LabIm using {Interp} interpolator\n')
        
        if ApplyPreResBlur:
            # Gaussian blur LabIm:
            BlurLabIm = gaussian_blur_im(Im=LabIm, Variance=PostResVariance)
            
            # Resample BlurLabIm using the chosen interpolator:
            ResLabIm = resample_im(Im=BlurLabIm, RefIm=TrgIm, sitkTx=sitkTx,
                                   Interp=Interp)
            
            msg = 'Image info for resampled blurred image:'
        else:
            # Resample LabIm using the chosen interpolator:
            ResLabIm = resample_im(Im=LabIm, RefIm=TrgIm, sitkTx=sitkTx,
                                   Interp=Interp)
            
            msg = 'Image info for resampled image:'
            
        if LogToConsole:
            print(msg)
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
        
        if ApplyPostResBlur:
            # Gaussian blur ResLabIm:
            ResLabIm = gaussian_blur_im(Im=ResLabIm, Variance=PostResVariance)
            
            if LogToConsole:
                print('Image info for blurred resampled image:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
            if LogToConsole:
                print('')
            
            # Binarise ResLabIm if required:
            if len(UniqueVals) != 2 or sum(UniqueVals) != 1:
                """
                ResLabIm is not binary. Find suitable threshold value that 
                approximately preserves the number of pre-blurred truth values
                scaled by VolumeRatio: 
                """
                Thresh = find_thresh(BinaryIm=LabIm, NonBinaryIm=ResLabIm,
                                     LogToConsole=LogToConsole)
                
                # Binary threshold ResLabIm:
                ResLabIm = binarise_im(Im=ResLabIm, Thresh=Thresh) 
                
                if LogToConsole:
                    print(f'\nImage info after binary thresholding at {Thresh}:')
                PixID, PixIDTypeAsStr, UniqueVals,\
                ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
                if LogToConsole:
                    print('')
        
        #print(f'\n   ResF2Sinds = {ResF2Sinds}')
        
        """ 
        Is ResF2Sinds empty and not expected to be? If so try the 
        "BlurThenLinear" approach.
        
        If F2Sinds isn't empty, ResF2Sinds shouldn't be empty either.
        
        F2Sinds will be empty if there were no segmentations/contours of 
        interest for the r^th ROI. In this case an empty ResF2Sinds is
        acceptable. 
        """
        
        if ResF2Sinds == []:
            print(f"There are {len(F2Sinds)} non-empty masks in the input",
                  f"label image but {len(ResF2Sinds)} non-empty frames in the",
                  f"resampled label image using {Interp}. Will Gaussian blur,",
                  "linearly resample and binarise...\n")
            
            Interp = 'BlurThenLinear'

    if Interp == 'BlurThenLinear':
        # Gaussian blur LabIm:
        BlurLabIm = gaussian_blur_im(Im=LabIm, Variance=PreResVariance)
        
        if LogToConsole:
            print('\nImage info for BlurLabIm:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = get_im_info(BlurLabIm, LogToConsole)
            print('\n\nBlurLabIm prior to resampling:')
            print(f'   BlurLabIm.GetSize() = {BlurLabIm.GetSize()}')
            print(f'   BlurLabIm.GetSpacing() = {BlurLabIm.GetSpacing()}')
            print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
            print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
        
        # Linearly resample BlurLabIm:
        ResLabIm = resample_im(Im=BlurLabIm, RefIm=TrgIm, sitkTx=sitkTx, 
                               Interp='Linear')
        
        if LogToConsole:
            print('\nImage info after resampling using linear interpolator:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
        
        if ApplyPostResBlur:
            # Gaussian blur ResLabIm:
            ResLabIm = gaussian_blur_im(Im=ResLabIm, Variance=PostResVariance)
            
        # Find suitable threshold value:
        Thresh = find_thresh(BinaryIm=LabIm, NonBinaryIm=ResLabIm,
                             LogToConsole=LogToConsole)
        
        # Binary threshold ResLabIm:
        ResLabIm = binarise_im(Im=ResLabIm, Thresh=Thresh) 
        
        if LogToConsole:
            print(f'\nImage info after binary thresholding {Thresh}:')
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
    
    # Ensure that ResLabIm is a 32-bit unsigned integer (PixID = 5):
    if PixID != 5: 
        if LogToConsole:
            print(f'\nResLabIm has PixelID = {PixID} ({PixIDTypeAsStr})).')
        
        # Convert ResLabIm from float to 32-bit unsigned integer:
        ResLabIm = change_im_dtype(Image=ResLabIm, NewPixelType='UInt32')
        
        if LogToConsole:
            print('\nImage info after converting to 32-bit unsigned int:')
        #print(f'\nThe metadata keys are:', ResLabIm.GetMetaDataKeys())
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = get_im_info(ResLabIm, LogToConsole)
        if LogToConsole:
            print('')
    
    # Convert ResLabIm to a pixel array:
    ResPixArr, ResF2Sinds = im_to_pixarr(ResLabIm)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    ResLabIm.SetMetaData("ResInterpUsed", Interp)
    
    if Interp == 'BlurThenLinear':
        # Store the threshold used as metadata:
        ResLabIm.SetMetaData("PostResThreshUsed", f"{Thresh}")
            
    if LogToConsole:
        # The number of frames before and after:
        N_before = len(F2Sinds)
        N_after = len(ResF2Sinds)
        
        print(f'\nThere were {N_before} frames in the label image')
        print(f'There are {N_after} frames in the resampled label image')
        print('After converting ResLabIm to a pixel array:')
        print(f'ResPixArr.shape = {ResPixArr.shape}')
        print(f'ResF2Sinds = {ResF2Sinds}')
        plot_two_ims(Im0=LabIm, Ind0=F2Sinds[0], 
                         PlotLabel0='Original label image', 
                         Im1=ResLabIm, Ind1=ResF2Sinds[0], 
                         PlotLabel1='Resampled label image')
        print('-'*120)
        
    return ResLabIm, ResPixArr, ResF2Sinds

def resample_labimByRoi(LabImByRoi, F2SindsByRoi, SrcIm, TrgIm, 
                        sitkTx=sitk.Transform(), Interp='NearestNeighbor', 
                        ApplyPreResBlur=False, PreResVariance=(1,1,1), 
                        ApplyPostResBlur=True, PostResVariance=(1,1,1), 
                        LogToConsole=False):
    """
    Resample a list 3D SimpleITK images representing binary label images. 
    
    Parameters
    ----------  
    LabImByRoi : list of SimpleITK Images
        A list (for each ROI) of 3D label images to be resampled.
    F2SindsByRoi : list of a list of ints
        List (for each ROI) of a list (for each frame) of slice numbers that 
        correspond to each frame in the label images.
    SrcIm : SimpleITK Image
        The Source 3D image whose gridspace matches the gridspace of the 
        label images.
    TrgIm : SimpleITK Image
        The Target 3D image whose gridspace the Source label images will be
        resampled to.
    sitkTx : SimpleITK Transform, optional (sitk.Transform() by default)
        The SimpleITK Transform to be used. The identity transform is default.
    Interp : str, optional ('NearestNeighbor' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    ApplyPreResBlur : bool, optional (False by default)
        If True, the pre-resampled LabIm will be Gaussian blurred. This 
        argument is only applicable if Interp is 'NearestNeighbor' or
        'LabelGaussian', since LabIm will be blurred prior to resampling if
        Interp is 'BlurThenLinear' regardless of the argument's value.
    PreResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
    ApplyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    PostResVariance : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        resampled labelmap image(s) is/are to be Gaussian blurred after  
        resampling.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    ResLabImByRoi : list of SimpleITK Images
        The list (for each ROI) of resampled 3D label image.
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of the pixel array representations of
        ResLabImByRoi.
    ResF2SindsByRoi : list of a list of ints
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.
    
    Note
    ----
    See Notes in resample_labim. 
    """
        
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of resample_labimByRoi():')
        print('\n\n', '-'*120)
        
    ResLabImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabImByRoi)):
        if LogToConsole:
            print(f'   Resampling of LabImByRoi[{r}]...')
        
        ResLabIm, ResPixArr,ResF2Sinds\
            = resample_labim(LabIm=LabImByRoi[r], F2Sinds=F2SindsByRoi[r], 
                             SrcIm=SrcIm, TrgIm=TrgIm, sitkTx=sitkTx,
                             Interp=Interp, ApplyPreResBlur=ApplyPreResBlur, 
                             PreResVariance=PreResVariance, 
                             ApplyPostResBlur=ApplyPostResBlur, 
                             PostResVariance=PostResVariance, 
                             LogToConsole=LogToConsole)
        
        ResLabImByRoi.append(ResLabIm)
        ResPixArrByRoi.append(ResPixArr)
        ResF2SindsByRoi.append(ResF2Sinds)
        
    if LogToConsole:
            print('-'*120)
        
    return ResLabImByRoi, ResPixArrByRoi, ResF2SindsByRoi

def align_ims_using_landmarks(Transform, FixFidsFpath, MovFidsFpath, 
                              FixIm, MovIm, FlipK=False, Buffer=3, 
                              LogToConsole=False):
    """
    Return the landmark-based aligned image and transform based on fiducials.
    
    Parameters
    ----------
    FixIm : SimpleITK Image 
        The 3D image that MovIm will be aligned to.
    MovIm : SimpleITK Image
        The 3D image that will be aligned to FixIm.
    FixFidsFpath : str, optional ('fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension.
    MovFidsFpath : str, optional ('moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension.
    FlipK : bool, optional (False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ.
    Buffer : int, optional (3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    Transform : str, optional ('affine' by default)
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    AlignedIm : SimpleITK Image
        The 3D aligned image.
    LandmarkTx : SimpleITK Transform
        The transform used to perform landmark-based alignment.
        
    Note
    ----
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(FixIm.GetDimension())
        sitk.BSplineTransform(FixIm.GetDimension())
        
    NumControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    LandmarkTx, FixPts, MovPts\
        = get_landmark_tx(Transform=Transform, 
                          FixFidsFpath=FixFidsFpath, 
                          MovFidsFpath=MovFidsFpath, 
                          FixIm=FixIm, MovIm=MovIm, FlipK=FlipK,
                          Buffer=Buffer, LogToConsole=LogToConsole)
    
    #print(f'LandmarkTx.GetName() = {LandmarkTx.GetName()}\n')
    
    AlignedIm = sitk.Resample(MovIm, FixIm, LandmarkTx, sitk.sitkLinear, 0.0, 
                              MovIm.GetPixelID())
    
    return AlignedIm, LandmarkTx