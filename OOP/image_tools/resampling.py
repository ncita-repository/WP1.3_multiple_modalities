# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 20:24:22 2021

@author: ctorti
"""


""" Resampling functions for SimpleITK images. """

from importlib import reload
import image_tools.operations
reload(image_tools.operations)
import image_tools.attrs_info
reload(image_tools.attrs_info)

import SimpleITK as sitk

from image_tools.attrs_info import get_im_info
from image_tools.operations import (
    change_im_dtype, find_thresh, binarise_im, gaussian_blur_im#, recursive_gaussian_blur_im
    )
from conversion_tools.pixarrs_ims import im_to_pixarr
from general_tools.fiducials import get_landmark_tx
from plotting_tools.general import plot_two_ims


def resample_im(im, refIm, sitkTx=sitk.Transform(3, sitk.sitkIdentity),
                #sitkTx=sitk.Transform(), 
                interp='Linear', p2c=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------   
    im : SimpleITK Image
        The 3D image to be resampled.
    refIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
    sitkTx : SimpleITK Transform, optional
        The SimpleITK Transform to be used. The default value is
        sitk.Transform(3, sitk.sitkIdentity) (identity transform).
    interp : str, optional
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
        The default value is 'Linear'.
    p2c : bool, optional
        Denotes whether some results will be logged to the console. The default
        value is False.
        
    Returns
    -------
    resIm : SimpleITK image
        The resampled 3D image.
    
    Note
    ----
    While a linear (or BSpline) interpolator is appropriate for intensity 
    images, only a NearestNeighbor interpolator is appropriate for binary 
    images (e.g. segmentations) so that no new labels are introduced.
    """
    
    # Define which interpolator to use:
    if interp == 'NearestNeighbor':    
        sitkInterp = sitk.sitkNearestNeighbor
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif interp == 'LabelGaussian':
        sitkInterp = sitk.sitkLabelGaussian
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif interp == 'Linear':
        sitkInterp = sitk.sitkLinear
        sitkPixType = sitk.sitkFloat32 # 20/08/21
        #sitkPixType = sitk.sitkUInt32 # 20/08/21
        
    elif interp == 'Bspline':
        sitkInterp = sitk.sitkBSpline
        sitkPixType = sitk.sitkFloat32
        
    else:
        msg = '"interp" must be "Linear", "BSpline" or "LabelGaussian".'
        
        raise Exception(msg)
    
    #print('\nUsing', Interpolation, 'interp\n')
    
    resampler = sitk.ResampleImageFilter()
    resampler.SetTransform(sitkTx)
    resampler.SetInterpolator(sitkInterp)
    
    resampler.SetReferenceImage(refIm)
    resampler.SetOutputPixelType(sitkPixType)
    resampler.SetDefaultPixelValue(0)
    
    resIm = resampler.Execute(im)
    
    #if p2c:
    #    resTx = resampler.GetTransform()
    #    print('\nResampling transform:\n', resTx)

    return resIm

def resample_labim_OLD(labim, f2sInds, srcIm, trgIm, interp, 
                   preResVar=(1,1,1), applyPostResBlur=True, 
                   postResVar=(1,1,1), p2c=False):
    """
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------  
    labim : SimpleITK image
        The 3D label image to be resampled.
    f2sInds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        labim.
    srcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of labim. 
    trgIm : SimpleITK image
        The Target 3D image whose gridspace the labim will be resampled to.
    interp : str
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    preResVar : tuple of floats, optional ((1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
    applyPostResBlur : bool, optional (True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    postResVar : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    resLabim : SimpleITK images
        The resampled 3D label image.
    resPixarr : Numpy arrays
        The resampled pixel array (converted from resLabim).
    resF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in resPixarr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, f2sInds = []) even if 
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
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print('\n\n', '-'*120)
        
    if not interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f'The chosen interpolation, {interp}, is not one of the '\
              + 'accepted inputs: \'NearestNeighbor\', \'LabelGaussian\', or '\
              + '\'BlurThenLinear\'.'
        
        raise Exception(msg)
    
    # Store the interpolation set as metadata:
    labim.SetMetaData("resInterpSet", interp)
    
    #srcSpacings = srcIm.GetSpacing()
    #trgSpacings = trgIm.GetSpacing()
    
    #spacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(srcSpacings, 
    #                                                      trgSpacings))
    
    #srcSpacings = np.array(srcIm.GetSpacing())
    #trgSpacings = np.array(trgIm.GetSpacing())
        
    #spacingsRatio = np.divide(trgSpacings, srcSpacings)
    
    #volumeRatio = np.prod(spacingsRatio)
    
    #""" 
    #The FWHM of the Gaussian is:
    #    fwhm = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
    #    
    #Due to the Nyquist-Shannon sampling theorem, the FWHM should be at
    #least 2x the voxel spacings of trgIm. For symmetry it makes sense for
    #one Source voxel width to be blurred to three Target voxels widths:
    #    sigma = (3*trgSpacings)/2.355
    #"""
    
    #if p2c:
    #    print(f'\nsrcSpacings = {srcSpacings}')
    #    print(f'trgSpacings = {trgSpacings}')
    #    print(f'\nspacingsRatio = {spacingsRatio}')
    
    if interp in ['NearestNeighbor', 'LabelGaussian']:
        if p2c:
            print('Attempting to resample labim using the',
                  f'{interp} interpolator...\n')
        
        """ Attempt to resample labim to the Target image's grid using the
        chosen labelmap interpolator. """
        resLabim = resample_im(im=labim, refIm=trgIm, interp=interp)
        
        if p2c:
            print('Image info for resLabim after resampling using',
                  f'{interp} interpolator:')
        pixID, pixIDTypeAsStr, uniqueVals,\
        resF2Sinds = get_im_info(resLabim, p2c)
        if p2c:
            print('')
        
        if applyPostResBlur:
            # Gaussian blur the image.
            resLabim = gaussian_blur_im(im=resLabim, var=postResVar)
            
            if p2c:
                print('Image info for resLabim after Gaussian blurring:')
            pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
                resLabim, p2c
                )
            if p2c:
                print('')
            
            #""" The original binary pixel array: """
            #PixArr_B, f2sInds = im_to_pixarr(labim)
            #
            #""" The non-binary pixel array following blur: """
            #pixArr_NB, f2sInds = im_to_pixarr(resLabim)
            #
            #ratio = GetVoxelVolumeRatio(labim, resLabim)
            #
            #if p2c:
            #    print(f'pixArr_B.shape = {pixArr_B.shape}')
            #    print(f'pixArr_NB.shape = {pixArr_NB.shape}\n')
            
            # Binarise the resampled labelmap if required.
            if len(uniqueVals) != 2 or sum(uniqueVals) != 1:
                """Find suitable threshold value that approximately preserves
                the number of pre-blurred truth values scaled by volumeRatio: 
                """
                thresh = find_thresh(
                    binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c
                    )
                
                # Binary threshold the image.
                resLabim = binarise_im(
                    im=resLabim, #thresh=postResThresh)
                    thresh=thresh) 
                
                if p2c:
                    print('\nImage info for resLabim after binary',
                          #f'thresholding at {postResThresh}:')
                          f'thresholding at {thresh}:')
                pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                    resLabim, p2c
                    )
                if p2c:
                    print('')
        
        #print(f'\n   resF2Sinds = {resF2Sinds}')
        
        """ Is resF2Sinds empty and not expected to be? If so, try the 
        "BlurThenLinear" approach.
        
        Note:
            If f2sInds isn't empty, resF2Sinds shouldn't be empty either.
        
            f2sInds will be empty if there were no segmentations/contours 
            of interest for the r^th ROI. In this case an empty resF2Sinds 
            is acceptable. """
        
        if resF2Sinds == []:
            print(f'There are {len(f2sInds)} non-empty masks in the input',
                  f'labelmap image but {len(resF2Sinds)} non-empty frames in',
                  f'the resampled labelmap image using {interp}.',
                  'Try using the \'BlurThenLinear\' approach.\n')
            
            interp = 'BlurThenLinear'

    if interp == 'BlurThenLinear':
        if p2c:
            print('\nResampling labim by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        #""" The original binary pixel array: """
        #pixArr_B, f2sInds = im_to_pixarr(labim)
        
        """ Convert labim from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        labim = change_im_dtype(im=labim, newPixType='Float32')
        
        if p2c:
            print('\nImage info for labim after converting to float:')
            pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                labim, p2c
                )
              
        # Gaussian blur the image.
        labim = gaussian_blur_im(im=labim, variance=preResVar)
        
        # Use the RecursiveGaussian image filter:
        #resLabim = recursive_gaussian_blur_im(
        #    resLabim, sigma=3, direction=2
        #    )
        
        if p2c:
            print('\nImage info for labim after applying Gaussian blur:')
            pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                labim, p2c
                )
        
        if p2c:
            print('\nlabim prior to resampling:')
            print(f'   labim.GetSize() = {labim.GetSize()}')
            print(f'   labim.GetSpacing() = {labim.GetSpacing()}')
            print(f'   trgIm.GetSize() = {trgIm.GetSize()}')
            print(f'   trgIm.GetSpacing() = {trgIm.GetSpacing()}')
        
        # Linearly resample labim to the Target image's grid.
        resLabim = resample_im(im=labim, refIm=trgIm, interp='Linear')
        
        if p2c:
            print('\nImage info for resLabim after resampling using a',
                  'linear interpolator:')
            
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
        
        if applyPostResBlur:
            # Gaussian blur the image.
            resLabim = gaussian_blur_im(
                im=resLabim, var=postResVar
                )
        
        #""" The non-binary pixel array following blur + linear resampling: """
        #pixArr_NB, f2sInds = im_to_pixarr(resLabim)
        #
        #if p2c:
        #        print(f'pixArr_B.shape = {pixArr_B.shape}')
        #        print(f'pixArr_NB.shape = {pixArr_NB.shape}\n')
            
        """Find suitable threshold value that approximately preserves the
        number of pre-blurred + linear resampled truth values scaled by
        volumeRatio: """
        thresh = find_thresh(
            binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c
            )
        
        # Binary threshold the image. 
        resLabim = binarise_im(
            im=resLabim, #thresh=postResThresh)
            thresh=thresh
            ) 
        
        if p2c:
            print('\nImage info for resLabim after binary thresholding',
                  #f'at {postResThresh}:')
                  f'at {thresh}:')
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
    
    # Ensure that resLabim is a 32-bit unsigned integer (pixID = 5).
    if pixID != 5: 
        if p2c:
            print(f'\nresLabim has PixelID = {pixID} ({pixIDTypeAsStr})).')
        
        # Convert resLabim from float to 32-bit unsigned integer. 
        resLabim = change_im_dtype(im=resLabim, newPixType='UInt32')
        
        if p2c:
            print('\nImage info for resLabim after converting to 32-bit',
                  'unsigned int:')
        #print(f'\nThe metadata keys are:', resLabim.GetMetaDataKeys())
        pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
            resLabim, p2c
            )
    
    # Convert resLabim to a pixel array. 
    resPixarr, resF2Sinds = im_to_pixarr(resLabim)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    resLabim.SetMetaData("resInterpUsed", interp)
    
    if p2c:
        print(f'\nThe key "resInterpUsed" with value "{interp}" was',
              'added to the metadata of resLabim. \nThe metadata keys are:',
              resLabim.GetMetaDataKeys())
    
    if interp == 'BlurThenLinear':
        """ Store the threshold used as metadata: """
        resLabim.SetMetaData("postResThreshUsed", f"{thresh}")
        
        if p2c:
            print(f'\nThe key "postResThreshUsed" with value "{thresh}" was',
                  'added to the metadata of resLabim. \nThe metadata keys:',
                  resLabim.GetMetaDataKeys())
            
    if p2c:
        print('\nAfter converting resLabim to a pixel array:')
        print(f'resPixarr.shape = {resPixarr.shape}')
        print(f'resF2Sinds = {resF2Sinds}')
        print('-'*120)
        
    return resLabim, resPixarr, resF2Sinds

def resample_labim(
        labim, f2sInds, im, refIm, sitkTx=sitk.Transform(3, sitk.sitkIdentity),
        #sitkTx=sitk.Transform(), 
        interp='NearestNeighbor', applyPreResBlur=False, preResVar=(1,1,1), 
        applyPostResBlur=True, postResVar=(1,1,1), p2c=False
        ):
    """
    Resample a 3D label image.
    
    Parameters
    ----------  
    labim : SimpleITK image
        The 3D label image to be resampled.
    f2sInds : list of ints
        List (for each frame) of slice numbers that correspond to each frame in 
        labim.
    im : SimpleITK image
        The 3D image whose gridspace matches the gridspace of labim. 
    refIm : SimpleITK image
        The 3D image whose gridspace the labim will be resampled to.
    sitkTx : SimpleITK Transform, optional
        The SimpleITK Transform to be used. The default value is 
        sitk.Transform(3, sitk.sitkIdentity) (identity transform).
    interp : str, optional
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
        The default value is 'NearestNeighbor'.
    applyPreResBlur : bool, optional
        If True, the pre-resampled labim will be Gaussian blurred. This 
        argument is only applicable if interp is 'NearestNeighbor' or
        'LabelGaussian', since labim will be blurred prior to resampling if
        interp is 'BlurThenLinear' regardless of the argument's value. The
        default value is False.
    preResVar : tuple of floats, optional
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling. The
        default value is (1,1,1).
    applyPostResBlur : bool, optional
        If True, the post-resampled labelmap image will be Gaussian blurred.
        The default value is True.
    postResVar : tuple of floats, optional
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling. The
        default value is (1,1,1).
    p2c : bool, optional
        Denotes whether some results will be logged to the console. The 
        default value is False.
    
    Returns
    -------
    resLabim : SimpleITK Images
        The resampled 3D label image.
    resPixarr : Numpy arrays
        The pixel array representation of resLabim.
    resF2Sinds : list of ints
        The list (for each frame) of the frame-to-slice indices in resPixarr.
    
    Note
    ----
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (label) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled label,
    and hence an empty list of frame-to-slice indices, f2sInds = []) even if 
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
    
    Rather than applying a Gaussian blur by default, might consider only 
    blurring if the labelim is very sparse (i.e. if there's only 1 segmentation
    as indicated by a length of 1 for f2sInds). Although blurring is needed to
    avoid aliasing effects common for resampling of single (or possibly more)
    framed segments, it's probably unneccessary for multi-framed ones.
    """
    
    if not interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f"The chosen interpolation, {interp}, is not one of the "\
              + "accepted arguments: 'NearestNeighbor', 'LabelGaussian', or "\
              + "'BlurThenLinear'."
        raise Exception(msg)
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of resample_labim():')
        print(f'Resampler (transform) Name = {sitkTx.GetName()}')
        print(f'                Parameters = {sitkTx.GetParameters()}')
        print(f'          Fixed Parameters = {sitkTx.GetFixedParameters()}')
        print(f'The chosen interpolation is {interp}.\n')
    
    """ 
    17/09/21: If interp = 'BlurThenLinear' but there are multiple indices in 
    f2sInds overwrite interp to 'NearestNeighbor' (or 'LabelGaussian'?) and
    overwrite applyPreResBlur and applyPostResBlur to False.
    """
    F = len(f2sInds)
    
    if 0: #F > 1 and interp == 'BlurThenLinear':
        #interp = 'NearestNeighbor'
        interp = 'LabelGaussian'
        print(f'*** Since there are {F} frames in this segment the',
              f'interpolation has been overwritten to {interp}.\n')
        
        if applyPreResBlur:
            applyPreResBlur = False
            print('*** The parameter applyPreResBlur has been overwritten',
                  f'to {applyPreResBlur}.\n')
        
        if applyPostResBlur:
            applyPostResBlur = False
            print('*** The parameter applyPostResBlur has been overwritten',
                  f'to {applyPostResBlur}.\n')
    
    # Store the interpolation set as metadata:
    labim.SetMetaData("resInterpSet", interp)
    
    if interp in ['NearestNeighbor', 'LabelGaussian']:
        if p2c:
            print(f'Attempting to resample labim using {interp} interpolator\n')
        
        if applyPreResBlur:
            # Gaussian blur labim:
            blurLabIm = gaussian_blur_im(im=labim, var=postResVar)
            
            # Resample blurLabIm using the chosen interpolator:
            resLabim = resample_im(
                im=blurLabIm, refIm=refIm, sitkTx=sitkTx, interp=interp
                )
            
            msg = 'Image info for resampled blurred image:'
        else:
            # Resample labim using the chosen interpolator:
            resLabim = resample_im(
                im=labim, refIm=refIm, sitkTx=sitkTx, interp=interp
                )
            
            msg = 'Image info for resampled image:'
            
        if p2c:
            print(msg)
        
        pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
            resLabim, p2c
            )
        if p2c:
            print('')
        
        if applyPostResBlur:
            # Gaussian blur resLabim:
            resLabim = gaussian_blur_im(im=resLabim, var=postResVar)
            
            if p2c:
                print('Image info for blurred resampled image:')
            pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
                resLabim, p2c
                )
            if p2c:
                print('')
            
            # Binarise resLabim if required:
            if len(uniqueVals) != 2 or sum(uniqueVals) != 1:
                """
                resLabim is not binary. Find suitable threshold value that 
                approximately preserves the number of pre-blurred truth values
                scaled by volumeRatio: 
                """
                thresh = find_thresh(
                    binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c
                    )
                
                # Binary threshold resLabim:
                resLabim = binarise_im(im=resLabim, thresh=thresh) 
                
                if p2c:
                    print(f'\nImage info after binary thresholding at {thresh}:')
                pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
                    resLabim, p2c
                    )
                if p2c:
                    print('')
        
        #print(f'\n   resF2Sinds = {resF2Sinds}')
        
        """ 
        Is resF2Sinds empty and not expected to be? If so try the 
        "BlurThenLinear" approach.
        
        If f2sInds isn't empty, resF2Sinds shouldn't be empty either.
        
        f2sInds will be empty if there were no segmentations/contours of 
        interest for the r^th ROI. In this case an empty resF2Sinds is
        acceptable. 
        """
        
        if resF2Sinds == []:
            print(f"There are {len(f2sInds)} non-empty masks in the input",
                  f"label image but {len(resF2Sinds)} non-empty frames in the",
                  f"resampled label image using {interp}. Will Gaussian blur,",
                  "linearly resample and binarise...\n")
            
            interp = 'BlurThenLinear'

    if interp == 'BlurThenLinear':
        # Gaussian blur labim:
        blurLabIm = gaussian_blur_im(im=labim, var=preResVar)
        
        if p2c:
            print('\nImage info for blurLabIm:')
            pixID, pixIDTypeAsStr, uniqueVals, f2sInds = get_im_info(
                blurLabIm, p2c
                )
            print('\n\nblurLabIm prior to resampling:')
            print(f'   blurLabIm.GetSize() = {blurLabIm.GetSize()}')
            print(f'   blurLabIm.GetSpacing() = {blurLabIm.GetSpacing()}')
            print(f'   refIm.GetSize() = {refIm.GetSize()}')
            print(f'   refIm.GetSpacing() = {refIm.GetSpacing()}')
        
        # Linearly resample blurLabIm:
        resLabim = resample_im(
            im=blurLabIm, refIm=refIm, sitkTx=sitkTx, interp='Linear'
            )
        
        """ 
        20/08/21: All zero value in resLabim, so instead try:
        """
        # TODO resolve this
        #print('\n\n*** Running sitk.Resample() rather than resample_im()..\n')
        #
        #resLabim = sitk.Resample(blurLabIm, refIm, sitkTx, 'Linear')
        """
        It was because sitkPixType was set to sitkUint32 instead of 
        sitkFloat32 in resample_im().
        """
        
        if p2c:
            print('\nImage info after resampling using linear interpolator:')
        pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
            resLabim, p2c
            )
        if p2c:
            print('')
        
        if applyPostResBlur:
            # Gaussian blur resLabim:
            resLabim = gaussian_blur_im(im=resLabim, var=postResVar)
            
        # Find suitable threshold value:
        thresh = find_thresh(binaryIm=labim, nonBinaryIm=resLabim, p2c=p2c)
        
        # Binary threshold resLabim:
        resLabim = binarise_im(im=resLabim, thresh=thresh) 
        
        if p2c:
            print(f'\nImage info after binary thresholding {thresh}:')
        pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
            resLabim, p2c
            )
        if p2c:
            print('')
    
    # Ensure that resLabim is a 32-bit unsigned integer (pixID = 5):
    if pixID != 5: 
        if p2c:
            print(f'\nresLabim has PixelID = {pixID} ({pixIDTypeAsStr})).')
        
        # Convert resLabim from float to 32-bit unsigned integer:
        resLabim = change_im_dtype(im=resLabim, newPixType='UInt32')
        
        if p2c:
            print('\nImage info after converting to 32-bit unsigned int:')
        #print(f'\nThe metadata keys are:', resLabim.GetMetaDataKeys())
        pixID, pixIDTypeAsStr, uniqueVals, resF2Sinds = get_im_info(
            resLabim, p2c
            )
        if p2c:
            print('')
    
    # Convert resLabim to a pixel array:
    resPixarr, resF2Sinds = im_to_pixarr(resLabim)
    
    # Store the interpolation used as metadata (which may be the same or 
    # different from the interpolation set): 
    resLabim.SetMetaData("resInterpUsed", interp)
    
    if interp == 'BlurThenLinear':
        # Store the threshold used as metadata:
        resLabim.SetMetaData("postResThreshUsed", f"{thresh}")
            
    if p2c:
        # The number of frames before and after:
        N_before = len(f2sInds)
        N_after = len(resF2Sinds)
        
        print(f'\nThere were {N_before} frames in the label image')
        print(f'There are {N_after} frames in the resampled label image')
        print('After converting resLabim to a pixel array:')
        print(f'resPixarr.shape = {resPixarr.shape}')
        print(f'resF2Sinds = {resF2Sinds}')
        plot_two_ims(
            im0=labim, ind0=f2sInds[0], plotLabel0='Original label image', 
            im1=resLabim, ind1=resF2Sinds[0], plotLabel1='Resampled label image')
        print('-'*120)
        
    return resLabim, resPixarr, resF2Sinds

def resample_labimBySeg(
        labimBySeg, f2sIndsBySeg, im, refIm, sitkTx=sitk.Transform(), 
        interp='NearestNeighbor', applyPreResBlur=False, preResVar=(1,1,1), 
        applyPostResBlur=True, postResVar=(1,1,1), p2c=False
        ):
    """
    Resample a list 3D SimpleITK images representing binary label images. 
    
    Parameters
    ----------  
    labimBySeg : list of SimpleITK Images
        A list (for each segment) of 3D label images to be resampled.
    f2sIndsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of slice numbers  
        that correspond to each frame in the label images.
    srcIm : SimpleITK Image
        The 3D image whose gridspace matches the gridspace of labimBySeg.
    refIm : SimpleITK Image
        The 3D image whose gridspace labimBySeg will be resampled to.
    sitkTx : SimpleITK Transform, optional (sitk.Transform() by default)
        The SimpleITK Transform to be used. The identity transform is default.
    interp : str, optional
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding). The default value
        is 'NearestNeighbor'.
    applyPreResBlur : bool, optional
        If True, the pre-resampled labim will be Gaussian blurred. This 
        argument is only applicable if interp is 'NearestNeighbor' or
        'LabelGaussian', since labim will be blurred prior to resampling if
        interp is 'BlurThenLinear' regardless of the argument's value. The
        default value is False.
    preResVar : tuple of floats, optional
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling. The default value is (1,1,1).
    applyPostResBlur : bool, optional
        If True, the post-resampled labelmap image will be Gaussian blurred.
        The default value is True.
    postResVar : tuple of floats, optional
        A tuple (for each dimension) of the variance to be applied if the 
        resampled labelmap image(s) is/are to be Gaussian blurred after  
        resampling. The default value is (1,1,1).
    p2c : bool, optional
        Denotes whether some results will be logged to the console. The
        default value is False.
    
    Returns
    -------
    resLabimBySeg : list of SimpleITK Images
        The list (for each segment) of resampled 3D label image.
    resPixarrBySeg : list of Numpy Arrays
        The list (for each segment) of the pixel array representations of
        resLabimBySeg.
    resF2SindsBySeg : list of a list of ints
        The list (for each segment) of a list (for each frame) of the
        frame-to-slice indices in resPixarrBySeg.
    
    Note
    ----
    See Notes in resample_labim. 
    """
        
    if p2c:
        print('\n\n', '-'*120)
        print('Running of resample_labimBySeg():')
        print('\n\n', '-'*120)
        
    resLabimBySeg = []
    resPixarrBySeg = []
    resF2SindsBySeg = []
    
    for r in range(len(labimBySeg)):
        if p2c:
            print(f'   Resampling of labimBySeg[{r}]...')
        
        resLabim, resPixarr, resF2Sinds\
            = resample_labim(
                labim=labimBySeg[r], f2sInds=f2sIndsBySeg[r], im=im, 
                refIm=refIm, sitkTx=sitkTx, interp=interp, 
                applyPreResBlur=applyPreResBlur, preResVar=preResVar, 
                applyPostResBlur=applyPostResBlur, postResVar=postResVar, 
                p2c=p2c
                )
        
        resLabimBySeg.append(resLabim)
        resPixarrBySeg.append(resPixarr)
        resF2SindsBySeg.append(resF2Sinds)
        
    if p2c:
        print('-'*120)
        
    return resLabimBySeg, resPixarrBySeg, resF2SindsBySeg

def align_ims_using_landmarks(
        transform, fixFidsFpath, movFidsFpath, fixIm, movIm, 
        flipK=False, buffer=3, p2c=False
        ):
    """
    Return the landmark-based aligned image and transform based on fiducials.
    
    Parameters
    ----------
    fixIm : SimpleITK Image 
        The 3D image that movIm will be aligned to.
    movIm : SimpleITK Image
        The 3D image that will be aligned to fixIm.
    fixFidsFpath : str, optional ('fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for fixIm. 
        The string need not contain the .txt extension.
    movFidsFpath : str, optional ('moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for movIm. 
        The string need not contain the .txt extension.
    flipK : bool, optional (False by default)
        Set to True if the fiducials specified in fixFidsFpath and movFidsFpath
        are indices obtained from ImageJ.
    buffer : int, optional (3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    transform : str, optional ('affine' by default)
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    alignedIm : SimpleITK Image
        The 3D aligned image.
    landmarkTx : SimpleITK Transform
        The transform used to perform landmark-based alignment.
        
    Note
    ----
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(fixIm.GetDimension())
        sitk.BSplineTransform(fixIm.GetDimension())
        
    numControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    landmarkTx, fixPts, movPts\
        = get_landmark_tx(
            transform=transform, fixFidsFpath=fixFidsFpath, 
            movFidsFpath=movFidsFpath, fixIm=fixIm, movIm=movIm, 
            flipK=flipK, buffer=buffer, p2c=p2c)
    
    #print(f'LandmarkTx.GetName() = {LandmarkTx.GetName()}\n')
    
    alignedIm = sitk.Resample(
        movIm, fixIm, landmarkTx, sitk.sitkLinear, 0.0, movIm.GetPixelID()
        )
    
    return alignedIm, landmarkTx