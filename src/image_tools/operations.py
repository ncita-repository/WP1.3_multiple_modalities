# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 14:20:01 2021

@author: ctorti
"""

""" Simple operations that apply to SimpleITK images. """

import time
import numpy as np
import SimpleITK as sitk
from conversion_tools.pixarrs_ims import im_to_pixarr
#from image_tools.attrs_info import get_im_info

def initialise_im(refIm, dtype='sameAsRefIm'):
    """
    Initialise an empty SimpleITK image with the same size, spacing, origin and
    direction as another image.  
    
    Parameters
    ----------
    refIm : SimpleITK Image
        The reference image whose attributes will be used for the initialised
        image.
    dtype : str, optional ('sameAsRefIm' by default)
        The pixel ID type as defined at the link in the Notes section. 
        Acceptable values are:
            - 'sameAsRefIm' (default)
            - 'uint8'
            - 'int8'
            - 'uint16'
            - 'int16'
            - 'uint32'
            - 'int32'
            - 'uint64'
            - 'int64'
            - 'float32'
            - 'float64'
        
    Returns
    -------
    newIm : SimpleITK Image
        An empty image with the same PixelID, Size, Spacing, Origin and 
        Direction as refIm.
    
    Notes
    -----
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    if dtype == 'sameAsRefIm':
        pixIDtype = refIm.GetPixelID()
    elif dtype == 'uint8':
        pixIDtype = sitk.sitkUInt8
    elif dtype == 'int8':
        pixIDtype = sitk.sitkInt8
    elif dtype == 'uint16':
        pixIDtype = sitk.sitkUInt16
    elif dtype == 'int16':
        pixIDtype = sitk.sitkInt16
    elif dtype == 'uint32':
        pixIDtype = sitk.sitkUInt32
    elif dtype == 'int32':
        pixIDtype = sitk.sitkInt32
    elif dtype == 'uint64':
        pixIDtype = sitk.sitkUInt64
    elif dtype == 'int64':
        pixIDtype = sitk.sitkInt64
    elif dtype == 'float32':
        pixIDtype = sitk.sitkFloat32
    elif dtype == 'float64':
        pixIDtype = sitk.sitkFloat64
    else:
        msg = f"dtype '{dtype}' is not a valid input."
        raise Exception(msg)
    
    #newIm = sitk.Image(refIm.GetSize(), refIm.GetPixelID())
    newIm = sitk.Image(refIm.GetSize(), pixIDtype)
    
    newIm.SetOrigin(refIm.GetOrigin())
    
    newIm.SetDirection(refIm.GetDirection())
    
    newIm.SetSpacing(refIm.GetSpacing())
    
    return newIm

def copy_im(im):
    """
    Copy a SimpleITK Image.  
    
    Parameters
    ----------
    im : SimpleITK Image
        The image to be copied.
        
    Returns
    -------
    imCopy : SimpleITK Image
        The copied image.
          
    Notes
    -----
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    imCopy = sitk.Image(im.GetSize(), im.GetPixelID())
    
    #imCopy.SetOrigin(im.GetOrigin())
    #imCopy.SetDirection(im.GetDirection())
    #imCopy.SetSpacing(im.GetSpacing())
    
    imCopy = sitk.Paste(im, imCopy)
    
    return imCopy

def im_max(im):
    """
    Find maximum of an SimpleITK Image.  
    
    Parameters
    ----------  
    im : SimpleITK Image
        
    Returns
    -------
    maximum : float or int
        Maximum value of im.
    """
    
    imFilt = sitk.MinimumMaximumImageFilter()
    imFilt.Execute(im)
    
    return imFilt.GetMaximum()

def im_min(im):
    """
    Find minimum of an SimpleITK Image.  
    
    Parameters
    ----------  
    im : SimpleITK Image
    
    Returns
    -------
    minimum : float or int
        Maximum value of im.
    """
    
    imFilt = sitk.MinimumMaximumImageFilter()
    imFilt.Execute(im)
    
    return imFilt.GetMinimum()

def get_im_stats(im, p2c=False):
    """
    Get various statistical values on a SimpleITK Image.

    Parameters
    ----------
    im : SimpleITK Image
    p2c : bool, optional
        If True stats will be printed to the console.

    Returns
    -------
    imMin : float
        Minimum value.
    imMax : float
        Maximum value.
    imMean : float
        Mean value.
    imSigma : float
        The standard deviation.
    imVariance : float
    imSum : float
        The sum of all voxels.
    """
    imFilt = sitk.StatisticsImageFilter()
    imFilt.Execute(im)
    
    imMin = imFilt.GetMinimum()
    imMax = imFilt.GetMaximum()
    imMean = imFilt.GetMean()
    imSigma = imFilt.GetSigma()
    imVariance = imFilt.GetVariance()
    imSum = imFilt.GetSum()
    
    if p2c:
        print(f'imMin = {imMin}')
        print(f'imMax = {imMax}')
        print(f'imMean = {imMean}')
        print(f'imSigma = {imSigma}')
        print(f'imVariance = {imVariance}')
        print(f'imSum = {imSum}')
    
    return imMin, imMax, imMean, imSigma, imVariance, imSum

def im_max_by_frame(im):
    """
    Return a list of the maximum value of each frame in a SimpleITK Image .  
    
    Parameters
    ----------
    im : SimpleITK Image
        
    Returns
    -------
    maxima : list of float or int
        A list of the maximum value of each frame in im.
    """
    
    F = im.GetSize()[2]
    
    maxima = []
    
    for f in range(F):
        maxima.append(im_max(im[:, :, f]))
    
    return maxima

def im_or(im0, im1):
    """
    Perform pixel-wise OR operation on two SimpleITK images.  
    
    Parameters
    ----------
    im0 : SimpleITK Image
    im1 : SimpleITK Image
        
    Returns
    -------
    imOr : SimpleITK Image
        The pixel-wise OR of im0 and im1.
    """
            
    imFilt = sitk.OrImageFilter()
    imOr = imFilt.Execute(im0, im1)
    
    return imOr

def im_add(im0, im1):
    """
    Perform pixel-wise addition of two SimpleITK images.  
    
    Parameters
    ----------
    im0 : SimpleITK Image
    im1 : SimpleITK Image
    
    Returns
    -------
    imSum : SimpleITK Image
        The pixel-wise sum of image0 and image1.
    """
    
    imFilt = sitk.AddImageFilter()
    imSum = imFilt.Execute(im0, im1)
    
    return imSum

def normalise_im(im, normBy='sum'):
    """
    Normalise a SimpleITK image to 1.0.  
    
    Parameters
    ----------
    im : SimpleITK Image
    normBy : str, optional
        Denotes how the image will be normalised. Acceptable inputs include:
            - 'sum' --> normalise by the sum of all voxels
            - 'max' --> normalise by the maximum value
        The default value is 'sum'.
    
    Returns
    -------
    normIm : SimpleITK Image
        im normalised to 1.0
    """
    
    if not isinstance(normBy, str):
        msg = "The input 'normBy' must be a string."
        raise Exception(msg)
    
    imMin, imMax, imMean, imSigma, imVariance, imSum = get_im_stats(im)
    
    if normBy == 'sum':
        normIm = im/imSum
    elif normBy == 'max':
        normIm = im/imMax
    else:
        msg = "The input 'normBy' must be 'sum' or 'max'. '{normBy}' is " +\
            "not an acceptable value."
        raise Exception(msg)
        
    return normIm
    
def pixarr_sum_of_ims(ims):
    """
    Return Numpy pixel array representation of the sum of a list of images.
    
    Useful for binary label images.

    Parameters
    ----------
    images : list of SimpleITK Images

    Returns
    -------
    sumPixarr : Numpy data array
        Pixel-wise sum of all images.
    """
    
    # Initialise sumIm:
    sumIm = initialise_im(ims[0])
    
    for i in range(1, len(ims)):
        sumIm = im_add(sumIm, ims[i])

    # Convert to numpy arrays:
    sumPixarr = sitk.GetArrayFromImage(sumIm)
    
    return sumPixarr

def get_pixel_area_ratio(im0, im1):
    """
    Get the pixel area ratio of SimpleITK images.
    
    Parameters
    ----------
    im0 : SimpleITK Image
    im1 : SimpleITK Image
    
    Returns
    -------
    ratio : float
        The ratio of the pixel area of im1 to im0.
    """
    
    area0 = np.prod(np.array(im0.GetSpacing())[0:2])
    
    area1 = np.prod(np.array(im1.GetSpacing())[0:2])
    
    return area1/area0

def get_voxel_volume_ratio(im0, im1):
    """
    Get the voxel volume ratio of SimpleITK images.
    
    Parameters
    ----------
    im0 : SimpleITK Image
    im1 : SimpleITK Image
    
    Returns
    -------
    ratio : float
        The ratio of the voxel volume of im1 to im0.
    """
    
    vol0 = np.prod(np.array(im0.GetSpacing()))
    
    vol1 = np.prod(np.array(im1.GetSpacing()))
    
    return vol1/vol0

def change_im_dtype(im, newPixType):
    """
    Convert a 3D SimpleITK Image to a new pixel type.
    
    Parameters
    ----------   
    im : SimpleITK Image
        The 3D image whose pixel type is to be converted.
    newPixType : str
        The pixel type to convert im to.  Acceptable inputs are a subset of
        the SimpleITK pixel types 
        (https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html):
         
        - 'UInt8' 	Unsigned 8-bit integer
        - 'UInt16' 	Unsigned 16-bit integer
        - 'UInt32' 	Unsigned 32-bit integer
        - 'UInt64' 	Unsigned 64-bit integer
        - 'Int8' 	    Signed 8-bit integer
        - 'Int16' 	Signed 16-bit integer
        - 'Int32' 	Signed 32-bit integer
        - 'Int64' 	Signed 64-bit integer
        - 'Float32' 	32-bit float
        - 'Float64' 	64-bit float
    
    Returns
    -------
    SimpleITK Image
        The 3D image with modified pixel type.
        
    Note
    ----
    While a linear (or BSpline) interpolator is appropriate for intensity 
    images, only a NearestNeighbor interpolator is appropriate for binary 
    images (e.g. segmentations) so that no new labels are introduced.
    """
    
    if newPixType == 'UInt8':
        pixTypeOut = sitk.sitkUInt8
    elif newPixType == 'UInt16':
        pixTypeOut = sitk.sitkUInt16
    elif newPixType == 'UInt32':
        pixTypeOut = sitk.sitkUInt32
    elif newPixType == 'UInt64':
        pixTypeOut = sitk.sitkUInt64
    elif newPixType == 'Int8':
        pixTypeOut = sitk.sitkInt8
    elif newPixType == 'Int16':
        pixTypeOut = sitk.sitkInt16
    elif newPixType == 'Int32':
        pixTypeOut = sitk.sitkInt32
    elif newPixType == 'Int64':
        pixTypeOut = sitk.sitkInt64
    elif newPixType == 'Float32':
        pixTypeOut = sitk.sitkFloat32
    elif newPixType == 'Float64':
        pixTypeOut = sitk.sitkFloat64
    else:
        msg = "'newPixType' must be either:\n 'UInt8'\n 'UInt16'\n 'UInt32'"\
              + "\n 'UInt64'\n 'Int8'\n 'Int16'\n 'Int32'\n 'Int64'"\
              + "\n 'Float32'\n 'Float64'."
        raise Exception(msg)
    
    Caster = sitk.CastImageFilter()
    
    Caster.SetOutputPixelType(pixTypeOut)
    
    return Caster.Execute(im)

def find_thresh(binaryIm, nonBinaryIm, p2c=False):
    """
    Find a suitable threshold level that if used to binary threshold a non-
    binary 3D SimpleITK image, nonBinaryIm, would result in a similar 
    distribution of zeros and ones in the binarised image as that of the
    original binary image, binaryIm.
    
    Parameters
    ---------- 
    binaryIm : SimpleITK Image 
        The binary image.
    nonBinaryIm : SimpleITK Image
        The non-binary image.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    thresh : float
        A threshold level that when applied to nonBinaryIm should result in 
        a similar distribution of zeros and ones as in binaryIm.
    """
    
    #from conversion_tools.pixarrs_ims import im_to_pixarr
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of find_thresh():')
        #print('\n\n', '-'*120)
    
    ratio = get_voxel_volume_ratio(binaryIm, nonBinaryIm)
    
    binPixarr, binf2sInds = im_to_pixarr(binaryIm)
    
    nonBinPixarr, nonBinf2sInds = im_to_pixarr(nonBinaryIm)
    
    numOfOnes_B = np.count_nonzero(binPixarr)
    
    # The expected number of 1s to retain in nonBinPixarr (i.e. expected
    # number of 1st in the binarised non-binary pixel array):
    expectedNumOfOnes_BNB = int(numOfOnes_B/ratio)
            
    vals_NB = nonBinPixarr.flatten()
    minVal = min(vals_NB)
    maxVal = max(vals_NB)
    
    freq, bins = np.histogram(vals_NB, bins=1000, range=[minVal, maxVal])
    
    """
    if p2c:
        print('\n      Distribution of values:')
        #for b, f in zip(bins[1:], freq):
        #    print(f'      {round(b, 8)} - {f}')
        for i in reversed(range(len(freq))):
                valRange = f'[{round(bins[i], 3)} - {round(bins[i+1], 3)}]'
                print(f'      {valRange} : {freq[i]}')
    """
    
    # The cumulative sum of the frequencies in reverse order:
    cumSumfreq = np.cumsum(np.flip(freq))
    
    # The absolute differences between the cumulative sum of the frequencies
    # and the number of ones in the binary image scaled by the ratio of the
    # voxel volumes:
    #absDiffs = np.abs(cumSumfreq - numOfOnes_B)
    #absDiffs = np.abs(cumSumfreq - numOfOnes_B/ratio) # 20/08/21
    absDiffs = np.abs(cumSumfreq - expectedNumOfOnes_BNB) # 20/08/21
    
    revInd = np.where(absDiffs == np.min(absDiffs))[0][0]
    
    ind = len(cumSumfreq) - revInd - 1
    
    thresh = bins[ind]
    
    if p2c:
        print(f'\nbinaryIm.GetSpacing() = {binaryIm.GetSpacing()}')
        print(f'\nnonBinaryIm.GetSpacing() = {nonBinaryIm.GetSpacing()}')
        print(f'\nVoxel volume ratio = {ratio}')
        print(f'\nNumber of 1s in binaryIm = {numOfOnes_B}')
        print(f'\nTarget number of ones = {expectedNumOfOnes_BNB}')
        #print(f'\ncumSumfreq[0:revInd+5] = {cumSumfreq[:revInd+5]}')
        #print(f'\nabsDiffs[0:revInd+5] = {absDiffs[:revInd+5]}')
        print(f'\nind = {ind}')
        print(f'bins[{ind}] = {bins[ind]}')
    
    binPixarr = (nonBinPixarr > thresh).astype(np.int_)
    
    numOfOnes_BNB = np.count_nonzero(binPixarr)
    
    if p2c:
        print(f'There are {numOfOnes_B} ones in the binary pixel array. ',
              f'\nThresholding the non-binary pixel array at {round(thresh,3)}',
              f'would result in {numOfOnes_BNB} ones (difference of',
              #f'{abs(numOfOnes_BNB - numOfOnes_B)}).')
              f'{abs(numOfOnes_BNB - expectedNumOfOnes_BNB)}).')
        print('-'*120)
    
    return thresh

def binarise_im(im, thresh=0.5):
    """
    Binary threshold a 3D SimpleITK image.
    
    Parameters
    ---------- 
    im : SimpleITK Image 
        The 3D image to be transformed.
    thresh : float, optional (0.5 by default)
        Lower limit threshold.
        
    Returns
    -------
    SimpleITK Image
        The 3D binary threshold image.
    
    Note
    ----
    If an empty image (all zeros) is binary thresholded, the resuling image
    has ones everywhere for some unknown reason.
    """
    
    imFilt = sitk.BinaryThresholdImageFilter()
    
    #imFilt.SetInput(im)
    
    """
    Setting UpperThreshold (e.g. to 1) removes all pixels that exceed 1, and
    for some reason, transformed binary labelmap images result in non-binary 
    values (including negative values).
    """
    
    imFilt.SetLowerThreshold(thresh)
    #imFilt.SetUpperThreshold(1)
    imFilt.SetOutsideValue(0)
    imFilt.SetInsideValue(1)
    #imFilt.SetOutputPixelType(sitk.sitkUInt32) # --> error:
    # 'BinaryThresholdImageFilter' object has no attribute 'SetOutputPixelType'
    #imFilt.DebugOn()
    
    #imFilt.Execute()
    #binIm = imFilt.GetOutput()
    
    return imFilt.Execute(im)

def gaussian_blur_im(im, var=(1,1,1), p2c=False):
    """
    Gaussian blur a 3D SimpleITK image.
    
    Parameters
    ---------- 
    im : SimpleITK Image 
        The 3D image to be blurred.
    var : tuple of floats, optional ((1,1,1) by default)
        The variance along all dimensions.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    blurredIm : SimpleITK Image
        The 3D blurred image.
    
    Note
    ----
    DiscreteGaussianImageFilter:
        Blurs an image by separable convolution with discrete gaussian kernels. 
        This filter performs Gaussian blurring by separable convolution of an 
        image and a discrete Gaussian operator (kernel).
        
    SmoothingRecursiveGaussianImageFilter:
        Computes the smoothing of an image by convolution with the Gaussian 
        kernels implemented as IIR filters
        
    The input image type must be float - otherwise odd results will arise. So
    check will be made prior to Gaussian blurring and the image type will be
    converted to float if applicable.
    """
    
    from image_tools.attrs_info import get_im_info
    
    if p2c:
        print('\n   Image info for im prior to Gaussian blurring:')
            
    pixID, pixIdtypeAsStr, uniqueValsPostTx, f2sInds = get_im_info(im, p2c)
    
    # Ensure im is a 32-bit float (pixID = 8).
    if pixID != 8: 
        if p2c:
            print(f'\n   im has PixelID = {pixID} ({pixIdtypeAsStr})). '\
                  'It will be converted to Float32.')
        
        # Convert image to Float32:
        im = change_im_dtype(im=im, newPixType='Float32')
        
        if p2c:
            print('\n   Image info for im after converting to Float32:')
            
        pixID, pixIdtypeAsStr, uniqueVals, f2sInds = get_im_info(im, p2c)
    
    imFilt = sitk.DiscreteGaussianImageFilter()
    #imFilt.SetMaximumKernelWidth(sigma)
    #print(f'   imFilt.GetMaximumKernelWidth() = {imFilt.GetMaximumKernelWidth()}')
    imFilt.SetVariance(var)
        
    blurredIm = imFilt.Execute(im)
    
    if p2c:
        print('\nApplying Gaussian blur...')
        print(f'   imFilt.GetMaximumError() = {imFilt.GetMaximumError()}')
        print(f'   imFilt.GetMaximumKernelWidth() = {imFilt.GetMaximumKernelWidth()}')
        print(f'   imFilt.GetUseImageSpacing() = {imFilt.GetUseImageSpacing()}')
        print(f'   imFilt.GetVariance() = {imFilt.GetVariance()}')
        print(f'\n   im.GetPixelID() = {im.GetPixelID()}')
        print(f'   im.GetPixelIDValue() = {im.GetPixelIDValue()}')
        print(f'   im.GetPixelIDTypeAsString() = {im.GetPixelIDTypeAsString()}')
        print(f'   blurredIm.GetPixelID() = {blurredIm.GetPixelID()}')
        print(f'   blurredIm.GetPixelIDValue() = {blurredIm.GetPixelIDValue()}')
        print(f'   blurredIm.GetPixelIDTypeAsString() = {blurredIm.GetPixelIDTypeAsString()}')
    
    return blurredIm

def recursive_gaussian_blur_im(im, sigma, direction, p2c=False):
    """
    Apply a recursive Gaussian blur to a 3D SimpleITK Image.
    
    Parameters
    ---------- 
    im : SimpleITK Image 
        The 3D image to be blurred.
    sigma : int
        The kernel width.
    direction : int
        The direction to apply blurring to. direction can be 0, 1, or 2 for the
        x, y or z direction.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    blurredIm : SimpleITK Image
        The 3D blurred image.
    
    Note
    ----
    http://itk-users.7.n7.nabble.com/ITK-users-sitk-DiscreteGaussianImageFilter-with-different-variances-per-direction-td38185.html
    
    The input image type must be float - otherwise odd results will arise. So
    check will be made prior to Gaussian blurring and the image type will be
    converted to float if applicable.
    """
    
    from image_tools.attrs_info import get_im_info
    
    if p2c:
        print('\n   Image info for im prior to Gaussian blurring:')
    
    pixID, pixIdtypeAsStr, uniqueValsPostTx, f2sInds = get_im_info(im, p2c)
    
    # Ensure im is a 32-bit float (pixID = 8).
    if pixID != 8: 
        if p2c:
            print(f'\n   im has PixelID = {pixID} ({pixIdtypeAsStr})). '\
                  'It will be converted to Float32.')
        
        # Convert im to Float32:
        im = change_im_dtype(im=im, newPixType='Float32')
        
        if p2c:
            print('\n   Image info for image after converting to Float32:')
            
        pixID, pixIdtypeAsStr, uniqueVals, f2sInds = get_im_info(im, p2c)
        
    imFilt = sitk.RecursiveGaussianImageFilter()
    imFilt.SetSigma(sigma)
    imFilt.SetDirection(direction)
    
    blurredIm = imFilt.Execute(im)
    
    if p2c:
        print('\nApplying Gaussian blur...')
        print(f'   imFilt.GetDirection() = {imFilt.GetDirection()}')
        print(f'   imFilt.GetNormalizeAcrossScale() = {imFilt.GetNormalizeAcrossScale()}')
        print(f'   imFilt.GetOrder() = {imFilt.GetOrder()}')
        print(f'   imFilt.GetSigma() = {imFilt.GetSigma()}')
        print(f'\n   im.GetPixelID() = {im.GetPixelID()}')
        print(f'   im.GetPixelIDValue() = {im.GetPixelIDValue()}')
        print(f'   im.GetPixelIDTypeAsString() = {im.GetPixelIDTypeAsString()}')
        print(f'   blurredIm.GetPixelID() = {blurredIm.GetPixelID()}')
        print(f'   blurredIm.GetPixelIDValue() = {blurredIm.GetPixelIDValue()}')
        print(f'   blurredIm.GetPixelIDTypeAsString() = {blurredIm.GetPixelIDTypeAsString()}')
    
    return blurredIm

def segment_im(
        im, closeHoles=True, threshFactor=0.25, kernelSize=3, 
        numClosingOps=1, p2c=False
        ):
    """
    Segment a 3D SimpleITK Image using binary thresholding and (optional) a 
    morphological closing.
    
    Parameters
    ---------- 
    im : SimpleITK Image 
        The 3D image to be segmented.
    closeHoles : bool, optional
        If True the binary thresholded image will be morphologically closed.
        The default value is True.
    threshFactor : float, optional
        The factor to apply to the mean value of im for binary thresholding,
        e.g. 0.25*meanVal where meanVal is the mean value of im. The default
        value is 0.25.
    kernelSize : int, optional
        The kernel size to use for morphological closing. The default value is
        3.
    numClosingOps : int, optional
        The number of times to perform the hole closing operation. The default
        value is 1.
    p2c : bool, optional
        Denotes whether some results will be logged to the console. The default
        value is False.
    
    Returns
    -------
    seg : SimpleITK Image
        The binary segmented image.
    """
    
    # Start timing:
    times = [time.time()]    
    
    imMin, imMax, imMean, imSigma, imVariance, imSum = get_im_stats(
        im, p2c=False
        )
    
    seg = im > threshFactor*imMean
    
    if closeHoles:
        closing = sitk.BinaryMorphologicalClosingImageFilter()
        for i in range(numClosingOps):
            closing.SetKernelRadius([kernelSize]*im.GetDimension())
            seg = closing.Execute(seg)
    
    times.append(time.time())
    dTime = times[-1] - times[-2]
    if p2c:
        print(f'*Took {dTime:.4f} s to segment the image.\n')
    
    return seg