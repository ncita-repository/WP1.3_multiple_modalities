# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:52:44 2021

@author: ctorti
"""

""" High-level conversions between pixel arrays and images and visa versa. """

import numpy as np
import SimpleITK as sitk
from conversion_tools.inds_pts_pixarrs import pixarr_to_labarr
from image_tools.attrs_info import get_im_info


def pixarr_to_im(pixarr, f2sInds, refIm):
    """ 
    Convert a 3D pixel array to a 3D SimpleITK image.  
    
    Parameters
    ----------
    pixarr : Numpy array
        PixelData as a Numpy array.
    f2sInds : list of ints
        The DICOM slice numbers that correspond to each frame in pixarr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number
    refIm : SimpleITK Image
        The 3D image that pixarr relates to.
    
    Returns
    -------
    labim : SimpleITK Image
        A SimpleITK image containing zero-padded indices from pixarr.
    """
    
    imSize = refIm.GetSize()
    
    #print(f'\nimSize = {imSize}')
    #print(f'pixarr.shape = {pixarr.shape}')
    #print(f'imSize[2] = {imSize[2]}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    
    labarr = pixarr_to_labarr(
        pixarr=pixarr, numOfSlices=imSize[2], f2sInds=f2sInds
        )
    
    labim = sitk.GetImageFromArray(labarr)
    
    # Set the Origin, Spacing and Direction of labim to that of refIm:
    labim.SetOrigin(refIm.GetOrigin())
    labim.SetSpacing(refIm.GetSpacing())
    labim.SetDirection(refIm.GetDirection())
    
    #print(f'\n\nThe maximum value in Labmap is {Labmap.max()}')
    #print(f'\n\nThe maximum value in labim is {ImageMax(labim)}')
        
    return labim

def pixarr_to_imsBySeg(pixarr, f2sIndsBySeg, frameNumsBySeg, refIm):
    """
    Convert a 3D pixel array to a list of 3D SimpleITK images - one per
    segment.  
    
    Parameters
    ----------
    pixarr : Numpy array
        PixelData as a Numpy array.
    f2sIndsBySeg : list of list of ints
        List of the DICOM slice numbers that correspond to each frame in pixarr 
        (e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number) - a list for each segment.
    frameNumsBySeg : list of list of ints
        List of the frame numbers that belong to each segment - one list per
        segment.
    refIm : SimpleITK Image
        The 3D image that pixarr relates to.
    
    Returns
    -------
    ims : List of SimpleITK Images
        A list of 3D zero-padded labelaps - one per segment.
    """
    
    ims = []
    
    #print(f'\nf2sIndsBySeg = {f2sIndsBySeg}\n')
    #print(f'frameNumsBySeg = {frameNumsBySeg}')
    
    for i in range(len(frameNumsBySeg)):
        # The starting (= a) and ending (= b) frame numbers:
        a = frameNumsBySeg[i][0]
        b = frameNumsBySeg[i][-1]
        
        #print(f'a = {a}, b = {b}')
        #print(f'pixarr.shape = {pixarr.shape}')
        #print(f'pixarr[{a}:{b+1}].shape = {pixarr[a:b+1].shape}')
        
        im = pixarr_to_im(pixarr=pixarr[a:b+1],
                          refIm=refIm,
                          f2sInds=f2sIndsBySeg[i])
        ims.append(im)
        
    return ims

def pixarrByRoi_to_imByRoi(pixarrByRoi, f2sIndsByRoi, refIm, p2c=False):
    """
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates for each
    pixel array in a list of pixel-arrays-by-ROI.
 
    Parameters
    ----------
    pixarrByRoi : list of Numpy arrays
        A list (for each ROI) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    f2sIndsByRoi : list of list of ints
        A list (for each ROI) of a list (for each frame) of the slice numbers 
        that correspond to each frame in each pixel array in pixarrByRoi.  This
        is equivalent to a list of Per-FrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds). 
    refIm : SimpleITK image
        The 3D image that pixarr relates to (e.g. SrcIm).
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    labimByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images representing the pixel
        arrays in pixarrByRoi.
    """
    
    #from ImageTools import GetImageInfo
    #from GeneralTools import UniqueItems, AreListsEqualToWithinEpsilon
    
    if p2c:
        print('\n\nStart', '-'*110)
        print('Running of pixarrByRoi_to_imByRoi():')
        print('-'*120)
        #print(f'   len(pixarrByRoi) = {len(pixarrByRoi)}')
        
    labimByRoi = []
    f2sIndsByRoi_new = []

    for r in range(len(pixarrByRoi)):
        labim = pixarr_to_im(
            pixarr=pixarrByRoi[r], f2sInds=f2sIndsByRoi[r], refIm=refIm
            )
        
        labimByRoi.append(labim)
    
        if p2c:
            print(f'\n   Image info for labimByRoi[{r}]:')
            
        pixID, pixIDtypeAsStr, uniqueVals, f2sInds = get_im_info(labim, p2c)
        
        f2sIndsByRoi_new.append(f2sInds)
    
    #""" Check that the number of unique f2sInds in f2sIndsByRoi_new matches the
    #number of unique f2sInds in f2sIndsByRoi (i.e. that no frames were lost). """
    #NumIn = []
    #NumOut = []
    #
    #for r in range(len(f2sIndsByRoi)):
    #    NumIn.append(len(UniqueItems(f2sIndsByRoi[r])))
    #    
    #    NumOut.append(len(UniqueItems(f2sIndsByRoi_new[r])))
    #    
    #print('\n', NumIn, NumOut)
    #
    #if not AreListsEqualToWithinEpsilon(NumIn, NumOut, epsilon=0.9):
    #    print('\nWARNING:  The number of unique f2sIndsByRoi used to generate',
    #          f'the labelmaps, {NumIn}, doesn\'t match that in the labelmap',
    #          f'images, {NumOut}.')
    
    if p2c:
        print('End', '-'*120)
        
    return labimByRoi

def im_to_pixarr(image, f2sInds=None):
    """
    Convert a 3D SimpleITK Image to a 3D pixel array.  
    
    Parameters
    ----------
    image : SimpleITK Image
        A image (e.g. labelmap).
    f2sInds : list of ints, optional
        Use PerFrameFunctionalGroupsSequence-to-slice indices 
        (PFFGStoSliceInds) if they are known. If this input is not specified a 
        list of the frame numbers that are non-zero will be determined. The
        default value is None.
        
    Returns
    -------
    pixarr : Numpy array
        The non-zero frames in Image as a Numpy array.
    f2sInds : list of ints
        List of the indices of non-zero frames in Image.  This should be
        equivalent to the PerFrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds).
    """
    
    # Convert image to a Numpy data array:
    labarr = sitk.GetArrayFromImage(image)
    
    if f2sInds == None:
        f2sInds = []
        
        # Loop through all frames in labarr and get the sum of all pixels:
        for i in range(labarr.shape[0]):
            labarrSum = np.amax(labarr[i])
            
            #print(f'i = {i}, labarrSum = {labarrSum}')
            
            if labarrSum:
                # This is a non-zero frame:
                f2sInds.append(i)
    
    # Initialise pixarr:
    pixarr = np.zeros((len(f2sInds), labarr.shape[1], labarr.shape[2]))
        
    # Use f2sInds to index the non-zero frames in labarr:
    for i in range(len(f2sInds)):
        pixarr[i] = labarr[f2sInds[i]]
        
    return pixarr, f2sInds

def im_to_binary_ones_im(image):
    """
    Create a 3D binary ones SimpleITK Image from a non-binary 3D SimpleITK 
    Image.  
    
    Parameters
    ----------
    image : SimpleITK Image
    
    Returns
    -------
    maskIm : SimpleITK Image
        A 3D ones SimpleITK image with the same image attributes as Image.
    
    Note
    ----
    This function is similar to pixarr_to_im so consider combining.
    """
    
    C, R, F = image.GetSize()
    
    pixarr = np.ones((F, R, C), dtype=np.uint8)
    
    maskIm = sitk.GetImageFromArray(pixarr)
    
    # Set the Origin, Spacing and Direction of maskIm to that of image:
    maskIm.SetOrigin(image.GetOrigin())
    maskIm.SetSpacing(image.GetSpacing())
    maskIm.SetDirection(image.GetDirection())
    
    return maskIm

def imBySeg_to_pixarrBySeg(imBySeg, p2c=False):
    """
    Convert a list of 3D SimpleITK Images to a list of 3D SEG pixel arrays.  
    
    Parameters
    ----------
    imBySeg : list of SimpleITK Images
        A list (e.g. for each segment) of images (e.g. labelmap).
    p2c : bool, optional
        If True results will be printed to the console. The default value is
        False.
        
    Returns
    -------
    pixarrBySeg : list of Numpy arrays
        A list (for each segment) of the non-zero frames in each image in
        ImagebySeg.
    f2sInds : list of a list of ints
        List (for each segment) of a list (for each frame) of the indices of 
        non-zero frames in each image in imBySeg.
    """
    
    if p2c:
        print('\n\n' + '-'*120)
        print('---> imBySeg_to_pixarrBySeg():\n')
        #print('-'*120)
    
    pixarrBySeg = []
    f2sIndsBySeg = []
    
    for i in range(len(imBySeg)):
        pixarr, f2sInds = im_to_pixarr(imBySeg[i])
        
        if p2c:
            print(f'{i}: f2sInds = {f2sInds}\n')
            
        pixarrBySeg.append(pixarr)
        
        f2sIndsBySeg.append(f2sInds)
    
    if p2c:
        print('<--- imBySeg_to_pixarrBySeg()')
        print('-'*120 + '\n')
    
    return pixarrBySeg, f2sIndsBySeg