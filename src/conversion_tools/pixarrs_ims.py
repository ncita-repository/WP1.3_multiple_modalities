# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:52:44 2021

@author: ctorti
"""

""" High-level conversions between pixel arrays and images and visa versa. """


import numpy as np
import SimpleITK as sitk
from conversion_tools.inds_pts_pixarrs import pixarr_to_labarr
#from image_tools.attrs_info import get_im_info
from general_tools.console_printing import (
    print_indsByRoi, print_ptsByCntByRoi, print_pixarrBySeg, print_labimBySeg
    )


def pixarr_to_im(pixarr, f2sInds, refIm):
    """ 
    Convert a 3D pixel array to a 3D SimpleITK image.  
    
    Parameters
    ----------
    pixarr : Numpy array
        PixelData as a FxRxC Numpy array, where F is the number of
    segmentations/masks.
    f2sInds : list of ints
        The DICOM slice numbers that correspond to each frame in pixarr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number
    refIm : SimpleITK Image
        The 3D image that pixarr relates to.
    
    Returns
    -------
    labim : SimpleITK Image
        A CxRxS SimpleITK image, where S is the number of slices in the DICOM
        series.
    """
    
    imSize = refIm.GetSize()
    
    #print(f'\nimSize = {imSize}')
    #print(f'pixarr.shape = {pixarr.shape}')
    #print(f'imSize[2] = {imSize[2]}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    
    # Convert from sparse pixel array to zero-padded label array:
    """
    pixarr is a FxRxC array, where F is the number of segmentations/masks
    labarr is a SxRxC array, where S is the number of slices in the DICOM 
    series
    """
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

def pixarrBySeg_to_labimBySeg(pixarrBySeg, f2sIndsBySeg, refIm, p2c=False):
    """
    Convert a 3D pixel array to a list (for each segment) of SimpleITK Image
    representations of the pixel arrays.
 
    Parameters
    ----------
    pixarrBySeg : list of Numpy arrays
        A list (for each ROI) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    f2sIndsBySeg : list of list of ints
        A list (for each ROI) of a list (for each frame) of the slice numbers 
        that correspond to each frame in each pixel array in pixarrBySeg.  This
        is equivalent to a list of Per-FrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds). 
    refIm : SimpleITK image
        The 3D image that pixarr relates to (e.g. SrcIm).
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    labimBySeg : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images representing the pixel
        arrays in pixarrBySeg.
    newF2SIndsBySeg list of list of ints
        As f2sIndsBySeg but the new list may have fewer items, since any 
        repeated index in f2sIndsBySeg (i.e. more than one segmentated object
        on any given slice) will be combined into the same frame in the label
        image. 
    """
    
    #from ImageTools import GetImageInfo
    #from GeneralTools import UniqueItems, AreListsEqualToWithinEpsilon
    from image_tools.attrs_info import get_im_info
    
    if p2c:
        print('\n\nStart', '-'*110)
        print('Running of pixarrBySeg_to_imByRoi():')
        print('-'*120)
        #print(f'   len(pixarrBySeg) = {len(pixarrBySeg)}')
        
    labimBySeg = []
    newF2SIndsBySeg = []

    for r in range(len(pixarrBySeg)):
        labim = pixarr_to_im(
            pixarr=pixarrBySeg[r], f2sInds=f2sIndsBySeg[r], refIm=refIm
            )
        
        labimBySeg.append(labim)
    
        if p2c:
            print(f'\n   Image info for labimBySeg[{r}]:')
            
        pixID, pixIDtypeAsStr, uniqueVals, f2sInds = get_im_info(labim, p2c)
        
        newF2SIndsBySeg.append(f2sInds)
    
    if p2c:
        print('\nThe (original) f2sIndsBySeg corresponding to the input',
              'pixarrBySeg:')
        print_indsByRoi(f2sIndsBySeg)
        print('\nThe new f2sIndsBySeg corresponding to the label image:')
        print_indsByRoi(newF2SIndsBySeg)
        print('End', '-'*120)
        
    return labimBySeg, newF2SIndsBySeg

def im_to_pixarr(image, f2sInds=None):
    """
    Convert an CxRxS SimpleITK Image to a sparse FxRxC pixel array, where C/R
    the number of columns/rows, S is the number of slices in the DICOM stack, 
    and F is the number of non-empty frames in image.
    
    This function can be used for any image type, but is specifically used for
    converting binary label images to sparse 3D pixel arrays (i.e. storing 
    masks/segmentations).  To convert an CxRxS image to a SxRxC pixel array use
    sitk.GetArrayViewFromImage(image).
    
    Parameters
    ----------
    image : SimpleITK Image
        A image (e.g. label image).
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