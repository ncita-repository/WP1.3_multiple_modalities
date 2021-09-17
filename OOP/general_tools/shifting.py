# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:24:52 2021

@author: ctorti
"""

import numpy as np
from copy import deepcopy
from general_tools.general import get_unique_items
#from io_tools.imports import import_im
from conversion_tools.inds_pts_cntdata import pts_to_inds, inds_to_pts
from general_tools.console_printing import print_indsByRoi, print_ptsByCntByRoi

def replace_ind_in_C2SindsByRoi(c2sIndsByRoi, indToReplace, replacementInd):
    
    #newC2SindsByRoi = list(c2sIndsByRoi) # this modifies c2sIndsByRoi despite
    # use of list() to make a copy
    newC2SindsByRoi = deepcopy(c2sIndsByRoi)
    
    for r in range(len(newC2SindsByRoi)):
        for c in range(len(newC2SindsByRoi[r])):
            if newC2SindsByRoi[r][c] == indToReplace:
                newC2SindsByRoi[r][c] = replacementInd
    
    return newC2SindsByRoi

def get_phys_shift_bt_slices(image0, sliceNum0, image1, sliceNum1):
    """
    Get the physical shift between sliceNum0 in image0 and sliceNum1 in image1
    based on the origins of the two slices (i.e. image0[0,0,sliceNum0] and
    image1[0,0,sliceNum1]).
    
    Parameters
    ----------
    image0 : SimpleITK image
    sliceNum0 : int
        The slice number in image0.
    image1 : SimpleITK image
    sliceNum1 : int
        The slice number in image1.
    
    Returns
    -------
    mmShift : Numpy array
        Shift in mm between image0[0,0,sliceNum0] and image1[0,0,sliceNum1].
    """
    
    mmShift = np.array(image1.TransformIndexToPhysicalPoint([0,0,sliceNum1])) \
            - np.array(image0.TransformIndexToPhysicalPoint([0,0,sliceNum0]))
            
    return mmShift

def get_voxel_shift_bt_slices(
        image0, sliceNum0, image1, sliceNum1, refImage, fractional=False
        ):
    """
    Get the shift in voxels between sliceNum0 in image0 and sliceNum1 in image1
    based on the origins of the two slices (i.e. image0[0,0,sliceNum0] and
    image1[0,0,sliceNum1]).
    
    Parameters
    ----------
    image0 : SimpleITK image
    sliceNum0 : int
        The slice number in image0.
    image1 : SimpleITK image
    sliceNum1 : int
        The slice number in image1.
    refImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert the
        physical shift to pixels.
    fractional : bool, optional (False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    
    Returns
    -------
    voxShift : Numpy array
        Shift in voxels between image0[0,0,sliceNum0] and image1[0,0,sliceNum1].
    """
    
    mmShift = get_phys_shift_bt_slices(image0, sliceNum0, image1, sliceNum1)
    
    voxSpacing = np.array(refImage.GetSpacing())
    
    voxShift = mmShift/voxSpacing
    
    if not fractional:
        voxShift = np.round(voxShift).astype('int')
            
    return voxShift

def shift_frame(frame, voxShift):
    """
    Shift the items in a 2D frame.
    
    Parameters
    ----------
    frame : Numpy array
        A 2D frame from PixelData with shape (1, R, C), where R = number of 
        rows, and C = number of columns.
    voxShift : Numpy array
        Voxel shift, e.g. [di, dj, dz].
        
    Returns
    -------
    shiftedFrame : Numpy array
        Shifted frame with shape (1, R, C). Only the x- and y-components are 
        shifted.
    """
    
    F, R, C = frame.shape
    
    if F > 1:
        msg = "'frame' must be a 2D frame with shape (1, R, C) but has shape"\
              + f" ({F}, {R}, {C})."
        
        raise Exception(msg)
    
    # Initialise shiftedFrame:
    shiftedFrame = np.zeros((1, R, C), dtype='uint')
    #shiftedFrame = np.empty_like(frame, dtype='uint') # this creates 42,932
    # unique values for some reason!
    
    #unique = get_unique_items(Nda=frame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in frame')
    #unique = get_unique_items(Nda=shiftedFrame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in the initialised',
    #      f'shiftedFrame: {unique[:11]}...')
    
    di, dj, dk = voxShift
    
    ##shiftedFrame[0, dj:, di:] = frame[0, :-(1+dj), :-(1+di)]
    ##shiftedFrame[0, :-(1+dj), :-(1+di)] = frame[0, dj:, di:]
    #shiftedFrame[0, :R-dj, :C-di] = frame[0, dj:, di:]
    
    if di > 0 and dj > 0:
        shiftedFrame[0, dj:, di:] = frame[0, :-dj, :-di]
        
    elif di < 0 and dj < 0:
        shiftedFrame[0, :dj, :di] = frame[0, -dj:, -di:]
        
    elif di > 0 and dj < 0:
        shiftedFrame[0, :dj, di:] = frame[0, -dj:, :-di]
        
    elif di < 0 and dj > 0:
        shiftedFrame[0, dj:, :di] = frame[0, :-dj, -di:]
        
    elif di == 0 and dj > 0:
        shiftedFrame[0, dj:, :] = frame[0, :-dj, :]
        
    elif di == 0 and dj < 0:
        shiftedFrame[0, :dj, :] = frame[0, -dj:, :]
        
    elif di > 0 and dj == 0:
        shiftedFrame[0, :, di:] = frame[0, :, :-di]
        
    elif di < 0 and dj == 0:
        shiftedFrame[0, :, :di] = frame[0, :, -di:]
        
    elif di == 0 and dj == 0:
        shiftedFrame[0] = frame[0]
        
    #unique = get_unique_items(Nda=shiftedFrame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in the shiftedFrame',
    #      'after shifting.')
            
    return shiftedFrame

def shift_frames_in_pixarrBySeg_OBSOLETE(PixArrBySeg, F2SindsBySeg, srcIm, 
                                srcSlcNum, trgIm, trgSlcNum, refImage, 
                                shiftInX=True, shiftInY=True, shiftInZ=True, 
                                fractional=False, p2c=False):
    """
    See shift_frames_in_pixarrBySeg in propagate.propagate.py
    
    Shift the items in each frame of each pixel array in a list of 
    pixel-arrays-by-segment.  Example uses: To compensate for differences in 
    z-positions of a single-framed Source pixel array to be "direct" copied to
    a Target slice at a different stack position; or compensating for 
    differences in the origin of Source and Target images for relationship-
    preserving copies of images with equal voxel spacings.
    
    Parameters
    ----------
    PixArrBySeg : list of Numpy arrays
        A list (for each segment) of pixel arrays.
    F2SindsBySeg : list of a list of ints
        List (for each segment) of a list (for each frame) of slice numbers that 
        correspond to each frame in the pixel arrays in PixArrBySeg.
    srcIm : SimpleITK image
        The Source 3D image.
    srcSlcNum : int
        The slice number within the Source image stack that the segment relates
        to.
    trgIm : SimpleITK image
        The Target 3D image.
    trgSlcNum : int
        The slice number within the Target image stack that the segment is to  
        be copied to following the shift in z-components.
    refImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert apply
        the shift in pixels (e.g. trgIm).
    shiftInX / shiftInY / shiftInZ : bool, optional (True by default)
        If True the pixel shift between the origins of srcIm and trgIm
        will be applied along the x / y / z direction.
    fractional : bool, optional (False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    PixArrBySeg : list of Numpy arrays
        As above but with shifted frame with shape (1, R, C).
    F2SindsBySeg : list of a list of ints
        As above but with shifted slice indices that correspond to the shifted 
        PixArrBySeg.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of shift_frames_in_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    # Get pixel shift between FromSliceNum in srcIm and ToSliceNum in 
    # trgIm in the Target image domain:
    voxShift = get_voxel_shift_bt_slices(
        image0=srcIm, sliceNum0=srcSlcNum, image1=trgIm, 
        sliceNum1=trgSlcNum, refImage=trgIm
        )
    
    if p2c:
        print(f'   \nvoxShift between slices = {voxShift}')
        print(f'   \nF2SindsBySeg prior to shifting = {F2SindsBySeg}')
        
        #print(f'\n   Prior to shifting frames:')
        ##print(f'PixArrBySeg = {PixArrBySeg}')
        #print(f'   len(PixArrBySeg) = {len(PixArrBySeg)}')
        #for r in range(len(PixArrBySeg)):
        #    print(f'      type(PixArrBySeg[{r}]) = {type(PixArrBySeg[r])}')
        #    print(f'      PixArrBySeg[{r}].shape = {PixArrBySeg[r].shape}')
    
    if not shiftInX:
        voxShift[0] = 0
    if not shiftInY:
        voxShift[1] = 0
    if not shiftInZ:
        voxShift[2] = 0
    
    if p2c:
        print(f'   \nvoxShift that will be applied = {voxShift}')
        
    for s in range(len(PixArrBySeg)):
        #if PixArrBySeg[s]:
        if PixArrBySeg[s].shape[0]:
            #print('   \nBefore shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            # Replace PixArrBySeg[s] with the result of shifting the in-plane 
            # elements, and replace F2SindsBySeg[s] with [trgSlcNum]:
            PixArrBySeg[s] = shift_frame(
                frame=PixArrBySeg[s], voxShift=voxShift
                )
            
            #print('   \nAfter shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            #F2SindsBySeg[s] = [trgSlcNum] # 11/02
            
            # Shift the frame-to-slice indices by voxShift[2]:
            F2SindsBySeg[s] = [F2SindsBySeg[s][i] + voxShift[2] for i in range(len(F2SindsBySeg[s]))]
            
            if p2c:
                unique = get_unique_items(Items=PixArrBySeg[s], IgnoreZero=False)
                
                print(f'\nThere are {len(unique)} unique items in',
                      f'PixArrBySeg[{s}] after shifting the frame.')
    
    if p2c:
        print(f'   F2SindsBySeg after shifting = {F2SindsBySeg}')
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg

def change_z_inds(indices, newZind):
    """
    Re-assign the z-components of a list of indeces.
    
    As required when making a direct copy of contour points.
    
    Parameters
    ----------
    indices : list of a list of ints
        List (for each point) of a list (for each dimension) of indices, 
        e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    
    Returns
    -------
    indices : list of a list of ints
        List (for each point) of a list (for each dimension) of the indices 
        with modified z-components, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    if not isinstance(newZind, int):
        raise Exception(f"newZind = {newZind} must be an integer.")
    
    for i in range(len(indices)):
        indices[i][2] = newZind
        
    return indices

def z_shift_ptsByCntByRoi(ptsByCntByRoi, newZind, dcmIm):
    """
    14/09/21: It seems for use case 1 or 2a, I previously ran:
        1. ReplaceIndInC2SindsByRoi() (now called change_z_inds())
        2. ZshiftPtsByCntByRoi (now z_shift_ptsByCntByRoi())
        3. ShiftPtsByCntByRoi() (now shift_ptsByCntByRoi()) 
    
    But #3 did what #1-2 did, so it seems there was redundancy (?).
    
    Shift the z-component of the points in ptsByCntByRoi (as required when
    making a direct copy of contour points).
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    dcmIm : SimpleITK Image
        SimpleITK Image representation of the DICOM series.
    
    Returns
    -------
    shiftedPtsByCntByRoi : list of list of a list of a list of floats
        As ptsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    #from ConversionTools import Points2Indices, Indices2Points
    
    # Import the image:
    #dcmIm = import_im(dicomDir)
    
    """ Convert ptsByCntByRoi to indices, modify the z-component of the
    indices, then convert back to physical points. """
    
    #indsByCntByRoi = []
    shiftedPtsByCntByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(ptsByCntByRoi)): 
        #indsByCnt = []
        shiftedPtsByCnt = []
        
        if ptsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(ptsByCntByRoi[r])):
                pts = ptsByCntByRoi[r][c]
                
                inds = pts_to_inds(points=pts, refIm=dcmIm, rounding=False)
                
                # Modify the z-component (k) of the indices to newZind:
                inds = change_z_inds(inds, newZind)
                
                #indsByCnt.append(inds)
                
                # Convert back to physical points:
                pts = inds_to_pts(inds, refIm=dcmIm)
                
                shiftedPtsByCnt.append(pts)
        
        #indsByCntByRoi.append(indsByCnt)
        shiftedPtsByCntByRoi.append(shiftedPtsByCnt)
        
    return shiftedPtsByCntByRoi

def shift_ptsByCntByRoi(
        ptsByCntByRoi, c2sIndsByRoi, voxShift, refIm, p2c=False
        ):
    """
    14/09/21: 
    
    Shift the points in each list of points in each list of contours in each
    list of ROIs given the required voxel shift.
    
    The points in ptsByCntByRoi will be converted to indices, modified based on
    required voxel shift, then convert back points. 
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for each contour) of a list (for each
        point) of a list (for each dimension) of coordinates.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of indices that denote
        the slice number that the contour relates to.
    voxShift : list of ints
        A list (for each dimension) of the required voxel shift.
    refIm : SimpleITK Image
        The image whose image domain the points will relate to (e.g. trgIm for
        copying of points from source to target domain).
    
    Returns
    -------
    shiftedPtsByCntByRoi : list of list of a list of a list of floats
        Modified ptsByCntByRoi.
    shiftedC2SindsByRoi : list of a list of ints
        Modified c2sIndsByRoi.
    
    Notes
    -----
    Example uses: 
        - to compensate for differences in z-positions of a single-contour 
        Source points to be "direct" copied to a Target slice at a different  
        stack position
        - compensating for differences in the origin of Source and Target 
        images for relationship-preserving copies of images with equal voxel 
        spacings
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of shift_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        print('\nPre-shift:')
        print_indsByRoi(c2sIndsByRoi)
        print_ptsByCntByRoi(ptsByCntByRoi)
    
    dim = len(voxShift)
    
    #indsByCntByRoi = []
    shiftedPtsByCntByRoi = []
    shiftedC2SindsByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(ptsByCntByRoi)): 
        #indsByCnt = []
        shiftedPtsByCnt = []
        shiftedC2Sinds = []
        
        # If the list of points-by-contour for this ROI is not empty..
        if ptsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(ptsByCntByRoi[r])):
                c2sInd = c2sIndsByRoi[r][c]
                pts = ptsByCntByRoi[r][c]
                
                # Convert from points to indices:
                #inds = pts_to_inds(points=pts, refIm=refIm, rounding=False)
                listOfInds = pts_to_inds(
                    points=pts, refIm=refIm, rounding=False
                    )
                
                #print(f'\nvoxShift = {voxShift}, \nlistOfInds = {listOfInds}')
                
                # Apply the required pixel shifts:
                #shiftedInds = [
                #    inds[i] + voxShift[i] for i in range(len(voxShift))
                #    ]
                shiftedInds = [
                    [inds[i] + voxShift[i] for i in range(dim)] for inds in listOfInds
                    ]
                
                #indsByCnt.append(inds)
                #print(f'\n\ntype(inds) = {type(inds)}')
                #print(f'\n\ntype(inds[0]) = {type(inds[0])}')
                
                # Convert back to physical points:
                shiftedPts = inds_to_pts(shiftedInds, refIm=refIm)
                
                shiftedPtsByCnt.append(shiftedPts)
                
                # Shift the contour-to-slice indices by voxShift[2]:
                #c2sIndsByRoi[r][c] += voxShift[2]
                shiftedC2Sinds.append(c2sInd + voxShift[2])
        
        #indsByCntByRoi.append(indsByCnt)
        shiftedPtsByCntByRoi.append(shiftedPtsByCnt)
        shiftedC2SindsByRoi.append(shiftedC2Sinds)
        
        if p2c:
            print('\nPost-shift:')
            print_indsByRoi(shiftedC2SindsByRoi)
            print_ptsByCntByRoi(shiftedPtsByCntByRoi)
            print('\nEnd of shift_ptsByCntByRoi\n')
            print('-'*120)
        
    return shiftedPtsByCntByRoi, shiftedC2SindsByRoi

def shift_ptsByCntByRoi_210914(
        ptsByCntByRoi, c2sIndsByRoi, srcIm, srcSlcNum, 
        trgIm, trgSlcNum, refImage, 
        shiftInX=True, shiftInY=True, shiftInZ=True, 
        fractional=False, p2c=False
        ):
    """
    14/09/21: 
    
    Shift the points in each list of points in each list of contours in each
    list of ROIs.  
    
    Example uses: 
        - to compensate for differences in z-positions of a single-contour 
        Source points to be "direct" copied to a Target slice at a different  
        stack position
        - compensating for differences in the origin of Source and Target 
        images for relationship-preserving copies of images with equal voxel 
        spacings
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of indices that denote
        the slice number that the contour relates to.
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    shiftedPtsByCntByRoi : list of list of a list of a list of floats
        As ptsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    #from ConversionTools import Points2Indices, Indices2Points
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of shift_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        
    """ Get voxel shift between FromSliceNum in srcIm and ToSliceNum in 
    trgIm in the Target image domain: """
    voxShift = get_voxel_shift_bt_slices(
        image0=srcIm, sliceNum0=srcSlcNum, 
        image1=trgIm, sliceNum1=trgSlcNum,
        refImage=trgIm
        )
    
    if p2c:
        print(f'   voxShift between slices = {voxShift}')
    
    if not shiftInX:
        voxShift[0] = 0
    if not shiftInY:
        voxShift[1] = 0
    if not shiftInZ:
        voxShift[2] = 0
    
    if p2c:
        print(f'   voxShift that will be applied = {voxShift}')
    
    
    """ Convert ptsByCntByRoi to indices, modify the indices, then convert back
    to physical points. """
    
    #indsByCntByRoi = []
    shiftedPtsByCntByRoi = []
    shiftedC2SindsByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(ptsByCntByRoi)): 
        #indsByCnt = []
        shiftedPtsByCnt = []
        shiftedC2Sinds = []
        
        if ptsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(ptsByCntByRoi[r])):
                pts = ptsByCntByRoi[r][c]
                
                """ inds are the list of relationship-preserving indices that
                correspond to pts: """
                inds = pts_to_inds(points=pts, refIm=trgIm, rounding=False)
                
                #print(f'\n\nvoxShift = {voxShift}, \ninds = {inds}')
                
                """ Apply the required pixel shifts: """
                inds = [inds[i] + voxShift[i] for i in range(len(voxShift))]
                
                #indsByCnt.append(inds)
                #print(f'\n\ntype(inds) = {type(inds)}')
                #print(f'\n\ntype(inds[0]) = {type(inds[0])}')
                
                """ Convert back to physical points: """
                pts = inds_to_pts(inds, refIm=trgIm)
                
                shiftedPtsByCnt.append(pts)
                
                """ Shift the contour-to-slice indices by voxShift[2]: """
                #c2sIndsByRoi[r][c] += voxShift[2]
                shiftedC2Sinds.append(c2sIndsByRoi[r][c] + voxShift[2])
        
        #indsByCntByRoi.append(indsByCnt)
        shiftedPtsByCntByRoi.append(shiftedPtsByCnt)
        shiftedC2SindsByRoi.append(shiftedC2Sinds)
        
    return shiftedPtsByCntByRoi, shiftedC2SindsByRoi

def get_srcF2SindsByRoi_to_trgF2SindsByRoi(
        srcF2SindsByRoi, srcIm, trgIm):
    """
    Shift the indices in srcF2SindsByRoi to account for any differences in the
    origins of srcIm and trgIm.

    Parameters
    ----------
    srcF2SindsByRoi : list of a list of ints
        List (for each segment/ROI) of a list (for each segmentation/contour) 
        of slice numbers that correspond to each Source segmentation/contour.
    srcIm : SimpleITK Image
        The Source 3D image.
    trgIm : SimpleITK Image
        The Target 3D image.

    Returns
    -------
    trgF2SindsByRoi : list of a list of ints
        List (for each segment/ROI) of a list (for each segmentation/contour) 
        of slice numbers that correspond to each Source segmentation/contour in 
        the Target image domain.
    """
    
    srcOrigin = srcIm.GetOrigin()
    trgOrigin = trgIm.GetOrigin()
    
    #dim = len(srcOrigin)
    #mmDiff = [trgOrigin[i] - srcOrigin[i] for i in range(dim)]
    
    dz_mm = trgOrigin[2] - srcOrigin[2]
    
    dk = trgIm.GetSpacing()[2]
    
    dz_pix = round(dz_mm/dk)
    
    trgF2SindsByRoi = []
    
    for r in range(len(srcF2SindsByRoi)):
        trgF2Sinds = []
        
        for c in range(len(srcF2SindsByRoi[r])):
            #trgF2Sinds.append(srcF2SindsByRoi[r][c] + dz_pix)
            trgF2Sinds.append(srcF2SindsByRoi[r][c] - dz_pix)
            
        trgF2SindsByRoi.append(trgF2Sinds)
    
    return trgF2SindsByRoi