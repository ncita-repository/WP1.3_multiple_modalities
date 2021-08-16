# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:24:52 2021

@author: ctorti
"""

import numpy as np
from general_tools.general import get_unique_items
from image_tools.imports import import_im
from conversion_tools.inds_pts import pts_to_inds, inds_to_pts

def replace_ind_in_C2SindsByRoi_OBSOLETE(C2SindsByRoi, IndToReplace, ReplacementInd):
    # See replace_ind_in_f2sIndsByRoi in propagate.propagate.py
    for r in range(len(C2SindsByRoi)):
        for c in range(len(C2SindsByRoi[r])):
            if C2SindsByRoi[r][c] == IndToReplace:
                C2SindsByRoi[r][c] = ReplacementInd
    
    return C2SindsByRoi

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

def get_voxel_shift_bt_slices(image0, sliceNum0, image1, sliceNum1, refImage,
                              fractional=False):
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

def shift_frames_in_pixarrBySeg_OBSOLETE(PixArrBySeg, F2SindsBySeg, SrcImage, 
                                SrcSliceNum, TrgImage, TrgSliceNum, refImage, 
                                ShiftInX=True, ShiftInY=True, ShiftInZ=True, 
                                fractional=False, LogToConsole=False):
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
    SrcImage : SimpleITK image
        The Source 3D image.
    SrcSliceNum : int
        The slice number within the Source image stack that the segment relates
        to.
    TrgImage : SimpleITK image
        The Target 3D image.
    TrgSliceNum : int
        The slice number within the Target image stack that the segment is to  
        be copied to following the shift in z-components.
    refImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert apply
        the shift in pixels (e.g. TrgImage).
    ShiftInX / ShiftInY / ShiftInZ : bool, optional (True by default)
        If True the pixel shift between the origins of SrcImage and TrgImage
        will be applied along the x / y / z direction.
    fractional : bool, optional (False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    PixArrBySeg : list of Numpy arrays
        As above but with shifted frame with shape (1, R, C).
    F2SindsBySeg : list of a list of ints
        As above but with shifted slice indices that correspond to the shifted 
        PixArrBySeg.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of shift_frames_in_pixarrBySeg():')
        print('\n\n', '-'*120)
        
    # Get pixel shift between FromSliceNum in SrcImage and ToSliceNum in 
    # TrgImage in the Target image domain:
    voxShift = get_voxel_shift_bt_slices(
        image0=SrcImage, sliceNum0=SrcSliceNum, image1=TrgImage, 
        sliceNum1=TrgSliceNum, refImage=TrgImage
        )
    
    if LogToConsole:
        print(f'   \npixShift between slices = {voxShift}')
        print(f'   \nF2SindsBySeg prior to shifting = {F2SindsBySeg}')
        
        #print(f'\n   Prior to shifting frames:')
        ##print(f'PixArrBySeg = {PixArrBySeg}')
        #print(f'   len(PixArrBySeg) = {len(PixArrBySeg)}')
        #for r in range(len(PixArrBySeg)):
        #    print(f'      type(PixArrBySeg[{r}]) = {type(PixArrBySeg[r])}')
        #    print(f'      PixArrBySeg[{r}].shape = {PixArrBySeg[r].shape}')
    
    if not ShiftInX:
        voxShift[0] = 0
    if not ShiftInY:
        voxShift[1] = 0
    if not ShiftInZ:
        voxShift[2] = 0
    
    if LogToConsole:
        print(f'   \nvoxShift that will be applied = {voxShift}')
        
    for s in range(len(PixArrBySeg)):
        #if PixArrBySeg[s]:
        if PixArrBySeg[s].shape[0]:
            #print('   \nBefore shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            # Replace PixArrBySeg[s] with the result of shifting the in-plane 
            # elements, and replace F2SindsBySeg[s] with [TrgSliceNum]:
            PixArrBySeg[s] = shift_frame(
                frame=PixArrBySeg[s], voxShift=voxShift
                )
            
            #print('   \nAfter shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            #F2SindsBySeg[s] = [TrgSliceNum] # 11/02
            
            # Shift the frame-to-slice indices by voxShift[2]:
            F2SindsBySeg[s] = [F2SindsBySeg[s][i] + voxShift[2] for i in range(len(F2SindsBySeg[s]))]
            
            if LogToConsole:
                unique = get_unique_items(Items=PixArrBySeg[s], IgnoreZero=False)
                
                print(f'\nThere are {len(unique)} unique items in',
                      f'PixArrBySeg[{s}] after shifting the frame.')
    
    if LogToConsole:
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

def z_shift_ptsByCntByRoi(PtsByCntByRoi, newZind, dicomDir):
    """
    Shift the z-component of the points in PtsByCntByRoi (as required when
    making a direct copy of contour points).
    
    Parameters
    ----------
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    NewPtsByCntByRoi : list of list of a list of a list of floats
        As PtsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    #from ConversionTools import Points2Indices, Indices2Points
    
    # Import the image:
    Im = import_im(dicomDir)
    
    """ Convert PtsByCntByRoi to indices, modify the z-component of the
    indices, then convert back to physical points. """
    
    #IndsByCntByRoi = []
    NewPtsByCntByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(PtsByCntByRoi)): 
        #IndsByCnt = []
        NewPtsByCnt = []
        
        if PtsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(PtsByCntByRoi[r])):
                pts = PtsByCntByRoi[r][c]
                
                inds = pts_to_inds(Points=pts, RefIm=Im, Rounding=False)
                
                # Modify the z-component (k) of the indices to newZind:
                inds = change_z_inds(inds, newZind)
                
                #IndsByCnt.append(inds)
                
                # Convert back to physical points:
                pts = inds_to_pts(inds, RefIm=Im)
                
                NewPtsByCnt.append(pts)
        
        #IndsByCntByRoi.append(IndsByCnt)
        NewPtsByCntByRoi.append(NewPtsByCnt)
        
    return NewPtsByCntByRoi

def shift_ptsByCntByRoi(PtsByCntByRoi, C2SindsByRoi, SrcImage, SrcSliceNum, 
                        TrgImage, TrgSliceNum, refImage, ShiftInX=True,  
                        ShiftInY=True, ShiftInZ=True, fractional=False, 
                        LogToConsole=False):
    """
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
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    newZind : int
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    dicomDir : str
        Directory containing the corresponding DICOMs.
    
    Returns
    -------
    NewPtsByCntByRoi : list of list of a list of a list of floats
        As PtsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    #from ConversionTools import Points2Indices, Indices2Points
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of shift_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        
    """ Get pixel shift between FromSliceNum in SrcImage and ToSliceNum in 
    TrgImage in the Target image domain: """
    pixShift = get_pix_shift_bt_slices(image0=SrcImage, 
                                       sliceNum0=SrcSliceNum, 
                                       image1=TrgImage, 
                                       sliceNum1=TrgSliceNum,
                                       refImage=TrgImage)
    
    if LogToConsole:
        print(f'   pixShift between slices = {pixShift}')
    
    if not ShiftInX:
        pixShift[0] = 0
    if not ShiftInY:
        pixShift[1] = 0
    if not ShiftInZ:
        pixShift[2] = 0
    
    if LogToConsole:
        print(f'   pixShift that will be applied = {pixShift}')
    
    
    """ Convert PtsByCntByRoi to indices, modify the indices, then convert back
    to physical points. """
    
    #IndsByCntByRoi = []
    NewPtsByCntByRoi = []
    NewC2SindsByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(PtsByCntByRoi)): 
        #IndsByCnt = []
        NewPtsByCnt = []
        NewC2Sinds = []
        
        if PtsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(PtsByCntByRoi[r])):
                Pts = PtsByCntByRoi[r][c]
                
                """ Inds are the list of relationship-preserving indices that
                correspond to Pts: """
                Inds = pts_to_inds(Points=Pts, RefIm=TrgImage, Rounding=False)
                
                #print(f'\n\npixShift = {pixShift}, \nInds = {Inds}')
                
                """ Apply the required pixel shifts: """
                Inds = [Inds[i] + pixShift[i] for i in range(len(pixShift))]
                
                #IndsByCnt.append(inds)
                #print(f'\n\ntype(Inds) = {type(Inds)}')
                #print(f'\n\ntype(Inds[0]) = {type(Inds[0])}')
                
                """ Convert back to physical points: """
                Pts = inds_to_pts(Inds, RefIm=TrgImage)
                
                NewPtsByCnt.append(Pts)
                
                """ Shift the contour-to-slice indices by pixShift[2]: """
                #C2SindsByRoi[r][c] += pixShift[2]
                NewC2Sinds.append(C2SindsByRoi[r][c] + pixShift[2])
        
        #IndsByCntByRoi.append(IndsByCnt)
        NewPtsByCntByRoi.append(NewPtsByCnt)
        NewC2SindsByRoi.append(NewC2Sinds)
        
    return NewPtsByCntByRoi, NewC2SindsByRoi

def get_srcF2SindsByRoi_to_trgF2SindsByRoi(
        SrcF2SindsByRoi, SrcImage, TrgImage):
    """
    Shift the indices in SrcF2SindsByRoi to account for any differences in the
    origins of SrcIm and TrgIm.

    Parameters
    ----------
    SrcF2SindsByRoi : list of a list of ints
        List (for each segment/ROI) of a list (for each segmentation/contour) 
        of slice numbers that correspond to each Source segmentation/contour.
    SrcImage : SimpleITK Image
        The Source 3D image.
    TrgImage : SimpleITK Image
        The Target 3D image.

    Returns
    -------
    TrgF2SindsByRoi : list of a list of ints
        List (for each segment/ROI) of a list (for each segmentation/contour) 
        of slice numbers that correspond to each Source segmentation/contour in 
        the Target image domain.
    """
    
    SrcOrigin = SrcImage.GetOrigin()
    TrgOrigin = TrgImage.GetOrigin()
    
    #dim = len(SrcOrigin)
    #mmDiff = [TrgOrigin[i] - SrcOrigin[i] for i in range(dim)]
    
    Dz_mm = TrgOrigin[2] - SrcOrigin[2]
    
    dk = TrgImage.GetSpacing()[2]
    
    Dz_pix = round(Dz_mm/dk)
    
    TrgF2SindsByRoi = []
    
    for r in range(len(SrcF2SindsByRoi)):
        TrgF2Sinds = []
        
        for c in range(len(SrcF2SindsByRoi[r])):
            #TrgF2Sinds.append(SrcF2SindsByRoi[r][c] + Dz_pix)
            TrgF2Sinds.append(SrcF2SindsByRoi[r][c] - Dz_pix)
            
        TrgF2SindsByRoi.append(TrgF2Sinds)
    
    return TrgF2SindsByRoi