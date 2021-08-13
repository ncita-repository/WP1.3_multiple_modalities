# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:34:59 2021

@author: ctorti
"""

""" Medium-level conversions between pixel arrays and indices or points. """

from itertools import product
import numpy as np
from shapely.geometry import Polygon, Point
from scipy.sparse import issparse
from skimage.measure import find_contours

from general_tools.console_printing import print_pts_by_cnt
from general_tools.console_printing import print_pixarr_by_seg
from image_tools.imports import import_im
from conversion_tools.inds_pts import inds_to_pts, pts_to_cntdata, pts_to_inds


def indsByCnt_to_pixarr(indsByCnt, refImage):
    """
    Convert a list of indices grouped by contour to a pixel array.
    
    Parameters
    ----------
    indsByCnt : list of lists of lists of ints
        A list (for each contour) of a list (for all polygons) of a list (for 
        each dimension) of indices.
    refImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
    Returns
    -------
    pixarr : Numpy data array
        A list of 2D masks of the indices - one mask per contour.
    
    Note
    ----
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py
    """
    
    #import numpy as np
    #from shapely.geometry import Polygon, Point
    #from itertools import product
    
    imSize = refImage.GetSize()
    
    Ncontours = len(indsByCnt)
    
    pixarr = np.zeros((Ncontours, imSize[1], imSize[0]))
    
    # Note that indsByCnt are a list of indices for each contour:
    for c in range(Ncontours):
        indices = indsByCnt[c]
        
        # Only proceed if there are at least 3 indices in this contour (the 
        # minimum number required to define a closed contour):
        if len(indices) > 2:
            # Convert list indices to a Shapely Polygon:
            poly = Polygon(indices)
                
            xMin, yMin, xMax, yMax = poly.bounds
         
            points = [Point(x, y) for x, y in
                          product(np.arange(int(xMin), np.ceil(xMax)),
                                  np.arange(int(yMin), np.ceil(yMax)))]
                      
            pointsInPoly = list(filter(poly.contains, points))
            
            for point in pointsInPoly:
                xx, yy = point.xy
                
                x = int(xx[0])
                y = int(yy[0])
                
                if 0 <= y < imSize[1] and 0 <= x < imSize[0]:
                    pixarr[c, y, x] = 1
    
    return pixarr

def pixarr_to_labarr(pixarr, numOfSlices, f2sInds):
    """
    Convert a sparse 3D pixel array (pixarr) to a zero-padded 3D pixel array
    (labarr).
    
    Parameters
    ----------
    pixarr : Numpy array
        SEG PixelData as a Numpy array.
    numOfSlices : int
        The total number of frames in the labelmap (e.g. the number of slices 
        in a image series).
    f2sInds : list of ints
        The DICOM slice numbers that correspond to each frame in pixarr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number.
        
    Returns
    -------
    labarr : Numpy array
        Zero-padded version of pixarr.
    
    Note
    ----
    labarr refers to a 3D pixel array that contains the non-empty frames in
    pixarr, padded with empty frames so that its size matches the size of the 
    3D DICOM image. This is not to be confused with a 'labmap image', which is 
    the SimpleITK image representation of a labmap array.
    """
    
    #import numpy as np
    
    numOfFrames, R, C = pixarr.shape
    
    """ Initialise labarr: """
    #Labmap = np.zeros((numOfSlices, R, C), dtype='bool')
    #Labmap = np.zeros((numOfSlices, R, C), dtype='uint')
    labarr = np.zeros((numOfSlices, R, C), dtype=pixarr.dtype)
    
    """ Note: numOfFrames is the number of non-zero frames that make up pixarr,
    whereas numOfSlices are the total frame count in 3D DICOM image. """
    
    #print(f'\npixarr.shape = {pixarr.shape}')
    ##print(f'\n\nThe maximum value in pixarr is {pixarr.max()}')
    #print(f'f2sInds = {f2sInds}')
    #print(f'labarr.shape = {labarr.shape}')
    
    if numOfFrames:
        """ Over-write all frames in labarr with non-zero 2D masks from 
        pixarr.  """
        for i in range(len(f2sInds)):
            """ The slice number for the i^th frame """
            s = f2sInds[i]
            
            #print(f'frame {i} corresponds to slice {s}')
            
            #labarr[s] = pixarr[i] # 13/02 This was over-writting frames for 
            # pixel arrays containing more than one frame for the same slice!!!
            
            labarr[s] += pixarr[i]
    
    #print(f'\n\nThe maximum value in labarr is {labarr.max()}')
    
    return labarr

def get_labarrByRoi(rts, dicomDir, c2sIndsByRoi, refIm):
    """ 
    Get a list of zero-padded 3D Numpy arrays - one per ROI.
    
    Parameters
    ----------
    rts : Pydicom object
        RTS object.
    dicomDir : str
        Directory containing the corresponding DICOMs.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of slice numbers  
        that corresponds to each contour.
    refIm : SimpleITK image
        The 3D image that pixarr relates to.
    
    Returns
    -------
    labarrByRoi : list of Numpy arrays
        A list of 3D zero-padded labelmap Numpy arrays - one per contour.
    
    Note
    ----
    While a pixel array (pixarr) contains only non-zero arrays, a "label
    array" (labarr) is a zero-padded version of pixarr, containing the same
    number of frames as in the SimpleITK image representation.
    """
    
    from RtsTools import GetIndsByRoi
    
    labarrByRoi = []
    
    indsByRoi = GetIndsByRoi(rts, dicomDir)
    
    for roiNum in range(len(indsByRoi)):
        pixarr = indsByCnt_to_pixarr(indsByCnt=indsByRoi[roiNum], 
                                     refImage=refIm)
    
        labarr = pixarr_to_labarr(pixarr=pixarr, 
                                  numOfSlices=refIm.GetSize()[2], 
                                  f2sInds=c2sIndsByRoi[roiNum])
        
        labarrByRoi.append(labarr)
    
    return labarrByRoi

def ptsByCnt_to_pixarr(ptsByCnt, refIm, p2c=False):
    """ 
    Convert a list of points by contour to a pixel array.
    
    Parameters
    ----------
    ptsByCnt : list of a list of a list of floats
        List (for all contours) of a list (for each point) of a list (for each 
        dimension) of coordinates.
    refIm : SimpleITK image
        The 3D image that ptsByCnt relate to.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarr : Numpy data array
        A pixel array containing as many frames as the number of contours in
        ptsByCnt.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of ptsByCnt_to_pixarr():')
        print('\n\n', '-'*120)
        
    indsByCnt = []
    
    if ptsByCnt == []:
        import numpy as np
        
        C, R, F = refIm.GetSize()
        
        #pixarr = np.zeros((F, R, C), dtype='uint')
        pixarr = np.zeros((F, R, C), dtype=np.uint8)
        
    else:
        for Pts in ptsByCnt:
            inds = pts_to_inds(Pts, refIm)
                
            indsByCnt.append(inds)
    
        """ Convert indsByCnt to a pixel array. """
        pixarr = indsByCnt_to_pixarr(indsByCnt, refIm)
    
    if p2c:
        print(f'   pixarr.shape = {pixarr.shape}')
        print('-'*120)
    
    return pixarr

def ptsByCntByRoi_to_pixarrByRoi(ptsByCntByRoi, refIm, p2c=False):
    """
    Convert a list of points in all contours in all ROIs to a list of pixel
    arrays grouped by ROI.
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    refIm : SimpleITK image
        The 3D image that relates to ptsByCntByRoi.
        
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
                       
    Returns
    -------
    pixarrByRoi : list of Numpy data arrays
        A list (for each ROI) of 3D pixel arrays. Each pixel array contains as
        many frames as the number of contours in the list of points-by-contour
        for that given ROI.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of ptsByCntByRoi_to_pixarrByRoi():')
        print('\n\n', '-'*120)
        
    pixarrByRoi = []
    
    R = len(ptsByCntByRoi)
    
    for r in range(R):        
        """ Convert to a pixel array. """
        pixarr = ptsByCnt_to_pixarr(ptsByCntByRoi[r], refIm, p2c)
        
        pixarrByRoi.append(pixarr)
        
    if p2c:
        shapes = [pixarrByRoi[r].shape for r in range(R)]
        print(f'   Shape of each pixarrByRoi = {shapes}')
        print('-'*120)
    
    return pixarrByRoi

def pixarr_to_indsByFrame(pixarr, f2sInds, thresh=0.5, p2c=False):
    """  
    Convert a 3D pixel array to a list (for each frame) of a list of indices
    that define the equivalent contour for the mask in each frame.
 
    Parameters
    ----------
    pixarr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        pixarr.
    f2sInds : list of list of ints
        A list (for each frame) of the slice numbers that correspond to each 
        frame in pixarr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    thresh : float, optional (0.5 by default)
        Threshold value used to binarise the labels in pixarr.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    indsByObjByFrame : list of a list of a list of floats
        A list (for each frame) of a list (for each contour) of a list of 
        [x, y, z] coordinates of the polygons converted from each frame in
        pixarr.
    
    Note
    ----
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py
    """
    
    #import numpy as np
    #from scipy.sparse import issparse
    #from skimage.measure import find_contours
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of pixarr_to_indsByFrame()')
        print('\n\n', '-'*120)
        
    indsByObjByFrame = []
    
    for frameNum, frame in enumerate(pixarr):
        if issparse(frame):
            frame = np.array(frame.astype('byte').todense())
 
        if (frame != 0).sum() == 0:
            """ If frame is empty, just skip it """
            continue
 
        """ Add an empty row and column at both sides of the frame to ensure 
        that any segmentation pixels near an edge are found. """
        expandedDims = (frame.shape[0] + 2, frame.shape[1] + 2)
        
        expandedFrame = np.zeros(expandedDims, dtype=float)
        
        expandedFrame[1:frame.shape[0] + 1, 1:frame.shape[1] + 1] = frame
 
        """ Get the indices that define contours in expandedFrame. """
        indsByObj_2D = find_contours(expandedFrame.T, thresh)
        
        if p2c:
            print(f'frameNum = {frameNum}')
            print(f'len(indsByObj_2D) = {len(indsByObj_2D)}')
            #print(f'indsByObj_2D = {indsByObj_2D}')
 
        """ Remove one row and column to shift the indices back to their 
        original locations. """
        indsByObj_2D = [np.subtract(x, 1).tolist() for x in indsByObj_2D]
        
        if p2c:
            #print(f'\n\n\nframeNum = {frameNum}')
            print(f'len(indsByObj_2D) = {len(indsByObj_2D)}')
            #print(f'indsByObj_2D = {indsByObj_2D}')
        
        sliceNum = float(f2sInds[frameNum])
            
        """ Add the index for the 3rd dimension. """
        indsByObj_3D = []
        
        for indsThisObj_2D in indsByObj_2D:
            indsThisObj_3D = [Ind_2D + [sliceNum] for Ind_2D in indsThisObj_2D]
            
            indsByObj_3D.append(indsThisObj_3D)
        
        if p2c:
            print(f'len(indsByObj_3D) = {len(indsByObj_3D)}')
            #print(f'\nindsByObj_3D = {indsByObj_3D}')
        
        indsByObjByFrame.append(indsByObj_3D)
        
    
    if p2c:
        print(f'\nlen(indsByObjByFrame) = {len(indsByObjByFrame)}')
        #print(f'\nindsByObjByFrame = {indsByObjByFrame}')
 
    return indsByObjByFrame

def pixarr_to_listOfInds(pixarr, f2sInds):
    """  
    Convert a binary pixel array to a list (for each of a list of indices of 
    non-zero elements).
    
    If pixarr is 2D, the function will return a list of length 1 of a list of
    indices.  If pixarr is 3D, the function will return a list of length N of a
    list of indices, where N is the size of pixarr along the third dimension. 
    
    Parameters
    ----------
    pixarr : Numpy array
        A binary pixel array representing a segmentation (2D) or segment (3D).
    f2sInds : List of ints
        The DICOM slice numbers that correspond to each frame in pixarr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number.
    
    Returns
    -------
    listOfInds : list of a list of ints
        A list (of length N, where N = number of frames in pixarr) of a list
        (for each dimension) of a list of integer coordinates for each non-zero
        element in each 2D frame in pixarr.
    """
    
    #import numpy as np
    
    dim = pixarr.ndim
    
    if dim == 2:
        R, C = pixarr.shape
        F = 1
    elif dim == 3:
        F, R, C = pixarr.shape
    else:
        msg = "'pixarr' must have 2 or 3 dimensions."
        raise Exception(msg)
    
    listOfInds = []
    
    for f in range(F):       
        yInds, xInds = np.where(pixarr[f] > 0)
        
        xInds = xInds.tolist()
        yInds = yInds.tolist()
        
        if dim == 2: 
            inds = [[yInds[i], xInds[i]] for i in range(len(xInds))]
        else:
            s = f2sInds[f]
            
            inds = [[yInds[i], xInds[i], s] for i in range(len(xInds))]
        
        listOfInds.append(inds)
        
    return listOfInds

def pixarr_to_ptsByCnt(pixarr, f2sInds, dicomDir, thresh=0.5, 
                       p2c=False):
    """
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates.
 
    Parameters
    ----------
    pixarr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        pixarr.
    f2sInds : list of list of ints
        A list (for each frame) of the slice numbers that correspond to each 
        frame in pixarr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    dicomDir : str
        Directory containing the DICOMs that relate to pixarr.
    thresh : float, optional (0.5 by default)
        Threshold value used to binarise the labels in pixarr.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    ptsByCnt : list of a list of floats
        A list (for each contour) of a list (for each dimension) of the 
        physical coordinates that define each distinct object in each frame of
        pixarr.
    cntdataByCnt : list of a list of strs
        A list (for each contour) of a flattened list of [x, y, z] physical 
        coordinates in ptsByCnt.
    c2sInds : list of ints
        List (for each contour) of slice numbers that correspond to each 
        contour in ptsByCnt.
    """
    
    #from general_tools.console_printing import print_pts_by_cnt
    #from general_tools.console_printing import print_pixarr_by_seg
    #from image_tools.imports import import_im
    #from conversion_tools.inds_pts import inds_to_pts, pts_to_cntdata
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of pixarr_to_ptsByCnt():')
        print('\n\n', '-'*120)
        print(f'Prior to conversion: \nf2sInds = {f2sInds}')
        print_pixarr_by_seg(pixarr)
        
    refIm = import_im(dicomDir)
    
    # Convert the pixel array to a list of indices-by-object-by-frame:
    """ A list of objects originating from the same frame will be treated as 
    distinct contours. """ 
    indsByObjByFrame = pixarr_to_indsByFrame(pixarr, f2sInds, thresh=0.5)
    
    if p2c:
        F = len(indsByObjByFrame)
        print(f'\nlen(indsByObjByFrame) = {F} frames')
        for f in range(F):
            O = len(indsByObjByFrame[f])
            print(f'   len(indsByObjByFrame[{f}]) = {O} objects')
            for o in range(O):
                N = len(indsByObjByFrame[f][o])
                print(f'      len(indsByObjByFrame[{f}][{o}]) = {N} indices')
    
    ptsByCnt = []
    cntdataByCnt = []
    c2sInds = []
    
    for f in range(len(indsByObjByFrame)):
        for o in range(len(indsByObjByFrame[f])):
            pts = inds_to_pts(inds=indsByObjByFrame[f][o], refIm=refIm)
            
            ptsByCnt.append(pts)
            
            cntdataByCnt.append(pts_to_cntdata(points=pts))
            
            c2sInds.append(f2sInds[f])
    
    #print(f'\n\n\nlen(ptsByCnt) = {len(ptsByCnt)}')
    #print(f'\nlen(cntdataByCnt) = {len(cntdataByCnt)}')
    #print(f'\nlen(c2sInds) = {len(c2sInds)}')
    #print(f'\nptsByCnt = {ptsByCnt}')
    #print(f'\ncntdataByCnt = {cntdataByCnt}')
    #print(f'\nc2sInds = {c2sInds}')
    
    if p2c:
        print(f'\nAfter conversion: \nc2sInds = {c2sInds}')
        print_pts_by_cnt(ptsByCnt)
        print('-'*120)
 
    return ptsByCnt, cntdataByCnt, c2sInds

def pixarrByRoi_to_ptsByCntByRoi(pixarrByRoi, f2sIndsByRoi, dicomDir, 
                                 thresh=0.5, p2c=False):
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
    dicomDir : str
        Directory containing the DICOMs that relate to pixarr.
    thresh : float, optional (0.5 by default)
        Threshold value used to binarise the labels in pixarr.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
 
    Returns
    -------
    ptsByCntByRoi : list of a list of a list of floats
        A list (for each ROI) of a list (for each contour) of a list (for each 
        dimension) of the physical coordinates that define each pixel array in
        pixarrByRoi.
    cntdataByObjByFrame : list of a list of a list of str
        A list (for each ROI) of a list (for each contour) of a flattened list 
        of [x, y, z] physical coordinates that define each pixel array in
        pixarrByRoi.
    c2sIndsByRoi : list of a list of ints
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour in ptsByCntByRoi.
    """
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of pixarrByRoi_to_ptsByCntByRoi():')
        print('\n\n', '-'*120)
        
    ptsByCntByRoi = []
    cntdataByCntByRoi = []
    c2sIndsByRoi = []
    
    #print(f'\n\n\npixarrByRoi = {pixarrByRoi}')
    #print(f'len(pixarrByRoi) = {len(pixarrByRoi)}')
    #for r in range(len(pixarrByRoi)):
    #    print(f'type(pixarrByRoi[{r}]) = {type(pixarrByRoi[r])}')
    #    print(f'pixarrByRoi[{r}].shape = {pixarrByRoi[r].shape}')
    
    for r in range(len(pixarrByRoi)):
        #if pixarrByRoi[r]:
        if pixarrByRoi[r].shape[0]:
            ptsByCnt, cntdataByCnt, c2sInds\
                = pixarr_to_ptsByCnt(pixarrByRoi[r], f2sIndsByRoi[r], 
                                     dicomDir, thresh, p2c)
            
            ptsByCntByRoi.append(ptsByCnt)
            cntdataByCntByRoi.append(cntdataByCnt)
            c2sIndsByRoi.append(c2sInds)
    
    #print(f'\n\n\nlen(ptsByCntByRoi) = {len(ptsByCntByRoi)}')
    #print(f'\nlen(cntdataByCntByRoi) = {len(cntdataByCntByRoi)}')
    #print(f'\nlen(c2sIndsByRoi) = {len(c2sIndsByRoi)}')
    #print(f'\nptsByCntByRoi = {ptsByCntByRoi}')
    #print(f'\ncntdataByCntByRoi = {cntdataByCntByRoi}')
    #print(f'\nc2sIndsByRoi = {c2sIndsByRoi}')
    
    if p2c:
        print('-'*120)
 
    return ptsByCntByRoi, cntdataByCntByRoi, c2sIndsByRoi

def inds_to_mask(inds, refImage):
    """
    Convert a list of indices to a (filled) mask (pixel array).
       
    Parameters
    ----------
    inds : list of a list of ints
        A list (for all points) of a list (for each dimension) of indices.
    refImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
    Returns
    -------
    mask : Numpy data array
        A 2D mask of the indices.
    
    Note
    ----
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py
    """
    
    # Only proceed if there are at least 3 indices in this contour (the 
    # minimum number required to define a closed contour):
    if len(inds) < 3:
        raise Exception(f"There are only {len(inds)} indices. A minimum of ",
                        "3 are required to define a closed contour.")
        
    C, R, S = refImage.GetSize()
    
    mask = np.zeros((1, R, C))
    
    #print(f'\nmask.shape = {mask.shape}')
    
    # Convert inds to a Shapely Polygon:
    poly = Polygon(inds)
    
    xMin, yMin, xMax, yMax = poly.bounds
 
    points = [Point(x, y) for x, y in 
              product(np.arange(int(xMin), np.ceil(xMax)),
                      np.arange(int(yMin), np.ceil(yMax)))]
              
    pointsInPoly = list(filter(poly.contains, points))
    
    for point in pointsInPoly:
        xx, yy = point.xy
        
        x = int(xx[0])
        y = int(yy[0])
        
        if 0 <= y < R and 0 <= x < C:
            #print(f'\n   mask.shape = {mask.shape}')
            
            mask[0, y, x] = 1
    
    return mask