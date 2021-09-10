# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:03:04 2021

@author: ctorti
"""


from shapely.geometry import Point, MultiPoint
#from image_tools.imports import import_im
from io_tools.imports import import_im
from conversion_tools.inds_pts_pixarrs import pixarr_to_ptsByCnt
from conversion_tools.inds_pts_pixarrs import pixarrByRoi_to_ptsByCntByRoi


def is_pt_in_poly(point, vertices):
    """
    Determine if a point lies within a polygon.
    
    Parameters
    ----------
    point : list of floats
        A list of length 3 of floats representing a 3D point.
    vertices : list of floats
        A list of length 8 of floats for all vertices that define a 3D volume.
        
    Returns
    -------
    pt_in_poly : bool
        True/False if Point is/isn't in Polygon.
    
    Note
    ----
    This hasn't been tested on 2D points and 2D polygons.
    """
    
    ##from shapely.geometry import Polygon, Point
    #from shapely.geometry import Point, MultiPoint
    
    """ Use of Polygon requires that the points be ordered specifically. Use
    of MultiPoint and the convex hull allows for the points to be in a random
    order. """
    
    #polygon = Polygon(vertices)
    polygon = MultiPoint(vertices).convex_hull
    
    point = Point(point)
    
    pt_in_poly = polygon.contains(point)
    
    return pt_in_poly

def get_im_verts(image):
    """
    Get the vertices that define the 3D extent of a 3D image.
    
    Parameters
    ----------
    image : SimpleITK image
        A 3D image whose physical extent is to be determined.
    
    Returns
    -------
    pts : list of floats
        A list of length 8 of floats for all vertices that define the 3D extent
        of image.
    """
    
    w = image.GetWidth()
    h = image.GetHeight()
    d = image.GetDepth()
    
    pts = [image.TransformIndexToPhysicalPoint((0, 0, 0)), 
           image.TransformIndexToPhysicalPoint((w, 0, 0)),
           image.TransformIndexToPhysicalPoint((w, h, 0)),
           image.TransformIndexToPhysicalPoint((0, h, 0)),
           image.TransformIndexToPhysicalPoint((0, h, d)),
           image.TransformIndexToPhysicalPoint((w, h, d)),
           image.TransformIndexToPhysicalPoint((w, 0, d)),
           image.TransformIndexToPhysicalPoint((0, 0, d))
           ]
    
    return pts

def get_im_extent(image):
    """
    Get the extent of a 3D image.
    
    Parameters
    ----------
    image : SimpleITK image
        A 3D image whose physical extent is to be determined.
    
    Returns
    -------
    extent : list of floats
        A list (for each dimension) of the max-to-min span along the direction.
    """
    
    sizes = image.GetSize()
    spacings = image.GetSpacing()
    
    extent = [sizes[i]*spacings[i] for i in range(len(sizes))]
    
    return extent

def get_slice_verts(image, slcNum):
    """
    Get the vertices that define the 3D extent of a slice in a 3D image. 
    
    Parameters
    ----------
    image : SimpleITK image
        A 3D image.
    slcNum : int
        The slice index of interest.
    
    Returns
    -------
    pts : list of floats
        A list of length 8 of floats for all vertices that define the 3D extent
        of the slice.
    """
    
    W = image.GetWidth()
    H = image.GetHeight()
    #D = image.GetDepth()
    
    pts = [image.TransformContinuousIndexToPhysicalPoint((0, 0, slcNum - 0.5)), 
           image.TransformContinuousIndexToPhysicalPoint((W, 0, slcNum - 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((W, H, slcNum - 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((0, H, slcNum - 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((0, H, slcNum + 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((W, H, slcNum + 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((W, 0, slcNum + 0.5)),
           image.TransformContinuousIndexToPhysicalPoint((0, 0, slcNum + 0.5))
           ]
    
    return pts

def prop_of_cnt_in_extent(points, trgIm, p2c=False):
    """
    Determine what proportion of points that make up a contour intersect the
    volume of a 3D image.
    
    Parameters
    ----------
    points : list of a list of floats
        A list (for each point) of a list (for each dimension) of coordinates, 
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    trgIm : SimpleITK Image
        The image whose grid/extent is to be checked for intersection with the
        points of interest (e.g. trgIm if copying contours from source to
        target).
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
        
    Returns
    -------
    fracProp : float
        The fractional proportion (normalised to 1) of the points in Contour 
        that intersect the extent of RefImage.
    """
    
    verts = get_im_verts(trgIm)
    
    isInside = []
    
    for pt in points:
        isInside.append(is_pt_in_poly(pt, verts))
    
    fracProp = sum(isInside)/len(isInside)
    
    if p2c:
        print('\n\n\nFunction prop_of_cnt_in_extent():')
        print(f'   Contour has {len(points)} points')
        
        print('\n   The vertices of the image grid:')
        for i in range(len(verts)):
            print(f'   {verts[i]}')
        N = len(isInside)
        limit = 10
        if N > limit:
            print(f'\n   The first {limit} points:')
            for i in range(limit):
                print(f'   {points[i]}')
        else:
            print(f'\n   The first {N} points:')
            for i in range(N):
                print(f'   {points[i]}')
        
        if fracProp < 1:
            print('\n   The vertices of the image grid:')
            for i in range(len(verts)):
                print(f'   {verts[i]}')
            N = len(isInside)
            limit = 3
            if N > limit:
                print(f'\n   The first {limit} points that lie outside of extent:')
                for i in range(limit):
                    print(f'   {points[i]}')
            else:
                print(f'\n   The first {N} points that lie outside of extent:')
                for i in range(N):
                    print(f'   {points[i]}')
                    
    return fracProp

def prop_of_rois_in_extent(ptsByCntByRoi, trgIm, p2c=False):
    """
    Determine what proportion of the points in a list of points by contour by
    ROI intersect with the volume of a 3D image.
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    trgIm : SimpleITK Image
        The image whose grid/extent is to be checked for intersection with the
        points of interest (e.g. trgIm if copying contours from source to
        target).
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
        
    Returns
    -------
    fracProp : float
        The fractional proportion (normalised to 1) of ptsByCntByRoi that
        intersect the extent of trgIm.
    """
    
    verts = get_im_verts(trgIm)
    
    isInside = []
    
    # Iterate through ROIs:
    for r in range(len(ptsByCntByRoi)):
        # Iterate through contours:
        for c in range(len(ptsByCntByRoi[r])):
            # Iterate through points:
            for p in range(len(ptsByCntByRoi[r][c])):
                isInside.append(is_pt_in_poly(ptsByCntByRoi[r][c][p], verts))
    
    fracProp = sum(isInside)/len(isInside)
                    
    return fracProp

def prop_of_pixarr_in_extent(pixarr, f2sInds, srcIm, trgIm, p2c=False):
    """
    Determine what proportion of voxels that make up a 3D pixel array intersect
    with the volume of a 3D image.
    
    Parameters
    ----------
    pixarr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    f2sInds : list of list of ints
        A list (for each frame) of the slice numbers that correspond to each 
        frame in pixarr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    srcIm : SimpleITK Image
        The image that relates to pixarr.
    trgIm : SimpleITK Image
        The image whose grid/extent is to be checked for intersection with the
        points of interest (e.g. trgIm if copying contours from source to
        target).
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
        
    Returns
    -------
    fracProp : float
        The fractional proportion (normalised to 1) of voxels in pixarr that 
        intersect the extent of RefImage.
    """
    
    #from conversion_tools.inds_pts_pixarrs import pixarr_to_ptsByCnt
    #from image_tools.imports import import_im
    
    F, R, C = pixarr.shape
    
    # Convert the pixel array to a list (for each frame in pixarr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    ptsByObjByFrame, cntdataByObjByFrame\
        = pixarr_to_ptsByCnt(pixarr, f2sInds, srcIm)
    
    verts = get_im_verts(trgIm)
    
    isInside = []
    
    for f in range(len(ptsByObjByFrame)):
        for o in range(len(ptsByObjByFrame[f])):
            for pt in ptsByObjByFrame[f][o]:
                isInside.append(is_pt_in_poly(pt, verts))
    
    fracProp = sum(isInside)/len(isInside)
    
    if p2c:
        print('\n\n\nFunction ProportionOfPixArrInExtent():')
        print(f'   pixarr has {F}x{R}x{C} (FxRxC)')
        print(f'   f2sInds = {f2sInds}')
        print(f'   len(ptsByObjByFrame) = {len(ptsByObjByFrame)}')
        print(f'   len(cntdataByObjByFrame) = {len(cntdataByObjByFrame)}')
        for f in range(len(ptsByObjByFrame)):
            print(f'   Frame {f} has {len(ptsByObjByFrame[f])} objects')
            for o in range(len(ptsByObjByFrame[f])):
                print(f'      Object {o} has {len(ptsByObjByFrame[f][o])} points')
                ##print(f'   len(CntDataByFrame[{i}]) = {len(CntDataByFrame[i])}')
        #print(f'   \nPtsByFrame = {PtsByFrame}')
        #print(f'   \nCntDataByFrame = {CntDataByFrame}')
        
        if False:
            print('\n   The vertices of the image grid:')
            for i in range(len(verts)):
                print(f'   {verts[i]}')
            
            print('\n   The points from the input indices:')
            for i in range(len(ptsByObjByFrame[f][o])):
                print(f'   {ptsByObjByFrame[f][o][i]}')
                
            print(f'\n   len(isInside) = {len(isInside)}')
            print(f'   sum(isInside) = {sum(isInside)}')
            
        if fracProp < 1:
            print('\n   The vertices of the image grid:')
            for i in range(len(verts)):
                print(f'   {verts[i]}')
            N = len(isInside)
            limit = 10
            if N > limit:
                print(f'\n   The first {limit} points that lie outside of extent:')
                #print([f'\n   {ptsByObjByFrame[0][0][i]}' for i in range(5)])
                for i in range(limit):
                    print(f'   {ptsByObjByFrame[0][0][i]}')
            else:
                print(f'\n   The first {N} points that lie outside of extent:')
                #print([f'\n   {ptsByObjByFrame[0][0][i]}' for i in range(N)])
                for i in range(N):
                    print(f'   {ptsByObjByFrame[0][0][i]}')
    
    return fracProp

def prop_of_segs_in_extent(
        pixarrBySeg, f2sIndsBySeg, srcIm, trgIm, p2c=False
        ):
    """
    Determine what proportion of voxels that make up a list of 3D pixel arrays 
    intersect with the volume of a 3D image.
    
    Parameters
    ----------
    pixarrBySeg : list of Numpy arrays
        A list (for each segment) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    f2sIndsBySeg : list of list of ints
        A list (for each segment) of a list (for each frame) of the slice 
        numbers that correspond to each frame in each pixel array in 
        pixarrBySeg.   
    srcIm : SimpleITK Image
        The image that relates to pixarrBySeg.
    trgIm : SimpleITK Image
        The image whose grid/extent is to be checked for intersection with the
        points of interest (e.g. trgIm if copying contours from source to
        target).
    p2c : bool, optional
        If True some results will be printed to the console. The default value
        is False.
        
    Returns
    -------
    fracProp : float
        The fractional proportion (normalised to 1) of voxels in pixarrBySeg that 
        intersect the extent of trgIm.
    """
    
    #from conversion_tools.inds_pts_pixarrs import pixarrByRoi_to_ptsByCntByRoi
    
    # Convert the list of pixel arrays-by-segment to a list of 
    # points-by-contour-by-ROI:
    ptsByCntByRoi, cntdataByCntByRoi, c2sIndsByRoi\
        = pixarrByRoi_to_ptsByCntByRoi(pixarrBySeg, f2sIndsBySeg, 
                                       srcDcmDir, Thresh=0.5)
    
    fracProp = prop_of_rois_in_extent(ptsByCntByRoi, trgIm, p2c)
    
    return fracProp