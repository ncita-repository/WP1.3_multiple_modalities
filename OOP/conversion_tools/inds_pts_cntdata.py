# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:47:07 2021

@author: ctorti
"""

""" Low-level conversions between indices and points and visa versa. """


import numpy
from general_tools.general import get_list_of_dtypes


def ind_to_pt(index, refIm):
    """
    Convert an index (in the Image Coordinate System) to a physical point (in 
    the Patient Coordinate System).
    
    Parameters
    ----------
    index : list of ints or floats
        A list (for each dimension) of an index, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that index belongs to.
        
    Returns
    -------
    point : list of floats
        A list (for each dimension) of the coordinates of index, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    # Check the data type of the first element in index:
    if isinstance(index[0], int):
        point = refIm.TransformIndexToPhysicalPoint(index)
    elif isinstance(index[0], float):
        point = refIm.TransformContinuousIndexToPhysicalPoint(index)
    else:
        msg = f"The data type of index is {type(index[0])}. It must be either"\
              + " 'int' or 'float'."
        raise Exception(msg)
         
    return list(point)

def ind_to_pt_210921(index, refIm):
    """
    14/09/21: This is unnecessarily complex.
    
    Convert an index (in the Image Coordinate System) to a physical point (in 
    the Patient Coordinate System).
    
    Parameters
    ----------
    index : list of ints or floats
        A list (for each dimension) of an index, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that index belongs to.
        
    Returns
    -------
    point : list of floats
        A list (for each dimension) of the coordinates of index, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    #import numpy
    #from general_tools.general import get_list_of_dtypes
    
    dtypes = get_list_of_dtypes(index)
    
    if len(dtypes) == 1:
        dtype = dtypes[0]
        
        """All items in index are of a common data type."""
        if dtype == int:
            point = refIm.TransformIndexToPhysicalPoint(index)
        
        elif dtype == float:
            point = refIm.TransformContinuousIndexToPhysicalPoint(index)
        
        elif dtype in [numpy.ndarray, numpy.float64, numpy.int32]:
            """ Convert from numpy array to list: """
            index = index.tolist()
            
            dtypes = get_list_of_dtypes(index)
            
            if len(dtypes) == 1:
                dtype = dtypes[0]
                
                """All items in index are of a common data type."""
                if dtype == int:
                    point = refIm.TransformIndexToPhysicalPoint(index)
                
                elif dtype == float:
                    point = refIm.TransformContinuousIndexToPhysicalPoint(index)
                
                else:
                    msg = f"The data type of index is {dtype}. It must be "\
                          + "'int', 'float', 'numpy.ndarray' or 'numpy.float64'."
            
                    raise Exception(msg)
                    
            elif len(dtypes) == 2 and int in dtypes and float in dtypes:
                """Allow for the possibility of a mixture of integers and 
                floats."""
                
                """ Convert integers to float: """
                index = [float(item) for item in index]
                
                point = refIm.TransformContinuousIndexToPhysicalPoint(index)
                
            else:
                msg = f"The data type of index is {dtypes}. It must be 'int',"\
                      + " 'float' or 'numpy.ndarray'."
            
        else:
            msg = f"The data type of index is {dtypes}. It must be 'int',"\
                  + " 'float', 'numpy.ndarray' or 'numpy.float64'."
            raise Exception(msg)
            
    elif len(dtypes) == 2 and int in dtypes and float in dtypes:
        """Allow for the possibility of a mixture of integers and floats."""
        
        """ Convert integers to float: """
        index = [float(item) for item in index]
        
        point = refIm.TransformContinuousIndexToPhysicalPoint(index)
        
    else:
        msg = f"The data type of index is {dtypes}. It must be either "\
              + "'int' or 'float'."
        raise Exception(msg)
         
    #return point # 29/01/2021
    return list(point) # 29/01/2021

def inds_to_pts(indices, refIm):
    """
    Convert a list of indices (in the Image Coordinate System) to a list of 
    physical points (in the Patient Coordinate System).
    
    Parameters
    ----------
    indices : list of a list of ints or floats
        A list (for each index) of a list (for each dimension) of indices, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that the indices belong to.
        
    Returns
    -------
    points : list of a list of floats
        A list (for each point) of a list (for each dimension) of the 
        coordinates of indices, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    points = []
    
    for index in indices:
        #print(f'\n\nindex = {index}')
        #index = [int(item) for item in index]
        
        point = ind_to_pt(index, refIm)
        
        points.append(point)
        
    return points

def pt_to_ind(point, refIm, rounding=True):
    """
    Convert a physical point (in the Patient Coordinate System) to an index
    (in the Image Coordinate System).
    
    Parameters
    ----------
    point : list of floats
        A list (for each dimension) of coordinates, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that point belongs to.
    rounding : bool, optional (True by default)
        If True, the integers in index will be rounded to the nearest integers.
        
    Returns
    -------
    index : list of ints or floats
        List (for each dimension) of the indices of points, e.g. 
            [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    if rounding:
        index = refIm.TransformPhysicalPointToIndex(point)
    else:
        index = refIm.TransformPhysicalPointToContinuousIndex(point)
    
    #return index # 29/01/2021
    return list(index) # 29/01/2021

def pts_to_inds(points, refIm, rounding=True):
    """
    Convert a list of physical points (in the Patient Coordinate System) to a
    list of indices (in the Image Coordinate System).
    
    Parameters
    ----------
    points : list of a list of floats
        A list (for each point) of a list (for each dimension) of coordinates, 
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    refIm : SimpleITK image
        A 3D image that occupies the grid that points belong to.
    rounding : bool optional (True by default)
        If True, the integers in index will be rounded to the nearest integers.
    
    Returns
    -------
    indices : list of a list of ints
        List (for each point) of a list (for each dimension) of the indices 
        of points, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    indices = []
    
    #print(f'\n\n\npoints = {points}\n\n')
    
    for point in points:
        #print(f'point = {point}')
        
        index = pt_to_ind(point, refIm, rounding)
        
        #print(f'index = {index}')
        
        indices.append(index)
        
    return indices

def cntdata_to_pts(cntdata):
    """
    Re-format a flat list of coordinates to a list of a list of [x, y, z]
    coordinates.
    
    Parameters
    ----------
    cntdata : list of strs
        Flat list of [x, y, z] coordinates as strings.
    
    Returns
    -------
    points : list of a list of floats
        List (for each point) of a list (for each dimension) of points,
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    cntdata = [float(item) for item in cntdata]
    
    points = []
            
    # Iterate for all points in threes:
    for p in range(0, len(cntdata), 3):
        point = [cntdata[p], cntdata[p+1], cntdata[p+2]]
        
        points.append(point)
        
    return points

def pts_to_cntdata(points):
    """
    Re-format a list of points to a flat list as required for the DICOM tag
    ContourData.
    
    Parameters
    ----------
    points : list of a list of floats
        List (for each point) of a list (for each dimension) of points,
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
    Returns
    -------
    cntdata : list of strs
        Flat list of coordinates in points converted from floats to strings.
    """
    
    cntdata = []
    
    for point in points:
        point = [str(item) for item in point]
        
        cntdata.extend(point)
        
    return cntdata

def ptsByCntByRoi_to_cntdataByCntByRoi(ptsByCntByRoi):
    """
    Convert a list of points-by-contour-by-ROI to a list of 
    cntdata-by-contour-by-ROI, where cntdata is flat list of points in
    string format (as required for ContourData DICOM tag).
    
    Parameters
    ----------
    ptsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
    Returns
    -------
    cntdataByCntByRoi : list of a list of a list of strs
        List (for each ROI) of a list (for all contours) of a flat list of 
        coordinates in ptsByCntByRoi converted from floats to strings.
    """
    
    cntdataByCntByRoi = []
        
    for r in range(len(ptsByCntByRoi)):
        cntdataByCnt = []
        
        if ptsByCntByRoi[r]:
            for c in range(len(ptsByCntByRoi[r])):
                pts = ptsByCntByRoi[r][c]
                
                cntdataByCnt.append(pts_to_cntdata(pts))
        
        cntdataByCntByRoi.append(cntdataByCnt)
        
    return cntdataByCntByRoi