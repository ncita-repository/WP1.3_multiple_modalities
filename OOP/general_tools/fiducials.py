# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 08:51:21 2021

@author: ctorti
"""

import SimpleITK as sitk
from conversion_tools.inds_pts_cntdata import inds_to_pts
from general_tools.general import flatten_list
    
    
def import_fiducials(filepath, flipK=False, refIm=None, p2c=False):
    """
    Import a list of indices or points from a suitably formatted text file.
    
    Parameters
    ----------  
    filepath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    flipK : bool, optional (False by default)
        If the fiducials are indices derived from ImageJ, set to True to flip
        the z-index, i.e. [i, j, S - k], where S = the number of slices in the
        3D image.
    refIm : SimpleITK image, optional unless flipK is True (None by default)
        The 3D SimpleITK image that corresponds to the fiducials.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------   
    fids : list of a list of ints or floats
        A list (for each index/point) of a list (for each dimension) of either
        indices (integers) or points (floats), 
        
        e.g. 1: [[i0, j0, k0], [i1, j1, k1], ...] 
        
        e.g. 2: [[x0, y0, z0], [x1, y1, z1], ...]
    fidType : str
        Denotes the fiducial type as 'index' or 'point'.
    
    Note
    ----
    For historical reasons, the text file must have the format expected of 
    Elastix/SimpleElastix, i.e.:
        
    <index or point>
    <number of indices/points>
    point1 x point1 y [point1 z]
    point2 x point2 y [point2 z]
    
    e.g.
    
    point
    3
    102.8 -33.4 57.0
    178.1 -10.9 14.5
    180.4 -18.1 78.9
    
    https://simpleelastix.readthedocs.io/PointBasedRegistration.html
    """
    
    if flipK == True and refIm == None:
        msg = "'flipK' is True so a reference image must be provided."
        
        raise Exception(msg)
    
    if not '.txt' in filepath:
        filepath += '.txt'
    
    fids = []
    
    with open(filepath, 'r') as file:
        contents = file.readlines()
            
        fidType = contents[0].replace('\n', '')
        
        if fidType not in ['index', 'point']:
            msg = f"The first line in {filepath} must be either 'index' or "\
                  + f"'point'.  '{fidType}' is not a valid entry."
            raise Exception(msg)
            
        n = int(contents[1].replace('\n', ''))
            
        for i in range(2, len(contents)):
            items = contents[i].replace('\n', '').split(' ')
            
            if fidType == 'index':
                fiducial = [int(item) for item in items]
                
                if flipK:
                    S = refIm.GetSize()[2]
                    
                    #fiducial[2] = S - fiducial[2] # 26/04/21
                    fiducial[2] = (S - 1) - fiducial[2] # 26/04/21
                    
            else:
                fiducial = [float(item) for item in items]
            
            fids.append(fiducial)
    
    N = len(fids)
    
    if N != n:
        msg = f"Warning: The second line in {filepath} specifies that there "\
              f"are {n} fiducials but {N} fiducials were parsed."
        
        print(msg)
    
    if p2c:
        print(f'Fiducial type: {fidType}')
        print(f'There are {len(fids)} fiducials: \n {fids}\n') 
    
    return fids, fidType

def get_fiducials_within_im_extents(
        fixFids, movFids, fixIm, movIm, buffer=3, p2c=False
        ):
    """
    Return fiducials belonging to two images which lie within its image extent, 
    and do not lie on the boundary or to within buffer pixels of it.  The 
    fiducials must be indices (not points).
    
    Parameters
    ----------
    fixIm : SimpleITK image 
        The 3D image that movIm will be registered to.
    fixFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for fixIm.
    movIm : SimpleITK image
        The 3D image that will be registered to fixIm.
    movFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for movIm.
    buffer : int, optional (3 by default)
        Rather than eliminating all fiducials at or beyond the boundary, also
        filter out any fiducials that lie within buffer pixels of any boundary.
    p2c : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    fixFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for fixIm.
    movFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for movIm.
    commonInds : list of ints
        A list of indices of the original fiducials that lie within both image
        boundaries.
    
    Note
    ----
    For BSpline Transforms it's not enough for the points to lie within the    
    extent of the image - they also must not lie on the boundary:

    https://itk.org/pipermail/insight-users/2012-June/044975.html
    
    It seems that a buffer of 2 (minimum) or 3 pixels is appropriate.
    """
    
    I_f, J_f, K_f = fixIm.GetSize()
    I_m, J_m, K_m = movIm.GetSize()
    
    # Store the indices of the fiducials that do not lie on any boundary of 
    # either image:
    commonInds = []
    
    for n in range(len(fixFids)):
        i_f, j_f, k_f = fixFids[n]
        i_m, j_m, k_m = movFids[n]
        
        # The n^th fiducial lies within the boundaries of both images if all of
        # the following are True:
        b0 = i_f > 0 + buffer
        b1 = i_f < I_f - 1 - buffer
        b2 = j_f > 0 + buffer
        b3 = j_f < J_f - 1 - buffer
        b2 = k_f > 0 + buffer
        b3 = k_f < K_f - 1 - buffer
        b4 = i_m > 0 + buffer
        b5 = i_m < I_m - 1 - buffer
        b6 = j_m > 0 + buffer
        b7 = j_m < J_m - 1 - buffer
        b8 = k_m > 0 + buffer
        b9 = k_m < K_m - 1 - buffer
        
        if b0 & b1 & b2 & b3 & b4 & b5 & b6 & b7 & b8 & b9:
            commonInds.append(n)
        
    fixFids = [fixFids[i] for i in commonInds]
    movFids = [movFids[i] for i in commonInds]
    
    if p2c:
        print(f'After rejecting fiducials on or to within {buffer} pixels of',
              f'the image boundaries there are now {len(fixFids)} fiducials.\n',
              f'fixFids:\n {fixFids}\nmovFids: {movFids}\n')
    
    return fixFids, movFids, commonInds

def get_landmark_pts_from_imageJ_inds(fidFpath, image):
    """ 
    Get a flat list of points from a list of indices from ImageJ's
    coordinate system.
    
    Parameters
    ----------
    fidFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    image : SimpleITK image
        The 3D DICOM image for which the fiducials relate to.
    
    Returns
    -------
    pts : list of floats
        Flat list of [x, y, z] coordinates converted from each index the text
        file, e.g. [x0, y0, z0, x1, y1, z1, ...]).
    """
    
    R, C, S = image.GetSize()
    
    inds, fidType = import_fiducials(fidFpath)
    
    if fidType != 'index':
        msg = "The text file must contain indices, not points."
        raise Exception(msg)
        
    """ Slice indexing in ImageJ is opposite to that of sitk. """
    inds = [[ind[0], ind[1], S - ind[2]] for ind in inds]
    
    pts = [image.TransformIndexToPhysicalPoint(ind) for ind in inds]
    
    """ Flatten pts """
    pts = [c for p in pts for c in p]
    
    return pts

def get_landmark_tx(
        txName, fixFidsFpath, movFidsFpath, fixIm, movIm,
        flipK=False, buffer=3, p2c=False
        ):
    """
    Return the landmark-based SimpleITK Transform based on fiducials.
    
    Parameters
    ----------
    txName : str
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    fixFidsFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for fixIm. 
        The string need not contain the .txt extension.
    movFidsFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for movIm. 
        The string need not contain the .txt extension.
    fixIm : SimpleITK Image 
        The 3D image that movIm will be aligned to.
    movIm : SimpleITK Image
        The 3D image that will be aligned to fixIm.
    flipK : bool, optional (False by default)
        Set to True if the fiducials specified in fixFidsFpath and movFidsFpath
        are indices obtained from ImageJ.
    buffer : int, optional (3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    landmarkTx : SimpleITK Transform
        The transform used to perform landmark-based alignment.
    fixPts : list of floats
        Flat list of landmark coordinates for fixIm.
    movPts : list of floats
        Flat list of landmark coordinates for movIm.
    
    Note
    ----
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(fixIm.GetDimension())
        sitk.BSplineTransform(fixIm.GetDimension())
        
    numControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    if not txName in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {txName}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    dim = fixIm.GetDimension()
    
    if txName == 'rigid':
        sitkTx = sitk.VersorRigid3DTransform()
    elif txName == 'affine':
        sitkTx = sitk.AffineTransform(dim)
    else:
        sitkTx = sitk.BSplineTransform(dim)
    
    # Import the fiducials from text files:
    fixFids, fixFidType = import_fiducials(fixFidsFpath, p2c=p2c)

    movFids, movFidType = import_fiducials(movFidsFpath, p2c=p2c)
    
    if txName == 'bspline':
        # Reject any fiducials at or near any image boundary:
        print(f'txName = {txName} so will reject fiducials at or near',
              'the boundary.\n')
        fixFids, movFids, commonInds = get_fiducials_within_im_extents(
            fixFids=fixFids, movFids=movFids, 
            fixIm=fixIm, movIm=movIm,
            buffer=buffer, p2c=p2c
            )
    
    # Convert from lists of indices to lists of points:
    if fixFidType == 'index':
        fixPts = inds_to_pts(indices=fixFids, refIm=fixIm)
    
    if movFidType == 'index':
        movPts = inds_to_pts(indices=movFids, refIm=movIm)
    
    if p2c:
        print(f'There are {len(fixPts)} fixed points: \n {fixPts}\n')
        print(f'There are {len(movPts)} moving points: \n {movPts}\n')
    
    # Flatten the list of points (as required by 
    # LandmarkBasedTransformInitializer):
    fixPts = flatten_list(fixPts)
    movPts = flatten_list(movPts)
    
    """ 
    Additional check (for troubleshooting) - are all points within image
    extents? 
    
    from GeneralTools import IsPointInPolygon
    
    fixVerts = GetImageVertices(image=fixIm)
    movVerts = GetImageVertices(image=movIm)
    
    fixPtsInExtent = []
    movPtsInExtent = []
    
    for i in range(len(fixPts)):
        fixPtsInExtent.append(IsPointIPolygon(fixPts[], fixVerts))
    """
    
    if txName == 'bspline':
        landmarkTx = sitk.LandmarkBasedTransformInitializer(
            transform=sitkTx, fixedLandmarks=fixPts, movingLandmarks=movPts,
            referenceImage=fixIm, numberOfControlPoints=4
            #BSplineNumberOfControlPoints=numControlPts
            )
    else:
        landmarkTx = sitk.LandmarkBasedTransformInitializer(
            transform=sitkTx, fixedLandmarks=fixPts, movingLandmarks=movPts,
            referenceImage=fixIm
            )
    
    # Cast landmarkTx from sitkTransform to the desired sitk transform:
    # (01/07/21)
    if txName == 'rigid':
        landmarkTx = sitk.VersorRigid3DTransform(landmarkTx)
    elif txName == 'affine':
        landmarkTx = sitk.AffineTransform(landmarkTx)
    else:
        landmarkTx = sitk.BSplineTransform(landmarkTx)
        
    if p2c:
        #print(f'Fixed fiducials (points): \n {fixPts}\n')
        #print(f'Moving fiducials (points): \n {movPts}\n')
        #print(f'Landmark-based initial transformation is: \n {landmarkTx.GetParameters()}')
        print("-------")
        print('Landmark-based initial transform:')
        print(landmarkTx)
        print("-------")
    
    #return landmarkTx # 01/07/21
    return landmarkTx, fixPts, movPts