# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 08:51:21 2021

@author: ctorti
"""

import SimpleITK as sitk
from conversion_tools.inds_pts import inds_to_pts
from general_tools.general import flatten_list
    
    
def import_fiducials(Filepath, FlipK=False, RefIm=None, LogToConsole=False):
    """
    Import a list of indices or points from a suitably formatted text file.
    
    Parameters
    ----------  
    Filepath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    FlipK : bool, optional (False by default)
        If the fiducials are indices derived from ImageJ, set to True to flip
        the z-index, i.e. [i, j, S - k], where S = the number of slices in the
        3D image.
    RefIm : SimpleITK image, optional unless FlipK is True (None by default)
        The 3D SimpleITK image that corresponds to the fiducials.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------   
    Fiducials : list of a list of ints or floats
        A list (for each index/point) of a list (for each dimension) of either
        indices (integers) or points (floats), 
        
        e.g. 1: [[i0, j0, k0], [i1, j1, k1], ...] 
        
        e.g. 2: [[x0, y0, z0], [x1, y1, z1], ...]
    
    FiducialType : str
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
    
    if FlipK == True and RefIm == None:
        msg = "'FlipK' is True so a reference image must be provided."
        
        raise Exception(msg)
    
    if not '.txt' in Filepath:
        Filepath += '.txt'
    
    Fiducials = []
    
    with open(Filepath, 'r') as file:
        contents = file.readlines()
            
        FiducialType = contents[0].replace('\n', '')
        
        if FiducialType not in ['index', 'point']:
            msg = f"The first line in {Filepath} must be either 'index' or "\
                  + f"'point'.  '{FiducialType}' is not a valid entry."
            raise Exception(msg)
            
        n = int(contents[1].replace('\n', ''))
            
        for i in range(2, len(contents)):
            items = contents[i].replace('\n', '').split(' ')
            
            if FiducialType == 'index':
                fiducial = [int(item) for item in items]
                
                if FlipK:
                    S = RefIm.GetSize()[2]
                    
                    #fiducial[2] = S - fiducial[2] # 26/04/21
                    fiducial[2] = (S - 1) - fiducial[2] # 26/04/21
                    
            else:
                fiducial = [float(item) for item in items]
            
            Fiducials.append(fiducial)
    
    N = len(Fiducials)
    
    if N != n:
        msg = f"Warning: The second line in {Filepath} specifies that there "\
              f"are {n} fiducials but {N} fiducials were parsed."
        
        print(msg)
    
    if LogToConsole:
        print(f'Fiducial type: {FiducialType}')
        print(f'There are {len(Fiducials)} fiducials: \n {Fiducials}\n') 
    
    return Fiducials, FiducialType

def get_fiducials_within_im_extents(FixFids, MovFids, FixIm, MovIm, Buffer=3,
                                    LogToConsole=False):
    """
    Return fiducials belonging to two images which lie within its image extent, 
    and do not lie on the boundary or to within Buffer pixels of it.  The 
    fiducials must be indices (not points).
    
    Parameters
    ----------
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
    FixFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for FixIm.
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    MovFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for MovIm.
    Buffer : int, optional (3 by default)
        Rather than eliminating all fiducials at or beyond the boundary, also
        filter out any fiducials that lie within Buffer pixels of any boundary.
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    FixFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for FixIm.
    MovFids : list of a list of ints
        A list (for each fiducial) of a list of [i, j, k] indices for MovIm.
    CommonInds : list of ints
        A list of indices of the original fiducials that lie within both image
        boundaries.
    
    Note
    ----
    For BSpline Transforms it's not enough for the points to lie within the    
    extent of the image - they also must not lie on the boundary:

    https://itk.org/pipermail/insight-users/2012-June/044975.html
    
    It seems that a buffer of 2 (minimum) or 3 pixels is appropriate.
    """
    
    I_f, J_f, K_f = FixIm.GetSize()
    I_m, J_m, K_m = MovIm.GetSize()
    
    # Store the indices of the fiducials that do not lie on any boundary of 
    # either image:
    CommonInds = []
    
    for n in range(len(FixFids)):
        i_f, j_f, k_f = FixFids[n]
        i_m, j_m, k_m = MovFids[n]
        
        # The n^th fiducial lies within the boundaries of both images if all of
        # the following are True:
        b0 = i_f > 0 + Buffer
        b1 = i_f < I_f - 1 - Buffer
        b2 = j_f > 0 + Buffer
        b3 = j_f < J_f - 1 - Buffer
        b2 = k_f > 0 + Buffer
        b3 = k_f < K_f - 1 - Buffer
        b4 = i_m > 0 + Buffer
        b5 = i_m < I_m - 1 - Buffer
        b6 = j_m > 0 + Buffer
        b7 = j_m < J_m - 1 - Buffer
        b8 = k_m > 0 + Buffer
        b9 = k_m < K_m - 1 - Buffer
        
        if b0 & b1 & b2 & b3 & b4 & b5 & b6 & b7 & b8 & b9:
            CommonInds.append(n)
        
    FixFids = [FixFids[i] for i in CommonInds]
    MovFids = [MovFids[i] for i in CommonInds]
    
    if LogToConsole:
        print(f'After rejecting fiducials on or to within {Buffer} pixels of',
              f'the image boundaries there are now {len(FixFids)} fiducials.\n',
              f'FixFids:\n {FixFids}\nMovFids: {MovFids}\n')
    
    return FixFids, MovFids, CommonInds

def get_landmark_pts_from_ImageJ_inds(FiducialFpath, Image):
    """ 
    Get a flat list of points from a list of indices from ImageJ's
    coordinate system.
    
    Parameters
    ----------
    FiducialFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    Image : SimpleITK image
        The 3D DICOM image for which the fiducials relate to.
    
    Returns
    -------
    Pts : list of floats
        Flat list of [x, y, z] coordinates converted from each index the text
        file, e.g. [x0, y0, z0, x1, y1, z1, ...]).
    """
    
    R, C, S = Image.GetSize()
    
    Inds, FiducialType = import_fiducials(FiducialFpath)
    
    if FiducialType != 'index':
        msg = "The text file must contain indices, not points."
        raise Exception(msg)
        
    """ Slice indexing in ImageJ is opposite to that of sitk. """
    Inds = [[ind[0], ind[1], S - ind[2]] for ind in Inds]
    
    Pts = [Image.TransformIndexToPhysicalPoint(ind) for ind in Inds]
    
    """ Flatten Pts """
    Pts = [c for p in Pts for c in p]
    
    return Pts

def get_landmark_tx(Transform, FixFidsFpath, MovFidsFpath, FixIm, MovIm,
                    FlipK=False, Buffer=3, LogToConsole=False):
    """
    Return the landmark-based SimpleITK Transform based on fiducials.
    
    Parameters
    ----------
    Transform : str
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    FixFidsFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension.
    MovFidsFpath : str
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension.
    FixIm : SimpleITK Image 
        The 3D image that MovIm will be aligned to.
    MovIm : SimpleITK Image
        The 3D image that will be aligned to FixIm.
    FlipK : bool, optional (False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ.
    Buffer : int, optional (3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    LandmarkTx : SimpleITK Transform
        The transform used to perform landmark-based alignment.
    FixPts : list of floats
        Flat list of landmark coordinates for FixIm.
    MovPts : list of floats
        Flat list of landmark coordinates for MovIm.
    
    Note
    ----
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(FixIm.GetDimension())
        sitk.BSplineTransform(FixIm.GetDimension())
        
    NumControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {Transform}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    dim = FixIm.GetDimension()
    
    if Transform == 'rigid':
        SitkTx = sitk.VersorRigid3DTransform()
    elif Transform == 'affine':
        SitkTx = sitk.AffineTransform(dim)
    else:
        SitkTx = sitk.BSplineTransform(dim)
    
    # Import the fiducials from text files:
    FixFids, FixFidType = import_fiducials(FixFidsFpath, 
                                           LogToConsole=LogToConsole)

    MovFids, MovFidType = import_fiducials(MovFidsFpath,
                                           LogToConsole=LogToConsole)
    
    if Transform == 'bspline':
        # Reject any fiducials at or near any image boundary:
        print(f'Transform = {Transform} so will reject fiducials at or near',
              'the boundary.\n')
        FixFids, MovFids, CommonInds\
            = get_fiducials_within_im_extents(FixFids=FixFids, MovFids=MovFids, 
                                              FixIm=FixIm, MovIm=MovIm,
                                              Buffer=Buffer,
                                              LogToConsole=LogToConsole)
    
    # Convert from lists of indices to lists of points:
    if FixFidType == 'index':
        FixPts = inds_to_pts(Indices=FixFids, RefIm=FixIm)
    
    if MovFidType == 'index':
        MovPts = inds_to_pts(Indices=MovFids, RefIm=MovIm)
    
    if LogToConsole:
        print(f'There are {len(FixPts)} fixed points: \n {FixPts}\n')
        print(f'There are {len(MovPts)} moving points: \n {MovPts}\n')
    
    # Flatten the list of points (as required by 
    # LandmarkBasedTransformInitializer):
    FixPts = flatten_list(FixPts)
    MovPts = flatten_list(MovPts)
    
    #""" Additional check (for troubleshooting) - are all points within image
    #extents? """
    #from GeneralTools import IsPointInPolygon
    #
    #FixVerts = GetImageVertices(Image=FixIm)
    #MovVerts = GetImageVertices(Image=MovIm)
    #
    #FixPtsInExtent = []
    #MovPtsInExtent = []
    #
    #for i in range(len(FixPts)):
    #    FixPtsInExtent.append(IsPointIPolygon(FixPts[], FixVerts))
        
    if Transform == 'bspline':
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=SitkTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm,
                                                            numberOfControlPoints=4)
                                                            #BSplineNumberOfControlPoints=NumControlPts)
    else:
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=SitkTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm)
    
    # Cast LandmarkTx from sitkTransform to the desired sitk transform:
    # (01/07/21)
    if Transform == 'rigid':
        LandmarkTx = sitk.VersorRigid3DTransform(LandmarkTx)
    elif Transform == 'affine':
        LandmarkTx = sitk.AffineTransform(LandmarkTx)
    else:
        LandmarkTx = sitk.BSplineTransform(LandmarkTx)
        
    if LogToConsole:
        #print(f'Fixed fiducials (points): \n {FixPts}\n')
        #print(f'Moving fiducials (points): \n {MovPts}\n')
        #print(f'Landmark-based initial transformation is: \n {LandmarkTx.GetParameters()}')
        print("-------")
        print('Landmark-based initial transform:')
        print(LandmarkTx)
        print("-------")
    
    #return LandmarkTx # 01/07/21
    return LandmarkTx, FixPts, MovPts