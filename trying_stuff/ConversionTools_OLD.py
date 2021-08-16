# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:29:58 2020

@author: ctorti
"""





def Inds2SparseLabmap_OLD(Indices, RefImage): # not useful since Sparse2Labelmap doesn't work
    """
    Note:
        
        The problem with this function is that by using the SimpleITK image 
        format any slice that contains more than one contour (e.g. belonging to 
        more than ROI) will exist on the same slice.  What I want is a Numpy 
        array with a different mask for each contour.  
        
        The other reason for using the SimpleITK
        image format was because I had intended to use sitk filters to convert 
        the sparse mask into a filled one (e.g. using 
        BinaryClosingByReconstructionImageFilter), but since that failed I have
        no reason to use the sitk image format.
        
    
    Convert a list of indices to a sparse 3D Numpy data array.  
    
    Inputs:
        Indices  - (List of a list of integers) List of list of point indices
        
        RefImage - (SimpleITK image) Reference image whose size will be used to
                   define the shape of the sparse data array
        
        
    Returns:
        #Matrix   - (Numpy array, dtype int) The point indices as a Numpy array
        Sparse   - (SimpleITK image) The point indices as a SimpleITK image
    """
    
    import numpy as np
    import SimpleITK as sitk
    
    ImSize = RefImage.GetSize()
    
    # Initialise the sparse matrix:
    """
    Note: .GetSize() provides the size as (Ncols, Nrows, Nslices) but for
    Numpy arrays the order will be (Nslices, Nrows, Ncols).
    """
    Matrix = np.zeros((ImSize[2], ImSize[1], ImSize[0]), dtype=int)
    
    #print(f'Matrix.shape = {Matrix.shape}')
    
    # Loop through all contours:
    for contour in Indices:
        # Ignore any contours that have less than 3 (the minimum number to
        # define a closed contour) indices (as may arise during manual 
        # contouring):
        if len(contour) > 2:
            # Each i,j,k index along x,y,z:
            for i,j,k in contour:
                #print(f'i = {i}, j = {j}, k = {k}')
                
                Matrix[k,j,i] = 1
                
    
    #return Matrix
    return sitk.GetImageFromArray(Matrix)
            
    






def Sparse2Labelmap_OLD(Sparse): # doesn't work
    """
    This doesn't work.
    
    
    Convert a sparse labelmap to a filled labelmap, i.e. close the contour
    indices to result in masks.  
    
    Inputs:
        #Matrix   - (Numpy array, dtype int) The point indices as a Numpy array
        Sparse   - (SimpleITK image) Sparse labelmap (containing only indices
                   that define the contours)
        
        
    Returns:
        Labmap   - (SimpleITK image) Labelmap from closing of contour indices
    """
    
    import SimpleITK as sitk
    
    # First convert from "labelmap" to binary:
    #Labmap2BinaryFilt = sitk.LabelMapToBinaryImageFilter()
    
    #Labmap = Labmap2BinaryFilt.Execute(Sparse)
    
    # Close labelmap:
    ClosingFilt = sitk.BinaryClosingByReconstructionImageFilter()
    ClosingFilt.FullyConnectedOn()
    
    
    #ClosingFilt = sitk.BinaryMorphologicalClosingImageFilter()
    
    
    
    #Labmap = ClosingFilt.Execute(Matrix)
    #Labmap = ClosingFilt.Execute(Labmap)
    Labmap = ClosingFilt.Execute(Sparse)
    
    
    return Labmap







def Contour2Indices_OLD(Points, Origin, Directions, Spacings, Rounding=True):
    """
    Convert a list of physical point (in Patient Coordinate System (PCS)) to a
    list of indeces (in Image Coordinate System (ICS).
    
    Inputs:
    ------
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points in the
        PCS, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
    Origin : list of floats
        A list (for each dimension) of the 3D image origin 
        (= ImagePositionPatient of the first slice in the DICOM series), 
        e.g. [x0, y0, z0].
        
    Directions : list of floats 
        List of the direction cosines along x (rows), y (columns) and z 
        (slices) (= the cross product of the x and y direction cosines from 
        ImageOrientationPatient appended to ImageOrientationPatient),
        e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz].
        
    Spacings : list of floats 
        List (for each dimension) of the pixel spacings along x, y and z, 
        e.g. [di, dj, dk].  The x and y spacings come from PixelSpacing.  The 
        z spacing is not the SliceThickness but instead is determined from 
        ImagePositionPatient and ImageOrientationPatient.
    
    Rounding : boolean (optional; True by default)
        If True, the integers in Index will be rounded to the nearest integers.
    
        
    Ouputs:
    ------
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        (in ICS) of Points, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    Indices = []
    
    print(f'\n\n\nPoints = {Points}\n\n')
    
    for Point in Points:
        print(f'Point = {Point}')
        
        Indices.append(Point2Index(Point, Origin, Directions, Spacings, 
                                   Rounding))
        
    return Indices






def ConvertContourDataToIndices_OLD(ContourData, DicomDir):
    """
    Not useful?..
    
    Convert contour data from a flat list of points in the Patient Coordinate
    System (PCS) to a list of [i, j, k] indices in the Image Coordinate System. 
    
    Inputs:
        ContourData - (List of floats) Flat list of points in the PCS,
                      e.g. [x0, y0, z0, x1, y1, z1, ...]
        
        DicomDir    - (String) Directory containing DICOMs
        
        
    Returns:
        Indices   - (List of floats) List of a list of [i, j, k] coordinates
                      of each point in ContourData converted to the ICS,
                      e.g. [[i0, j0, k0], [i1, j1, k1], ...]
    """
    
    from ImageTools import GetImageAttributes
    
    # Get the image attributes:
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    # Initialise the array of contour points:
    Indices = []
    
    # Iterate for all contour points in threes:
    for p in range(0, len(ContourData), 3):
        # The point in Patient Coordinate System:
        pt = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        # Convert pt to Image Coordinate System:
        ind = Point2Index(pt, IPPs[0], Dirs, Spacings)
        
        Indices.append(ind)
        
    return Indices







def Point2Index_OLD(Point, Origin, Directions, Spacings, Rounding=True):
    """
    WARNING:  Due to simplifications (see Equations below), this function is 
    only valid for the case where the cosine directions are:
        X = [1, 0, 0]
        Y = [0, 1, 0]
        Z = [0, 0, 1]
        
    COMMENT:
        This function was unnecessary since there's a built in function in 
        SimpleITK that does this:
        Index = RefIm.TransformPhysicalPointToIndex(Point)
        
    
                      
    Convert a physical point (in Patient Coordinate System (PCS)) to an index
    (in Image Coordinate System (ICS).
    
    Inputs:
    ------
    
    Point : list of floats
        A list (for each dimension) of a point in PCS, e.g. [x, y, z].
        
    Origin : list of floats
        A list (for each dimension) of the 3D image origin 
        (= ImagePositionPatient of the first slice in the DICOM series), 
        e.g. [x0, y0, z0].
        
    Directions : list of floats 
        List of the direction cosines along x (rows), y (columns) and z 
        (slices) (= the cross product of the x and y direction cosines from 
        ImageOrientationPatient appended to ImageOrientationPatient),
        e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz].
        
    Spacings : list of floats 
        List (for each dimension) of the pixel spacings along x, y and z, 
        e.g. [di, dj, dk].  The x and y spacings come from PixelSpacing.  The 
        z spacing is not the SliceThickness but instead is determined from 
        ImagePositionPatient and ImageOrientationPatient.
        
    Rounding : boolean (optional; True by default)
        If True, the integers in Index will be rounded to the nearest integers.
    
        
    Outputs:
    -------
    
    Index : list of integers
        List (for each dimension) of the index (i.e. in ICS) of Point, 
        e.g. [i, j, k].
        
    
    Equations:
    ---------
    
    P = S + (di*i).X + (dj*j).Y + (dk*k).Z
    
    where P = (Px, Py, Pz) is the point in the Patient Coordinate System
    
          S = (Sx, Sy, Sz) is the origin in the Patient Coordinate System (i.e.
          the ImagePositionPatient)
          
          X = (Xx, Xy, Xz) is the direction cosine vector along rows (x) 
          
          Y = (Yx, Yy, Yz) is the direction cosine vector along columns (y)
          
          Z = (Zx, Zy, Zz) is the direction cosine vector along slices (z)
          
          di, dj and dk are the pixel spacings along x, y and z
          
          i, j, and k are the row (x), column (y) and slice (z) indices
          
          
    Solve for i, j and k, i.e.:
        
        Vx = Xx*di*i + Yx*dj*j + Zx*dk*k
        
        Vy = Xy*di*i + Yy*dj*j + Zy*dk*k
        
        Vz = Xz*di*i + Yz*dk*k + Zz*dk*k
        
    where Vx = Px - Sx
    
          Vy = Py - Sy
          
          Vz = Pz - Sz
          
    Solving for i, j and k results in the expressions:
        
        i = (1/di)*d/e
            
        j = (1/dj)*( (a/b)*di*i + c/b )
        
        k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
        or
        k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
        
        (either expression for k should be fine)
        
    In the case where X = [1, 0, 0]
                      Y = [0, 1, 0]
                      Z = [0, 0, 1]
                      
    The expressions simplify to:
        
        a = 0
        
        b = 1
        
        c = Vy
        
        d = Vx
        
        e = 1
        
        i = Vx/di
        
        j = Vy/dj
        
        k = Vz/dk
    """
    
    #print(f'\nPoint = {Point}')
    #print(f'Origin = {Origin}')
    #print(f'Directions = {Directions}')
    #print(f'Spacings = {Spacings}')
    #print(f'Rounding = {Rounding}')
    
    # Define S, X, Y and Z:
    S = Origin # the origin
    
    X = Directions[0:3] # the row (x) direction cosine vector
    Y = Directions[3:6] # the column (y) direction cosine vector
    Z = Directions[6:] # the slice (z) direction cosine vector
    
    # The indices of the largest direction cosines along rows, columns and 
    # slices:
    ind_i = X.index(max(X)) # typically = 0
    ind_j = Y.index(max(Y)) # typically = 1
    ind_k = Z.index(max(Z)) # typically = 2

    # The pixel spacings:
    di = Spacings[0]
    dj = Spacings[1]
    dk = Spacings[2] 
    
    # Simplifying expressions:
    Xx = X[ind_i]
    Xy = X[ind_j]
    Xz = X[ind_k]
    
    Yx = Y[ind_i]
    Yy = Y[ind_j]
    Yz = Y[ind_k]
    
    Zx = Z[ind_i]
    Zy = Z[ind_j]
    Zz = Z[ind_k]
    
    # Define simplifying expressions:
    Vx = Point[0] - S[0]
    Vy = Point[1] - S[1]
    Vz = Point[2] - S[2]
    
    a = Xz*Zy/Zz - Xy
    
    b = Yy - Yz*Zy/Zz
    
    c = Vy - Vz*Zy/Zz
    
    d = Vx - (c/b)*Yx - (Zx/Zz)*(Vz - (c/b)*Yz)
    
    e = Xx + (a/b)*Yx - (Zx/Zz)*(Xz + (a/b)*Yz)
    
    # Solve for i, j and k:
    i = (1/di)*d/e
    
    j = (1/dj)*( (a/b)*di*i + c/b )
    
    k = (1/dk)*(1/Zz)*(Vz - Xz*di*i - Yz*dj*j)
    #k = (1/dk)*(1/Zz)*( Vz - (c/b)*Yz - (Xz + (a/b)*Yz)*di*i )
    
    #print(f'i = {i}, j = {j}, k = {k}')
    
    if Rounding: # 8/12/20
        # Convert fractional indices to integers:
        i = int(round(i))
        j = int(round(j))
        k = int(round(k))
    
    return [i, j, k]







def Index2Point_OLD(Index, Origin, Directions, Spacings):
    """
    COMMENT:
        This function was unnecessary since there's a built in function in 
        SimpleITK that does this:
        Point = RefIm.TransformIndexToPhysicalPoint(Index)
    
    Convert an index (in Image Coordinate System (ICS)) to a physical point 
    (in Patient Coordinate System (PCS)).
    
    Inputs:
    ------
    
    Index : list of integers
        List (for each dimension) of the indices (in ICS) of a point, 
        e.g. [i, j, k].
        
    Origin : list of floats
        A list (for each dimension) of the 3D image origin 
        (= ImagePositionPatient of the first slice in the DICOM series), 
        e.g. [x0, y0, z0].
        
    Directions : list of floats 
        List of the direction cosines along x (rows), y (columns) and z 
        (slices) (= the cross product of the x and y direction cosines from 
        ImageOrientationPatient appended to ImageOrientationPatient),
        e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz].
        
    Spacings : list of floats 
        List (for each dimension) of the pixel spacings along x, y and z, 
        e.g. [di, dj, dk].  The x and y spacings come from PixelSpacing.  The 
        z spacing is not the SliceThickness but instead is determined from 
        ImagePositionPatient and ImageOrientationPatient.
    
        
    Ouputs:
    ------
    
    Point : list of floats
        List (for each dimension) of the coordinates of a point, 
        e.g. [x, y, z].
        
        
    Equations:
    ---------
    
    P = S + (di*i).X + (dj*j).Y + (dk*k).Z
    
    where P = (Px, Py, Pz) is the point in the Patient Coordinate System
    
          S = (Sx, Sy, Sz) is the origin in the Patient Coordinate System (i.e.
          the ImagePositionPatient)
          
          X = (Xx, Xy, Xz) is the direction cosine vector along rows (x) 
          
          Y = (Yx, Yy, Yz) is the direction cosine vector along columns (y)
          
          Z = (Zx, Zy, Zz) is the direction cosine vector along slices (z)
          
          di, dj and dk are the pixel spacings along x, y and z
          
          i, j, and k are the row (x), column (y) and slice (z) indices
    """
    
    if False: # 08/12/20
        x = Origin[0] + Index[0]*Spacings[0]*Directions[0] \
                      + Index[0]*Spacings[0]*Directions[1] \
                      + Index[0]*Spacings[0]*Directions[2] \
        
        y = Origin[1] + Index[1]*Spacings[1]*Directions[3] \
                      + Index[1]*Spacings[1]*Directions[4] \
                      + Index[1]*Spacings[1]*Directions[5] \
            
        z = Origin[2] + Index[2]*Spacings[2]*Directions[6] \
                      + Index[2]*Spacings[2]*Directions[7] \
                      + Index[2]*Spacings[2]*Directions[8] \
                  
    x = Origin[0] + Index[0]*Spacings[0]*Directions[0] \
                  + Index[1]*Spacings[1]*Directions[3] \
                  + Index[2]*Spacings[2]*Directions[6] \
    
    y = Origin[1] + Index[0]*Spacings[0]*Directions[1] \
                  + Index[1]*Spacings[1]*Directions[4] \
                  + Index[2]*Spacings[2]*Directions[7] \
        
    z = Origin[2] + Index[0]*Spacings[0]*Directions[2] \
                  + Index[1]*Spacings[1]*Directions[5] \
                  + Index[2]*Spacings[2]*Directions[8] \
        
    return [x, y, z]






def Indices2Points_OLD(Indices, Origin, Directions, Spacings):
    """
    Convert a list of indices (in Image Coordinate System (ICS)) to a list of
    physical points (in Patient Coordinate System (PCS)).
    
    Inputs:
    ------
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        (in ICS) of Points, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
        
    Origin : list of floats
        A list (for each dimension) of the 3D image origin 
        (= ImagePositionPatient of the first slice in the DICOM series), 
        e.g. [x0, y0, z0].
        
    Directions : list of floats 
        List of the direction cosines along x (rows), y (columns) and z 
        (slices) (= the cross product of the x and y direction cosines from 
        ImageOrientationPatient appended to ImageOrientationPatient),
        e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz].
        
    Spacings : list of floats 
        List (for each dimension) of the pixel spacings along x, y and z, 
        e.g. [di, dj, dk].  The x and y spacings come from PixelSpacing.  The 
        z spacing is not the SliceThickness but instead is determined from 
        ImagePositionPatient and ImageOrientationPatient.
    
        
    Ouputs:
    ------
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points in the
        PCS, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    Points = []
    
    for index in Indices:
        Points.append(Index2Point(index, Origin, Directions, Spacings))
        
    return Points






def PixArr2IndsByFrame_OLD(PixArr, FrameToSliceInds, Thresh=0.5):
    """
    Convert a 3D pixel array to a list (for each frame) of a list of indices
    that define the equivalent contour for the mask in each frame.
    
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
 
    Inputs:
    ------
    
    PixArr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        PixArr.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
 
    
    Outputs:
    -------
    
    IndicesByFrame : list of floats
        A list (for each contour/frame) of a flat list of [x, y, z] coordinates
        of the polygons that define each frame in PixArr.
        
    PtsByCnt : list of a list of floats
        A list (for each contour/frame) of a list (for each dimension) of the 
        polygons that define each frame in PixArr.
        
    #CntToSliceInds : list of integers
    #    List of slice numbers that correspond to each contour.
    """
    
    import numpy as np
    from scipy.sparse import issparse
    from skimage.measure import find_contours
    
    IndsByFrame = []
    
    for FrameNum, Mask in enumerate(PixArr):
        if issparse(Mask):
            Mask = np.array(Mask.astype('byte').todense())
 
        if (Mask != 0).sum() == 0:
            # If Mask is empty, just skip it
            continue
 
        # Add an empty row and column at both sides of the mask to ensure that
        # any masks near an edge are found:
        ExpandedDims = (Mask.shape[0] + 2, Mask.shape[1] + 2)
        
        ExpandedMask = np.zeros(ExpandedDims, dtype=float)
        
        ExpandedMask[1:Mask.shape[0] + 1, 1:Mask.shape[1] + 1] = Mask
 
        # Get the indices that define contours in ExpandedMask:
        """
        I'm not sure how this will behave if there is more than one label (i.e.
        more than one closed mask on any given frame).  Since the returned
        item is a list of length 1, presummably a frame with multiple masks 
        will result in a list of length > 1 for each mask.  This needs to be
        tested and confirmed.
        """
        IndsByObj_2D = find_contours(ExpandedMask.T, Thresh)
        
        print(f'\n\n\nFrameNum = {FrameNum}')
        #print(f'\nlen(IndsByObj_2D) = {len(IndsByObj_2D)}')
        #print(f'\nIndsByObj_2D = {IndsByObj_2D}')
 
        # Remove one row and column to shift the indices back to their original
        # locations:
        IndsByObj_2D = [np.subtract(x, 1).tolist() for x in IndsByObj_2D]
        
        if len(IndsByObj_2D) > 1:
            print(f'\nWarning: {len(IndsByObj_2D)} objects were found in',
                  f'frame {FrameNum} in PixArr.')
            
            print(f'\nIndsByObj_2D = {IndsByObj_2D}')
            
            """
            Consider what else needs to be done. 
            e.g.
            
            for IndsThisObj_2D in IndsByObj_2D:
                ...
            """
            
        else:
            """
            Only one object was found. Reduce the level of nested lists to one
            list of 2D indices.
            """
            Inds_2D = IndsByObj_2D[0]
        
        #print(f'\n\n\nFrameNum = {FrameNum}')
        print(f'\nlen(Inds_2D) = {len(Inds_2D)}')
        print(f'\nInds_2D = {Inds_2D}')
        
        SliceNum = FrameToSliceInds[FrameNum]
            
        # Add the index for the 3rd dimension:
        Inds_3D = [Ind_2D + [SliceNum] for Ind_2D in Inds_2D]
        
        print(f'\nlen(Inds_3D) = {len(Inds_3D)}')
        print(f'\nInds_3D = {Inds_3D}')
        
        IndsByFrame.append(Inds_3D)
        
        
    print(f'\n\n\n\nlen(IndsByFrame) = {len(IndsByFrame)}')
    print(f'\nIndsByFrame = {IndsByFrame}')
 
    return IndsByFrame







def PixArr2PtsByContour_OLD(PixArr, FrameToSliceInds, DicomDir, Thresh=0.5):
    """
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates.
 
    Inputs:
    ------
    
    PixArr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        PixArr.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
 
    
    Outputs:
    -------
    
    PtsByCnt : list of a list of floats
        A list (for each contour/frame) of a list (for each dimension) of the 
        physical coordinates that define the mask in each frame in PixArr.
        
    CntDataByCnt : list of a list of strings
        A list (for each contour/frame) of a flattened list of [x, y, z]
        physical coordinates that define the mask in each frame in PixArr as
        strings (format of ContourData tag in RTS).
    """
    
    from ImageTools import GetImageAttributes
    
    # Convert the pixel array to a list of indices-by-frame:
    IndsByFrame = PixArr2IndsByFrame(PixArr, FrameToSliceInds, Thresh=0.5)
    
    Size, Spacings, ST,\
    IPPs, Dirs = GetImageAttributes(DicomDir, Package='pydicom')
    
    PtsByCnt = []
    CntDataByCnt = []
    
    for Inds in IndsByFrame:
        Pts = Indices2Points(Indices=Inds, Origin=IPPs[0], Directions=Dirs, 
                             Spacings=Spacings)
        
        PtsByCnt.append(Pts)
        
        CntDataByCnt.append(Points2ContourData(Points=Pts))
    
    
    #print(f'\n\n\nlen(PtsByCnt) = {len(PtsByCnt)}')
    #print(f'\nlen(CntDataByCnt) = {len(CntDataByCnt)}')
    #print(f'\nPtsByCnt = {PtsByCnt}')
    #print(f'\nCntDataByCnt = {CntDataByCnt}')
 
    return PtsByCnt, CntDataByCnt









def PixArr2PtsByContour_OLD2(PixArr, FrameToSliceInds, DicomDir, Thresh=0.5):
    """
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates.
 
    Inputs:
    ------
    
    PixArr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        PixArr.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
 
    
    Outputs:
    -------
    
    PtsByCnt : list of a list of floats
        A list (for each contour/frame) of a list (for each dimension) of the 
        physical coordinates that define the mask in each frame in PixArr.
        
    CntDataByCnt : list of a list of strings
        A list (for each contour/frame) of a flattened list of [x, y, z]
        physical coordinates that define the mask in each frame in PixArr as
        strings (format of ContourData tag in RTS).
    """
    
    #from ImageTools import GetImageAttributes
    from ImageTools import ImportImage
    
    RefIm = ImportImage(DicomDir)
    
    # Convert the pixel array to a list of indices-by-frame:
    IndsByFrame = PixArr2IndsByFrame(PixArr, FrameToSliceInds, Thresh=0.5)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir, Package='pydicom')
    
    PtsByCnt = []
    CntDataByCnt = []
    
    for Inds in IndsByFrame:
        Pts = Indices2Points(Indices=Inds, RefIm=RefIm)
        
        PtsByCnt.append(Pts)
        
        CntDataByCnt.append(Points2ContourData(Points=Pts))
    
    
    #print(f'\n\n\nlen(PtsByCnt) = {len(PtsByCnt)}')
    #print(f'\nlen(CntDataByCnt) = {len(CntDataByCnt)}')
    #print(f'\nPtsByCnt = {PtsByCnt}')
    #print(f'\nCntDataByCnt = {CntDataByCnt}')
 
    return PtsByCnt, CntDataByCnt










def PtsByCntByRoi2PixArrByRoi_OLD(PtsByCntByRoi, RefIm, LogToConsole=False):
    """
    Convert RTS data (e.g. from GetRtsDataOfInterest()) to a list of pixel
    arrays grouped by ROI.
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates to be copied from 
        Rts.
    
    RefIm : SimpleITK image
        The 3D image that relates to the RTS data in PtsByCntByRoi.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrByRoi : list of Numpy data arrays
        A list (for each ROI) of 3D pixel arrays.
    """
    
    PixArrByRoi = []
    
    # Convert points-by-contour-by-ROI to indices-by-contour-by-ROI:
    IndsByCntByRoi = []
    
    R = len(PtsByCntByRoi)
    
    for r in range(R):
        IndsByCnt = []
        
        if PtsByCntByRoi[r]:
            for Pts in PtsByCntByRoi[r]:
                Inds = Points2Indices(Pts, RefIm)
                    
                IndsByCnt.append(Inds)
            
        IndsByCntByRoi.append(IndsByCnt)
        
        # Convert IndsByCnt to PixArr:
        PixArr = IndsByContour2PixArr(IndsByCnt, RefIm)
        
        PixArrByRoi.append(PixArr)
        
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Results from PtsByCntByRoi2PixArrByRoi():')
        shapes = [PixArrByRoi[r].shape for r in range(R)]
        print(f'   PixArrByRoi.shape = {shapes}')
        print('-'*120)
    
    return PixArrByRoi






def PixArr2PtsByContour_OLD(PixArr, F2Sinds, DicomDir, Thresh=0.5):
    """
    04/02/21:
        This was modified (see PixArr2PtsByCnt()) since the function below 
        creates a list for each object in any given frame.  The revised 
        function treats each object as a separate contour, rather than nesting
        the contours.
        
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates.
 
    Inputs:
    ******
    
    PixArr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        PixArr.
    
    F2Sinds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
    
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
 
    
    Outputs:
    *******
    
    PtsByObjByFrame : list of a list of a list of floats
        A list (for each frame) of a list (for each object, i.e. contour) of a 
        list (for each dimension) of the physical coordinates that define the 
        mask in each frame in PixArr.
        
    CntDataByObjByFrame : list of a list of a list of strings
        A list (for each frame) of a list (for each object, i.e. contour) of a 
        flattened list of [x, y, z] physical coordinates that define the mask 
        in each frame in PixArr as strings (format of ContourData tag in RTS).
    """
    
    #from ImageTools import GetImageAttributes
    from ImageTools import ImportImage
    
    RefIm = ImportImage(DicomDir)
    
    # Convert the pixel array to a list of indices-by-object-by-frame:
    IndsByObjByFrame = PixArr2IndsByFrame(PixArr, F2Sinds, Thresh=0.5)
    
    #Size, Spacings, ST, IPPs, Dirs,\
    #ListOfWarnings = GetImageAttributes(DicomDir, Package='pydicom')
    
    PtsByObjByFrame = []
    CntDataByObjByFrame = []
    
    for f in range(len(IndsByObjByFrame)):
        PtsByObj = []
        CntDataByObj = []
    
        for o in range(len(IndsByObjByFrame[f])):
            Pts = Indices2Points(Indices=IndsByObjByFrame[f][o], 
                                 RefIm=RefIm)
            
            PtsByObj.append(Pts)
            
            CntDataByObj.append(Points2ContourData(Points=Pts))
        
        PtsByObjByFrame.append(PtsByObj)
        CntDataByObjByFrame.append(CntDataByObj)
    
    
    #print(f'\n\n\nlen(PtsByCnt) = {len(PtsByCnt)}')
    #print(f'\nlen(CntDataByCnt) = {len(CntDataByCnt)}')
    #print(f'\nPtsByCnt = {PtsByCnt}')
    #print(f'\nCntDataByCnt = {CntDataByCnt}')
 
    return PtsByObjByFrame, CntDataByObjByFrame







