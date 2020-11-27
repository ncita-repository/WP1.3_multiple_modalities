# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:38:31 2020

@author: ctorti

Modified from:
    https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
"""



# Import packages and functions:
#import numpy as np
#import SimpleITK as sitk
#from pydicom import dcmread
#from itertools import product
#from scipy.sparse import lil_matrix, issparse
#from skimage.measure import find_contours
#from shapely.geometry import MultiPolygon, Polygon, Point
##from warnings import warn

#from GeneralTools import UnpackPoints
#from GeneralTools import UniqueItems

##from DicomTools import ImportDicoms
#from DicomTools import GetDicomSOPuids

#from ImageTools import GetImageAttributes
##from ImageTools import InitialiseImage
##from ImageTools import AddImages

##from RtsTools import GetCIStoSliceInds
#from RtsTools import GetCStoSliceInds

##from SegTools import GetPFFGStoSliceInds





"""
******************************************************************************
******************************************************************************
CONVERSION FUNCTIONS
******************************************************************************
******************************************************************************
"""       


"""
Note function Point2Index is somewhat redundant since a point can be converted
to an index with a single line of code using SimpleITK:
    
    SitkIm.TransformPhysicalPointToIndex(Point)
"""


def GetPtsInContour(RtsFpath, DicomDir, RoiNum, SliceNum):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in a 
    contour in a ROI.
    
    Inputs:
    ------
    
    RtsFpath : string
        Full path of the DICOM-RTSTRUCT file.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
        
    RoiNum : integer
        The ROI number (counting from 0) that contains the contour of interest.
        
    SliceNum : integer
        The slice number (counting from 0) that corresponds to the contour of
        interest.
                       
                            
    Outputs:
    -------
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points. 
    """
    
    from pydicom import dcmread
    from DicomTools import GetDicomSOPuids
    from RtsTools import GetCStoSliceInds
    
    # Import the RTSTRUCT ROI:
    Rts = dcmread(RtsFpath)
    
    # Get the DICOM SOP UIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourSequence-to-slice indices:
    CStoSliceInds = GetCStoSliceInds(Rts, SopUids)
    
    # The contour number of interest:
    ContourNum = CStoSliceInds.index(SliceNum)
    
    # The ROI Contour Sequence:
    RoiContourSequence = Rts.ROIContourSequence
    
    # The number of ROIs:
    NumOfRois = len(RoiContourSequence)
    
    if RoiNum > NumOfRois - 1:
        raise Exception(f"There are only {NumOfRois} ROIs in the RTS.\n",
                        f"RoiNum = {RoiNum} > {NumOfRois}.")
    
    
    # The Contour Sequences in the ROI of interest:
    ContourSequence = Rts.ROIContourSequence[RoiNum].ContourSequence
    
    # The number of contours in this ROI:
    NumOfContours = len(ContourSequence)
    
    if ContourNum > NumOfContours - 1:
        raise Exception(f"There are only {NumOfContours} contours in the RTS.",
                        f"\nContourNum = {ContourNum} > {NumOfContours}.")
        
    ContourData = ContourSequence[ContourNum].ContourData
    
    ContourData = [float(item) for item in ContourData]
    
    Points = []
            
    # Iterate for all points in threes:
    for p in range(0, len(ContourData), 3):
        point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        Points.append(point)
        
    return Points


        




    

def GetPtsByRoi(RtsFpath, DicomDir):
    """
    Get the physical points (i.e. in the Patient Coordinate System) in all 
    contours for all ROIs.
    
    Inputs:
    ------
    
    RtsFpath : string
        Full path of the DICOM-RTSTRUCT file.
        
    DicomDir : string 
        Directory containing the corresponding DICOMs.
                       
                            
    Outputs:
    -------
    
    PtsByRoi : list of list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each 
        dimension) of points.
    """
    
    from pydicom import dcmread
    from DicomTools import GetDicomSOPuids
    from RtsTools import GetCStoSliceInds
    
    # Import the RTSTRUCT ROI:
    Roi = dcmread(RtsFpath)
    
    # Get the DICOM SOP UIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourImageSequence-to-slice indices:
    #CIStoSliceInds = GetCIStoSliceInds(Roi, SopUids)
    
    # Get the ContourSequence-to-slice indices:
    CStoSliceInds = GetCStoSliceInds(Roi, SopUids)
    
    Nrois = len(CStoSliceInds)
    
    PtsByRoi = []
    
    for r in range(Nrois):
        # Number of contours in this ROI:
        Ncontours = len(CStoSliceInds[r])
        
        # Get a list of the Contour Sequences in this ROI:
        ContourSequences = Roi.ROIContourSequence[r].ContourSequence
        
        # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
        #RefSopUids = [sequence.ContourImageSequence[0]\
        #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
        
        
        PtsByContour = []
        
        
        # Loop through each CStoSliceInds:
        for c in range(Ncontours):
            
            # Get the indeces of all matching RefSopUids:
            #ind = RefSopUids.index(SopUid) # this will only find the first 
            # match, and there may be more than one!
            #inds = [i for i, e in enumerate(RefSopUids) if e==SopUid]
            
            # Iterate for each index in inds:
            #for ind in inds:
            
            ContourSequence = ContourSequences[c]
            
            ContourData = [float(item) for item in ContourSequence.ContourData]
    
            # Initialise the array of points for this contour:
            points = []
            
            # Iterate for all points in threes:
            for p in range(0, len(ContourData), 3):
                point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
                
                points.append(point)
        
        
            PtsByContour.append(points)
            
            
        PtsByRoi.append(PtsByContour)

    return PtsByRoi









def Point2Index(Point, Origin, Directions, Spacings):
    """
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
    
        
    Ouputs:
    ------
    
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
    
    # Convert fractional indices to integers:
    i = int(round(i))
    j = int(round(j))
    k = int(round(k))
    
    return [i, j, k]






def ContourToIndices(Points, Origin, Directions, Spacings):
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
    
        
    Ouputs:
    ------
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        (in ICS) of Points, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    Indices = []
    
    for point in Points:
        Indices.append(Point2Index(point, Origin, Directions, Spacings))
        
    return Indices




    
    
    
    
def ConvertContourDataToIndices(ContourData, DicomDir):
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






def GetIndsByRoi(RtsFpath, DicomDir):
    """
    Get the indices of the points (i.e. in the Image Coordinate System) in all 
    contours for all ROIs.
    
    
    Inputs:
        RtsFpath   - (String) Full path of the DICOM-RTSTRUCT file 
        
        DicomDir   - (String) Directory containing the corresponding DICOMs 
                       
                            
    Returns:
        IndsInRois - (List of list of a list of floats) List (for each ROI) of
                     a list (for all contours) of a list (for each dimension)
                     of indices 
    
    """
    
    from ImageTools import GetImageAttributes
    
    PtsByRoi = GetPtsByRoi(RtsFpath, DicomDir)
    
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir)
    
    
    Nrois = len(PtsByRoi)
    
    IndsByRoi = []
    
    for r in range(Nrois):
        IndsByContour = []
        
        # The list of points in the contours for this ROI:
        PtsByContour = PtsByRoi[r]
        
        for contour in PtsByContour:
            Indices = ContourToIndices(contour, IPPs[0], Dirs, Spacings)
                
            IndsByContour.append(Indices)
            
        IndsByRoi.append(IndsByContour)
        
    
    return IndsByRoi








def ChangeZinds(Indices, NewZind):
    """
    Re-assign the z-components of a list of indeces.
    
    Inputs:
    ------
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of indices, 
        e.g. [[i0, j0, k0], [i1, j1, k1], ...].
        
    NewZind : integer
        The value to re-assign to the z-components (i.e. k) of Indeces.
    
    
        
    Ouputs:
    ------
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        with modified z-components, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    if not isinstance(NewZind, int):
        raise Exception(f"NewZind = {NewZind} must be an integer.")
    
    for i in range(len(Indices)):
        Indices[i][2] = NewZind
        
    return Indices








def Ind2Pt(Index, Origin, Directions, Spacings):
    """
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
    
    x = Origin[0] + Spacings[0]*Index[0]*Directions[0] \
        + Spacings[0]*Index[0]*Directions[1] \
        + Spacings[0]*Index[0]*Directions[2] \
    
    y = Origin[1] + Spacings[1]*Index[1]*Directions[3] \
        + Spacings[1]*Index[1]*Directions[4] \
        + Spacings[1]*Index[1]*Directions[5] \
        
    z = Origin[2] + Spacings[2]*Index[2]*Directions[6] \
        + Spacings[2]*Index[2]*Directions[7] \
        + Spacings[2]*Index[2]*Directions[8] \
        
    return [x, y, z]







def Inds2Pts(Indices, Origin, Directions, Spacings):
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
        Points.append(Ind2Pt(index, Origin, Directions, Spacings))
        
    return Points







def Pts2ContourData(Points):
    """
    Format a list of points to a flat list as required for ContourData.
    
    Inputs:
    ------
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points in the
        PCS, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
        
    Outputs:
    -------
    
    ContourData : list of strings
        Flat list of coordinates in Points converted from floats to strings.
    """
    
    ContourData = []
    
    for point in Points:
        point = [str(item) for item in point]
        
        ContourData.extend(point)
        
    return ContourData







def IndsByContour2PixArr(IndsByContour, RefImage):
    """
    Convert a list of indices grouped by contour to a pixel array.
       
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    Inputs:
    ------
    
    IndsByContour : list of lists of lists of integers
        A list (for each contour) of a list (for all polygons) of a list (for 
        each dimension) of indices.
                     
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
        
    Outputs:
    -------
    
    PixArr : Numpy data array
        A list of 2D masks of the indices - one mask per contour.
    """
    
    import numpy as np
    from shapely.geometry import Polygon, Point
    from itertools import product
    
    ImSize = RefImage.GetSize()
    
    Ncontours = len(IndsByContour)
    
    PixArr = np.zeros((Ncontours, ImSize[1], ImSize[0]))
    
    # Note that IndsByContour are a list of indices for each contour:
    for c in range(Ncontours):
        Indices = IndsByContour[c]
        
        # Only proceed if there are at least 3 indices in this contour (the 
        # minimum number required to define a closed contour):
        if len(Indices) > 2:
            # Convert list Indices to a Shapely Polygon:
            Poly = Polygon(Indices)
                
            xMin, yMin, xMax, yMax = Poly.bounds
         
            points = [Point(x, y) for x, y in
                      product(np.arange(int(xMin), np.ceil(xMax)),
                              np.arange(int(yMin), np.ceil(yMax)))]
                      
            PointsInPoly = list(filter(Poly.contains, points))
            
            for point in PointsInPoly:
                xx, yy = point.xy
                
                x = int(xx[0])
                y = int(yy[0])
                
                if 0 <= y < ImSize[1] and 0 <= x < ImSize[0]:
                    PixArr[c, y, x] = 1
    
        
    return PixArr








def GetLabmapImsByRoi(RtsFpath, DicomDir, CStoSliceIndsByRoi, RefIm):
    """
    Get a list of 3D labelmap (SimpleITK) images - one per ROI.  
    
    Inputs:
    ------
    
    RtsFpath : string
        Full path of the DICOM-RTSTRUCT file.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    CStoSliceIndsByRoi : list of a list of integers
        List of slice numbers that correspond to each Referenced SOP instance 
        UID in the Contour Sequence - one list per ROI.
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    -------
        
    LabmapIms : List of SimpleITK images
        A list of 3D zero-padded labelmap (SimpleITK) images - one per ROI.
    """
    
    LabmapIms = []
    
    IndsByRoi = GetIndsByRoi(RtsFpath, DicomDir)
    
    for RoiNum in range(len(IndsByRoi)):
        PixArr = IndsByContour2PixArr(IndsByContour=IndsByRoi[RoiNum], 
                                      RefImage=RefIm)
    
        LabmapIm = PixArr2Labmap(PixArr=PixArr, 
                                 NumOfSlices=RefIm.GetSize()[2], 
                                 FrameToSliceInds=CStoSliceIndsByRoi[RoiNum])
        
        LabmapIms.append(LabmapIm)
    
    return LabmapIms








def PixArr2Labmap(PixArr, NumOfSlices, FrameToSliceInds):
    """
    Convert a 3D SEG pixel array to a 3D SEG labelmap.  The labelmap will 
    contain the 2D pixel arrays at the same frame position as the slice
    position of the corresponding DICOM series (with 2D zero masks in-between).
    
    Inputs:
    ------
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
        
    NumOfSlices : integer
        The total number of frames in the labelmap (e.g. the number of slices 
        in a image series).
        
    FrameToSliceInds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number.
        
    Outputs:
    -------
    
    Labmap : Numpy array
        Zero-padded version of PixArr.
    """
    
    import numpy as np
    
    NumOfFrames, NumOfRows, NumOfCols = PixArr.shape
    
    # Initialise Labmap:
    #Labmap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='bool')
    Labmap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='uint')
    
    #print(f'\nPixArr.shape = {PixArr.shape}')
    ##print(f'\n\nThe maximum value in PixArr is {PixArr.max()}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    #print(f'Labmap.shape = {Labmap.shape}')
    
    # Over-write all frames in Labmap with non-zero 2D masks from PixArr:
    for i in range(len(FrameToSliceInds)):
        # The slice number for the i^th frame:
        s = FrameToSliceInds[i]
        
        #print(f'\nThe slice number for the {i}^th frame is {s}')
        
        Labmap[s] = PixArr[i]
        
    
    #print(f'\n\nThe maximum value in Labmap is {Labmap.max()}')
    
    return Labmap





def PixArr2Image(PixArr, FrameToSliceInds, RefIm):
    """
    Convert a 3D pixel array to a 3D SimpleITK image.  
    
    Inputs:
    ------
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
    
    FrameToSliceInds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    -------
        
    LabmapIm : SimpleITK image
        A zero-padded version of PixArr.
    """
    
    import SimpleITK as sitk
    
    ImSize = RefIm.GetSize()
    
    #print(f'\nImSize = {ImSize}')
    #print(f'PixArr.shape = {PixArr.shape}')
    #print(f'ImSize[2] = {ImSize[2]}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    
    nda = PixArr2Labmap(PixArr, ImSize[2], FrameToSliceInds)
    
    LabmapIm = sitk.GetImageFromArray(nda)
    
    # Set the Origin, Spacing and Direction of LabmapIm to that of RefIm:
    LabmapIm.SetOrigin(RefIm.GetOrigin())
    LabmapIm.SetSpacing(RefIm.GetSpacing())
    LabmapIm.SetDirection(RefIm.GetDirection())
    
    #print(f'\n\nThe maximum value in nda is {nda.max()}')
    #print(f'\n\nThe maximum value in LabmapIm is {ImageMax(LabmapIm)}')
        
    return LabmapIm






def PixArr2ImagesBySeg(PixArr, FrameToSliceIndsBySeg, FrameNumsBySeg, RefIm):
    """
    Convert a 3D pixel array to a list of 3D SimpleITK images - one per
    segment.  
    
    Inputs:
    ------
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
    
    FrameToSliceIndsBySeg : List of list of integers
        List of the DICOM slice numbers that correspond to each frame in PixArr 
        (e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number) - a list for each segment.
    
    FrameNumsBySeg : List of list of integers
        List of the frame numbers that belong to each segment - one list per
        segment.
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    -------
        
    LabmapIms : List of SimpleITK images
        A list of 3D zero-padded labelaps - one per segment.
    """
    
    LabmapIms = []
    
    #print(f'\nFrameToSliceIndsBySeg = {FrameToSliceIndsBySeg}\n')
    #print(f'FrameNumsBySeg = {FrameNumsBySeg}')
    
    for i in range(len(FrameNumsBySeg)):
        # The starting (= a) and ending (= b) frame numbers:
        a = FrameNumsBySeg[i][0]
        b = FrameNumsBySeg[i][-1]
        
        
        #print(f'a = {a}, b = {b}')
        #print(f'PixArr.shape = {PixArr.shape}')
        #print(f'PixArr[{a}:{b+1}].shape = {PixArr[a:b+1].shape}')
        
        LabmapIm = PixArr2Image(PixArr=PixArr[a:b+1],
                                RefIm=RefIm,
                                FrameToSliceInds=FrameToSliceIndsBySeg[i])
        
        LabmapIms.append(LabmapIm)
        
    return LabmapIms






def Image2PixArr(LabmapIm, Non0FrameInds=None):
    """
    Convert a 3D SimpleITK image to a 3D SEG pixel array.  
    
    Inputs:
    ------
    
    LabelmapIm : SimpleITK image
        A zero-padded version of PixelData.
        
    Non0FrameInds : List of integers (optional; None by default)  
        Use PerFrameFunctionalGroupsSequence-to-slice indices 
        (PFFGStoSliceInds) if they are known. If this input is not specified a 
        list of the frame numbers that are non-zero will be determined.
        
        
    Outputs:
    -------
        
    PixArr : Numpy array
        The non-zero frames in LabmapIm as a Numpy array.
                           
    Non0FrameInds : List of integers
        List of the indices of non-zero frames in LabelmapIm.  This should be
        equivalent to the PerFrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds).
    """
    
    import SimpleITK as sitk
    import numpy as np
    
    # Convert from SimpleITK image to a Numpy array:
    Nda = sitk.GetArrayFromImage(LabmapIm)
    
    if not Non0FrameInds:
        Non0FrameInds = []
        
        #print('')
        
        # Loop through all frames in Nda and get the sum of all pixels:
        for i in range(Nda.shape[0]):
            Sum = np.amax(Nda[i])
            
            #print(f'i = {i}, Sum = {Sum}')
            
            if Sum:
                # This is a non-zero frame:
                Non0FrameInds.append(i)
    
    # Initialise PixArr:
    PixArr = np.zeros((len(Non0FrameInds), Nda.shape[1], Nda.shape[2]))
        
    # Use Non0FrameInds to index the non-zero frames in Nda:
    for i in range(len(Non0FrameInds)):
        PixArr[i] = Nda[Non0FrameInds[i]]
        
    
    return PixArr, Non0FrameInds
        





def GetFrameFromPixArr(PixArr, FrameNum):
    """
    Extract a single (2D) frame from a 3D pixel array.  
    
    Inputs:
    ------
    
    PixArr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
        
    FrameNum : integer
        The number of slices in the DICOM series.
        
    
    Outputs:
    -------
    
    Frame : Numpy array
        Frame from PixArr.

    """
    
    import numpy as np
    
    F, R, C = PixArr.shape
    
    # Initialise Frame:
    Frame = np.zeros((1, R, C), dtype='uint')
    
    Frame[0] = PixArr[FrameNum]
            
    return Frame







def GetPhysicalShiftBetweenSlices(Image, SliceNum0, SliceNum1):
    """
    Get the physical shift between Slice0 and Slice1 in an image.
    
    Inputs:
    ------
    
    Image : SimpleITK image
    
    SliceNum0 : integer
    
    SliceNum1 : integer
    
    
    Outputs:
    -------
    
    mmShift : Numpy array
        Shift in mm between index [0,0,SliceNum1] and index [0,0,SliceNum0].
    """
    
    import numpy as np
    
    mmShift = np.array(Image.TransformIndexToPhysicalPoint([0,0,SliceNum1])) \
            - np.array(Image.TransformIndexToPhysicalPoint([0,0,SliceNum0]))
            
    return mmShift





def GetPixelShiftBetweenSlices(Image, SliceNum0, SliceNum1, Fractional=False):
    """
    Get the shift between Slice0 and Slice1 in an image in pixels.
    
    Inputs:
    ------
    
    Image : SimpleITK image
    
    SliceNum0 : integer
    
    SliceNum1 : integer
    
    Fractional : boolean (optional; False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    
    
    Outputs:
    -------
    
    PixShift : Numpy array
        Shift in pixels between index [0, 0, SliceNum1] and index 
        [0, 0, SliceNum0].
    """
    
    import numpy as np
    
    mmShift = GetPhysicalShiftBetweenSlices(Image, SliceNum0, SliceNum1)
    
    PixSpacing = np.array(Image.GetSpacing())
    
    PixShift = mmShift/PixSpacing
    
    if not Fractional:
        PixShift = np.round(PixShift).astype('int')
            
    return PixShift






    
def ShiftFrame(Frame, PixShift):
    """
    Shift the items in a 2D frame.
    
    Inputs:
    ------
    
    Frame : Numpy array
        A 2D frame from PixelData with shape (1, R, C), where R = number of 
        rows, and C = number of columns.
        
    PixShift : Numpy array
        Shift in pixels (e.g. [di, dj, dz].
        
    
    Outputs:
    -------
    
    ShiftedFrame : Numpy array
        Shifted frame with shape (1, R, C). Only the x- and y-components are 
        shifted.

    """
    
    import numpy as np
    
    F, R, C = Frame.shape
    
    if F > 1:
        msg = f"'Frame' must be a 2D frame with shape (1, R, C) but has shape"\
              + f" ({F}, {R}, {C})."
        
        raise Exception(msg)
    
    # Initialise ShiftedFrame:
    ShiftedFrame = np.zeros((1, R, C), dtype='uint')
    #ShiftedFrame = np.empty_like(Frame, dtype='uint') # this creates 42,932
    # unique values for some reason!
    
    #unique = UniqueItems(Nda=Frame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in Frame')
    #unique = UniqueItems(Nda=ShiftedFrame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in the initialised',
    #      f'ShiftedFrame: {unique[:11]}...')
    
    di, dj, dk = PixShift
    
    ##ShiftedFrame[0, dj:, di:] = Frame[0, :-(1+dj), :-(1+di)]
    ##ShiftedFrame[0, :-(1+dj), :-(1+di)] = Frame[0, dj:, di:]
    #ShiftedFrame[0, :R-dj, :C-di] = Frame[0, dj:, di:]
    
    if di > 0 and dj > 0:
        ShiftedFrame[0, dj:, di:] = Frame[0, :-dj, :-di]
        
    elif di < 0 and dj < 0:
        ShiftedFrame[0, :dj, :di] = Frame[0, -dj:, -di:]
        
    elif di > 0 and dj < 0:
        ShiftedFrame[0, :dj, di:] = Frame[0, -dj:, :-di]
        
    elif di < 0 and dj > 0:
        ShiftedFrame[0, dj:, :di] = Frame[0, :-dj, -di:]
        
    elif di == 0 and dj > 0:
        ShiftedFrame[0, dj:, :] = Frame[0, :-dj, :]
        
    elif di == 0 and dj < 0:
        ShiftedFrame[0, :dj, :] = Frame[0, -dj:, :]
        
    elif di > 0 and dj == 0:
        ShiftedFrame[0, :, di:] = Frame[0, :, :-di]
        
    elif di < 0 and dj == 0:
        ShiftedFrame[0, :, :di] = Frame[0, :, -di:]
        
    elif di == 0 and dj == 0:
        ShiftedFrame[0] = Frame[0]
        
    #unique = UniqueItems(Nda=ShiftedFrame, NonZero=False)
    #print(f'\n---> There are {len(unique)} unique items in the ShiftedFrame',
    #      'after shifting.')
            
    return ShiftedFrame







def Distance(item0, item1):
    """
    Compute the distance between points or indices.  
    
    Inputs:
        item0    - (List of floats/integers) A 2D/3D point or point indices
        
        item1    - (List of floats/integers) A 2D/3D point or point indices
        
    Returns:
        distance - (Scalar float/integer) Distance between item0 and item1 in
                   Euclidean space or in pixels

    """
    
    SumOfDims = len(item0) + len(item1)
    
    if SumOfDims == 4:
        # 2D points/indices.
        distance = ((item0[0] - item1[0])**2 \
                    + (item0[1] - item1[1])**2)**0.5
                    
        return distance
        
    elif SumOfDims == 6:
        # 3D points/indices.
        distance = ((item0[0] - item1[0])**2 \
                    + (item0[1] - item1[1])**2 \
                    + (item0[2] - item1[2])**2)**0.5
                    
        return distance
        
        
    else:
        msg = "The inputs must both be 2D or 3D lists of points/indices."
        
        raise Exception(msg)







def CentroidOfContour(Contour):
    """
    Compute the centroid of a list of points (or indices).  
    
    Inputs:
        Contour  - (List of list of floats/integers) List of a list of points 
                   or indices
        
    Returns:
        Centroid - (List of floats/integers) Centroid of Contour

    """
    
    import numpy as np
    
    Contour = np.array(Contour)
    
    return Contour.mean(axis=0)






def LabelContoursByCentroid(Contours): # incomplete
    """
    
    This function is incomplete 
    
    Given a list of points (or indices) assign an integer label to each contour 
    using the contour centroids as a best guess of which ROI they belong to.  
    
    Inputs:
        Contours - (List of list of a list of floats/integers) 
                   List of a list of points/indices for each contour
        
    Returns:
        RoiNums  - (List of integers) Integer label assigning each contour to
                   a ROI

    """
    
    Centroids = []
    
    for Contour in Contours:
        Centroids.append(CentroidOfContour(Contour))
        
        
    RoiNums = []
    RoiCentroids = []
    
    # Assign the first contour the number 0, and the first RoiCentroid 
    # Centroids[0]:
    RoiNums.append(0)
    RoiCentroids.append(Centroids[0])
    
    
    # Loop through Centroids starting from the 2nd one:
    for c in range(1, len(Centroids)):
        distances = []
        
        for r in range(len(RoiCentroids)):
            # Get the distance between the centroids but ignore the z-component:
            distance = Distance(Centroids[c][:2], RoiCentroids[r][:2])
            
            distances.append(distance)
            
            
        # Find index of minimum distance:
        ind = distances.index(min(distances))
        
        #if ind in RoiNums:
            
            
        
        
    return RoiNums
        





def Labelmap2Contours(Labelmap, threshold=0.5):
    """
    
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    Takes a mask and returns a MultiPolygon
 
    Parameters
    ----------
    mask : array
        Sparse or dense array to identify polygon contours within.
    threshold : float, optional
        Threshold value used to separate points in and out of resulting
        polygons. 0.5 will partition a boolean mask, for an arbitrary value
        binary mask choose the midpoint of the low and high values.
 
    Output
    ------
    MultiPolygon
        Returns a MultiPolygon of all masked regions.
 
    """
    
    import numpy as np
    from scipy.sparse import issparse
    from skimage.measure import find_contours
    
    ContourData = []
    ContourToSliceNums = []
    
    for z, Mask in enumerate(Labelmap):
        if issparse(Mask):
            Mask = np.array(Mask.astype('byte').todense())
 
        if (Mask != 0).sum() == 0:
            # If Mask is empty, just skip it
            continue
 
        # Add an empty row and column around the mask to make sure edge masks
        # are correctly determined
        ExpandedDims = (Mask.shape[0] + 2, Mask.shape[1] + 2)
        
        ExpandedMask = np.zeros(ExpandedDims, dtype=float)
        
        ExpandedMask[1:Mask.shape[0] + 1, 1:Mask.shape[1] + 1] = Mask
 
        Contours = find_contours(ExpandedMask.T, threshold)
 
        # Subtract off 1 to shift coords back to their real space:
        Contours = [np.subtract(x, 1).tolist() for x in Contours]
 
        v = []
        
        for poly in Contours:
            NewPoly = [point + [z] for point in poly]
            
            v.append(NewPoly)
            
        ContourData.extend(v)
        
        ContourToSliceNums.append(z)
 
    return ContourData, ContourToSliceNums









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










def Inds2Mask(Inds, RefImage):
    """
    Convert a list of indices to a (filled) mask (pixel array).
       
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    Inputs:
    ------
    
    Inds : list of a list of integers
        A list (for all points) of a list (for each dimension) of indices.
                     
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
        
    Outputs:
    -------
    
    Mask : Numpy data array
        A 2D mask of the indices.
    """
    
    import numpy as np
    from itertools import product
    from shapely.geometry import Polygon, Point
    
    # Only proceed if there are at least 3 indices in this contour (the 
    # minimum number required to define a closed contour):
    if len(Inds) < 3:
        raise Exception(f"There are only {len(Inds)} indices. A minimum of ",
                        "3 are required to define a closed contour.")
        
    
    ImSize = RefImage.GetSize()
    
    Mask = np.zeros((1, ImSize[1], ImSize[0]))
    
    # Convert Inds to a Shapely Polygon:
    Poly = Polygon(Inds)
        
    xMin, yMin, xMax, yMax = Poly.bounds
 
    points = [Point(x, y) for x, y in 
              product(np.arange(int(xMin), np.ceil(xMax)),
                      np.arange(int(yMin), np.ceil(yMax)))]
              
    PointsInPoly = list(filter(Poly.contains, points))
    
    for point in PointsInPoly:
        xx, yy = point.xy
        
        x = int(xx[0])
        y = int(yy[0])
        
        if 0 <= y < ImSize[1] and 0 <= x < ImSize[0]:
            Mask[1, y, x] = 1
    
        
    return Mask