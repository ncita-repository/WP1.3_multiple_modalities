# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:38:31 2020

@author: ctorti

Modified from:
    https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
"""



# Import packages and functions:
import numpy as np
import SimpleITK as sitk
from pydicom import dcmread
from itertools import product
from scipy.sparse import lil_matrix, issparse
from skimage.measure import find_contours
from shapely.geometry import MultiPolygon, Polygon, Point
#from warnings import warn

from GeneralTools import UnpackPoints

#from DicomTools import ImportDicoms
from DicomTools import GetDicomSOPuids

from ImageTools import GetImageAttributes
#from ImageTools import InitialiseImage
#from ImageTools import AddImages

#from RtsTools import GetCIStoDcmInds
from RtsTools import GetCStoDcmInds

#from SegTools import GetPFFGStoDcmInds




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

def Point2Index(Point, Origin, Directions, Spacings):
    """
    Convert a physical point (in Patient Coordinate System (PCS)) to an index
    (in Image Coordinate System (ICS).
    
    Input:
        Point      - (List of floats) Point in PCS, e.g. [x, y, z]
        
        Origin     - (List of floats) The 3D image origin 
                     (= ImagePositionPatient of the first slice in the DICOM 
                     series), e.g. [x0, y0, z0]
        
        Directions - (List of floats) The direction cosine along x (rows), 
                     y (columns) and z (slices) (= the cross product of the 
                     x and y direction cosines from ImageOrientationPatient 
                     appended to ImageOrientationPatient),
                     e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
        
        Spacings   - (List of floats) The pixel spacings along x, y and z 
                     (= SliceThickness appended to PixelSpacing), 
                     e.g. [di, dj, dk]
    
        
    Returns:
        Index      - (List of integers) Index (i.e. in ICS) of PointPCS,
                     e.g. [i, j, k]
        
    
    Equations:
    
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
    
    return [i, j, k]






def ConvertContourDataToIcsPoints(ContourData, DicomDir):
    """
    Convert contour data from a flat list of points in the Patient Coordinate
    System (PCS) to a list of [x, y, z] points in the Image Coordinate System. 
    
    Inputs:
        ContourData - (List of floats) Flat list of points in the PCS,
                      e.g. [x0, y0, z0, x1, y1, z1, ...]
        
        DicomDir    - (String) Directory containing DICOMs
        
        
    Returns:
        Indices   - (List of floats) List of a list of [x, y, z] coordinates
                      of each point in ContourData converted to the ICS,
                      e.g. [[x0, y0, z0], [x1, y1, z1], ...]
    """
    
    
    # Get the image attributes:
    Size, Spacings, SliceThick, IPPs, Directions = GetImageAttributes(DicomDir)
    
    # Initialise the array of contour points:
    Indices = []
    
    # Iterate for all contour points in threes:
    for p in range(0, len(ContourData), 3):
        # The point in Patient Coordinate System:
        pt = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        # Convert pt to Image Coordinate System:
        ind = Point2Index(pt, IPPs[0], Directions, Spacings)
        
        Indices.append(ind)
        
    return Indices







def PixArr2Labmap(PixelArray, NumOfSlices, FrameToSliceInds):
    """
    Convert a 3D SEG pixel array to a 3D SEG labelmap.  The labelmap will 
    contain the 2D pixel arrays at the same frame position as the slice
    position of the corresponding DICOM series (with 2D zero masks in-between).
    
    Inputs:
        PixelArray       - (Numpy array) The SEG's PixelData loaded as a 
                           Numpy array with dtype uint8
        
        NumOfSlices      - (Integer) The number of slices in the DICOM series
        
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number
        
    Returns:
        LabelMap         - (Numpy array) A zero-padded version of PixelArray
                           with dtype uint

    """
    
    #print('\nRunning PixArr2Labmap()...')
    
    NumOfFrames, NumOfRows, NumOfCols = PixelArray.shape
    
    #print(f'\n\nThe maximum value in PixelArray is {PixelArray.max()}')
    
    # Initialise LabelMap:
    #LabelMap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='bool')
    LabelMap = np.zeros((NumOfSlices, NumOfRows, NumOfCols), dtype='uint')
    
    # Over-write all frames in LabelMap with non-zero 2D masks from PixelArray:
    for i in range(len(FrameToSliceInds)):
        # The DICOM slice number for the i^th frame:
        s = FrameToSliceInds[i]
        
        #print(f'\nThe DICOM slice number for the {i}^th frame is {s}')
        
        LabelMap[s] = PixelArray[i]
        
    
    #print(f'\n\nThe maximum value in LabelMap is {LabelMap.max()}')
    
    #print('\nCompleted running PixArr2Labmap.')
    
    return LabelMap





def PixArr2Image(PixelArray, Image, FrameToSliceInds):
    """
    Convert a 3D SEG pixel array to a 3D SimpleITK image using an intermediate
    conversion from 3D SEG pixel array to a Numpy array (using PixArr2Labmap).  
    
    Inputs:
        PixelArray       - (Numpy array) The SEG's PixelData loaded as a 
                           Numpy array with dtype uint8
                           
        Image            - (SimpleITK image) The 3D image that PixelArray
                           relates to
        
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number
        
    Returns:
        LabelMapImage    - (SimpleITK image) A zero-padded version of 
                           PixelArray

    """
    
    #print('\nRunning PixArr2Image()...')
    
    ImSize = Image.GetSize()
    
    nda = PixArr2Labmap(PixelArray, ImSize[2], FrameToSliceInds)
    
    #print('\nContinuing with running PixArr2Image()..')
    
    LabelMapImage = sitk.GetImageFromArray(nda)
    
    # Set the Origin, Spacing and Direction of LabelMapImage to that of Image:
    LabelMapImage.SetOrigin(Image.GetOrigin())
    LabelMapImage.SetSpacing(Image.GetSpacing())
    LabelMapImage.SetDirection(Image.GetDirection())
    
    #print(f'\n\nThe maximum value in nda is {nda.max()}')
    #print(f'\n\nThe maximum value in LabelMapImage is {ImageMax(LabelMapImage)}')
    
    #print('\nCompleted running PixArr2Image().')
    
        
    return LabelMapImage




def Image2PixArr(LabmapImage):
    """
    Convert a 3D SimpleITK image to a 3D SEG pixel array.  
    
    Inputs:
        LabelMapImage    - (SimpleITK image) A zero-padded version of 
                           PixelArray
        
        
    Returns:
        PixelArray       - (Numpy array) The SEG's PixelData as a Numpy array
                           with dtype uint8
                           
        FrameToSliceInds - (List of integers) The DICOM slice numbers that
                           correspond to each frame in PixelArray,
                           e.g. PFFGStoDcmInds 
                           = PerFrameFunctionalGroupsSequence-to-DICOM slice
                           number

    """
    
    # Convert from SimpleITK image to a Numpy array:
    LabmapNda = sitk.GetArrayFromImage(LabmapImage)
    
    FrameToSliceInds = []
    
    # Loop through all frames in LabmapNda:
    for i in range(LabmapNda.shape[0]):
        Sum = np.amax(LabmapNda[i])
        
        if Sum:
            # This is a non-zero frame:
            FrameToSliceInds.append(i)
            
            
    
    # Initialise PixelArray:
    PixelArray = np.zeros((len(FrameToSliceInds), 
                           LabmapNda.shape[1], 
                           LabmapNda.shape[2]))
        
    # Use FrameToSliceInds to index the non-zero frames in LabmapNda:
    for i in range(len(FrameToSliceInds)):
        PixelArray[i] = LabmapNda[FrameToSliceInds[i]]
        
    
    return PixelArray, FrameToSliceInds
        




def PixArrForOneFrame(PixelArray, FrameNum):
    """
    Modify a 3D SEG pixel array so that all but one frame has non-zero 2D masks
    (i.e. set all other frames to zero).  
    
    Inputs:
        PixelArray    - (Numpy array) The SEG's PixelData loaded as a Numpy
                        array with dtype uint8
        
        FrameNum      - (Integer) The number of slices in the DICOM series
        
    Returns:
        NewPixelArray - (Numpy array) New PixelArray

    """
    
    NumOfFrames, NumOfRows, NumOfCols = PixelArray.shape
    
    # Initialise NewPixelArray:
    NewPixelArray = np.zeros((1, NumOfRows, NumOfCols), dtype='uint')
    
    NewPixelArray[0] = PixelArray[FrameNum]
            
    return NewPixelArray






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
        raise Exception("The inputs must both be 2D or 3D lists of points/indices")







def CentroidOfContour(Contour):
    """
    Compute the centroid of a list of points (or indices).  
    
    Inputs:
        Contour  - (List of list of floats/integers) List of a list of points 
                   or indices
        
    Returns:
        Centroid - (List of floats/integers) Centroid of Contour

    """
    
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








def GetPtsInContours(RtsFpath, DicomDir):
    
    # Import the RTSTRUCT ROI:
    Roi = dcmread(RtsFpath)
    
    # Get the DICOM SOP UIDs:
    SopUids = GetDicomSOPuids(DicomDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices:
    #CIStoDcmInds = GetCIStoDcmInds(Roi, SopUids)
    
    # Get the ContourSequence-to-DICOM slice indices:
    CStoDcmInds = GetCStoDcmInds(Roi, SopUids)
    
    #Ncontours = len(CIStoDcmInds)
    Ncontours = len(CStoDcmInds)
    
    
    # Get a list of the Contour Sequences in the Roi:
    ContourSequences = Roi.ROIContourSequence[0].ContourSequence
    
    # Get a list of all Referenced SOP Instance UIDs from the ContourSequences:
    #RefSopUids = [sequence.ContourImageSequence[0]\
    #              .ReferencedSOPInstanceUID for sequence in ContourSequences]
    
    
    
    Contours = []
    
    
    # Loop through each CStoDcmInds:
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
    
    
        Contours.append(points)

    return Contours





def GetIndsInContours(RtsFpath, DicomDir):
    
    PtsInContours = GetPtsInContours(RtsFpath, DicomDir)
    
    Size, Spacings, SliceThick,\
    IPPs, Dirs = GetImageAttributes(DicomDir=DicomDir, Package='pydicom')
    
    IndsInContours = []
    
    for contour in PtsInContours:
        Indices = []
        
        for pointPCS in contour:
            inds = Point2Index(Point=pointPCS, Origin=IPPs[0], 
                               Directions=Dirs, Spacings=Spacings)
            
            # Convert fractional indices to integers:
            inds = [int(round(item)) for item in inds]
            
            Indices.append(inds)
            
        
        IndsInContours.append(Indices)
        
    
    return IndsInContours







def Inds2SparseLabmap(Indices, RefImage): # not useful since Sparse2Labelmap doesn't work
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
            
    






def Sparse2Labelmap(Sparse): # doesn't work
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










def IndsInContours2PixArr(IndsInContours, RefImage):
    """
    Poly2Mask(polygons, im_size):
    
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    
    Converts polygons to a sparse binary mask.
 
    >>> from sima.ROI import poly2mask
    >>> poly1 = [[0,0], [0,1], [1,1], [1,0]]
    >>> poly2 = [[0,1], [0,2], [2,2], [2,1]]
    >>> mask = poly2mask([poly1, poly2], (3, 3))
    >>> mask[0].todense()
    matrix([[ True, False, False],
            [ True,  True, False],
            [False, False, False]], dtype=bool)
 
    Parameters
    ----------
    polygons : sequence of coordinates or sequence of Polygons
        A sequence of polygons where each is either a sequence of (x,y) or
        (x,y,z) coordinate pairs, an Nx2 or Nx3 numpy array, or a Polygon
        object.
    im_size : tuple
        Final size of the resulting mask
 
    Output
    ------
    mask
        A list of sparse binary masks of the points contained within the
        polygons, one mask per plane
 
    """
    #if len(im_size) == 2:
    #    im_size = (1,) + im_size
 
    #polygons = _reformat_polygons(polygons)
    
    ImSize = RefImage.GetSize()
    
    Ncontours = len(IndsInContours)
    
    #masks = []
    PixArr = np.zeros((Ncontours, ImSize[1], ImSize[0]))
    
    # Note that IndsInContours are a list of indices for each contour:
    for c in range(Ncontours):
        Indices = IndsInContours[c]
        
        # Only proceed if there are at least 3 indices in this contour:
        if len(Indices) > 2:
            # The z index should be the same for all indices:
            #z = Indices[0][2]
            
            #mask = np.zeros((1, ImSize[1], ImSize[0]))
            
            # Convert list Indices to a Shapely Polygon:
            poly = Polygon(Indices)
                
            x_min, y_min, x_max, y_max = poly.bounds
         
            points = [Point(x, y) for x, y in
                      product(np.arange(int(x_min), np.ceil(x_max)),
                              np.arange(int(y_min), np.ceil(y_max)))]
                      
            points_in_poly = list(filter(poly.contains, points))
            
            for point in points_in_poly:
                xx, yy = point.xy
                
                x = int(xx[0])
                y = int(yy[0])
                
                if 0 <= y < ImSize[1] and 0 <= x < ImSize[0]:
                    #mask[z, y, x] = 1
                    #mask[c, y, x] = 1
                    PixArr[c, y, x] = 1
               
        #masks.append(mask)
        #PixArr[c] = mask
        
    #for z_coord in np.arange(mask.shape[0]):
    #    masks.append(lil_matrix(mask[z_coord, :, :]))
        
    return PixArr