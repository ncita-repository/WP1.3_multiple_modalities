# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:38:31 2020

@author: ctorti

    
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




def Point2Index(Point, RefIm, Rounding=True):
    """
    Convert a physical point (in the Patient Coordinate System) to an index
    (in the Image Coordinate System).
    
    Inputs:
    ******
    
    Point : list of floats
        A list (for each dimension) of coordinates, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
        
    RefIm : SimpleITK image
        A 3D image that occupies the grid that Point belongs to.
        
    Rounding : boolean (optional; True by default)
        If True, the integers in Index will be rounded to the nearest integers.
        
        
    Ouputs:
    ******
    
    Index : list of integers or floats
        List (for each dimension) of the indices of Points, e.g. 
            [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    if Rounding:
        Index = RefIm.TransformPhysicalPointToIndex(Point)
    else:
        Index = RefIm.TransformPhysicalPointToContinuousIndex(Point)
    
    #return Index # 29/01/2021
    return list(Index) # 29/01/2021







def Points2Indices(Points, RefIm, Rounding=True):
    """
    Convert a list of physical points (in the Patient Coordinate System) to a
    list of indices (in the Image Coordinate System).
    
    Inputs:
    ******
    
    Points : list of a list of floats
        A list (for each point) of a list (for each dimension) of coordinates, 
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
    RefIm : SimpleITK image
        A 3D image that occupies the grid that Points belong to.
        
    Rounding : boolean (optional; True by default)
        If True, the integers in Index will be rounded to the nearest integers.
    
        
    Ouputs:
    ******
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        of Points, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    Indices = []
    
    #print(f'\n\n\nPoints = {Points}\n\n')
    
    for Point in Points:
        #print(f'Point = {Point}')
        
        Index = Point2Index(Point, RefIm, Rounding)
        
        #print(f'Index = {Index}')
        
        Indices.append(Index)
        
    return Indices






def Index2Point(Index, RefIm):
    """
    Convert an index (in the Image Coordinate System) to a physical point (in 
    the Patient Coordinate System).
    .
    
    Inputs:
    ******
    
    Index : list of integers or floats
        A list (for each dimension) of an index, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    
    RefIm : SimpleITK image
        A 3D image that occupies the grid that Index belongs to.
        
        
    Ouputs:
    ******
    
    Point : list of floats
        A list (for each dimension) of the coordinates of Index, e.g.
        [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    from GeneralTools import DataTypeOfElementsInList
    
    Dtypes = DataTypeOfElementsInList(Index)
    
    if len(Dtypes) == 1:
        Dtype = Dtypes[0]
        
        """All items in Index are of a common data type."""
        if Dtype == int:
            Point = RefIm.TransformIndexToPhysicalPoint(Index)
        
        elif Dtype == float:
            Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
            
        else:
            msg = f"The data type of Index is {Dtype}. It must be either "\
                  + "'int' or 'float'."
            
            raise Exception(msg)
            
    elif len(Dtypes) == 2 and int in Dtypes and float in Dtypes:
        """Allow for the possibility of a mixture of integers and floats."""
        
        # Convert integers to float:
        Index = [float(item) for item in Index]
        
        Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
        
    else:
        msg = f"The data type of Index is {Dtypes}. It must be either "\
              + "'int' or 'float'."
         
        
    #return Point # 29/01/2021
    return list(Point) # 29/01/2021






def Indices2Points(Indices, RefIm):
    """
    Convert a list of indices (in the Image Coordinate System) to a list of 
    physical points (in the Patient Coordinate System).
    
    Inputs:
    ******
    
    Indices : list of a list of integers or floats
        A list (for each index) of a list (for each dimension) of indices, e.g. 
        [[i0, j0, k0], [i1, j1, k1], ...].
    
    RefIm : SimpleITK image
        A 3D image that occupies the grid that the Indices belong to.
        
        
    Ouputs:
    ******
    
    Points : list of a list of floats
        A list (for each point) of a list (for each dimension) of the 
        coordinates of Indices, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    Points = []
    
    for Index in Indices:
        #print(f'\n\nIndex = {Index}')
        
        #Index = [int(item) for item in Index]
        
        Point = Index2Point(Index, RefIm)
        
        Points.append(Point)
        
    return Points
    




def ContourData2Points(ContourData):
    """
    Re-format a flat list of coordinates to a list of a list of [x, y, z]
    coordinates.
    
    Inputs:
    ******
    
    ContourData : list of strings
        Flat list of [x, y, z] coordinates as strings.
    
        
    Outputs:
    *******
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points,
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
    """
    
    ContourData = [float(item) for item in ContourData]
    
    Points = []
            
    # Iterate for all points in threes:
    for p in range(0, len(ContourData), 3):
        point = [ContourData[p], ContourData[p+1], ContourData[p+2]]
        
        Points.append(point)
        
    return Points








def Points2ContourData(Points):
    """
    Re-format a list of points to a flat list as required for ContourData.
    
    Inputs:
    ******
    
    Points : list of a list of floats
        List (for each point) of a list (for each dimension) of points,
        e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        
        
    Outputs:
    *******
    
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
    ******
    
    IndsByContour : list of lists of lists of integers
        A list (for each contour) of a list (for all polygons) of a list (for 
        each dimension) of indices.
                     
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
        
    Outputs:
    *******
    
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








def GetLabmapImsByRoi(Rts, DicomDir, CStoSliceIndsByRoi, RefIm):
    """
    Get a list of 3D labelmap (SimpleITK) images - one per ROI.  
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    CStoSliceIndsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each Referenced SOP instance UID in the Contour Sequence.
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    *******
        
    LabmapIms : List of SimpleITK images
        A list of 3D zero-padded labelmap (SimpleITK) images - one per ROI.
    """
    
    from RtsTools import GetIndsByRoi
    
    LabmapIms = []
    
    IndsByRoi = GetIndsByRoi(Rts, DicomDir)
    
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
    ******
    
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
    *******
    
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
    ******
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
    
    FrameToSliceInds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    *******
        
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
    ******
    
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
    *******
        
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
    ******
    
    LabelmapIm : SimpleITK image
        A zero-padded version of PixelData.
        
    Non0FrameInds : List of integers (optional; None by default)  
        Use PerFrameFunctionalGroupsSequence-to-slice indices 
        (PFFGStoSliceInds) if they are known. If this input is not specified a 
        list of the frame numbers that are non-zero will be determined.
        
        
    Outputs:
    *******
        
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
        








def PixArr2IndsByFrame(PixArr, FrameToSliceInds, Thresh=0.5):
    """
    Convert a 3D pixel array to a list (for each frame) of a list of indices
    that define the equivalent contour for the mask in each frame.
    
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
 
    Inputs:
    ******
    
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
    *******
    
    IndsByObjByFrame : list of a list of a list of floats
        A list (for each frame) of a list (for each contour) of a list of 
        [x, y, z] coordinates of the polygons converted from each frame in
        PixArr.
    """
    
    import numpy as np
    from scipy.sparse import issparse
    from skimage.measure import find_contours
    
    IndsByObjByFrame = []
    
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
        IndsByObj_2D = find_contours(ExpandedMask.T, Thresh)
        
        #print(f'\n\n\nFrameNum = {FrameNum}')
        #print(f'\nlen(IndsByObj_2D) = {len(IndsByObj_2D)}')
        #print(f'\nIndsByObj_2D = {IndsByObj_2D}')
 
        # Remove one row and column to shift the indices back to their original
        # locations:
        IndsByObj_2D = [np.subtract(x, 1).tolist() for x in IndsByObj_2D]
        
        #print(f'\n\n\nFrameNum = {FrameNum}')
        #print(f'\nlen(IndsByObj_2D) = {len(IndsByObj_2D)}')
        #print(f'\nIndsByObj_2D = {IndsByObj_2D}')
        
        SliceNum = float(FrameToSliceInds[FrameNum])
            
        # Add the index for the 3rd dimension:
        IndsByObj_3D = []
        
        for IndsThisObj_2D in IndsByObj_2D:
            IndsThisObj_3D = [Ind_2D + [SliceNum] for Ind_2D in IndsThisObj_2D]
            
            IndsByObj_3D.append(IndsThisObj_3D)
        
        #print(f'\nlen(IndsByObj_3D) = {len(IndsByObj_3D)}')
        #print(f'\nIndsByObj_3D = {IndsByObj_3D}')
        
        IndsByObjByFrame.append(IndsByObj_3D)
        
        
    #print(f'\nlen(IndsByObjByFrame) = {len(IndsByObjByFrame)}')
    #print(f'\nIndsByObjByFrame = {IndsByObjByFrame}')
 
    return IndsByObjByFrame









def PixArr2PtsByContour(PixArr, FrameToSliceInds, DicomDir, Thresh=0.5):
    """
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates.
 
    Inputs:
    ******
    
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
    IndsByObjByFrame = PixArr2IndsByFrame(PixArr, FrameToSliceInds, Thresh=0.5)
    
    #Size, Spacings, ST,\
    #IPPs, Dirs = GetImageAttributes(DicomDir, Package='pydicom')
    
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









def NumOfListsAtDepthTwo(NestedList):
    """
    Return the number of lists in a nested list at depth two.
    
    Example uses:  Determining the number of objects/contours in 
    PtsByObjByFrame or CntDataByObjByFrame.  
    
    N.B.:
        
    PtsByObjByFrame is a list (for each frame) of a list (for each object,
    i.e. contour) of a list (for each dimension) of the physical coordinates 
    that define the mask in each frame in a pixel array.
    
    CntDataByObjByFrame is a list (for each frame) of a list (for each object, 
    i.e. contour) of a flattened list of [x, y, z] physical coordinates that 
    define the mask in each frame in PixArr as strings.
        
    So the number of contours is the total sum of all objects in the nested 
    list in the first two levels.
    
    
    Inputs:
    ******
    
    NestedList : list of a list of a list of items
        Example: A list (for each frame) of a list (for each object, i.e. 
        contour) of a list (for each dimension) of the physical coordinates 
        that define the mask in each frame in a pixel array.
 
    
    Outputs:
    *******
    
    N : integer
        The total number of lists in a nested list up to two levels down.
    """
    
    N = 0
    
    for i in range(len(NestedList)):
        N = N + len(NestedList[i])
        
    return N




            
    

def Indices2Mask(Inds, RefImage):
    """
    Convert a list of indices to a (filled) mask (pixel array).
       
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    Inputs:
    ******
    
    Inds : list of a list of integers
        A list (for all points) of a list (for each dimension) of indices.
                     
    RefImage : SimpleITK image
        Reference image whose size will be used to define the shape of the 
        pixel array.
        
        
    Outputs:
    *******
    
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
        
    
    C, R, S = RefImage.GetSize()
    
    Mask = np.zeros((1, R, C))
    
    #print(f'\nMask.shape = {Mask.shape}')
    
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
        
        if 0 <= y < R and 0 <= x < C:
            #print(f'\n   Mask.shape = {Mask.shape}')
            
            Mask[0, y, x] = 1
    
        
    return Mask






def ConvertImagePixelType(Image, NewPixelType):
    """
    Convert a 3D SimpleITK image to a new pixel type.
    
    Inputs:
    ******                      
        
    Image : SimpleITK image
        The 3D image whose pixel type is to be converted.
        
    NewPixelType : string
        The pixel type to convert Image to.  Acceptable inputs are a subset of
        the SimpleITK pixel types 
        (https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html):
            
        - 'sitkUInt32' 	Unsigned 32 bit integer
        - 'sitkFloat32' 	32 bit float

    
    
    Outputs:
    *******
    
    Image : SimpleITK image
        The 3D image with modified pixel type.

    
    Note:
    ----
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """
    
    import SimpleITK as sitk
    
    if NewPixelType == 'sitkUInt32':
        OutputPixelType = sitk.sitkUInt32
        
    elif NewPixelType == 'sitkFloat32':
        OutputPixelType = sitk.sitkFloat32
    
    else:
        msg = '"NewPixelType" must be either "sitkUInt32" or "sitkFloat32".'
        
        raise Exception(msg)
        
    
    Caster = sitk.CastImageFilter()
    
    Caster.SetOutputPixelType(OutputPixelType)
    
    Image = Caster.Execute(Image)
    
    return Image
    