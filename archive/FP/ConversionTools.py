# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:38:31 2020

@author: ctorti

    
"""



"""
******************************************************************************
******************************************************************************
CONVERSION FUNCTIONS
******************************************************************************
******************************************************************************
"""       


def Point2Index(Point, RefIm, Rounding=True):
    """
    07/07/21:
        See pt_to_ind in conversion_tools.inds_pts.py
        
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
    07/07/21:
        See pts_to_inds in conversion_tools.inds_pts.py
        
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
    07/07/21:
        See ind_to_pt in conversion_tools.inds_pts.py
        
    Convert an index (in the Image Coordinate System) to a physical point (in 
    the Patient Coordinate System).
    
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
    
    import numpy
    from GeneralTools import GetListOfDataTypes
    
    Dtypes = GetListOfDataTypes(Index)
    
    if len(Dtypes) == 1:
        Dtype = Dtypes[0]
        
        """All items in Index are of a common data type."""
        if Dtype == int:
            Point = RefIm.TransformIndexToPhysicalPoint(Index)
        
        elif Dtype == float:
            Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
        
        elif Dtype in [numpy.ndarray, numpy.float64, numpy.int32]:
            """ Convert from numpy array to list: """
            Index = Index.tolist()
            
            Dtypes = GetListOfDataTypes(Index)
            
            if len(Dtypes) == 1:
                Dtype = Dtypes[0]
                
                """All items in Index are of a common data type."""
                if Dtype == int:
                    Point = RefIm.TransformIndexToPhysicalPoint(Index)
                
                elif Dtype == float:
                    Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
                
                else:
                    msg = f"The data type of Index is {Dtype}. It must be "\
                          + "'int', 'float', 'numpy.ndarray' or 'numpy.float64'."
            
                    raise Exception(msg)
                    
            elif len(Dtypes) == 2 and int in Dtypes and float in Dtypes:
                """Allow for the possibility of a mixture of integers and 
                floats."""
                
                """ Convert integers to float: """
                Index = [float(item) for item in Index]
                
                Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
                
            else:
                msg = f"The data type of Index is {Dtypes}. It must be 'int',"\
                      + " 'float' or 'numpy.ndarray'."
            
        else:
            msg = f"The data type of Index is {Dtypes}. It must be 'int',"\
                  + " 'float', 'numpy.ndarray' or 'numpy.float64'."
            
            raise Exception(msg)
            
    elif len(Dtypes) == 2 and int in Dtypes and float in Dtypes:
        """Allow for the possibility of a mixture of integers and floats."""
        
        """ Convert integers to float: """
        Index = [float(item) for item in Index]
        
        Point = RefIm.TransformContinuousIndexToPhysicalPoint(Index)
        
    else:
        msg = f"The data type of Index is {Dtypes}. It must be either "\
              + "'int' or 'float'."
              
        raise Exception(msg)
         
        
    #return Point # 29/01/2021
    return list(Point) # 29/01/2021


def Indices2Points(Indices, RefIm):
    """
    07/07/21:
        See inds_to_pts in conversion_tools.inds_pts.py
        
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
    07/07/21:
        See cntdata_to_pts in conversion_tools.inds_pts.py
        
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
    07/07/21:
        See pts_to_cntdata in conversion_tools.inds_pts.py
       
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


def PtsByCntByRoi2CntDataByCntByRoi(PtsByCntByRoi):
    """
    07/07/21:
        See ptsByCntByRoi_to_cntdataByCntByRoi in 
        conversion_tools.inds_pts.py
        
    Convert a list of points-by-contour-by-ROI to a list of 
    ContourData-by-contour-by-ROI, where ContourData is flat list of points in
    string format (as required for ContourData tag).
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
        
    Outputs:
    *******
    
    CntDataByCntByRoi : list of a list of a list of strings
        List (for each ROI) of a list (for all contours) of a flat list of 
        coordinates in PtsByCntByRoi converted from floats to strings.
    """
    
    CntDataByCntByRoi = []
        
    for r in range(len(PtsByCntByRoi)):
        CntDataByCnt = []
        
        if PtsByCntByRoi[r]:
            for c in range(len(PtsByCntByRoi[r])):
                pts = PtsByCntByRoi[r][c]
                
                CntDataByCnt.append(Points2ContourData(pts))
        
        CntDataByCntByRoi.append(CntDataByCnt)
        
    return CntDataByCntByRoi



def IndsByContour2PixArr(IndsByCnt, RefImage):
    """
    07/07/21:
        See indsByCnt_to_pixarr in conversion_tools.inds_pts_pixarrs.py
        
    Convert a list of indices grouped by contour to a pixel array.
       
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
    
    Inputs:
    ******
    
    IndsByCnt : list of lists of lists of integers
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
    
    Ncontours = len(IndsByCnt)
    
    PixArr = np.zeros((Ncontours, ImSize[1], ImSize[0]))
    
    # Note that IndsByContour are a list of indices for each contour:
    for c in range(Ncontours):
        Indices = IndsByCnt[c]
        
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


def PtsByCnt2PixArr(PtsByCnt, RefIm, LogToConsole=False):
    """
    07/07/21:
        See ptsByCnt_to_pixarr in conversion_tools.inds_pts_pixarrs.py
      
    Convert a list of points by contour to a pixel array.
    
    Inputs:
    ******
    
    PtsByCnt : list of a list of a list of floats
        List (for all contours) of a list (for each point) of a list (for each 
        dimension) of coordinates.
    
    RefIm : SimpleITK image
        The 3D image that PtsByCnt relate to.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArr : Numpy data arrays
        A pixel array containing as many frames as the number of contours in
        PtsByCnt.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of PtsByCnt2PixArr():')
        print('\n\n', '-'*120)
        
    IndsByCnt = []
    
    if PtsByCnt == []:
        import numpy as np
        
        C, R, F = RefIm.GetSize()
        
        #PixArr = np.zeros((F, R, C), dtype='uint')
        PixArr = np.zeros((F, R, C), dtype=np.uint8)
        
    else:
        for Pts in PtsByCnt:
            Inds = Points2Indices(Pts, RefIm)
                
            IndsByCnt.append(Inds)
    
        """ Convert IndsByCnt to a pixel array. """
        PixArr = IndsByContour2PixArr(IndsByCnt, RefIm)
    
    
    if LogToConsole:
        print(f'   PixArr.shape = {PixArr.shape}')
        print('-'*120)
    
    return PixArr


def PixArr2ListOfInds(PixArr, F2Sinds):
    """
    07/07/21:
        See pixarr_to_listOfInds in conversion_tools.inds_pts_pixarrs.py
      
    Convert a binary pixel array to a list (for eachof a list of indices of non-zero 
    elements.  If PixArr is 2D, the function will return a list of length 1 of 
    a list of indices.  If PixArr is 3D, the function will return a list of
    length N of a list of indices, where N is the size of PixArr along the 
    third dimension. 
    
    Inputs:
    ******
    
    PixArr : Numpy array
        A binary pixel array representing a segmentation (2D) or segment (3D).
    
    F2Sinds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number.
    
    Outputs:
    *******
    
    ListOfInds : list of a list of integers
        A list (of length N, where N = number of frames in PixArr) of a list
        (for each dimension) of a list of integer coordinates for each non-zero
        element in each 2D frame in PixArr.
    """
    
    import numpy as np
    
    dim = PixArr.ndim
    
    if dim == 2:
        R, C = PixArr.shape
        F = 1
    elif dim == 3:
        F, R, C = PixArr.shape
    else:
        msg = "'PixArr' must have 2 or 3 dimensions."
        raise Exception(msg)
    
    ListOfInds = []
    
    for f in range(F):       
        yInds, xInds = np.where(PixArr[f] > 0)
        
        xInds = xInds.tolist()
        yInds = yInds.tolist()
        
        if dim == 2: 
            Inds = [[yInds[i], xInds[i]] for i in range(len(xInds))]
        else:
            s = F2Sinds[f]
            
            Inds = [[yInds[i], xInds[i], s] for i in range(len(xInds))]
        
        ListOfInds.append(Inds)
        
    return ListOfInds


def PtsByCntByRoi2PixArrByRoi(PtsByCntByRoi, RefIm, LogToConsole=False):
    """
    07/07/21:
        See ptsByCntByRoi_to_pixarrByRoi in conversion_tools.inds_pts_pixarrs.py
     
    Convert a list of points in all contours in all ROIs to a list of pixel
    arrays grouped by ROI.
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
    
    RefIm : SimpleITK image
        The 3D image that relates to PtsByCntByRoi.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
                       
                            
    Outputs:
    *******
    
    PixArrByRoi : list of Numpy data arrays
        A list (for each ROI) of 3D pixel arrays. Each pixel array contains as
        many frames as the number of contours in the list of points-by-contour
        for that given ROI.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of PtsByCntByRoi2PixArrByRoi():')
        print('\n\n', '-'*120)
        
    PixArrByRoi = []
    
    R = len(PtsByCntByRoi)
    
    for r in range(R):        
        """ Convert to a pixel array. """
        PixArr = PtsByCnt2PixArr(PtsByCntByRoi[r], RefIm, LogToConsole)
        
        PixArrByRoi.append(PixArr)
        
    if LogToConsole:
        shapes = [PixArrByRoi[r].shape for r in range(R)]
        print(f'   Shape of each PixArrByRoi = {shapes}')
        print('-'*120)
    
    return PixArrByRoi


def GetLabmapArrByRoi(Rts, DicomDir, C2SindsByRoi, RefIm):
    """
    07/07/21:
        See get_labarrByRoi in conversion_tools.inds_pts_pixarrs.py
     
    Get a list of 3D labelmap Numpy arrays - one per ROI.  
    
    Inputs:
    ******
    
    Rts : Pydicom object
        RTS object.
    
    DicomDir : string
        Directory containing the corresponding DICOMs.
    
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers  
        that corresponds to each contour.
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    *******
        
    LabmapArrByRoi : List of Numpy arrays
        A list of 3D zero-padded labelmap Numpy arrays - one per contour.
    """
    
    from RtsTools import GetIndsByRoi
    
    LabmapArrByRoi = []
    
    IndsByRoi = GetIndsByRoi(Rts, DicomDir)
    
    for RoiNum in range(len(IndsByRoi)):
        PixArr = IndsByContour2PixArr(IndsByContour=IndsByRoi[RoiNum], 
                                      RefImage=RefIm)
    
        LabmapArr = PixArr2LabmapArr(PixArr=PixArr, 
                                     NumOfSlices=RefIm.GetSize()[2], 
                                     FrameToSliceInds=C2SindsByRoi[RoiNum])
        
        LabmapArrByRoi.append(LabmapArr)
    
    return LabmapArrByRoi


def PixArr2LabmapArr(PixArr, NumOfSlices, F2Sinds):
    """
    07/07/21:
        See pixarr_to_labarr in conversion_tools.inds_pts_pixarrs.py
     
    Convert a 3D SEG pixel array to a 3D SEG labelmap array.  The array will 
    contain the 2D pixel arrays at the same frame position as the slice
    position of the corresponding DICOM series (with 2D zero masks in-between).
    
    
    
    Inputs:
    ******
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
        
    NumOfSlices : integer
        The total number of frames in the labelmap (e.g. the number of slices 
        in a image series).
        
    F2Sinds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number.
        
    Outputs:
    *******
    
    Labmap : Numpy array
        Zero-padded version of PixArr.
    
    
    Notes:
    *****
    
    'Labmap array' refers to a 3D pixel array that contains the non-empty 
    frames in PixArr, padded with empty frames so that the size of the labelmap
    array matches the size of the 3D DICOM image. This is not to be confused 
    with 'Labmap image', which is the SimpleITK image representation of a 
    labmap array.
    """
    
    import numpy as np
    
    NumOfFrames, R, C = PixArr.shape
    
    """ Initialise LabmapArr: """
    #Labmap = np.zeros((NumOfSlices, R, C), dtype='bool')
    #Labmap = np.zeros((NumOfSlices, R, C), dtype='uint')
    LabmapArr = np.zeros((NumOfSlices, R, C), dtype=PixArr.dtype)
    
    """ Note: NumOfFrames is the number of non-zero frames that make up PixArr,
    whereas NumOfSlices are the total frame count in 3D DICOM image. """
    
    #print(f'\nPixArr.shape = {PixArr.shape}')
    ##print(f'\n\nThe maximum value in PixArr is {PixArr.max()}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    #print(f'LabmapArr.shape = {LabmapArr.shape}')
    
    if NumOfFrames:
        """ Over-write all frames in LabmapArr with non-zero 2D masks from 
        PixArr.  """
        for i in range(len(F2Sinds)):
            """ The slice number for the i^th frame """
            s = F2Sinds[i]
            
            #print(f'Frame {i} corresponds to slice {s}')
            
            #LabmapArr[s] = PixArr[i] # 13/02 This was over-writting frames for 
            # pixel arrays containing more than one frame for the same slice!!!
            
            LabmapArr[s] += PixArr[i]
        
    
    #print(f'\n\nThe maximum value in LabmapArr is {LabmapArr.max()}')
    
    return LabmapArr





def PixArr2Image(PixArr, F2Sinds, RefIm):
    """
    07/07/21:
        See pixarr_to_im in conversion_tools.pixarrs_ims.py
        
    Convert a 3D pixel array to a 3D SimpleITK image.  
    
    Inputs:
    ******
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
    
    F2Sinds : List of integers
        The DICOM slice numbers that correspond to each frame in PixArr,
        e.g. PFFGStoSliceInds = PerFrameFunctionalGroupsSequence-to-slice
        number
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to.
    
        
    Outputs:
    *******
        
    LabIm : SimpleITK image
        A SimpleITK image containing zero-padded indices from PixArr.
    """
    
    import SimpleITK as sitk
    
    ImSize = RefIm.GetSize()
    
    #print(f'\nImSize = {ImSize}')
    #print(f'PixArr.shape = {PixArr.shape}')
    #print(f'ImSize[2] = {ImSize[2]}')
    #print(f'FrameToSliceInds = {FrameToSliceInds}')
    
    LabmapArr = PixArr2LabmapArr(PixArr=PixArr, NumOfSlices=ImSize[2], 
                                 F2Sinds=F2Sinds)
    
    LabIm = sitk.GetImageFromArray(LabmapArr)
    
    # Set the Origin, Spacing and Direction of LabIm to that of RefIm:
    LabIm.SetOrigin(RefIm.GetOrigin())
    LabIm.SetSpacing(RefIm.GetSpacing())
    LabIm.SetDirection(RefIm.GetDirection())
    
    #print(f'\n\nThe maximum value in Labmap is {Labmap.max()}')
    #print(f'\n\nThe maximum value in LabIm is {ImageMax(LabIm)}')
        
    return LabIm


def Image2BinaryOnesImage(Image):
    """
    07/07/21:
        See im_to_binary_ones_im in conversion_tools.pixarrs_ims.py
        
    Create a 3D binary ones SimpleITK image from a non-binary 3D SimpleITK image.  
    
    Inputs:
    ******
    
    Image : SimpleITK image
    
        
    Outputs:
    *******
        
    MaskIm : SimpleITK image
        A 3D ones SimpleITK image with the same image attributes as Image.
    
    
    Notes:
    *****
    
    This function is similar to PixArr2Image.
    """
    
    import numpy as np
    import SimpleITK as sitk
    
    C, R, F = Image.GetSize()
    
    PixArr = np.ones((F, R, C), dtype=np.uint8)
    
    MaskIm = sitk.GetImageFromArray(PixArr)
    
    # Set the Origin, Spacing and Direction of LabIm to that of Image:
    MaskIm.SetOrigin(Image.GetOrigin())
    MaskIm.SetSpacing(Image.GetSpacing())
    MaskIm.SetDirection(Image.GetDirection())
    
    return MaskIm


def PixArr2ImagesBySeg(PixArr, F2SindsBySeg, FrameNumsBySeg, RefIm):
    """
    07/07/21:
        See pixarr_to_imsBySeg in conversion_tools.pixarrs_ims.py
    
    Convert a 3D pixel array to a list of 3D SimpleITK images - one per
    segment.  
    
    Inputs:
    ******
    
    PixArr : Numpy array
        SEG PixelData as a Numpy array.
    
    F2SindsBySeg : List of list of integers
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
        
    LabIms : List of SimpleITK images
        A list of 3D zero-padded labelaps - one per segment.
    """
    
    LabIms = []
    
    #print(f'\nF2SindsBySeg = {F2SindsBySeg}\n')
    #print(f'FrameNumsBySeg = {FrameNumsBySeg}')
    
    for i in range(len(FrameNumsBySeg)):
        # The starting (= a) and ending (= b) frame numbers:
        a = FrameNumsBySeg[i][0]
        b = FrameNumsBySeg[i][-1]
        
        
        #print(f'a = {a}, b = {b}')
        #print(f'PixArr.shape = {PixArr.shape}')
        #print(f'PixArr[{a}:{b+1}].shape = {PixArr[a:b+1].shape}')
        
        LabIm = PixArr2Image(PixArr=PixArr[a:b+1],
                             RefIm=RefIm,
                             F2Sinds=F2SindsBySeg[i])
        
        LabIms.append(LabIm)
        
    return LabIms


def Image2PixArr(Image, Non0FrameInds=None):
    """
    06/07/21:
        See conversion_tools.im_pixarr.py
        
    Convert a 3D SimpleITK image to a 3D SEG pixel array.  
    
    Inputs:
    ******
    
    Image : SimpleITK image
        A image (e.g. labelmap).
        
    Non0FrameInds : List of integers (optional; None by default)  
        Use PerFrameFunctionalGroupsSequence-to-slice indices 
        (PFFGStoSliceInds) if they are known. If this input is not specified a 
        list of the frame numbers that are non-zero will be determined.
        
        
    Outputs:
    *******
        
    PixArr : Numpy array
        The non-zero frames in Image as a Numpy array.
                           
    Non0FrameInds : List of integers
        List of the indices of non-zero frames in Image.  This should be
        equivalent to the PerFrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds).
    """
    
    import SimpleITK as sitk
    import numpy as np
    
    """ Convert from SimpleITK image to a Numpy data array: """
    Nda = sitk.GetArrayFromImage(Image)
    
    if Non0FrameInds == None:
        Non0FrameInds = []
        
        #print('')
        
        """ Loop through all frames in Nda and get the sum of all pixels: """
        for i in range(Nda.shape[0]):
            Sum = np.amax(Nda[i])
            
            #print(f'i = {i}, Sum = {Sum}')
            
            if Sum:
                """ This is a non-zero frame: """
                Non0FrameInds.append(i)
    
    """ Initialise PixArr: """
    PixArr = np.zeros((len(Non0FrameInds), Nda.shape[1], Nda.shape[2]))
        
    """ Use Non0FrameInds to index the non-zero frames in Nda: """
    for i in range(len(Non0FrameInds)):
        PixArr[i] = Nda[Non0FrameInds[i]]
        
    
    return PixArr, Non0FrameInds


def ImageBySeg2PixArrBySeg_OLD(ImageBySeg, Non0FrameIndsBySeg=None): # 08/06/21
    """
    Convert a list of 3D SimpleITK images to a list of 3D SEG pixel arrays.  
    
    Inputs:
    ******
    
    ImageBySeg : List of SimpleITK images
        A list (e.g. for each segment) of images (e.g. labelmap).
        
    Non0FrameIndsBySeg : List of a list of integers (optional; None by default)  
        A list (for each segment) of a list of the indices of all non-zero
        frames in each image. If this input is not specified it will be
        determined.
        
        
    Outputs:
    *******
        
    PixArrBySeg : List of Numpy arrays
        A list (for each segment) of the non-zero frames in each image in
        ImagebySeg.
                           
    Non0FrameIndsBySeg : List of a list of integers
        List (for each segment) of a list (for each frame) of the indices of 
        non-zero frames in each image in ImageBySeg.
    """
    
    #import SimpleITK as sitk
    #import numpy as np
    
    PixArrBySeg = []
    #Non0FrameIndsBySeg = [] # 08/06/21
    
    if Non0FrameIndsBySeg == None:
        Non0FrameIndsBySeg_out = []
    
    for i in range(len(ImageBySeg)):
        #if Non0FrameInds == None: # 08/06/21
        if Non0FrameIndsBySeg == None: # 08/06/21
            PixArr, Non0FrameInds = Image2PixArr(ImageBySeg[i])
        else:
            #PixArr, Non0FrameInds = Image2PixArr(ImageBySeg[i], Non0FrameInds) # 08/06/21
            PixArr, Non0FrameInds = Image2PixArr(ImageBySeg[i], 
                                                 Non0FrameIndsBySeg[i]) # 08/06/21
            
        PixArrBySeg.append(PixArr)
        
        #Non0FrameIndsBySeg.append(Non0FrameInds) # 08/06/21
        Non0FrameIndsBySeg_out.append(Non0FrameInds) # 08/06/21
    
    
    #return PixArrBySeg, Non0FrameIndsBySeg # 08/06/21
    
    # 08/06/21:
    if Non0FrameIndsBySeg != None:
        Non0FrameIndsBySeg_out = Non0FrameIndsBySeg 
        
    return PixArrBySeg, Non0FrameIndsBySeg_out # 08/06/21


def ImageBySeg2PixArrBySeg(ImageBySeg):
    """
    07/07/21:
        See imBySeg_to_pixarrBySeg in conversion_tools.pixarrs_ims.py
    
    Convert a list of 3D SimpleITK images to a list of 3D SEG pixel arrays.  
    
    Inputs:
    ******
    
    ImageBySeg : List of SimpleITK images
        A list (e.g. for each segment) of images (e.g. labelmap).
        
        
    Outputs:
    *******
        
    PixArrBySeg : List of Numpy arrays
        A list (for each segment) of the non-zero frames in each image in
        ImagebySeg.
                           
    Non0FrameIndsBySeg : List of a list of integers
        List (for each segment) of a list (for each frame) of the indices of 
        non-zero frames in each image in ImageBySeg.
    """
    
    #import SimpleITK as sitk
    #import numpy as np
    
    PixArrBySeg = []
    Non0FrameIndsBySeg = []
    
    for i in range(len(ImageBySeg)):
        PixArr, Non0FrameInds = Image2PixArr(ImageBySeg[i])
            
        PixArrBySeg.append(PixArr)
        
        Non0FrameIndsBySeg.append(Non0FrameInds)
    
    
    return PixArrBySeg, Non0FrameIndsBySeg


def PixArr2IndsByFrame(PixArr, F2Sinds, Thresh=0.5, LogToConsole=False):
    """
    07/07/21:
        See pixarr_to_inds_by_frame in conversion_tools.inds_pts_pixarrs
        
    Convert a 3D pixel array to a list (for each frame) of a list of indices
    that define the equivalent contour for the mask in each frame.
    
    Adapted from https://programtalk.com/vs2/python/7636/sima/sima/ROI.py/
 
    Inputs:
    ******
    
    PixArr : Numpy array
        A FxRxC (frames x rows x cols) Numpy array containing F RxC masks in 
        PixArr.
    
    F2Sinds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
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
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of PixArr2IndsByFrame()')
        print('\n\n', '-'*120)
        
    IndsByObjByFrame = []
    
    for FrameNum, Frame in enumerate(PixArr):
        if issparse(Frame):
            Frame = np.array(Frame.astype('byte').todense())
 
        if (Frame != 0).sum() == 0:
            """ If Mask is empty, just skip it """
            continue
 
        """ Add an empty row and column at both sides of the frame to ensure 
        that any segmentation pixels near an edge are found. """
        ExpandedDims = (Frame.shape[0] + 2, Frame.shape[1] + 2)
        
        ExpandedFrame = np.zeros(ExpandedDims, dtype=float)
        
        ExpandedFrame[1:Frame.shape[0] + 1, 1:Frame.shape[1] + 1] = Frame
 
        """ Get the indices that define contours in ExpandedFrame. """
        IndsByObj_2D = find_contours(ExpandedFrame.T, Thresh)
        
        if LogToConsole:
            print(f'FrameNum = {FrameNum}')
            print(f'len(IndsByObj_2D) = {len(IndsByObj_2D)}')
            #print(f'IndsByObj_2D = {IndsByObj_2D}')
 
        """ Remove one row and column to shift the indices back to their 
        original locations. """
        IndsByObj_2D = [np.subtract(x, 1).tolist() for x in IndsByObj_2D]
        
        if LogToConsole:
            #print(f'\n\n\nFrameNum = {FrameNum}')
            print(f'len(IndsByObj_2D) = {len(IndsByObj_2D)}')
            #print(f'IndsByObj_2D = {IndsByObj_2D}')
        
        SliceNum = float(F2Sinds[FrameNum])
            
        """ Add the index for the 3rd dimension. """
        IndsByObj_3D = []
        
        for IndsThisObj_2D in IndsByObj_2D:
            IndsThisObj_3D = [Ind_2D + [SliceNum] for Ind_2D in IndsThisObj_2D]
            
            IndsByObj_3D.append(IndsThisObj_3D)
        
        if LogToConsole:
            print(f'len(IndsByObj_3D) = {len(IndsByObj_3D)}')
            #print(f'\nIndsByObj_3D = {IndsByObj_3D}')
        
        IndsByObjByFrame.append(IndsByObj_3D)
        
    
    if LogToConsole:
        print(f'\nlen(IndsByObjByFrame) = {len(IndsByObjByFrame)}')
        #print(f'\nIndsByObjByFrame = {IndsByObjByFrame}')
 
    return IndsByObjByFrame



def PixArrByRoi2LabImByRoi(PixArrByRoi, F2SindsByRoi, RefIm, 
                           LogToConsole=False):
    """
    07/07/21:
        See pixarrByRoi_to_imByRoi in conversion_tools.pixarrs_ims.py
    
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates for each
    pixel array in a list of pixel-arrays-by-ROI.
 
    Inputs:
    ******
    
    PixArrByRoi : list of Numpy arrays
        A list (for each ROI) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    
    F2SindsByRoi : List of list of integers
        A list (for each ROI) of a list (for each frame) of the slice numbers 
        that correspond to each frame in each pixel array in PixArrByRoi.  This
        is equivalent to a list of Per-FrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds). 
    
    RefIm : SimpleITK image
        The 3D image that PixArr relates to (e.g. SrcIm).
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    LabImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images representing the pixel
        arrays in PixArrByRoi.
    """
    
    from ImageTools import GetImageInfo
    #from GeneralTools import UniqueItems, AreListsEqualToWithinEpsilon
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Running of PixArrByRoi2LabImByRoi():')
        print('\n\n', '-'*120)
        #print(f'   len(PixArrByRoi) = {len(PixArrByRoi)}')
        
    LabImByRoi = []
    NewF2SindsByRoi = []

    for r in range(len(PixArrByRoi)):
        LabIm = PixArr2Image(PixArr=PixArrByRoi[r],
                             F2Sinds=F2SindsByRoi[r],
                             RefIm=RefIm)
        
        LabImByRoi.append(LabIm)
    
        if LogToConsole:
            print(f'\n   Image info for LabImByRoi[{r}]:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(LabIm, LogToConsole)
        
        NewF2SindsByRoi.append(F2Sinds)
    
    
    #""" Check that the number of unique F2Sinds in NewF2SindsByRoi matches the
    #number of unique F2Sinds in F2SindsByRoi (i.e. that no frames were lost). """
    #NumIn = []
    #NumOut = []
    #
    #for r in range(len(F2SindsByRoi)):
    #    NumIn.append(len(UniqueItems(F2SindsByRoi[r])))
    #    
    #    NumOut.append(len(UniqueItems(NewF2SindsByRoi[r])))
    #    
    #print('\n', NumIn, NumOut)
    #
    #if not AreListsEqualToWithinEpsilon(NumIn, NumOut, epsilon=0.9):
    #    print('\nWARNING:  The number of unique F2SindsByRoi used to generate',
    #          f'the labelmaps, {NumIn}, doesn\'t match that in the labelmap',
    #          f'images, {NumOut}.')
    
    if LogToConsole:
        print('-'*120)
        
    return LabImByRoi



def PixArr2PtsByCnt(PixArr, F2Sinds, DicomDir, Thresh=0.5, LogToConsole=False):
    """
    07/07/21:
        See pixarr_to_ptsByCnt in conversion_tools.inds_pts_pixarrs
        
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
 
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
        
    Outputs:
    *******
    
    PtsByCnt : list of a list of floats
        A list (for each contour) of a list (for each dimension) of the 
        physical coordinates that define each distinct object in each frame of
        PixArr.
        
    CntDataByCnt : list of a list of strings
        A list (for each contour) of a flattened list of [x, y, z] physical 
        coordinates in PtsByCnt.
        
    C2Sinds : list of integers
        List (for each contour) of slice numbers that correspond to each 
        contour in PtsByCnt.
    """
    
    #from ImageTools import GetImageAttributes
    from ImageTools import ImportImage
    
    if LogToConsole:
        from GeneralTools import PrintPixArrBySeg, PrintPtsByCnt
        print('\n\n', '-'*120)
        print(f'Running of PixArr2PtsByCnt():')
        print('\n\n', '-'*120)
        print(f'Prior to conversion: \nF2Sinds = {F2Sinds}')
        PrintPixArrBySeg(PixArr)
        
    RefIm = ImportImage(DicomDir)
    
    # Convert the pixel array to a list of indices-by-object-by-frame:
    """ A list of objects originating from the same frame will be treated as 
    distinct contours. """ 
    IndsByObjByFrame = PixArr2IndsByFrame(PixArr, F2Sinds, Thresh=0.5)
    
    if LogToConsole:
        F = len(IndsByObjByFrame)
        print(f'\nlen(IndsByObjByFrame) = {F} frames')
        for f in range(F):
            O = len(IndsByObjByFrame[f])
            print(f'   len(IndsByObjByFrame[{f}]) = {O} objects')
            for o in range(O):
                N = len(IndsByObjByFrame[f][o])
                print(f'      len(IndsByObjByFrame[{f}][{o}]) = {N} indices')
    
    
    PtsByCnt = []
    CntDataByCnt = []
    C2Sinds = []
    
    for f in range(len(IndsByObjByFrame)):
        for o in range(len(IndsByObjByFrame[f])):
            Pts = Indices2Points(Indices=IndsByObjByFrame[f][o], 
                                 RefIm=RefIm)
            
            PtsByCnt.append(Pts)
            
            CntDataByCnt.append(Points2ContourData(Points=Pts))
            
            C2Sinds.append(F2Sinds[f])
    
    
    #print(f'\n\n\nlen(PtsByCnt) = {len(PtsByCnt)}')
    #print(f'\nlen(CntDataByCnt) = {len(CntDataByCnt)}')
    #print(f'\nlen(C2Sinds) = {len(C2Sinds)}')
    #print(f'\nPtsByCnt = {PtsByCnt}')
    #print(f'\nCntDataByCnt = {CntDataByCnt}')
    #print(f'\nC2Sinds = {C2Sinds}')
    
    if LogToConsole:
        print(f'\nAfter conversion: \nC2Sinds = {C2Sinds}')
        PrintPtsByCnt(PtsByCnt)
        print('-'*120)
 
    return PtsByCnt, CntDataByCnt, C2Sinds



def PixArrByRoi2PtsByCntByRoi(PixArrByRoi, F2SindsByRoi, DicomDir, Thresh=0.5, 
                              LogToConsole=False):
    """
    07/07/21:
        See pixarrByRoi_to_ptsByCntByRoi in conversion_tools.inds_pts_pixarrs
        
    Convert a 3D pixel array to a list (for each frame/contour) of a list (for 
    each point) of a list (for each dimension) of physical coordinates for each
    pixel array in a list of pixel-arrays-by-ROI.
 
    Inputs:
    ******
    
    PixArrByRoi : list of Numpy arrays
        A list (for each ROI) of a FxRxC (frames x rows x cols) Numpy array 
        containing F RxC masks in each pixel array.
    
    F2SindsByRoi : List of list of integers
        A list (for each ROI) of a list (for each frame) of the slice numbers 
        that correspond to each frame in each pixel array in PixArrByRoi.  This
        is equivalent to a list of Per-FrameFunctionalGroupsSequence-to-slice 
        indices (PFFGStoSliceInds).
    
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
    Thresh : float (optional; 0.5 by default)
        Threshold value used to binarise the labels in PixArr.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
 
    
    Outputs:
    *******
    
    PtsByCntByRoi : list of a list of a list of floats
        A list (for each ROI) of a list (for each contour) of a list (for each 
        dimension) of the physical coordinates that define each pixel array in
        PixArrByRoi.
        
    CntDataByObjByFrame : list of a list of a list of strings
        A list (for each ROI) of a list (for each contour) of a flattened list 
        of [x, y, z] physical coordinates that define each pixel array in
        PixArrByRoi.
    
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour in PtsByCntByRoi.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of PixArrByRoi2PtsByCntByRoi():')
        print('\n\n', '-'*120)
        
    PtsByCntByRoi = []
    CntDataByCntByRoi = []
    C2SindsByRoi = []
    
    #print(f'\n\n\nPixArrByRoi = {PixArrByRoi}')
    #print(f'len(PixArrByRoi) = {len(PixArrByRoi)}')
    #for r in range(len(PixArrByRoi)):
    #    print(f'type(PixArrByRoi[{r}]) = {type(PixArrByRoi[r])}')
    #    print(f'PixArrByRoi[{r}].shape = {PixArrByRoi[r].shape}')
    
    for r in range(len(PixArrByRoi)):
        #if PixArrByRoi[r]:
        if PixArrByRoi[r].shape[0]:
            PtsByCnt, CntDataByCnt,\
            C2Sinds = PixArr2PtsByCnt(PixArrByRoi[r], F2SindsByRoi[r], 
                                      DicomDir, Thresh, LogToConsole)
            
            PtsByCntByRoi.append(PtsByCnt)
            CntDataByCntByRoi.append(CntDataByCnt)
            C2SindsByRoi.append(C2Sinds)
    
    
    #print(f'\n\n\nlen(PtsByCntByRoi) = {len(PtsByCntByRoi)}')
    #print(f'\nlen(CntDataByCntByRoi) = {len(CntDataByCntByRoi)}')
    #print(f'\nlen(C2SindsByRoi) = {len(C2SindsByRoi)}')
    #print(f'\nPtsByCntByRoi = {PtsByCntByRoi}')
    #print(f'\nCntDataByCntByRoi = {CntDataByCntByRoi}')
    #print(f'\nC2SindsByRoi = {C2SindsByRoi}')
    
    if LogToConsole:
        print('-'*120)
 
    return PtsByCntByRoi, CntDataByCntByRoi, C2SindsByRoi


def Indices2Mask(Inds, RefImage):
    """
    07/07/21:
        See inds_to_mask in conversion_tools.inds_pts_pixarrs
      
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



def BinaryIm2Labelmap(BinaryIm):
    import SimpleITK as sitk

    ImFilt = sitk.BinaryImageToLabelMapFilter()
    
    Labelmap = ImFilt.Execute(BinaryIm)
    
    return Labelmap



def Labelmap2BinaryIm(Labelmap):
    import SimpleITK as sitk

    ImFilt = sitk.LabelMapToBinaryImageFilter()
    
    BinaryIm = ImFilt.Execute(Labelmap)
    
    return BinaryIm