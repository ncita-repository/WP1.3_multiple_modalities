# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:57:05 2020

@author: ctorti
"""


# Import packages:
#from copy import deepcopy
#import numpy as np
#import matplotlib.pyplot as plt

##import importlib
##import DicomTools
##importlib.reload(DicomTools)

#from DicomTools import GetDicomSOPuids

#from SegTools import GetPFFGStoSliceInds



"""
******************************************************************************
******************************************************************************
GENERAL UTILITY FUNCTIONS
******************************************************************************
******************************************************************************
"""


def ItemsUniqueToWithin(List, epsilon=1e-5):
    
    UniqueList = []
    
    L = len(List)
    
    if L < 2:
        print(f'The list has {L} items (it must have at least 2 items).')
        
        return UniqueList
    
    else:
        UniqueList.append(List[0])
        
        for i in range(L - 1):
            if abs(List[i] - List[0]) > epsilon:
                UniqueList.append(List[i])
                
        return UniqueList
    






def UniqueItems(Items, IgnoreZero=False):
    """
    Get list of unique items in a list or Numpy data array.
    
    Inputs:
    ------
    
    Items : list or Numpy data array
    
    IgnoreZero : boolean (optional; False by default)
        If IgnoreZero=True only unique non-zero items will be returned.
        
    
    Outputs:
    -------
    
    UniqueItems : list or Numpy data array
        List or Numpy array of unique items.
    """
    
    import numpy as np
    from copy import deepcopy
    
    OrigItems = deepcopy(Items)
    
    # If Items is a list convert to a Numpy array:
    if isinstance(Items, list):
        #Items = np.ndarray.flatten(Items)
        Items = np.array(Items)
    
    if IgnoreZero:
        inds = np.nonzero(Items)
        
        Items = Items[inds]
        
    #print(f'\nItems.shape = {Items.shape}')
    
    #UniqueItems = list(set(Items.tolist()))
    UniqueItems = np.unique(Items)
    
    # If Items was a list convert back to a list:
    if isinstance(OrigItems, list):
        UniqueItems = list(UniqueItems) # 18/12/2020
        #UniqueItems = sorted(list(UniqueItems)) # 18/12/2020
        

    
    return UniqueItems
        
    
        
def UniqueItems_new(Items, IgnoreZero=False, MaintainOrder=False):
    """
    04/09/21:
        Copied from FP/GeneralTools.py
        
    Get list of unique items in a list or Numpy data array.
    
    Inputs:
    ******
    
    Items : list of integers/floats or a Numpy data array
    
    IgnoreZero : boolean (optional; False by default)
        If IgnoreZero=True only unique non-zero items will be returned.
        
    MaintainOrder : boolean (optional; False by default)
        If MaintainOrder=True the unique items will be returned in their
        original order.
        
    
    Outputs:
    *******
    
    UniqueItems : list or Numpy data array or None
        List or Numpy array of unique items, or None if Items is empty.
    """
    
    import numpy as np
    from copy import deepcopy
    
    
    if not isinstance(Items, list) and not isinstance(Items, np.ndarray):
        msg = f'The input "Items" is data type {type(Items)}.  Acceptable '\
              + 'data types are "list" and "numpy.ndarray".'
        
        raise Exception(msg)
    
    
    OrigItems = deepcopy(Items)
    
    if isinstance(Items, list):
        #Items = np.array(Items)
    
        ##Items = list(np.concatenate(Items).flat)
        #Items = np.concatenate(Items)
        
        Items = np.array(FlattenList(Items))
            

    if Items.size == 0:
        return None
    
    
    if IgnoreZero:
        inds = np.nonzero(Items)
        
        Items = [Items[ind] for ind in inds[0]]
    
    if MaintainOrder:
        UniqueItems = []
        
        for item in Items:
            if not item in UniqueItems:
                UniqueItems.append(item)
        
    else:
        #UniqueItems = list(set(Items.tolist()))
        UniqueItems = np.unique(Items)
    
    # If Items was a list convert back to a list:
    if isinstance(OrigItems, list):
        UniqueItems = list(UniqueItems) # 18/12/2020
        #UniqueItems = sorted(list(UniqueItems)) # 18/12/2020
        

    
    return UniqueItems
    
    
    

def AppendItemToListAndSort(OrigList, ItemToAppend):
    """
    Append an item to a list and return the sorted list. 
    
    Inputs:
    ------
    
    OrigList : list
        The list to append to. 
        
    ItemToAppend : string, integer or float
        The item to append to OrigList.
        
        
    Outputs:
    -------
    
    NewList : list
        OrigList with ItemToAppend appended to it and sorted. 
    """
    
    from copy import deepcopy
    
    NewList = deepcopy(OrigList)

    NewList.append(ItemToAppend)
    
    NewList.sort()
    
    return NewList









def Unpack(Items):
    """
    Unpack list of 2D or 3D lists (e.g. [x, y] or [x, y, z] coordinates, or 
    [i, j] or [i, j, k] indices) into separate lists for each dimension.
    
    Input:
        Items - (List of floats/indices) List of a list of items,
                 e.g. [[x0, y0, z0], [x1, y1, z1], ...]
        
    Returns:
        X     - (List of floats/indices) List of all x coordinates in Items,
                 e.g. [x0, x1, x2, ...]
        
        Y     - (List of floats/indices) List of all y coordinates in Items,
                 e.g. [y0, y1, y2, ...]
        
        Z     - (List of floats/indices) List of all z coordinates in Items,
                 e.g. [z0, z1, z2, ...]
    """
    
    import numpy as np
    
    nda = np.array(Items)
    
    if len(nda.shape) < 2:
        msg = "The input argument is a list of length 1. There is nothing to" \
              + " unpack."
        
        raise Exception(msg)
        
    dim = nda.shape[1]
    
    # Initialise unpacked lists:
    X = []
    Y = []
    Z = []
    
    if dim == 2:
        # Unpack tuple and store in X, Y:
        for x, y in Items:
            X.append(x)
            Y.append(y)
            
        return X, Y
    
    elif dim == 3:
        # Unpack tuple and store in X, Y, Z:
        for x, y, z in Items:
            X.append(x)
            Y.append(y)
            Z.append(z)
            
        return X, Y, Z
    
    else:
        msg = f"The input argument has dim = {dim}. Only dim = 2 or dim = 3 "\
              + "is allowed."
        
        raise Exception(msg)








def AreItemsEqualToWithinEpsilon(Item0, Item1, epsilon=1e-06):
    """
    Are two items equal to within epsilon?
    
    Default value of epsilon is 1e-06.
    """
    
    AbsDiff = abs(Item0 - Item1)
    
    if AbsDiff < epsilon:
        IsEqual = True
    else:
        IsEqual = False
        
    return IsEqual






def AreListsEqualToWithinEpsilon(List0, List1, epsilon=1e-06):
    """
    Are two lists equal to within epsilon?
    
    Default value of epsilon is 1e-06.
    """
    
    import numpy as np
    
    # Get the maximum value of their absolute differences:
    AbsDiffs = abs(np.array(List0) - np.array(List1))
    MaxAbsDiff = max(AbsDiffs)
    
    if MaxAbsDiff < epsilon:
        IsEqual = True
    else:
        IsEqual = False
        
    return IsEqual






def PrintTitle(Title):
    
    ul = '*' * len(Title)
    print('\n' + Title)
    print(ul + '\n')
    
    return





def GetExtremePoints(Image):
    """
    Get the vertices that define the 3D extent of a 3D image.
    
    Inputs:
    ------
    
    Image : SimpleITK image
        A 3D image whose physical extent is to be determined.
    
        
    Outputs:
    -------
    
    pts : list of floats
        A list of length 8 of floats for all vertices that define the 3D extent
        of Image.
    """
    pts = [Image.TransformIndexToPhysicalPoint((0, 0, 0)), 
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), 0, 0)),
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), Image.GetHeight(), 0)),
           Image.TransformIndexToPhysicalPoint((0, Image.GetHeight(), 0)),
           Image.TransformIndexToPhysicalPoint((0, 0, Image.GetDepth())), 
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), 0, Image.GetDepth())),
           Image.TransformIndexToPhysicalPoint((Image.GetWidth(), Image.GetHeight(), Image.GetDepth())),
           Image.TransformIndexToPhysicalPoint((0, Image.GetHeight(), Image.GetDepth()))]
    
    return pts





def IsPointInPolygon(point, vertices):
    """
    Determine if a 3D point is in a 3D polygon.
    
    Inputs:
    ------
    
    Point : list of floats
        A list of length 3 of floats representing a 3D point.
                     
    Vertices : list of floats
        A list of length 8 of floats for all vertices that define a 3D volume.
        
        
    Outputs:
    -------
    
    PtInPoly : boolean
        True/False if Point is/isn't in Polygon.
    """
    
    from shapely.geometry import Polygon, Point
    
    poly = Polygon(vertices)
    
    pt = Point(point)
    
    return poly.contains(pt)







def ProportionOfPixArrInExtent_OLD(PixArr, FrameToSliceInds, DicomDir):
    """
    Determine what proportion of voxels that make up a 3D pixel array reside
    within the 3D volume of a 3D image.
    
    Inputs:
    ------
    
    PixArr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
        
    Outputs:
    -------
    
    p : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArr2PtsByContour
    from ImageTools import ImportImage
    
    # Convert the pixel array to a list (for each frame in PixArr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    PtsByFrame, CntDataByFrame = PixArr2PtsByContour(PixArr, 
                                                     FrameToSliceInds, 
                                                     DicomDir)
    
    Image = ImportImage(DicomDir)
    
    Vertices = GetExtremePoints(Image)
    
    IsInside = []
    
    for PtsInFrame in PtsByFrame:
        for pt in PtsInFrame:
            IsInside.append(IsPointInPolygon(pt, Vertices))
            
    return sum(IsInside)/len(IsInside)






def ProportionOfPixArrInExtent(PixArr, FrameToSliceInds, DicomDir):
    """
    COMMENT 19/01/2021:
        I was previously mostly getting 100% proportions and the odd 0%.  Now
        most or all of the results are 0%.  This needs to be checked.
        
    
    Determine what proportion of voxels that make up a 3D pixel array reside
    within the 3D volume of a 3D image.
    
    Inputs:
    ------
    
    PixArr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
        
    Outputs:
    -------
    
    p : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArr2PtsByContour
    from ImageTools import ImportImage
    
    # Convert the pixel array to a list (for each frame in PixArr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    PtsByObjByFrame,\
    CntDataByObjByFrame = PixArr2PtsByContour(PixArr, 
                                              FrameToSliceInds, 
                                              DicomDir)
    
    print('\n\n\nFunction ProportionOfPixArrInExtent():')
    print(f'   FrameToSliceInds = {FrameToSliceInds}')
    print(f'   len(PtsByObjByFrame) = {len(PtsByObjByFrame)}')
    print(f'   len(CntDataByObjByFrame) = {len(CntDataByObjByFrame)}')
    for f in range(len(PtsByObjByFrame)):
        print(f'   Frame {f} has {len(PtsByObjByFrame[f])} objects')
        for o in range(len(PtsByObjByFrame[f])):
            print(f'      Object {o} has {len(PtsByObjByFrame[f][o])} points')
            #print(f'   len(CntDataByFrame[{i}]) = {len(CntDataByFrame[i])}')
    #print(f'   \nPtsByFrame = {PtsByFrame}')
    #print(f'   \nCntDataByFrame = {CntDataByFrame}')
    
    
    Image = ImportImage(DicomDir)
    
    Vertices = GetExtremePoints(Image)
    
    IsInside = []
    
    for f in range(len(PtsByObjByFrame)):
        for o in range(len(PtsByObjByFrame[f])):
            for pt in PtsByObjByFrame[f][o]:
                IsInside.append(IsPointInPolygon(pt, Vertices))
            
    return sum(IsInside)/len(IsInside)







def CropNonZerosIn2dMask(Orig2dMask):
    
    import numpy as np
    
    #print('\nShape of Orig2dMask =', Orig2dMask.shape)
    
    # The non-zero elements in the Orig2dMask:
    non0 = np.nonzero(Orig2dMask)
    
    # If there are non-zero elements:
    if non0[0].size:
    
        a = non0[0][0]
        b = non0[0][-1]
        
        c = non0[1][0]
        d = non0[1][-1]
    
        # The non-zero (cropped) array:
        Cropped2dMask = Orig2dMask[a:b, c:d]
        
        #print('\nShape of Cropped2dMask =', Cropped2dMask.shape)
        
        return Cropped2dMask
    
    else:
        #print('There were no non-zero elements in the 2D mask.')
        
        return non0[0]
    




def Compare2dMasksFrom3dMasks(OrigSegRoi, NewSegRoi, OrigDicomDir, NewDicomDir):
    """ Compare cropped masks of non-zero elements only. """ 
    
    from DicomTools import GetDicomSOPuids
    from SegTools import GetPFFGStoSliceInds
    import matplotlib.pyplot as plt
    
    # Get the DICOM SOP UIDs:
    OrigSOPuids = GetDicomSOPuids(DicomDir=OrigDicomDir)
    NewSOPuids = GetDicomSOPuids(DicomDir=NewDicomDir)
    
    # Get the Per-frameFunctionalGroupsSequence-to-slice indices:
    OrigPFFGStoSliceInds = GetPFFGStoSliceInds(OrigSegRoi, OrigSOPuids)
    NewPFFGStoSliceInds = GetPFFGStoSliceInds(NewSegRoi, NewSOPuids)
    
    # Combined indices from OrigPFFGStoSliceInds and NewPFFGStoSliceInds:
    AllInds = OrigPFFGStoSliceInds + NewPFFGStoSliceInds
    
    # Remove duplicates:
    AllInds = list(set(AllInds))
    
    
    # Get the 3D SEG masks:
    Orig3dMask = OrigSegRoi.pixel_array
    New3dMask = NewSegRoi.pixel_array
    
    OrigShape = Orig3dMask.shape
    NewShape = New3dMask.shape
    
    print(f'Segments exist in OrigSegRoi on slices {OrigPFFGStoSliceInds}')
    print(f'Shape of Orig3dMask = {OrigShape}')
    print(f'\nSegments exist in NewSegRoi on slices {NewPFFGStoSliceInds}')
    print(f'Shape of New3dMask = {NewShape}\n')
    
    # Initialise the 3D cropped SEG masks:
    Orig3dMaskCropped = []
    New3dMaskCropped = []
    
    for i in range(OrigShape[0]):
        cropped = CropNonZerosIn2dMask(Orig3dMask[i])
        
        Orig3dMaskCropped.append(cropped)
        
    for i in range(NewShape[0]):
        cropped = CropNonZerosIn2dMask(New3dMask[i])
        
        New3dMaskCropped.append(cropped)
    
    
    Nrows = len(AllInds)
    
    Ncols = 2
    
    n = 1 # initialised sub-plot number
    
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 6*Nrows))
    
    for i in range(len(AllInds)):
        SliceNum = AllInds[i]
        
        # Does slice SliceNum have a segment in OrigSegRoi or NewSegRoi?
        if SliceNum in OrigPFFGStoSliceInds:
            OrigFrameNum = OrigPFFGStoSliceInds.index(SliceNum)
        
            ax = plt.subplot(Nrows, Ncols, n, aspect='equal')
            ax.imshow(Orig3dMaskCropped[OrigFrameNum])
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'Orig slice {SliceNum}')
            
        n += 1 # increment sub-plot number
        
        if SliceNum in NewPFFGStoSliceInds:
            NewFrameNum = NewPFFGStoSliceInds.index(SliceNum)
    
            ax = plt.subplot(Nrows, Ncols, n, aspect='equal')
            ax.imshow(New3dMaskCropped[NewFrameNum])
            ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
            ax.set_title(f'New slice {SliceNum}')
            
        n += 1 # increment sub-plot number
    
    return









def GetPhysicalShiftBetweenSlices(Image0, SliceNum0, Image1, SliceNum1):
    """
    Get the physical shift between SliceNum0 in Image0 and SliceNum1 in Image1
    based on the origins of the two slices (i.e. Image0[0,0,SliceNum0] and
    Image1[0,0,SliceNum1]).
    
    Inputs:
    ------
    
    Image0 : SimpleITK image
    
    SliceNum0 : integer
        The slice number in Image0.
    
    Image1 : SimpleITK image
    
    SliceNum1 : integer
        The slice number in Image1.
    
    
    Outputs:
    -------
    
    mmShift : Numpy array
        Shift in mm between Image0[0,0,SliceNum0] and Image1[0,0,SliceNum1].
    """
    
    import numpy as np
    
    mmShift = np.array(Image1.TransformIndexToPhysicalPoint([0,0,SliceNum1])) \
            - np.array(Image0.TransformIndexToPhysicalPoint([0,0,SliceNum0]))
            
    return mmShift







def GetPixelShiftBetweenSlices(Image0, SliceNum0, Image1, SliceNum1, RefImage,
                               Fractional=False):
    """
    Get the shift in pixels between SliceNum0 in Image0 and SliceNum1 in Image1
    based on the origins of the two slices (i.e. Image0[0,0,SliceNum0] and
    Image1[0,0,SliceNum1]).
    
    Inputs:
    ------
    
    Image0 : SimpleITK image
    
    SliceNum0 : integer
        The slice number in Image0.
    
    Image1 : SimpleITK image
    
    SliceNum1 : integer
        The slice number in Image1.
        
    RefImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert the
        physical shift to pixels.
    
    Fractional : boolean (optional; False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    
    
    Outputs:
    -------
    
    PixShift : Numpy array
        Shift in pixels between Image0[0,0,SliceNum0] and Image1[0,0,SliceNum1].
    """
    
    import numpy as np
    
    mmShift = GetPhysicalShiftBetweenSlices(Image0, SliceNum0, 
                                            Image1, SliceNum1)
    
    PixSpacing = np.array(RefImage.GetSpacing())
    
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








def MeanPixArr(PixArr, binary=True):
    """
    Perform pixel-by-pixel mean of all frames in a pixel array. Output a 
    binary pixel array if binary=True (or undefined).
    
    Inputs:
    ------
    
    PixArr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
        
    binary : boolean (optional; True by default)
        If binary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by 0.5.
        
        
    Outputs:
    -------
    
    MeanPixArr : Numpy array
        Mean pixel array.
    """
    
    import numpy as np
    
    F, R, C = PixArr.shape
    
    # Initialise MeanPixArr:
    MeanPixArr = np.zeros((1, R, C), dtype='uint')
    
    
    result = np.mean(PixArr, axis=0)
            
    if binary:
        MeanPixArr[0] = (result >= 0.5) * result
    else:
        MeanPixArr[0] = result
        
    return MeanPixArr







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
        
        
        
        
def SplitString(String, SearchString, Separator='_'):
    """
    Split a string to the string of characters on either side of a search
    string bounded by a separator.  
    
    Inputs:
    *******
    
    String : string
        The input string to be separated.
        
    SearchString : string
        The string that defines the split point.
        
    Separator : string (optional; '_' by default)
        The string that separates words within a string.
        
        e.g. The string "Series-4-Slice-9" has separator '-'; The string
        "Series_4_Slice_9" has separator '_'.
        
        
    Outputs:
    *******
        
    StringList : list of strings
        A list of length two of strings (e.g. ['string0', 'string1']).
    """
    
    StringSplit = String.split(Separator)
    
    ind = StringSplit.index(SearchString)
    
    StringList = []
    
    if ind - 1 < 0:
        StringList.append('')
    else:
        StringList.append(StringSplit[ind-1])
        
    if ind + 1 >= len(StringSplit):
        StringList.append('')
    else:
        StringList.append(StringSplit[ind+1])
        
    return StringList





def SliceNoFromFpath(Fpath, SearchString, Separator='_'):
    """
    Parse a filepath for a slice number based on a search string and (optional)
    separator.  This will only work for specific filepath formats.
    
    Inputs:
    *******
    
    Fpath : string
        The filepath to be parsed.
        
    SearchString : string
        The string that identifies the slice number.
        
        e.g. The string "S4_slice9" would require the search string "slice";
        The string "S4_s9" would require "s".
        
    Separator : string (optional; '_' by default)
        The string that separates words within a string.
        
        e.g. The string "S4-s9" has separator '-'; The string "S4_s9" has 
        separator '_'.
        
        
    Outputs:
    *******
        
    SliceNo : int
        The slice number parsed from Fpath.
    """
    
    from os.path import split, splitext
    
    Fname = split(splitext(Fpath)[0])[1]
    
    IntString = Fname.split(SearchString)[1]
    
    try:
        SliceNo = int(IntString)
        
        return SliceNo
        
    except ValueError:
        msg = f'Cannot convert {IntString} into int. \nFpath = {Fpath}, \n' \
              + f'Fname = {Fname}, \nIntString = {IntString}.'
        
        print(msg)
        
        return None
    
def FlattenList(ListOfLists):
    """
    04/09/21:
        Copied from FP/GeneralTools.py
    """
    if len(ListOfLists) == 0:
        return ListOfLists
    if isinstance(ListOfLists[0], list):
        return FlattenList(ListOfLists[0]) + FlattenList(ListOfLists[1:])
    return ListOfLists[:1] + FlattenList(ListOfLists[1:])

def ReduceListOfStringFloatsTo16(OrigList):
    """ 
    04/09/21:
        Copied from FP/GeneralTools.py
    
    Decimal strings (VR DS in DICOM standard) must be limited to 16 
    characters. Reduce a list of string floats to the required limit.
    
    Inputs:
    ******
        
    OrigList : list of strings
        The original list of string floats.
    
    
    Outputs:
    *******
        
    NewList : list of strings
        The list of strings with characters limited to 16.
    """
    
    #NewList = [item[:16] for item in OrigList]
    
    """
    print(f'\n\n\nOrigList = {OrigList}')
    print(f'type(OrigList) = {type(OrigList)}')
    print(f'type(OrigList[0]) = {type(OrigList[0])}')
    print(f'OrigList[0] = {OrigList[0]}')
    print(f'len(OrigList[0]) = {len(OrigList[0])}')
    print('\n\n\n')
    """
    
    charLimit = 16
    
    newList = []
    
    for item in OrigList:
        if len(item) > charLimit:
            newList.append(item[:charLimit])
        else:
            newList.append(item)
    
    return newList