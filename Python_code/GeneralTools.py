# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:57:05 2020

@author: ctorti
"""


"""
******************************************************************************
******************************************************************************
GENERAL UTILITY FUNCTIONS
******************************************************************************
******************************************************************************
"""

def GetIndOfNearestSlice(RefIm, RefInd, Im, UseCentre=True):
    """ 
    Get the index of the slice in Im that is nearest in z-position to the
    RefInd^th slice in RefIm.
    
    Inputs:
    ******
    
    RefIm : SimpleITK image
    
    RefInd : integer
        The slice index for RefIm.
    
    Im : SimpleITK image
        The image whose slice index nearest to RefIm[RefInd] is to be 
        determined.
    
    UseCentre : boolean (optional; True by default)
        If True, the z-position of the central pixel in each slice will be 
        considered. If False, the origin (0, 0) of each slice will.
    
    
    Outputs:
    *******
    
    Ind : integer
        The slice index in Im whose z-position is nearest to RefIm[RefInd].
    
    MinDiff : float
        The difference in z-position between RefIm[RefInd] and Im[Ind] in mm.
    """
    
    R_r, C_r, S_r = RefIm.GetSize()
    R_i, C_i, S_i = Im.GetSize()
        
    if UseCentre:
        i_r = R_r//2
        j_r = C_r//2
        
        i_i = R_i//2
        j_i = C_i//2
    else:
        i_r = 0
        j_r = 0
        
        i_i = 0
        j_i = 0
    
    """ The z-position of RefIm at RefInd: """
    RefZ = RefIm.TransformIndexToPhysicalPoint([i_r,j_r,RefInd])[2]
    
    """ All z-positions in Im: """
    Zs = [Im.TransformIndexToPhysicalPoint([i_i,j_i,k])[2] for k in range(S_i)]
    
    Diffs = [abs(RefZ - Zs[k]) for k in range(S_i)]
    
    MinDiff = min(Diffs)
    
    Ind = Diffs.index(MinDiff)
    
    return Ind, MinDiff





    
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
    






def FlattenList(ListOfLists):
    if len(ListOfLists) == 0:
        return ListOfLists
    if isinstance(ListOfLists[0], list):
        return FlattenList(ListOfLists[0]) + FlattenList(ListOfLists[1:])
    return ListOfLists[:1] + FlattenList(ListOfLists[1:])







def UniqueItems(Items, IgnoreZero=False, MaintainOrder=False):
    """
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
    




def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def xor(a, b):
    # use the symmetric difference operator (i.e. XOR operator)
    return list(set(a) ^ set(b))





def AppendItemToListAndSort(OrigList, ItemToAppend):
    """
    Append an item to a list and return the sorted list. 
    
    Inputs:
    ******
    
    OrigList : list
        The list to append to. 
        
    ItemToAppend : string, integer or float
        The item to append to OrigList.
        
        
    Outputs:
    *******
    
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
    
    Inputs:
    ******
    
    Items : list of floats/indices
        List of a list of items, e.g. [[x0, y0, z0], [x1, y1, z1], ...].
     
    
    Outputs:
    *******
    X : list of floats/indices
        List of all x coordinates in Items, e.g. [x0, x1, x2, ...].
    
    Y : list of floats/indices
        List of all y coordinates in Items, e.g. [y0, y1, y2, ...].
    
    Z : list of floats/indices
        List of all z coordinates in Items, e.g. [z0, z1, z2, ...].
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





def NumOfListsAtDepthTwo(NestedList):
    """
    Return the number of lists in a nested list at depth two.
    
    Inputs:
    ******
    
    NestedList : list of a list of a list of items
        A nested list of items.
 
    
    Outputs:
    *******
    
    N : integer
        The total number of lists in a nested list up to two levels down.
    """
    
    N = 0
    
    for i in range(len(NestedList)):
        N = N + len(NestedList[i])
        
    return N






def NumOfItemsInAListByDepth(List):
    """
    Return the number of items in a list at each depth.
    
    Inputs:
    ******
    
    List : a list 
        A list of items or nested list of any number of lists.
 
    
    Outputs:
    *******
    
    NumOfItemsByDepth : list of integers
        The number of items in List at each depth.
    """
    
    NumOfItemsByDepth = []
    
    while isinstance(List, list):
        NumOfItemsByDepth.append(len(List))
        
        NumOfItemsNextDepth = []
        
        for i in range(len(List)): 
            NumOfItemsNextDepth.append(len(List[i]))
        
        """ need to continue... """
        
    return NumOfItemsByDepth








def PrintTitle(Title):
    
    ul = '*' * len(Title)
    print('\n' + Title)
    print(ul + '\n')
    
    return





def PrintIndsByRoi(IndsByRoi):
    """ Print the list of indices in each ROI, e.g. from a list of C2SindsByRoi
    or F2SindsBySeg. """
    
    R = len(IndsByRoi)
    print(f'There are {R} ROIs/segments in IndsByRoi')
    
    C = [len(IndsByRoi[r]) for r in range(R)]
    print(f'There are {C} contours/frames in each ROI/segment')
    
    print('IndsByRoi:')
    
    for r in range(R):
        if r == 0:
            msg = f'   [{IndsByRoi[r]}'
        else:
            msg = f'    {IndsByRoi[r]}'
        
        if r < R - 1:
            msg += ','
        
        if r == R - 1:
            msg += ']'
                
        print(msg)
        
    return






def PrintPtsByCnt(PtsByCnt):
    """ Print the number of points in each contour in PtsByCnt. """
    
    C = len(PtsByCnt)
    print(f'There are {C} contours')
    
    for c in range(C):
        P = len(PtsByCnt[c])
        print(f'   There are {P} points in the {c}^th contour')
    
    return







def PrintPtsByCntByRoi(PtsByCntByRoi):
    """ Print the number of points in each contour of each ROI in 
    PtsByCntByRoi. """
    
    R = len(PtsByCntByRoi)
    print(f'There are two ROIs/segments in IndsByRoi')
    
    #C = [len(PtsByCntByRoi[r]) for r in range(R)]
    #print(f'There are {C} contours/frames in each ROI/segment')
    
    for r in range(R):
        C = len(PtsByCntByRoi[r])
        print(f'There are {C} contours in the {r}^th ROI')
        
        for c in range(C):
            P = len(PtsByCntByRoi[r][c])
            print(f'   There are {P} points in the {c}^th contour')
    
    return







def PrintPixArrBySeg(PixArrBySeg):
    """ If PixArrBySeg is a list of pixel arrays (as the input variable name
    suggests):
        - print the number of pixel arrays in the list
        - print the shape of each pixel array
    
    If PixArrBySeg is a pixel array (not a list as the variable name suggests):
        - print the shape of the pixel array """
    
    if isinstance(PixArrBySeg, list):
        S = len(PixArrBySeg)
        print(f'There are {S} pixel arrays in PixArrBySeg')
        
        for s in range(S):
            print(f'   The {s}^th pixel array has shape {PixArrBySeg[s].shape}')
    else:
        print(f'The pixel array has shape {PixArrBySeg.shape}')
    
    return





def PrintPixArrShapeBySeg(PixArrBySeg):
    #import numpy as np
    
    S = len(PixArrBySeg)
    
    print('\n   Shape of PixArrBySeg:')
    
    for s in range(S):
        if s == 0:
            msg = f'   [{PixArrBySeg[s].shape}'
        else:
            msg = f'    {PixArrBySeg[s].shape}'
        
        if s < S - 1:
            msg += ','
        
        if s == S - 1:
            msg += ']'
                
        print(msg)

    F = [PixArrBySeg[s].shape[0] for s in range(S)]
    
    print(f'   Number of frames in each segment = {F}')
    
    return






def IsPointInPolygon(point, vertices):
    """
    Determine if a point lies within a polygon.
    
    Inputs:
    ******
    
    Point : list of floats
        A list of length 3 of floats representing a 3D point.
                     
    Vertices : list of floats
        A list of length 8 of floats for all vertices that define a 3D volume.
        
        
    Outputs:
    *******
    
    PtInPoly : boolean
        True/False if Point is/isn't in Polygon.
        
        
    Note:
    ----
    
    This hasn't been tested on 2D points and 2D polygons.
    """
    
    #from shapely.geometry import Polygon, Point
    from shapely.geometry import Point, MultiPoint
    
    """ Use of Polygon requires that the points be ordered specifically. Use
    of MultiPoint and the convex hull allows for the points to be in a random
    order. """
    
    #polygon = Polygon(vertices)
    polygon = MultiPoint(vertices).convex_hull
    
    point = Point(point)
    
    PtInPoly = polygon.contains(point)
    
    return PtInPoly









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
    






def CropPixArr(PixArr):
    
    import numpy as np
    
    mask = PixArr.any(axis=0)
    
    Rinds, Cinds = np.where(mask)
    
    rMin = Rinds.min()
    rMax = Rinds.max()
    cMin = Cinds.min()
    cMax = Cinds.max()

    IndsExtent = [rMin, rMax, cMin, cMax]
    
    CropPixArr = PixArr[:, rMin:rMax+1, cMin:cMax+1]

    return CropPixArr, IndsExtent







def CropLabmapArr(LabmapArr):
    """
    Crop a 3D binary pixel array (e.g. "labelmap array") to a non-zero pixel
    array.
    
    Inputs:
    ******
    
    LabmapArr : Numpy array
        A binary pixel array containing binary masks along the 3rd dimension.
    
    
    Outputs:
    *******
    
    PixArr : Numpy array
        LabmapArr cropped within the bounds of all non-zero elements.    
    """
    
    import numpy as np
    
    """ Get the slices indices of non-zero frames: """
    Inds = []
        
    """ Loop through all frames in Nda and get the sum of all pixels: """
    for i in range(LabmapArr.shape[0]):
        Sum = np.amax(LabmapArr[i])
        
        if Sum:
            Inds.append(i)
            
    a = Inds[0]
    b = Inds[-1] + 1
    
    """ Get the indices that surround all non-zero elements along the other
    axes: """
    
    mask = LabmapArr.any(axis=0)
    
    Rinds, Cinds = np.where(mask)
    
    c = Rinds.min()
    d = Rinds.max() + 1
    e = Cinds.min()
    f = Cinds.max() + 1
    
    

    return LabmapArr[a:b, c:d, e:f]








def CropLabmapIm(LabmapIm, KeepAspect=False):
    """
    Crop a 3D binary pixel array from a labelmap image to a non-zero pixel
    array.
    
    Inputs:
    ******
    
    LabmapIm : SimpleITK image
        A binary image containing binary masks along the 3rd dimension.
    
    KeepAspect : boolean (optional; False by default)
        KeepAspect will determine the shape of PixArr.  If KeepAspect=False,
        the shape of PixArr will be the smallest cube that encloses all non-
        zero elements in LabmapArr.  If True the shape of PixArr will be such 
        that the relative spacings are maintained. 
        
        e.g. If Spacings=(0.5, 0.5, 1), PixArr will have shape (N, N, 2N), 
        where N = the maximum number of rows or columns that encloses all non-
        zero elements in LabmapArr. 
    
    Outputs:
    *******
    
    PixArr : Numpy array
        LabmapArr cropped within the bounds of all non-zero elements.    
    """
    
    import numpy as np
    import SimpleITK as sitk
    from ConversionTools import Image2PixArr
    
    LabmapArr = sitk.GetArrayViewFromImage(LabmapIm)
    
    PixArr, F2Sinds = Image2PixArr(LabmapIm)
    
    F, C, R = PixArr.shape
    
    R, C, S = LabmapIm.GetSize()
    
    di, dj, dk = LabmapIm.GetSpacing()
    
    
            
    imin = F2Sinds[0]
    imax = F2Sinds[-1] + 1
    
    """ Get the indices that surround all non-zero elements along the other
    axes: """
    
    mask = PixArr.any(axis=0)
    
    Rinds, Cinds = np.where(mask)
    
    jmin = Rinds.min()
    jmax = Rinds.max() + 1
    kmin = Cinds.min()
    kmax = Cinds.max() + 1
    
    
    if KeepAspect:
        # The number of non-zero rows and cols:
        N_j = jmax - jmin
        N_k = kmax - kmin
        
        if N_j > N_k:
            # Increase the number of columns:
            kmin = kmin - (N_j - N_k)//2
            kmax = kmin + N_j
        else:
            # Increase the number of rows:
            jmin = jmin - (N_k - N_j)//2
            jmax = jmin + N_k
            
        
        # Ratio of z-to-x spacing:
        dk_di = dk/di
        
        # The number of pixels along x or y:
        N = jmax - jmin
        
        # Scale by dk_di:
        imax = imin + int(round(N/dk_di))
        

    return LabmapArr[imin:imax, jmin:jmax, kmin:kmax]







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
    ******
    
    Image0 : SimpleITK image
    
    SliceNum0 : integer
        The slice number in Image0.
    
    Image1 : SimpleITK image
    
    SliceNum1 : integer
        The slice number in Image1.
    
    
    Outputs:
    *******
    
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
    ******
    
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
    *******
    
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






def ReplaceIndInC2SindsByRoi(C2SindsByRoi, IndToReplace, ReplacementInd):
    
    for r in range(len(C2SindsByRoi)):
        for c in range(len(C2SindsByRoi[r])):
            if C2SindsByRoi[r][c] == IndToReplace:
                C2SindsByRoi[r][c] = ReplacementInd
    
    return C2SindsByRoi






def GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1(SrcC2SindsByRoi, SrcImage, TrgImage):
    """ This doesn't seem to work for some reason. """
    
    TrgC2SindsByRoi = []
    
    for r in range(len(SrcC2SindsByRoi)):
        TrgC2Sinds = []
        
        for c in range(len(SrcC2SindsByRoi[r])):
            s = SrcC2SindsByRoi[r][c]
            
            """ Get the Source IPP for this slice: """
            IPP = SrcImage.TransformIndexToPhysicalPoint([0,0,s])
            
            """ Convert to an index in TrgIm: """
            Ind = TrgImage.TransformPhysicalPointToIndex(IPP)
            
            TrgC2Sinds.append(Ind[2])
            
        TrgC2SindsByRoi.append(TrgC2Sinds)
    
    return TrgC2SindsByRoi






def GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcC2SindsByRoi, SrcImage, TrgImage):
    
    """ Compare origins: """
    SrcOrigin = SrcImage.GetOrigin()
    TrgOrigin = TrgImage.GetOrigin()
    
    #dim = len(SrcOrigin)
    
    #mmDiff = [TrgOrigin[i] - SrcOrigin[i] for i in range(dim)]
    
    Dz_mm = TrgOrigin[2] - SrcOrigin[2]
    
    dk = TrgImage.GetSpacing()[2]
    
    Dz_pix = round(Dz_mm/dk)
    
    TrgC2SindsByRoi = []
    
    for r in range(len(SrcC2SindsByRoi)):
        TrgC2Sinds = []
        
        for c in range(len(SrcC2SindsByRoi[r])):
            #TrgC2Sinds.append(SrcC2SindsByRoi[r][c] + Dz_pix)
            TrgC2Sinds.append(SrcC2SindsByRoi[r][c] - Dz_pix)
            
        TrgC2SindsByRoi.append(TrgC2Sinds)
    
    return TrgC2SindsByRoi







def ShiftFrame(Frame, PixShift):
    """
    Shift the items in a 2D frame.
    
    Inputs:
    ******
    
    Frame : Numpy array
        A 2D frame from PixelData with shape (1, R, C), where R = number of 
        rows, and C = number of columns.
        
    PixShift : Numpy array
        Shift in pixels, e.g. [di, dj, dz].
        
    
    Outputs:
    *******
    
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







def ShiftFramesInPixArrBySeg(PixArrBySeg, F2SindsBySeg, SrcImage, SrcSliceNum,
                             TrgImage, TrgSliceNum, RefImage, ShiftInX=True,
                             ShiftInY=True, ShiftInZ=True, Fractional=False,
                             LogToConsole=False):
    """
    Shift the items in each frame of each pixel array in a list of 
    pixel-arrays-by-segment.  Example uses: To compensate for differences in 
    z-positions of a single-framed Source pixel array to be "direct" copied to
    a Target slice at a different stack position; or compensating for 
    differences in the origin of Source and Target images for relationship-
    preserving copies of images with equal voxel spacings.
    
    Inputs:
    ******
    
    PixArrBySeg : list of Numpy arrays
        A list (for each segment) of pixel arrays.
    
    F2SindsBySeg : list of a list of integers
        List (for each segment) of a list (for each frame) of slice numbers that 
        correspond to each frame in the pixel arrays in PixArrBySeg.
        
    SrcImage : SimpleITK image
        The Source 3D image.
    
    SrcSliceNum : integer
        The slice number within the Source image stack that the segment relates
        to.
    
    TrgImage : SimpleITK image
        The Target 3D image.
    
    TrgSliceNum : integer
        The slice number within the Target image stack that the segment is to  
        be copied to following the shift in z-components.
        
    RefImage : SimpleITK image
        The reference image whose voxel spacings will be used to convert apply
        the shift in pixels (e.g. TrgImage).
    
    ShiftInX / ShiftInY / ShiftInZ : boolean (optional; True by default)
        If True the pixel shift between the origins of SrcImage and TrgImage
        will be applied along the x / y / z direction.
    
    Fractional : boolean (optional; False by default)
        If True the pixel shift will be rounded to the nearest integer value.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    PixArrBySeg : list of Numpy arrays
        As above but with shifted frame with shape (1, R, C).
    
    F2SindsBySeg : list of a list of integers
        As above but with shifted slice indices that correspond to the shifted 
        PixArrBySeg.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ShiftFramesInPixArrBySeg():')
        
    # Get pixel shift between FromSliceNum in SrcImage and ToSliceNum in 
    # TrgImage in the Target image domain:
    PixShift = GetPixelShiftBetweenSlices(Image0=SrcImage, 
                                          SliceNum0=SrcSliceNum, 
                                          Image1=TrgImage, 
                                          SliceNum1=TrgSliceNum,
                                          RefImage=TrgImage)
    
    if LogToConsole:
        print(f'   \nPixShift between slices = {PixShift}')
        print(f'   \nF2SindsBySeg prior to shifting = {F2SindsBySeg}')
        
        #print(f'\n   Prior to shifting frames:')
        ##print(f'PixArrBySeg = {PixArrBySeg}')
        #print(f'   len(PixArrBySeg) = {len(PixArrBySeg)}')
        #for r in range(len(PixArrBySeg)):
        #    print(f'      type(PixArrBySeg[{r}]) = {type(PixArrBySeg[r])}')
        #    print(f'      PixArrBySeg[{r}].shape = {PixArrBySeg[r].shape}')
    
    if not ShiftInX:
        PixShift[0] = 0
    if not ShiftInY:
        PixShift[1] = 0
    if not ShiftInZ:
        PixShift[2] = 0
    
    if LogToConsole:
        print(f'   \nPixShift that will be applied = {PixShift}')
        
    for s in range(len(PixArrBySeg)):
        #if PixArrBySeg[s]:
        if PixArrBySeg[s].shape[0]:
            #print('   \nBefore shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            # Replace PixArrBySeg[s] with the result of shifting the in-plane 
            # elements, and replace F2SindsBySeg[s] with [TrgSliceNum]:
            PixArrBySeg[s] = ShiftFrame(Frame=PixArrBySeg[s], PixShift=PixShift)
            
            #print('   \nAfter shifting frame:')
            #print(f'   type(PixArrBySeg[{s}]) = {type(PixArrBySeg[s])}')
            #print(f'   PixArrBySeg[{s}].shape = {PixArrBySeg[s].shape}')
            
            #F2SindsBySeg[s] = [TrgSliceNum] # 11/02
            
            # Shift the frame-to-slice indices by PixShift[2]:
            F2SindsBySeg[s] = [F2SindsBySeg[s][i] + PixShift[2] for i in range(len(F2SindsBySeg[s]))]
            
            if LogToConsole:
                unique = UniqueItems(Items=PixArrBySeg[s], IgnoreZero=False)
                
                print(f'\nThere are {len(unique)} unique items in',
                      f'PixArrBySeg[{s}] after shifting the frame.')
    
    if LogToConsole:
        print(f'   F2SindsBySeg after shifting = {F2SindsBySeg}')
        print('-'*120)
        
    return PixArrBySeg, F2SindsBySeg
    
    





def ChangeZinds(Indices, NewZind):
    """
    Re-assign the z-components of a list of indeces (as required when
    making a direct copy of contour points).
    
    Inputs:
    ******
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of indices, 
        e.g. [[i0, j0, k0], [i1, j1, k1], ...].
        
    NewZind : integer
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    
    
        
    Ouputs:
    ******
    
    Indices : list of a list of integers
        List (for each point) of a list (for each dimension) of the indices 
        with modified z-components, e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    """
    
    if not isinstance(NewZind, int):
        raise Exception(f"NewZind = {NewZind} must be an integer.")
    
    for i in range(len(Indices)):
        Indices[i][2] = NewZind
        
    return Indices







def ZshiftPtsByCntByRoi(PtsByCntByRoi, NewZind, DicomDir):
    """
    Shift the z-component of the points in PtsByCntByRoi (as required when
    making a direct copy of contour points).
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
    NewZind : integer
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
        
    Ouputs:
    ******
    
    NewPtsByCntByRoi : list of list of a list of a list of floats
        As PtsByCntByRoi but with modified z-components.
    """
    
    from ImageTools import ImportImage
    from ConversionTools import Points2Indices, Indices2Points
    
    # Import the image:
    Im = ImportImage(DicomDir)
    
    """ Convert PtsByCntByRoi to indices, modify the z-component of the
    indices, then convert back to physical points. """
    
    #IndsByCntByRoi = []
    NewPtsByCntByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(PtsByCntByRoi)): 
        #IndsByCnt = []
        NewPtsByCnt = []
        
        if PtsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(PtsByCntByRoi[r])):
                pts = PtsByCntByRoi[r][c]
                
                inds = Points2Indices(Points=pts, RefIm=Im, Rounding=False)
                
                # Modify the z-component (k) of the indices to NewZind:
                inds = ChangeZinds(inds, NewZind)
                
                #IndsByCnt.append(inds)
                
                # Convert back to physical points:
                pts = Indices2Points(inds, RefIm=Im)
                
                NewPtsByCnt.append(pts)
        
        #IndsByCntByRoi.append(IndsByCnt)
        NewPtsByCntByRoi.append(NewPtsByCnt)
        
    return NewPtsByCntByRoi





def ShiftPtsByCntByRoi(PtsByCntByRoi, C2SindsByRoi, SrcImage, SrcSliceNum, 
                       TrgImage, TrgSliceNum, RefImage, ShiftInX=True,  
                       ShiftInY=True, ShiftInZ=True, Fractional=False, 
                       LogToConsole=False):
    
    #NewZind, DicomDir, 
    
    """
    Shift the points in each list of points in each list of contours in each
    list of ROIs.  Example uses: To compensate for differences in 
    z-positions of a single-contour Source points to be "direct" copied to
    a Target slice at a different stack position; or compensating for 
    differences in the origin of Source and Target images for relationship-
    preserving copies of images with equal voxel spacings.
    
    Inputs:
    ******
    
    PtsByCntByRoi : list of list of a list of a list of floats
        List (for each ROI) of a list (for all contours) of a list (for each
        point) of a list (for each dimension) of coordinates.
        
    NewZind : integer
        The value to re-assign to the z-components (i.e. k) of Indeces (e.g.
        ToSliceNum).
    
    DicomDir : string 
        Directory containing the corresponding DICOMs.
    
        
    Ouputs:
    ******
    
    NewPtsByCntByRoi : list of list of a list of a list of floats
        As PtsByCntByRoi but with modified z-components.
    """
    
    #from ImageTools import ImportImage
    from ConversionTools import Points2Indices, Indices2Points
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ShiftPtsByCntByRoi():')
        
    """ Get pixel shift between FromSliceNum in SrcImage and ToSliceNum in 
    TrgImage in the Target image domain: """
    PixShift = GetPixelShiftBetweenSlices(Image0=SrcImage, 
                                          SliceNum0=SrcSliceNum, 
                                          Image1=TrgImage, 
                                          SliceNum1=TrgSliceNum,
                                          RefImage=TrgImage)
    
    if LogToConsole:
        print(f'   PixShift between slices = {PixShift}')
    
    if not ShiftInX:
        PixShift[0] = 0
    if not ShiftInY:
        PixShift[1] = 0
    if not ShiftInZ:
        PixShift[2] = 0
    
    if LogToConsole:
        print(f'   PixShift that will be applied = {PixShift}')
    
    
    """ Convert PtsByCntByRoi to indices, modify the indices, then convert back
    to physical points. """
    
    #IndsByCntByRoi = []
    NewPtsByCntByRoi = []
    NewC2SindsByRoi = []
    
    # Iterate through ROIs:
    for r in range(len(PtsByCntByRoi)): 
        #IndsByCnt = []
        NewPtsByCnt = []
        NewC2Sinds = []
        
        if PtsByCntByRoi[r]:
            # Iterate through contours:
            for c in range(len(PtsByCntByRoi[r])):
                Pts = PtsByCntByRoi[r][c]
                
                """ Inds are the list of relationship-preserving indices that
                correspond to Pts: """
                Inds = Points2Indices(Points=Pts, RefIm=TrgImage, 
                                      Rounding=False)
                
                #print(f'\n\nPixShift = {PixShift}, \nInds = {Inds}')
                
                """ Apply the required pixel shifts: """
                Inds = [Inds[i] + PixShift[i] for i in range(len(PixShift))]
                
                #IndsByCnt.append(inds)
                #print(f'\n\ntype(Inds) = {type(Inds)}')
                #print(f'\n\ntype(Inds[0]) = {type(Inds[0])}')
                
                """ Convert back to physical points: """
                Pts = Indices2Points(Inds, RefIm=TrgImage)
                
                NewPtsByCnt.append(Pts)
                
                """ Shift the contour-to-slice indices by PixShift[2]: """
                #C2SindsByRoi[r][c] += PixShift[2]
                NewC2Sinds.append(C2SindsByRoi[r][c] + PixShift[2])
        
        #IndsByCntByRoi.append(IndsByCnt)
        NewPtsByCntByRoi.append(NewPtsByCnt)
        NewC2SindsByRoi.append(NewC2Sinds)
        
    return NewPtsByCntByRoi, NewC2SindsByRoi






def MeanFrameInPixArr(PixArr, F2Sinds, MakeBinary=True, BinaryThresh=0.5,
                      LogToConsole=False):
    """
    Perform pixel-by-pixel mean of all frames in a pixel array. Output a 
    binary pixel array if binary=True (or undefined).
    
    Inputs:
    ******
    
    PixArr : Numpy array
        PixelData from a SEG file loaded as a Numpy array
        
    F2Sinds : List of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.
        
    MakeBinary : boolean (optional; True by default)
        If MakeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by BinaryThreshold.
    
    BinaryThreshold : float (optional; 0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if MakeBinary = True.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
        
    Outputs:
    *******
    
    MeanPixArr : Numpy array
        Mean pixel array.
    
    MeanF2Sind : integer or float
        The integer or fractional frame-to-slice index corresponding to the 
        mean pixel array.  MeanF2Sind will be an integer if MeanF2Sind % 1 = 0,
        and a fractional float otherwise.
    """
    
    import numpy as np
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of MeanFrameInPixArr():')
        
    F, R, C = PixArr.shape
    
    if LogToConsole:
        print(f'\nPixArr.shape = {PixArr.shape}')
        #print(f'np.amax(PixArr, axis=0) = {np.amax(PixArr, axis=0)}')
        #print(f'np.max(PixArr, axis=0) = {np.max(PixArr, axis=0)}')
        print(f'Max along each frame = {[np.amax(PixArr[f]) for f in range(PixArr.shape[0])]}')
    
    MeanPixArr = np.mean(PixArr, axis=0)
    
    if LogToConsole:
        print(f'Max of MeanPixArr after averaging = {np.amax(MeanPixArr)}')
            
    if MakeBinary:
        #MeanPixArr = (MeanPixArr >= BinaryThresh) * MeanPixArr
        MeanPixArr = (MeanPixArr >= BinaryThresh) * 1
        
        if LogToConsole:
            print('Max of MeanPixArr after binary thresholding =',
                  f'{np.amax(MeanPixArr)}')
    
    MeanPixArr = np.reshape(MeanPixArr, (1, R, C))
    
    if LogToConsole:
        print(f'Max of MeanPixArr after reshape = {np.amax(MeanPixArr)}')
    
    MeanPixArr = MeanPixArr.astype(dtype='uint8')
    
    if LogToConsole:
        print(f'Max of MeanPixArr after change to uint8 = {np.amax(MeanPixArr)}')
    
    """ Find the mean slice index: """
    MeanF2Sind = sum(F2Sinds)/len(F2Sinds)
    
    if MeanF2Sind % 1 == 0:
        MeanF2Sind = int(MeanF2Sind)
    
    if LogToConsole:
        print('-'*120)
        
    return MeanPixArr, MeanF2Sind






def MeanFrameInPixArrBySeg(PixArrBySeg, F2SindsBySeg, MakeBinary=True, 
                           BinaryThresh=0.5, LogToConsole=False):
    """
    Perform pixel-by-pixel mean of all frames in all pixel arrays in a list
    of pixel-arrays-by-ROI which have more than 1 frame. No changes will be
    made to empty or single-framed pixel arrays.  If there are less than 2
    frames in any given pixel array the pixel array will be returned unchanged.
    Likewise for the F2SindsBySeg.
    
    Inputs:
    ******
    
    PixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    
    F2SindsBySeg : list of a list of integers
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in PixArrBySeg.
        
    MakeBinary : boolean (optional; True by default)
        If MakeBinary = True the mean pixel array will be converted to a binary
        pixel array by thresholding pixel values by BinaryThreshold.
    
    BinaryThreshold : float (optional; 0.5 by default)
        The threshold that will be used when converting to a binary pixel array
        if MakeBinary = True.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    MeanPixArrByRoi : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
        
    MeanF2SindByRoi : list of a list (of length 1) of an integer or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the mean pixel array.  MeanF2Sind will be 
        an integer if MeanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of MeanFrameInPixArrBySeg():')
        
    MeanPixArrBySeg = []
    MeanF2SindBySeg = []
    
    for s in range(len(PixArrBySeg)):
        PixArr = PixArrBySeg[s]
        F2Sinds = F2SindsBySeg[s]
        
        """ The number of frames in PixArr """
        F = PixArr.shape[0]
        
        if F > 1: 
            if LogToConsole:
                print(f'\nPixArrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be averaged.')
            
            MeanPixArr,\
            MeanF2Sind = MeanFrameInPixArr(PixArr, F2Sinds, MakeBinary, 
                                           BinaryThresh, LogToConsole)
            
            MeanPixArrBySeg.append(MeanPixArr)
            MeanF2SindBySeg.append([MeanF2Sind])
            
            if LogToConsole:
                print(f'\nPixArr had {F} frames but now has',
                      f'{MeanPixArr.shape[0]}.')
        else:
            MeanPixArrBySeg.append(PixArr)
            
            MeanF2SindBySeg.append(F2Sinds)
            
            if LogToConsole:
                print(f'\nPixArr has {F} frames so no frame averaging required.')
    
    if LogToConsole:
        print('-'*120)
        
    return MeanPixArrBySeg, MeanF2SindBySeg







def OrFrameOfPixArrBySeg(PixArrBySeg, F2SindsBySeg, LogToConsole=False):
    """
    Perform pixel-by-pixel logical OR of all frames in all pixel arrays in a 
    list of pixel-arrays-by-segment which have more than 1 frame. No changes 
    will be made to empty or single-framed pixel arrays. If there are less than
    2 frames in any given pixel array the pixel array will be returned 
    unchanged.  Likewise for the F2SindsBySeg.  
    
    Inputs:
    ******
    
    PixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays.
    
    F2SindsBySeg : list of a list of integers
        A list (for each ROI/segment) of a list (for each contour/frame) of the 
        slice numbers that correspond to each contour/frame in each pixel array
        in PixArrBySeg.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    OrPixArrBySeg : list of Numpy arrays
        A list (for each ROI/segment) of pixel arrays averaged to a single-
        framed pixel array.
    
    MeanF2SindBySeg : list of a list (of length 1) of an integer or float
        A list (for each ROI/segment) of a list (of length 1, since there is 
        only 1 frame per pixel array) of the integer or fractional frame-to-
        slice index corresponding to the OR pixel array.  MeanF2Sind will be 
        an integer if MeanF2Sind % 1 = 0, and a fractional float otherwise.
    """
    
    import numpy as np
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of OrFrameOfPixArrBySeg():')
        
    # Get the shape of any pixel array:
    F, R, C = PixArrBySeg[0].shape
    
    OrPixArrBySeg = []
    MeanF2SindBySeg = []
    
    for s in range(len(PixArrBySeg)):
        PixArr = PixArrBySeg[s]
        F2Sinds = F2SindsBySeg[s]
        
        """ The number of frames in PixArr """
        F = PixArr.shape[0]
        
        if F > 1: 
            if LogToConsole:
                print(f'\nPixArrBySeg[{s}] has {F} frames so the pixel arrays',
                      'will be logically ORed.')
            
            OrPixArr = PixArr.any(axis=0)
            
            """ Reshape to a (1, R, C) pixel array: """
            OrPixArr = np.reshape(OrPixArr, (1, R, C))
            
            """ Convert from boolean to uint8: """
            OrPixArr = OrPixArr.astype('uint8')
            
            OrPixArrBySeg.append(OrPixArr)
            
            """ Find the mean slice index: """
            MeanF2Sind = sum(F2Sinds)/len(F2Sinds)
            
            """
            14/02: It may be useful to return the fractional index but for now
            I'll simply round it to the nearest integer...
            if MeanF2Sind % 1 == 0:
                MeanF2Sind = int(MeanF2Sind)
            """
            MeanF2Sind = round(MeanF2Sind)
            
            MeanF2SindBySeg.append([MeanF2Sind])
            
            if LogToConsole:
                print(f'\nPixArr had {F} frames but now has,',
                      f'{OrPixArr.shape[0]}.')
        else:
            OrPixArrBySeg.append(PixArr)
            
            MeanF2SindBySeg.append(F2Sinds)
            
            if LogToConsole:
                print(f'\nPixArr has {F} frames so no logical ORing required.')
    
    if LogToConsole:
        print(f'MeanF2SindBySeg = {MeanF2SindBySeg}')
        print('-'*120)
        
    return OrPixArrBySeg, MeanF2SindBySeg







def Distance(item0, item1):
    """
    Compute the distance between points or indices.  
    
    Inputs:
    ******
    
    item0 : list of floats/integers
        A 2D/3D point or point indices.
        
    item1 : list of floats/integers
        A 2D/3D point or point indices.
    
    
    Outputs:
    *******
    
    distance : float/integer
        Distance between item0 and item1 in Euclidean space or in pixels.
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





def GetListOfDataTypes(ListOfItems, Unique=True):
    """
    Return a list of the data type of the elements in the list.  If the input
    Unique is True or omitted, the set of data types will be returned.  Hence
    a single-itemed list containing a single data type is a permitted output.
    
    Inputs:
    *******
    
    ListOfItems : list of items
        A list of items (e.g. integers, floats, strings, etc.).
        
    Unique : boolean (optional; True by default)
        If True the set of data types will be returned (i.e. only unique data
        types).
        
        
    Outputs:
    *******
        
    Dtypes : list of data types
        A list of the data type of each item in ListOfItems.
    """
    
    Dtypes = [type(item) for item in ListOfItems]
    
    if Unique:
        Dtypes = list(set(Dtypes))
    
    return Dtypes








def GetVectorLength(Vector):
    """
    Get vector length.
    
    Inputs:
    ******
    
    Vector : list of floats
        A list (for each dimension) of vector components along each axis,
        e.g. [Vx, Vy, Vz] or [Vx, Vy].
        
        
    Outputs:
    *******
    
    VectorL : float
        The length of Vector.
    """
    
    dim = len(Vector)
        
    if dim == 2:
        VectorL = ( Vector[0]**2 + Vector[1]**2 ) ** 0.5
        
    elif dim == 3:
        
        VectorL = ( Vector[0]**2 + Vector[1]**2 + Vector[2]**2 ) ** 0.5
        
    else:
        msg = f'The vector has dim {dim}. It must have 2 or 3 dimensions.'
        
        raise Exception(msg)
    
        
    return VectorL






def SumOfSquaresOfMatrix(Matrix):
    """
    Compute the root-sum-of-squares of a matrix.
    
    Inputs:
    ******
    
    Matrix : Numpy array
        A Numpy array (including matrices).
        
        
    Outputs:
    *******
    
    RSS : float
        The root-sum-of-squares of Matrix.
    """
    
    from numpy import sqrt, sum, square
    
    RSS = sqrt(sum(square(Matrix)))
    
    return RSS








def IsMatrixOrthogonal(Matrix, LogToConsole=False):
    """
    Determine if a matrix is orthogonal.
    
    Inputs:
    ******
    
    Matrix : Numpy array of square shape
        A Numpy matrix of square shape.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    IsOrthogonal : boolean
        True if Matrix is orthogonal, False otherwise.
        
        
    Notes:
    -----
    
    The following tests for orthogonality of a matrix M: 
        
        M is orthogonal if M*M.T = M.T*M, 
        i.e. if np.dot(M, M.T) - np.dot(M.T, M) = 0
    """
    
    import numpy as np
    
    M = Matrix
    
    if len(M.shape) != 2:
        msg = f'The matrix has {len(M.shape)} dimensions. It must have 2.'
        
        raise Exception(msg)
        
    if M.shape[0] != M.shape[1]:
        msg = f'The matrix has shape {M.shape}. It must be a square matrix.'
        
        raise Exception(msg)
        
    
    result =  np.dot(M, M.T) - np.dot(M.T, M)
    
    SSE = SumOfSquaresOfMatrix(result)
    
    IsOrthogonal = AreItemsEqualToWithinEpsilon(SSE, 0, epsilon=1e-06)
    
    if LogToConsole:
        print(f'\nM*M.T - M.T*M = \n{result}')
        print(f'\nSSE = {SSE}')
        
        if IsOrthogonal:
            print(f'\nM is orthogonal')
        else:
            print(f'\nM is not orthogonal')
        
        
    return IsOrthogonal






def IsMatrixOrthonormal(Matrix, LogToConsole=False):
    """
    Determine if a matrix is orthonormal.
    
    Inputs:
    ******
    
    Matrix : Numpy array of square shape
        A Numpy matrix of square shape.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
        
    Outputs:
    *******
    
    IsOrthonormal : boolean
        True if Matrix is orthonormal, False otherwise.
        
        
    Notes:
    -----
    
    The following tests for orthonomality of a matrix M: 
        
        M is orthonormal if its determinant (np.linalg.det(M)) = +/- 1.
        
        M is orthonormal if M*M.T = I, where I is the identity matrix, 
        i.e. if np.dot(M, M.T) - np.identity(M.shape[0]) = 0
        
        M is orthonormal if M.T = M^-1, 
        i.e. if M.T - np.linalg.inv(M) = 0
    """
    
    import numpy as np
    
    M = Matrix
    
    #print(type(M))
    print(M)
    
    if len(M.shape) != 2:
        msg = f'The matrix has {len(M.shape)} dimensions. It must have 2.'
        
        raise Exception(msg)
        
    if M.shape[0] != M.shape[1]:
        msg = f'The matrix has shape {M.shape}. It must be a square matrix.'
        
        raise Exception(msg)
        
    
    """ Either of the three tests below (Det = 1?, M.M_T = 0?, 
    SSE of M.T - M^-1 = 0?) would be suitable. Determinant test is chosen for
    no particular reason. """
    Det = np.linalg.det(M)
    IsOrthonormal = AreItemsEqualToWithinEpsilon(abs(Det), 1, epsilon=1e-06)
    
    MdotM_T =  np.dot(M, M.T) - np.identity(M.shape[0])
    SSEofMdotM_T = SumOfSquaresOfMatrix(MdotM_T)
    #IsOrthonormal = AreItemsEqualToWithinEpsilon(SSEofMdotM_T, 0, 
    #                                             epsilon=1e-06)
    
    TequalsInv = M.T - np.linalg.inv(M)
    SSEofTequalsInv = SumOfSquaresOfMatrix(TequalsInv)
    #IsOrthonormal = AreItemsEqualToWithinEpsilon(SSEofTequalsInv, 0, 
    #                                             epsilon=1e-06)
    
    
    if LogToConsole:
        print(f'\nThe determinant of M = {Det}')
        print(f'\nM*M.T - I = \n{MdotM_T}')
        print(f'\nSSE = {SSEofMdotM_T}')
        print(f'\nM.T - M^-1 = \n{SSEofTequalsInv}')
        print(f'\nSSE = {SSEofTequalsInv}')
        
        if IsOrthonormal:
            print(f'\nM is orthonormal')
        else:
            print(f'\nM is not orthonormal')
        
    return IsOrthonormal






def ReduceDictToKeysMatchingKeyword(Dict, Keyword='', KeepDirectories=True,
                                    KeepFilepaths=True):
    """ Reduce a dictionary to a dictionary containing only the elements in
    Dict whose key matches a keyword. If the keyword is empty (''), the 
    elements won't be filtered by key.  Elements can be filtered further based
    on whether their keys are paths to directories or files. """
    
    import os
    
    NewDict = {}
    
    for Key, Val in Dict.items():
        if Keyword:
            """ Filter by keys that match Keyword: """
            if Keyword in Key:
                """ Filter by keys that are directory paths: """
                if os.path.isdir(Val) and KeepDirectories:
                    NewDict[Key] = Val
                """ Filter by keys that are file paths: """
                if os.path.isfile(Val) and KeepFilepaths:
                    NewDict[Key] = Val
        else:
            """ Filter by keys that are directory paths: """
            if os.path.isdir(Val) and KeepDirectories:
                NewDict[Key] = Val
            """ Filter by keys that are file paths: """
            if os.path.isfile(Val) and KeepFilepaths:
                NewDict[Key] = Val
    
    return NewDict




def ReduceDictToKeysMatchingKeywords(Dict, ListOfKeywords, 
                                     KeepDirectories=True, KeepFilepaths=True):
    """ Reduce a dictionary to a dictionary containing only the elements in
    Dict whose key matches all of the keywords in a list of keyword. """
    
    
    for Keyword in ListOfKeywords:
        Dict = ReduceDictToKeysMatchingKeyword(Dict, Keyword, KeepDirectories,
                                               KeepFilepaths)
    
    return Dict





def DisplayAttributesDict(Dict):
    """ Display dictionary containing image attributes as a pandas data frame
    and return the data frame. 
    
    See https://stackoverflow.com/questions/49661018/displaying-embedded-newlines-in-a-text-column-of-a-pandas-dataframe
    
    for two options below.
    """
    
    import pandas as pd
    from IPython.display import display, HTML
    
    df = pd.DataFrame.from_dict(Dict, orient ='index')
    
    #option = 1
    option = 2
    
    if option == 1:
        display(HTML(df.to_html().replace("\\n","<br>")))
    else:
        display(df.style.set_properties(**{'text-align': 'center',
                                           'white-space': 'pre-wrap'}))
    
    return df




def DisplayAttributesOfAllData(Keyword='Dcm', RemoveFromKeys='_DcmDir',
                               Package='pydicom'):
    """ Display the image attributes of all DICOMs whose directories are listed
    in CreateDictOfPaths() as a pandas dataframe. Optional inputs allow for:
    - filtering of paths in the dictionary (containing root directories of 
    sessions, directories containing DICOM series, and filepaths of RTS/SEG
    files;
    - removing of character strings from the keys in the dictionary of paths
    when creating the dictionary of image attributes;
    - selection of the package to use when obtaining the image attributes. """
    
    import time
    from RoiCopyTools import CreateDictOfPaths
    from ImageTools import GetImageAttributesFromDict
    #from GeneralTools import DisplayAttributesDict
    
    Dict = CreateDictOfPaths()
    
    print(f'There are {len(list(Dict.keys()))} keys in Dict\n')
    
    ##Package = 'sitk'
    #Package = 'pydicom'
    
    """ Start timing. """
    times = []
    times.append(time.time())
        
    Dict, AttDict = GetImageAttributesFromDict(Dict, Keyword, Package, 
                                               RemoveFromKeys)
    
    print(f'\nThere are {len(list(AttDict.keys()))} keys in AttDict\n')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\n*Took {Dtime} s using {Package}.')
    
    AttDf = DisplayAttributesDict(AttDict)
    
    return AttDf






def ExportDictionaryToJson(Dictionary, FileName, ExportDir='cwd'):
    """
    Export a dictionary to a JSON file.
    
    Inputs:
    ******
    
    Dictionary : dictionary
        The dictionary to be exported.
        
    FileName : string
        The file name to assign to the exported file. If FileName doesn't 
        include the '.json' extension it will be added automatically.
    
    ExportDir : string (optional; 'cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    import json
    import os
    from pathlib import Path
    
    if not '.json' in FileName:
        FileName += '.json'
    
    if ExportDir == 'cwd':
        FilePath = FileName
    else:
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        FilePath = os.path.join(ExportDir, FileName)
    
    with open(FilePath, 'w') as JsonFile:
        json.dump(Dictionary, JsonFile)
    
    return




def ImportDictionaryFromJson(FileName):
    """
    Import a dictionary from a JSON file from the current working directory.
    
    Inputs:
    ******
    FileName : string
        The file name of the JSON file (which must be in the current working 
        directory).  FileName does not need to include the .json extension.
        
        
    Outputs:
    *******
        
    Dictionary : dictionary
    """
    
    import json
    
    if not '.json' in FileName:
        FileName += '.json'
    
    with open(FileName, 'r') as JsonFile:
        Dictionary = json.load(JsonFile)
    
    return Dictionary





def ExportListToTxt(Lst, FileName, ExportDir='cwd'):
    """
    Export a list to a txt file.
    
    Inputs:
    ******
    
    Lst : list
        The list to be exported.
        
    FileName : string
        The file name to assign to the exported file. If FileName doesn't 
        include the '.txt' extension it will be added automatically.
    
    ExportDir : string (optional; 'cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    import os
    from pathlib import Path
    
    if not '.txt' in FileName:
        FileName += '.txt'
    
    if ExportDir == 'cwd':
        FilePath = FileName
    else:
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        FilePath = os.path.join(ExportDir, FileName)
    
    with open(FilePath, 'w') as f:
        for line in Lst:
            f.write(line)
    
    return





def ImportListFromTxt(FileName):
    """
    Import a list from a text file from the current working directory.
    
    Inputs:
    ******
        
    FileName : string
        The file name of the text file (which must be in the current working 
        directory).  FileName does not need to include the .txt extension.
    
    
    Outputs:
    *******
        
    Lst : list
    """
    
    #import os
    #from pathlib import Path
    
    if not '.txt' in FileName:
        FileName += '.txt'
    
    Lst = []
    
    with open(FileName, 'r') as file:
        for line in file:
            Lst.append(line.replace('\n', '').split(' '))
    
    return Lst




def DeleteFile(Filepath):
    import os
    
    try:
        os.remove(Filepath)
        print(f'{Filepath} deleted')
    except PermissionError as error:
        print(error)
        print(f'{Filepath} cannot be deleted\n')



def DeleteDirectory(Dirpath):
    import shutil
    
    try:
        shutil.rmtree(Dirpath)
        print(f'{Dirpath} deleted')
    except PermissionError as error:
        print(error)
        print(f'{Dirpath} cannot be deleted\n')



def ReduceListOfStringFloatsTo16(OrigList):
    """ Decimal strings (VR DS in DICOM standard) must be limited to 16 
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
    
    NewList = [item[:16] for item in OrigList]
    
    return NewList