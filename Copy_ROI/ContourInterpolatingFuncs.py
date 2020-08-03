# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:50:23 2020

@author: ctorti
"""



""" 

******************************************************************************
******************************************************************************
Contour interpolating functions
******************************************************************************
******************************************************************************


The following dictionaries are used as inputs to some of the following 
functions and are appended to by others:
   
***************
FixContourData:
***************

A dictionary of the contour points belonging to the Fixed (input) image, 
containing the keys:

    
  'SliceNo'     - List of integers (for original contour points) or floats (for
                  interpolated contours that correspond to non-imaged planes,
                  e.g. slice interpolated at 2.5 in between slices 2 and 3).
                  
  'ContourType' - List of integers for each DICOM slice; 
                  The integer value denotes whether there is or is not contour 
                  data for any given slice, and if there is, whether it is an
                  original contour, an interpolated contour, or transformed 
                  contour.
                  0 -> No contour points 
                  1 -> Original contour points
                  2 -> Interpolated contour points
                  3 -> Transformed contour points
                  
  'PointPCS'    - List of input points in the Patient Coordinate System if 
                  contour data exists for any given slice, e.g.:
                  [[x0, y0, z0], [x1, y1, z1], ...], 
                  i.e. a list of lists of lists), or [] if not.
                
  'PointICS'    - Same as 'PointPCS' but with points in the Image Coordinate
                  System.

              
**********                          
PointData:
**********        

A dictionary containing the keys:
    
  'PointNo'      - List of integers of all points in ContourData.
  
  'InSliceNo'    - List of integers (for original contour points) or floats 
                   (for interpolated contours that correspond to non-imaged
                   planes, e.g. slice interpolated at 2.5 in between slices
                   2 and 3).
                  
  'InContourNo'  - List of integers of the contour number each point belongs to
                   from the input image.
  
  'ContourType'  - List of integers for each DICOM slice; 
                   The integer value denotes whether there is or is not contour
                   data for any given slice, and if there is, whether it is an 
                   original contour, an interpolated contour, or transformed 
                   contour.
                   0 -> No contour points 
                   1 -> Original contour points
                   2 -> Interpolated contour points
                   3 -> Transformed contour points
                  
  'InPointIndex' - List of indeces [x, y, z] for each point in the input image.
                   
                  
  'InPointPCS'   - List of points in the Patient Coordinate System if contour 
                   data exists for any given slice, e.g.:
                   [[x0, y0, z0], [x1, y1, z1], ...], 
                   i.e. a list of lists of lists), or [] if not.
                
  'InPointICS'   - Same as 'InPointPCS' but with points in the Image 
                   Coordinate System.
                   
  'OutPointIndex - Like 'InPointIndex' but for the output/transformed image.
                   
  'OutPointPCS'  - Like 'InPointPCS' but for the output/transformed image.
                
  'OutPointICS' - Same as 'OutPointPCS' but with points in the Image
                  Coordinate System.
           
            
***********
InterpData:
***********
    
A dictionary containing the keys:
    
  'InterpSliceInd'       - A list of the indices of the interpolated slices 
                          (can be floats if interpolating between two 
                          neighbouring slices), e.g.: [13.5, 14.5, ...]
                          
  'BoundingSliceInds'    - A list of a list of the indices of the nearest 
                           slices that surround the plane to be interpolated, 
                           e.g.: [ [13, 14], [14, 15], ...]
                           
  'FixOSIsOrigNode1'     - A list of booleans whose value denotes whether the 
                           point/node in the list FixOSContour1Pts were 
                           original (True) or added (False).
                           
  'FixOSIsOrigNode2'     - Like FixOSIsOrigNode1 but for FixOSContour2Pts.
  
  'FixOSContour1Inds'    - A list of a list of [x, y, z] indices for each point
                           in FixOSContour1Pts.
                           
  'FixOSContour2Inds'    - Like FixOSContour1Inds but for FixOSContour2Pts.
  
  'FixInterpContourInds' - Like FixOSContour1Inds but for FixInterpContourPts.
  
  'FixOSContour1Pts'     - A list of a list of [x, y, z] coordinates for each
                           point in the over-sampled Contour1 for the Fixed 
                           image.
                           
  'FixOSContour2Pts'     - Like FixOSContour2Pts but for Contour2.
  
  'FixInterpContourPts'  - Like FixOSContour2Pts but for the contour 
                           interpolated from FixOSContour1Pts and 
                           FixOSContour2Pts. 
                           
  'MovOSContour1Inds'    - Like FixOSContour1Inds but for MovOSContour1Pts.
  
  'MovOSContour2Inds'    - Like FixOSContour1Inds but for MovOSContour2Pts.
  
  'MovInterpContourInds' - Like FixOSContour1Inds but for MovInterpContourPts.
  
  'MovOSContour1Pts'     - A list of a list of [x, y, z] coordinates for each
                           point in FixOSContour1Pts transformed to the Moving
                           image domain.
                           
  'MovOSContour2Pts'     - Like MovOSContour1Pts but for the transformation of
                           FixOSContour2Pts.
                           
  'MovInterpContourPts'  - Like MovOSContour1Pts but for the transformation of
                           FixInterpContourPts.
                           


Note:
*****
    PointData is a flat dictionary, that contains contour data on a 
    point-by-point basis.
    
    InContourData is a multi-list dictionary, that contains contour data 
    arranged on a slice-by-slice and contour-by-contour basis for all slices
    within the input DICOM image (including slices that do not contain any
    contour data).
    
    Many of the functions below were adapted from a JavaScript conversion of
    Matt Orton's Matlab code:
    
    https://bitbucket.org/icrimaginginformatics/ohif-viewer-xnat/src/master/Packages/icr-peppermint-tools/client/lib/util/freehandInterpolate/
    
"""





def GetIndsOfSlicesWithContours(ContourData, ContourType):
    """
    Get the indices of the rows (slices) within FixContourData that contain a 
    single list of contour points (i.e. only one contour per slice).
    
    Inputs:
        ContourData        - Dictionary containing the key 'PointPCS', whose
                             values are a list of contour points [x, y, z] 
                             arranged along rows for each DICOM slice, and 
                             along columns for different contours (i.e. a list 
                             of lists of lists). Slices without contours have 
                             an empty list ([]).
                             
        ContourType        - Integer that denotes what type of contour to 
                             include in the list SlicesWithContours.
                             Acceptable values are:
                                 1  -> Original contours only
                                 2  -> Interpolated contours only
                                 3  -> Interpolated, transformed contours that
                                       intersect the Moving image planes only
                                 -1 -> Any contours
                      
    Returns:
        SlicesWithContours - List of integers of the indices (row numbers) in
                             FixContourData that contains a single list of 
                             contour points.
    """
    
    # Initialise list of indices of the slices that contain contours:
    SlicesWithContours = []
    
    for i in range(len(ContourData['ContourType'])):
        if ContourType == -1 and ContourData['ContourType'][i] > 0:
            SlicesWithContours.append(i)
            
        else:
            if ContourData['ContourType'][i] == ContourType:
                SlicesWithContours.append(i)
        
    
    if SlicesWithContours:
        # Reduce to list of unique indices:
        SlicesWithContours = list(set(SlicesWithContours))
        
    else:
        # Print warning to console if SlicesWithContours = []:
        print(f'No slices in ContourData found with ContourType = {ContourType}')
        
            
    return SlicesWithContours





def GetBoundingSliceInds(FixContourData, PointData, InterpSliceInd):
    """ 
    Get the slice indices (row numbers) in FixContourData of the nearest slices 
    that bound the slice to be interpolated.
    
    Inputs:
        FixContourData      - Dictionary containing the key 'PointPCS', 
                              whose values are a list of contour points 
                              [x, y, z] arranged along rows for each DICOM 
                              slice, and along columns for different contours 
                              (i.e. a list of lists of lists). 
                              Slices without contours have empty list ([]).
                              
        PointData          - Dictionary containing the key 'InSliceNo' whose
                             values are a list of integers of the indices of
                             the DICOM slices that contain contour data.
                             
        InterpSliceInd     - The slice index within the DICOM image where a 
                             contour is to be interpolated; May be an integer  
                             (i.e. interpolation at an imaging plane) or a  
                             float (i.e. interpolation between imaging planes).
                      
    Returns:
        BoundingSliceInds  - List of integers of the indices (row numbers) in
                             FixContourData that contains a single list of 
                             contour points.
    
        
    Notes:
        
        1. The JavaScript function _getBoundingPair checks if a contour is an 
        interpolated contour, and if so, it is excluded from the list of 
        possible contours to interpolate between (hence the index of an 
        interpolated contour cannot be included in ContourPair). 
        This has not yet been implemented in GetBoundingSliceInds.
        
        2. Like _getBoundingPair, GetBoundingSliceInds will exclude any slices 
        that have multiple contours.  A more robust version will be required to 
        interpolate contours on slices with multiple contours present.
    """
    
    # Get the indices of the slices that contain contours:
    #SlicesWithContours = GetIndsOfSlicesWithContours(ContourData=ContourData)
    SlicesWithContours = list(set(PointData['InSliceNo']))
    
    #print('SlicesWithContours =', SlicesWithContours)
    
    
    # Get extent of region of contours by slice number (equivalent to function 
    # _getExtentOfRegion in generateInterpolationData.js):
    #Extent = [SlicesWithContours[0], SlicesWithContours[-1]]
    Extent = [min(SlicesWithContours), max(SlicesWithContours)]
    
    #print('Extent =', Extent)

    # Initialise BoundingSliceInds:
    BoundingSliceInds = []
    
    # Initialise booleans:
    FoundLowerSlice = False
    FoundUpperSlice = False
    
    # Cannot interpolate a slice whose index is <= Extent[0], or >= Extent[1]:
    if InterpSliceInd <= Extent[0] or InterpSliceInd >= Extent[1]:
        print(f'Cannot interpolate slice number {InterpSliceInd}')
        return []
    
    else:
        #print(f'Looking for the nearest lower slice to {InterpSliceInd}...')
        
        # Iterate backwards from (InterpSliceInd - 1) to SlicesWithContours[0] 
        # to determine the nearest lowest slice index that contains a single 
        # contour:
        if False:
            for i in range(InterpSliceInd - 1, min(SlicesWithContours) - 1, -1):
                #print(f'\nFixContourData['InPointPCS'][{i}] =', FixContourData['InPointPCS'][i])
                #print(f'\nlen(FixContourData['InPointPCS'][{i}]) = {len(FixContourData['InPointPCS'][i])}')
                
                if len(FixContourData['PointPCS'][i]) == 1:
                    #print(f'\nlen(FixContourData['PointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundLowerSlice = True
                    
                    break
            
            
        """ Allow for the possibility of InterpInd to be a fraction (not an
        integer) by taking floor of InterpInd using (InterpInd//1): """
        FloorInterpInd = int(InterpSliceInd//1)
        
        for i in range(FloorInterpInd, min(SlicesWithContours) - 1, -1):
            if len(FixContourData['PointPCS'][i]) == 1:
                # Append i to BoundingSliceInds if i != InterpSliceInd (which 
                # could arise if InterpSliceInd is an integer):
                if i != InterpSliceInd:
                    BoundingSliceInds.append(i)
                    
                    FoundLowerSlice = True
                    
                    break
                
        
        
                
                
        #print(f'Looking for the nearest upper slice to {InterpSliceInd}...')
        
        # Iterate forwards from (InterpSliceInd + 1) to SlicesWithContours[-1] 
        # to determine the nearest upper slice index that contains a single 
        # contour:
        if False:
            for i in range(InterpSliceInd + 1, max(SlicesWithContours) + 1):
                #print(f'\nFixContourData['PointPCS'][{i}] =', FixContourData['PointPCS'][i])
                #print(f'\nlen(FixContourData['PointPCS'][{i}]) = {len(FixContourData['PointPCS'][i])}')
                #print(f'Appending {i} to BoundingSliceInds')
                
                if len(FixContourData['PointPCS'][i]) == 1:
                    #print(f'\nlen(FixContourData['PointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundUpperSlice = True
                    
                    break
                
        """ Allow for the possibility of InterpSliceInd to be a fraction (not 
        an integer) by taking ceiling of InterpSliceInd using 
        # -(-InterpSliceInd//1): """
        CeilInterpInd = int(-(-InterpSliceInd//1))
        
        for i in range(CeilInterpInd, max(SlicesWithContours) + 1):
            if len(FixContourData['PointPCS'][i]) == 1:
                 # Append i to BoundingSliceInds if i != InterpSliceInd (which 
                 # could arise if InterpInd is an integer):
                if i != InterpSliceInd:
                    BoundingSliceInds.append(i)
                    
                    FoundUpperSlice = True
                    
                    break
        
        # Deal with case where no lower, or upper, or both upper nor lower 
        # slice index can be used for interpolation:
        if not (FoundLowerSlice and FoundUpperSlice):
            if not FoundLowerSlice and FoundUpperSlice:
                print('Did not find suitable lower slice to',
                      f'interpolate slice number {InterpSliceInd}.')
                
            if not FoundUpperSlice and FoundLowerSlice:
                print('Did not find suitable upper slice to',
                      f'interpolate slice number {InterpSliceInd}.')
                
            if not FoundLowerSlice and not FoundUpperSlice:
                print('Did not find suitable lower or upper slice',
                      f'to interpolate slice number {InterpSliceInd}.')
            
            return []
        
        else:
            return BoundingSliceInds 
    




def GenerateClosedContour(Contour):
    """
    Close an open contour by appending the first point to the end of the list.
    
    Inputs:
        Contour - An open list of contour points.
        
    Returns:
        Contour - A closed list of contour points.
    """
    
    # Import packages:
    #import copy
    
    #Contour = copy.deepcopy(Points)
    
    # Append the first point to create a closed contour:
    #Contour.append(Points[0])
    Contour.append(Contour[0])
    
    return Contour


def ReverseIfAntiClockwise(Contour):
    """
    Ensures that a contour's points are arranged in a clockwise order.
    
    Inputs:
        Contour - A list of contour points.
        
    Returns:
        Contour - A list of contour points whose points have been reversed if
                  they were counter-clockwise, or the original list if not.
    """
    
    
    P = len(Contour)
    
    # Get mean x position of all points:
    #MeanX = sum([Contour[i][0] for i in range(P)]) / P # correct but might not
    # be clear outside Python
    
    MeanX = 0
    
    for i in range(P):
        MeanX += Contour[i][0] / P
        
    
    # Check the direction of the points:
    CheckSum = 0
    j = 1
    k = 2
    
    for i in range(P):
        CheckSum += (Contour[j][0] - MeanX) * (Contour[k][1] - Contour[i][1])
        
        j += 1
        k += 1
        
        if j >= P:
            j = 0
        if k >= P:
            k = 0
            
            
    if CheckSum > 0:
        Contour.reverse()
        
        print('Contour was reversed.')
        
    return Contour
            

""" Generate list of cumulative perimeters of a contour at each node """

def GetCumulativePerimeters(Contour):
    """
    Generate a list of cumulative "perimeters" of a contour.
    
    Inputs:
        Contour  - A list of contour points.
        
    Returns:
        CumPerim - A list of cumulative "perimeters" of the contour.
    """
    
    #print(f'len(Contour) = {len(Contour)}')
    
    # Initialise the cumulative perimeter and list of cumulative perimeters:
    CumP = 0
    CumPerim = [0]
    
    # Loop through each point:
    for i in range(1, len(Contour)):
        # The segment length between every two points:
        SegmentL = (\
                   (Contour[i-1][0] - Contour[i][0])**2 \
                   + (Contour[i-1][1] - Contour[i][1])**2 \
                   + (Contour[i-1][2] - Contour[i][2])**2 \
                   )**(1/2)
        
        CumP = CumP + SegmentL
              
        CumPerim.append(CumP + SegmentL)
        
    return CumPerim





def NormaliseCumulativePerimeters(CumPerimeters): 
    """
    Normalised list of cumulative contour "perimeters".
    
    Inputs:
        CumPerimeters - List of cumulative contour "perimeters".
    
    Returns:
        CumPerimNorm  - List of normalised cumulative contour "perimeters".
    """
    
    # List comprehension is compact but not reader-friendly outside Python:
    #return [CumPerimeters[i]/CumPerimeters[-1] for i in range(len(CumPerimeters))]
    
    # Initialise the list:
    CumPerimNorm = []
    
    for i in range(len(CumPerimeters)):
        CumPerimNorm.append(CumPerimeters[i]/CumPerimeters[-1])
        
    return CumPerimNorm



def GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm, NumNodesToAdd):
    """ 
    Super-sample a list of normalised cumulative "perimeters".
    
    Inputs:
        CumPerimetersNorm - List of normalised cumulative contour "perimeters".
        
        NumNodesToAdd     - List of integers whose values are the number of 
                            nodes to be added within each segment along the
                            list of normalised cumulative contour perimeters.
        
    Returns:
        CumSSPerimNorm    - List of super-sampled normalised cumulative contour
                            perimeters.    
    
    Notes: 
        
        1. The first (NumNodesToAdd - 2) items in CumSSPerimNorm are the 
        cumulative contour perimeters that ensure that a minimum desired inter-
        node spacing is satisfied.  The -2 is due to skipping of the 0 and 1
        since they are already covered by the first and final items in the 
        input list CumPerimetersNorm.
        
        The next len(CumPerimetersNorm) items in CumSSPerimNorm are the items
        in the input list CumPerimetersNorm (concatenated to the end of the 
        list CumSSPerimNorm following the initial step described above).
        
        2. This function is called _getInterpolatedPerim in the JavaScript code. 
    """
    
    # The inter-node spacing is:
    dP = 1 / (NumNodesToAdd - 1)
    
    """ Exclude the first (0) and last item (1) normalised cumulative perimeter
    since they are already covered by the original normalised cumulative 
    perimeters (i.e. iterate from 1 to NumNodesToAdd - 1)
    """
    # Initialise the normalised cumulative super-sampled perimeters:
    CumSSPerimNorm = [dP]
    
    # First append the cumulative perimeters for the nodes that ensure the
    # minimum inter-node spacing required:
    #for i in range(1, NumNodesToAdd - 1):
    """ 
    Not sure why the upper limit needs to be (NumNodesToAdd - 2) instead of
    (NumNodesToAdd - 1), since the lower limit of 1 (rather than 0) already
    accounts for one, and the -1 accounts for the other. The JavaScript code 
    also has -2 so it's probably correct.
    """
    for i in range(1, NumNodesToAdd - 2):
        CumSSPerimNorm.append(CumSSPerimNorm[-1] + dP)
        
    # Next concatenate the normalised cumulative perimeters:
    CumSSPerimNorm += CumPerimetersNorm
    
    return CumSSPerimNorm
    


def GetIsOrigNode(Contour, NumNodesToAdd):
    """
    Generate an array of booleans that indicate whether the nodes are original
    or interpolated (i.e. super-sampled).
    
    Inputs:
        Contour       - A list of contour points.
        
        NumNodesToAdd - List of integers whose values are the number of 
                        nodes to be added within each segment along the list
                        of normalised cumulative contour perimeters.
                        
    Returns:
        IsOrigNode    - List of booleans whose values denote whether the 
                        corresponding contour point is an original node (True)
                        or an added node (False).
    
    
    Note:
        
        This is called _getIndicatorArray in the JavaScript code.
    """
    
    IsOrigNode = []
    
    # The first NumNodesToAdd - 2 nodes were additional nodes added to achieve 
    # the required minimum inter-node spacing:
    for i in range(NumNodesToAdd - 2):
        IsOrigNode.append(False)
        
    # The next len(Contour) nodes are original nodes from Contour:
    for i in range(len(Contour)):
        IsOrigNode.append(True)
        
        
    return IsOrigNode
    
    

def GetNodesToAddPerSegment(CumSSPerimNorm, IsOrigNode):
    """
    Generate list of the number of nodes to be added for each segment of a
    contour.
    
    Inputs:
        CumSSPerimNorm       - List of super-sampled normalised cumulative 
                               contour perimeters.
        
        IsOrigNode           - List of booleans whose values denote whether the 
                               corresponding contour point is an original node
                               (True) or an added node (False).
        
    Returns:
        Inds                 - List of indices that sort CumSSPerimNorm in
                               ascending order.
                               
        IsOrigNodeSorted     - List of IsOrigNode sorted by Inds.
        
        IndsOrigNodes        - List of indices of original nodes in 
                               CumSSPerimNorm sorted by Inds.
                               
        NodesToAddPerSegment - A list of the integer number of nodes to be
                               added to each original node/point.
    """
    
    # Get the indices that sort the normalised cumulative super-sampled 
    # perimeters:
    Inds = [i[0] for i in sorted(enumerate(CumSSPerimNorm), key=lambda x:x[1])]
    
    # Order the boolean list:
    IsOrigNodeSorted = []
    
    for ind in Inds:
        IsOrigNodeSorted.append(IsOrigNode[ind])
    
    # Get ordered indices of original nodes:
    """
    Seems this isn't needed (I've likely misunderstood the JavaScript code).
    """
    #OrigIndsSorted = []
    #
    #for i in range(len(Inds)):
    #    if IsOrigNodeSorted[i]:
    #        OrigIndsSorted.append(Inds[i])
            
    # Get the indicies of the original nodes in IsOrigNodeSorted:
    IndsOrigNodes = []
    
    for i in range(len(IsOrigNodeSorted)):
        if IsOrigNodeSorted[i]:
            IndsOrigNodes.append(i)
        

    # Order the normalised cumulative super-sampled perimeters:
    """ This is not in the JavaScript code """
    #CumSSPerimNormSorted = [CumSSPerimNorm[i] for i in Inds]

    # Get the list of normailsed cumulative super-sampled perimeters for the
    # original nodes:
    """ This is not in the JavaScript code """
    #CumSSPerimNormSortedOrig = [CumSSPerimNormSorted[i] for i in OrigIndsSorted]

    # The segment lengths are:
    """ This is not in the JavaScript code """
    #SegmentL = [CumSSPerimNormSortedOrig[i+1] - CumSSPerimNormSortedOrig[i] for i in range(len(CumSSPerimNormSortedOrig) - 1)]
    
    # Initialise the list of the number of nodes to add to each segment:
    NodesToAddPerSegment = []
    
    #for i in range(len(OrigIndsSorted) - 1):
        #NodesToAddPerSegment.append(OrigIndsSorted[i + 1] - OrigIndsSorted[i])
        
    for i in range(len(IndsOrigNodes) - 1):
        NodesToAddPerSegment.append(IndsOrigNodes[i + 1] - IndsOrigNodes[i])
    
        
    #return NodesToAddPerSegment
    #return Inds, IsOrigNodeSorted, OrigIndsSorted, NodesToAddPerSegment
    return Inds, IsOrigNodeSorted, IndsOrigNodes, NodesToAddPerSegment



    
    
def SuperSampleContour(Contour, NodesToAddPerSegment):
    """
    Super-sample a contour.
    
    Inputs:
        Contour              - A list of contours points.
        
        NodesToAddPerSegment - A list of the integer number of nodes to be 
                               added to each original node/point in Contour.
        
    Returns:
        SSContour            - A list of super-sampled (with added nodes) 
                               contour points.
        
        SSIsOrigNode         - A list of booleans whose values denote whether  
                               the corresponding contour point in SSContour is 
                               an original node (True) or an added node (False).
    """
    
    # Initialise list of super-sampled contour points and the list of booleans 
    # that indicate if the point/node in the super-sampled list is original 
    # (False if the node was interpolated):
    SSContour = []
    SSIsOrigNode = []
    
    #print(f'len(Contour)         = {len(Contour)}')
    #print(f'len(NodesToAddPerSegment) = {len(NodesToAddPerSegment)}')
    
    # Iterate through each contour point without returning to the 0th point 
    # (which is the final node in the list), i.e. the super-sampled contour 
    # will not be closed:
    for i in range(len(Contour) - 1):
        #print(i)
        
        # Add the original node:
        SSContour.append(Contour[i])
        SSIsOrigNode.append(True)
        
        # Add the linearly interpolated nodes.  Get the spacings along x and y
        # for this segment:
        dx = (Contour[i + 1][0] - Contour[i][0]) / (NodesToAddPerSegment[i] + 1)
        dy = (Contour[i + 1][1] - Contour[i][1]) / (NodesToAddPerSegment[i] + 1)
        
        # Add dx and dy to the x and y components of the last item in SSContour
        # (keep the z component unchanged):
        """
        Not sure why j is iterated to NodesToAddPerSegment[i] - 1, instead of
        NodesToAddPerSegment[i] so deviating from the JavaScript code:
        """
        for j in range(NodesToAddPerSegment[i] - 1):
        #for j in range(NodesToAddPerSegment[i]):
            SSContour.append([
                              SSContour[-1][0] + dx, 
                              SSContour[-1][1] + dy, 
                              SSContour[-1][2]
                              ])
            SSIsOrigNode.append(False)
            
    return SSContour, SSIsOrigNode





#def IntegrateLineLengths(Contour1, Contour2, StartInd2):
#    
#    # Initialise the cummulative line lengths:
#    LineLength = 0
#    
#    for i in range(len(Contour1)):
#        Point1 = Contour1[i]
#        
#        for j in range(len(Contour2)):
#            # The index of the point for Contour2 starts
#            # at StartInd2 and increments to N, then is 
#            # "wrapped" back to 0 and so on until N 
#            # iterations:
#            if StartInd2 + j < len(Contour2):
#                Ind2 = StartInd2 + j
#            else:
#                Ind2 = len(Contour2) - (StartInd2 + j)
#                
#            Point2 = Contour2[Ind2]
#            
#            # The vector length between Point1 and Point2:
#            VectorL = ((Point1[0] - Point2[0])**2 \
#                       + (Point1[1] - Point2[1])**2 \
#                       + (Point1[2] - Point2[2])**2)**(1/2)
#            
#            # Add the line length between Point1 and Point2:
#            LineLength = LineLength + VectorL
#            
#            
#    return LineLength



def IntegrateLineLengths(Contour1, Contour2):
    """ 
    Integrate the line lengths that connect every point in Contour1 with the
    corresponding point in Contour2.  
    
    Inputs:
        Contour1   - A list of contour points.
        
        Contour2   - A list of contour points.
        
    Returns:
        LineLength - The sum of all line lengths that join the points in 
                     Contour1 and Contour2.
                     
    Note: 
        Contour1 and Contour2 must have equal lengths. 
        
    """
    
    # Initialise the cummulative line lengths:
    LineLength = 0
    
    for i in range(len(Contour1)):
            # The vector length between Point1 and Point2:
            VectorL = (
                       (Contour1[i][0] - Contour2[i][0])**2 \
                       + (Contour1[i][1] - Contour2[i][1])**2 \
                       + (Contour1[i][2] - Contour2[i][2])**2 \
                       )**(1/2)
            
            LineLength += VectorL
            
            
    return LineLength





def FindIndForMinCumLength(Contour1, Contour2):
    """ 
    Determine the index for Contour2 that would minimise the integrated line
    lengths obtained by joining the points in Contour1 to the corresponding 
    points in Contour2 with the index shift applied to Contour2.
    
    Inputs:
        Contour1 - A list of contour points.
        
        Contour2 - A list of contour points.
        
    Returns:
        Shift    - The starting index for Contour2 that minimises the 
                   integrated line lengths.
    
    
    Note:
        The Shift index is determined by looping through each of N possible
        combinations of Shift index (for N points in Contour2) and applying the 
        function IntegrateLineLengths() with those shifts. The starting index 
        that yields the minimum cumulative line length is equivalent to 
        minimising the surface area (i.e. catenoid) that would be formed by 
        joining each i^th point in Contour1 to each (i + Shift)^th point in 
        Contour 2.

        https://en.wikipedia.org/wiki/Minimal_surface_of_revolution

    """
    
    # Create list to store N cumulative line lengths for all N combinations of 
    # starting index for Contour2:
    LineLengths = []

    # Loop through all N possible starting indeces for Contour2:
    for i in range(len(Contour2)):
        #CumLineLengths.append(IntegrateLineLengths(Contour1=Contour1, 
        #                                           Contour2=Contour2, 
        #                                           StartInd2=i))
        
        # Apply shift to Contour2:
        Contour2Shifted = ShiftContourPoints(Contour=Contour2, Shift=i)
        
        # Get the cumumlative line length for this shift:
        LineLength = IntegrateLineLengths(Contour1=Contour1, 
                                          Contour2=Contour2Shifted)
        
        LineLengths.append(LineLength)

    
    # The index that yields the minimum cumulative line length:
    Shift = LineLengths.index(min(LineLengths))
    
    return Shift






def ShiftContourPoints(Contour, Shift):
    """
    Generate a new list of contour points by shifting the points in the 
    original list of contour points by an integer number.
    
    Inputs:
        Contour        - A list of contour points.
        
        Shift          - A integer value.
        
    Returns:
        ShiftedContour - A list of the contour points in Contour shifted by
                         Shift integer positions.
    """
    
    # Initialise the shifted contour:
    ShiftedContour = []
    
    for i in range(len(Contour)):
        j = i + Shift
        # Wrap j if (i + Shift) >= len(Contour):
        if i + Shift >= len(Contour):
            j = j - (len(Contour) - 1) 

        ShiftedContour.append(Contour[j])
        
    return ShiftedContour






def ReduceNodesOfContours(SuperSampledContour1, SuperSampledContour2, 
                          SuperSampledIsOrigNode1, SuperSampledIsOrigNode2):
    """
    Reduce the number of nodes in a list of super-sampled contour points to a
    "over-sampled" number of nodes.
    
    Inputs:
        SuperSampledContour1    - A list of super-sampled contour points.
        
        SuperSampledContour2    - A list of super-sampled contour points.
        
        SuperSampledIsOrigNode1 - A list of booleans whose values denote   
                                  whether the corresponding contour point in 
                                  SuperSampledContour1 is an original node 
                                  (True) or an added node (False).
        
        SuperSampledIsOrigNode2 - A list of booleans whose values denote   
                                  whether the corresponding contour point in 
                                  SuperSampledContour2 is an original node 
                                  (True) or an added node (False).
        
    Returns:
        OverSampledContour1     - A list of over-sampled contour points after 
                                  removing of any node in SuperSampledContour1
                                  if both that node and the corresponding node
                                  in SuperSampledContour2 was added (i.e. not
                                  original).
        
        OverSampledContour2     - A list of over-sampled contour points after 
                                  removing of any node in SuperSampledContour2
                                  if both that node and the corresponding node
                                  in SuperSampledContour1 was added (i.e. not
                                  original).
        
        OverSampledIsOrigNode1  - A list of booleans whose values denote   
                                  whether the corresponding contour point in 
                                  OverSampledContour1 is an original node 
                                  (True) or an added node (False).
        
        OverSampledIsOrigNode2  - A list of booleans whose values denote   
                                  whether the corresponding contour point in 
                                  OverSampledContour2 is an original node 
                                  (True) or an added node (False).
    """
    
    # Initialise the list of points in the over-sampled contours with some
    # nodes removed, and the corresponding IsOrigNode lists:
    OverSampledContour1 = []
    OverSampledContour2 = []
    OverSampledIsOrigNode1 = []
    OverSampledIsOrigNode2 = []
    
    for i in range(len(SuperSampledContour1)):
        if SuperSampledIsOrigNode1[i] or SuperSampledIsOrigNode2[i]:
            OverSampledContour1.append(SuperSampledContour1[i])
            OverSampledContour2.append(SuperSampledContour2[i])
            
            OverSampledIsOrigNode1.append(SuperSampledIsOrigNode1[i])
            OverSampledIsOrigNode2.append(SuperSampledIsOrigNode2[i])
            
    return OverSampledContour1, OverSampledContour2, \
           OverSampledIsOrigNode1, OverSampledIsOrigNode2
           
           


def InterpolateBetweenContours(Contour1, Contour2, IsOrigNode1, IsOrigNode2,
                               FracInterpSliceInd, Contour1HasMorePts,
                               InterpolateAllPts):
    """ 
    Generate a list of intepolated contour points between two contours.
    
    Inputs:
        Contour1            - A list of contour points.
        
        Contour2            - A list of contour points.
        
        IsOrigNode1         - A list of booleans whose values denote whether
                              the corresponding contour point in Contour1 is an 
                              original node (True) or an added node (False).
                              
        IsOrigNode2         - A list of booleans whose values denote whether 
                              the corresponding contour point in Contour1 is an 
                              original node (True) or an added node (False).
                              
        FracInterpSliceInd  - A integer or float that denotes the fractional 
                              distance of the interpolated contour's plane 
                              within the region bounded by the planes for 
                              Contour1 and Contour2.
                             
        Contour1HasMorePts  - A boolean that denotes whether or not Contour1
                              has more points than Contour2.
                              
        InterpolateAllPts   - A boolean that determines whether or not 
                              interpolation will be performed only between a
                              pair of contour points for which the corresponding
                              node in the longer of the original contours was
                              an original node.  See Note 4.
    
    Returns:
        InterpolatedContour - A list of interpolated contour points.
    
    Notes:
        1. Contour1 and Contour2 would normally be over-sampled contours (i.e. 
        not the original lists of contours).
        
        2. Contour1 and Contour2 must have equal lengths. 
        
        3. The fraction FracInterpSliceInd would be 0 if the interpolated 
        contour were to coincide with the plane containing Contour1, and 1 if 
        it were to coincide with the plane containing Contour2.
        
        4. In the JavaScript code the interpolation is performed between the 
        pair of contour points only if the corresponding contour point in the 
        longer list of contour points was an original node. This results in an 
        interpolated contour with the same number of points in it as the number
        of True values in the longer of IsOrigNode1 and IsOrigNode2, and hence,
        the length of interpolated contour points is shorter than the length of
        Contour1 or Contour2.
        
        In this function I've allowed for either using the truncation used in
        the JavaScript code (InterpolateAllPts = False) or to interpolate all 
        pairs of points, so that the length of Contour1, Contour2 and 
        InterpolatedContour are equal (InterpolateAllPts = True).
    """
    
    N1 = len(Contour1)
    N2 = len(Contour2)
    
    if not N1 == N2:
        print(f'Contour1 (N = {N1}) and Contour2 (N = {N2}) do not have',
              'equal lengths.')
        
        return
    
    # Initialise the interpolated contour:
    InterpolatedContour = []
    
    # If the original Contour1 has more points than the original Contour2..
    if Contour1HasMorePts:
        IsOrigNode = IsOrigNode1
    else:
        IsOrigNode = IsOrigNode2
        
    for i in range(len(Contour1)):
        if InterpolateAllPts:
            # Interpolate all points:
            InterpX = (1 - FracInterpSliceInd) * Contour1[i][0] \
                       + FracInterpSliceInd * Contour2[i][0]
                
            InterpY = (1 - FracInterpSliceInd) * Contour1[i][1] \
                      + FracInterpSliceInd * Contour2[i][1]
                      
            InterpZ = (1 - FracInterpSliceInd) * Contour1[i][2] \
                      + FracInterpSliceInd * Contour2[i][2]
        
            InterpolatedContour.append([InterpX, InterpY, InterpZ])
            
            
        else:
            # Only interpolate between points for which IsOrigNode is True:
            if IsOrigNode[i]:
                InterpX = (1 - FracInterpSliceInd) * Contour1[i][0] \
                          + FracInterpSliceInd * Contour2[i][0]
                
                InterpY = (1 - FracInterpSliceInd) * Contour1[i][1] \
                          + FracInterpSliceInd * Contour2[i][1]
                          
                InterpZ = (1 - FracInterpSliceInd) * Contour1[i][2] \
                          + FracInterpSliceInd * Contour2[i][2]
                          
                InterpolatedContour.append([InterpX, InterpY, InterpZ])
            
            
    return InterpolatedContour






def InterpolateContours(FixContourData, PointData, InterpSliceInd, dP, 
                        InterpolateAllPts):
    """
    Main function for contour interpolation.
    
    Inputs:
        FixContourData    - Dictionary containing the key 'PointPCS', whose
                            values are a list of contour points [x, y, z] 
                            arranged along rows for each DICOM slice, and along
                            columns for different contours (i.e. a list of  
                            lists of lists). Slices without contours have empty 
                            list ([]).
                              
        PointData         - Dictionary containing the key 'InSliceNo' whose 
                            values are a list of integers of the indices of the
                            DICOM slices that contain contour data.
                             
        InterpSliceInd    - The index within the DICOM image where a contour is
                            to be interpolated; May be an integer (i.e. 
                            interpolation at an imaging plane) or a float (i.e.  
                            interpolation between imaging planes).
                             
        dP                - A float denoting the desired inter-node spacing to 
                            use when super-sampling the contours used for 
                            interpolation.
                         
        InterpolateAllPts - A boolean that determines whether or not 
                            interpolation will be performed only between a pair
                            of contour points for which the corresponding node
                            in the longer of the original contours was an
                            original node.  See Note 4 in the function 
                            InterpolateBetweenContours.
        
    Returns:
        InterpContour     - A list of interpolated contour points.
     
    """
    
    import copy
    
    # Get bounding slice indices:
    BoundingSliceInds = GetBoundingSliceInds(FixContourData=FixContourData,
                                             PointData=PointData,
                                             InterpSliceInd=InterpSliceInd)
    
    
    # Escape if None was returned for BoundingSliceInds (i.e. if two slice 
    # indices were not found):
    if not BoundingSliceInds:
        return
    
    print(f'Interpolating between slice {BoundingSliceInds[0]} and',
          f'{BoundingSliceInds[1]} for slice at {InterpSliceInd}')
    
    # Define the contours:
    """
    Note: 
    
    Need to add [0] to the end of FixContourData['PointPCS'][BoundingSliceInds[0]] 
    since the contour points are a list within the list for each slice.
    """
    Contour1 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[0]][0])
    Contour2 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[1]][0])
    
    # Close the contours:
    Contour1.append(Contour1[0])
    Contour2.append(Contour2[0])
    
    # Get the segment lengths (surplus to interpolation algorithm):
    #SegmentLengths1 = GetSegmentLengths(Contour1)
    #SegmentLengths2 = GetSegmentLengths(Contour2)
    
    # Ensure that the points in both contours run clockwise:
    Contour1 = ReverseIfAntiClockwise(Contour=Contour1)
    Contour2 = ReverseIfAntiClockwise(Contour=Contour2)
    
    # Get the cumulative perimeters of the two bounding contours:
    CumPerim1 = GetCumulativePerimeters(Contour=Contour1)
    CumPerim2 = GetCumulativePerimeters(Contour=Contour2)
    
    # Normalise the cumulative perimeters:
    CumPerim1Norm = NormaliseCumulativePerimeters(CumPerimeters=CumPerim1)
    CumPerim2Norm = NormaliseCumulativePerimeters(CumPerimeters=CumPerim2)
    
    # The number of nodes required for Contour1 and Contour2 to achieve the 
    # minimum inter-node spacing of dP:
    """ 
    Note:
        a//b yields floor(a/b)
        
        - (-a//b) yields ceil(a/b)
    """
    MinNumNodes1 = round(- (- CumPerim1[-1] // dP))
    MinNumNodes2 = round(- (- CumPerim2[-1] // dP))
    
    # The minimum number of nodes that will need to be added:
    MinNumNodes = max(MinNumNodes1, MinNumNodes2)
    
    # The number of nodes that will have to be added to Contour1 and Contour2 
    # are MinNumNodes plus the number of original nodes in Contour2 and Contour1
    # respectively (i.e. the nodes not common to either)?:
    NumNodesToAdd1 = MinNumNodes + len(Contour2)
    NumNodesToAdd2 = MinNumNodes + len(Contour1)
    
    # Get the cumulative super-sampled perimeters:
    CumSSPerim1 = GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm=CumPerim1Norm,
                                                          NumNodesToAdd=NumNodesToAdd1)
    CumSSPerim2 = GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm=CumPerim2Norm,
                                                          NumNodesToAdd=NumNodesToAdd2)
    
    # Get list of booleans for original (True) and added (False) nodes:
    IsOrigNode1 = GetIsOrigNode(Contour=Contour1, NumNodesToAdd=NumNodesToAdd1)
    IsOrigNode2 = GetIsOrigNode(Contour=Contour2, NumNodesToAdd=NumNodesToAdd2)
    
    # Get the nodes per segment:
    Inds1,\
    IsOrigNodeSorted1,\
    IndsOrigNodes1,\
    NodesToAddPerSegment1 = GetNodesToAddPerSegment(CumSSPerimNorm=CumSSPerim1, 
                                                    IsOrigNode=IsOrigNode1)
    
    Inds2,\
    IsOrigNodeSorted2,\
    IndsOrigNodes2,\
    NodesToAddPerSegment2 = GetNodesToAddPerSegment(CumSSPerimNorm=CumSSPerim2, 
                                                    IsOrigNode=IsOrigNode2)
    
    # Super-sample the contours:
    SSContour1,\
    SSIsOrigNode1 = SuperSampleContour(Contour=Contour1, 
                                       NodesToAddPerSegment=NodesToAddPerSegment1)
    
    SSContour2,\
    SSIsOrigNode2 = SuperSampleContour(Contour=Contour2, 
                                       NodesToAddPerSegment=NodesToAddPerSegment2)
    
    # Find the starting index for SSContour2 that minimises the integrated line 
    # lengths (i.e. minimal area of the surfaces created by linking points in 
    # SSContour1 and SSContour2 N different ways for all N points in 
    # SSContour2):
    IndMinCumL = FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)
    
    # Shift the points in SSContour2 (and also SSIsOrigNode2):
    SSContour2Shifted = ShiftContourPoints(Contour=SSContour2, Shift=IndMinCumL)
    SSIsOrigNode2Shifted = ShiftContourPoints(Contour=SSIsOrigNode2, Shift=IndMinCumL)
    
    # Remove added nodes uncommon to both original contours:
    OSContour1,\
    OSContour2,\
    OSIsOrigNode1,\
    OSIsOrigNode2 = ReduceNodesOfContours(SuperSampledContour1=SSContour1,
                                          SuperSampledContour2=SSContour2Shifted,
                                          SuperSampledIsOrigNode1=SSIsOrigNode1, 
                                          SuperSampledIsOrigNode2=SSIsOrigNode2Shifted)
    
    # Define Contour1HasMorePts:
    if len(Contour1) > len(Contour2):
        Contour1HasMorePts = True
    else:
        Contour1HasMorePts = False
    
    
    """ 
    Note:  What if len(Contour1) = len(Contour2)?...
    """
    
    # Define the fractional interpolated slice number:
    FracInterpInd = (InterpSliceInd - BoundingSliceInds[0]) / (BoundingSliceInds[1] - BoundingSliceInds[0])
    
    # Linearly interpolate between the two over-sampled contours:
    InterpContour = InterpolateBetweenContours(Contour1=OSContour1,
                                               Contour2=OSContour2,
                                               IsOrigNode1=OSIsOrigNode1,
                                               IsOrigNode2=OSIsOrigNode2,
                                               FracInterpSliceInd=FracInterpInd,
                                               Contour1HasMorePts=Contour1HasMorePts,
                                               InterpolateAllPts=InterpolateAllPts)
    
    
    
    #return InterpContour
    return BoundingSliceInds, OSContour1, OSIsOrigNode1,\
           OSContour2, OSIsOrigNode2, InterpContour


    
    """ Deviate from the JavaScript code - try interpolating using the 
    super-sampled contours (i.e. before reducing the num of nodes): """
    if False:
        InterpSSContour = InterpolateBetweenContours(Contour1=SSContour1,
                                                     Contour2=SSContour2,
                                                     IsOrigNode1=SSIsOrigNode1,
                                                     IsOrigNode2=SSIsOrigNode2,
                                                     FracInterpSliceInd=FracInterpInd,
                                                     Contour1HasMorePts=Contour1HasMorePts)
        
        # Define boolean list of indices to keep in the super-sampled 
        # interpolated contour:
        IsOrigNode1or2 = IsOrigNode1 or IsOrigNode2
        
        # Reduce the number of nodes of InterpSSContour:
        ReducedInterpSSContour = ReduceNodesOfContour(Contour=InterpSSContour, 
                                                      IndsToKeep=IsOrigNode1or2)
        
        return InterpContour, InterpSSContour, ReducedInterpSSContour
    
    


   



#def TransformInterpolatedContours(OSContour1, OSContour2, InterpContour, 
#                                  InterpData, FixedDicomDir, MovingDicomDir):
def TransformInterpolatedContours(InterpData, FixedDicomDir, MovingDicomDir):
    """
    Transform the interpolated contours and the over-sampled contours that were
    used to interpolate, using the transformation filter from Elastix, within
    the dictionary InterpData.
    
    Inputs:
        InterpData         - A dictionary containing interpolated data.
        
        FixedDicomDir      - Directory containing the DICOMs for the Fixed 
                             image domain.
        
        MovingDicomDir     - Directory containing the DICOMs for the Moving 
                             image domain.
        
        
    Returns:
        InterpData         - As above but with added entries for the 
                             transformed contours.
                             
        FixIm              - 3D Fixed image as a SimpleITK object
        
        MovIm              - 3D Moving image as a SimpleITK object
        
        RegIm              - 3D Registered image as a SimpleITK object
        
        ElastixImFilt      - Image transformation filter used by Elastix to
                             transform MovingImage to the Fixed image domain.
    
    """

    # Import packages and functions:
    import SimpleITK as sitk
    from GetImageAttributes import GetImageAttributes
    from RunSimpleElastixReg import RunSimpleElastixReg
    from CreateInputFileForElastix import CreateInputFileForElastix
    from TransformPoints import TransformPoints
    from ParseTransformixOutput import ParseTransformixOutput
                                                                     
    # Chose which package to use to get the Image Plane Attributes:
    #package = 'sitk'
    package = 'pydicom'

    # Get the Image Attributes for the images:
    FixOrigin, FixDirs,\
    FixSpacings, FixDims = GetImageAttributes(DicomDir=FixedDicomDir, 
                                              Package=package)
    
    MovOrigin, MovDirs,\
    MovSpacings, MovDims = GetImageAttributes(DicomDir=MovingDicomDir, 
                                              Package=package)

    # Read in the 3D stack of Fixed DICOMs:
    FixReader = sitk.ImageSeriesReader()
    FixNames = FixReader.GetGDCMSeriesFileNames(FixedDicomDir)
    FixReader.SetFileNames(FixNames)
    FixIm = FixReader.Execute()
    FixIm = sitk.Cast(FixIm, sitk.sitkFloat32)
    
    # Read in the 3D stack of Moving DICOMs:
    MovReader = sitk.ImageSeriesReader()
    MovNames = MovReader.GetGDCMSeriesFileNames(MovingDicomDir)
    MovReader.SetFileNames(MovNames)
    MovIm = MovReader.Execute()
    MovIm = sitk.Cast(MovIm, sitk.sitkFloat32)
    
    # Register MovIm to FixIm to obtain the transformation filter:
    RegIm, ElastixImFilt = RunSimpleElastixReg(FixedIm=FixIm, MovingIm=MovIm, 
                                               LogToConsole=False)

                                                         
    # Create lists of keys in InterpData whose values are to be transformed 
    # (will loop through each key):
    BaseKeys = ['OSContour1', 'OSContour2', 'InterpContour']
    
    FixIndsKeys = ['Fix' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
    FixPtsKeys = ['Fix' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]
    
    MovIndsKeys = ['Mov' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
    MovPtsKeys = ['Mov' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]
    
    # Initialise the additional keys to be added to InterpData (this isn't
    # required but ensures that the keys are added in the same order for
    # Mov as for Fix):
    InterpData.update({
                       'MovOSContour1Inds':[],
                       'MovOSContour2Inds':[],
                       'MovInterpContourInds':[],
                       'MovOSContour1Pts':[],
                       'MovOSContour2Pts':[],
                       'MovInterpContourPts':[]
                       })
    
    # Loop through each row in InterpData and transform the points in FixPtsKeys:
    for i in range(len(InterpData['InterpSliceInd'])):
        
        
        for k in range(len(FixPtsKeys)):
            Points = InterpData[FixPtsKeys[k]][i]
            
            # Create inputpoints.txt for Elastix, containing the contour points to be
            # transformed:
            CreateInputFileForElastix(Points=Points)
    
            # Transform MovingContourPts:
            TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)
    
            # Parse outputpoints.txt:
            PtNos, FixInds, FixPts_PCS,\
            FixOutInds, MovPts_PCS,\
            Def_PCS, MovInds = ParseTransformixOutput()
            
            # Initialise the keys to be added to InterpData:
            #InterpData.update({MovIndsKeys[k] : [],
            #                   MovPtsKeys[k]  : []
            #                   })
            
            # Append the transformed results to InterpData:
            InterpData[FixIndsKeys[k]].append(FixInds) # the indices of the input (fixed) points
            #InterpData[MovIndsKeys[k]].append(MovInds) # the indices of the output (moving/transformed) points
            #InterpData[MovPtsKeys[k]].append(MovPts_PCS) # the output (moving/transformed) points
            
            # Get list of keys in InterpData:
            keys = list(InterpData.keys())
            
            if MovIndsKeys[k] in keys:
                InterpData[MovIndsKeys[k]].append(MovInds) # the indices of the output (moving/transformed) points
            else:
                InterpData.update({MovIndsKeys[k] : [MovInds]}) # the indices of the output (moving/transformed) points
                
            if MovPtsKeys[k] in keys:
                InterpData[MovPtsKeys[k]].append(MovPts_PCS) # the output (moving/transformed) points
            else:
                InterpData.update({MovPtsKeys[k] : [MovPts_PCS]}) # the output (moving/transformed) points
    
    
    return InterpData, FixIm, MovIm, RegIm, ElastixImFilt

    
    

    
def GetLinePlaneIntersection(PlaneN, PlaneP, LineP0, LineP1, OnLineSegment):
    """ 
    Get the intersection between a line and a plane.
    
    Inputs:
        PlaneN        - A list of the [x, y, z] components of the plane's 
                        normal
                        
        PlaneP        - A list of the [x, y z] components of a point lying on 
                        the plane

        LineP0        - A list of the [x, y, z] components of a point lying on
                        the line
                        
        LineP1        - A list of the [x, y z] components of a point lying on 
                        the line
                        
        OnLineSegment - Boolean value that determines whether the intersection
                        point must lie on the line segment (True) or not (False)
                        
    Returns:
        IntersP       - A list of the [x, y, z] components of the point that
                        intersects the line and plane; or None
    
    Code adapted from:
    https://stackoverflow.com/questions/5666222/3d-line-plane-intersection
    """
    
    epsilon=1e-6
        
    LineV = [LineP1[0] - LineP0[0],
             LineP1[1] - LineP0[1],
             LineP1[2] - LineP0[2]
             ]
    
    PlaneN_dot_LineV = PlaneN[0] * LineV[0] + \
                       PlaneN[1] * LineV[1] + \
                       PlaneN[2] * LineV[2]
    
    if abs(PlaneN_dot_LineV) < epsilon:
        raise RuntimeError("The line is parallel to the plane")
        return None
        
    PlaneLineV = [LineP0[0] - PlaneP[0],
                  LineP0[1] - PlaneP[1],
                  LineP0[2] - PlaneP[2]
                  ]
    
    PlaneN_dot_PlaneLineV = PlaneN[0] * PlaneLineV[0] + \
                            PlaneN[1] * PlaneLineV[1] + \
                            PlaneN[2] * PlaneLineV[2]
    
    factor = - PlaneN_dot_PlaneLineV / PlaneN_dot_LineV
    
    IntersP = [LineP0[0] + factor * LineV[0],
               LineP0[1] + factor * LineV[1],
               LineP0[2] + factor * LineV[2]
               ]
    
    if OnLineSegment and ( factor < 0 ) or ( factor > 1 ):
        # IntersP does not lie on the line segment.
        return None
        
    else:
        return IntersP
            
    



def GetIntersectingContoursInMovingPlanes(InterpData, MovingDicomDir, 
                                          UseInterp, OnLineSegment):
    """
    Generate contours in the moving planes by finding the intersection of the 
    lines that join each over-sampled-to-interpolated contour #1 to each 
    imaging plane, and likewise for each interpolated-to-over-sampled contour 
    #2.
    
    Inputs:
        InterpData     - A dictionary containing interpolated data.
        
        MovingDicomDir - Directory containing the DICOMs for the Moving image
                         domain.
                         
        UseInterp      - Boolean that determines whether interpolated contours
                         will be used (True) when finding intersecting points
                         or not (False).
                         
        OnLineSegment  - Boolean value that determines whether the intersection
                         point must lie on the line segment (True) or not (False)
        
        
    Returns:
        MovContourData - A dictionary containing the contour data belonging to
                         the Moving image domain.
    """
    
    # Import packages and functions:
    from GetImageAttributes import GetImageAttributes
    from PCStoICS import PCStoICS
    
    # Chose which package to use to get the Image Plane Attributes:
    #package = 'sitk'
    package = 'pydicom'
    
    # Get the Image Attributes for the Moving image:
    MovOrigin, MovDirs,\
    MovSpacings, MovDims = GetImageAttributes(DicomDir=MovingDicomDir, 
                                              Package=package)
    
    # The imaging plane normal:
    PlaneNormal = MovDirs[6:9]
    
    # Initialise a list of contours consisting of the intersecting points 
    # for each imaging plane:
    IntersContours = [[] for i in range(MovDims[2])]
    
    
    # Loop through each set of interpolation data in InterpData:
    for i in range(len(InterpData['InterpSliceInd'])):
    
        # Loop through each point in Contour1, InterpContour and Contour2:
        for j in range(len(InterpData['MovOSContour1Pts'][i])):
    
            Point0 = InterpData['MovOSContour1Pts'][i][j]
            Point1 = InterpData['MovInterpContourPts'][i][j]
            Point2 = InterpData['MovOSContour2Pts'][i][j]
    
            # Loop through all imaging planes in Moving Image and get the 
            # intersecting points of Vector1 and Vector2 with each plane:
            for k in range(MovDims[2]):
    
                # Need a point lying on the plane - use the plane's origin:
                PlaneOrigin = [MovOrigin[0] + k*MovSpacings[2]*PlaneNormal[0], 
                               MovOrigin[1] + k*MovSpacings[2]*PlaneNormal[1], 
                               MovOrigin[2] + k*MovSpacings[2]*PlaneNormal[2]
                               ]
    
                #print(f'Length of PlaneNormal is {GetVectorLength(PlaneNormal)}')
                
                if UseInterp:
                    IntersPoint0 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                            PlaneP=PlaneOrigin,
                                                            LineP0=Point0,
                                                            LineP1=Point1,
                                                            OnLineSegment=OnLineSegment)
                
                    IntersPoint1 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                            PlaneP=PlaneOrigin,
                                                            LineP0=Point1,
                                                            LineP1=Point2,
                                                            OnLineSegment=OnLineSegment)
                    
                else:
                    IntersPoint0 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                            PlaneP=PlaneOrigin,
                                                            LineP0=Point0,
                                                            LineP1=Point2,
                                                            OnLineSegment=OnLineSegment)
                
                
                # Append IntersPoint0 and IntersPoint1 to IntersContours if 
                # the intersecting points were found: 
                if IntersPoint0:
                    IntersContours[k].append(IntersPoint0)
                if UseInterp and IntersPoint1:
                    IntersContours[k].append(IntersPoint1)
                    
                
                if False:
                    print('The intersection of Line0 joining point', Point0, 'in slice', 
                          InterpData['BoundingSliceInds'][i][0], 'and \npoint', Point1,
                          'in slice', InterpData['InterpSliceInd'][i], 'and the', 
                          f'{i}^th image plane with origin\n', PlaneOrigin, 'is', 
                          IntersPoint0, '.\n\n')
    
                    print('The intersection of Line1 joining point', Point1, 'in slice', 
                          InterpData['InterpSliceInd'][i], 'and \npoint', Point2,
                          'in slice', InterpData['BoundingSliceInds'][i][1], 'and the',  
                          f'{i}^th image plane with origin\n', PlaneOrigin, 'is', 
                          IntersPoint1, '.\n\n')
            
    """
    Note:
        IntersContours has lists of points by contour. Modify it so that each 
        slice has a list of a list of points (whereby the first list represents 
        the list of points that belong to a contour, and the list of points 
        follow), so that it's consistent with the structure of FixContourData.
    """
    
    # Initialise IntersContoursPCS:
    IntersContoursPCS = []
    
    for PtsThisSlice in IntersContours:    
        if PtsThisSlice: # i.e. if not empty
            IntersContoursPCS.append([PtsThisSlice])
        else:
            IntersContoursPCS.append([])
    
    
    
    # Create a dictionary to store contour points for Moving image (obtained by
    # transforming the input/Fixed points, interpolating and finding 
    # intersection points with the image planes in Moving image) arranged by 
    # slices and contours:
    """
    Note:
        Assign the number 3 for ContourType for all slices that have contours 
        and 0 for all slices that have no contours. 
        
        N.B. ContourType 0 -> No contours
                         1 -> Original contour
                         2 -> Interpolated contour
                         3 -> Contour obtained by transforming original or 
                              interpolated points in the Fixed image domain and 
                              finding the intersection with the Moving image 
                              planes
        
    """
    
    # Create a list of slice numbers:
    SliceNo = list(range(MovDims[2]))
    
    # Initialise ContourType:
    ContourType = []
    
    for contour in IntersContours:
        if contour:
            ContourType.append(3)
            
        else:
            ContourType.append(0)
    
    # Initialise IntersContoursICS:
    IntersContoursICS = []
    
    for i in range(len(IntersContoursPCS)):
        PtsThisSlicePCS = IntersContoursPCS[i]
        
        if PtsThisSlicePCS:
            # Convert points from PCS to ICS:
            """
            Note: There should only be one contour since the contour
            interpolation code can currently only deal with one contour
            per slice, but to future-proof this loop through all contours.
            """ 
            PtsThisSliceICS = []
            
            for c in range(len(PtsThisSlicePCS)):
                PtsThisSliceThisContourPCS = PtsThisSlicePCS[c]
                
                PtsThisSliceThisContourICS = PCStoICS(Pts_PCS=PtsThisSliceThisContourPCS, 
                                                      Origin=MovOrigin, 
                                                      Directions=MovDirs, 
                                                      Spacings=MovSpacings)
                                                      
                PtsThisSliceICS.append(PtsThisSliceThisContourICS)
                                                      
                
            IntersContoursICS.append(PtsThisSliceICS)
                                     
            
        else:
            IntersContoursICS.append([])
    
    
    MovContourData = {
                      'SliceNo'     : SliceNo, 
                      'ContourType' : ContourType,
                      'PointPCS'    : IntersContoursPCS,
                      'PointICS'    : IntersContoursICS
                       }

    return MovContourData






def CopyRois(FixedDicomDir, MovingDicomDir, FixedRoiFpath, dP, 
             InterpolateAllPts, LogToConsole, PlotResults, AnnotatePtNums, 
             ExportResults):
    """
    Copy ROIs from the Fixed to Moving image domain.  
    
    This is a main function that runs the following functions:
        1. InterpolateContours
        
        2. TransformInterpolatedContours
        
        3. GetIntersectingContoursInMovingPlanes
        
        
    Interpolation will be midway between all contours in a ROI collection for
    the Fixed image domain by running the main function InterpolateContours() 
    iteratively.
    
    The interpolated contours and the contours used to interpolate them (the
    over-sampled original contours) are then transformed from the Fixed to 
    Moving image domain.
    
    Then the intersecting points between the lines that link each point in the
    first over-sampled contour and the interpolated contour, and each point in
    the interpolated contour and the second over-sampled contour, with the 
    imaging planes is found. These points form contours that lie in the Moving
    image planes.
    
    
    Inputs:
        FixedDicomDir     - Directory containing the DICOMs for the Fixed 
                            image domain.
        
        MovingDicomDir    - Directory containing the DICOMs for the Moving 
                            image domain.
                             
        FixedRoiFpath     - Full path to the ROI Collection RT-STRUCT file.
                             
        dP                - A float denoting the desired inter-node spacing to 
                            use when super-sampling the contours used for 
                            interpolation.
                         
        InterpolateAllPts - A boolean that determines whether or not 
                            interpolation will be performed only between a pair
                            of contour points for which the corresponding node
                            in the longer of the original contours was an
                            original node.  See Note 4 in the function 
                            InterpolateBetweenContours.
                            
        LogToConsole      - Log some results to the console.
                            
        PlotResults       - Boolean that denotes whether or not results will be
                            plotted.
        
        ExportResults     - Boolean that denotes whether or not plotted results
                            are to be exported.
        
    Returns:
        PointData         - Dictionary containing the data on all points in all
                            contours for the Fixed image and with added entries 
                            for interpolated points.
        
        FixContourData    - Dictionary containing the data on all points 
                            arranged slices, including the interpolated points.
        
        InterpData        - A dictionary containing interpolated data.
    
    """
    
    # Import packages and functions:
    import time
    from GetImageAttributes import GetImageAttributes
    from GetInputPoints import GetInputPoints
    from PCStoICS import PCStoICS
    
    # Force ExportResults to False if PlotResults is False:
    if not PlotResults:
        ExportResults = False
    
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Chose which package to use to get the Image Plane Attributes:
    #package = 'sitk'
    package = 'pydicom'
    
    # Get the Image Attributes for the Fixed image:
    FixOrigin, FixDirs,\
    FixSpacings, FixDims = GetImageAttributes(DicomDir=FixedDicomDir, 
                                              Package=package)
    
    
    # Get contour points into necessary arrays:
    FixPtsPCS, FixPtsBySliceAndContourPCS,\
    FixPtsICS, FixPtsBySliceAndContourICS,\
    LUT, PointData, FixContourData = GetInputPoints(FixedDicomDir=FixedDicomDir, 
                                                    FixedRoiFpath=FixedRoiFpath)
    
    
    # Get the indices of the slices that contain contours in the Fixed image:
    FixSlicesWithContours = GetIndsOfSlicesWithContours(ContourData=FixContourData,
                                                        ContourType=1)
    
    print('Contours exist on slices', FixSlicesWithContours)
    
    # Initialise a dictionary to store the interpolation results:
    InterpData = {'InterpSliceInd'       : [],
                  'BoundingSliceInds'    : [],
                  'FixOSIsOrigNode1'     : [],
                  'FixOSIsOrigNode2'     : [],
                  'FixOSContour1Inds'    : [],
                  'FixOSContour2Inds'    : [],
                  'FixInterpContourInds' : [],
                  'FixOSContour1Pts'     : [],
                  'FixOSContour2Pts'     : [],
                  'FixInterpContourPts'  : []
                  }
    
    # Iterate through all FixSlicesWithContours and interpolate between each
    # two slices:
    for i in range(len(FixSlicesWithContours) - 1):
        InterpSliceInd = FixSlicesWithContours[i] + 0.5
        
        if LogToConsole:
            print(f'Attempting to interpolate at slice {InterpSliceInd}...')
    
    
        #InterpContour = InterpolateContours(FixContourData, PointData, InterpSliceInd, dP)
        
        BoundingSliceInds,\
        OSContour1,\
        OSIsOrigNode1,\
        OSContour2,\
        OSIsOrigNode2,\
        InterpContourPCS = InterpolateContours(FixContourData, PointData, 
                                               InterpSliceInd, dP,
                                               InterpolateAllPts)
    
        # Continue only if None was not returned for InterpContourPCS:
        if InterpContourPCS:
            #print(f'len(InterpContourPCS) = {len(InterpContourPCS)}')
    
            # The maximum contour number in the original contour data:
            LastCntNo = PointData['InContourNo'][-1]
    
            # Add interpolated contour to PointData:
            for i in range(len(InterpContourPCS)):
                PointData['PointNo'].append(PointData['PointNo'][-1] + 1)
                PointData['InSliceNo'].append(InterpSliceInd)
                PointData['InContourNo'].append(LastCntNo + 1)
                PointData['InContourType'].append(2)
                PointData['InPointIndex'].append([])
                PointData['InPointPCS'].append(InterpContourPCS[i])
                PointData['InPointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS[i], 
                                                        Origin=FixOrigin, 
                                                        Directions=FixDirs, 
                                                        Spacings=FixSpacings))
    
            # Add interpolated contour to FixContourData:
            FixContourData['SliceNo'].append(InterpSliceInd)
            FixContourData['ContourType'].append(2)
            #FixContourData['PointPCS'].append(InterpContourPCS)
            #FixContourData['PointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS, 
            #                                           Origin=FixOrigin, 
            #                                           Directions=FixDirs, 
            #                                           Spacings=FixSpacings))
            """
            The list of contour points needs to be in a list structure to be 
            consistent with the other contour data in the dictionary.  
            Should verify this (July 23).
            """
            FixContourData['PointPCS'].append([InterpContourPCS])
            InterpContourICS = PCStoICS(Pts_PCS=InterpContourPCS, 
                                        Origin=FixOrigin, 
                                        Directions=FixDirs, 
                                        Spacings=FixSpacings)
            FixContourData['PointICS'].append([InterpContourICS])
        
            
            # Store the interpolation output in InterpData:
            InterpData['InterpSliceInd'].append(InterpSliceInd)
            InterpData['BoundingSliceInds'].append(BoundingSliceInds)
            InterpData['FixOSIsOrigNode1'].append(OSIsOrigNode1)
            InterpData['FixOSIsOrigNode2'].append(OSIsOrigNode2)
            InterpData['FixOSContour1Pts'].append(OSContour1)
            InterpData['FixOSContour2Pts'].append(OSContour2)
            InterpData['FixInterpContourPts'].append(InterpContourPCS)
                         
   
                
     
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'Took {Dtime} s to interpolate {len(FixSlicesWithContours)}',
              'contours.\n\n')
        
    
    # Transform the interpolated contours and over-sampled contours used to
    # interpolate:
    if LogToConsole:
        print('Registering images...')
 
    InterpData, FixIm, MovIm, RegIm,\
    ElastixImFilt = TransformInterpolatedContours(InterpData, FixedDicomDir, 
                                                  MovingDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'Took {Dtime} s to register the Moving image stack to the',
              'Fixed image domain and transform the contours.\n\n')
    
    # Generate contours in the Moving image planes:
    UseInterp = True
    #UseInterp = False
    
    OnLineSegment = True
    
    MovContourData = GetIntersectingContoursInMovingPlanes(InterpData, 
                                                           MovingDicomDir,
                                                           UseInterp,
                                                           OnLineSegment)
    
    # Get the indices of the slices that contain contours in the Moving image:
    MovSlicesWithContours = GetIndsOfSlicesWithContours(ContourData=MovContourData,
                                                        ContourType=3)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    if LogToConsole:
        print(f'Took {Dtime} s to copy {len(FixSlicesWithContours)} contours',
              f'in the Fixed image domain to {len(MovSlicesWithContours)}',
              'contours in the Moving image domain.')
    
    
    """ Plot various results: """
    if PlotResults:
        import RegUtilityFuncs as ruf
        import time
        
        #""" Get additional data for plots: """
        #import copy
        
        # Get the bounding slice indices for the slice at InterpSliceInd:
        #BoundingSliceInds = GetBoundingSliceInds(FixContourData, PointData, InterpSliceInd)

        #Contour1 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[0]][0])
        #Contour2 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[1]][0])

        #if InterpSliceInd.is_integer():
        #    # The actual contour that exists at slice InterpSliceInd:
        #    ActualContour = copy.deepcopy(FixContourData['PointPCS'][InterpSliceInd][0])


        # Combine data into lists to plot:
        #Contours = [Contour1, Contour2, InterpContourPCS]
        #Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
        #          f'Contour at slice {BoundingSliceInds[1]}', 
        #          f'Interpolated contour at {InterpSliceInd}']
        #Colours = ['b', 'g', 'r', 'y']
        ##Shifts = [0, 2, 4, 6] # to help visualise
        #Shifts = [0, 0, 0, 0] # no shift
        #if InterpSliceInd.is_integer():
        #    Contours.append(ActualContour)
        #    Labels.append(f'Actual contour at {InterpSliceInd}')

        # Plot the interpolated and original contours:
        PlotInterpContours2D(InterpData=InterpData, 
                             FixedOrMoving='Fixed', 
                             dP=dP, 
                             AnnotatePtNums=False, 
                             SubPlots=True, 
                             ExportPlot=ExportResults)
                             
        # Plot the transformed interpolated and original contours:
        PlotInterpContours2D(InterpData=InterpData,
                             FixedOrMoving='Moving', 
                             dP=dP, 
                             AnnotatePtNums=False, 
                             SubPlots=True, 
                             ExportPlot=ExportResults)
                             
        # Plot the intersecting points of the transformed interpolated
        # and original contours with the image planes in the Moving
        # domain:
        PlotIntersectingPts2D(MovContourData=MovContourData, 
                              dP=dP, 
                              AnnotatePtNums=False, 
                              SubPlots=True, 
                              ExportPlot=ExportResults)
                
                
                
        """ Plot DICOM images containing contours overlaid: """        
        
        # Plot all slices containing contours in either image domains.
        # Concatenate the lists:
        FixOrMovSlicesWithContours = FixSlicesWithContours + MovSlicesWithContours
        
        # Reduce to set of unique indices:
        FixOrMovSlicesWithContours = list(set(FixOrMovSlicesWithContours))
        
        # Choose the perspective:
        Perspective = 'axial'
        #Perspective = 'sagittal'
        #Perspective = 'coronal'
        
        # Choose whether to plot contour points as lines or dots:
        ContoursAs = 'lines'
        #ContoursAs = 'dots'

        # Create the filename for the exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                   + f'_Result_of_registration_and_interpolation.png'
        
        ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, 
                                                                         mov_im=MovIm, 
                                                                         reg_im=RegIm,
                                                                         fix_pts=FixContourData['PointICS'], 
                                                                         mov_pts=MovContourData['PointICS'],
                                                                         export_fig=ExportResults,
                                                                         export_fname=FigFname,
                                                                         plot_slices=FixOrMovSlicesWithContours,
                                                                         perspective=Perspective,
                                                                         contours_as=ContoursAs,
                                                                         LogToConsole=False)
    
    
    return PointData, FixContourData, InterpData, MovContourData
    
    



    

"""
******************************************************************************
Supplementary functions : 
******************************************************************************
"""



def ReduceNodesOfContour(Contour, IndsToKeep):
    """
    Reduce the number of nodes in a list of contour points.
    
    Inputs:
        Contour        - A list of super-sampled contour points.
        
        IndsToKeep     - A list of booleans whose values denote whether the   
                         corresponding contour point in Contour is to remain
                         (True) or not (False).
        
    Returns:
        ReducedContour - A reduced list of contour points.
        
        
    Note:
        This function is similar to ReduceNodesOfContours() but it operates on
        a single set of contour points.
    """
    
    ReducedContour = []
    
    for i in range(len(Contour)):
        if IndsToKeep[i]:
            ReducedContour.append(Contour[i])
            
    return ReducedContour



def GetSegmentLengths(Contour):
    """
    Generate a list of segment lengths for each segment in a contour.
    
    Inputs:
        Contour   - A list of contour points.
        
    Returns:
        SegmentLs - A list of segment lengths.
    
    """
    
    # Initialise the list of segment lengths:
    SegmentLs = []
    
    # Loop through each point:
    for i in range(len(Contour) - 1):
        # The segment length between every two points:
        SegmentL = (\
                   (Contour[i+1][0] - Contour[i][0])**2 \
                   + (Contour[i+1][1] - Contour[i][1])**2 \
                   + (Contour[i+1][2] - Contour[i][2])**2 \
                   )**(1/2)
              
        SegmentLs.append(SegmentL)
        
    return SegmentLs





def GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo):
    ContourTypes = ContourData['ContourType']
    
    inds = []
    
    for i in range(len(ContourTypes)):
        if ContourTypes[i] == ContourTypeNo:
            inds.append(i)
            
    return inds



def GetPointsInImagePlanes(ContourData, ContourTypeNo, DicomDir):
    """
    Create a closed "contour" of the 4 (+ 1 to close the contour) points that
    define the extent of the imaging planes.
    """
    
    from GetImageAttributes import GetImageAttributes
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
        
    # Get the slice numbers of the imaging planes that contain contours:
    SliceNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    # Create "contours" that mark out corners that define each plane in 
    # SliceNums:
    Planes = []
    
    for i in range(len(SliceNums)):
        s = SliceNums[i]
        
        # Let Point0 be equivalent to the origin moved along s*Spacings[2] in 
        # the z-direction (Dirs[6:8]):
        Point0 = [Origin[0] + s*Spacings[2]*Dirs[6], 
                  Origin[1] + s*Spacings[2]*Dirs[7], 
                  Origin[2] + s*Spacings[2]*Dirs[8]
                  ]
        
        # Let Point1 be equivalent to Point0 moved along Dims[0]*Spacings[0] in
        # the x-direction (Dirs[0:2]):
        Point1 = [Point0[0] + Dims[0]*Spacings[0]*Dirs[0],
                  Point0[1] + Dims[0]*Spacings[0]*Dirs[1],
                  Point0[2] + Dims[0]*Spacings[0]*Dirs[2]
                  ]
        
        # Let Point2 be equivalent to Point1 moved along Dims[1]*Spacings[1] in 
        # the y-direction (Dirs[3:5]):
        Point2 = [Point1[0] + Dims[1]*Spacings[1]*Dirs[3],
                  Point1[1] + Dims[1]*Spacings[1]*Dirs[4],
                  Point1[2] + Dims[1]*Spacings[1]*Dirs[5]
                  ]
        
        # Let Point3 be equivalent to Point0 moved along Dims[1]*Spacings[1] in
        # the y-direction (Dirs[3:5]):
        Point3 = [Point0[0] + Dims[1]*Spacings[1]*Dirs[3],
                  Point0[1] + Dims[1]*Spacings[1]*Dirs[4],
                  Point0[2] + Dims[1]*Spacings[1]*Dirs[5]
                  ]
        
        
        Planes.append([Point0, Point1, Point2, Point3, Point0])
        
    
    return Planes
            


def GetGridOfPointsInImagePlanes(ContourData, ContourTypeNo, DicomDir):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes.
    """
    
    from GetImageAttributes import GetImageAttributes
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
        
    # Get the slice numbers of the imaging planes that contain contours:
    SliceNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    # Create "contours" that mark out corners that define each plane in 
    # SliceNums:
    Planes = []
    
    for s in SliceNums:
        
        Points = []
        
        for i in range(Dims[0]):
            for j in range(Dims[1]):
                x = Origin[0] + i*Spacings[0]*Dirs[0] + j*Spacings[1]*Dirs[3] + s*Spacings[2]*Dirs[6]
                
                y = Origin[1] + i*Spacings[0]*Dirs[1] + j*Spacings[1]*Dirs[4] + s*Spacings[2]*Dirs[7]
                
                z = Origin[2] + i*Spacings[0]*Dirs[2] + j*Spacings[1]*Dirs[5] + s*Spacings[2]*Dirs[8]
                
                Points.append([x, y, z])
          
        Planes.append(Points)
        

    return Planes
 

def ExportContourPtsToCsv(ContourData, Label):
    """
    ContourData - e.g. FixContourData, MovContourData
    
    Label       - e.g. 'Fix', 'Mov'
    """
    
    import csv

    for i in range(len(ContourData['PointPCS'])):
        Points = ContourData['PointPCS'][i]
        
        if Points:
            Points = Points[0]
            
            SliceNo = ContourData['SliceNo'][i]
            
            with open(Label + f"ContourPts_slice{SliceNo}.csv", "w") as output:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerows(Points)
                
                #for r in range(len(Points)):
                #    writer.writerow(Points[r])
                
    return
    

    
    

def PlotInterpContours2D_OLD(Contours, Labels, Colours, Shifts, 
                             InterpSliceInd, BoundingSliceInds, dP, 
                             AnnotatePtNums, ExportPlot):
    """
    Plot interpolation results.
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    MarkerSize = 5
    AnnotationPeriod = 5 # annotate every AnnotationPeriod^th point
    #AnnotationPeriod = 10 # annotate every AnnotationPeriod^th point
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14), dpi=300)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    
    for i in range(len(Contours)):
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
    
        for x, y, z in Contours[i]:
            #X.append(x)
            X.append(x + Shifts[i])
            #Y.append(y)
            Y.append(y + Shifts[i])
    
        # Define linestyle:
        if 'nterp' in Labels[i]: # i.e. interpolated contour
            if 'uper' in Labels[i]: # super-sampled
                LineStyle='dashed'
            else: # reduced nodes
                LineStyle=(0, (5, 10)) # loosely dashed   
        else:
            LineStyle='solid'
            
        # Plot line:
        ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[i], label=Labels[i]);    
        ax.legend(loc='upper left', fontsize='large')
        # Plot dots:
        if True:
            plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[i]);
        
        # Plot the first and last points with different markers to help
        # identify them:
        #plt.plot(X[0], Y[0], '>', markersize=MarkerSize + 5, c=Colours[i]);
        #plt.plot(X[-1], Y[-1], 's', markersize=MarkerSize + 5, c=Colours[i]);
        # Annotate with point numbers:
        if AnnotatePtNums:
            P = list(range(len(X))) # list of point numbers
            # Annotate every point:
            #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
            #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
            if False:
                plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5, c=Colours[i]);
            # Annotate every AnnotationPeriod^th point with the point number:
            for p in range(0, len(X), AnnotationPeriod):
                plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[i]);
        # Also annotate the last point:
        #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
        
    plt.axis('equal')
    
    plt.title(f'Interpolation of contour at slice {InterpSliceInd} between ' \
              + f'slices {BoundingSliceInds[0]} and {BoundingSliceInds[1]}')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_at_slice_{InterpSliceInd}_of_contours_at_' \
                + f'{BoundingSliceInds[0]}_and_{BoundingSliceInds[1]}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    

def PlotInterpContours2D(InterpData, FixedOrMoving, dP, AnnotatePtNums, 
                         SubPlots, ExportPlot):
    """
    Plot interpolation results.
    
    Inputs:
        InterpData     - Dictionary containing interpolated data
        
        FixedOrMoving  - String; Acceptable values are "Fixed" or "Moving"
        
        dP             - Float value; Minimum inter-node spacing used for 
                         interpolation
        
        AnnotatePtNums - Boolean value determines whether points are annotated with
                         point numbers
                     
        SubPlots       - Boolean value determines whether a single figure with 
                         sub-plots or individual plots are generated.
                    
    Returns:
        None; Plot exported if ExportPlot = True
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    MarkerSize = 5
    
    Colours = ['b', 'g', 'r']
    
    #AnnotationPeriod = 5 # annotate every AnnotationPeriod^th point
    AnnotationPeriod = 10 # annotate every AnnotationPeriod^th point
    
    # Get number of interp data sets in InterpData:
    N = len(InterpData['InterpSliceInd'])
    
    # Get the number of rows and columns required of the subplot:
    if N == 1:
        Nrows = 1
        Ncols = 1
    elif N > 1:
        Ncols = 2
        Nrows = - (- N//Ncols) # i.e. floor of N / Ncols
    else:
        return
        
    # Prepare the figure:
    if ExportPlot:
        DPI = 300
    else:
        DPI = 100
        
    if SubPlots:
        fig = plt.figure(figsize=(5*Ncols, 5*Nrows), dpi=DPI)

    
    for i in range(N):
        # Set up the axes for this sub-plot:
        if SubPlots:
            ax = fig.add_subplot(Nrows, Ncols, i+1)
        else:
            fig = plt.figure(figsize=(10, 10), dpi=DPI)
            ax = fig.add_subplot(111)
        
        if 'ix' in FixedOrMoving:
            Contour1 = InterpData['FixOSContour1Pts'][i]
            Contour2 = InterpData['FixOSContour2Pts'][i]
            ContourI = InterpData['FixInterpContourPts'][i]
        elif 'ov' in FixedOrMoving:
            Contour1 = InterpData['MovOSContour1Pts'][i]
            Contour2 = InterpData['MovOSContour2Pts'][i]
            ContourI = InterpData['MovInterpContourPts'][i]
        else:
            print('Input "FixedOrMoving" must be "Fixed" or "Moving".')
            return
        
        AllContours = [Contour1, Contour2, ContourI]
        
        SliceInd1 = InterpData['BoundingSliceInds'][i][0]
        SliceInd2 = InterpData['BoundingSliceInds'][i][1]
        SliceIndI = InterpData['InterpSliceInd'][i]
        
        AllLabels = [f'{SliceInd1}', 
                     f'{SliceInd2}', 
                     f'{SliceIndI}']
        
        for j in range(len(AllContours)):
            # Unpack tuple and store each x,y tuple in arrays X and Y:
            X = []
            Y = []
        
            for x, y, z in AllContours[j]:
                X.append(x)
                #X.append(x + Shifts[j])
                Y.append(y)
                #Y.append(y + Shifts[j])
        
            # Define linestyle:
            if 'nterp' in AllLabels[j]: # i.e. interpolated contour
                if 'uper' in AllLabels[j]: # super-sampled
                    LineStyle='dashed'
                else: # reduced nodes
                    LineStyle=(0, (5, 10)) # loosely dashed   
            else:
                LineStyle='solid'
                
            # Plot line:
            ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]);    
            ax.legend(loc='upper left', fontsize='small')
            # Plot dots:
            if True:
                plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[j]);
            
            # Plot the first and last points with different markers to help
            # identify them:
            #plt.plot(X[0], Y[0], '>', markersize=MarkerSize + 5, c=Colours[i]);
            #plt.plot(X[-1], Y[-1], 's', markersize=MarkerSize + 5, c=Colours[i]);
            # Annotate with point numbers:
            if AnnotatePtNums:
                P = list(range(len(X))) # list of point numbers
                # Annotate every point:
                #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
                #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
                if False:
                    plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5, c=Colours[j]);
                # Annotate every AnnotationPeriod^th point with the point number:
                for p in range(0, len(X), AnnotationPeriod):
                    plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[j]);
            # Also annotate the last point:
            #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
        
        plt.axis('equal')
    
        plt.title('Interpolation between ' + FixedOrMoving \
                  + f' slices {SliceInd1} and {SliceInd2}')
        
        if not SubPlots:
            # Create filename for exported figure:
            FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                       + '_Interpolation_between_' + FixedOrMoving \
                       + f'_slices_{SliceInd1}_and_{SliceInd2}_dP_{dP}.png'
            
            if ExportPlot:
                plt.savefig(FigFname, bbox_inches='tight')
    
    FirstSlice = InterpData['BoundingSliceInds'][0][0]
    LastSlice = InterpData['BoundingSliceInds'][-1][1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                   + '_Interpolation_between_' + FixedOrMoving \
                   + f'_slices_{FirstSlice}_and_{LastSlice}_dP_{dP}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return


    
    
def PlotInterpolatedContours3D_OLD(InterpData, dP, ExportPlot):
    """
    Plot interpolation results.
    
    """
    
    import matplotlib.pyplot as plt
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import time
    
    FirstSliceInd = InterpData['BoundingSliceInds'][0][0]
    LastSliceInd = InterpData['BoundingSliceInds'][-1][1]
    
    
    NumPairs = len(InterpData['InterpSliceInd'])
    
    # Group all original and interpolated contours:
    FixContours = []
    MovContours = []
    
    # Store the ave Z values for each contour:
    FixAveZs = []
    MovAveZs = []
    
    for i in range(NumPairs):
        # Append Contour1:
        FixContour = InterpData['FixOSContour1Pts'][i]
        FixContours.append(FixContour)
        
        MovContour = InterpData['MovOSContour1Pts'][i]
        MovContours.append(MovContour)
        
        F = len(FixContour)
        M = len(MovContour)
        
        FixAveZ = sum([FixContour[i][2] for i in range(F)]) / F
        MovAveZ = sum([MovContour[i][2] for i in range(M)]) / M
        
        FixAveZs.append( [FixAveZ] * F )
        MovAveZs.append( [MovAveZ] * M )
        
        # Append the interpolated contour:
        FixContour = InterpData['FixInterpContourPts'][i]
        FixContours.append(FixContour)
        
        MovContour = InterpData['MovInterpContourPts'][i]
        MovContours.append(MovContour)
        
        F = len(FixContour)
        M = len(MovContour)
        
        FixAveZ = sum([FixContour[i][2] for i in range(F)]) / F
        MovAveZ = sum([MovContour[i][2] for i in range(M)]) / M
        
        FixAveZs.append( [FixAveZ] * F )
        MovAveZs.append( [MovAveZ] * M )
     
    # Append the last of Contour2:
    FixContour = InterpData['FixOSContour2Pts'][-1]
    FixContours.append(FixContour)
    
    MovContour = InterpData['MovOSContour2Pts'][-1]
    MovContours.append(MovContour)
    
    F = len(FixContour)
    M = len(MovContour)
    
    FixAveZ = sum([FixContour[i][2] for i in range(F)]) / F
    MovAveZ = sum([MovContour[i][2] for i in range(M)]) / M
    
    FixAveZs.append( [FixAveZ] * F )
    MovAveZs.append( [MovAveZ] * M )
    
    # Group all contours:
    AllContours = [FixContours, MovContours]
    
    AllAveZs = [FixAveZs, MovAveZs]
    
    #Colours = ['b', 'g', 'r', 'm', 'y']
    #Shifts = [0, 2, 4, 6, 8] # to help visualise
    #Shifts = [0, 0, 0, 0, 0] # no shift
    
    #MarkerSize = 5
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        #fig = plt.subplots(1, 2, figsize=(14, 14), dpi=300)
        fig = plt.figure(figsize=(10, 10), dpi=300)
    else:
        #fig = plt.subplots(1, 2, figsize=(14, 14))
        fig = plt.figure(figsize=(10, 10))
    #ax = plt.axes(projection="3d")


    
    # Loop through each list of contours for each sub-plot:
    for i in range(len(AllContours)):
        Contours = AllContours[i]
        
        AveZs = AllAveZs[i]
        
        # Set up the axes for this sub-plot:
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        
        # Loop through each contour:
        for j in range(len(Contours)):
            Contour = Contours[j]
            
            AveZ = AveZs[j]
            
            #print(len(Contour))
                
            # Store each x,y,z coordinate in arrays X, Y and Z:
            X = []
            Y = []
            Z = []
            
            # Unpack tuple and append to arrays X, Y and Z:
            for x, y, z in Contour:
                X.append(x)
                Y.append(y)
                Z.append(z)
        
            # Define linestyle:
            if (i % 2): # i.e. if i is even, it's an interpolated contour
                #LineStyle='dashed'
                LineStyle=(0, (5, 10)) # loosely dashed   
            else: # if i is odd, it's not an interpolated contour
                LineStyle='solid'
                
            # Plot line:
            ax.plot3D(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray'); 
            #Axes3D.plot(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray'); 
            
            #ax.legend(loc='upper left', fontsize='large')
            
            #print(len(X), len(Y), len(Z), len(AveZ))
            
            # Plot dots:
            #ax.scatter3D(X, Y, Z, markersize=MarkerSize, c=Z, cmap='hsv');
            #ax.scatter3D(X, Y, Z, c=Z, cmap='hsv');
            ax.scatter3D(X, Y, Z, c=AveZ, cmap='hsv');
            
        
        plt.title(f'Fixed image contours between slice {FirstSliceInd} and ' \
                  + f'{LastSliceInd} \nand interpolation inbetween slices')
        
        # Increment the sub-plot number:
        #i += 1

        
    
    
    #plt.title(f'Interpolation of contours between slice {FirstSliceInd} and ' \
    #          + f'slice {LastSliceInd}')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_bt_slice_{FirstSliceInd}_and_' \
                + f'{LastSliceInd}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    
    
    
    
def axisEqual3D(ax):
    """
    Source:
    https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio    
    """
    
    import numpy as np
    
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
    
            
def PlotInterpolatedContours3D(InterpData, FixContourData, MovContourData, dP,
                               PlotImagingPlanes, FixDicomDir, MovDicomDir, 
                               CombinePlots, ExportPlot):
    """
    Plot interpolation results.
    
    """
    
    import matplotlib.pyplot as plt
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    
    FirstSliceInd = InterpData['BoundingSliceInds'][0][0]
    LastSliceInd = InterpData['BoundingSliceInds'][-1][1]
    
    
    NumPairs = len(InterpData['InterpSliceInd'])
        
    
    # Group all original and interpolated contours:
    FixContours = []
    MovContours = []
    
    
    for i in range(NumPairs):
        # Append Contour1:
        FixContour = InterpData['FixOSContour1Pts'][i]
        FixContours.append(FixContour)
        
        MovContour = InterpData['MovOSContour1Pts'][i]
        MovContours.append(MovContour)
        
        # Append the interpolated contour:
        FixContour = InterpData['FixInterpContourPts'][i]
        FixContours.append(FixContour)
        
        MovContour = InterpData['MovInterpContourPts'][i]
        MovContours.append(MovContour)
        
     
    # Append the last of Contour2:
    FixContour = InterpData['FixOSContour2Pts'][-1]
    FixContours.append(FixContour)
    
    MovContour = InterpData['MovOSContour2Pts'][-1]
    MovContours.append(MovContour)
    
    # Group all contours:
    AllContours = [FixContours, MovContours]
    
    
    if PlotImagingPlanes:
        # Get points that define each imaging plane that contains contours:
        if True:
            FixPlanes = GetPointsInImagePlanes(ContourData=FixContourData, 
                                               ContourTypeNo=1, 
                                               DicomDir=FixDicomDir)
            
            MovPlanes = GetPointsInImagePlanes(ContourData=MovContourData, 
                                               ContourTypeNo=3, 
                                               DicomDir=MovDicomDir)
        if False:
            FixPlanes = GetGridOfPointsInImagePlanes(ContourData=FixContourData, 
                                                     ContourTypeNo=1, 
                                                     DicomDir=FixDicomDir)
            
            MovPlanes = GetGridOfPointsInImagePlanes(ContourData=MovContourData, 
                                                     ContourTypeNo=3, 
                                                     DicomDir=MovDicomDir)
        
        # Group all points belonging to the planes:
        AllPlanes = [FixPlanes, MovPlanes]                    
    
    
    Colours = ['b', 'r', 'g', 'm', 'y']
    MarkerSize = 3
    
    # Create a figure with two subplots and the specified size:
    #if ExportPlot:
    #    #fig = plt.subplots(1, 2, figsize=(14, 14), dpi=300)
    #    #fig = plt.figure(figsize=(10, 10), dpi=300)
    #    fig = plt.figure(figsize=plt.figaspect(0.5)*1.5) # Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
    #else:
    #    #fig = plt.subplots(1, 2, figsize=(14, 14))
    #    #fig = plt.figure(figsize=(10, 10))
    #    fig = plt.figure(figsize=plt.figaspect(0.5)*1.5) # Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
    #ax = plt.axes(projection="3d")

    #if CombinePlots:
    #    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    # Loop through each list of contours for each sub-plot:
    for i in range(len(AllContours)):
        Contours = AllContours[i]
        
        # Set up the axes for this sub-plot:
        #if not CombinePlots:
        #    ax = fig.add_subplot(2, 1, i+1, projection='3d')
        
        # Store all x,y,z coordinates in arrays X, Y and Z:
        X = []
        Y = []
        Z = []
            
        # Loop through each contour:
        for j in range(len(Contours)):
            Contour = Contours[j]
            
            #print(len(Contour))
                
            # Unpack tuple and append to arrays X, Y and Z:
            for x, y, z in Contour:
                X.append(x)
                Y.append(y)
                Z.append(z)
        
        # Define linestyle:
        if (i % 2): # i.e. if i is even, it's an interpolated contour
            #LineStyle='dashed'
            LineStyle=(0, (5, 10)) # loosely dashed   
        else: # if i is odd, it's not an interpolated contour
            LineStyle='solid'
            
            
        #fig = plt.figure(figsize=plt.figaspect(0.5)*1.5, projection='3d') # Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
        #fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
        fig = plt.figure(figsize=(14, 14))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot line:
        ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[i]); 
        #Axes3D.plot(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray'); 
        
        #ax.legend(loc='upper left', fontsize='large')
        
        #print(len(X), len(Y), len(Z), len(AveZ))
        
        # Plot dots:
        ##ax.scatter3D(X, Y, Z, markersize=MarkerSize, c=Z, cmap='hsv');
        #ax.scatter3D(X, Y, Z, c=Z, cmap='hsv');
        ax.scatter3D(X, Y, Z, s=MarkerSize, c=Colours[i], cmap='hsv');
        
        
        
        if PlotImagingPlanes:  
            #X = PlanesX[i]
            #Y = PlanesY[i]
            #Z = PlanesZ[i]
            #
            ##ax.plot_surface(X, Y, Z)
            #
            #XX, YY = np.meshgrid(X, Y)
            #   
            #ax.plot_surface(XX, YY, Z)
            
            Planes = AllPlanes[i]
            
            # Store all x,y,z coordinates in arrays X, Y and Z:
            X = []
            Y = []
            Z = []
                
            # Loop through each contour:
            for j in range(len(Planes)):
                Plane = Planes[j]
                
                #print(len(Contour))
                    
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in Plane:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
            
            # Plot line:
            ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c='gray'); 
            
        
            #plt.axis('equal') # doesn't work
            #ax.set_aspect('equal')
            #plt.axis('square')
            
            axisEqual3D(ax)
            
            if i == 0:
                txt = 'Fixed'
            else:
                txt = 'Moving'
                
            plt.title(txt + f' image contours between slice {FirstSliceInd}' \
                      + f' and {LastSliceInd} and interpolation inbetween' \
                      + ' slices')    
        
        if CombinePlots:
            plt.title('Fixed and Moving image contours between slice ' \
                      + f'{FirstSliceInd} and {LastSliceInd} \nand ' \
                      + 'interpolation inbetween slices')
        #else:
        #    plt.title(txt + f' image contours between slice {FirstSliceInd}' \
        #              + f' and {LastSliceInd} \nand interpolation inbetween' \
        #              + ' slices')
        
        # Increment the sub-plot number:
        i += 1

    
    #plt.title(f'Interpolation of contours between slice {FirstSliceInd} and ' \
    #          + f'slice {LastSliceInd}')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_bt_slice_{FirstSliceInd}_and_' \
                + f'{LastSliceInd}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    
    

def PlotIntersectingPoints2D_OLD(MovContourData, SliceNum, dP, AnnotatePtNums,
                                 ExportPlot):
    """
    Plot intersecting points.
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    MarkerSize = 5
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14), dpi=300)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    
    
    Points = MovContourData['PointPCS'][SliceNum]
    
    for i in range(len(Points)):
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
    
        for x, y, z in Points[i]:
            X.append(x)
            Y.append(y)

                
    # Plot line:
    ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c='b');    
    #ax.legend(loc='upper left', fontsize='large')
    # Plot dots:
    if True:
        plt.plot(X, Y, '.', markersize=MarkerSize, c='k');
    # Annotate with point numbers:
    if AnnotatePtNums:
        P = list(range(len(X))) # list of point numbers
        # Annotate every point:
        #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
        #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
        if False:
            plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5);
        # Annotate every 5th point with the point number:
        for p in range(0, len(X), 5):
            plt.text(X[p], Y[p], p, fontsize=MarkerSize+5);
        
    # Also annotate the last point:
    #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5);
        
    plt.axis('equal')
    
    plt.title(f'Intersecting points on Moving slice {SliceNum}')

    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Intersecting_points_on_slice_{SliceNum}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    return



def PlotIntersectingPts2D(MovContourData, dP, AnnotatePtNums, SubPlots, 
                          ExportPlot):
    """
    Plot intersecting points.
    
    AnnotatePtNums - Boolean value determines whether points are annotated with
                     point numbers
                     
    SubPlots       - Boolean value determines whether a single figure with 
                     sub-plots or individual plots are generated.
                     
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    #AnnotationPeriod = 5 # annotate every AnnotationPeriod^th point
    AnnotationPeriod = 10 # annotate every AnnotationPeriod^th point
    
    # Get the indices of the slices in ContourData that have the desired 
    # ContourType value:
    Inds = GetIndsOfSliceNumsOfContourType(ContourData=MovContourData, 
                                           ContourTypeNo=3)
    
    N = len(Inds)
    
    # Get the number of rows and columns required of the subplot:
    if N == 1:
        Nrows = 1
        Ncols = 1
    elif N > 1:
        Ncols = 2
        Nrows = - (- N//Ncols) # i.e. floor of N / Ncols
    else:
        return
        
    # Prepare the figure:
    if ExportPlot:
        DPI = 300
    else:
        DPI = 100
        
    if SubPlots:
        fig = plt.figure(figsize=(5*Ncols, 5*Nrows), dpi=DPI)
        MarkerSize = 2
    else:
        MarkerSize = 5

    
    for i in range(N):
        # Set up the axes for this sub-plot:
        if SubPlots:
            ax = fig.add_subplot(Nrows, Ncols, i+1)
        else:
            fig = plt.figure(figsize=(10, 10), dpi=DPI)
            ax = fig.add_subplot(111)
        
        SliceNum = Inds[i]
        
        Points = MovContourData['PointPCS'][SliceNum]
        
        for j in range(len(Points)):
            # Unpack tuple and store each x,y tuple in arrays X and Y:
            X = []
            Y = []
        
            for x, y, z in Points[j]:
                X.append(x)
                Y.append(y)
                
            # Plot line:
            ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c='b');    
            # Plot dots:
            if True:
                plt.plot(X, Y, '.', markersize=MarkerSize, c='k');
            
            # Plot the first and last points with different markers to help
            # identify them:
            #plt.plot(X[0], Y[0], '>', markersize=MarkerSize + 5, c=Colours[i]);
            #plt.plot(X[-1], Y[-1], 's', markersize=MarkerSize + 5, c=Colours[i]);
            # Annotate with point numbers:
            if AnnotatePtNums:
                P = list(range(len(X))) # list of point numbers
                # Annotate every point:
                #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
                #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
                if False:
                    plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5, c='k');
                # Annotate every AnnotationPeriod^th point with the point number:
                for p in range(0, len(X), AnnotationPeriod):
                    plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c='k');
            # Also annotate the last point:
            #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
        
        #plt.axis('tight') 
        plt.axis('equal')
        
    
        plt.title(f'Intersecting points on Moving slice {SliceNum}')
        
        if not SubPlots:
            # Create filename for exported figure:
            FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                       + f'_Intersecting_points_on_slice_{SliceNum}_dP_{dP}.png'
            
            if ExportPlot:
                plt.savefig(FigFname, bbox_inches='tight')
    
    FirstSlice = Inds[0]
    LastSlice = Inds[-1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_Intersecting' \
                   + f'_points_between_slices_{FirstSlice}_and_{LastSlice}_dP_{dP}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return
    

    
def PlotJoinedPoints2D(PointsJoined, ExportPlot, SliceNum):
    """
    Plot joined points.
    
    """
    
    import matplotlib.pyplot as plt
    import time
    import numpy as np
    
    MarkerSize = 5
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    
    # Create a figure with two subplots and the specified size:
    if ExportPlot:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14), dpi=300)
        #fig, ax = plt.subplots(1, 1, figsize=(14, 14), dpi=150)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    
    
    # Loop through each collection of points joined to every other point:
    for i in range(len(PointsJoined)):
        Points = PointsJoined[i]
        
        # Unpack tuple and store each x,y tuple in arrays X and Y:
        X = []
        Y = []
        
        for j in range(len(Points)):
            # Unpack tuple and store each x,y tuple in arrays X and Y:
            #X = []
            #Y = []
        
            for x, y, z in Points[j]:
                X.append(x)
                Y.append(y)
    
        
        # Generate a random colour for this contour:
        colour = [item/255 for item in list(np.random.choice(range(256), size=3))]
                            
                            
        # Plot line:
        ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=colour);    
        #ax.legend(loc='upper left', fontsize='large')
        # Plot dots:
        plt.plot(X, Y, '.', markersize=MarkerSize, c='k');
            
            
        if False:
            for i in range(len(Points)):
                # Plot line:
                ax.plot(Points[i][0], Points[i][1], linestyle=LineStyle, linewidth=1, c='b');    
                #ax.legend(loc='upper left', fontsize='large')
                # Plot dots:
                plt.plot(Points[i][0], Points[i][1], '.', markersize=MarkerSize, c='k');
                # Annotate with point numbers:
                if False:#AnnotatePtNums:
                    P = list(range(len(X))) # list of point numbers
                    # Annotate every point:
                    #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
                    #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
                    if False:
                        plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5);
                    # Annotate every 5th point with the point number:
                    for p in range(0, len(X), 5):
                        plt.text(X[p], Y[p], p, fontsize=MarkerSize+5);
                    
                # Also annotate the last point:
                #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5);
        
    plt.axis('equal')
    
    plt.title(f'Every point in intersecting points joined to every other' \
              + f'point for slice no {SliceNum}')

    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Intersecting_pts_for_slice_{SliceNum}_joined_to_every_' \
                + 'other_point.png'
    
    if ExportPlot:
        print('exporting plot...')
        plt.savefig(FigFname, bbox_inches='tight')
    