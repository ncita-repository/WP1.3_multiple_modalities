# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:50:23 2020

@author: ctorti
"""



""" 
Contour interpolating functions


The following dictionaries are used as inputs to some of the following 
functions and are appended to by others:
    
FixContourData - Dictionary of the contour points belonging to the Fixed 
(input) image, containing the keys:
    
              'SliceNo'     - List of integers (for original contour points) or
                              floats (for interpolated contours that correspond
                              to non-imaged planes, e.g. slice interpolated at
                              2.5 in between slices 2 and 3).
                              
              'ContourType' - List of integers for each DICOM slice; 
                              The integer value denotes whether there is or is
                              not contour data for any given slice, and if 
                              there is, whether it is an original contour,
                              an interpolated contour, or transformed contour.
                              0 -> No contour points 
                              1 -> Original contour points
                              2 -> Interpolated contour points
                              3 -> Transformed contour points
                              
              'PointPCS'    - List of input points in the Patient Coordinate 
                              System if contour data exists for any given slice 
                              (e.g. [[x0, y0, z0], [x1, y1, z1], ...], i.e. a 
                              list of lists of lists), or [] if not.
                            
              'PointICS'    - Same as 'PointPCS' but with points in the
                              Image Coordinate System.

              
                              
PointData   - Dictionary containing the keys:
    
              'PointNo'      - List of integers of all points in ContourData.
              
              'InSliceNo'    - List of integers (for original contour points) 
                               or floats (for interpolated contours that 
                               correspond to non-imaged planes, e.g. slice 
                               interpolated at 2.5 in between slices 2 and 3).
                              
              'InContourNo'  - List of integers of the contour number each 
                               point belongs to from the input image.
              
              'ContourType'  - List of integers for each DICOM slice; 
                               The integer value denotes whether there is or is
                               not contour data for any given slice, and if 
                               there is, whether it is an original contour,
                               an interpolated contour, or transformed contour.
                               0 -> No contour points 
                               1 -> Original contour points
                               2 -> Interpolated contour points
                               3 -> Transformed contour points
                              
              'InPointIndex' - List of indeces [x, y, z] for each point in the
                               input image.
                              
              'InPointPCS'   - List of points in the Patient Coordinate System 
                               if contour data exists for any given slice (e.g. 
                               [[x0, y0, z0], [x1, y1, z1], ...], i.e. a list  
                               of lists of lists), or [] if not.
                            
              'InPointICS'   - Same as 'InPointPCS' but with points in the
                               Image Coordinate System.
                               
              'OutPointIndex - Like 'InPointIndex' but for the output/
                               transformed image.
                               
              'OutPointPCS'  - Like 'InPointPCS' but for the output/transformed
                              image.
                            
              'OutPointICS' - Same as 'OutPointPCS' but with points in the 
                              Image Coordinate System.


Note:
    PointData is a flat dictionary, that contains contour data on a 
    point-by-point basis.
    
    InContourData is a multi-list dictionary, that contains contour data 
    arranged on a slice-by-slice and contour-by-contour basis for all slices
    within the input DICOM image (including slices that do not contain any
    contour data).
"""


def GetIndsOfSlicesWithContours(ContourData):
    """
    Get the indices of the rows (slices) within ContourData that contain a 
    single list of contour points (i.e. only one contour per slice).
    
    Inputs:
        ContourData         - Dictionary containing the key 'PointPCS', 
                              whose values are a list of contour points 
                              [x, y, z] arranged along rows for each DICOM 
                              slice, and along columns for different contours 
                              (i.e. a list of lists of lists). 
                              Slices without contours have empty list ([]).
                      
    Returns:
        SlicesWithContours - List of integers of the indices (row numbers) in
                             ContourData that contains a single list of contour
                             points.
    """
    
    # Initialise list of indices of the slices that contain contours:
    SlicesWithContours = []
    
    for i in range(len(ContourData['ContourType'])):
        if ContourData['ContourType'][i] == 1:
            SlicesWithContours.append(i)
            
    # Reduce to list of unique indices:
    SlicesWithContours = list(set(SlicesWithContours))
            
    return SlicesWithContours


def GetBoundingSliceInds(ContourData, PointData, InterpSliceInd):
    """ 
    Get the slice indices (row numbers) in ContourData of the nearest slices 
    that bound the slice to be interpolated.
    
    Inputs:
        ContourData         - Dictionary containing the key 'PointPCS', 
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
                             ContourData that contains a single list of contour
                             points.
    
        
    Notes:
        
        1. Matt Orton's function _getBoundingPair checks if a contour is an 
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
                #print(f'\nContourData['InPointPCS'][{i}] =', ContourData['InPointPCS'][i])
                #print(f'\nlen(ContourData['InPointPCS'][{i}]) = {len(ContourData['InPointPCS'][i])}')
                
                if len(ContourData['PointPCS'][i]) == 1:
                    #print(f'\nlen(ContourData['PointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundLowerSlice = True
                    
                    break
            
            
        """ Allow for the possibility of InterpInd to be a fraction (not an
        integer) by taking floor of InterpInd using (InterpInd//1): """
        FloorInterpInd = int(InterpSliceInd//1)
        
        for i in range(FloorInterpInd, min(SlicesWithContours) - 1, -1):
            if len(ContourData['PointPCS'][i]) == 1:
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
                #print(f'\nContourData['PointPCS'][{i}] =', ContourData['PointPCS'][i])
                #print(f'\nlen(ContourData['PointPCS'][{i}]) = {len(ContourData['PointPCS'][i])}')
                #print(f'Appending {i} to BoundingSliceInds')
                
                if len(ContourData['PointPCS'][i]) == 1:
                    #print(f'\nlen(ContourData['PointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundUpperSlice = True
                    
                    break
                
        """ Allow for the possibility of InterpSliceInd to be a fraction (not 
        an integer) by taking ceiling of InterpSliceInd using 
        # -(-InterpSliceInd//1): """
        CeilInterpInd = int(-(-InterpSliceInd//1))
        
        for i in range(CeilInterpInd, max(SlicesWithContours) + 1):
            if len(ContourData['PointPCS'][i]) == 1:
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
    #MeanX = sum([Contour[i][0] for i in range(P)])
    
    MeanX = 0
    
    for i in range(P):
        MeanX += Contour[i][0]
    
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
        
        2. This function is called _getInterpolatedPerim in Matt Orton's code. 
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
    accounts for one, and the -1 accounts for the other. Matt's code also has
    -2 so it's probably correct.
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
        
        This is called _getIndicatorArray in Matt Orton's code
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
    Seems this isn't needed (I've likely misunderstood Matt's code).
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
    """ This is not in Matt's code """
    #CumSSPerimNormSorted = [CumSSPerimNorm[i] for i in Inds]

    # Get the list of normailsed cumulative super-sampled perimeters for the
    # original nodes:
    """ This is not in Matt's code """
    #CumSSPerimNormSortedOrig = [CumSSPerimNormSorted[i] for i in OrigIndsSorted]

    # The segment lengths are:
    """ This is not in Matt's code """
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
        NodesToAddPerSegment[i] so deviating from Matt Orton's code:
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
        
        4. In Matt Orton's code the interpolation is performed between the pair
        of contour points only if the corresponding contour point in the longer
        list of contour points was an original node. This results in an 
        interpolated contour with the same number of points in it as the number
        of True values in the longer of IsOrigNode1 and IsOrigNode2, and hence,
        the length of interpolated contour points is shorter than the length of
        Contour1 or Contour2.
        
        In this function I've allowed for either using the truncation that Matt
        uses (InterpolateAllPts = False) or to interpolate all pairs of points, 
        so that the length of Contour1, Contour2 and InterpolatedContour are 
        equal (InterpolateAllPts = True).
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



""" Extra (supplementary) functions : """



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



def InterpolateContours(ContourData, PointData, InterpSliceInd, dP, 
                        InterpolateAllPts):
    """
    Main function for contour interpolation.
    
    Inputs:
        ContourData       - Dictionary containing the key 'PointPCS', whose
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
    BoundingSliceInds = GetBoundingSliceInds(ContourData=ContourData,
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
    
    Need to add [0] to the end of ContourData['PointPCS'][BoundingSliceInds[0]] 
    since the contour points are a list within the list for each slice.
    """
    Contour1 = copy.deepcopy(ContourData['PointPCS'][BoundingSliceInds[0]][0])
    Contour2 = copy.deepcopy(ContourData['PointPCS'][BoundingSliceInds[1]][0])
    
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


    
    """ Deviate from Matt's code - try interpolating using the super-sampled 
    contours (i.e. before reducing the num of nodes): """
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
    
    


    

def PlotInterpolationResults2D(Contours, Labels, Colours, Shifts, InterpSliceInd, 
                               BoundingSliceInds, dP, ExportPlot):
    """
    Plot interpolation results.
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    #Colours = ['b', 'g', 'r', 'm', 'y']
    #Shifts = [0, 2, 4, 6, 8] # to help visualise
    #Shifts = [0, 0, 0, 0, 0] # no shift
    
    MarkerSize = 5
    
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
        
        #P = list(range(len(X))) # list of point numbers
        #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
        #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
        #plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5, c=Colours[i]);
        
        # Annotate every 5th point with the point number:
        if False:
            for p in range(0, len(X), 5):
                plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[i]);
            
        # Also annotate the last point:
        #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
        
    
    plt.title(f'Interpolation of contour at slice {InterpSliceInd} between ' \
              + f'slices {BoundingSliceInds[0]} and {BoundingSliceInds[1]}')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_at_slice_{InterpSliceInd}_of_contours_at_' \
                + f'{BoundingSliceInds[0]}_and_{BoundingSliceInds[1]}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    
    
    
    
def PlotInterpolationResults3D_OLD(InterpData, dP, ExportPlot):
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
        i += 1

        
    
    
    #plt.title(f'Interpolation of contours between slice {FirstSliceInd} and ' \
    #          + f'slice {LastSliceInd}')
    
    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + f'_Interpolation_bt_slice_{FirstSliceInd}_and_' \
                + f'{LastSliceInd}_dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    
    
def PlotInterpolationResults3D(InterpData, dP, ExportPlot):
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
        
        # Set up the axes for this sub-plot:
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        
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
            
        # Plot line:
        ax.plot3D(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray'); 
        #Axes3D.plot(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray'); 
        
        #ax.legend(loc='upper left', fontsize='large')
        
        #print(len(X), len(Y), len(Z), len(AveZ))
        
        # Plot dots:
        #ax.scatter3D(X, Y, Z, markersize=MarkerSize, c=Z, cmap='hsv');
        ax.scatter3D(X, Y, Z, c=Z, cmap='hsv');
            
        
        if i == 0:
            txt = 'Fixed'
        else:
            txt = 'Moving'
            
        plt.title(txt + f' image contours between slice {FirstSliceInd} and' \
                  + f' {LastSliceInd} \nand interpolation inbetween slices')
        
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
    