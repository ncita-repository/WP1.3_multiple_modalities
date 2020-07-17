# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:50:23 2020

@author: ctorti
"""



""" Contour interpolating functions """


def GetIndsOfSlicesWithContours(ContourData):
    # Initialise list of indices of the slices that contain contours:
    SlicesWithContours = []
    
    for i in range(len(ContourData['ContourType'])):
        if ContourData['ContourType'][i] == 1:
            SlicesWithContours.append(i)
            
    # Reduce to list of unique indices:
    SlicesWithContours = list(set(SlicesWithContours))
            
    return SlicesWithContours


def GetBoundingSliceInds(ContourData, PointData, InterpInd):
    """ 
    Note 1: 
        ContourData is a dictionary containing the key ['InPointPCS'], whose
        values are a list of contour points arranged along rows for 
        different slices and along columns for different contours (i.e. 
        ContourData is a list of lists). Slices without contours have empty 
        list ([]).
        
    Note 2:
        _getBoundingPair checks if a contour is an interpolated contour, and if 
        so, it is excluded from the list of possible contours to interpolate
        between (hence the index of an interpolated contour cannot be included 
        in ContourPair). This has not yet been implemented in 
        GetBoundingSliceInds.
        
    Note 3: 
        Like _getBoundingPair, GetBoundingSliceInds will exclude any slices 
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
    if InterpInd <= Extent[0] or InterpInd >= Extent[1]:
        print(f'Cannot interpolate slice number {InterpInd}')
        return []
    
    else:
        #print(f'Looking for the nearest lower slice to {InterpInd}...')
        
        # Iterate backwards from (InterpInd - 1) to SlicesWithContours[0] to determine 
        # the nearest lowest slice index that contains a single contour:
        if False:
            for i in range(InterpInd - 1, min(SlicesWithContours) - 1, -1):
                #print(f'\nContourData['InPointPCS'][{i}] =', ContourData['InPointPCS'][i])
                #print(f'\nlen(ContourData['InPointPCS'][{i}]) = {len(ContourData['InPointPCS'][i])}')
                
                if len(ContourData['InPointPCS'][i]) == 1:
                    #print(f'\nlen(ContourData['InPointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundLowerSlice = True
                    
                    break
            
        """ Allow for the possibility of InterpInd to be a fraction (not an
        integer) by taking floor of InterpInd using (InterpInd//1): """
        FloorInterpInd = int(InterpInd//1)
        
        for i in range(FloorInterpInd, min(SlicesWithContours) - 1, -1):
            if len(ContourData['InPointPCS'][i]) == 1:
                # Append i to BoundingSliceInds if i != InterpInd (which could
                # arise if InterpInd is an integer):
                if i != InterpInd:
                    BoundingSliceInds.append(i)
                    
                    FoundLowerSlice = True
                    
                    break
                
        
        
                
                
        #print(f'Looking for the nearest upper slice to {InterpInd}...')
        
        # Iterate forwards from (InterpInd + 1) to SlicesWithContours[-1] to determine 
        # the nearest upper slice index that contains a single contour:
        if False:
            for i in range(InterpInd + 1, max(SlicesWithContours) + 1):
                #print(f'\nContourData['InPointPCS'][{i}] =', ContourData['InPointPCS'][i])
                #print(f'\nlen(ContourData['InPointPCS'][{i}]) = {len(ContourData['InPointPCS'][i])}')
                #print(f'Appending {i} to BoundingSliceInds')
                
                if len(ContourData['InPointPCS'][i]) == 1:
                    #print(f'\nlen(ContourData['InPointPCS'][{i}]) == 1')
                    #print(f'Appending {i} to BoundingSliceInds')
                    
                    BoundingSliceInds.append(i)
                    
                    FoundUpperSlice = True
                    
                    break
                
        """ Allow for the possibility of InterpInd to be a fraction (not an
        integer) by taking ceiling of InterpInd using -(-InterpInd//1): """
        CeilInterpInd = int(-(-InterpInd//1))
        
        for i in range(CeilInterpInd, max(SlicesWithContours) + 1):
            if len(ContourData['InPointPCS'][i]) == 1:
                 # Append i to BoundingSliceInds if i != InterpInd (which could
                # arise if InterpInd is an integer):
                if i != InterpInd:
                    BoundingSliceInds.append(i)
                    
                    FoundUpperSlice = True
                    
                    break
        
        # Deal with case where no lower, or upper, or both upper nor lower 
        # slice index can be used for interpolation:
        if not (FoundLowerSlice and FoundUpperSlice):
            if not FoundLowerSlice and FoundUpperSlice:
                print('Did not find suitable lower slice to',
                      f'interpolate slice number {InterpInd}.')
                
            if not FoundUpperSlice and FoundLowerSlice:
                print('Did not find suitable upper slice to',
                      f'interpolate slice number {InterpInd}.')
                
            if not FoundLowerSlice and not FoundUpperSlice:
                print('Did not find suitable lower or upper slice',
                      f'to interpolate slice number {InterpInd}.')
            
            return []
        
        else:
            return BoundingSliceInds 
    




def GenerateClosedContour(Contour):
    # Import packages:
    #import copy
    
    #Contour = copy.deepcopy(Points)
    
    # Append the first point to create a closed contour:
    #Contour.append(Points[0])
    Contour.append(Contour[0])
    
    return Contour


def ReverseIfAntiClockwise(Contour):
    
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
    # List comprehension is compact but not reader-friendly outside Python:
    #return [CumPerimeters[i]/CumPerimeters[-1] for i in range(len(CumPerimeters))]
    
    # Initialise the list:
    CumPerimNorm = []
    
    for i in range(len(CumPerimeters)):
        CumPerimNorm.append(CumPerimeters[i]/CumPerimeters[-1])
        
    return CumPerimNorm



def GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm, NumNodesToAdd):
    """ Note:
        This function is called _getInterpolatedPerim in Matt Orton's code. 
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
    
    

def GetNodesToAddPerSegment(CumSuperSampledPerim, IsOrigNode):
    
    # Get the indices that sort the cumulative super-sampled perimeters:
    Inds = [i[0] for i in sorted(enumerate(CumSuperSampledPerim), key=lambda x:x[1])]
    
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
        

    # Order the cumulative super-sampled perimeters:
    """ This is not in Matt's code """
    #CumSSPerimSorted = [CumSuperSampledPerim[i] for i in Inds]

    # Get the list of cumulative super-sampled perimeters for the original 
    # nodes:
    """ This is not in Matt's code """
    #CumSSPerimSortedOrig = [CumSSPerimSorted[i] for i in OrigIndsSorted]

    # The segment lengths are:
    """ This is not in Matt's code """
    #SegmentL = [CumSSPerimSortedOrig[i+1] - CumSSPerimSortedOrig[i] for i in range(len(CumSSPerimSortedOrig) - 1)]
    
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



""" 
Estimate the area of the surface created by linking points from Contour1 to
points in Contour2, starting with the first point in Contour1 and a point in 
Contour2 given by the index StartInd2, and adding up all line lengths.

Note: There are N points in Contour1 and Contour2.
"""

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
    """ Note: Contour1 and Contour2 must have equal lengths. """
    
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



""" 
Using the function AddLineLengths(), loop through each N possible combinations 
of starting index StartInd2, and work out which starting index yields the 
minimum cumulative line length.
"""

def FindIndForMinCumLength(Contour1, Contour2):
    
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

    
    # Return the index that yields the minimum cumulative line length:
    return LineLengths.index(min(LineLengths))





""" 
Defining Contour1 and Contour2 as the two contours, shift the points in
Contour2 by StartInd (the amount that yields the minimum cumulative line length).
Then the points in the Contour1 and ROContour2, when joined point-by-point, 
will create the minimal surface area (i.e. catenoid) of the two contours:

https://en.wikipedia.org/wiki/Minimal_surface_of_revolution
"""


def ShiftContourPoints(Contour, Shift):
    
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
                               FracInterpSliceInd, Contour1HasMorePts):
    """ Note Contour1 and Contour2 must have equal lengths. """
    
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
        if True:
            if IsOrigNode[i]:
                InterpX = (1 - FracInterpSliceInd) * Contour1[i][0] \
                          + FracInterpSliceInd * Contour2[i][0]
                
                InterpY = (1 - FracInterpSliceInd) * Contour1[i][1] \
                          + FracInterpSliceInd * Contour2[i][1]
                          
                InterpZ = (1 - FracInterpSliceInd) * Contour1[i][2] \
                          + FracInterpSliceInd * Contour2[i][2]
                          
                InterpolatedContour.append([InterpX, InterpY, InterpZ])
            
        if False:
            """ Try something different from Matt's code: """
            InterpX = (1 - FracInterpSliceInd) * Contour1[i][0] \
                       + FracInterpSliceInd * Contour2[i][0]
                
            InterpY = (1 - FracInterpSliceInd) * Contour1[i][1] \
                      + FracInterpSliceInd * Contour2[i][1]
                      
            InterpZ = (1 - FracInterpSliceInd) * Contour1[i][2] \
                      + FracInterpSliceInd * Contour2[i][2]
        
            InterpolatedContour.append([InterpX, InterpY, InterpZ])
            
            
    return InterpolatedContour



""" Extra (supplementary) functions : """


""" New function, similar to but slight modification of ReduceNodesOfContours 
so that it operates on a single set of contour points.
"""
def ReduceNodesOfContour(Contour, IndsToKeep):
    
    ReducedContour = []
    
    for i in range(len(Contour)):
        if IndsToKeep[i]:
            ReducedContour.append(Contour[i])
            
    return ReducedContour



def GetSegmentLengths(Contour):
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



def InterpolateContours(ContourData, PointData, InterpSliceInd, dP):
    import copy
    
    # Get bounding slice indices:
    BoundingSliceInds = GetBoundingSliceInds(ContourData=ContourData,
                                             PointData=PointData,
                                             InterpInd=InterpSliceInd)
    
    
    # Escape if None was returned for BoundingSliceInds (i.e. if two slice 
    # indices were not found):
    if not BoundingSliceInds:
        return
    
    print(f'Interpolating between slice {BoundingSliceInds[0]} and',
          f'{BoundingSliceInds[1]} for slice at {InterpSliceInd}')
    
    # Define the contours:
    """
    Note: 
    
    Need to add [0] to the end of ContourData['InPointPCS'][BoundingSliceInds[0]] 
    since the contour points are a list within the list for each slice.
    """
    Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
    Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])
    
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
    Inds1, IsOrigNodeSorted1, IndsOrigNodes1, NodesToAddPerSegment1\
    = GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim1, 
                              IsOrigNode=IsOrigNode1)
    Inds2, IsOrigNodeSorted2, IndsOrigNodes2, NodesToAddPerSegment2\
    = GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim2, 
                              IsOrigNode=IsOrigNode2)
    
    # Super-sample the contours:
    SSContour1, SSIsOrigNode1 = SuperSampleContour(Contour=Contour1, 
                                                   NodesToAddPerSegment=NodesToAddPerSegment1)
    SSContour2, SSIsOrigNode2 = SuperSampleContour(Contour=Contour2, 
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
    OSContour1, \
    OSContour2, \
    OSIsOrigNode1, \
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
                                               Contour1HasMorePts=Contour1HasMorePts)
    
    """ 
    Note:
        
    It would be sensible to add the interpolated contour points to PointData
    and ContourData, but along wih the PCS points it would be wise to add the
    ICS points at the same time, but that requires the origin, cosine 
    directions and pixel spacings of the corresponding image...
    """
    
    
    
    return InterpContour


    
    """ Deviate from Matt's code - try interpolating using the super-sampled 
    contours (i.e. before reducing the num of nodes): """
    if False:
        InterpSSContour = InterpolateBetweenContours(Contour1=SSContour1,
                                                     Contour2=SSContour2,
                                                     IsOrigNode1=SSIsOrigNode1,
                                                     IsOrigNode2=SSIsOrigNode2,
                                                     FracInterpSliceInd=FracInterpInd,
                                                     Contour1HasMorePts=Contour1HasMorePts)
        
        # Define boolean list of indices to keep in the super-sampled interpolated
        # contour:
        IsOrigNode1or2 = IsOrigNode1 or IsOrigNode2
        
        # Reduce the number of nodes of InterpSSContour:
        ReducedInterpSSContour = ReduceNodesOfContour(Contour=InterpSSContour, IndsToKeep=IsOrigNode1or2)
        
        return InterpContour, InterpSSContour, ReducedInterpSSContour
    
    


    

def PlotInterpolationResults(Contours, Labels, Colours, Shifts, InterpSliceInd, 
                             BoundingSliceInds, dP, ExportPlot):
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
    