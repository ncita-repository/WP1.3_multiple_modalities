# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:50:23 2020

@author: ctorti
"""



""" Contour interpolating functions """



#def GetBoundingSliceInds(InterpInd, Extent, ContourData):
def GetBoundingSliceInds(InterpInd, SliceInds, ContourData):
    """ 
    Note 1: 
        ContourData is a list of contour points 
        arranged along rows for different slices and along 
        columns for different contours (i.e. ContourData 
        is a list of lists). Slices without contours have
        empty list ([]).
        
    Note 2:
        _getBoundingPair checks if a contour is an
        interpolated contour, and if so, it is excluded
        from the list of possible contours to interpolate
        between (hence the index of an interpolated contour
        cannot be included in ContourPair). This has not yet
        been implemented in GetBoundingSliceInds.
        
    Note 3: 
        Like _getBoundingPair, GetBoundingSliceInds will 
        exclude any slices that have multiple contours.
        A more robust version will be required to interpolate
        contours on slices with multiple contours present.
    """
    
    # Get extent of region of contours by slice number 
    # (equivalent to function _getExtentOfRegion in 
    # generateInterpolationData.js):
    Extent = [SliceInds[0], SliceInds[-1]]
    
    print('Extent =', Extent)


    # Initialise BoundingSliceInds:
    BoundingSliceInds = []
    
    # Initialise booleans:
    FoundLowerSlice = False
    FoundUpperSlice = False
    
    # Cannot interpolate a slice whose index is <= Extent[0], or
    # >= Extent[1]:
    if InterpInd <= Extent[0] or InterpInd >= Extent[1]:
        print(f'Cannot interpolate slice number {InterpInd}')
        return []
    
    else:
        #print(f'Looking for the nearest lower slice to {InterpInd}...')
        
        # Iterate backwards from (InterpInd - 1) to SliceInds[0] 
        # to determine the nearest lowest slice index that 
        # contains a single contour:
        for i in range(InterpInd - 1, SliceInds[0] - 1, -1):
            #print(f'\nContourData[{i}] =', ContourData[i])
            #print(f'\nlen(ContourData[{i}]) = {len(ContourData[i])}')
            
            if len(ContourData[i]) == 1:
                #print(f'\nlen(ContourData[{i}]) == 1')
                #print(f'Appending {i} to BoundingSliceInds')
                
                BoundingSliceInds.append(i)
                
                FoundLowerSlice = True
                
                break
                #continue
                
        #print(f'Looking for the nearest upper slice to {InterpInd}...')
        
        # Iterate forwards from (InterpInd + 1) to SliceInds[-1]
        # to determine the nearest upper slice index that contains
        # a single contour:
        for i in range(InterpInd + 1, SliceInds[-1] + 1):
            #print(f'\nContourData[{i}] =', ContourData[i])
            #print(f'\nlen(ContourData[{i}]) = {len(ContourData[i])}')
            #print(f'Appending {i} to BoundingSliceInds')
            
            if len(ContourData[i]) == 1:
                #print(f'\nlen(ContourData[{i}]) == 1')
                #print(f'Appending {i} to BoundingSliceInds')
                
                BoundingSliceInds.append(i)
                
                FoundUpperSlice = True
                
                break
        
        # Deal with case where no lower, or upper, or both upper nor 
        # lower slice index can be used for interpolation:
        #if len(BoundingSliceInds) < 2:
        #    print('Cannot find suitable lower and/or upper slice to ',
        #          f'interpolate slice number {InterpInd}.')
        #    return []
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
    
    print(f'len(Contour) = {len(Contour)}')
    
    # Initialise the cumulative perimeter and list of cumulative perimeters:
    CumP = 0
    CumPerim = []
    
    # Loop through each point:
    for i in range(len(Contour) - 1):
        Point1 = Contour[i]
        Point2 = Contour[i+1]
        
         # The segment length between Point1 and Point2:
        SegmentL = (\
                   (Point1[0] - Point2[0])**2 \
                   + (Point1[1] - Point2[1])**2 \
                   + (Point1[2] - Point2[2])**2 \
                   )**(1/2)
        
        CumP = CumP + SegmentL
              
        CumPerim.append(CumP + SegmentL)
        
    return CumPerim





def NormaliseCumulativePerimeters(CumPerimeters): 
    # Below is compact but might not be understood by a non-Python user:
    #return [CumPerimeters[i]/CumPerimeters[-1] for i in range(len(CumPerimeters))]
    
    # Initialise the list:
    CumPerimNorm = []
    
    for i in range(len(CumPerimeters)):
        CumPerimNorm.append(CumPerimeters[i]/CumPerimeters[-1])
        
    return CumPerimNorm



def GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm, NumNodesToAdd):
    # The inter-node spacing is:
    dP = 1 / (NumNodesToAdd - 1)
    
    """ Exclude the first (0) and last item (1) normalised cumulative perimeter
    since they are already covered by the original normalised cumulative 
    perimeters (i.e. iterate from 1 to NumNodesToAdd - 1)
    """
    # Initialise the normalised cumulative super-sampled perimeters:
    CumSSPerimNorm = [dP]
    
    # First append the cumulative perimeters for the nodes that ensure
    # the minimum inter-node spacing required:
    for i in range(1, NumNodesToAdd - 1):
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
    
    # The first NumNodesToAdd - 2 nodes were additional nodes
    # added to achieve the required minimum inter-node spacing:
    for i in range(NumNodesToAdd - 2):
        IsOrigNode.append(False)
        
    # The next len(Contour) nodes are original nodes from Contour:
    for i in range(len(Contour)):
        IsOrigNode.append(True)
        
        
    return IsOrigNode
    
    

def GetNodesPerSegment(CumSuperSampledPerim, IsOrigNode):
    # Get the indices that sort the cumulative super-sampled perimeters:
    Inds = [i[0] for i in sorted(enumerate(CumSuperSampledPerim), key=lambda x:x[1])]
    
    # Order the boolean list:
    IsOrigNodeSorted = [IsOrigNode[i] for i in Inds]
    
    # Get ordered indices of original nodes:
    OrigIndsSorted = []
    
    for i in range(len(Inds)):
        if IsOrigNodeSorted[i]:
            OrigIndsSorted.append(Inds[i])
        

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
    
    # Initialise the list of the number of nodes in each segment:
    NodesPerSegment = []
    
    for i in range(len(OrigIndsSorted) - 1):
        NodesPerSegment.append(OrigIndsSorted[i + 1] - OrigIndsSorted[i])
        
        
    return NodesPerSegment



    
    
def SuperSampleContour(Contour, NodesPerSegment):
    # Initialise list of super-sampled contour points and the list of booleans 
    # that indicate if the point/node in the super-sampled list is original 
    # (False if the node was interpolated):
    SSContour = []
    SSIsOrigNode = []
    
    #print(f'len(Contour)         = {len(Contour)}')
    #print(f'len(NodesPerSegment) = {len(NodesPerSegment)}')
    
    # Iterate through each contour point without returning to the 0th point 
    # (which is the final node in the list):
    for i in range(len(Contour) - 1):
        #print(i)
        
        # Add the original node:
        SSContour.append(Contour[i])
        SSIsOrigNode.append(True)
        
        # Add the linearly interpolated nodes:
        dx = (Contour[i + 1][0] - Contour[i][0]) / (NodesPerSegment[i] + 1)
        dy = (Contour[i + 1][1] - Contour[i][1]) / (NodesPerSegment[i] + 1)
        
        # Add dx and dy to the x and y components of the last item in SSContour
        # (keep the z component unchanged):
        for j in range(NodesPerSegment[i] - 1):
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

def IntegrateLineLengths(Contour1, Contour2, StartInd2):
    
    # Initialise the cummulative line lengths:
    LineLength = 0
    
    for i in range(len(Contour1)):
        Point1 = Contour1[i]
        
        for j in range(len(Contour2)):
            # The index of the point for Contour2 starts
            # at StartInd2 and increments to N, then is 
            # "wrapped" back to 0 and so on until N 
            # iterations:
            if StartInd2 + j < len(Contour2):
                Ind2 = StartInd2 + j
            else:
                Ind2 = len(Contour2) - (StartInd2 + j)
                
            Point2 = Contour2[Ind2]
            
            # The vector length between Point1 and Point2:
            VectorL = ((Point1[0] - Point2[0])**2 \
                       + (Point1[1] - Point2[1])**2 \
                       + (Point1[2] - Point2[2])**2)**(1/2)
            
            # Add the line length between Point1 and Point2:
            LineLength = LineLength + VectorL
            
            
    return LineLength




""" 
Using the function AddLineLengths(), loop through each N possible 
combinations of starting index StartInd2, and work out which
starting index yields the minimum cumulative line length.
"""

def FindIndForMinCumLength(Contour1, Contour2):
    
    # Create list to store N cumulative line lengths for all
    # N combinations of starting index for Contour2:
    CumLineLengths = []

    # Loop through all N possible starting indeces for Contour2:
    for i in range(len(Contour2)):
        CumLineLengths.append(IntegrateLineLengths(Contour1=Contour1, 
                                                   Contour2=Contour2, 
                                                   StartInd2=i))

    
    # Return the index that yields the minimum cumulative line length:
    return CumLineLengths.index(min(CumLineLengths))





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
    
    # Shift Contour by starting at the index StartInd then iterating through
    # all N subsequent points, wrapping the indeces past N:
    for i in range(len(Contour)):
        # The index of the point for Contour starts at StartInd and increments 
        # to N, then is "wrapped" back to 0 and so on until N iterations:
        if Shift + i < len(Contour):
            Ind = Shift + i
        else:
            Ind = len(Contour) - (Shift + i)

        ShiftedContour.append(Contour[Ind])
        
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


def InterpolateBetweenContours(OverSampledContour1, OverSampledContour2, 
                               OverSampledIsOrigNode1, OverSampledIsOrigNode2,
                               FracInterpSliceInd, Contour1HasMorePts):
    # Initialise the interpolated contour:
    InterpolatedContour = []
    
    if Contour1HasMorePts:
        inds = OverSampledIsOrigNode1
    else:
        inds = OverSampledIsOrigNode2
        
    for i in range(len(OverSampledContour1)):
        if inds[i]:
            InterpX = (1 - FracInterpSliceInd) * OverSampledContour1[i][0] \
                      + FracInterpSliceInd * OverSampledContour2[i][0]
            
            InterpY = (1 - FracInterpSliceInd) * OverSampledContour1[i][1] \
                      + FracInterpSliceInd * OverSampledContour2[i][1]
                      
            InterpZ = (1 - FracInterpSliceInd) * OverSampledContour1[i][2] \
                      + FracInterpSliceInd * OverSampledContour2[i][2]
            
            InterpolatedContour.append([InterpX, InterpY, InterpZ])
            
    return InterpolatedContour