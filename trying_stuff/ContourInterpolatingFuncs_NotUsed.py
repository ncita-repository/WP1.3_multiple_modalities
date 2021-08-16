# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 14:25:16 2020

@author: ctorti



Functions from ContourInterpolatingFuncs.py that are not used (including older
versions of functions)

"""





#def TransformInterpolatedContours(OSContour1, OSContour2, InterpContour, 
#                                  InterpData, FixedDicomDir, MovingDicomDir):
def TransformInterpolatedContours_v1(InterpData, FixedDicomDir, MovingDicomDir):
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
    
    
    Note:
        This function loops through the three contours Contour1, InterpContour 
        and Contour2 for each row in InterpData (i.e. 3*4 = 12 iterations), 
        each time creating inputpoints.txt, transforming the points, then
        parsing outputpoints.txt.  
        
        This is very inefficient due to the numerous file I/O operations.  So
        it should be modified to gather all points into one list, transform 
        all of them, parse all of them, then re-group them into the original
        list indices.
        
        Also, Contour2 on row 0 is equivalent to Contour1 on row 1, 
        Contour2 on row 1 is equivalent to Contour1 on row 2, etc. so there's 
        redundancy.
        
        Could instead gather all points in all Contour1 and InterpContour and
        Contour2 on the last row of InterpData.
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
            #TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)
            TransformPoints(Points=Points, Image=MovIm, 
                            TransformFilter=ElastixImFilt)
    
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








def GetIntersectingPointsInMovingPlanes_v1(InterpData, MovingDicomDir, 
                                           UseInterp, MustBeOnLineSeg):
    """
    Get the list of intersection points of the lines that join each 
    over-sampled-to-interpolated contour #1 to each imaging plane, and likewise 
    for each interpolated-to-over-sampled contour #2.
    
    Inputs:
        InterpData      - A dictionary containing interpolated data.
        
        MovingDicomDir  - Directory containing the DICOMs for the Moving image
                          domain.
                         
        UseInterp       - Boolean that determines whether interpolated contours
                          will be used (True) when finding intersecting points
                          or not (False).
                         
        MustBeOnLineSeg - Boolean value that determines whether the 
                          intersection point must lie on the line segment 
                          (True) or not (False)
        
        
    Returns:
        MovContourData  - A dictionary containing the contour data belonging to
                          the Moving image domain.
                          
                          
    Note:
        If UseInterp = True, there may be two intersection points for any given
        imaging plane for each of P points in the contours:  One intersection
        from the line segment created by joining a point from Contour1 to the
        interpolated contour, and another intersection from the line segment
        created by joining a point from the interpolated contour to Contour2.
        
        This probably explains some odd looking features in the collection of
        intersecting points.
        
        In v2 of this function of the two intersecting points, the one which is
        furthest to the imaging plane is rejected.
        
        The other big difference with v2 is that v1 returns MovContourData, 
        a list of a list of points for every slice in the Moving image domain.
        By contrast, v2 adds additional keys to InterpData.
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
            # intersecting points of Point0->Point1 and Point1->Point2 with 
            # each plane:
            for k in range(MovDims[2]):
    
                # Need a point lying on the plane - use the plane's origin:
                PlaneOrigin = [MovOrigin[0] + k*MovSpacings[2]*PlaneNormal[0], 
                               MovOrigin[1] + k*MovSpacings[2]*PlaneNormal[1], 
                               MovOrigin[2] + k*MovSpacings[2]*PlaneNormal[2]
                               ]
    
                #print(f'Length of PlaneNormal is {GetVectorLength(PlaneNormal)}')
                
                if UseInterp:
                    IntersPoint0, \
                    OnLineSeg0 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                          PlaneP=PlaneOrigin,
                                                          LineP0=Point0,
                                                          LineP1=Point1,
                                                          MustBeOnLineSeg=MustBeOnLineSeg)
                
                    IntersPoint1, \
                    OnLineSeg1 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                          PlaneP=PlaneOrigin,
                                                          LineP0=Point1,
                                                          LineP1=Point2,
                                                          MustBeOnLineSeg=MustBeOnLineSeg)
                    
                else:
                    IntersPoint0, \
                    OnLineSeg0 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                          PlaneP=PlaneOrigin,
                                                          LineP0=Point0,
                                                          LineP1=Point2,
                                                          MustBeOnLineSeg=MustBeOnLineSeg)
                
                
                #print('PlaneNormal =', PlaneNormal)
                #print('PlaneOrigin =', PlaneOrigin)
                #print('Point0 =', Point0)
                #print('Point1 =', Point1)
                #print('Point2 =', Point2)
                #print('IntersPoint0 =', IntersPoint0, 'OnLineSeg0 =', OnLineSeg0)
                #print('IntersPoint1 =', IntersPoint1, 'OnLineSeg1 =', OnLineSeg1, '\n')
                
                if 0:
                    if IntersPoint0 or IntersPoint1:
                        print(f'i = {i}, j = {j}, k = {k} \n')
                        
                        if IntersPoint0:
                            print('IntersPoint0 =', IntersPoint0)
                            print('OnLineSeg0 =', OnLineSeg0, '\n')
                            
                        if IntersPoint1:
                            print('IntersPoint1 =', IntersPoint1)
                            print('OnLineSeg1 =', OnLineSeg1, '\n')
                            
                if 0:
                    if IntersPoint0 and IntersPoint1:
                        print(f'i = {i}, j = {j}, k = {k} \n')
                        
                        print('IntersPoint0 =', IntersPoint0)
                        print('OnLineSeg0 =', OnLineSeg0, '\n')
                        
                        print('IntersPoint1 =', IntersPoint1)
                        print('OnLineSeg1 =', OnLineSeg1, '\n')
                    
                    
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
            
    #print('IntersContours', IntersContours, '\n')
    
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






def GetIntersectingPointsInMovingPlanes_v2(InterpData, MovingDicomDir, 
                                           UseInterp, MustBeOnLineSeg,
                                           OnlyKeepClosestIntersPt):
    """
    Get the list of intersection points of the lines that join each 
    over-sampled-to-interpolated contour #1 to each imaging plane, and likewise 
    for each interpolated-to-over-sampled contour #2.
    
    Inputs:
        InterpData      - A dictionary containing interpolated data.
        
        MovingDicomDir  - Directory containing the DICOMs for the Moving image
                          domain.
                         
        UseInterp       - Boolean that determines whether interpolated contours
                          will be used (True) when finding intersecting points
                          or not (False).
                         
        MustBeOnLineSeg - Boolean value that determines whether the 
                          intersection point must lie on the line segment 
                          (True) or not (False)
                         
        OnlyKeepClosestIntersPt - Boolean value that determines whether or not the
                              intersection point whose midpoint is furthest
                              from the imaging plane will be rejected.
        
        
    Returns:
        InterpData      - Updated with additional keys 'AllIntPtsAllSlices' and
                          'AllIntDtsAllSlices'.
                          
    
    
    
    Note:
        v2 is meant to overcome a limitation of v1 whereby for any given
        triplet of points from Contour1, the interpolated contour and Contour2,
        there may be two intersection points with any given plane - one from
        the intersection of the line segment created by joining the point from
        Contour1 and the interpolated contour, and one from the intersection of
        the line segment created by joining the point from the interpolated 
        contour and Contour2.
        
        In this version, of the two intersection points, the one furthest from
        the imaging plane is rejected.
        
        The other big difference with v1 is that v1 returns MovContourData, 
        a list of a list of points for every slice in the Moving image domain.
        By contrast, v2 adds additional keys to InterpData: 
        'AllIntPtsAllSlices' and 'AllIntDtsAllSlices'.  
        For P points within the contours Contour1, InterpContour and Contour2
        for any given row in InterpData, 'AllIntPtsAllSlices' contains a list
        of length P (number of contour points) of a list of length S (number of
        slices) of the list of [x, y, z] coordinates.
            
        Likewise, 'AllIntDtsAllSlices' contains a list of a list of a list of
        distances between the plane and the midway point of the two points 
        whose line segment's intersection with the plane defined the point in
        the corresponding location in 'AllIntPtsAllSlices'.
        
        
    Comment 01/09/2020:
        
        I'm not sure it makes any sense to concern myself with which of two
        intersection points is closer to the midpoints of the line segments 
        because I no longer believe that two intersection points will arise
        that are both on their respective line segments.  So should revert to
        v1 of this function.
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
    
    # Initialise a list of points and distances for each imaging plane:
    AllIntersPoints = [[] for i in range(MovDims[2])]
    AllIntersDistances = [[] for i in range(MovDims[2])]
    
    #AllIntPtsAllSlices = []
    #AllIntDtsAllSlices = []
    
    InterpData.update({'AllIntPtsAllSlices' : [],
                       'AllIntDtsAllSlices' : []
                       })
    
    
    # Loop through each set of interpolation data in InterpData:
    for i in range(len(InterpData['InterpSliceInd'])):
        
        AllIntPtsAllSlices = []
        AllIntDtsAllSlices = []
    
        # Loop through each point in Contour1, InterpContour and Contour2:
        for p in range(len(InterpData['MovOSContour1Pts'][i])):
    
            Point0 = InterpData['MovOSContour1Pts'][i][p]
            Point1 = InterpData['MovInterpContourPts'][i][p]
            Point2 = InterpData['MovOSContour2Pts'][i][p]
            
            # Initialise the intersection points at all slices:
            IntPtsAllSlices = []
            IntDtsAllSlices = []
    
            # Loop through all imaging planes in Moving Image and get the 
            # intersecting points of Point0->Point1 and Point1->Point2 with 
            # each plane:
            for s in range(MovDims[2]):
    
                # Need a point lying on the plane - use the plane's origin:
                PlaneOrigin = [MovOrigin[0] + s*MovSpacings[2]*PlaneNormal[0], 
                               MovOrigin[1] + s*MovSpacings[2]*PlaneNormal[1], 
                               MovOrigin[2] + s*MovSpacings[2]*PlaneNormal[2]
                               ]
    
                #print(f'Length of PlaneNormal is {GetVectorLength(PlaneNormal)}')
                
                if UseInterp:
                    IntersPoint0, \
                    OnLineSeg0 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                          PlaneP=PlaneOrigin,
                                                          LineP0=Point0,
                                                          LineP1=Point1,
                                                          MustBeOnLineSeg=MustBeOnLineSeg)
                
                    IntersPoint1, \
                    OnLineSeg1 = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                          PlaneP=PlaneOrigin,
                                                          LineP0=Point1,
                                                          LineP1=Point2,
                                                          MustBeOnLineSeg=MustBeOnLineSeg)
                    
                    #print('PlaneNormal =', PlaneNormal)
                    #print('PlaneOrigin =', PlaneOrigin)
                    #print('Point0 =', Point0)
                    #print('Point1 =', Point1)
                    #print('Point2 =', Point2)
                    #print('IntersPoint0 =', IntersPoint0, 'OnLineSeg0 =', OnLineSeg0)
                    #print('IntersPoint1 =', IntersPoint1, 'OnLineSeg1 =', OnLineSeg1, '\n')
                        

                    # Get the absolute z-distance between the z-coordinate of
                    # IntersPoint0 and the midway z-position of the points used
                    # to define that line (i.e. the z-components of Point0 and 
                    # Point1):
                    midZ0 = ( Point1[2] + Point0[2] ) / 2
                    midZ1 = ( Point2[2] + Point1[2] ) / 2
                            
                    if IntersPoint0:
                        d0 = abs( midZ0 - IntersPoint0[2] )
                        
                    if IntersPoint1:
                        d1 = abs( midZ1 - IntersPoint1[2] )
                        
                    
                    """ 
                    Note:
                        I could use the boolean values OnLineSeg0 and 
                        OnLineSeg1 to retain or reject IntersPoint0 and
                        IntersPoint1, but one will only be rejected if 
                        either OnLineSeg0 or OnLineSeg1 is False.
                        
                        The benefit of using d0 and d1 is that one of these
                        absolute distances will almost certainly be smaller
                        than the other.
                    """
                        
                    if IntersPoint0 and IntersPoint1:
                        if d0 < d1:
                            # IntersPoint0 is closer to the imaging plane, so
                            # use this point:
                            IntersPoint = IntersPoint0
                            IntersDistance = d0
                        else:
                            # IntersPoint1 is closer to the imaging plane:
                            IntersPoint = IntersPoint1
                            IntersDistance = d1
                    
                    elif IntersPoint0 and not IntersPoint1:
                        IntersPoint = IntersPoint0
                        #IntersDistance = abs( midZ0 - IntersPoint0[2] )
                        IntersDistance = d0
                        
                    elif not IntersPoint0 and IntersPoint1:
                        IntersPoint = IntersPoint1
                        #IntersDistance = abs( midZ1 - IntersPoint1[2] )
                        IntersDistance = d1
                        
                    else:
                        IntersPoint = []
                        #IntersDistance = None
                        IntersDistance = []
                            
                        
                    
                else:
                    IntersPoint, \
                    OnLineSeg = GetLinePlaneIntersection(PlaneN=PlaneNormal, 
                                                         PlaneP=PlaneOrigin,
                                                         LineP0=Point0,
                                                         LineP1=Point2,
                                                         MustBeOnLineSeg=MustBeOnLineSeg)
                    
                    if IntersPoint:
                        midZ = ( Point2[2] + Point0[2] ) / 2
                    
                        IntersDistance = abs( midZ - IntersPoint[2] )
                        
                    else:
                        IntersPoint = []
                        #IntersDistance = None
                        IntersDistance = []

                
                #print(f'IntersPoint = {IntersPoint}')
                #print(f'IntersDistance = {IntersDistance}')
                if OnlyKeepClosestIntersPt or ( not UseInterp ):
                    AllIntersPoints[s].append(IntersPoint)
                    AllIntersDistances[s].append(IntersDistance)
                    
                    IntPtsAllSlices.append(IntersPoint)
                    IntDtsAllSlices.append(IntersDistance)
                    
                else:
                    #print('THIS WAS TRUE')
                    # Append both IntersPoint0 and IntersPoint1 if they exist: 
                    if IntersPoint0:
                        #print('FIRST IF')
                        AllIntersPoints[s].append(IntersPoint0)
                        AllIntersDistances[s].append(d0)
                        
                        IntPtsAllSlices.append(IntersPoint0)
                        IntDtsAllSlices.append(d0)
                        
                    if UseInterp and IntersPoint1:
                        #print('SECOND IF')
                        AllIntersPoints[s].append(IntersPoint1)
                        AllIntersDistances[s].append(d1)
                        
                        IntPtsAllSlices.append(IntersPoint1)
                        IntDtsAllSlices.append(d1)
                    
                    
                    
            AllIntPtsAllSlices.append(IntPtsAllSlices)
            AllIntDtsAllSlices.append(IntDtsAllSlices)
        
        InterpData['AllIntPtsAllSlices'].append(AllIntPtsAllSlices)
        InterpData['AllIntDtsAllSlices'].append(AllIntDtsAllSlices)  
        
        #import pandas as pd
        #pd.DataFrame(InterpData)
        #print('\n\n', InterpData['AllIntPtsAllSlices'], '\n\n')
        #print('\n\n', InterpData['AllIntDtsAllSlices'], '\n\n')

    return InterpData




def GetClosestIntersectionPoints(InterpData, MaxDistToImagePlane):
    
    #InterpData.update({'MinIntPts' : [],
    #                   'MinIntDts' : [],
    #                   'MinIntSlice' :
    #                   })
    
    # Initialise the list of the closest intersecting points, their distances 
    # to the nearest imaging plane, and the list slice indices of the nearest 
    # imaging plane for all point pairs:
    AllClosestPts = []
    AllClosestDts = []
    AllClosestInds = []
    
    # Loop through each set of interpolation data in InterpData:
    for i in range(len(InterpData['InterpSliceInd'])):
        
        AllIntPtsAllSlices = InterpData['AllIntPtsAllSlices'][i]
        AllIntDtsAllSlices = InterpData['AllIntDtsAllSlices'][i]
    
        # Initialise the list of the closest intersecting points, their 
        # distances to the nearest imaging plane, and the list slice indices of 
        # the nearest imaging plane for this point pair:
        ClosestDts = []
        ClosestPts = []
        ClosestInds = []
            
        # Loop through each list of distances for each point pair:
        for p in range(len(AllIntDtsAllSlices)):
            # Get the intersection point and distances for all slices for this 
            # point pair:
            IntPtsAllSlices = AllIntPtsAllSlices[p]
            IntDtsAllSlices = AllIntDtsAllSlices[p]
            
            # Initialise the list of distances for all slices for this point 
            # pair:
            Distances = []
            
            # Loop through the list of distances for each slice:
            for s in range(len(IntDtsAllSlices)):
                # Append to Distances if not [], and 100 (arbitrary) if empty 
                # (to maintain the relationship between the distances and the 
                # slice num they relate to):
                if IntDtsAllSlices[s]:
                    Distances.append(IntDtsAllSlices[s])
                else:
                    Distances.append(100)
                    
            # The index of the closest point:
            ind = Distances.index(min(Distances))
            
            ClosestDt = Distances[ind] 
            
            # Only append if the closest distance is < MaxDistToImagePlane:
            if ClosestDt < MaxDistToImagePlane:
                ClosestDts.append(Distances[ind])
                ClosestPts.append(IntPtsAllSlices[ind])
                ClosestInds.append(ind)
                
        AllClosestDts.append(ClosestDts)
        AllClosestPts.append(ClosestPts)
        AllClosestInds.append(ClosestInds)
    
    InterpData.update({'ClosestIntersPts'  : AllClosestPts,
                       'ClosestIntersDts'  : AllClosestDts,
                       'ClosestIntersInds' : AllClosestInds
                       })    
        
    return InterpData




def GroupIntersPtsBySliceIndAndContourNum(InterpData, MovingDicomDir):
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
    
    # Initialise a list of interesting points and distances for each imaging 
    # plane:
    #IntersPts = [ [] for i in range(MovDims[2]) ]
    #IntersDts = [ [] for i in range(MovDims[2]) ]
    #IntersPts = [ [[]] for i in range(MovDims[2]) ]
    #IntersDts = [ [[]] for i in range(MovDims[2]) ]
    PtsPCSBySliceByContour = [ [] for i in range(MovDims[2]) ]
    DtsBySliceByContour = [ [] for i in range(MovDims[2]) ]
    ISIBySliceByContour = [ [] for i in range(MovDims[2]) ]
    PtsICSBySliceByContour = [ [] for i in range(MovDims[2]) ]
    CTBySliceByContour = [ [] for i in range(MovDims[2]) ]

    
    # Loop through each set of interpolation data in InterpData:
    for i in range(len(InterpData['InterpSliceInd'])):
        
        Points = InterpData['ClosestIntersPts'][i]
        Distances = InterpData['ClosestIntersDts'][i]
        Inds = InterpData['ClosestIntersInds'][i]
        
        # Get the interpolated slice index for this i (to track where the
        # re-grouped point originated from): 
        InterpSliceInd = InterpData['InterpSliceInd'][i]
        
        PtsBySlice = [ [] for i in range(MovDims[2]) ]
        DtsBySlice = [ [] for i in range(MovDims[2]) ]
        ISIBySlice = [ [] for i in range(MovDims[2]) ]
        
        # Loop through each slice index:
        for j in range(len(Inds)):
            #IntersPts[Inds[j]].append(Points[j])
            #IntersDts[Inds[j]].append(Distances[j])
            
            #IntersPts[Inds[j]][i].append(Points[j])
            #IntersDts[Inds[j]][i].append(Distances[j])
            
            PtsBySlice[Inds[j]].append(Points[j])
            DtsBySlice[Inds[j]].append(Distances[j])
            ISIBySlice[Inds[j]].append(InterpSliceInd)
         
        
        #PtsBySliceByContour[i].append(PtsBySlice)
        #DtsBySliceByContour[i].append(DtsBySlice)
        #PtsBySliceByContour.append(PtsBySlice)
        #DtsBySliceByContour.append(DtsBySlice)
        
        
        
        # Loop through each slice in PtsBySlice:
        for s in range(MovDims[2]):
            # If the list of points for this slice is not empty:
            if PtsBySlice[s]:
                PtsPCSBySliceByContour[s].append(PtsBySlice[s])
                DtsBySliceByContour[s].append(DtsBySlice[s])
                ISIBySliceByContour[s].append(ISIBySlice[s])
                CTBySliceByContour[s].append(3)
                
                PtsICS = PCStoICS(Pts_PCS=PtsBySlice[s], 
                                  Origin=MovOrigin, 
                                  Directions=MovDirs, 
                                  Spacings=MovSpacings)
        
                PtsICSBySliceByContour[s].append(PtsICS)
                
            #else:
            #    CTBySliceByContour[s].append(0)
    
    
    # Create a list of slice numbers:
    SliceNo = list(range(MovDims[2]))
    
    MovContourData = {
                      'SliceNo'        : SliceNo, 
                      'ContourType'    : CTBySliceByContour,
                      'InterpSliceInd' : ISIBySliceByContour,
                      'PointPCS'       : PtsPCSBySliceByContour,
                      'PointICS'       : PtsICSBySliceByContour,
                      'Distance2Plane' : DtsBySliceByContour
                       }
    
    
    #return IntersPts, IntersDts
    return MovContourData
    #return MovContourData, CTBySliceByContour, PtsPCSBySliceByContour, \
    #       PtsICSBySliceByContour, DtsBySliceByContour, ISIBySliceByContour
    
    
    
    
    
    
    
def TravellingSalesman(Points):
    """
    Find the shortest route to join all points in Points by bruteforce.
    
    Input:
        Points  - A list of [x, y, z] (or [x, y]) coordinates of a list of 
                 points
                 
    Returns:
        Ordered - Points sorted to minimise the contour perimeter
    
    Code modified from:
    https://codereview.stackexchange.com/questions/81865/travelling-salesman-using-brute-force-and-heuristics
    
    
    Note:
        This functions is inappropriate since the time cost is O(N!).
        The permutations for N = 10 points took 1.5 s
        With N = 15 took I had to interrupt the kernel.
    """
    from itertools import permutations
    
    start = Points[0]
    
    Ordered = min([perm for perm in permutations(Points) if perm[0] == start], key=GetContourPerimeter)
    
    return Ordered



