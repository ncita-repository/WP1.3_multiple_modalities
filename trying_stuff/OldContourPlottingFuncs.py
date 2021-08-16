# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:29:22 2020

@author: ctorti
"""



def PlotInterpContours2D_OLD(Contours, Labels, Colours, Shifts, 
                             InterpSliceInd, BoundingSliceInds, dP, 
                             MaxDistToImagePlane, AnnotatePtNums, ExportPlot):
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
        X = []; Y = []
    
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
                + f'_Interp_at_slice_{InterpSliceInd}_of_contours_at_' \
                + f'{BoundingSliceInds[0]}_and_{BoundingSliceInds[1]}_dP_{dP}' \
                + f'_MaxDistToPlane_{MaxDistToImagePlane}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    return





def GetPointsInImagePlanesWithContoursOLD(ContourData, ContourTypeNo, DicomDir):
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






def GetMeshGridsForImagePlanesOLD(ContourData, ContourTypeNo, DicomDir):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes.
    """
    
    from GetImageAttributes import GetImageAttributes
    import numpy as np
    #from scipy.interpolate import griddata
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
        
    # Get the slice numbers of the imaging planes that contain contours:
    SliceNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    Meshgrids = []
    
    for s in SliceNums:
        xMin = Origin[0] + s*Spacings[2]*Dirs[6]
        yMin = Origin[1] + s*Spacings[2]*Dirs[7]
        #zMin = Origin[2] + s*Spacings[2]*Dirs[8]
        
        xMax = xMin + Dims[0]*Spacings[0]*Dirs[0] + Dims[1]*Spacings[1]*Dirs[3]
        yMax = yMin + Dims[0]*Spacings[0]*Dirs[1] + Dims[1]*Spacings[1]*Dirs[4]
        #zMax = zMin + Dims[0]*Spacings[0]*Dirs[2] + Dims[1]*Spacings[1]*Dirs[5]
        
        #X = np.linspace(xMin, xMax, 100)
        X = np.linspace(xMin, xMax, Dims[0])
        #Y = np.linspace(yMin, yMax, 100)
        Y = np.linspace(yMin, yMax, Dims[1])
        #Z = np.linspace(zMin, zMax, 100)
        
        #XX, YY, ZZ = np.meshgrid(X, Y, Z)
        #Meshgrids.append( np.meshgrid(X, Y, Z) )
        
        XX, YY = np.meshgrid(X, Y)
        
        #YZ, ZZ = np.meshgrid(Y, Z)
        #ZZ1, ZZ2 = np.meshgrid(Z, Z)
        #ZZ = Z.reshape(XX.shape)
        #ZZ = griddata((X, Y), Z, (XX, YY), method='cubic')
        
        #ZZ = np.zeros(XX.shape, dtype=float)
        #ZZ = np.zeros((Dims[0], Dims[1]), dtype=float)
        ZZ = np.zeros((Dims[1], Dims[0]), dtype=float)
        
        for i in range(Dims[0]):
            for j in range(Dims[1]):
                #ZZ[i,j] = Origin[2] + i*Spacings[0]*Dirs[2] + j*Spacings[1]*Dirs[5] + s*Spacings[2]*Dirs[8] 
                ZZ[j,i] = Origin[2] + i*Spacings[0]*Dirs[2] + j*Spacings[1]*Dirs[5] + s*Spacings[2]*Dirs[8] 
                
                
        #print('XX.shape =', XX.shape)
        #print('YY.shape =', YY.shape)
        #print('ZZ.shape =', ZZ.shape)
        
        Meshgrids.append([XX, YY, ZZ])
        #Meshgrids.append([XX, YY, ZZ1 + ZZ2])
        #Meshgrids.append([X, Y, Z])
        
        
    return Meshgrids




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
            X = []; Y = []; Z = []
            
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
        
    return
    





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
        X = []; Y = []
    
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
        X = []; Y = []
        
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
        
    return
    




