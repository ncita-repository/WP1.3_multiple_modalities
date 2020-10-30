# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:28:58 2020

@author: ctorti
"""





def GetMeshGridsForImagePlanes(PlaneNums, DicomDir):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes for select plane numbers.
    
    USE GetMeshGridsForImagePlanesInContourData INSTEAD.
    """
    
    from GetImageAttributes import GetImageAttributes
    import numpy as np
    #from scipy.interpolate import griddata
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
    
    Meshgrids = []
    
    for s in PlaneNums:
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


def GetMeshGridsForImagePlanesNEW(DicomDir, Convert2ICS):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes for all image planes.
    
    Note that the points are in the Patient Coordinate System but will be 
    converted to the Image Coordinate System if Convert2ICS=True.
    
    07/09/2020
    """
    
    from GetImageAttributes import GetImageAttributes
    import numpy as np
    #from scipy.interpolate import griddata
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
    
    Meshgrids = []
    
    for s in range(Dims[2]):
        print(f'\ns = {s}')
        print('Origin =', Origin)
        print('s*Spacings[2]*Dirs[6] =', s*Spacings[2]*Dirs[6])
        print('s*Spacings[2]*Dirs[7] =', Origin)
        
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
        
        if Convert2ICS:
            # Import function:
            from PCStoICS import PCStoICS
            
            """ This is not an elegant approach (and it's VERY INEFFICIENT) but 
            a quick means of getting a result. """
            for i in range(Dims[0]):
                for j in range(Dims[1]):
                    PtPCS = [XX[j,i], YY[j,i], ZZ[j,i]]
                    
                    PtICS = PCStoICS(Pts_PCS=PtPCS, Origin=Origin, 
                                     Directions=Dirs, Spacings=Spacings)
                    
                    # Over-write the [j,i]^th elements in XX, YY and ZZ:
                    XX[j,i] = PtICS[0]
                    YY[j,i] = PtICS[1]
                    ZZ[j,i] = PtICS[2]
        
        
        
        Meshgrids.append([XX, YY, ZZ])
        #Meshgrids.append([XX, YY, ZZ1 + ZZ2])
        #Meshgrids.append([X, Y, Z])
        
        
    return Meshgrids






def GetMeshGridsForImagePlanesInContours(Contours, OnlyPlanesWithContours, 
                                         DicomDir):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes from a list of contours.
    """
    
    if OnlyPlanesWithContours:
        # Get the indices of the imaging planes that contain contours:
        PlaneNums = GetIndsOfSliceNumsWithContours(Contours)
    else:
        # Create a list from 0 to len(Contours):
        PlaneNums = list(range(len(Contours)))
    
    Meshgrids = GetMeshGridsForImagePlanes(PlaneNums, DicomDir)
          
    return Meshgrids


def GetMeshGridsForImagePlanesInContourData(ContourData, ContourTypeNo, 
                                            DicomDir):
    """
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes.
    """
        
    # Get the indices of the imaging planes that contain contours of the 
    # required ContourTypeNo:
    PlaneNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    Meshgrids = GetMeshGridsForImagePlanes(PlaneNums, DicomDir)
        
    return Meshgrids





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
    
    
    Note:  
        I wasn't able to get retain a function that can either plot all results
        as subplots or as individual figures whilst maintaining equal, shared
        axes for all subplots.  So the function PlotInterpContours2D was split
        into two functions:  PlotInterpContours2D for plotting individual 
        figures, and PlotInterpContours2DSubPlots for plotting subplots.
        
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
        fig = plt.figure(figsize=(5*Ncols, 5*Nrows), dpi=DPI) # works
        ##fig = plt.figure(figsize=(5*Ncols, 5*Nrows), dpi=DPI, sharex=True, sharey=True)


    # Store the axes:
    #axs = []
    
    for i in range(N):
        # Set up the axes for this sub-plot:
        if SubPlots:
            ax = fig.add_subplot(Nrows, Ncols, i+1) # works
            ##ax = fig.add_subplot(Nrows, Ncols, i+1, sharex=True, sharey=True)
            ##axs.append(fig.add_subplot(Nrows, Ncols, i+1))
            #fig, ax = plt.subplots(Nrows, Ncols, sharex=True, sharey=True, figsize=(5*Ncols, 5*Nrows), dpi=DPI)
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
            X = []; Y = []
        
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
            #ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]);    
            plt.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]); # works
            ax.legend(loc='upper left', fontsize='small') # works
            # Plot dots:
            if True:
                plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[j]); # works
            
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
                    plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[j]); # works
            # Also annotate the last point:
            #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
        
        plt.axis('equal')
        
        # Equalise the axes range for all plots:
        #xRangesMin = xRangesMax = yRangesMin = yRangesMax = []
        
        # Equalise the axes ticks for all plots:
        #OldxTicksMins = OldxTicksMaxs = OldyTicksMins = OldyTicksMaxs = []

        #for ax in axs:
            #print('xlim =', ax.get_xlim())
            #print('ylim =', ax.get_ylim(), '\n\n')
            #xRangesMin.append(ax.get_xlim()[0])
            #xRangesMax.append(ax.get_xlim()[1])
            #yRangesMin.append(ax.get_ylim()[0])
            #yRangesMax.append(ax.get_ylim()[1])
            #print('xticks =', ax.get_xticks())
            #print('yticks =', ax.get_yticks(), '\n\n')
            #OldxTicksMins.append(ax.get_xticks()[0])
            #OldxTicksMaxs.append(ax.get_xticks()[-1])
            #OldyTicksMins.append(ax.get_yticks()[0])
            #OldyTicksMaxs.append(ax.get_yticks()[-1])

            
        #print('OldxTicks =', OldxTicks, '\n\n')
        #print('OldxTicks[0] =', OldxTicks[0], '\n\n')
        #print('OldxTicks[0][0] =', OldxTicks[0][0], '\n\n')
        #print('OldxTicks[0][-1] =', OldxTicks[0][-1], '\n\n')
        
        #xRangeMinMax = [min(xRangesMin), max(xRangesMax)]
        #yRangeMinMax = [min(yRangesMin), max(yRangesMax)]
        #OldxTicksMin = min(OldxTicksMins)
        #OldxTicksMax = max(OldxTicksMaxs)
        #OldyTicksMin = min(OldyTicksMins)
        #OldyTicksMax = max(OldyTicksMaxs)
        
        #print('OldxTicksMin =', OldxTicksMin)
        #print('OldxTicksMax =', OldxTicksMax)
        #print('OldyTicksMin =', OldyTicksMin)
        #print('OldyTicksMax =', OldyTicksMax, '\n\n')
        
        #print('xRangesMin =', xRangesMin)
        #print('xRangesMax =', xRangesMax)
        #print('yRangesMin =', yRangesMin)
        #print('yRangesMin =', yRangesMin)
        #print('xRangeMinMax =', xRangeMinMax)
        #print('yRangeMinMax =', yRangeMinMax)
        
        # Create a new set of ticks:
        #NewxTicks = list(range(int(OldxTicksMin), int(OldxTicksMax), 5))
        #NewyTicks = list(range(int(OldyTicksMin), int(OldyTicksMax), 5))
        
        #for ax in axs:
        #    ax.set_xlim(xRangeMinMax)
        #    ax.set_ylim(yRangeMinMax)
            #ax.set_xticks(NewxTicks)
            #ax.set_yticks(NewyTicks)
    
        plt.title('Interpolation between ' + FixedOrMoving \
                  + f' slices {SliceInd1} and {SliceInd2}')
        
        if not SubPlots:
            # Create filename for exported figure:
            FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                       + '_Interp_bt_' + FixedOrMoving \
                       + f'_slices_{SliceInd1}_and_{SliceInd2}__dP_{dP}.png'
            
            if ExportPlot:
                plt.savefig(FigFname, bbox_inches='tight')
    
    FirstSlice = InterpData['BoundingSliceInds'][0][0]
    LastSlice = InterpData['BoundingSliceInds'][-1][1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                   + '_Interp_bt_' + FixedOrMoving \
                   + f'_slices_{FirstSlice}_and_{LastSlice}__dP_{dP}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return



def PlotInterpContours2DSubPlots(InterpData, FixedOrMoving, dP, 
                                 AnnotatePtNums, SubPlots, ExportPlot):
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
    
    
    Note:  
        I wasn't able to get retain a function that can either plot all results
        as subplots or as individual figures whilst maintaining equal, shared
        axes for all subplots.  So the function PlotInterpContours2D was split
        into two functions:  PlotInterpContours2D for plotting individual 
        figures, and PlotInterpContours2DSubPlots for plotting subplots.
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
        
    fig, axs = plt.subplots(Nrows, Ncols, sharex=True, sharey=True, figsize=(5*Ncols, 5*Nrows), dpi=DPI)
    
    # Initialise sub-plot iterator:
    i = 0
    
    for r in range(Nrows):
        for c in range(Ncols):
            # Set up the axes for this sub-plot:
            if False:
                if SubPlots:
                    #ax = fig.add_subplot(Nrows, Ncols, i+1) # works
                    ##ax = fig.add_subplot(Nrows, Ncols, i+1, sharex=True, sharey=True)
                    ##axs.append(fig.add_subplot(Nrows, Ncols, i+1))
                    #fig, ax = plt.subplots(Nrows, Ncols, sharex=True, sharey=True, figsize=(5*Ncols, 5*Nrows), dpi=DPI)
                    continue
                else:
                    fig = plt.figure(figsize=(10, 10), dpi=DPI)
                    ax = fig.add_subplot(111)
                    #axs.append(fig.add_subplot(111))
            
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
                X = []; Y = []
            
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
                #ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]);    
                #plt.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]); # works
                axs[r, c].plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[j], label=AllLabels[j]);   
                #ax.legend(loc='upper left', fontsize='small') # works
                axs[r, c].legend(loc='upper left', fontsize='small') # works
                #axs[-1].legend(loc='upper left', fontsize='small')
                # Plot dots:
                if True:
                    #plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[j]); # works
                    axs[r, c].plot(X, Y, '.', markersize=MarkerSize, c=Colours[j]);
                
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
                        #plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[j]); # works
                        axs[r, c].text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[j]);
                # Also annotate the last point:
                #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
            
            #plt.axis('equal')
            #axs[r, c].axis('equal')
            axs[r, c].set_aspect('equal', adjustable='box')
            #axs[r, c].set_aspect(1.0/axs[r, c].get_data_ratio(), adjustable='box')
        
            #plt.title('Interpolation between ' + FixedOrMoving \
            #          + f' slices {SliceInd1} and {SliceInd2}')
            
            axs[r, c].set_title('Interpolation between ' + FixedOrMoving \
                                + f' slices {SliceInd1} and {SliceInd2}')
            
            if not SubPlots:
                # Create filename for exported figure:
                FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                           + '_Interp_bt_' + FixedOrMoving \
                           + f'_slices_{SliceInd1}_and_{SliceInd2}__dP_{dP}.png'
                
                if ExportPlot:
                    plt.savefig(FigFname, bbox_inches='tight')
                   
            # Increment sub-plot iterator:
            i += 1 
    
    FirstSlice = InterpData['BoundingSliceInds'][0][0]
    LastSlice = InterpData['BoundingSliceInds'][-1][1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                   + '_Interp_bt_' + FixedOrMoving \
                   + f'_slices_{FirstSlice}_and_{LastSlice}__dP_{dP}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return




def GetPointsInImagePlanes(DicomDir):
    """
    Create a closed "contour" of the 4 (+ 1 to close the contour) points that
    define the extent of the imaging planes.
    """
    
    from GetImageAttributes import GetImageAttributes
        
    package = 'pydicom'
    
    # Get the image attributes:
    Origin, Dirs, Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                                      Package=package)
    
    # Create "contours" that mark out corners that define each plane in 
    # SliceNums:
    Planes = []
    
    for s in range(Dims[2]):
        
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




def GetPointsInImagePlanes_v2(Origin, Dirs, Spacings, Dims):
    """
    Create a closed "contour" of the 4 (+ 1 to close the contour) points that
    define the extent of the imaging planes.
    
    Adapted from GetPointsInImagePlanes but uses image attributes as inputs
    rather than DicomDir.
    """
    
    # Create "contours" that mark out corners that define each plane in 
    # SliceNums:
    Planes = []
    
    for s in range(Dims[2]):
        
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





def GetPointsInImagePlanesWithContours(ContourData, ContourTypeNo, DicomDir):
    """
    Create a closed "contour" of the 4 (+ 1 to close the contour) points that
    define the extent of the imaging planes that contain contour points only.
    """
    
    # Get the list of the list of the four points at the corners of the image
    # plane for all imaging planes:
    AllPlanes = GetPointsInImagePlanes(ContourData, ContourTypeNo, DicomDir)
        
    # Get the slice numbers of the imaging planes that contain contours:
    SliceNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    # Create "contours" that mark out corners that define each plane in 
    # SliceNums:
    Planes = []
    
    for i in range(len(SliceNums)):
        s = SliceNums[i]
        
        Planes.append(AllPlanes[s])
        
    return Planes





#def GetGridOfPointsInImagePlanes(ContourData, ContourTypeNo, DicomDir):
def GetGridOfPointsInImagePlanes(DicomDir):
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
    #SliceNums = GetIndsOfSliceNumsOfContourType(ContourData, ContourTypeNo)
    
    SliceNums = list(range(Dims[2]))
    
    # Initialise the list grid of points for each plane:
    GridPtsPerPlane = []
    #XCoordsPerPlane = []
    #YCoordsPerPlane = []
    #ZCoordsPerPlane = []
    
    for s in SliceNums:
        
        Points = []
        #X = []
        #Y = []
        #Z = []
        
        for i in range(Dims[0]):
            for j in range(Dims[1]):
                x = Origin[0] + i*Spacings[0]*Dirs[0] + j*Spacings[1]*Dirs[3] + s*Spacings[2]*Dirs[6]
                
                y = Origin[1] + i*Spacings[0]*Dirs[1] + j*Spacings[1]*Dirs[4] + s*Spacings[2]*Dirs[7]
                
                z = Origin[2] + i*Spacings[0]*Dirs[2] + j*Spacings[1]*Dirs[5] + s*Spacings[2]*Dirs[8]
                
                Points.append([x, y, z])
                #X.append(x)
                #Y.append(y)
                #Z.append(z)
          
        GridPtsPerPlane.append(Points)
        #XCoordsPerPlane.append(X)
        #YCoordsPerPlane.append(Y)
        #ZCoordsPerPlane.append(Z)
        

    return GridPtsPerPlane#, XCoordsPerPlane, YCoordsPerPlane, ZCoordsPerPlane



def GroupPointsIntoNearestPlanesTOOSLOW(Points, DicomDir):
    """ 
    Group a list of points into the nearest plane. 
    
    
    NOTE:
        This is taking too long.  It takes ~5 s to process one point, it will 
        take 21 min to process the 253 points in that one contour.  Scaling up 
        to multiple contours let along multiple contours with a higher number 
        of points is unrealistic.
    """
    
    import time
    
    # Start timing:
    times = []
    times.append(time.time())
    
    #GridPtsPerPlane, \
    #XCoordsPerPlane, \
    #YCoordsPerPlane, \
    #ZCoordsPerPlane = GetGridOfPointsInImagePlanes(DicomDir)
    GridPtsPerPlane = GetGridOfPointsInImagePlanes(DicomDir)
    
    S = len(GridPtsPerPlane)
    
    # Initialise the grouped points:
    GroupedPts = [ [] for i in range(S) ]
    
    for i in range(len(Points)):
        #pt_x, pt_y, pt_z = Points[i]
            
        # Initialise the list of plane distances:
        dPlanes = []
            
        # Loop through all planes in GridPtsPerPlane, find the indices of the 
        # grid point closest in [x, y], then calculate the distance:
        for j in range(len(GridPtsPerPlane)):
            GridPts = GridPtsPerPlane[j]
            
            # Initialise the distances between Points[i] and each :
            #dists_x = []
            #dists_y = []
            #dists_z = []
            
            ## Initialise the list of point distances:
            #dPoints = []
            
            ## Loop through each grid point:
            #for k in range(len(GridPts)):
            #    #GridPt = GridPts[k]
            #    
            #    dPoints.append(GetVectorLength(Points[i], GridPts[k]))
            
            dPoints = GetVectorLengths(Points[i], GridPts)
            
            # The index of the nearest point in GridPts to Points[i]:
            PtInd = dPoints.index(min(dPoints))
            
            #print(f'Points[{i}] = {Points[i]}')
            #print(f'GridPts[{PtInd}] = {GridPts[PtInd]}')
            
            dPlanes.append(GetVectorLength(Points[i], GridPts[PtInd]))
            
        # The index of the nearest plane in GridPtsPerPlane to Points[i]:
        PlInd = dPlanes.index(min(dPlanes))
        
        # Append Points[i] to GroupedPts
        GroupedPts[PlInd].extend(Points[i])
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        
        print(f'Took {Dtime} s to find the nearest plane index ({PlInd}) for',
              f'point {i} out of {len(Points)} points.\n\n')
        
    return GroupedPts




        
    
            
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
        X = []; Y = []; Z = []
            
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
            X = []; Y = []; Z = []
                
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
                + f'{LastSliceInd}__dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return
    
    
   
"""


"""


def PlotInterpolatedContours3DNew(InterpData, FixContourData, MovContourData, 
                                  Convert2ICS, dP, PlotImagingPlanes, 
                                  FixDicomDir, MovDicomDir, CombinePlots, 
                                  ExportPlot):
    """
    Plot interpolation results.
    
    
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    
    """ Aug 8 """
    if Convert2ICS:
        from PCStoICS import PCStoICS
        from GetImageAttributes import GetImageAttributes
    
        # Get the image attributes:
        FixOrigin, FixDirs, \
        FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                  Package='pydicom')
        
        MovOrigin, MovDirs, \
        MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                  Package='pydicom')
    
    
    FirstSliceInd = InterpData['BoundingSliceInds'][0][0]
    LastSliceInd = InterpData['BoundingSliceInds'][-1][1]
    
    
    NumOfPairs = len(InterpData['InterpSliceInd'])
        
    
    # Group all original and interpolated contours:
    FixContours = []
    MovContours = []
    
    """
    For each contour pair (Contour1 and Contour2) in InterpData, Contour2 of
    contour pair i is repeated as Contour1 for contour pair (i+1).  So only 
    need to gather Contour1 and InterpContour in each contour pair + Contour2 
    in the last pair.
    """
    
    # First gather all Contour1 and InterpContour by looping through each 
    # contour pair in InterpData.  The keys of the contours to gather are:
    keys = ['FixOSContour1Pts', 'FixInterpContourPts', 'FixOSContour2Pts',
            'MovOSContour1Pts', 'MovInterpContourPts', 'MovOSContour2Pts']
    
    
    for i in range(NumOfPairs):
        for key in keys:
            Contour = InterpData[key][i]
            
            """ Aug 7 """
            if Convert2ICS:
                if 'Fix' in key:
                    Origin=FixOrigin
                    Directions=FixDirs
                    Spacings=FixSpacings
                    
                if 'Mov' in key:
                    Origin=MovOrigin
                    Directions=MovDirs
                    Spacings=MovSpacings
                    
                # Convert from PCS to ICS:
                Contour = PCStoICS(Pts_PCS=Contour, Origin=Origin, 
                                   Directions=Directions, Spacings=Spacings)
                
            # Append Contour to FixContours or MovContours only append 
            # 'FixOSContour2Pts' or 'MovOSContour2Pts' if i = NumOfPairs - 1:
            if 'Fix' in key:
                if ( not '2' in key ) or ( ( '2' in key ) and ( i == NumOfPairs - 1 ) ):
                    FixContours.append(Contour)
                
            if 'Mov' in key:
                if ( not '2' in key ) or ( ( '2' in key ) and ( i == NumOfPairs - 1 ) ):
                    MovContours.append(Contour)
        
     
    # Now gather Contour2 in the last contour pair:
    #FixContour = InterpData['FixOSContour2Pts'][-1]
    #FixContours.append(FixContour)
    
    #MovContour = InterpData['MovOSContour2Pts'][-1]
    #MovContours.append(MovContour)
    
    # Group all contours:
    AllContours = [FixContours, MovContours]
    
    
    if PlotImagingPlanes:
        # Get meshgrids:
        FixMeshgrids = GetMeshGridForImagePlanes(ContourData=FixContourData, 
                                                 ContourTypeNo=1, 
                                                 DicomDir=FixDicomDir)
        
        MovMeshgrids = GetMeshGridForImagePlanes(ContourData=MovContourData, 
                                                 ContourTypeNo=3, 
                                                 DicomDir=MovDicomDir)

        
        # Group all meshgrids:
        AllMeshgrids = [FixMeshgrids, MovMeshgrids]                    
    
    
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
        fig = plt.figure(figsize=(14, 14))
        ax = fig.add_subplot(111, projection='3d')
    
        Contours = AllContours[i]
        
        # Set up the axes for this sub-plot:
        #if not CombinePlots:
        #    ax = fig.add_subplot(2, 1, i+1, projection='3d')
        
        # Store all x,y,z coordinates in arrays X, Y and Z:
        allX = []; allY = []; allZ = []
            
        # Loop through each contour:
        for j in range(len(Contours)):
            Contour = Contours[j]
            
            X = []; Y = []; Z = []
            
            # Unpack tuple and append to arrays X, Y and Z:
            for x, y, z in Contour:
                X.append(x)
                Y.append(y)
                Z.append(z)
                allX.append(x)
                allY.append(y)
                allZ.append(z)
                
            # Plot line:
            ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[i]); 
        
        # Define linestyle:
        #if (i % 2): # i.e. if i is even, it's an interpolated contour
        #    #LineStyle='dashed'
        #    LineStyle=(0, (5, 10)) # loosely dashed   
        #else: # if i is odd, it's not an interpolated contour
        #    LineStyle='solid'
            
            
        ##fig = plt.figure(figsize=plt.figaspect(0.5)*1.5, projection='3d') # Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
        ##fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
        #fig = plt.figure(figsize=(14, 14))
        #ax = fig.add_subplot(111, projection='3d')
        
        # Plot line:
        #ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[i]); 
        ##Axes3D.plot(X, Y, Z, linestyle=LineStyle, linewidth=1, c='gray');
        
        #ax.legend(loc='upper left', fontsize='large')
        
        #print(len(X), len(Y), len(Z), len(AveZ))
        
        # Plot dots:
        ##ax.scatter3D(X, Y, Z, markersize=MarkerSize, c=Z, cmap='hsv');
        #ax.scatter3D(X, Y, Z, c=Z, cmap='hsv');
        ax.scatter3D(allX, allY, allZ, s=MarkerSize, c=Colours[i], cmap='hsv');
        
        #print('FixMeshgrids shape =', FixMeshgrids.shape)
        
        if PlotImagingPlanes:  
            Meshgrids = AllMeshgrids[i]
            
            # Plot surface:
            for j in range(len(Meshgrids)):
            #for j in range(3):
                ax.plot_surface(Meshgrids[j][0], Meshgrids[j][1], Meshgrids[j][2], alpha=0.2)
                #ax.plot_trisurf(Meshgrids[j][0], Meshgrids[j][1], Meshgrids[j][2], 
                #                cmap=cm.jet, linewidth=0.1)
                #ax.plot_surface(Meshgrids[j][0], Meshgrids[j][1], Meshgrids[j][2], 
                #                rstride=1, cstride=1, cmap=cm.coolwarm,
                #                linewidth=0, antialiased=False)
            
            #plt.axis('equal') # doesn't work
            #ax.set_aspect('equal')
            #plt.axis('square')
            
            #axisEqual3D(ax)
            
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
                + f'{LastSliceInd}__dP_{dP}.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return






def PlotIntersectingPts2D(MovContourData, dP, UseInterp, AnnotatePtNums, 
                          SubPlots, ExportPlot):
    """
    Plot intersecting points.
    
    AnnotatePtNums - Boolean value determines whether points are annotated with
                     point numbers
                     
    SubPlots       - Boolean value determines whether a single figure with 
                     sub-plots or individual plots are generated.
                     
    
    Note:  
        I wasn't able to get retain a function that can either plot all results
        as subplots or as individual figures whilst maintaining equal, shared
        axes for all subplots.  So the function PlotIntersectingPts2D was split
        into two functions:  PlotIntersectingPts2D for plotting individual 
        figures, and PlotIntersectingPts2DSubPlots for plotting subplots.
        
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
            X = []; Y = []
        
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
                       + f'_Inters_pts_on_slice_{SliceNum}__dP_{dP}__' \
                       + f'UseInterp_{UseInterp}.png'
            
            if ExportPlot:
                plt.savefig(FigFname, bbox_inches='tight')
    
    FirstSlice = Inds[0]
    LastSlice = Inds[-1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                    + f'_Inters_pts_bt_slices_{FirstSlice}_and_{LastSlice}__' \
                    + f'dP_{dP}__UseInterp_{UseInterp}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return
    



def PlotIntersectingPts2DSubPlots(MovContourData, dP, UseInterp, 
                                  AnnotatePtNums, SubPlots, ExportPlot):
    """
    Plot intersecting points.
    
    AnnotatePtNums - Boolean value determines whether points are annotated with
                     point numbers
                     
    SubPlots       - Boolean value determines whether a single figure with 
                     sub-plots or individual plots are generated.
                     
    
    Note:  
        I wasn't able to get retain a function that can either plot all results
        as subplots or as individual figures whilst maintaining equal, shared
        axes for all subplots.  So the function PlotIntersectingPts2D was split
        into two functions:  PlotIntersectingPts2D for plotting individual 
        figures, and PlotIntersectingPts2DSubPlots for plotting subplots.
                     
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    colours = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    
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
        #MarkerSize = 2
        MarkerSize = 3
    else:
        MarkerSize = 5

    fig, axs = plt.subplots(Nrows, Ncols, sharex=True, sharey=True, figsize=(5*Ncols, 5*Nrows), dpi=DPI)
    
    # Initialise sub-plot iterator:
    i = 0
    
    for r in range(Nrows):
        for c in range(Ncols):
            if False:
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
                #ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c='b');
                axs[r, c].plot(X, Y, linestyle=LineStyle, linewidth=1, c=colours[j]);    
                # Plot dots:
                if True:
                    #plt.plot(X, Y, '.', markersize=MarkerSize, c='k');
                    axs[r, c].plot(X, Y, '.', markersize=MarkerSize, c=colours[j]);
                
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
                        #plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c='k');
                        axs[r, c].text(X[p], Y[p], p, fontsize=MarkerSize+5, c='k');
                # Also annotate the last point:
                #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
            
            #plt.axis('tight') 
            #plt.axis('equal')
            #axs[r, c].axis('equal')
            axs[r, c].set_aspect('equal', adjustable='box')
            #axs[r, c].set_aspect(1.0/axs[r, c].get_data_ratio(), adjustable='box')
        
            #plt.title(f'Intersecting points on Moving slice {SliceNum}')
            axs[r, c].set_title(f'Intersecting points on Moving slice {SliceNum}')
            
            if not SubPlots:
                # Create filename for exported figure:
                FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                           + f'_Inters_pts_on_slice_{SliceNum}__dP_{dP}__' \
                           + f'UseInterp_{UseInterp}.png'
                
                if ExportPlot:
                    plt.savefig(FigFname, bbox_inches='tight')
                    
            # Increment sub-plot iterator:
            i += 1 
    
    FirstSlice = Inds[0]
    LastSlice = Inds[-1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                   + f'_Inters_pts_bt_slices_{FirstSlice}_and_{LastSlice}__' \
                   + f'dP_{dP}__UseInterp_{UseInterp}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return



def Plot2DSubPlotsOfMovContourData(MovContourData, dP, AnnotatePtNums, SubPlots, 
                                   ExportPlot):
    """
    Plot intersecting points.
    
    AnnotatePtNums - Boolean value determines whether points are annotated with
                     point numbers
                     
    SubPlots       - Boolean value determines whether a single figure with 
                     sub-plots or individual plots are generated.
                     
    
    Note:  
        I wasn't able to get retain a function that can either plot all results
        as subplots or as individual figures whilst maintaining equal, shared
        axes for all subplots.  So the function PlotIntersectingPts2D was split
        into two functions:  PlotIntersectingPts2D for plotting individual 
        figures, and PlotIntersectingPts2DSubPlots for plotting subplots.
                     
    
    """
    
    import matplotlib.pyplot as plt
    import time
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    Colours = ['b', 'g', 'r', 'c', 'm', 'k', 'y']
    
    #AnnotationPeriod = 5 # annotate every AnnotationPeriod^th point
    AnnotationPeriod = 10 # annotate every AnnotationPeriod^th point
    
    # Get the indices of the slices in ContourData that have the desired 
    # ContourType value:
    #Inds = GetIndsOfSliceNumsOfContourType(ContourData=MovContourData, 
    #                                       ContourTypeNo=3)
    
    Inds = []
    
    ContourType = MovContourData['ContourType']
    
    for s in range(len(ContourType)):
        # If ContourType[s] is not empty:
        if ContourType[s]:
            Inds.append(s)
    
    #print('Inds =', Inds, '\n\n')
    
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
        
    #print(f'N = {N}, Nrows = {Nrows}, Ncols = {Ncols}')
    
    # Prepare the figure:
    if ExportPlot:
        DPI = 300
    else:
        DPI = 100
        
    if SubPlots:
        MarkerSize = 2
    else:
        MarkerSize = 5

    fig, axs = plt.subplots(Nrows, Ncols, sharex=True, sharey=True, figsize=(5*Ncols, 5*Nrows), dpi=DPI)
    
    # Initialise sub-plot iterator:
    i = 0
    
    for r in range(Nrows):
        for c in range(Ncols):
            ## Set up the axes for this sub-plot:
            #if SubPlots:
            #    ax = fig.add_subplot(Nrows, Ncols, i+1)
            #else:
            #    fig = plt.figure(figsize=(10, 10), dpi=DPI)
            #    ax = fig.add_subplot(111)
            
            #print(f'r = {r}, c = {c}, i = {i}, Inds[{i}] = {Inds[i]}')
            
            SliceNum = Inds[i]
            
            Contours = MovContourData['PointPCS'][SliceNum]
            
            for j in range(len(Contours)):
                Contour = Contours[j]
                Colour = Colours[j]
                
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
            
                for x, y, z in Contour:
                    X.append(x)
                    Y.append(y)
                    
                # Plot line:
                #ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c='b');
                axs[r, c].plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colour);    
                # Plot dots:
                if True:
                    #plt.plot(X, Y, '.', markersize=MarkerSize, c='k');
                    axs[r, c].plot(X, Y, '.', markersize=MarkerSize, c=Colour);
                
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
                        #plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c='k');
                        axs[r, c].text(X[p], Y[p], p, fontsize=MarkerSize+5, c='k');
                # Also annotate the last point:
                #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
            
            #plt.axis('tight') 
            #plt.axis('equal')
            #axs[r, c].axis('equal')
            axs[r, c].set_aspect('equal', adjustable='box')
            #axs[r, c].set_aspect(1.0/axs[r, c].get_data_ratio(), adjustable='box')
        
            #plt.title(f'Intersecting points on Moving slice {SliceNum}')
            axs[r, c].set_title(f'Intersecting points on Moving slice {SliceNum}')
            
            if not SubPlots:
                # Create filename for exported figure:
                FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                           + f'_Intersecting_points_on_slice_{SliceNum}__dP_{dP}.png'
                
                if ExportPlot:
                    plt.savefig(FigFname, bbox_inches='tight')
                    
            # Increment sub-plot iterator:
            i += 1 
            
            if i == N:
                break
    
    FirstSlice = Inds[0]
    LastSlice = Inds[-1]
    
    if SubPlots:
        # Create filename for exported figure:
        FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_Intersecting' \
                   + f'_points_between_slices_{FirstSlice}_and_{LastSlice}__dP_{dP}.png'
        
        if ExportPlot:
            plt.savefig(FigFname, bbox_inches='tight')
        
    return
    





def PlotPoints(Points, AnnotatePtNums, AnnotateEveryNthSegment, PlotTitle, 
               ExportPlot):
    """ Plot a list of points with option of annotating the point numbers """
    
    import matplotlib.pyplot as plt
    import time
    
    MarkerSize = 5
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    
    # Create a figure with two subplots and the specified size:
    fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    
    
    # Unpack tuple and store each x,y tuple in arrays X and Y:
    X = []
    Y = []

    for x, y, z in Points:
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

        for p in range(0, len(X) - 1, AnnotateEveryNthSegment):
            plt.text(X[p], Y[p], p, fontsize=MarkerSize+5);
        
    # Shift the position of the annotation for the last point so it doesn't overlap the
    # annotation for the first point:
    plt.text(1.005*X[-1], 1.005*Y[-1], len(X) - 1, fontsize=MarkerSize+5);
        
    plt.axis('equal')
    
    #plt.title(f'Every point in intersecting points joined to every other' \
    #          + f'point for slice no {SliceNum}')
    plt.title(PlotTitle)

    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        print('exporting plot...')
        plt.savefig(FigFname, bbox_inches='tight')
    
    return




def PlotAllContours(ContoursPerSlice, AnnotatePtNums, PlotLineSegments, 
                    EveryNth, PlotTitle, ExportPlot):
    """ Plot a list of contours for every slice with option of annotating the 
    point numbers """
    
    import matplotlib.pyplot as plt
    import time
    
    MarkerSize = 5
    
    #LineStyle='dashed'
    #LineStyle=(0, (5, 10)) # loosely dashed  
    LineStyle='solid'
    
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']
    
    # Create a figure with two subplots and the specified size:
    fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    
    
    # Number of slices:
    S = len(ContoursPerSlice)
    
    # Store all contours (skipping slices that have no contours):
    Contours = []
    PtsPerContour = []
    
    i = 0
    
    for s in range(S):
        ContoursThisSlice = ContoursPerSlice[s]
        
        if ContoursThisSlice:
            C = len(ContoursThisSlice)
            
            for c in range(C):
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
            
                for x, y, z in ContoursThisSlice[c]:
                    X.append(x)
                    Y.append(y)
            
                            
                # Plot this contour:
                ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[i]);    
                #ax.legend(loc='upper left', fontsize='large')
                
                # Plot dots:
                plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[i]);
                
                # Annotate with point numbers:
                if AnnotatePtNums:
                    P = list(range(len(X))) # list of point numbers
                    # Annotate every point:
                    #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
                    #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
                    
                    #plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5);
            
                    for p in range(0, len(X) - 1, EveryNth,):
                        plt.text(X[p], Y[p], p, fontsize=MarkerSize+5,
                                 c=Colours[i]);
                    
                    # Shift the position of the annotation for the last point so it doesn't overlap the
                    # annotation for the first point:
                    plt.text(1.005*X[-1], 1.005*Y[-1], len(X) - 1, 
                             fontsize=MarkerSize+5, c=Colours[i]);
                             
                Contours.append(ContoursThisSlice[c])
                
                PtsPerContour.append(len(ContoursThisSlice[c]))
                     
                i += 1
                
        
    # Plot line segments:
    if PlotLineSegments:
        # Continue only if all contours have the same number of points:
        if len(list(set(PtsPerContour))) == 1:
            # The number of points:
            P = PtsPerContour[0]
            
            # The number of contours:
            C = len(Contours)

            
            for p in range(0, P, EveryNth):
                # Initialise the list of points that sweep down the contours
                # for the p^th point in all contours:
                SweepPts = []
                
                for c in range(C):
                    SweepPts.append(Contours[c][p])
                
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
            
                for x, y, z in SweepPts:
                    X.append(x)
                    Y.append(y)
            
                # Plot this curve:
                ax.plot(X, Y, linestyle='dashed', linewidth=1, c='gray');
                
            
        else:
            print('\nCannot plot the contour sweep because the contours do',
                  'not have an equal number of points:', PtsPerContour)
        
        
    plt.axis('equal')
    
    #plt.title(f'Every point in intersecting points joined to every other' \
    #          + f'point for slice no {SliceNum}')
    plt.title(PlotTitle)

    # Create filename for exported figure:
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        print('exporting plot...')
        plt.savefig(FigFname, bbox_inches='tight')
    
    return





def PlotContoursAndPlanes3D(Contours, Convert2ICS, PlotImagingPlanes, 
                            PlotJoiningSegments, EveryNthSegment,
                            SegmentColoursRandom, DicomDir, PlotTitle, 
                            ExportPlot):
    """
    Plot contours (list of a list of points) and imaging planes (if desired)
    from a list of contours.
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    
    if Convert2ICS:
        from PCStoICS import PCStoICS
        from GetImageAttributes import GetImageAttributes
    
        # Get the image attributes:
        Origin, Directions, \
        Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                            Package='pydicom')
    
    #print(f'len(Meshgrids) = {len(Meshgrids)}')
    
    C = len(Contours)
    
    #print(f'len(Contours) = {C}')
    
    # Create a list of strings for colours, repeating 6 colours N times so that 
    # there are sufficient colours for each contour:
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']
    n = len(Colours)
    
    # Take the ceiling of 
    N = - ( - C // n )
    
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']*N
    
    # Initialise the number of points in each slice containing contours, and
    # the indices of the slices containing contours:
    NumPtsPerSliceContainingContours = []
    IndsSlicesContainingContours = []
    
    MarkerSize = 3
    
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection='3d')
    
        
    # Loop through all slices:    
    for s in range(C):
        contours = Contours[s]
        #meshgrid = Meshgrids[s]
        
        # Get the shape of the contour data for this slice:
        DataShape = np.array(contours).shape
        
        #LenDataShape = len(DataShape)
        
        Ncontours = DataShape[0]
        
        # If Ncontours != 0, get the Npoints in each contour:
        if Ncontours:
            # Initialise Npoints for each contour:
            #Npoints = []
            Npoints = 0
            
            # Loop through all contours:
            for c in range(Ncontours):
                #Npoints.append(len(contours[c]))
                Npoints += len(contours[c])
                
            NumPtsPerSliceContainingContours.append(Npoints)
            
            IndsSlicesContainingContours.append(s)
        else:
            Npoints = 0
            
        
        
        if False: #LogToConsole:
            print('\n\ntxt =', txt)
        
            print('\nind =', ind)
            
            print(f'\nDataShape = {DataShape}')
            #print(f'LenDataShape = {LenDataShape}')
            
            #print(f'\nNcts = {Ncts}')
            print(f'\nNcontours[ind={ind}] = {Ncontours}')
            
            #print(f'\npts[ind={ind}] =\n\n', pts[ind])
            print(f'\nNpoints[ind={ind}] = {Npoints}')
        
        
        
        
        # Continue if points exist for this slice:
        if contours:
            # Loop through each array of contour points:
            for c in range(len(contours)):
                contour = contours[c]
                
                if Convert2ICS:
                    # Convert from PCS to ICS:
                    contour = PCStoICS(Pts_PCS=contour, Origin=Origin, 
                                       Directions=Directions, Spacings=Spacings)
                    
    
                X = []; Y = []; Z = []
                
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in contour:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                    
                # Plot line:
                ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[s]);
                
                # Plot dots:
                ax.scatter3D(X, Y, Z, s=MarkerSize, c=Colours[s], cmap='hsv');
        
        #plt.title('Contours and imaging planes')
        #plt.title(PlotTitle)
        ax.set_title(PlotTitle)
        
        
    if PlotImagingPlanes:
        # Get meshgrids:
        #Meshgrids = GetMeshGridsForImagePlanesInContours(Contours=Contours, 
        #                                                 OnlyPlanesWithContours=OnlyPlanesWithContours, 
        #                                                 DicomDir=DicomDir)  
        ConvertMeshgrids2ICS = False
        Meshgrids = GetMeshGridsForImagePlanesNEW(DicomDir=DicomDir,
                                                  Convert2ICS=ConvertMeshgrids2ICS) 
        
        PlotAllImagePlanes = True
        PlotAllImagePlanes = False
        
        if PlotAllImagePlanes:
            # Loop through all meshgrids:
            for s in range(len(Meshgrids)):
                # Plot surface:
                ax.plot_surface(Meshgrids[s][0], Meshgrids[s][1], Meshgrids[s][2], 
                                #alpha=0.2)
                                alpha=0.1)
                                #alpha=1)
        
        else:
            # Get indices of slices that contain contours:
            inds = []
            
            for c in range(C):
                if Contours[c]:
                    inds.append(c)
            
            #for i in inds: # 14/09
            for i in range(inds[0] - 2, inds[-1] + 8): # 14/09
                # Plot surface:
                ax.plot_surface(Meshgrids[i][0], Meshgrids[i][1], Meshgrids[i][2], 
                                #alpha=0.2)
                                #alpha=0.1)
                                alpha=1)
    
    
    if PlotJoiningSegments:
        #print('NumPtsPerSliceContainingContours =', 
        #      NumPtsPerSliceContainingContours)
        
        print('IndsSlicesContainingContours =', IndsSlicesContainingContours)
        
        # Get the list of unique number of points per slice containing contours:
        SetOfNumPts = list(set(NumPtsPerSliceContainingContours))
                    
        print('SetOfNumPts =', SetOfNumPts) 
        
        # If len(SetOfNumPts) = 1, all contours have the same number of points,
        # so can be linked along each point in each contour:
        if len(SetOfNumPts) == 1:
            # The number of points:
            P = SetOfNumPts[0]
            
            print(f'P = {P}')
            
            # Loop through each point:
            for p in range(0, P, EveryNthSegment):
                # Initialise the list of points to be linked:
                LinkedPts = []
                
                # Loop through each contour:
                for i in IndsSlicesContainingContours:
                    LinkedPts.append(Contours[i][0][p])
                    
                if Convert2ICS:
                    # Convert from PCS to ICS:
                    LinkedPts = PCStoICS(Pts_PCS=LinkedPts, Origin=Origin, 
                                         Directions=Directions, 
                                         Spacings=Spacings)
                    
                X = []; Y = []; Z = []
                
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in LinkedPts:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                 
                if SegmentColoursRandom:
                    # Generate a random colour for this contour:
                    colour = [item/255 for item in list(np.random.choice(range(256), size=3))]
                else:
                    colour = 'gray'
        
                # Plot line:
                ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=colour);
                
        
        
    # Create filename for exported figure:
    #FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
    #            + f'_Contours_and_image_planes.png'
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + '_' + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return






def PlotContoursAndPlanes3DNEW(ContourData, PlotInterpContours, 
                               PlotImagingPlanes, Convert2ICS,
                               DicomDir, PlotTitle, ExportPlot):
    """
    Plot contours (list of a list of points) and imaging planes (if desired)
    from a dictionary.
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    
    #print(f'len(Meshgrids) = {len(Meshgrids)}')
    
    if PlotInterpContours:
        #C = len(Contours)
        C = len(ContourData['PointPCS'])
        
        Inds = list(range(C))
    else:
        # Get Inds of slices that have only original contours:
        Inds = GetIndsOfSliceNumsOfContourType(ContourData=ContourData, 
                                               ContourTypeNo=1)
        
        C = len(Inds)
    
    #print(f'len(Contours) = {C}')
    
    # Create a list of strings for colours, repeating 6 colours N times so that 
    # there are sufficient colours for each contour:
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']
    n = len(Colours)
    
    # Take the ceiling of 
    N = - ( - C // n )
    
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']*N
    
    MarkerSize = 3
    
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection='3d')
    
        
    # Loop through all slices:    
    for s in range(C):
        i = Inds[s]

        #contours = Contours[i]
        contours = ContourData['PointPCS'][i] # <-- dontours parallel to planes but not within bounds
        #contours = ContourData['PointICS'][i] # <-- contours not parallel to planes
        #meshgrid = Meshgrids[s]
        
        # Get the shape of the contour data for this slice:
        DataShape = np.array(contours).shape
        
        #LenDataShape = len(DataShape)
        
        Ncontours = DataShape[0]
        
        # If Ncontours != 0, get the Npoints in each contour:
        if Ncontours:
            # Initialise Npoints for each contour:
            Npoints = []
            
            # Loop through all contours:
            for c in range(Ncontours):
                Npoints.append(len(contours[c]))
        else:
            Npoints = 0
        
        if False: #LogToConsole:
            print('\n\ntxt =', txt)
        
            print('\nind =', ind)
            
            print(f'\nDataShape = {DataShape}')
            #print(f'LenDataShape = {LenDataShape}')
            
            #print(f'\nNcts = {Ncts}')
            print(f'\nNcontours[ind={ind}] = {Ncontours}')
            
            #print(f'\npts[ind={ind}] =\n\n', pts[ind])
            print(f'\nNpoints[ind={ind}] = {Npoints}')
        
        
        
        
        # Continue if points exist for this slice:
        if contours:
            # Loop through each array of contour points:
            for c in range(len(contours)):
                contour = contours[c]
                
                #if Convert2ICS:
                #    # Convert from PCS to ICS:
                #    contour = PCStoICS(Pts_PCS=contour, Origin=Origin, 
                #                       Directions=Directions, Spacings=Spacings)
                    
    
                X = []; Y = []; Z = []
                
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in contour:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                    
                # Plot line:
                ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[s]);
                
                # Plot dots:
                ax.scatter3D(X, Y, Z, s=MarkerSize, c=Colours[s], cmap='hsv');
        
        #plt.title('Contours and imaging planes')
        #plt.title(PlotTitle)
        ax.set_title(PlotTitle)
        
        
    if PlotImagingPlanes:
        # Get meshgrids:
        #Meshgrids = GetMeshGridsForImagePlanesInContours(Contours=Contours, 
        #                                                 OnlyPlanesWithContours=OnlyPlanesWithContours, 
        #                                                 DicomDir=DicomDir)   
        Meshgrids = GetMeshGridsForImagePlanesNEW(DicomDir=DicomDir, 
                                                  Convert2ICS=Convert2ICS) 
        
        PlotAllImagePlanes = True
        PlotAllImagePlanes = False
        
        if PlotAllImagePlanes:
            # Loop through all meshgrids:
            for s in range(len(Meshgrids)):
                # Plot surface:
                ax.plot_surface(Meshgrids[s][0], Meshgrids[s][1], Meshgrids[s][2], 
                                alpha=0.2)
        
        else:
            # Get Inds of slices that have only original contours:
            Inds = GetIndsOfSliceNumsOfContourType(ContourData=ContourData, 
                                               ContourTypeNo=1)
        
            for i in Inds:
                #print(f'i = {i}')
                
                # Plot surface:
                ax.plot_surface(Meshgrids[i][0], Meshgrids[i][1], Meshgrids[i][2], 
                                alpha=0.2)

    
    # Create filename for exported figure:
    #FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
    #            + f'_Contours_and_image_planes.png'
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + '_' + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return





def PlotContoursAndPlanes3D_NEW_2(InterpData, ContourData, FixOrMov,
                                  PlotImagingPlanes, PlotAllImagePlanes,
                                  PlotLineSegments, PlotEveryNthSegment,
                                  ConvertMeshgrids2ICS, DicomDir, 
                                  PlotTitle, ExportPlot):
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    #from PCStoICS import PCStoICS
    from GetImageAttributes import GetImageAttributes
    
    # Get the image attributes:
    Origin, Directions, \
    Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                        Package='pydicom')
    
    
    #Colours = ['b', 'r', 'g', 'm', 'y', 'k']
    #n = len(Colours)
    
    # Take the ceiling of 
    #N = - ( - C // n )
    
    #Colours = ['b', 'r', 'g', 'm', 'y', 'k']*N
    
    MarkerSize = 3
    
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection='3d')
    
    
    # Loop through each row in InterpData and plot the contours:
    for i in range(len(InterpData['BoundingSliceInds'])):
        
        Contour1 = InterpData[FixOrMov + 'OSContour1Pts'][i]
        Contour2 = InterpData[FixOrMov + 'OSContour2Pts'][i]
        
        # Unpack the arrays of x-, y- and z- components for Contour1 and 
        # Contour2: 
        X1 = []; Y1 = []; Z1 = []
                
        # Unpack tuple and append to arrays X, Y and Z:
        for x, y, z in Contour1:
            X1.append(x)
            Y1.append(y)
            Z1.append(z)
            
        X2 = []; Y2 = []; Z2 = []
                
        # Unpack tuple and append to arrays X, Y and Z:
        for x, y, z in Contour2:
            X2.append(x)
            Y2.append(y)
            Z2.append(z)
            
        # Plot contours:
        ax.plot3D(X1, Y1, Z1, linestyle='solid', linewidth=1, c='b');
        ax.plot3D(X2, Y2, Z2, linestyle='solid', linewidth=1, c='b');
        
        # Plot points:
        ax.scatter3D(X1, Y1, Z1, s=MarkerSize, c='b', cmap='hsv');
        ax.scatter3D(X2, Y2, Z2, s=MarkerSize, c='b', cmap='hsv');
        
        
        if PlotLineSegments:
            # Plot the line segments that connect the points in Contour1 and 
            # Contour2:
            for p in range(0, len(X1), PlotEveryNthSegment):
                ax.plot3D([X1[p], X2[p]], 
                          [Y1[p], Y2[p]],
                          [Z1[p], Z2[p]],
                          linestyle='dashed', linewidth=1, c='gray');
        
    
    # Plot the imaging planes:
    if PlotImagingPlanes:
        # Get meshgrids:  
        Meshgrids = GetMeshGridsForImagePlanesNEW(DicomDir=DicomDir,
                                                  Convert2ICS=ConvertMeshgrids2ICS) 
        
        #PlotAllImagePlanes = True
        #PlotAllImagePlanes = False
        
        if PlotAllImagePlanes:
            meshgrids = Meshgrids
        else:
            # Get Inds of slices that have only original contours:
            if FixOrMov == 'Fix':
                ContourTypeNo = 1
            else:
                ContourTypeNo = 3
                
            Inds = GetIndsOfSliceNumsOfContourType(ContourData=ContourData, 
                                                   ContourTypeNo=ContourTypeNo)
            
            #print('FixOrMov =', FixOrMov, ', ContourTypeNo =', ContourTypeNo,
            #      ', Inds =', Inds)
            
            #print('\nContourData["ContourType"] =', ContourData['ContourType'])
        
            meshgrids = []
            
            for i in Inds:
                meshgrids.append(Meshgrids[i])
                
                
        # Loop through all meshgrids:
        for s in range(len(meshgrids)):
            # Plot surface:
            ax.plot_surface(meshgrids[s][0], meshgrids[s][1], meshgrids[s][2], 
                            alpha=0.2)
                
           
    
    # Plot the intersection points of the transformed points on the Moving 
    # image planes if FixOrMov = 'Mov':
    if FixOrMov == 'Mov':
        Contours = ContourData['PointPCS']
        #Contours = ContourData['PointICS']
        
        colours = ['r', 'g', 'm', 'c', 'y', 'k']
        
        # Loop through each slice:
        for s in range(len(Contours)):
            # Contours for this slice:
            contours = Contours[s]
            
            # If contours != []:
            if contours:
                # Loop through each array of contour points:
                for c in range(len(contours)):
                    # Initialise lists of x and y coordinates: 
                    XI = []
                    YI = []
                    ZI = []
                    
                    # Unpack tuple of x,y,z coordinates and store in
                    # lists XI, YI and ZI:
                    for x, y, z in contours[c]:
                        XI.append(x)
                        YI.append(y)
                        ZI.append(z)
                        
                        # Plot the contour points as lines and points:
                        ax.plot3D(XI, YI, ZI, linestyle='solid', linewidth=1, c=colours[c]);
                        ax.scatter3D(XI, YI, ZI, s=MarkerSize, c=colours[c], cmap='hsv');
        
    
    # Create filename for exported figure:
    #FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
    #            + f'_Contours_and_image_planes.png'
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + '_' + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
    
    
    
    
    return






def PlotContourSweepAndPlanes3D(ContourSweep, Contours, Convert2ICS, 
                                PlotImagingPlanes, EveryNthSegment,
                                DicomDir, PlotTitle, ExportPlot):
    """
    Plot contour sweep (list of list of points that link the same node across
    contours, a new list for each node) and imaging planes (if desired).
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    
    if Convert2ICS:
        from PCStoICS import PCStoICS
        from GetImageAttributes import GetImageAttributes
    
        # Get the image attributes:
        Origin, Directions, \
        Spacings, Dims = GetImageAttributes(DicomDir=DicomDir, 
                                            Package='pydicom')
    
    #print(f'len(Meshgrids) = {len(Meshgrids)}')
    
    N = len(ContourSweep)
    
    #print(f'N = {N}')
    
    
    
    MarkerSize = 3
    
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection='3d')
    
        
    # Loop through all nodes:    
    for n in range(0, N, EveryNthSegment):
        SweepThisNode = ContourSweep[n]
        #meshgrid = Meshgrids[s]
        
        if Convert2ICS:
            # Convert from PCS to ICS:
            SweepThisNode = PCStoICS(Pts_PCS=SweepThisNode, Origin=Origin, 
                                     Directions=Directions, Spacings=Spacings)
            

        X = []; Y = []; Z = []
        
        # Unpack tuple and append to arrays X, Y and Z:
        for x, y, z in SweepThisNode:
            X.append(x)
            Y.append(y)
            Z.append(z)
            
        # Generate a random colour for this curve:
        colour = [item/255 for item in list(np.random.choice(range(256), size=3))]
        
        # Plot contour sweep curve:
        ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=colour);
        
        # Plot dots:
        ax.scatter3D(X, Y, Z, s=MarkerSize, c='g', cmap='hsv');
        
        
    if PlotImagingPlanes:
        # Get meshgrids:
        #Meshgrids = GetMeshGridsForImagePlanesInContours(Contours=Contours, 
        #                                                 OnlyPlanesWithContours=OnlyPlanesWithContours, 
        #                                                 DicomDir=DicomDir)  
        ConvertMeshgrids2ICS = False
        Meshgrids = GetMeshGridsForImagePlanesNEW(DicomDir=DicomDir,
                                                  Convert2ICS=ConvertMeshgrids2ICS) 
        
        PlotAllImagePlanes = True
        PlotAllImagePlanes = False
        
        if PlotAllImagePlanes:
            # Loop through all meshgrids:
            for s in range(len(Meshgrids)):
                # Plot surface:
                ax.plot_surface(Meshgrids[s][0], Meshgrids[s][1], Meshgrids[s][2], 
                                alpha=0.2)
        
        else:
            C = len(Contours)
            
            # Get indices of slices that contain contours:
            inds = []
            
            for c in range(C):
                if Contours[c]:
                    inds.append(c)
            
            for i in inds:
                """ 10/09: Fudge correction to deal with displacement: """
                i += 12
                
                # Plot surface:
                ax.plot_surface(Meshgrids[i][0], Meshgrids[i][1], Meshgrids[i][2], 
                                alpha=0.2)
    
    #plt.title('Contours and imaging planes')
    #plt.title(PlotTitle)
    ax.set_title(PlotTitle)
        
    # Create filename for exported figure:
    #FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
    #            + f'_Contours_and_image_planes.png'
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + '_' + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return






def GetMeshGridsForImagePlanes_v2(Origin, Dirs, Spacings, Dims, Convert2ICS,
                                  LogToConsole):
    """
    17/09/2020
    Create a list of points that define the grid of pixels that
    define the extent of the imaging planes for all image planes.
    
    Note that the points are in the Patient Coordinate System but will be 
    converted to the Image Coordinate System if Convert2ICS=True.
    
    Adapted from GetMeshGridsForImagePlanesNEW to avoid need to use
    SimpleElastix.
    """
    
    import numpy as np
    #from scipy.interpolate import griddata
    
    
    Meshgrids = []
    
    for s in range(Dims[2]):
        if LogToConsole:
            print(f'\ns = {s}')
            print('Origin =', Origin)
            print('s*Spacings[2]*Dirs[6] =', s*Spacings[2]*Dirs[6])
            print('s*Spacings[2]*Dirs[7] =', Origin)
        
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
        
        if Convert2ICS:
            # Import function:
            from PCStoICS import PCStoICS
            
            """ This is not an elegant approach (and it's VERY INEFFICIENT) but 
            a quick means of getting a result. """
            for i in range(Dims[0]):
                for j in range(Dims[1]):
                    PtPCS = [XX[j,i], YY[j,i], ZZ[j,i]]
                    
                    PtICS = PCStoICS(Pts_PCS=PtPCS, Origin=Origin, 
                                     Directions=Dirs, Spacings=Spacings)
                    
                    # Over-write the [j,i]^th elements in XX, YY and ZZ:
                    XX[j,i] = PtICS[0]
                    YY[j,i] = PtICS[1]
                    ZZ[j,i] = PtICS[2]
        
        
        
        Meshgrids.append([XX, YY, ZZ])
        #Meshgrids.append([XX, YY, ZZ1 + ZZ2])
        #Meshgrids.append([X, Y, Z])
        
        
    return Meshgrids







def PlotContoursAndPlanes3D_v2(Contours, Origin, Directions, Spacings, 
                               Dimensions, Convert2ICS, PlotImagingPlanes, 
                               PlotJoiningSegments, EveryNthSegment,
                               SegmentColoursRandom, PlotTitle, LogToConsole,
                               ExportPlot):
    """
    Plot contours (list of a list of points) and imaging planes (if desired)
    from a list of contours.
    
    Adapted from PlotContoursAndPlanes3D but avoids use of SimpleElastix.
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #from mpl_toolkits import mplot3d
    #from mpl_toolkits import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import time
    from PCStoICS import PCStoICS

    
    #print(f'len(Meshgrids) = {len(Meshgrids)}')
    
    C = len(Contours)
    
    #print(f'len(Contours) = {C}')
    
    # Create a list of strings for colours, repeating 6 colours N times so that 
    # there are sufficient colours for each contour:
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']
    n = len(Colours)
    
    # Take the ceiling of 
    N = - ( - C // n )
    
    Colours = ['b', 'r', 'g', 'm', 'y', 'k']*N
    
    # Initialise the number of points in each slice containing contours, and
    # the indices of the slices containing contours:
    NumPtsPerSliceContainingContours = []
    IndsSlicesContainingContours = []
    
    MarkerSize = 3
    
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection='3d')
    
        
    # Loop through all slices:    
    for s in range(C):
        contours = Contours[s]
        #meshgrid = Meshgrids[s]
        
        # Get the shape of the contour data for this slice:
        DataShape = np.array(contours).shape
        
        #LenDataShape = len(DataShape)
        
        Ncontours = DataShape[0]
        
        # If Ncontours != 0, get the Npoints in each contour:
        if Ncontours:
            # Initialise Npoints for each contour:
            #Npoints = []
            Npoints = 0
            
            # Loop through all contours:
            for c in range(Ncontours):
                #Npoints.append(len(contours[c]))
                Npoints += len(contours[c])
                
            NumPtsPerSliceContainingContours.append(Npoints)
            
            IndsSlicesContainingContours.append(s)
        else:
            Npoints = 0
            
        
        
        if False: #LogToConsole:
            print('\n\ntxt =', txt)
        
            print('\nind =', ind)
            
            print(f'\nDataShape = {DataShape}')
            #print(f'LenDataShape = {LenDataShape}')
            
            #print(f'\nNcts = {Ncts}')
            print(f'\nNcontours[ind={ind}] = {Ncontours}')
            
            #print(f'\npts[ind={ind}] =\n\n', pts[ind])
            print(f'\nNpoints[ind={ind}] = {Npoints}')
        
        
        
        
        # Continue if points exist for this slice:
        if contours:
            # Loop through each array of contour points:
            for c in range(len(contours)):
                contour = contours[c]
                
                if Convert2ICS:
                    # Convert from PCS to ICS:
                    contour = PCStoICS(Pts_PCS=contour, Origin=Origin, 
                                       Directions=Directions, Spacings=Spacings)
                    
    
                X = []; Y = []; Z = []
                
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in contour:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                    
                # Plot line:
                ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=Colours[s]);
                
                # Plot dots:
                ax.scatter3D(X, Y, Z, s=MarkerSize, c=Colours[s], cmap='hsv');
        
        #plt.title('Contours and imaging planes')
        #plt.title(PlotTitle)
        ax.set_title(PlotTitle)
        
        
    if PlotImagingPlanes: 
        #ConvertMeshgrids2ICS = False
        
        Meshgrids = GetMeshGridsForImagePlanes_v2(Origin=Origin, 
                                                  Dirs=Directions, 
                                                  Spacings=Spacings, 
                                                  Dims=Dimensions,
                                                  Convert2ICS=Convert2ICS,
                                                  LogToConsole=LogToConsole) 
        
        
        #PlotAllImagePlanes = True
        PlotAllImagePlanes = False
        
        if PlotAllImagePlanes:
            # Loop through all meshgrids:
            for s in range(len(Meshgrids)):
                # Plot surface:
                ax.plot_surface(Meshgrids[s][0], Meshgrids[s][1], Meshgrids[s][2], 
                                #alpha=0.2)
                                alpha=0.1)
                                #alpha=1)
        
        else:
            # Get indices of slices that contain contours:
            inds = []
            
            for c in range(C):
                if Contours[c]:
                    inds.append(c)
            
            #for i in inds: # 14/09
            for i in range(inds[0] - 2, inds[-1] + 8): # 14/09
                # Plot surface:
                ax.plot_surface(Meshgrids[i][0], Meshgrids[i][1], Meshgrids[i][2], 
                                #alpha=0.2)
                                alpha=0.1)
                                #alpha=1)
    
    
    if PlotJoiningSegments:
        if LogToConsole:
            print('NumPtsPerSliceContainingContours =', 
                  NumPtsPerSliceContainingContours)
            
            print('IndsSlicesContainingContours =', 
                  IndsSlicesContainingContours)
        
        # Get the list of unique number of points per slice containing contours:
        SetOfNumPts = list(set(NumPtsPerSliceContainingContours))
             
        if LogToConsole:
            print('SetOfNumPts =', SetOfNumPts) 
        
        # If len(SetOfNumPts) = 1, all contours have the same number of points,
        # so can be linked along each point in each contour:
        if len(SetOfNumPts) == 1:
            # The number of points:
            P = SetOfNumPts[0]
            
            if LogToConsole:
                print(f'P = {P}')
            
            # Loop through each point:
            for p in range(0, P, EveryNthSegment):
                # Initialise the list of points to be linked:
                LinkedPts = []
                
                # Loop through each contour:
                for i in IndsSlicesContainingContours:
                    LinkedPts.append(Contours[i][0][p])
                    
                if Convert2ICS:
                    # Convert from PCS to ICS:
                    LinkedPts = PCStoICS(Pts_PCS=LinkedPts, Origin=Origin, 
                                         Directions=Directions, 
                                         Spacings=Spacings)
                    
                X = []; Y = []; Z = []
                
                # Unpack tuple and append to arrays X, Y and Z:
                for x, y, z in LinkedPts:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
                 
                if SegmentColoursRandom:
                    # Generate a random colour for this contour:
                    colour = [item/255 for item in list(np.random.choice(range(256), size=3))]
                else:
                    colour = 'gray'
        
                # Plot line:
                ax.plot3D(X, Y, Z, linestyle='solid', linewidth=1, c=colour);
                
        
        
    # Create filename for exported figure:
    #FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
    #            + f'_Contours_and_image_planes.png'
    FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
                + '_' + PlotTitle.replace(' ', '_') + '.png'
    
    if ExportPlot:
        plt.savefig(FigFname, bbox_inches='tight')
        
    return






