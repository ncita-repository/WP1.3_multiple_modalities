# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 20:41:36 2020

@author: ctorti
"""



def PlotDicomsAndRois_v2(DicomDir, RoiFpath, SlicesToPlot, LogToConsole, 
                         ExportPlot, ExportFname):
    
    
    
    # Import packages and functions:
    from GetDicoms import GetDicoms
    #import time
    #from GetImageAttributes import GetImageAttributes
    #import importlib
    #import GetInputPoints
    #importlib.reload(GetInputPoints)
    from GetInputPoints import GetInputPoints
    #import PCStoICS
    #importlib.reload(PCStoICS)
    #from PCStoICS import PCStoICS
    #import copy
    import numpy as np
    from ContourInterpolatingFuncs import GetIndsOfNonEmptyItems
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    
    # Get all DICOMs in DicomDir:
    DicomFpaths, Dicoms = GetDicoms(DirPath=DicomDir, SortMethod='slices', 
                                    Debug=False)
    
    # Get contour points from RT-STRUCT file:
    if RoiFpath:
        PtsPCS, PtsBySliceAndContourPCS,\
        PtsICS, PtsBySliceAndContourICS,\
        LUT, PointData, ContourData,\
        NotUsed = GetInputPoints(FixedDicomDir=DicomDir, 
                                 FixedRoiFpath=RoiFpath)
    
    
    # Prepare the indices to plot:
    if SlicesToPlot == -1:
        # All slices will be plotted:
        SliceInds = list(range(len(DicomFpaths)))
        
    elif SlicesToPlot == -2:
        # All slices containing contours will be plotted.  Get the indices of 
        # the slices that contain contours if RoiFpath is not None, if it is
        # None, plot all slices:
        if RoiFpath:
            SliceInds = GetIndsOfNonEmptyItems(List=PtsBySliceAndContourICS)
            
        else:
            SliceInds = list(range(len(DicomFpaths)))
        
    elif isinstance(SlicesToPlot, list):
        # SlicesToPlot is a list so use it:
        SliceInds = SlicesToPlot
        
    else:
        print('SlicesToPlot must either be a list of integers, or -1 or -2.')
        
        return
        
    
    Nslices = len(SliceInds)
    
    # Set the number of subplot rows and columns:
    #Nrows = np.int8(np.round(np.sqrt(Nslices)))
    #Ncols = np.int8(np.round(np.sqrt(Nslices)))
    Ncols = np.int8(np.ceil(np.sqrt(Nslices)))
    #Nrows = np.int8(np.ceil(np.sqrt(Nslices)))
    Nrows = np.int8(np.ceil(Nslices/Ncols))
    
    # Limit the number of rows to 5 otherwise it's difficult to make sense of
    # the images:
    if Ncols > 5:
        Ncols = 5
        Nrows = np.int8(np.ceil(Nslices/Ncols))
    
    
    # Initiate the figure:
    if ExportPlot:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 7*Nrows), dpi=300) # 14/09
    else:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(7*Ncols, 9*Nrows)) # 14/09
    
    #colours = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    
    #LineWidth = 1
    LineWidth = 0.5
    
    MarkerSize = 1
    
    ContoursAs = 'lines'
    
    n = 0 # sub-plot number
    
    ## Loop through each slice: 14/09
    for s in range(Nslices):
        if LogToConsole:
                print('\ns =', s)
                
        n += 1 # increment sub-plot number
        
        #plt.subplot(Nrows, Ncols, n) # 14/09
        #plt.axis('off') # 14/09
        ax = plt.subplot(Nrows, Ncols, n) # 14/09

        #ax = plt.subplot(Nrows, Ncols, n, aspect=AR) # 14/09
          
        # Plot the pixel array:
        ax.imshow(Dicoms[SliceInds[s]].pixel_array, cmap=plt.cm.Greys_r)
        #ax.set_aspect(ar)
        ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
        
        ax.set_title(f'Slice {SliceInds[s]}')
        #plt.axis('off')
        
        
        
        # Plot contours if they exist for this slice:
        if RoiFpath:
            # Get the shape of the contour data for this slice:
            DataShape = np.array(PtsBySliceAndContourICS[SliceInds[s]]).shape
            
            #LenDataShape = len(DataShape)
            
            Ncontours = DataShape[0]
            
            # If Ncontours != 0, get the Npoints in each contour:
            if Ncontours:
                # Initialise Npoints for each contour:
                Npoints = []
                
                for c in range(Ncontours):
                    Npoints.append(len(PtsBySliceAndContourICS[SliceInds[s]][c]))
            else:
                Npoints = 0
            
            if LogToConsole:
                print(f'\nDataShape = {DataShape}')
                #print(f'LenDataShape = {LenDataShape}')
                
                #print(f'\nNcts = {Ncts}')
                print(f'\nNcontours[s={s}] = {Ncontours}')
                
                #print(f'\npts[ind={ind}] =\n\n', pts[ind])
                print(f'\nNpoints[s={s}] = {Npoints}')
            
            # Plot contours if Ncontours != 0:
            if Ncontours:
                # Loop through each array of contour points:
                for c in range(Ncontours):
                    # Number of points in this contour:
                    Npts = len(PtsBySliceAndContourICS[SliceInds[s]][c])
                    #Npts = Npoints[c] # this line or above is fine
                    
                    # Initialise lists of x and y coordinates: 
                    X = []
                    Y = []
                    
                    if LogToConsole:
                        print(f'\nNpts[s={s}][c={c}] = {Npts}')
                    
                        #print(f'\ncontours[ind={ind}][c={c}] =', contours[ind][c])
                        
                    
                    # Unpack tuple of x,y,z coordinates and store in
                    # lists X and Y:
                    #for x, y, z in contours[ind][c]:
                    for x, y, z in PtsBySliceAndContourICS[SliceInds[s]][c]:
                        X.append(x)
                        Y.append(y)
        
                    # Flip Y pixel coordinates?:
                    if False:
                        N_r, N_c = np.shape(Dicoms[SliceInds[s]].pixel_array)
            
                        Y = N_c - np.array(Y)
                    
                    if ContoursAs=='lines':
                        ax.plot(X, Y, linewidth=LineWidth, c='r')
                        ax.plot(X, Y, '.', markersize=MarkerSize, c='r')
                        
                        
                    elif ContoursAs=='dots':
                        plt.plot(X, Y, '.', markersize=MarkerSize, c='r')
                        
                    else:
                        print('\nContoursAs argument must be either',
                              '"lines" or "dots".')
                            
                                
            #ax = plt.gca()
            #ax.set_aspect(AR)
            #ax.set_aspect('equal')
            #plt.tight_layout()
            
            
    if ExportPlot:
        plt.savefig(ExportFname, bbox_inches='tight')
        
        
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
        
    return
        
    
    
    
        