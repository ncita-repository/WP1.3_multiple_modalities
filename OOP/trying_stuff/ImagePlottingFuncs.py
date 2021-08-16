# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:22:45 2020

@author: ctorti
"""




# Import packages:
import matplotlib.pyplot as plt
import SimpleITK as sitk
#import time
import numpy as np
#import os
#import re
#from plotly.subplots import make_subplots 
#import plotly.express as px
#import plotly.graph_objects as go



def PrepDataForFixMovRegPlotWithContours(FixImage, MovImage, RegImage, 
                                         FixContours, MovContours,
                                         SlicesToPlot,
                                         #FixSlicesToPlot, MovSlicesToPlot,
                                         Perspective, LogToConsole):
                                         
    
    #import ContourInterpolatingFuncs as cif
    
    FixNpa = sitk.GetArrayFromImage(FixImage)
    MovNpa = sitk.GetArrayFromImage(MovImage)
    RegNpa = sitk.GetArrayFromImage(RegImage)
    
    FixShape = FixNpa.shape # e.g. (30, 256, 212) => (z, x, y)
    MovShape = MovNpa.shape
    #RegShape = RegNpa.shape
    
    # Get dimensions in [x, y, z] order:
    FixDims = [FixShape[i] for i in [1, 2, 0]]
    MovDims = [MovShape[i] for i in [1, 2, 0]]
    #RegDims = [RegShape[i] for i in [1, 2, 0]]
    
    if LogToConsole: 
        #print(f'\nFixed dims = {FixDims}\nMoving dims = {MovDims}\nRegistered dims = {RegDims}')
        print(f'Fixed image dims      = {FixDims}')
        print(f'Moving image dims     = {MovDims}')
        print(f'Registered image dims = {FixDims}')
    
    FixOrigin = FixImage.GetOrigin()
    MovOrigin = MovImage.GetOrigin()
    #RegOrigin = RegImage.GetOrigin()
    
    FixSpacings = FixImage.GetSpacing() # e.g. (0.8984375, 0.8984375, 5.0) => (x, y, z)
    MovSpacings = MovImage.GetSpacing()
    #RegSpacings = RegImage.GetSpacing()
    
    FixOrient = [FixImage.GetDirection()[i] for i in range(len(FixImage.GetDirection()))]
    MovOrient = [MovImage.GetDirection()[i] for i in range(len(MovImage.GetDirection()))]
    
    # Reverse the sign of the cross terms so that the vectors are defined in
    # the same way as Pydicom:
    for i in [1, 2, 3, 5, 6, 7]:
        FixOrient[i] = - FixOrient[i]
        MovOrient[i] = - MovOrient[i]
    
    # Number of columns to plot:
    Ncols = 3
    
    
    # Get the number of slices and aspect ratio (depends on the chosen 
    # perspective):
    if Perspective == 'axial':
        FixNumSlices = FixShape[0] # z
        MovNumSlices = MovShape[0]
        #RegNumSlices = RegShape[0]
        
        FixPixAR = FixSpacings[1]/FixSpacings[0] # y/x
        MovPixAR = MovSpacings[1]/MovSpacings[0]
        #RegAR = RegSpacings[1]/RegSpacings[0]
    
        
        if True: # 11/09
            if SlicesToPlot == -1: # plot all slices
                # Number of rows to plot:
                Nrows = max(FixNumSlices, MovNumSlices)
                
                # Create list of slice indices:
                FixInds = list(range(FixNumSlices))
                MovInds = list(range(MovNumSlices))
             
            else:
                FixInds = SlicesToPlot
                MovInds = SlicesToPlot
                
                Nrows = len(SlicesToPlot)
                
        # 11/09:
        if False:
            if FixSlicesToPlot == -1: # plot all slices
                # Create list of slice indices:
                FixInds = list(range(FixNumSlices))
            else:
                FixInds = FixSlicesToPlot
                
            if MovSlicesToPlot == -1: # plot all slices
                MovInds = list(range(MovNumSlices))
            else:
                MovInds = MovSlicesToPlot
                
            Nrows = max([len(FixInds), len(MovInds)])
            
            
        # Get the extent of the x and y axes in the form 
        # [left, right, bottom, top]:
        for FixInd in FixInds:
            FixXmin = FixOrigin[0] + FixInd*FixSpacings[2]*FixOrient[6]
            
            FixYmin = FixOrigin[1] + FixInd*FixSpacings[2]*FixOrient[7]
            
            FixXmax = FixXmin + FixDims[0]*FixSpacings[0]*FixOrient[0] \
                              + FixDims[1]*FixSpacings[1]*FixOrient[3]
                              
            FixYmax = FixYmin + FixDims[0]*FixSpacings[0]*FixOrient[1] \
                              + FixDims[1]*FixSpacings[1]*FixOrient[4]
        
            FixAxisExtent = [FixXmin, FixXmax, FixYmax, FixYmin]
            
            FixDx = (FixXmax - FixXmin) / FixDims[0]
            
            FixDy = (FixYmax - FixYmin) / FixDims[1]
            
            FixAxisAR = FixDx / FixDy
            
        for MovInd in MovInds:
            MovXmin = MovOrigin[0] + MovInd*MovSpacings[2]*MovOrient[6]
            
            MovYmin = MovOrigin[1] + MovInd*MovSpacings[2]*MovOrient[7]
            
            MovXmax = MovXmin + MovDims[0]*MovSpacings[0]*MovOrient[0] \
                              + MovDims[1]*MovSpacings[1]*MovOrient[3]
                              
            MovYmax = MovYmin + MovDims[0]*MovSpacings[0]*MovOrient[1] \
                              + MovDims[1]*MovSpacings[1]*MovOrient[4]
        
            MovAxisExtent = [MovXmin, MovXmax, MovYmax, MovYmin]
            
            MovDx = (MovXmax - MovXmin) / MovDims[0]
            
            MovDy = (MovYmax - MovYmin) / MovDims[1]
            
            MovAxisAR = MovDx / MovDy    
            
        
    elif Perspective == 'coronal':
        FixNumSlices = FixShape[1] # x
        MovNumSlices = MovShape[1]
        #RegNumSlices = RegShape[1]
        
        #FixAR = FixSpacings[1]/FixSpacings[2] # y/z
        #MovAR = MovSpacings[1]/MovSpacings[2]
        #RegAR = RegSpacings[1]/RegSpacings[2]
        FixPixAR = FixSpacings[2]/FixSpacings[1] # z/y
        MovPixAR = MovSpacings[2]/MovSpacings[1]
        #RegAR = RegSpacings[2]/RegSpacings[1]
        
        """ It's not practical to plot all slices along x because there 
        will be hundreds of slices, so create a list of indices with a
        more manageable number. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nslices = min(FixShape[0], MovShape[0])
        
        FixInds = np.linspace(0, FixNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        MovInds = np.linspace(0, MovNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        Nrows = len(FixInds)
          
        if True: # 11/09
            if SlicesToPlot != -1: # NOT plot all slices
                # Plot the indices in plot_slices normalised by the total number
                # of rows and as a proportion of the total number of slices:
                FixInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(SlicesToPlot))]
                MovInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(SlicesToPlot))]
                
                Nrows = len(FixInds)
            
        # 11/09:
        if False:
            if FixSlicesToPlot != -1: # NOT plot all slices
                # Plot the indices in plot_slices normalised by the total number
                # of rows and as a proportion of the total number of slices:
                FixInds = [int( ( (FixSlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(FixSlicesToPlot))]
            if MovSlicesToPlot != -1:
                MovInds = [int( ( (MovSlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(MovSlicesToPlot))]
                
            Nrows = max([len(FixInds), len(MovInds)])
        
        
    elif Perspective == 'sagittal':
        FixNumSlices = FixShape[2] # y
        MovNumSlices = MovShape[2]
        #RegNumSlices = RegShape[2]
        
        FixPixAR = FixSpacings[2]/FixSpacings[0] # z/x
        MovPixAR = MovSpacings[2]/MovSpacings[0]
        #RegAR = RegSpacings[2]/RegSpacings[0]
        
        """ Same comment as above applies for coronal perspective. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nslices = min(FixShape[0], MovShape[0])
  
        FixInds = np.linspace(0, FixNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        MovInds = np.linspace(0, MovNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        Nrows = len(FixInds)
            
        if True: # 11/09
            if SlicesToPlot != -1: # NOT plot all slices
                # Plot the indices in plot_slices normalised by the total number
                # of rows and as a proportion of the total number of slices:
                FixInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(SlicesToPlot))]
                MovInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(SlicesToPlot))]
                
                Nrows = len(MovInds)
            
        # 11/09:
        if False:
            if FixSlicesToPlot != -1: # NOT plot all slices
                # Plot the indices in plot_slices normalised by the total number
                # of rows and as a proportion of the total number of slices:
                FixInds = [int( ( (FixSlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(FixSlicesToPlot))]
            if MovSlicesToPlot != -1:
                MovInds = [int( ( (MovSlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(MovSlicesToPlot))]
                
            Nrows = max([len(FixInds), len(MovInds)])
        
    else:
        print('Input "perspective" must be "axial", "sagital" or "coronal".')
        
        return
    
    if LogToConsole:
        print('\nFixInds =', FixInds)
        print('\nMovIinds =', MovInds)
        print('\nFixPixAR =', FixAR)
        print('\nMovPixAR =', MovAR)
    
    # Combine the data for fixed, moving and registered images in a list:
    Dims = [FixDims, MovDims, FixDims]
    Inds = [FixInds, MovInds, FixInds]
    Npas = [FixNpa, MovNpa, RegNpa]
    Origins = [FixOrigin, MovOrigin, FixOrigin]
    Spacings = [FixSpacings, MovSpacings, FixSpacings]
    Orients = [FixOrient, MovOrient, FixOrient]
    PixARs = [FixPixAR, MovPixAR, FixPixAR]
    NumSlices = [FixNumSlices, MovNumSlices, FixNumSlices]
    Txts = ['Fixed slice', 'Moving slice', 'Registered slice']
    Contours = [FixContours, MovContours, FixContours]
    AxisExtents = [FixAxisExtent, MovAxisExtent, FixAxisExtent]
    AxisARs = [FixAxisAR, MovAxisAR, FixAxisAR]
    
    return Dims, Inds, Npas, Origins, Spacings, Orients, PixARs, AxisExtents, \
           AxisARs, NumSlices, Txts, Contours, Nrows, Ncols
           



def PlotFixMovRegImagesAndContours(FixImage, MovImage, RegImage, 
                                   FixContours, MovContours,
                                   SlicesToPlot,
                                   #FixSlicesToPlot, MovSlicesToPlot,
                                   Perspective, UsePCSforAxes,
                                   ContoursAs, LogToConsole, ExportFig, 
                                   ExportFname):
    
    # Import packages:
    #import matplotlib.pyplot as plt
    #import numpy as np
    #import time
    
    """ New 14/09:
        Convert axes from ICS to PCS
        
        If UsePCSforAxes = True, FixContours and MovContours must be in PCS!
    """
    
    
    # Prepare data required for plotting:
    Dims, Inds, Npas, Origins, Spacings, Orients, PixARs, AxisExtents, \
    AxisARs, NumSlices, Txts, Contours, Nrows, Ncols \
    = PrepDataForFixMovRegPlotWithContours(FixImage, MovImage, RegImage, 
                                           FixContours, MovContours,
                                           SlicesToPlot,
                                           #FixSlicesToPlot, 
                                           #MovSlicesToPlot,
                                           Perspective,
                                           LogToConsole)
    
    
    #return Npas, NumSlices, Inds

#if False:
    
    # Create a figure with two subplots and the specified size:
    if ExportFig:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 7*Nrows), dpi=300) # 14/09
    else:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 7*Nrows)) # 14/09
    
    colours = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    
    LineWidth = 1
    LineWidth = 0.5
    
    MarkerSize = 1
    
    n = 0 # sub-plot number
    
    ## Loop through each slice: 14/09
    #for s in range(Nrows):
    # Loop through each row: 14/09
    for r in range(Nrows):
        # Loop through each data type:
        for i in range(Ncols):
            dims = Dims[i]
            npa = Npas[i]
            origin = Origins[i]
            spacings = Spacings[i]
            orient = Orients[i]
            PixAR = PixARs[i]
            AxisExtent = AxisExtents[i]
            AxisAR = AxisARs[i]
            Nslices = NumSlices[i]
            txt = Txts[i]
            contours = Contours[i]
            inds = Inds[i]
            ind = inds[r] # 14/09
            
            
            if LogToConsole:
                print('\n\ntxt =', txt)
             
                print('\ninds =', inds)
                
                print('\nNslices =', Nslices)
                
                print('\nPixAR =', PixAR)
            
            n += 1 # increment sub-plot number
            
            # Plot the image:
            #if s <= len(Inds[i]) - 1: # 11/09
            #if s + 1 <= Nslices: # 11/09
            if ind + 1 <= Nslices: # 14/09
                #ind = Inds[i][s]
                #ind = inds[s] # 14/09
                
                if LogToConsole:
                        print('\nind =', ind)
                
                #plt.subplot(Nrows, Ncols, n) # 14/09
                ax = plt.subplot(Nrows, Ncols, n) # 14/09
                #plt.axis('off') # 14/09
                #ax = plt.subplot(Nrows, Ncols, n, aspect=AR) # 14/09
                    
                if Perspective=='axial':
                    arr = npa[ind,:,:]
                    
                    SliceAxis = 'z'
                    #SliceLoc = origin[2] + spacings[2]*s # 14/09
                    SliceLoc = origin[2] + spacings[2]*ind
                    
                    SlicePlane = '(X, Y)'
                    SliceOrient = orient[0:2]
                    
                elif Perspective=='sagittal':
                    arr = np.flipud(npa[:,:,ind]) # flipud to correct orientation
                    
                    SliceAxis = 'y'
                    #SliceLoc = origin[1] + spacings[1]*s # 14/09
                    SliceLoc = origin[1] + spacings[1]*ind # 14/09
                    
                    SlicePlane = '(X, Z)'
                    SliceOrient = [orient[0], orient[2]]
                    
                elif Perspective=='coronal':
                    arr = np.flipud(npa[:,ind,:]) # flipud to correct orientation
                    
                    SliceAxis = 'x'
                    #SliceLoc = origin[0] + spacings[0]*s # 14/09
                    SliceLoc = origin[0] + spacings[0]*ind # 14/09
                    
                    SlicePlane = '(Y, Z)'
                    SliceOrient = orient[1:3]
                    
                SliceLoc = round(SliceLoc, 1)
                
                SliceOrient = [round(SliceOrient[i], 3) for i in range(len(SliceOrient))]
                
                #xMin = origin[0] + ind*spacings[2]*orient[6]
                #yMin = origin[1] + ind*spacings[2]*orient[7]
                #
                #xMax = xMin + dims[0]*spacings[0]*orient[0] + dims[1]*spacings[1]*orient[3]
                #yMax = yMin + dims[0]*spacings[0]*orient[1] + dims[1]*spacings[1]*orient[4]
                #
                #AxisExtent = [xMin, xMax, yMax, yMin]
                #
                #dx = (xMax - xMin) / dims[0]
                #dy = (yMax - yMin) / dims[1]
                #ar = dx / dy
                
                #plt.imshow(arr, cmap=plt.cm.Greys_r, aspect=AR); # 14/09
                #plt.imshow(arr, cmap=plt.cm.Greys_r, extent=AxisExtent,
                #           aspect=AR); # 14/09
                #ax = plt.gca()
                #ax.set_aspect(AR) 
                #ax.imshow(arr, cmap=plt.cm.Greys_r, extent=AxisExtent,
                #           aspect=AR); # 14/09
                #ax.imshow(arr, cmap=plt.cm.Greys_r, extent=AxisExtent); # 14/09
                if UsePCSforAxes:
                    ax.imshow(arr, cmap=plt.cm.Greys_r, extent=AxisExtent,
                               aspect=PixAR); # 14/09
                               #aspect=AxisAR); # 14/09
                    #ax.set_aspect(ar)
                else:
                    ax.imshow(arr, cmap=plt.cm.Greys_r, aspect=AxisAR); # 14/09
                
                #plt.title(txt + f' {ind + 1}/{Nslices}' \
                ax.set_title(txt + f' {ind + 1}/{Nslices}' \
                          + '\n' + SliceAxis + f'0 = {SliceLoc} mm' \
                          + '\n' + SlicePlane + f' = {SliceOrient}')
                #plt.axis('off')
                
                
                # Plot contours (currently for axial view only) if they exist
                # (i.e. if contours[ind] != []):
                if Perspective == 'axial' and contours[ind]:
                    # Get the shape of the contour data for this slice:
                    DataShape = np.array(contours[ind]).shape
                    
                    #LenDataShape = len(DataShape)
                    
                    Ncontours = DataShape[0]
                    
                    # If Ncontours != 0, get the Npoints in each contour:
                    if Ncontours:
                        # Initialise Npoints for each contour:
                        Npoints = []
                        
                        for c in range(Ncontours):
                            Npoints.append(len(contours[ind][c]))
                    else:
                        Npoints = 0
                    
                    if LogToConsole:
                        #print('\n\ntxt =', txt)
                    
                        #print('\nind =', ind)
                        
                        print('\nAxisAR =', AxisAR)
                        
                        print(f'\nDataShape = {DataShape}')
                        #print(f'LenDataShape = {LenDataShape}')
                        
                        #print(f'\nNcts = {Ncts}')
                        print(f'\nNcontours[ind={ind}] = {Ncontours}')
                        
                        #print(f'\npts[ind={ind}] =\n\n', pts[ind])
                        print(f'\nNpoints[ind={ind}] = {Npoints}')
                    
                    # Plot contours if Ncontours != 0:
                    if Ncontours:
                        # Loop through each array of contour points:
                        for c in range(Ncontours):
                            # Generate a random colour for this contour:
                            #FixRegColour = [item/255 for item in list(np.random.choice(range(256), size=3))]
                            #MovColour = [item/255 for item in list(np.random.choice(range(256), size=3))]
                            
                            #FixRegColour = [1, 0, 0] # red
                            #MovColour = [1, 0, 0] # red
        
                            # Number of points in this contour:
                            Npts = len(contours[ind][c])
                            #Npts = Npoints[c] # this line or above is fine
                            
                            # Initialise lists of x and y coordinates: 
                            X = []
                            Y = []
                            
                            if LogToConsole:
                                print(f'\nNpts[ind={ind}][c={c}] = {Npts}')
                            
                                #print(f'\ncontours[ind={ind}][c={c}] =', contours[ind][c])
                            
                            # Unpack tuple of x,y,z coordinates and store in
                            # lists X and Y:
                            for x, y, z in contours[ind][c]:
                                X.append(x)
                                Y.append(y)
                
                            # Flip Y pixel coordinates?:
                            if False:
                                N_r, N_c = np.shape(npa[ind,:,:])
                    
                                Y = N_c - np.array(Y)
                            
                            if ContoursAs=='lines':
                                
                                #if i == 1: # i.e. moving
                                #    plt.plot(X, Y, linewidth=1, c=MovColour)
                                #    
                                #else: # i.e. fixed or registered
                                #    plt.plot(X, Y, linewidth=1, c=FixRegColour)
                                
                                #print(f'i = {i}, n = {n}')
                                
                                #plt.plot(X, Y, linewidth=LineWidth, c=colours[c])
                                #plt.plot(X, Y, '.', markersize=MarkerSize, c=colours[c])
                                #print('')
                                ax.plot(X, Y, linewidth=LineWidth, c=colours[c])
                                ax.plot(X, Y, '.', markersize=MarkerSize, c=colours[c])
                                
                                
                            elif ContoursAs=='dots':
                                
                                #if i == 1: # i.e. moving
                                #    plt.plot(X, Y, '.', markersize=1, c=MovColour)
                                #    
                                #else: # i.e. fixed or registered
                                #    plt.plot(X, Y, '.', markersize=1, c=FixRegColour)
                                    
                                plt.plot(X, Y, '.', markersize=MarkerSize, c=colours[c])
                                
                            else:
                                print('\nContoursAs argument must be either',
                                      '"lines" or "dots".')
                                
                                
                #ax = plt.gca()
                #ax.set_aspect(AR)
                #ax.set_aspect('equal')
                #plt.tight_layout()
            
            
    if ExportFig:
        plt.savefig(ExportFname, bbox_inches='tight')
        
        
    print(f'PixAR = {PixAR}')
    print(f'AxisAR = {AxisAR}')
        
    return








def PrepDataForImageAndContoursPlot(FixContours, MovContours, SlicesToPlot, 
                                    Perspective, LogToConsole):
                                         
    """ 
    Adapted from PrepDataForFixMovRegPlotWithContours; 
    Doesn't require SimpleElastix.
    
    Ensure that the files FixIm.npy, MovIm.npy, RegIm.npy, 
    FixOrigin.json, FixDirs.json, FixSpacings.json, FixDims.json, 
    MovOrigin.json, MovDirs.json, MovSpacings.json and MovDims.json 
    have been copied to the working directory.
    """
    
    from ImportNumpyArray import ImportNumpyArray as ImportNpa
    from ImportImageAttributesFromJson import ImportImageAttributesFromJson as ImportImAtt
    from ContourInterpolatingFuncs import GetIndsOfNonEmptyItems
    
    
    # Import the data:

    FixNpa = ImportNpa('FixIm')
    MovNpa = ImportNpa('MovIm')
    RegNpa = ImportNpa('RegIm')
    
    FixOrig, FixDirs, FixSpacings, FixDims = ImportImAtt(FilePrefix='Fix')
    MovOrig, MovDirs, MovSpacings, MovDims = ImportImAtt(FilePrefix='Mov')

    FixShape = FixNpa.shape # e.g. (30, 256, 212) => (z, x, y)
    MovShape = MovNpa.shape
    #RegShape = RegNpa.shape
    
    # Get dimensions in [x, y, z] order:
    FixDims = [FixShape[i] for i in [1, 2, 0]]
    MovDims = [MovShape[i] for i in [1, 2, 0]]
    #RegDims = [RegShape[i] for i in [1, 2, 0]]
    
    if LogToConsole: 
        #print(f'\nFixed dims = {FixDims}\nMoving dims = {MovDims}\nRegistered dims = {RegDims}')
        print(f'Fixed image dims      = {FixDims}')
        print(f'Moving image dims     = {MovDims}')
        print(f'Registered image dims = {FixDims}')
    
    
    # Number of columns to plot:
    Ncols = 3
    
    
    # Get the number of slices and aspect ratio (depends on the chosen 
    # perspective):
    if Perspective == 'axial':
        FixNumSlices = FixShape[0] # z
        MovNumSlices = MovShape[0]
        #RegNumSlices = RegShape[0]
        
        FixPixAR = FixSpacings[1]/FixSpacings[0] # y/x
        MovPixAR = MovSpacings[1]/MovSpacings[0]
        #RegAR = RegSpacings[1]/RegSpacings[0]
    
 
        if SlicesToPlot == -1: # plot all slices
            # Number of rows to plot:
            Nrows = max(FixNumSlices, MovNumSlices)
            
            # Create list of slice indices:
            FixInds = list(range(FixNumSlices))
            MovInds = list(range(MovNumSlices))
            
            Nrows = len(FixInds)
            
        elif SlicesToPlot == -2: # plot all slices containing contours
            # Get the indices of the slices that contain contours for the 
            # Fixed and Moving image stacks:
            FixSlicesWithContours = GetIndsOfNonEmptyItems(List=FixContours)
            
            MovSlicesWithContours = GetIndsOfNonEmptyItems(List=MovContours)
            
            # Plot all slices containing contours in either image domains.
            # Concatenate the lists:
            FixOrMovSlicesWithContours = FixSlicesWithContours + MovSlicesWithContours
            
            # Reduce to set of unique indices:
            FixInds = list(set(FixOrMovSlicesWithContours)) 
            MovInds = list(set(FixOrMovSlicesWithContours)) 
            
            Nrows = len(FixInds)
         
        else:
            FixInds = SlicesToPlot
            MovInds = SlicesToPlot
            
            Nrows = len(SlicesToPlot)
                
            
            
        # Get the extent of the x and y axes in the form 
        # [left, right, bottom, top]:
        for FixInd in FixInds:
            FixXmin = FixOrig[0] + FixInd*FixSpacings[2]*FixDirs[6]
            
            FixYmin = FixOrig[1] + FixInd*FixSpacings[2]*FixDirs[7]
            
            FixXmax = FixXmin + FixDims[0]*FixSpacings[0]*FixDirs[0] \
                              + FixDims[1]*FixSpacings[1]*FixDirs[3]
                              
            FixYmax = FixYmin + FixDims[0]*FixSpacings[0]*FixDirs[1] \
                              + FixDims[1]*FixSpacings[1]*FixDirs[4]
        
            FixAxisExtent = [FixXmin, FixXmax, FixYmax, FixYmin]
            
            FixDx = (FixXmax - FixXmin) / FixDims[0]
            
            FixDy = (FixYmax - FixYmin) / FixDims[1]
            
            FixAxisAR = FixDx / FixDy
            
        for MovInd in MovInds:
            MovXmin = MovOrig[0] + MovInd*MovSpacings[2]*MovDirs[6]
            
            MovYmin = MovOrig[1] + MovInd*MovSpacings[2]*MovDirs[7]
            
            MovXmax = MovXmin + MovDims[0]*MovSpacings[0]*MovDirs[0] \
                              + MovDims[1]*MovSpacings[1]*MovDirs[3]
                              
            MovYmax = MovYmin + MovDims[0]*MovSpacings[0]*MovDirs[1] \
                              + MovDims[1]*MovSpacings[1]*MovDirs[4]
        
            MovAxisExtent = [MovXmin, MovXmax, MovYmax, MovYmin]
            
            MovDx = (MovXmax - MovXmin) / MovDims[0]
            
            MovDy = (MovYmax - MovYmin) / MovDims[1]
            
            MovAxisAR = MovDx / MovDy    
            
        
    elif Perspective == 'coronal':
        FixNumSlices = FixShape[1] # x
        MovNumSlices = MovShape[1]
        #RegNumSlices = RegShape[1]
        
        #FixAR = FixSpacings[1]/FixSpacings[2] # y/z
        #MovAR = MovSpacings[1]/MovSpacings[2]
        #RegAR = RegSpacings[1]/RegSpacings[2]
        FixPixAR = FixSpacings[2]/FixSpacings[1] # z/y
        MovPixAR = MovSpacings[2]/MovSpacings[1]
        #RegAR = RegSpacings[2]/RegSpacings[1]
        
        """ It's not practical to plot all slices along x because there 
        will be hundreds of slices, so create a list of indices with a
        more manageable number. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nslices = min(FixShape[0], MovShape[0])
        
        FixInds = np.linspace(0, FixNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        MovInds = np.linspace(0, MovNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        Nrows = len(FixInds)
          
        # if NOT plot all slices or all slices containing contours:
        if (SlicesToPlot != -1) and (SlicesToPlot != -2): 
            # Plot the indices in plot_slices normalised by the total number
            # of rows and as a proportion of the total number of slices:
            FixInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(SlicesToPlot))]
            MovInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(SlicesToPlot))]
            
            Nrows = len(FixInds)
            
        
        
    elif Perspective == 'sagittal':
        FixNumSlices = FixShape[2] # y
        MovNumSlices = MovShape[2]
        #RegNumSlices = RegShape[2]
        
        FixPixAR = FixSpacings[2]/FixSpacings[0] # z/x
        MovPixAR = MovSpacings[2]/MovSpacings[0]
        #RegAR = RegSpacings[2]/RegSpacings[0]
        
        """ Same comment as above applies for coronal perspective. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nslices = min(FixShape[0], MovShape[0])
  
        FixInds = np.linspace(0, FixNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        MovInds = np.linspace(0, MovNumSlices-1, num=Nslices, endpoint=True,
                              retstep=False, dtype=int)
        
        Nrows = len(FixInds)
            
        # if NOT plot all slices or all slices containing contours:
        if (SlicesToPlot != -1) and (SlicesToPlot != -2): 
            # Plot the indices in plot_slices normalised by the total number
            # of rows and as a proportion of the total number of slices:
            FixInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*FixNumSlices ) for i in range(len(SlicesToPlot))]
            MovInds = [int( ( (SlicesToPlot[i] + 1) )/Nslices*MovNumSlices ) for i in range(len(SlicesToPlot))]
            
            Nrows = len(MovInds)
            
        
    else:
        print('Input "perspective" must be "axial", "sagital" or "coronal".')
        
        return
    
    if LogToConsole:
        print('\nFixInds =', FixInds)
        print('\nMovIinds =', MovInds)
        print('\nFixPixAR =', FixPixAR)
        print('\nMovPixAR =', MovPixAR)
    
    # Combine the data for fixed, moving and registered images in a list:
    Dims = [FixDims, MovDims, FixDims]
    Inds = [FixInds, MovInds, FixInds]
    Npas = [FixNpa, MovNpa, RegNpa]
    Origins = [FixOrig, MovOrig, FixOrig]
    Spacings = [FixSpacings, MovSpacings, FixSpacings]
    Dirs = [FixDirs, MovDirs, FixDirs]
    PixARs = [FixPixAR, MovPixAR, FixPixAR]
    NumSlices = [FixNumSlices, MovNumSlices, FixNumSlices]
    Txts = ['Fixed slice', 'Moving slice', 'Registered slice']
    Contours = [FixContours, MovContours, FixContours]
    AxisExtents = [FixAxisExtent, MovAxisExtent, FixAxisExtent]
    AxisARs = [FixAxisAR, MovAxisAR, FixAxisAR]
    
    return Dims, Inds, Npas, Origins, Spacings, Dirs, PixARs, AxisExtents, \
           AxisARs, NumSlices, Txts, Contours, Nrows, Ncols
           
           
           
           
           


def PlotImagesAndContours(FixContours, MovContours,
                          SlicesToPlot, Perspective, UsePCSforAxes,
                          ContoursAs, LogToConsole, ExportFig, ExportFname):
    
    # Import packages:
    #import matplotlib.pyplot as plt
    #import numpy as np
    #import time
    from PCStoICS import PCStoICS
    from PCStoICS import PCStoICSnested
    
    """ 
    Adapted from PlotFixMovRegImagesAndContours; Doesn't require SimpleElastix.
    """
    
    
    # Prepare data required for plotting:
    Dims, Inds, Npas, Origins, Spacings, Orients, PixARs, AxisExtents, \
    AxisARs, NumSlices, Txts, Contours, Nrows, Ncols \
    = PrepDataForImageAndContoursPlot(FixContours, MovContours, SlicesToPlot,
                                      Perspective, LogToConsole)
    
    
    
    # Create a figure with two subplots and the specified size:
    if ExportFig:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows), dpi=300) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows), dpi=300) # 14/09
    else:
        #plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        #fig, ax = plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows)) # 14/09
        fig, ax = plt.subplots(Nrows, Ncols, figsize=(15, 9*Nrows)) # 14/09
    
    colours = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
    
    LineWidth = 1
    LineWidth = 0.5
    
    MarkerSize = 1
    
    n = 0 # sub-plot number
    
    ## Loop through each slice: 14/09
    #for s in range(Nrows):
    # Loop through each row: 14/09
    for r in range(Nrows):
        # Loop through each data type:
        for i in range(Ncols):
            #dims = Dims[i]
            npa = Npas[i]
            origin = Origins[i]
            spacings = Spacings[i]
            orient = Orients[i]
            PixAR = PixARs[i]
            AxisExtent = AxisExtents[i]
            AxisAR = AxisARs[i]
            Nslices = NumSlices[i]
            txt = Txts[i]
            contours = Contours[i]
            inds = Inds[i]
            ind = inds[r] # 14/09
            
            
            if LogToConsole:
                print('\n\ntxt =', txt)
             
                print('\ninds =', inds)
                
                print('\nNslices =', Nslices)
                
                print('\nPixAR =', PixAR)
            
            n += 1 # increment sub-plot number
            
            # Plot the image:
            #if s <= len(Inds[i]) - 1: # 11/09
            #if s + 1 <= Nslices: # 11/09
            if ind + 1 <= Nslices: # 14/09
                #ind = Inds[i][s]
                #ind = inds[s] # 14/09
                
                if LogToConsole:
                        print('\nind =', ind)
                
                #plt.subplot(Nrows, Ncols, n) # 14/09
                #plt.axis('off') # 14/09
                ax = plt.subplot(Nrows, Ncols, n) # 14/09

                #ax = plt.subplot(Nrows, Ncols, n, aspect=AR) # 14/09
                    
                if Perspective=='axial':
                    arr = npa[ind,:,:]
                    
                    SliceAxis = 'z'
                    #SliceLoc = origin[2] + spacings[2]*s # 14/09
                    SliceLoc = origin[2] + spacings[2]*ind
                    
                    SlicePlane = '(X, Y)'
                    SliceOrient = orient[0:2]
                    
                elif Perspective=='sagittal':
                    arr = np.flipud(npa[:,:,ind]) # flipud to correct orientation
                    
                    SliceAxis = 'y'
                    #SliceLoc = origin[1] + spacings[1]*s # 14/09
                    SliceLoc = origin[1] + spacings[1]*ind # 14/09
                    
                    SlicePlane = '(X, Z)'
                    SliceOrient = [orient[0], orient[2]]
                    
                elif Perspective=='coronal':
                    arr = np.flipud(npa[:,ind,:]) # flipud to correct orientation
                    
                    SliceAxis = 'x'
                    #SliceLoc = origin[0] + spacings[0]*s # 14/09
                    SliceLoc = origin[0] + spacings[0]*ind # 14/09
                    
                    SlicePlane = '(Y, Z)'
                    SliceOrient = orient[1:3]
                    
                SliceLoc = round(SliceLoc, 1)
                
                SliceOrient = [round(SliceOrient[i], 3) for i in range(len(SliceOrient))]
                
                
                if UsePCSforAxes:
                    ax.imshow(arr, cmap=plt.cm.Greys_r, extent=AxisExtent,
                               #aspect=PixAR); # 14/09
                               aspect=AxisAR); # 14/09
                    #ax.set_aspect(ar)
                    ax.set_xlabel('Pixels'); ax.set_ylabel('Pixels')
                    
                else:
                    ax.imshow(arr, cmap=plt.cm.Greys_r, aspect=AxisAR); # 14/09
                    ax.set_xlabel('mm'); ax.set_ylabel('mm')
                
                ax.set_title(txt + f' {ind + 1}/{Nslices}' \
                          + '\n' + SliceAxis + f'0 = {SliceLoc} mm' \
                          + '\n' + SlicePlane + f' = {SliceOrient}')
                #plt.axis('off')
                
                
                # Plot contours (currently for axial view only) if they exist
                # (i.e. if contours[ind] != []):
                if Perspective == 'axial' and contours[ind]:
                    # Get the shape of the contour data for this slice:
                    DataShape = np.array(contours[ind]).shape
                    
                    #LenDataShape = len(DataShape)
                    
                    Ncontours = DataShape[0]
                    
                    # If Ncontours != 0, get the Npoints in each contour:
                    if Ncontours:
                        # Initialise Npoints for each contour:
                        Npoints = []
                        
                        for c in range(Ncontours):
                            Npoints.append(len(contours[ind][c]))
                    else:
                        Npoints = 0
                    
                    if LogToConsole:
                        #print('\n\ntxt =', txt)
                    
                        #print('\nind =', ind)
                        
                        print('\nAxisAR =', AxisAR)
                        
                        print(f'\nDataShape = {DataShape}')
                        #print(f'LenDataShape = {LenDataShape}')
                        
                        #print(f'\nNcts = {Ncts}')
                        print(f'\nNcontours[ind={ind}] = {Ncontours}')
                        
                        #print(f'\npts[ind={ind}] =\n\n', pts[ind])
                        print(f'\nNpoints[ind={ind}] = {Npoints}')
                    
                    # Plot contours if Ncontours != 0:
                    if Ncontours:
                        # Loop through each array of contour points:
                        for c in range(Ncontours):
                            # Number of points in this contour:
                            Npts = len(contours[ind][c])
                            #Npts = Npoints[c] # this line or above is fine
                            
                            # Initialise lists of x and y coordinates: 
                            X = []
                            Y = []
                            
                            if LogToConsole:
                                print(f'\nNpts[ind={ind}][c={c}] = {Npts}')
                            
                                #print(f'\ncontours[ind={ind}][c={c}] =', contours[ind][c])
                                
                                
                            if UsePCSforAxes:
                                # Convert the points from PCS to ICS:
                                contour = PCStoICS(Pts_PCS=contours[ind][c],
                                #contour = PCStoICSnested(Pts_PCS=contours[ind][c],
                                                   Origin=origin,
                                                   Directions=orient,
                                                   Spacings=spacings)
                                
                            else:
                                contour = contours[ind][c]
                                
                            
                            # Unpack tuple of x,y,z coordinates and store in
                            # lists X and Y:
                            #for x, y, z in contours[ind][c]:
                            for x, y, z in contour:
                                X.append(x)
                                Y.append(y)
                
                            # Flip Y pixel coordinates?:
                            if False:
                                N_r, N_c = np.shape(npa[ind,:,:])
                    
                                Y = N_c - np.array(Y)
                            
                            if ContoursAs=='lines':
                                ax.plot(X, Y, linewidth=LineWidth, c=colours[c])
                                ax.plot(X, Y, '.', markersize=MarkerSize, c=colours[c])
                                
                                
                            elif ContoursAs=='dots':
                                plt.plot(X, Y, '.', markersize=MarkerSize, c=colours[c])
                                
                            else:
                                print('\nContoursAs argument must be either',
                                      '"lines" or "dots".')
                                
                                
                #ax = plt.gca()
                #ax.set_aspect(AR)
                #ax.set_aspect('equal')
                #plt.tight_layout()
            
            
    if ExportFig:
        plt.savefig(ExportFname, bbox_inches='tight')
        
        
    #print(f'PixAR = {PixAR}')
    #print(f'AxisAR = {AxisAR}')
        
    return









def display_image_and_contours(image, 
                               points, 
                               plot_slices,
                               perspective,
                               contours_as,
                               export_fig,
                               export_fname):
    
    
    
    npa = sitk.GetArrayFromImage(image)
    
    shape = npa.shape # e.g. (30, 256, 212) => (z, x, y)
    
    # Get dimensions in [x, y, z] order:
    dims = [shape[i] for i in [1, 2, 0]]
    
    print(f'Image dims      = {dims}')

    origin = image.GetOrigin()
    
    spacing = image.GetSpacing() # e.g. (0.8984375, 0.8984375, 5.0) => (x, y, z)
    
    orient = [image.GetDirection()[i] for i in range(len(image.GetDirection()))]
    
    # Reverse the sign of the cross terms so that the vectors are defined in
    # the same way as Pydicom:
    for i in [1, 2, 3, 5, 6, 7]:
        orient[i] = - orient[i]
    
    # The number of contours in points:
    Ncontours = 0
    
    for s in range(len(points)):
        Ncontours = Ncontours + len(points[s])
    
    
    # Number of columns to plot:
    Ncols = 1
    
    
    # Get the number of slices and aspect ratio (depends on the chosen 
    # perspective):
    if perspective == 'axial':
        Nslices = shape[0] # z

        AR = spacing[1]/spacing[0] # y/x
    
        
        if plot_slices == -1: # plot all slices
            # Number of rows to plot:
            Nrows = Nslices
            
            # Create list of slice indices:
            inds = list(range(Nslices))

        
        else:
            inds = plot_slices
            
            Nrows = len(plot_slices)
        
    elif perspective == 'coronal':
        Nslices = shape[1] # x
        
        AR = spacing[2]/spacing[1] # z/y
        
        """ It's not practical to plot all slices along x because there 
        will be hundreds of slices, so create a list of indices with a
        more manageable number. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nrows = shape[0]
        
        inds = np.linspace(0, Nslices-1, num=Nrows, endpoint=True,
                           retstep=False, dtype=int)

            
        if plot_slices != -1: # NOT plot all slices
            # Plot the indices in plot_slices normalised by the total number
            # of rows and as a proportion of the total number of slices:
            inds = [int( ( (plot_slices[i] + 1) )/Nslices*Nrows ) for i in range(len(plot_slices))]
            
            Nrows = len(inds)
        
        
    elif perspective == 'sagittal':
        Nslices = shape[2] # y
        
        AR = spacing[2]/spacing[0] # z/x

        
        """ Same comment as above applies for coronal perspective. """
        
        # Arbitrarily divide by the same number of axial slices:
        Nrows = shape[0]
  
        inds = np.linspace(0, Nslices-1, num=Nrows, endpoint=True,
                           retstep=False, dtype=int)
            
        if plot_slices != -1: # NOT plot all slices
            # Plot the indices in plot_slices normalised by the total number
            # of rows and as a proportion of the total number of slices:
            inds = [int( ( (plot_slices[i] + 1) )/Nslices*Nrows ) for i in range(len(plot_slices))]
            
        
    else:
        print('Input "perspective" must be "axial", "sagital" or "coronal".')
        
        return
    
    #print('\ninds =', inds)

    
    
    i = 0 # sub-plot number
    
    ##print(f'\nfix_AR = {fix_AR}\nmov_AR = {mov_AR}\nreg_AR = {reg_AR}')
    #print(f'\nfix_inds (N = {len(fix_inds)}) =', fix_inds)
    #print(f'\nmov_inds (N = {len(mov_inds)}) =', mov_inds)

    # Create a figure with two subplots and the specified size:
    plt.subplots(Nrows, Ncols, figsize=(5*Ncols, 5*Nrows))
    #if use_matplotlib:
    #    plt.subplots(Nrows, Ncols, figsize=(14, 5*Nrows))
    #    
    #if use_plotly:
    #    fig = make_subplots(rows=Nrows, cols=Ncols)
    
    # Contour counters:
    #fix_c = 0 # same counter will be used for registered
    #mov_c = 0
    
    # Loop through each slice:
    for r in range(Nrows):
        # Loop through each data type:
        i = i + 1 # increment sub-plot number
        
        # Plot the image:
        if r <= len(inds) - 1:
            ind = inds[r]
            
            plt.subplot(Nrows, Ncols, i)
            #if use_matplotlib:
            #    plt.subplot(Nrows, Ncols, i)
                
            if perspective=='axial':
                arr = npa[ind,:,:]
                
                sliceAxis = 'z'
                sliceLoc = origin[2] + spacing[2]*r
                
                slicePlane = '(X, Y)'
                sliceOrient = orient[0:2]
                
            elif perspective=='sagittal':
                arr = np.flipud(npa[:,:,ind]) # flipud to correct orientation
                
                sliceAxis = 'y'
                sliceLoc = origin[1] + spacing[1]*r
                
                slicePlane = '(X, Z)'
                sliceOrient = [orient[0], orient[2]]
                
            elif perspective=='coronal':
                arr = np.flipud(npa[:,ind,:]) # flipud to correct orientation
                
                sliceAxis = 'x'
                sliceLoc = origin[0] + spacing[0]*r
                
                slicePlane = '(Y, Z)'
                sliceOrient = orient[1:3]
                
            sliceLoc = round(sliceLoc, 1)
            
            sliceOrient = [round(sliceOrient[i], 3) for i in range(len(sliceOrient))]
            
            plt.imshow(arr, cmap=plt.cm.Greys_r, aspect=AR);
            #plt.title(f'Fixed slice {ind + 1}/{S}')
            plt.title(f'Slice {ind + 1}/{Nslices}' \
                      + '\n' + sliceAxis + f'0 = {sliceLoc} mm' \
                      + '\n' + slicePlane + f' = {sliceOrient}')
            plt.axis('off')
            
            
            # Plot contours (currently for axial view only):
            if perspective == 'axial':
                # Get the number of contours for this slice:
                Ncts = len(points[ind])

                
                print(f'\nSlice (ind =) {ind} has (Ncts =) {Ncts} contours')
                
                
                #print(f'\npoints[ind={ind}] =\n\n', points[ind])
                
                # Plot contours if Ncts is not 0:
                if Ncts:
                    # Loop through each array of contour points:
                    for c in range(Ncts):
                        # Generate a random colour for this contour:
                        colour = [item/255 for item in list(np.random.choice(range(256), size=3))]
    
                        # Number of points in this contour:
                        Npts = len(points[ind][c])
                        
                        # Unpack tuple and store each x,y tuple in arrays 
                        # X and Y:
                        X = []
                        Y = []
                        
                        print(f'   Contour (c =) {c} has (Npts =) {Npts} points')
                        
                        #print(f'\npoints[ind={ind}][c={c}] =\n\n', points[ind][c])
                        
                        for x, y, z in points[ind][c]:
                            X.append(x)
                            Y.append(y)
            
                        # Flip Y pixel coordinates?:
                        if False:
                            N_r, N_c = np.shape(npa[ind,:,:])
                
                            Y = N_c - np.array(Y)
                        
                        if contours_as=='lines':
                            plt.plot(X, Y, linewidth=1, c='red');
                            #plt.plot(X, Y, linewidth=1, c=colour);
                            
                            
                        elif contours_as=='dots':
                            plt.plot(X, Y, '.r', markersize=1);
                            #plt.plot(X, Y, '.', markersize=1, c=colour);
                                
                        else:
                            print('\ncontours_as argument must be either',
                                  '"lines" or "dots".')
    
            
            
    if export_fig:
        #plt.savefig('Images_and_reg_result_sitk.png', bbox_inches='tight')
        plt.savefig(export_fname, bbox_inches='tight')
        
    return






