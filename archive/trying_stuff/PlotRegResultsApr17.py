# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 17:32:03 2020

@author: ctorti
"""


"""
Function:
    PlotRegResults()
    
Purpose:
    Plot registration results.  The plot will contain an R x C array of 
    sub-plots for all R DICOMs contained in the arrays of DICOMs FixedDicoms,
    MovingDicoms and RegisteredDicoms (C = 3), with contours overlaid from 
    their corresponding ROI objects FixedRois, MovingRois and RegisteredRois.

Input:
    FixedDicoms      - array of DICOM objects used as "fixed" DICOMs in the
                       registration
                        
    MovingDicoms     - array of DICOM objects used as "moving" DICOMs in the
                       registration
    
    RegisteredDicoms - array of registered DICOM objects
    
    FixedRois        - ROI object for FixedDicoms
    
    MovingRois       - ROI object for MovingDicoms
    
    RegisteredRois   - ROI object for RegisteredDicoms
    
    DataDict         - dictionary of variables passed from previous functions
    
    FixedKey         - the key in DataDict that contains info of the data used
                       as fixed input to the registration
    
    MovingKey        - the key in DataDict that contains info of the data used
                       as moving input to the registration
    
    RegisteredKey    - the key in DataDict that contains info of the registered
                       data
    

    
Returns:
    DataDict         - updated dictionary
    
"""



def PlotRegResultsApr17(FixedDicoms, MovingDicoms, RegisteredDicoms, \
                   FixedRois, MovingRois, RegisteredRois, DataDict, \
                   TransformParametersFpath, RowToPlot, xShift, yShift, \
                   CombineBy):
    
    # Import packages:
    #import pydicom
    import matplotlib.pyplot as plt
    import numpy as np
    import time, os
    import importlib
    import GetContourPoints
    importlib.reload(GetContourPoints)
    from GetContourPoints import GetContourPoints
    import textwrap
    from ParseTransformParameters import ParseTransformParameters
    import copy
    
    # Get the keys for the fixed, moving and registered DICOMs:
    FixedKey = DataDict['FixedKey']
    MovingKey = DataDict['MovingKey']
    RegisteredKey = DataDict['RegisteredKey']
    
    
    # Get the Image Position Patient for the Fixed, Moving DICOMs:
    """ Note: The IPP for the Registered DICOMs are the same as that of the 
    Moving DICOMs since the latter were used as a template for the former, 
    copying most DICOM tags unchanged, including IPP. """
    FixedImPos = FixedDicoms[0].ImagePositionPatient
    MovingImPos = MovingDicoms[0].ImagePositionPatient
    #RegisteredImPos = RegisteredDicoms[0].ImagePositionPatient
    
    # Get the Pixel Spacing for the Fixed, Moving DICOMs:
    FixedPixSpacing = FixedDicoms[0].PixelSpacing
    MovingPixSpacing = MovingDicoms[0].PixelSpacing
    
    # The ratio of the Pixel Spacing of the Moving to Fixed DICOMs:
    XPixSpacingRatio = MovingPixSpacing[0]/FixedPixSpacing[0]
    YPixSpacingRatio = MovingPixSpacing[1]/FixedPixSpacing[1]
    
    # Get the Xshift and Yshift between the origin of the Fixed and Moving IPP
    # and scale by the ratio of their Pixel Spacings:
    XShift = float((FixedImPos[0] - MovingImPos[0])*XPixSpacingRatio)
    YShift = float((FixedImPos[1] - MovingImPos[1])*YPixSpacingRatio)
    
    """ Note:  The two lines above overwrite the xShift and yShift that were
    manually provided as arguments to the function. To use those remember to
    comment out the above two lines. """
    
    # Combine all DICOMs into an array:
    AllDicoms = [FixedDicoms, MovingDicoms, RegisteredDicoms]
    
    # Combine all ROIs into an array:
    AllRois = [FixedRois, MovingRois, RegisteredRois]
    
    # Get the slice numbers:
    FixedSliceNos = DataDict[FixedKey]['SliceNos']
    MovingSliceNos = DataDict[MovingKey]['SliceNos']
    RegisteredSliceNos = DataDict[RegisteredKey]['SliceNos']
            
    # Combine all SliceNos into an array:
    AllSliceNos = [FixedSliceNos, MovingSliceNos, RegisteredSliceNos]
    
    # Get the scan positions:
    FixedScanPos = DataDict[FixedKey]['ScanPos']
    MovingScanPos = DataDict[MovingKey]['ScanPos']
    RegisteredScanPos = DataDict[RegisteredKey]['ScanPos']
            
    # Combine all scan positions into an array:
    AllScanPos = [FixedScanPos, MovingScanPos, RegisteredScanPos]
    
    # Get the scan position differences:
    FixedToMovingScanPosDiffs = DataDict['FixedToMovingScanPosDiffs']
    
    # Create array of string for plotting purposes:
    AllNames = ['Fixed', 'Moving', 'Registered']
    
    
    # Create the filepath for the exported plot:
    # Get the current date-time:
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%H", time.gmtime())
    
    # Get the feature that was used to manually determine the required shift
    # that was applied to the target DICOMs:
    ShiftFeature = DataDict['ShiftFeature']
    
    # Create the filename for the exported plot:
    #PlotFname = 'Registration results (' + CurrentDateTime + ').jpg'
    PlotFname = 'Reg results with fixed and registered images ' \
                'combined using ' + CombineBy \
                + f', XShift {round(XShift,2)}, YShift {round(YShift,2)} ' \
                + 'using ' + ShiftFeature + ' ' + CurrentDateTime + '.jpg'
    
    # Define the directory to export the plot:
    PlotDir = DataDict[RegisteredKey]['RoiDir']
    
    # Create the filepath for the exported plot:
    PlotFpath = os.path.join(PlotDir, PlotFname)
    
    # Get the origin from the TransformParameters text file:
    Origin = ParseTransformParameters(TransformParametersFpath, 'Origin')
            
    # Convert to float:
    """ April 20:  This is now done within ParseTransformParameters() """
    #Origin = [float(item) for item in Origin]
    
    #print(f'\nOrigin = {Origin}')
    #print(f'Origin[0] = {Origin[0]}')
    #print(f'type(Origin[0]) =', str(type(Origin[0])))
    #print(f'Origin[1] = {Origin[1]}')
    #print(f'type(Origin[1]) =', str(type(Origin[1])))


    # Configure plot:
    
    # Set the title font:
    fontSize=11
    
    # Set the number of subplot rows and columns:
    Ncols = len(AllDicoms)
    Ncols = len(AllDicoms) + 1 # +1 for additional col showing combined image
    cols = list(range(Ncols))
    
    
    #rows = math.ceil(len(FixedDicoms)/cols)
    #rows = np.ceil(len(FixedDicoms)/cols).astype('int')
    if RowToPlot=='all':
        Nrows = len(FixedDicoms) # all DICOMs objects have the same length
        rows = list(range(Nrows))
    else:
        Nrows = 1
        rows = [RowToPlot]   
    
    #print(f'\nrows = {rows}, cols = {cols}')
    #print(f'\nlen(FixedDicoms) = {len(FixedDicoms)}')
    #print(f'\nlen(MovingDicoms) = {len(MovingDicoms)}')
    #print(f'\nlen(RegisteredDicoms) = {len(RegisteredDicoms)}')
    
    # Decide whether to flip the images up-down so that the orientation is
    # the same as shown in the OHIF-Viewer:
    #flip = False
    flip = True
    
    # Define the character limit for text wrapping of the figure titles:
    CharLim = 40
    
    
    # Plot results:
    
    #plt.figure(figsize=(7,15), dpi=300);
    #plt.subplots(nrows=rows, ncols=cols, figsize=(5*3,5*rows));
    # Make the row heights a bit larger to account for wrapped text:
    plt.subplots(nrows=Nrows, ncols=Ncols, figsize=(5*3,6*Nrows), dpi=300);
    
    i = 0 # for subplot pos
    
    
    # Iterate through each slice:
    for row in rows:
        #print(f'r = {r}')
        
        # Iterate through each of three DICOMs:
        for col in cols[0:3]:
            #print(f'   c = {c}')
            print(f'\n' + AllNames[col] + ':\n-----------')
            #print(f'   This Dicom has {len(AllDicoms[c])} slices.')
            #print(f'\nPlotting {r}th row in {c}th col:' \
            #      + AllNames[c] + ' objects..')
            
            # Get the DICOM, slice number, scan position and scan positional 
            # difference for this row and column:
            Dicom = AllDicoms[col][row]
            Roi = AllRois[col]
            SliceNo = AllSliceNos[col][row]
            ScanPos = AllScanPos[col][row]
            ScanPosDiff = FixedToMovingScanPosDiffs[row]
            """ Note:
                If col = 1 (MovingDicoms), the scan positional difference is 
                FixedToMovingScanPosDiffs. If col = 0 (FixedDicoms), the scan
                positional difference is -FixedToMovingScanPosDiffs. """
            if col == 0:
                ScanPosDiff = - ScanPosDiff
            
            # Get the array of flattened Contour Points:
            ContourPoints = GetContourPoints(Dicom, Roi)
            
            # The number of contours:
            Ncontours = len(ContourPoints)
            
            # Get the pixel array:
            PixArray = Dicom.pixel_array
            
            # Get the Series Description:
            SeriesDesc = Dicom.SeriesDescription
            
            # Get the z-component of the Image Position Patient:
            #zPos = str(round(float(Dicom.ImagePositionPatient[2]), 2))
            
      
            # Increment sub-plot number:
            i = i + 1  
            
            #print(f'\ni = {i}')
            
            plt.subplot(Nrows,Ncols,i, aspect='equal')
            
            
            if flip:
                plt.pcolormesh(np.flipud(PixArray));
            else:
                plt.pcolormesh(PixArray);
            
            plt.title(f'Slice number = {SliceNos[row]} in\n' 
                      #+ SeriesDesc \
                      + '\n'.join(textwrap.wrap(SeriesDesc, CharLim)) \
                      + f'\n(z = {round(ScanPos,2)} mm, dz = {round(ScanPosDiff,2)} mm)', \
                      size=fontSize);
            plt.axis('off');
            
            
            """ April 7:  Need to check this """
            
            # Plot the contour(s) for the fixed and registered images only:
            print(f'\nThere are {Ncontours} contours.')
            if Ncontours > 0 and col in [0, 2]:   
                for c in range(Ncontours):
                    print(f'\nContour number {c}:')
                    # Unpack tuple and store each x,y tuple in arrays X and Y:
                    X = []
                    Y = []
                    print(f'\nThere are {len(ContourPoints[c])} points in this contour.')
                    for x, y in ContourPoints[c]:
                        X.append(x)
                        Y.append(y)
                        
                    
                    #print(f'\nOrigin before sign reversal = {Origin}')
                    #Origin = [-item for item in Origin]
                    #print(f'Origin after sign reversal = {Origin}')
                    #print(f'\n(X[0], Y[0]) before flip = ({X[0]}, {Y[0]})')
                    #print(f'Origin before flip = {Origin}')
                    
                    if False:
                        # Copy Origin:
                        OriginCopy = copy.deepcopy(Origin)
                        
                        # Reserve sign:
                        OriginCopy = [-item for item in OriginCopy]
                        
                        print(f'TransformParameter Origin before sign change = {Origin}')
                        print(f'TransformParameter Origin after sign change  = {OriginCopy}')
            
                    # Flip Y pixel coordinates if flip=True:
                    if flip:
                        Nx, Ny = np.shape(PixArray)
                        
                        if False:
                            #print('\nY =', Y)
                            #print('\nN_c - Y =', [N_c - y for y in Y])
                            #print(f'N_c = {N_c}')
                            print(f'Pixel array dims = {np.shape(PixelArray)}')
                            print(f'First contour point before up-down flip       = ({X[0]}, {Y[0]})')
                            print(f'First contour point after up-down flip        = ({X[0]}, {N_y - Y[0]})')
                            print(f'TransformParameter Origin before up-down flip = {Origin}')
                            print(f'TransformParameter Origin after up-down flip  = ({Origin[0]}, {N_y - Origin[1]})')
                    
                        Y = Ny - np.array(Y)
                        
                        if False:
                            OriginCopy[1] = Ny - OriginCopy[1]

            
                    """ April 17:
                        Shift (X,Y) by Origin from TransformParameters for the
                        registration plot:
                    """
                    if col==0: # Plot contours for fixed image
                        plt.plot(X, Y, linewidth=1, c='red');
                    elif col==2: # Plot contours for registered image
                        #print('\nContourPoints:\n', ContourPoints)
                        #print('\nX:\n', X)
                        #print('\nY:\n', Y)
                        #print('\nOrigin =', Origin)
                        
                        plt.plot(X, Y, linewidth=1, c='red');
                        
                        # Subtract Origin from X,Y:
                        if False:#True:
                            Xdiff = [x - OriginCopy[0] for x in X]
                            Ydiff = [y - OriginCopy[1] for y in Y]
                            
                            #print(f'\n(X[0] - OriginCopy[0], Y[0] - OriginCopy[1]) = ({X[0] - OriginCopy[0]}, {Y[0] - OriginCopy[1]})')
                            #print(f'Xdiff[0] = {Xdiff[0]}')
                            #print(f'Ydiff[0] = {Ydiff[0]}')
                        if False:
                            Xdiff = [OriginCopy[0] - x for x in X]
                            Ydiff = [OriginCopy[1] - y for y in Y]
                            
                            #print(f'\n(OriginCopy[0] - X[0], OriginCopy[1] - Y[0]) = ({OriginCopy[0] - X[0]}, {OriginCopy[1] - Y[0]})')
                            #print(f'Xdiff[0] = {Xdiff[0]}')
                            #print(f'Ydiff[0] = {Ydiff[0]}')
                       
                        if False:
                            print('\ntype(X[0])                                    =', type(X[0]))
                            print(f'len(X)                                        = {len(X)}')
                            print('type(Xdiff[0])                                =', type(Xdiff[0]))
                            print(f'len(Xdiff)                                    = {len(Xdiff)}')
                            
                            
                            # Apply (Xdiff,Ydiff) to (X,Y):
                            #Xshifted = [x - Xdiff for x in X]
                            #Xshifted = np.array(X) - Xdiff
                            #Xshifted = X - Xdiff # <-- why doesn't this work?
                            XShifted = [X[i] - Xdiff[i] for i in range(len(X))]
                            #Yshifted = [y - Ydiff for y in Y]
                            #Yshifted = np.array(Y) - Ydiff
                            #Yshifted = Y - Ydiff # <-- why doesn't this work?
                            YShifted = [Y[i] - Ydiff[i] for i in range(len(Y))]
                            
                            print(f'\n(X[0], Y[0])                                  = ({X[0]}, {Y[0]})')
                            print(f'(Xdiff[0], Ydiff[0])                          = ({Xdiff[0]}, {Ydiff[0]})')
                            print(f'(XShifted[0], YShifted[0])                    = ({XShifted[0]}, {YShifted[0]})')
                            print(f'len(XShifted)                                 = {len(XShifted)}')
                            
                            #print(f'\n(X, Y)               = {list(zip(X,Y))}')
                            #print(f'\n(Xdiff, Ydiff)       = {list(zip(Xdiff,Ydiff))}')
                            #print(f'\n(Xshifted, Yshifted) = {list(zip(Xshifted,Yshifted))}')
                        
                            plt.plot(Xshifted, Yshifted, linewidth=1, c='red');
                            #plt.plot(0, 0, 'wo');
                            
                        
                        # Apply XShift and YShift to X and Y:
                        XShifted = [X[i] + XShift for i in range(len(X))]
                        YShifted = [Y[i] + YShift for i in range(len(Y))]
                        
                        plt.plot(XShifted, YShifted, linewidth=1, c='white');
                        
            #else:
            #    print(f'\nNo countours available for slice no. {r}')
            
            
        # Plot the 4th column (combined checkerboard or vertical split):
        # Increment sub-plot number:
        i = i + 1  
        
        plt.subplot(Nrows,Ncols,i, aspect='equal')
        
        # Get the pixel array of the fixed and registered DICOMs:
        FixedPixArray = FixedDicoms[RowToPlot].pixel_array
        RegisteredPixArray = RegisteredDicoms[RowToPlot].pixel_array
        
        # Get the max pixel values:
        MaxFixed = np.max(FixedPixArray)
        MaxRegistered = np.max(RegisteredPixArray)
        
        print(f'\nmax of Fixed      = {MaxFixed}')
        print(f'max of Registered = {MaxRegistered}')
    
        
        # Equate the contrast in the two pixel arrays:
        if MaxFixed > MaxRegistered:
            scale = int(MaxFixed/MaxRegistered)
            #print(f'scale             = {scale}')
            #print(f'int(scale)        = {int(scale)}')
            #print(f'round(scale)      = {round(scale)}')
            #print(f'int(round(scale)) = {int(round(scale))}')
            
            #RegisteredPixArray = RegisteredPixArray*int(round(scale))
            RegisteredPixArray = RegisteredPixArray*scale
        else:
            scale = int(MaxRegistered/MaxFixed)
            #print(f'scale             = {scale}')
            #print(f'int(scale)        = {int(scale)}')
            #print(f'round(scale)      = {round(scale)}')
            #print(f'int(round(scale)) = {int(round(scale))}')
            
            #FixedPixArray = FixedPixArray*int(round(scale))
            FixedPixArray = FixedPixArray*scale
        
        print(f'\nmax of Fixed      = {MaxFixed}')
        print(f'max of Registered = {MaxRegistered}')
        
        if CombineBy=='checkerboard':
            # Create the checkboard mask:
            #Mask = np.kron([[1, 0] * 4, [0, 1] * 4] * 4, np.ones((32, 32)))
            Mask = np.kron([[1, 0] * 8, [0, 1] * 8] * 8, np.ones((16, 16)))
            #Mask = np.kron([[1, 0] * 16, [0, 1] * 16] * 16, np.ones((8, 8)))
        
        elif CombineBy=='split':
            # Start by creating a mask of ones:
            Mask = np.ones((Nx,Ny))
            
            # Set columns Nx/2 onwards = 0:
            Mask[:,int(Nx/2):] = 0
            
    
        
        Mask = np.array(Mask, dtype=bool)
        
        # Create the inverse mask:
        InvMask = np.logical_not(Mask)
        
        # Combine FixedPixArray and RegisteredPixArray:
        CombinedPixArray = FixedPixArray*Mask + RegisteredPixArray*InvMask
        
        
        if flip:
            plt.pcolormesh(np.flipud(CombinedPixArray));
        else:
            plt.pcolormesh(CombinedPixArray);
        
        plt.title('Combined fixed and registered \nimages with copied contours' \
                  + '\nin red and shifted by' \
                  + f'\n({round(XShift,2)}, {round(YShift,2)}) mm in white', \
                  size=fontSize);
        plt.axis('off');
            
        
        # Plot contours for the combined image:
        plt.plot(X, Y, linewidth=1, c='red');
                        
        plt.plot(XShifted, YShifted, linewidth=1, c='white');
            
            
       
    # Export the plot:
    if True:#False:
        plt.savefig(PlotFpath, bbox_inches='tight')
        
        print('Plot exported to:\n\n', PlotFpath)
    
    
    # Update the dictionary:
    DataDict.update({'PlotFpath':PlotFpath})
    
    #return DataDict
    #return X, Y, Xdiff, Ydiff, Xshifted, Yshifted, Origin, OriginCopy
    return Origin, X, Y, XShifted, XShifted
 
            