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



def PlotRegResults(FixedDicoms, MovingDicoms, RegisteredDicoms, \
                   FixedRois, MovingRois, RegisteredRois, DataDict):
    
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
    
    # Get the keys for the fixed, moving and registered DICOMs:
    FixedKey = DataDict['FixedKey']
    MovingKey = DataDict['MovingKey']
    RegisteredKey = DataDict['RegisteredKey']
    
    
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
    PlotFname = 'Registration results using ' + ShiftFeature + ' ' \
                + CurrentDateTime + '.jpg'
    
    # Define the directory to export the plot:
    PlotDir = DataDict[RegisteredKey]['RoiDir']
    
    # Create the filepath for the exported plot:
    PlotFpath = os.path.join(PlotDir, PlotFname)
            
    
    # Configure plot:
    
    # Set the title font:
    fontSize=11
    
    # Set the number of subplot rows and columns:
    cols = len(AllDicoms)
    #rows = math.ceil(len(FixedDicoms)/cols)
    #rows = np.ceil(len(FixedDicoms)/cols).astype('int')
    rows = len(FixedDicoms) # all DICOMs objects have the same length
    
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
    plt.subplots(nrows=rows, ncols=cols, figsize=(5*3,6*rows));
    
    i = 0 # for subplot pos
    
    
    # Iterate through each slice:
    for r in range(rows):
        #print(f'r = {r}')
        
        # Iterate through each of three DICOMs:
        for c in range(cols):
            #print(f'   c = {c}')
            #print(f'   This Dicom has {len(AllDicoms[c])} slices.')
            #print(f'\nPlotting {r}th row in {c}th col:' \
            #      + AllNames[c] + ' objects..')
                        
            # Get the DICOM, ROI, slice number, scan position and scan  
            # positionaldifference for this row and column:
            Dicom = AllDicoms[c][r]
            Roi = AllRois[c]
            SliceNo = AllSliceNos[c][r]
            ScanPos = AllScanPos[c][r]
            ScanPosDiff = FixedToMovingScanPosDiffs[r]
            """ Note:
                If c = 1 (MovingDicoms), the scan positional difference is 
                FixedToMovingScanPosDiffs. If c = 0 (FixedDicoms), the scan
                positional difference is -FixedToMovingScanPosDiffs. """
            if c == 0:
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
            
            plt.subplot(rows,cols,i, aspect='equal')
            
            
            if flip:
                plt.pcolormesh(np.flipud(PixArray));
            else:
                plt.pcolormesh(PixArray);
            
            plt.title(f'Slice number = {SliceNo} in\n' 
                      #+ SeriesDesc \
                      + '\n'.join(textwrap.wrap(SeriesDesc, CharLim)) \
                      + f'\n(z = {round(ScanPos,2)} mm, dz = {round(ScanPosDiff,2)} mm)', \
                      size=fontSize);
            plt.axis('off');
            
            
            """ April 7:  Need to check this """
            
            # Plot the contour(s):
            """ Only plot the contours over the Fixed DICOMs and their copies
            over the registered DICOMs (so only for c=0 or c=2). """
            if Ncontours > 0 and c in [0, 2]:   
                for c in range(Ncontours):
                    # Unpack tuple and store each x,y tuple in arrays X and Y:
                    X = []
                    Y = []
                    for x, y in ContourPoints[c]:
                        X.append(x)
                        Y.append(y)
            
                    # Flip Y pixel coordinates if flip=True:
                    if flip:
                        N_r, N_c = np.shape(PixArray)
            
                        Y = N_c - np.array(Y)
            
                    plt.plot(X, Y, linewidth=0.5, c='red');
            #else:
            #    print(f'\nNo countours available for slice no. {r}')
            
       
    # Export the plot:
    plt.savefig(PlotFpath, bbox_inches='tight')
    
    print('Plot exported to:\n\n', PlotFpath)
    
    
    # Update the dictionary:
    DataDict.update({'PlotFpath':PlotFpath})
    
    return DataDict
 
            