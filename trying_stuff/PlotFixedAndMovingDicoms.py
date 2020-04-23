# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:21:53 2020

@author: ctorti
"""




"""
Function:
    PlotFixedAndMovingDicoms()
    
Purpose:
    Plot the DICOMs that will be inputted to the function RegisterDicoms().  
    The plot will contain an R x C array of sub-plots for all R DICOMs 
    contained in the arrays of DICOMs FixedDicoms and MovingDicoms (C = 2).

Input:
    FixedDicoms      - array of DICOM objects used as "fixed" DICOMs in the
                       registration
                        
    MovingDicoms     - array of DICOM objects used as "moving" DICOMs in the
                       registration
    
    DataDict         - dictionary of variables passed from previous functions
    
    FixedKey         - the key in DataDict that contains info of the data used
                       as fixed input to the registration
    
    MovingKey        - the key in DataDict that contains info of the data used
                       as moving input to the registration
                       
    RegisteredKey    - the key in DataDict that contains info of the registered
                       data
    
    
Returns:
    N/A
    
"""



def PlotFixedAndMovingDicoms(FixedDicoms, MovingDicoms, DataDict):
    
    # Import packages:
    #import pydicom
    import matplotlib.pyplot as plt
    import numpy as np
    import time, os
    #import importlib
    #import GetContourPoints
    #importlib.reload(GetContourPoints)
    #from GetContourPoints import GetContourPoints
    import textwrap
    
    
    # Get the keys for the fixed, moving and registered DICOMs:
    FixedKey = DataDict['FixedKey']
    MovingKey = DataDict['MovingKey']
    RegisteredKey = DataDict['RegisteredKey']
    
    # Combine all DICOMs into an array:
    AllDicoms = [FixedDicoms, MovingDicoms]
    
    # Get the slice numbers:
    FixedSliceNos = DataDict[FixedKey]['SliceNos']
    MovingSliceNos = DataDict[MovingKey]['SliceNos']
            
    # Combine all SliceNos into an array:
    AllSliceNos = [FixedSliceNos, MovingSliceNos]
    
    # Get the scan positions:
    FixedScanPos = DataDict[FixedKey]['ScanPos']
    MovingScanPos = DataDict[MovingKey]['ScanPos']
    
    # Get the scan position differences:
    FixedToMovingScanPosDiffs = DataDict['FixedToMovingScanPosDiffs']
            
    # Combine all scan positions into an array:
    AllScanPos = [FixedScanPos, MovingScanPos]
    
    # Create array of string for plotting purposes:
    AllNames = ['Fixed', 'Moving']
    
    
    # Create the filepath for the exported plot:
    # Get the current date-time:
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%H", time.gmtime())
    
    # Get the feature that was used to manually determine the required shift
    # that was applied to the target DICOMs:
    ShiftFeature = DataDict['ShiftFeature']
    
    # Create the filename for the exported plot:
    PlotFname = 'Fixed and Moving DICOMs using ' + ShiftFeature + ' ' \
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
            
            # Get the DICOM, slice number, scan position and scan positional 
            # difference for this row and column:
            Dicom = AllDicoms[c][r]
            SliceNo = AllSliceNos[c][r]
            ScanPos = AllScanPos[c][r]
            ScanPosDiff = FixedToMovingScanPosDiffs[r]
            """ Note:
                If c = 1 (MovingDicoms), the scan positional difference is 
                FixedToMovingScanPosDiffs. If c = 0 (FixedDicoms), the scan
                positional difference is -FixedToMovingScanPosDiffs. """
            if c == 0:
                ScanPosDiff = - ScanPosDiff
            
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
            
            
       
    # Export the plot:
    plt.savefig(PlotFpath, bbox_inches='tight')
    
    print('Plot exported to:\n\n', PlotFpath)
    
    
    # Update the dictionary:
    #DataDict.update({'PlotFpath':PlotFpath})
    
    #return DataDict
    return