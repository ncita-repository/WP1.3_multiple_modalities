# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:41:21 2020

@author: ctorti
"""

""" Function:
    PlotDICOMsAndROIs()
    
Purpose:
    Plot DICOMs and overlay ROIs from DICOM-RTSRTUCT file

Input:
    dicomsPath - path to directory containing DICOM files
    
    roiPath    - path to directory containing DICOM-RTSTRUCT file
    
    exportPlot - (binary) If =1 export the plot
    
    exportDir  - path to directory where plot should be exported to
    
    
Returns:
    dataDict   - dictionary containing the DICOM pixel arrays and contour data
                from the DICOM-RTSTRUCT file; the keys of the dictionary are
                the SOP Instance UIDs of the DICOMs (and the Reference SOP
                Instance UIDs of the ROIs)
"""



import numpy as np
import matplotlib.pyplot as plt
import pydicom
from DicomHelperFuncs import GetDicomFpaths
import copy
import os

def PlotDICOMsAndROIs(dicomsPath, roiPath, exportPlot, exportDir):
    # First get the filepaths of the DICOMs:
    dicomsFpaths = GetDicomFpaths(dicomsPath)
    
    
    # Create a dictionary, dicomsDict to store the pixel array, image origin,
    # and pixel spacing, using the SOP Instance UIDs as keys:
    dicomsDict = {}

    for f in range(len(dicomsFpaths)):
        fpath = dicomsFpaths[f]
        dicom = pydicom.read_file(fpath)
        uid = dicom.SOPInstanceUID
        #position = float(dicom.ImagePositionPatient) # TypeError: float() 
        # argument must be a string or a number, not 'MultiValue'
        position = np.array(dicom.ImagePositionPatient).astype(np.float)
        orientation = np.array(dicom.ImageOrientationPatient).astype(np.float)
        pixSpacing = np.array(dicom.PixelSpacing).astype(np.float)
        pixelArray = dicom.pixel_array
    
        dicomsDict[uid] = {'Frame number':f+1, \
                           'Position':position, \
                           'Orientation':orientation, \
                           'Pixel spacing':pixSpacing, \
                           'Pixel array':pixelArray
                          }
    
    
    # Load the DICOM-RTSTRUCT file:
    roi = pydicom.dcmread(roiPath)
    
    
    # Create a dictionary, contoursDict, to store the contour data and number 
    # of contour points, using the Referenced SOP Instance UIDs as keys:
    
    # Note:  
    # Since the contour number, number of contour points and the contour data 
    # elements are all strings, convert the contour number and number of 
    # contour points to integers, and store the contour data array as a numpy 
    # array so it can easily be converted to float.
    
    contoursDict = {}
    
    # For any given UID...
    # count the number of contours:
    # Ncontours = 0 # initial value
    # create a list of the contour number(s):
    # contourNo = [] 
    # create a list of the number of contour points:
    # NcontourPts = []
    # create a list of contour data:
    # contourData = []
    
    # Loop for each contour sequence:
    for sequence in roi.ROIContourSequence[0].ContourSequence:
        uid = sequence.ContourImageSequence[0].ReferencedSOPInstanceUID
        
        # Get the keys (UIDs) in contoursDict so far:
        keys = list(contoursDict.keys())
        
        # Check if this uid is in keys. If it isn't this is the first time
        # for this uid.  Initialise the values:
        if uid not in keys:
            Ncontours = 0
            contourNo = [] 
            NcontourPts = []
            contourData = []
            
        # Append the values for this sequence:
        Ncontours = Ncontours + 1 # increment
        contourNo.append(int(sequence.ContourNumber))
        NcontourPts.append(int(sequence.NumberOfContourPoints))
        contourData.append(np.array(sequence.ContourData).astype(np.float))
         
        # Store the values for this uid. Values will be over-written if there's more than
        # one contour for this uid, but the lists have been appended already:
        contoursDict[uid] = {'Number of contours':Ncontours, \
                             'Contour number':contourNo, \
                             'Number of contour points':NcontourPts, \
                             'Contour data':contourData
                            }
        
        
    # Combine contoursDict and dicomsDict into a new dictionary called 
    # dataDict:
    Ckeys = list(contoursDict.keys())
    Dkeys = list(dicomsDict.keys())
    
    # Create a new dictionary called dataDict that combines both dictionaries. 
    # Start by copying dicomsDict:
    dataDict = copy.deepcopy(dicomsDict)
    
    for Dkey in Dkeys: # cycle through each key in dicomsDict
        # Check if Dkey is in Ckeys:
        if Dkey in Ckeys:
            # There is contour data for this Dkey, so update dataDict[Dkey]
            # with contoursDict[Ckey]:
            ind = Ckeys.index(Dkey) # index of matching key in contoursDict
            dataDict[Dkey].update(contoursDict[Dkey]) # Note: Dkey=Ckey
        
        else:
            # There is no contour data for this Dkey, so fill in the 
            # missing values:
            dataDict[Dkey].update({'Number of contours':0, \
                                  'Contour number':[], \
                                  'Number of contour points':0, \
                                  'Contour data':[]
                                  })
        
            
    # Create a new array "contourPts" containing the (x,y) image 
    # coordinates of each contour point and save it in dataDict:
    
    # Note:
    # Conversion of contour points from patient coordinates to image coordinates
    # adapted from James Petts' (ICR) Java code: 
    
    for uid, values in dataDict.items():
        # Get the values for this uid:
        # Note that for any given uid, NcontourPts and contourData could
        # be lists of values, if multiple contours existed for that uid.  
        # Even if there was only one contour, the single value (or single 
        # list) will be enclosed in [].
        position = values['Position'] # e.g. [x, y, z]
        orientation = values['Orientation'] # e.g. [row_x, row_y, row_z, col_x, col_y, col_z]
        pixSpacing = values['Pixel spacing'] # e.g. [dx, dy]
        Ncontours = values['Number of contours'] # e.g. 1 or 2, etc.
        contourData = values['Contour data'] # e.g. [[x1,y1,z1,...]] or [[x1,y1,z1,...], [x1,y1,z1,...]], etc.
    
        # The patient coordinates:
        x0 = position[0]
        y0 = position[1]
        z0 = position[2]
        
        # The direction cosines along rows and columns:
        #row_x, row_y, row_z, col_x, col_y, col_z = orientation
        theta_r = orientation[0:3]
        theta_c = orientation[3:6]
        
        # The indeces of the largest direction cosines along rows and columns:
        ind_r = np.where(theta_r==max(theta_r))[0][0]
        ind_c = np.where(theta_c==max(theta_c))[0][0]
        
        # The pixel spacings:
        dx = pixSpacing[0]
        dy = pixSpacing[1]
        
        # Create an array of (x,y) points from the flattened contour data:
        contourPts = []
        
        # Iterate for each contour for this uid:
        for c in range(Ncontours):
            # Store the contour points for each contour within any given uid:
            contourPts_c = []
        
            # Iterate for all countour points in threes:
            for p in range(0, len(contourData[c]), 3):
                # Define vector v as the the patient position subtracted from 
                # the contour coordinates:
                v = contourData[c][p] - x0, \
                contourData[c][p+1] - y0, \
                contourData[c][p+2] - z0
                
                # Convert from patient coord system to image coord system:          
                i = (v[ind_r] - theta_c[ind_r]*v[ind_c]/theta_c[ind_c]) / \
                (theta_r[ind_r]*dx * (1 - (theta_c[ind_r]*theta_r[ind_c])/(theta_r[ind_r]*theta_c[ind_c])))
                
                j = (v[ind_c] - theta_r[ind_c]*i*dx) / (theta_c[ind_c]*dy)
    
                contourPts_c.append([i,j])
                
            # Append contourPts with contourPts_c:
            contourPts.append(contourPts_c)
    
        # Add contourPts to dictData for this uid (key):
        dataDict[uid].update({'Contour points':contourPts})
    
    
    # Plot the DICOM scans with contours overlayed:
    # Configure plot:

    # Set the title font:
    fontSize=4
    
    Nuids = len(dataDict.keys())
    
    # Set the number of subplot rows and columns:
    rows = np.int8(np.round(np.sqrt(Nuids)))
    cols = np.int8(np.round(np.sqrt(Nuids)))
    
    # Decide whether to flip the images up-down so that the orientation is
    # the same as shown in the OHIF-Viewer:
    #flip = False
    flip = True
    
    # Plot results:
    
    plt.figure(figsize=(15,15), dpi=300);
    
    i = 0 # for subplot pos
    
    for uid, values in dataDict.items(): 
        frameNo = values['Frame number']
        pixArray = values['Pixel array']
        Ncontours = values['Number of contours']
        contourPts = values['Contour points']
        
        # Plot the pixel array:
        i = i + 1    
        plt.subplot(rows,cols,i, aspect='equal')
        if flip:
            plt.pcolormesh(np.flipud(pixArray));
        else:
            plt.pcolormesh(pixArray);
        
        plt.title(f'Frame {frameNo}', size=fontSize);
        plt.axis('off');
        
        # Plot the contours:
        if Ncontours > 0:   
            for c in range(Ncontours):
                # Unpack tuple and store each x,y tuple in arrays X and Y:
                X = []
                Y = []
                for x, y in contourPts[c]:
                    X.append(x)
                    Y.append(y)
    
                # Create 2D arrays for the X and Y values using meshgrid:
                #X, Y = np.meshgrid(xList, yList)
                
                # Flip Y pixel coordinates if flip=True:
                if flip:
                    N_r, N_c = np.shape(pixArray)
                    
                    Y = N_c - np.array(Y)
    
                plt.plot(X, Y, linewidth=0.5, c='red');
    
    
    
    if exportPlot:
        # Set up details for exporting of plot:
        #seriesDescription = dicom.SeriesDescription
        
        
        print('Study description =', dicom.StudyDescription)
        print('Patient ID =', dicom.PatientID)
        print('Series description =', dicom.SeriesDescription)
        
        
        exportFname = dicom.StudyDescription + ' ' + dicom.PatientID + ' ' \
        + dicom.SeriesDescription + '.jpg'
        
        exportFpath = os.path.join(exportDir, exportFname)
        
        # Export the plot:
        plt.savefig(exportFpath, bbox_inches='tight')
        
        print('Plot exported to:\n\n', exportFpath)
    
    return dataDict
    #return