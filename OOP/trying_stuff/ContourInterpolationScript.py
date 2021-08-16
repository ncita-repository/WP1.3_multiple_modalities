# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:00:16 2020

@author: ctorti
"""


PlotResults2D = True
PlotResults2D = False

PlotResults3D = True
#PlotResults3D = False



# Import packages and functions
import SimpleITK as sitk
print(sitk.Version())
import numpy as np
import time, os, copy
import importlib
import pydicom
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib inline
#%matplotlib notebook
#from ipywidgets import interact, fixed
#from IPython.display import clear_output

import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import GetInputPoints
importlib.reload(GetInputPoints)
import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)

#import GetImageAttributes
#importlib.reload(GetImageAttributes)
from GetImageAttributes import GetImageAttributes
from GetInputPoints import GetInputPoints
#import CreateInputFileForElastix
#importlib.reload(CreateInputFileForElastix)
from CreateInputFileForElastix import CreateInputFileForElastix
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
from ParseTransformixOutput import ParseTransformixOutput
from GetOutputPoints import GetOutputPoints
from PCStoICS import PCStoICS
import RegUtilityFuncs as ruf
import ContourInterpolatingFuncs as cif 

from GetVectorLength import GetVectorLength
from NormaliseVector import NormaliseVector


# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm' # tumour
#FixRoiFname = r'AIM_20200626_104631.dcm' # ventricles
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose which package to use to get the Image Plane Attributes:
#package = 'sitk'
package = 'pydicom'


# Get the Image Attributes for the images:
FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                              Package=package)

MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                              Package=package)


# Read in the 3D stack of Fixed DICOMs:
FixReader = sitk.ImageSeriesReader()
FixNames = FixReader.GetGDCMSeriesFileNames(FixDicomDir)
FixReader.SetFileNames(FixNames)
FixIm = FixReader.Execute()
FixIm = sitk.Cast(FixIm, sitk.sitkFloat32)

# Read in the 3D stack of Moving DICOMs:
MovReader = sitk.ImageSeriesReader()
MovNames = MovReader.GetGDCMSeriesFileNames(MovDicomDir)
MovReader.SetFileNames(MovNames)
MovIm = MovReader.Execute()
MovIm = sitk.Cast(MovIm, sitk.sitkFloat32)

# Get some info on FixIm and MovIm:
#ruf.ShowImagesInfo(FixIm, MovIm)

# Register MovIm to FixIm:
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Get contour points into necessary arrays:
FixPtsPCS, FixPtsBySliceAndContourPCS,\
FixPtsICS, FixPtsBySliceAndContourICS,\
LUT, PointData, ContourData = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                                             Origin=FixOrigin, Directions=FixDirs,
                                             Spacings=FixSpacings)




import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

# Define the slice(s) to be interpolated:
#InterpSliceInd = 14
#InterpSliceInd = 13.5

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Get the indices of the slices that contain contours:
SlicesWithContours = cif.GetIndsOfSlicesWithContours(ContourData=ContourData)

print('Contours exist on slices', SlicesWithContours)

# Initialise a dictionary to store the interpolation results:
InterpData = {'InterpSliceInd'       : [],
              'BoundingSliceInds'    : [],
              'FixOSIsOrigNode1'     : [],
              'FixOSIsOrigNode2'     : [],
              'FixOSContour1Inds'    : [],
              'FixOSContour2Inds'    : [],
              'FixInterpContourInds' : [],
              'FixOSContour1Pts'     : [],
              'FixOSContour2Pts'     : [],
              'FixInterpContourPts'  : []
              }

# Iterate through all SlicesWithContours and interpolate between each
# two slices:
for i in range(len(SlicesWithContours) - 1):
    InterpSliceInd = SlicesWithContours[i] + 0.5
    
    print(f'Attempting to interpolate at slice {InterpSliceInd}...')


    #InterpContour = cif.InterpolateContours(ContourData, PointData, InterpSliceInd, dP)
    
    BoundingSliceInds,\
    OSContour1,\
    OSIsOrigNode1,\
    OSContour2,\
    OSIsOrigNode2,\
    InterpContour = cif.InterpolateContours(ContourData, PointData, InterpSliceInd, dP)

    # Continue only if None was not returned for InterpContour:
    if InterpContour:
        #print(f'len(InterpContour) = {len(InterpContour)}')

        # The maximum contour number in the original contour data:
        LastCntNo = PointData['InContourNo'][-1]

        # Add interpolated contour to PointData:
        for i in range(len(InterpContour)):
            PointData['PointNo'].append(PointData['PointNo'][-1] + 1)
            PointData['InSliceNo'].append(InterpSliceInd)
            PointData['InContourNo'].append(LastCntNo + 1)
            PointData['ContourType'].append(2)
            PointData['InPointIndex'].append([])
            PointData['InPointPCS'].append(InterpContour[i])
            PointData['InPointICS'].append(PCStoICS(Pts_PCS=InterpContour[i], 
                                                    Origin=FixOrigin, 
                                                    Directions=FixDirs, 
                                                    Spacings=FixSpacings))

        # Add interpolated contour to ContourData:
        ContourData['InSliceNo'].append(InterpSliceInd)
        ContourData['ContourType'].append(2)
        ContourData['InPointPCS'].append(InterpContour)
        ContourData['InPointICS'].append(PCStoICS(Pts_PCS=InterpContour, 
                                                  Origin=FixOrigin, 
                                                  Directions=FixDirs, 
                                                  Spacings=FixSpacings))
        
        # Store the interpolation output in InterpData:
        InterpData['InterpSliceInd'].append(InterpSliceInd)
        InterpData['BoundingSliceInds'].append(BoundingSliceInds)
        InterpData['FixOSIsOrigNode1'].append(OSIsOrigNode1)
        InterpData['FixOSIsOrigNode2'].append(OSIsOrigNode2)
        InterpData['FixOSContour1Pts'].append(OSContour1)
        InterpData['FixOSContour2Pts'].append(OSContour2)
        InterpData['FixInterpContourPts'].append(InterpContour)


        """ GET ADDITIONAL DATA FOR PLOTTING """


        # Get the bounding slice indices for the slice at InterpSliceInd:
        BoundingSliceInds = cif.GetBoundingSliceInds(ContourData, PointData, InterpSliceInd)

        Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
        Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])

        if InterpSliceInd.is_integer():
            # The actual contour that exists at slice InterpSliceInd:
            ActualContour = copy.deepcopy(ContourData['InPointPCS'][InterpSliceInd][0])




        """ PLOT RESULTS """
        
        if PlotResults2D:

            ExportPlot = False
            #ExportPlot = True
    
            Contours = [Contour1, Contour2, InterpContour]
            Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
                      f'Contour at slice {BoundingSliceInds[1]}', 
                      f'Interpolated contour at {InterpSliceInd}']
            Colours = ['b', 'g', 'r', 'y']
            #Shifts = [0, 2, 4, 6] # to help visualise
            Shifts = [0, 0, 0, 0] # no shift
            if InterpSliceInd.is_integer():
                Contours.append(ActualContour)
                Labels.append(f'Actual contour at {InterpSliceInd}')
    
            cif.PlotInterpolationResults(Contours=Contours, Labels=Labels, 
                                         Colours=Colours, Shifts=Shifts,
                                         InterpSliceInd=InterpSliceInd,
                                         BoundingSliceInds=BoundingSliceInds,
                                         dP=dP, ExportPlot=ExportPlot)
            
            
            
            
# Transform the contour points in InterpData:
# Create lists of keys that can be looped through:
BaseKeys = ['OSContour1', 'OSContour2', 'InterpContour']

FixIndsKeys = ['Fix' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
FixPtsKeys = ['Fix' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

MovIndsKeys = ['Mov' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
MovPtsKeys = ['Mov' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

# Loop through each row in InterpData and transform the points in FixPtsKeys:
for i in range(len(InterpData['InterpSliceInd'])):
    
    
    for k in range(len(FixPtsKeys)):
        Points = InterpData[FixPtsKeys[k]][i]
        
        # Create inputpoints.txt for Elastix, containing the contour points to be
        # transformed:
        CreateInputFileForElastix(Points=Points)

        # Transform MovingContourPts:
        TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

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
            


if PlotResults3D:
    import ContourInterpolatingFuncs
    importlib.reload(ContourInterpolatingFuncs)
    import ContourInterpolatingFuncs as cif
    
    ExportPlot = False
    #ExportPlot = True
    
    cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=ExportPlot)

