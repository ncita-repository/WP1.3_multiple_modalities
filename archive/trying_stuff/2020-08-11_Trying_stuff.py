#!/usr/bin/env python
# coding: utf-8

# In[64]:


# Import packages and functions
import SimpleITK as sitk
print(sitk.Version())
#import numpy as np
import time, os, copy
import os
import importlib
#import pydicom
#import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')
#from ipywidgets import interact, fixed
#from IPython.display import clear_output

#import RegUtilityFuncs
#importlib.reload(RegUtilityFuncs)
#import GetInputPoints
#importlib.reload(GetInputPoints)
import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)

#import GetImageAttributes
#importlib.reload(GetImageAttributes)
#from GetImageAttributes import GetImageAttributes
#from GetInputPoints import GetInputPoints
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

# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm' # tumour
#FixRoiFname = r'AIM_20200626_104631.dcm' # ventricles
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)


# Get a list of a list of points of the four points at the corners of the image plane
# for all imaging planes in FixImage:
FixPlanes = cif.GetPointsInImagePlanes(DicomDir=FixDicomDir)

#FixPlanes

# Create meshgrids for the image planes:
#FixMeshgrids = cif.GetMeshGridsForImagePlanesInContours(Contours=FixPlanes, 
#                                                        OnlyPlanesWithContours=True, 
#                                                        DicomDir=FixDicomDir)

#Meshgrids

Contours = FixPlanes
Convert2ICS = False
PlotImagingPlanes = True
DicomDir = FixDicomDir
PlotTitle = 'Points on Fixed image planes and Fixed image planes'
ExportPlot = False
#ExportPlot = True

cif.PlotContoursAndPlanes3D(Contours=Contours, Convert2ICS=Convert2ICS, 
                            PlotImagingPlanes=PlotImagingPlanes, DicomDir=DicomDir,
                            PlotTitle=PlotTitle, ExportPlot=ExportPlot)


# In[24]:


FixPlanes


# In[25]:


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

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Register:
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm, LogToConsole=True)


# In[28]:


FixPlanes
#cif.GetPointsInImagePlanes(DicomDir=FixDicomDir)


# ### Transform points slice (contour) by slice (contour) (VERY INEFFICIENT):

# In[31]:


#import TransformPoints
#importlib.reload(TransformPoints)
#from TransformPoints import TransformPoints
import time

# Start timing:
times = []
times.append(time.time())

# Initialise the transformed FixPlanes:
FixPlanesT = []

# Useful info:
NumOfPts = 0
NumOfFiles = 0

# Loop through each list of points in FixPlanes:
for s in range(len(FixPlanes)):
    
    Points = FixPlanes[s]
    
    # If Points is not []:
    if Points:

        # Create inputpoints.txt for Elastix, containing the contour points to be
        # transformed:
        #CreateInputFileForElastix(Points=Points)

        # Transform the input points:
        #TransformInputPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

        # Transform the points:
        TransformPoints(Points=Points, Image=MovIm, TransformFilter=ElastixImFilt)

        # Parse outputpoints.txt:
        PtNos, InInds, InPts_PCS, FixOutInds,        OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

        #PtNos, InInds, InPts_PCS, InPts_ICS,\
        #FixOutInds, OutPts_PCS, OutPts_ICS,\
        #Defs_PCS, Defs_ICS, MovOutInds,\
        #OutPtsArr_PCS, OutPtsArr_ICS = GetOutputPoints(Origin=MovOrigin, Directions=MovDirs,
        #                                               Spacings=MovSpacings, Dimensions=MovDims)
        
        FixPlanesT.append(OutPts_PCS)
        
        NumOfPts += len(PtNos)
        
        NumOfFiles += 1
    
    else:
        FixPlanes.T.append([])


        
times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)

print(f'Took {Dtime} s to transform {NumOfPts} points in {NumOfFiles} files.\n\n')

""" Took 28.7 s to transform 150 points in 30 files. """

FixPlanesT


# ### Transform all points in a single operation (MUCH MORE EFFICIENT):

# In[44]:


#import TransformPoints
#importlib.reload(TransformPoints)
#from TransformPoints import TransformPoints
import time

# Start timing:
times = []
times.append(time.time())

# Initialise a list of all points (along rows, no grouping by slice index):
AllPoints = []

# Initialise list of the number of points contained per slice (so that 
# transformed points can be re-grouped by the original slice index):
NumPtsPerSlice = []

# Loop through each list of points in FixPlanes:
for s in range(len(FixPlanes)):
    
    Points = FixPlanes[s]
    
    # If Points is not []:
    if Points:
        #AllPoints.append(Points)
        AllPoints.extend(Points)
        
        NumPtsPerSlice.append(len(Points))

        
# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
#CreateInputFileForElastix(Points=Points)

# Transform the input points:
#TransformInputPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Transform the points:
TransformPoints(Points=AllPoints, Image=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
PtNos, InInds, InPts_PCS, FixOutInds,OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

#PtNos, InInds, InPts_PCS, InPts_ICS,\
#FixOutInds, OutPts_PCS, OutPts_ICS,\
#Defs_PCS, Defs_ICS, MovOutInds,\
#OutPtsArr_PCS, OutPtsArr_ICS = GetOutputPoints(Origin=MovOrigin, Directions=MovDirs,
#                                               Spacings=MovSpacings, Dimensions=MovDims)

#print(f'{len(AllPoints)} points were transformed to {len(OutPts_PCS)} points\n\n')

# Initialise the transformed FixPlanes:
FixPlanesT = []

## Keep track of the number of points that are re-grouped
## at each loop:
#NumPtsRegrouped = 0

# Useful info:
NumOfPts = 0

for s in range(len(FixPlanes)):
    # If NumPtsPerSlice[s] is not 0:
    if NumPtsPerSlice[s]:
        PointsThisSlice = []
        
        #print(f's = {s}, NumPtsPerSlice[{s}] = {NumPtsPerSlice[s]}')
        #print(f'p will loop from {NumOfPts} to {NumOfPts + NumPtsPerSlice[s]}\n\n')
        
        for p in range(NumOfPts, NumOfPts + NumPtsPerSlice[s]):
            PointsThisSlice.append(OutPts_PCS[p])
        
        FixPlanesT.append(PointsThisSlice)
        
        NumOfPts += len(PointsThisSlice)
    
    else:
        FixPlanes.T.append([])


        
times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)

print(f'Took {Dtime} s to transform {NumOfPts} points in 1 file.\n\n')

""" 
Took 28.7 s to transform 150 points in 30 files. 
Took 0.7 s to transform 150 points in 1 file.
"""

FixPlanesT


# # => All points should be transformed in one operation due to the inefficiency of file I/O operations!!!

# In[65]:


Contours = FixPlanesT
Convert2ICS = False
PlotImagingPlanes = True
DicomDir = MovDicomDir
PlotTitle = 'Transformed points on Fixed image planes and Moving image planes'
ExportPlot = False
#ExportPlot = True

cif.PlotContoursAndPlanes3D(Contours=Contours, Convert2ICS=Convert2ICS, 
                            PlotImagingPlanes=PlotImagingPlanes, DicomDir=DicomDir,
                            PlotTitle=PlotTitle, ExportPlot=ExportPlot)


# In[66]:


# Plot points from the Fixed image planes onto the Moving image planes:
Contours = FixPlanes
Convert2ICS = False
PlotImagingPlanes = True
DicomDir = MovDicomDir
PlotTitle = 'Points on Fixed image planes and Moving image planes'
ExportPlot = False
#ExportPlot = True

cif.PlotContoursAndPlanes3D(Contours=Contours, Convert2ICS=Convert2ICS, 
                            PlotImagingPlanes=PlotImagingPlanes, DicomDir=DicomDir,
                            PlotTitle=PlotTitle, ExportPlot=ExportPlot)


# In[67]:


# Plot points from the Moving image planes onto the Fixed image planes:
Contours = FixPlanesT
Convert2ICS = False
PlotImagingPlanes = True
DicomDir = FixDicomDir
PlotTitle = 'Transformed points on Fixed image planes and Fixed image planes'
ExportPlot = False
#ExportPlot = True

cif.PlotContoursAndPlanes3D(Contours=Contours, Convert2ICS=Convert2ICS, 
                            PlotImagingPlanes=PlotImagingPlanes, DicomDir=DicomDir,
                            PlotTitle=PlotTitle, ExportPlot=ExportPlot)


# In[ ]:




