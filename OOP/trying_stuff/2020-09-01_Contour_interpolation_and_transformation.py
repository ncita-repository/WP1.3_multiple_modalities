#!/usr/bin/env python
# coding: utf-8

# In[22]:


# Import packages and functions
#import SimpleITK as sitk
#print(sitk.Version())
#import numpy as np
#import time, os, copy
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
#from CreateInputFileForElastix import CreateInputFileForElastix
#from RunSimpleElastixReg import RunSimpleElastixReg
#from TransformPoints import TransformPoints
#from ParseTransformixOutput import ParseTransformixOutput
#from GetOutputPoints import GetOutputPoints
#from PCStoICS import PCStoICS
#import RegUtilityFuncs as ruf
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

""" Parameters for interpolation: """

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Choose whether or not to interpolate between all contour point pairs or only 
# those for which the node of the longest list of contour points was an 
# original node (see Note 4 in the function InterpolateBetweenContours):
InterpolateAllPts = True

""" Parameters for intersecting points: """

# Choose whether or not to only consider intersection points that line on the
# line segment:
OnLineSegment = True

# Choose whether or not to use the interpolated contour points when searching
# for intersection points:
UseInterp = True
UseInterp = False


""" Parameters for logging to console and plotting: """

# Choose whether to log some results to the console:
LogToConsole = True

# Choose whether to plot results to the console:
PlotResults = True
#PlotResults = False

# Choose whether to annotate points with their numbers:
AnnotatePtNums = True
AnnotatePtNums = False

# Choose whether to export plotted results to file:
ExportResults = True
ExportResults = False


PointData, FixContourData,InterpData, MovContourData, FixIm, MovIm, RegIm = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, dP,
                                   InterpolateAllPts, UseInterp, LogToConsole,
                                   PlotResults, AnnotatePtNums, ExportResults)

"""
Previously took 18 s to register and transform points row-by-row and column-by-column in 
InterpData (i.e. multiple file I/Os).

Now takes 9 s to register the images, and <1 s to transform all points.
"""


# import RegUtilityFuncs
# importlib.reload(RegUtilityFuncs)
# import RegUtilityFuncs as ruf
# import time
# 
# Perspective = 'axial'
# 
# PlotImagingPlanes = False
# #PlotImagingPlanes = True
# 
# ExportFig = True
# ExportFig = False
# 
# ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) \
#               + f'_Reg_results_with_Fixed_contours_only__dP_{dP}__' \
#               + f'UseInterp_{UseInterp}__' + Perspective + '.png'
# 
# ruf.PlotFixMovRegImagesAndContours(FixImage=FixIm, MovImage=MovIm, RegImage=RegIm, 
#                                    FixContours=FixContourData['PointICS'], 
#                                    MovContours=MovContourData['PointICS'],
#                                    SlicesToPlot=[13,14,15,16,17,18], 
#                                    Perspective=Perspective, ContoursAs='lines', 
#                                    LogToConsole=False, ExportFig=ExportFig, 
#                                    ExportFname=ExportFname)

# In[53]:


FixContourData['PointICS']

#pd.DataFrame(FixContourData)


# import ContourInterpolatingFuncs
# importlib.reload(ContourInterpolatingFuncs)
# import ContourInterpolatingFuncs as cif
# 
# %matplotlib notebook
# 
# Convert2ICS = True
# Convert2ICS = False
# 
# PlotInterpContours = True
# #PlotInterpContours = False
# 
# PlotImagingPlanes = True
# #PlotImagingPlanes = False
# 
# CombinePlots = True
# CombinePlots = False
# 
# ExportPlot = False
# #ExportPlot = True
# 
# PlotTitle = 'Fixed'
# 
# ## Ignore interpolated contours
# #Contours = FixContourData['PointICS'][i] for i in range(29)
# 
# ##cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
# ##cif.PlotInterpolatedContours3D(InterpData=InterpData,
# #cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
# #                                  FixContourData=FixContourData, 
# #                                  MovContourData=MovContourData,
# #                                  Convert2ICS=Convert2ICS,
# #                                  dP=dP,
# #                                  PlotImagingPlanes=PlotImagingPlanes,
# #                                  FixDicomDir=FixDicomDir,
# #                                  MovDicomDir=MovDicomDir,
# #                                  CombinePlots=CombinePlots, 
# #                                  ExportPlot=ExportPlot)
# 
# ##cif.PlotContoursAndPlanes3D(Contours=FixContourData['PointICS'], 
# #cif.PlotContoursAndPlanes3D(Contours=FixContourData['PointPCS'], 
# #                            Convert2ICS=Convert2ICS, 
# #                            PlotImagingPlanes=PlotImagingPlanes, 
# #                            DicomDir=FixDicomDir, 
# #                            PlotTitle=PlotTitle, 
# #                            ExportPlot=ExportPlot)
# 
# cif.PlotContoursAndPlanes3DNEW(ContourData=FixContourData,
#                                PlotInterpContours=PlotInterpContours,
#                                PlotImagingPlanes=PlotImagingPlanes,
#                                Convert2ICS=Convert2ICS,
#                                DicomDir=FixDicomDir, 
#                                PlotTitle=PlotTitle, 
#                                ExportPlot=ExportPlot)

# In[126]:


import pandas as pd

pd.DataFrame(FixContourData)


# In[180]:


import pandas as pd

pd.DataFrame(MovContourData)


# In[181]:


pd.DataFrame(InterpData)


# In[145]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif


FixOrMov = 'Fix'
#FixOrMov = 'Mov'

PlotImagingPlanes = True
#PlotImagingPlanes = False

PlotAllImagePlanes = True
PlotAllImagePlanes = False

ConvertMeshgrids2ICS = True
ConvertMeshgrids2ICS = False

ExportPlot = False
#ExportPlot = True

if FixOrMov == 'Fix':
    PlotTitle = 'Fixed'
    ContourData = FixContourData
    DicomDir = FixDicomDir
    
if FixOrMov == 'Mov':
    PlotTitle = 'Moving'
    ContourData = MovContourData
    DicomDir = MovDicomDir

cif.PlotContoursAndPlanes3D_NEW_2(InterpData=InterpData, 
                                  ContourData=ContourData, 
                                  FixOrMov=FixOrMov,
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  PlotAllImagePlanes=PlotAllImagePlanes,
                                  ConvertMeshgrids2ICS=ConvertMeshgrids2ICS,
                                  DicomDir=DicomDir,
                                  PlotTitle=PlotTitle, 
                                  ExportPlot=ExportPlot)


# In[144]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif


FixOrMov = 'Fix'
FixOrMov = 'Mov'

PlotImagingPlanes = True
#PlotImagingPlanes = False

PlotAllImagePlanes = True
PlotAllImagePlanes = False

ConvertMeshgrids2ICS = True
ConvertMeshgrids2ICS = False

ExportPlot = False
#ExportPlot = True

if FixOrMov == 'Fix':
    PlotTitle = 'Fixed'
    ContourData = FixContourData
    DicomDir = FixDicomDir
    
if FixOrMov == 'Mov':
    PlotTitle = 'Moving'
    ContourData = MovContourData
    DicomDir = MovDicomDir

cif.PlotContoursAndPlanes3D_NEW_2(InterpData=InterpData, 
                                  ContourData=ContourData, 
                                  FixOrMov=FixOrMov,
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  PlotAllImagePlanes=PlotAllImagePlanes,
                                  ConvertMeshgrids2ICS=ConvertMeshgrids2ICS,
                                  DicomDir=DicomDir,
                                  PlotTitle=PlotTitle, 
                                  ExportPlot=ExportPlot)


# In[ ]:




