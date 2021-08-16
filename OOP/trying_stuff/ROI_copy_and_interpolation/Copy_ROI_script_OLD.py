#!/usr/bin/env python
# coding: utf-8

# ### Main script:

# In[ ]:


# Import packages and functions
import os
#import importlib
#import ContourInterpolatingFuncs
#importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

# Define the root directory of the data:
RootDataDir = r'C:\Code\WP1.3_multiple_modalities\Contour_interpolation\data\ACRIN-FMISO-Brain-011'

# Define FixedDicomDir and MovingDicomDir:
FixDicomDir = os.path.join(RootDataDir, 
                           r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(RootDataDir, 
                           r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = os.path.join(RootDataDir, 
                           r'04-10-1960-MRI Brain wwo Contrast-69626\ROI_Collection')
FixRoiFname = r'AIM_20200511_073405.dcm'
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Decide whether or not to interpolate between all contour point pairs or only 
# those for which the node of the longest list of contour points was an 
# original node (see Note 4 in the function InterpolateBetweenContours):
InterpolateAllPts = True

# Decide whether or not to only consider intersection points that line on the
# line segment:
OnLineSegment = True

# Decide whether to log some results to the console:
LogToConsole = True

# Decide whether to plot results to the console:
PlotResults = True
#PlotResults = False

# Decide whether to annotate points with their numbers:
AnnotatePtNums = True
AnnotatePtNums = False

# Decide whether to export plotted results to file:
ExportResults = True
ExportResults = False


PointData, FixContourData,InterpData, MovContourData = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                          dP, InterpolateAllPts, LogToConsole, 
                                          PlotResults, AnnotatePtNums, ExportResults)


# ### 3D plots:

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'notebook')

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True

cif.PlotInterpolatedContours3D(InterpData=InterpData, 
                               FixContourData=FixContourData, 
                               MovContourData=MovContourData,
                               FixDicomDir=FixDicomDir,
                               MovDicomDir=MovDicomDir,
                               dP=dP, 
                               PlotImagingPlanes=PlotImagingPlanes,
                               CombinePlots=CombinePlots, ExportPlot=ExportPlot)

