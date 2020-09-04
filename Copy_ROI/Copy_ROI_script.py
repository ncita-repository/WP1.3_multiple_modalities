#!/usr/bin/env python
# coding: utf-8

# # Guidance
# 
# 
# ## Data
# 
# The data used is from Subject 011 from the TCIA ACRIN-FMISO-Brain data set obtained here:
# 
# https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948305#a6b63241cb4843a291c45ef8a62d2073
# 
# You can selectively download data for that subject only.  Specifically two series from two sessions were used:
# 
# 8-T1 SE AXIAL POST FS FC-59362 (Series 8) within 04-10-1960-MRI Brain wwo Contrast-69626
# 
# and
# 
# 8-T1 SE AXIAL POST FC-81428 (Series 8) within 06-11-1961-MRI Brain wwo Contrast-79433
# 
# An ROI Collection was manually created within the OHIF-Viewer for the series from 1960.  
# 
# All above data is located here:
# 
# WP1.3_multiple_modalities\Copy_ROI\data\ACRIN-FMISO-Brain-011
# 
# in separate sub-folders for the "fixed" and "moving" DICOMs, and for the ROI Collection that relates to the "fixed" DICOMs.
# 
# In the Main script below you may have to modify the directory path for "RootDataDir".
# 
# 
# 
# ## SimpleElastix
# 
# I recommend trying to install SimpleElastix using a pre-compiled Python module obtained here:
# 
# https://www.francescosantini.com/wp/2020/03/06/simpleelastix-precompiled-python-modules/
# 
# If it doesn't work you will have to compile it from source:
# 
# https://simpleelastix.readthedocs.io/GettingStarted.html
# 
# There's a "more official" wheel provided here:
# 
# https://sourceforge.net/projects/simpleelastix/
# 
# but it might be for an older version of SimpleElastix (and I think I had trouble with this wheel and ended up having to compile it - I hadn't yet discovered Francesco Santini's page).
# 
# 
# 
# ## Other packages
# 
# Most of the other packages can be pip installed straightforwardly.  
# 
# 
# 
# ## Description of approach
# 
# The contours are read in and stored in a dictionary containing S lists, for each slice in the fixed DICOM stack. Slices without contours have empty lists ([]).  
# 
# Every two adjacent contours are interpolated - firstly by interpolating within each contour (i.e. up-sampling) so that a minimum inter-node threshold is met.  Then each up-sampled contour is down-sampled accordingly so that both pairs of contours have an equal number of points (i.e. "over-sampled" contours).
# 
# Each pair of over-sampled contours are interpolated between them.  Hence there are a trio of contours all containing the same number of points:  A pair of over-sampled contours and their interpolated-between contour.  This process is repeated for every pair of adjacent contours.
# 
# The result is a dictionary ("InterpData") that contains the following groups of points:
# 
# Contour1, Contour2-1, Interp1-2
# Contour2-2, Contour3-1, Interp2-3
# Contour3-2, Contour4-1, Interp3-4
# ...
# 
# where Interp1-2 is the contour interpolated between Contour1 and Contour2, and Contour2-1 and Contour2-2 represent the same original contour (Contour2), but with different sampling such that Contour2-1 has the same number of contour points as Contour1, and likewise Contour2-2 has the same number of points as Contour3-1.  
# 
# The moving image stack is registered to the fixed moving stack.  All points from all trios of contours are transformed using the registration transformation.
# 
# The automatically generated outputpoints.txt file is parsed to collect the indices and coordinates of the transformed points (and also the indices of the input points).
# 
# For each pair of transformed contours, each point in the pair define a line segment, whose intersection with all planes in the moving image domain are found.  Only intersections that lie on the line segments are considered.
# 
# The result is a list of intersection points, one list for each slice in the moving image domain.  Each list is comprised of a sub-list of points, whereby the points in any given sub-list are the intersections that came about from different pairs of contours.  
# 
# In the next version, all contours will be re-sampled such that they all have an equal number of points (rather than an equal number amongst two adjacent contours only).  And rather than find the intersection between the planes and straight line segments that join two adjacent contours, instead the N^th point in all contours will be joined using a spline (or similar) fitting, and the intersection of that curve with the planes will be determined.

# # Main script:

# In[2]:


# Import packages and functions
import os
#import importlib
#import ContourInterpolatingFuncs
#importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')

# Define the root directory of the data:
RootDataDir = r'C:\Code\WP1.3_multiple_modalities\Copy_ROI\data\ACRIN-FMISO-Brain-011'

# Define FixedDicomDir and MovingDicomDir:
FixDicomDir = os.path.join(RootDataDir, 'Fixed_DICOMs')
MovDicomDir = os.path.join(RootDataDir, 'Moving_DICOMs')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = os.path.join(RootDataDir, 'Fixed_ROI_Collection')
FixRoiFname = r'AIM_20200511_073405.dcm'
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


PointData, FixContourData,InterpData, MovContourData = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                            dP, InterpolateAllPts, UseInterp,
                                            LogToConsole, PlotResults, 
                                            AnnotatePtNums, ExportResults)


# In[ ]:




