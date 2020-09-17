#!/usr/bin/env python
# coding: utf-8

# In[1]:


from ImportDictionaryFromJson import ImportDictionaryFromJson as ImportDict
#import importlib
#import ContourInterpolatingFuncs
#importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 
import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline
#%matplotlib notebook
#from ImportNumpyArray import ImportNumpyArray as ImportNpa
#from ImportImageAttributesFromJson import ImportImageAttributesFromJson as ImportImAtt

# First copy the files Simple_ContourData_pre.json, Simple_ContourData_post.json 
# and Simple_ContourSweepData.json to the working directory.

# Import the contour data dictionaries:

ContourDataPre = ImportDict(FileName='Simple_ContourData_pre')
ContourDataPost = ImportDict(FileName='Simple_ContourData_post')
ContourSweepData = ImportDict(FileName='Simple_ContourSweepData')

# Decide on the following parameters:

""" Parameters for intersecting points: """

# Choose whether or not to only consider intersection points that line on the
# line segment:
MustBeOnLineSeg = True
#MustBeOnLineSeg = False

# Set a limit on the maximum distance between the intersection point and 
# contours points that define the line segment:
MaxDistToPts = 5
MaxDistToPts = 2.5
MaxDistToPts = 1


""" Parameters for logging to console and plotting: """

# Choose whether to log some results to the console:
LogToConsole = True
#LogToConsole = False

# Choose whether to plot results to the console:
PlotResults = True
PlotResults = False

# Choose whether to annotate points with their numbers:
AnnotatePtNums = True
AnnotatePtNums = False

# Choose whether to export plotted results to file:
ExportResults = True
ExportResults = False

""" end of list of parameters """


# Assign Method variable (not yet implemented in the function below):
Method = ''

# Get the intersecting points of the contours on the moving planes: 
ContourData = cif.GetIntersectingPtsFromContourSweepData2(ContourSweepData=ContourSweepData, 
                                                          MustBeOnLineSeg=MustBeOnLineSeg, 
                                                          MaxDistToPts=MaxDistToPts,
                                                          ContourData=ContourDataPre)


# In[2]:


# Display ContourData as a pandas DataFrame:

pd.DataFrame(ContourData)


# In[2]:


# Plot images and contours:

# Ensure that the files FixIm.npy, MovIm.npy, RegIm.npy, 
# FixOrigin.json, FixDirs.json, FixSpacings.json, FixDims.json, 
# MovOrigin.json, MovDirs.json, MovSpacings.json and MovDims.json 
# have been copied to the working directory.

import ImagePlottingFuncs
import importlib
importlib.reload(ImagePlottingFuncs)
from ImagePlottingFuncs import PlotImagesAndContours
from ContourInterpolatingFuncs import GetIndsOfNonEmptyItems
from PCStoICS import PCStoICSnested
from ImportImageAttributesFromJson import ImportImageAttributesFromJson as ImportImAtt
import time

# Decide on the following parameters:

""" Parameters for plotting: """

# Choose which set of Moving (intersection) points to plot:
#MovKey = 'MovSSPtsPCS' # super-sampled points
MovKey = 'MovOSPtsPCS' # over-sampled points

# Choose which slices to plot:
#SlicesToPlot = [17, 18, 19]
#SlicesToPlot = -1 # plot all slices
SlicesToPlot = -2 # plot only slices containing contours

# Choose the perspective:
Perspective = 'axial'
#Perspective = 'sagittal'
#Perspective = 'coronal'

# Chose whether to keep plot axes in units of pixels or mm:
UsePCSforAxes = False # pixels
#UsePCSforAxes = True # mm <-- this isn't working properly yet

# Choose whether to plot contour points as lines or dots:
ContoursAs = 'lines'
#ContoursAs = 'dots'

LogToConsole = True
LogToConsole = False

ExportResults = True
ExportResults = False

""" end of list of parameters """



if UsePCSforAxes:
    FixContours = ContourData['FixPtsPCS']
    
    MovContours = ContourData[MovKey]
    
else:
    # Import the image attributes data:
    FixOrig, FixDirs, FixSpacings, FixDims = ImportImAtt(FilePrefix='Fix')
    MovOrig, MovDirs, MovSpacings, MovDims = ImportImAtt(FilePrefix='Mov')
    
    FixContours = PCStoICSnested(Pts_PCS=ContourData['FixPtsPCS'], 
                                 Origin=FixOrig, 
                                 Directions=FixDirs, 
                                 Spacings=FixSpacings)
    
    MovContours = PCStoICSnested(Pts_PCS=ContourData[MovKey], 
                                 Origin=MovOrig, 
                                 Directions=MovDirs, 
                                 Spacings=MovSpacings)
    


# Create filename for exported figure:
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())       + f'_Dicoms_and_contours__MustBeOnLineSeg_{MustBeOnLineSeg}__'       + f'MaxDistToPts_{MaxDistToPts}__' + Perspective + '.png'

PlotImagesAndContours(FixContours=FixContours, 
                      MovContours=MovContours, 
                      SlicesToPlot=SlicesToPlot, 
                      Perspective=Perspective, 
                      UsePCSforAxes=UsePCSforAxes, 
                      ContoursAs=ContoursAs, 
                      LogToConsole=LogToConsole, 
                      ExportFig=ExportResults, 
                      ExportFname=ExportFname)


# In[19]:


import ContourPlottingFuncs
importlib.reload(ContourPlottingFuncs)
from ContourPlottingFuncs import PlotContoursAndPlanes3D_v2 as Plot3D
#%matplotlib inline
#%matplotlib notebook
   
# Decide on the following parameters:

""" Parameters for plotting: """

# Key of contours to plot:
key = 'FixPtsPCS'; 
key = 'FixSSPtsPCS'
key = 'FixOSPtsPCS'
key = 'TraPtsPCS'; 
key = 'TraSSPtsPCS'
key = 'TraOSPtsPCS'
#key = 'MovSSPtsPCS'
#key = 'MovOSPtsPCS'

# Choose whether to convert points from PCS to ICS:
Convert2ICS = False
#Convert2ICS = True

#PlotImagingPlanes=False
PlotImagingPlanes=True

#PlotJoiningSegments = False
PlotJoiningSegments = True

# Choose how many segments to plot:
"""
Depending on whether plotting original list of points, or
an over-sampled or super-sampled list, it is easier to
visualise if not all line segments are displayed.
"""
EveryNthSegment = 5 # plot every 5th
EveryNthSegment = 20
#EveryNthSegment = 50

# Choose whether to use random colours for each segment:
SegmentColoursRandom = False

LogToConsole = True
LogToConsole = False

ExportResults = True
ExportResults = False

""" end of list of parameters """


# Depending on which contours are plotted determines which
# image attributes will be used to plot the planes:
if 'Fix' in key:
    # Plot the Fixed image planes:
    Origin = FixOrig
    Dirs = FixDirs
    Spacings = FixSpacings
    Dims = FixDims

if ('Tra' in key) or ('Mov' in key):
    # Plot the Moving planes when plotting either transformed or
    # moving contours:
    Origin = MovOrig
    Dirs = MovDirs
    Spacings = MovSpacings
    Dims = MovDims



Plot3D(Contours=ContourData[key], Origin=Origin, Directions=Dirs, Spacings=Spacings, 
       Dimensions=Dims, Convert2ICS=Convert2ICS, PlotImagingPlanes=PlotImagingPlanes, 
       PlotJoiningSegments=PlotJoiningSegments, EveryNthSegment=EveryNthSegment, 
       SegmentColoursRandom=SegmentColoursRandom, PlotTitle=key, LogToConsole=LogToConsole,
       ExportPlot=ExportResults)
       


# In[ ]:




