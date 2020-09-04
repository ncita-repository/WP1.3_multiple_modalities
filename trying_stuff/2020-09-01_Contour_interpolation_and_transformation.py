#!/usr/bin/env python
# coding: utf-8

# In[6]:


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

# Decide whether or not to interpolate between all contour point pairs or only 
# those for which the node of the longest list of contour points was an 
# original node (see Note 4 in the function InterpolateBetweenContours):
InterpolateAllPts = True

""" Parameters for intersecting points: """

# Decide whether or not to only consider intersection points that line on the
# line segment:
OnLineSegment = True

# Decide whether or not to only consider the intersection point that is closest to
# the imaging plane:
OnlyKeepClosestIntersPt = True
OnlyKeepClosestIntersPt = False

# Decide on the maximum distance for an intersecting point from the imaging planes:
MaxDistToImagePlane = 5
MaxDistToImagePlane = 2.5
#MaxDistToImagePlane = 0.5
#MaxDistToImagePlane = 10

""" Parameters for logging to console and plotting: """

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
                                            dP, InterpolateAllPts, OnlyKeepClosestIntersPt,
                                            MaxDistToImagePlane, LogToConsole, PlotResults, 
                                            AnnotatePtNums, ExportResults)

"""
Previously took 18 s to register and transform points row-by-row and column-by-column in 
InterpData (i.e. multiple file I/Os).

Now takes 9 s to register the images, and <1 s to transform all points.
"""


# In[180]:


import pandas as pd

pd.DataFrame(MovContourData)


# In[181]:


pd.DataFrame(InterpData)


# In[ ]:





# In[ ]:





# In[ ]:




