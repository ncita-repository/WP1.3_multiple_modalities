#!/usr/bin/env python
# coding: utf-8

# In[9]:


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

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Decide whether or not to interpolate between all contour point pairs
# or only those for which the node of the longest list of contour points
# was an original node (see Note 4 in the function 
# InterpolateBetweenContours):
InterpolateAllPts = True

# Decide whether to log some results to the console:
LogToConsole = True

# Decide whether to plot results to the console:
PlotResults = True
#PlotResults = False

# Decide whether to export plotted results to file:
ExportResults = True
ExportResults = False

# Force ExportResults to False if PlotResults is False:
if not PlotResults:
    ExportResults = False



PointData, FixContourData,InterpData, MovContourData = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                          dP, InterpolateAllPts, LogToConsole, 
                                          PlotResults, ExportResults)



# In[18]:


import pandas as pd

PointDataDF = pd.DataFrame(data=PointData)
FixContourDataDF = pd.DataFrame(data=FixContourData)
InterpDataDF = pd.DataFrame(data=InterpData)
MovContourDataDF = pd.DataFrame(data=MovContourData)


# In[19]:


PointDataDF


# In[20]:


FixContourDataDF


# In[21]:


InterpDataDF


# In[22]:


MovContourDataDF


# In[17]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        cif.PlotIntersectingPoints2D(MovContourData, i, dP, ExportResults)


# In[ ]:





# In[ ]:




