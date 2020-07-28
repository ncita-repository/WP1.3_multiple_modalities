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


# In[152]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        cif.PlotIntersectingPoints2D(MovContourData, i, dP, AnnotatePtNums=False, ExportPlot=False)


# ### July 27: Tidy up the closed "contours".  Need to:
# 
# ### 1. Re-order the points so they follow a clock-wise direction
# 
# ### 2. Close the contours either by:
# ###    a) Interpolating if possible
# ###    b) Joining the gaps if not
# 
# ### The tricky part will be in identifying the gaps - what will constitute a gap, and how to determine whether to use interpolation or join with a straight segment?
# 
# ### Here are some ideas:
# 
# ### What consistitues a gap - Determine the min or average segment length. If any given segment length exceeds the min or ave by a certain threshold it is considered a gap.
# 
# ### How to close - if the gap length exceeds a threshold value join the nodes with a line.  If not, use interpolation.  The threshold could be either set in absolute terms or as a length relative to the "perimeter" of the "contour".

# In[47]:


import numpy as np

#for i in range(len(MovContourData['PointPCS'])):
#    Points = MovContourData['PointPCS'][i]
    

    
for s in range(len(FixContourData['PointPCS'])):
    Points = FixContourData['PointPCS'][s]
    
    if Points:
        print(f's = {s}\n\n')
        
        Points = Points[0] # assumes only one contour per slice on contour points
        
        P = len(Points)
        
        MeanX = 0
        
        for i in range(P):
            MeanX += Points[i][0] / P
            
        # Check the direction of the points:
        CheckSum = 0
        j = 1
        k = 2

        for i in range(P):
            Result = (Points[j][0] - MeanX) * (Points[k][1] - Points[i][1])
            
            print(f'Result = {Result}')
            
            CheckSum += Result

            j += 1
            k += 1

            if j >= P:
                j = 0
            if k >= P:
                k = 0


        if CheckSum > 0:
            #Contour.reverse()

            print(f'CheckSum = {CheckSum}')
            
            
        print(f'CheckSum = {CheckSum}\n\n')


# In[ ]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

AnnotatePtNums = True

for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        cif.PlotIntersectingPoints2D(MovContourData, i, dP, AnnotatePtNums, ExportPlot=True)


# ### Export contour data to csv for Matt Orton to compare with his contour interpolation algorithm:

# In[71]:


import csv

for i in range(len(FixContourData['PointPCS'])):
    Points = FixContourData['PointPCS'][i]
    
    if Points:
        Points = Points[0]
        
        SliceNo = FixContourData['SliceNo'][i]
        
        with open(f"FixContourPts_slice{SliceNo}.csv", "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(Points)
            
            #for r in range(len(Points)):
            #    writer.writerow(Points[r])
    

for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        Points = Points[0]
        
        SliceNo = MovContourData['SliceNo'][i]
        
        with open(f"MovContourPts_slice{SliceNo}.csv", "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(Points)
            


# In[72]:


# Work with one collection of points:

i = 12

Points = MovContourData['PointPCS'][i]

cif.PlotIntersectingPoints2D(MovContourData, i, dP, AnnotatePtNums, ExportResults)


# In[104]:


i = 12

Points = MovContourData['PointPCS'][i][0]

P = len(Points)

# Join every point to every other point:

PointsJoined = []

for i in range(P):
    PointsJoinedThisPoint = []
    
    for j in range(P):
        if j != i:
            PointsJoinedThisPoint.append([Points[i], Points[j]])
        
    PointsJoined.append(PointsJoinedThisPoint)
    #PointsJoined.append([PointsJoinedThisPoint])
    
    
PointsJoined


# In[105]:


print(len(Points))
print(len(PointsJoined))


# In[106]:


PointsJoined[0]


# In[107]:


len(PointsJoined[0])


# In[122]:


PointsJoined[0][0]


# In[124]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 


cif.PlotJoinedPoints2D(PointsJoined=PointsJoined, ExportPlot=False, SliceNum=12)


# In[148]:


from GetImageAttributes import GetImageAttributes

GetImageAttributes(MovDicomDir, 'pydicom')

print(cif.GetIndsOfSliceNumsOfContourType(FixContourData, 1))
print(cif.GetIndsOfSliceNumsOfContourType(MovContourData, 3))


# In[168]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True


#cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
cif.PlotInterpolatedContours3D(InterpData=InterpData, 
                               FixContourData=FixContourData, 
                               MovContourData=MovContourData,
                               FixDicomDir=FixDicomDir,
                               MovDicomDir=MovDicomDir,
                               dP=dP, 
                               PlotImagingPlanes=PlotImagingPlanes,
                               CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[134]:


InterpDataDF


# In[138]:


FixContourDataDF
MovContourDataDF


# In[ ]:




