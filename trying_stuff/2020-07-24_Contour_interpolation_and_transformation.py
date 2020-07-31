#!/usr/bin/env python
# coding: utf-8

# In[3]:


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

# Decide whether or not to only consider intersection points that line
# on the line segment:
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

# Force ExportResults to False if PlotResults is False:
if not PlotResults:
    ExportResults = False



PointData, FixContourData,InterpData, MovContourData = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                          dP, InterpolateAllPts, LogToConsole, 
                                          PlotResults, AnnotatePtNums, ExportResults)



# ### Export contour data to csv for Matt Orton to compare with his contour interpolation algorithm:

# In[71]:


if False:
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
            

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 
            
ExportContourPtsToCsv(ContourData=FixContourData, Label='Fix')
ExportContourPtsToCsv(ContourData=MovContourData, Label='Mov')


# In[215]:


- (-3//2)


# In[180]:


# Re-plot the interpolated contours with equal aspect ratio:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

for i in range(len(InterpData['InterpSliceInd'])):
    Contour1 = InterpData['FixOSContour1Pts'][i]
    Contour2 = InterpData['FixOSContour2Pts'][i]
    ContourI = InterpData['FixInterpContourPts'][i]
    Contours = [Contour1, Contour2, ContourI]
    #Contours = [Contour1[0], Contour2[0], ContourI[0]]
    
    SliceInd1 = InterpData['BoundingSliceInds'][i][0]
    SliceInd2 = InterpData['BoundingSliceInds'][i][1]
    SliceIndI = InterpData['InterpSliceInd'][i]
    
    Labels = [f'Contour at slice {SliceInd1}', 
              f'Contour at slice {SliceInd2}', 
              f'Interpolated contour at {SliceIndI}']
    Colours = ['b', 'g', 'r']
    Shifts = [0, 0, 0] # no shift
    
    AnnotatePtNums = True
    AnnotatePtNums = False
    
    ExportPlot = True
    ExportPlot = False
    
    cif.PlotInterpolatedContours2D(Contours=Contours, Labels=Labels, 
                                   Colours=Colours, Shifts=Shifts,
                                   InterpSliceInd=SliceIndI,
                                   BoundingSliceInds=[SliceInd1, SliceInd2],
                                   dP=dP, AnnotatePtNums=AnnotatePtNums, 
                                   ExportPlot=ExportPlot)


# In[ ]:


# Re-plot the interpolated contours with equal aspect ratio:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

FixedOrMoving = 'Fixed'
FixedOrMoving = 'Moving'

AnnotatePtNums = True
AnnotatePtNums = False

SubPlots = True
SubPlots = False

ExportPlot = True
#ExportPlot = False

cif.PlotInterpContours2D(InterpData=InterpData, 
                         FixedOrMoving=FixedOrMoving,
                         dP=dP, 
                         AnnotatePtNums=AnnotatePtNums, 
                         SubPlots=SubPlots,
                         ExportPlot=ExportPlot)


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


# In[262]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

AnnotatePtNums=True
#AnnotatePtNums=False

SubPlots = True
#SubPlots = False

ExportPlot=True
ExportPlot=False

#for i in range(len(MovContourData['PointPCS'])):
#    Points = MovContourData['PointPCS'][i]
#    
#    if Points:
#        cif.PlotIntersectingPoints2D(MovContourData, i, dP, 
#                                     AnnotatePtNums=AnnotatePtNums, 
#                                     ExportPlot=ExportPlot)


cif.PlotIntersectingPts2D(MovContourData, dP, AnnotatePtNums, SubPlots, ExportPlot)


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


# In[257]:


# Work with one collection of points:

i = 12

Points = MovContourData['PointPCS'][i]

AnnotatePtNums = True
#AnnotatePtNums = False

ExportPlot = True
ExportPlot = False

cif.PlotIntersectingPoints2D_OLD(MovContourData, i, dP, AnnotatePtNums, ExportPlot)


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


# In[208]:


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


# ### July 29:
# 
# ### The intersection points for slice 13 look very odd around x = 10 mm, y = 45 mm, and it might be down to the fact that points of intersection from tranlated interpolated contours inherently have more uncertainty than if I were finding the intersection between translated original contours.

# In[8]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

UseInterp = True
UseInterp = False

OnLineSegment = True
    
MovContourData2 = cif.GetIntersectingContoursInMovingPlanes(InterpData, 
                                                            MovDicomDir,
                                                            UseInterp,
                                                            OnLineSegment)

AnnotatePtNums=True
AnnotatePtNums=False

SubPlots = True
#SubPlots = False

ExportPlot=True
ExportPlot=False

#for i in range(len(MovContourData2['PointPCS'])):
#    Points = MovContourData2['PointPCS'][i]
#    
#    if Points:
#        cif.PlotIntersectingPoints2D(MovContourData2, i, dP, 
#                                     AnnotatePtNums=AnnotatePtNums, 
#                                     ExportPlot=ExportPlot)

cif.PlotIntersectingPts2D(MovContourData2, dP, AnnotatePtNums, SubPlots, ExportPlot)


# ### The results look surprisingly similar.  
# 
# ### It's still unclear is why contours are not closed.

# ### July 30:
# 
# ### Might have found and fixed a bug in the function GetIntersectingContoursInMovingPlanes.

# In[304]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

UseInterp = True
#UseInterp = False

# Decide whether to only consider intersection points that lie 
# between the two contour points used to find the intersection 
# points (based on their z-components):
BetweenLinePts = True
#BetweenLinePts = False
    
MovContourData2 = cif.GetIntersectingContoursInMovingPlanes(InterpData, 
                                                          MovDicomDir,
                                                          UseInterp,
                                                          BetweenLinePts)
AnnotatePtNums=True
AnnotatePtNums=False

SubPlots = True
#SubPlots = False

ExportPlot=True
ExportPlot=False

cif.PlotIntersectingPts2D(MovContourData2, dP, AnnotatePtNums, SubPlots, ExportPlot)


# In[297]:


MovContourData2DF = pd.DataFrame(data=MovContourData2)

MovContourData2DF


# In[305]:


from ExportDictionaryToJson import ExportDictionaryToJson
from LoadDictionaryFromJson import LoadDictionaryFromJson

# Export InterpData to a JSON file:
ExportDictionaryToJson(Dictionary=InterpData, FileName='InterpData')

# Load the dictionary from the JSON file:
test = LoadDictionaryFromJson(FileName='InterpData')

test = pd.DataFrame(data=test)

test[['InterpSliceInd', 'BoundingSliceInds', 'MovOSContour1Pts', 'MovOSContour2Pts', 'MovInterpContourPts']]


# In[ ]:




