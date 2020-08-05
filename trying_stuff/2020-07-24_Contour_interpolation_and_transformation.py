#!/usr/bin/env python
# coding: utf-8

# In[62]:


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
PlotResults = False

# Decide whether to annotate points with their numbers:
AnnotatePtNums = True
AnnotatePtNums = False

# Decide whether to export plotted results to file:
ExportResults = True
ExportResults = False



PointData, FixContourData,InterpData, MovContourData = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                          dP, InterpolateAllPts, LogToConsole, 
                                          PlotResults, AnnotatePtNums, ExportResults)


"""
Timings to interpolate 5 contours using sqrt:
4.2 s
4.1 s
4.2 s

Timings to interpolate 5 contours without using sqrt:
4.0 s
4.0 s
4.2 s
"""


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


# In[ ]:


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


# In[ ]:


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


# ### Aug 3: I might have mixed up points in different coordinate systems when getting the intersection points - the points were in the PCS but the planes' normals and points are defined in the ICS.
# 
# ### I've copied GetIntersectingContoursInMovingPlanes and have a new version GetIntersectingContoursInMovingPlanes_v2 that converts the points to ICS first.
# 
# ### See if this helps...

# In[11]:


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
    
MovContourData3 = cif.GetIntersectingContoursInMovingPlanes_v2(InterpData, 
                                                               MovDicomDir,
                                                               UseInterp,
                                                               BetweenLinePts)
AnnotatePtNums=True
AnnotatePtNums=False

SubPlots = True
#SubPlots = False

ExportPlot=True
ExportPlot=False

cif.PlotIntersectingPts2D(MovContourData3, dP, AnnotatePtNums, SubPlots, ExportPlot)


# In[13]:


import pandas as pd

pd.DataFrame(data=MovContourData3)


# ### On second thought the above must be wrong since the image attributes are in the PCS!

# In[ ]:


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


# ### 04/08:  Re-order points clockwise

# In[6]:


import copy

# Work with the 17th slice in MovContourData:

i = 17

Points = copy.deepcopy(MovContourData['PointPCS'][i][0])

AnnotatePtNums = True
#AnnotatePtNums = False

ExportPlot = True
ExportPlot = False

cif.PlotIntersectingPoints2D_OLD(MovContourData, i, dP, AnnotatePtNums, ExportPlot)


# i = 17
# 
# Points = copy.deepcopy(MovContourData['PointPCS'][i][0])
# 
# 
# # Define the origin as the point with minimal y-coordinate.
# 
# # Get list of y-coords of Points:
# Points_y = [Points[i][1] for i in range(len(Points))]
# 
# # Get the index of the point with minimum y-coordinate: 
# IndMinY = Points_y.index(min(Points_y))
# 
# OriginPt = Points[IndMinY]
# 
# AllInds = list(range(len(Points)))
# #OtherInds = AllInds.remove(IndMinY)
# #OtherInds = list(range(len(Points))).remove(IndMinY)
# 
# OtherInds = []
# 
# for i in AllInds:
#     if AllInds[i] != IndMinY:
#         OtherInds.append(i)
#         
# OtherPts = [Points[i] for i in OtherInds]
# 
# import math
# 
# 
# for i in range(len(OtherPts)):
#     angle = math.atan2(OtherPts[i][1] - OriginPt[1], OtherPts[i][0] - OriginPt[0])
#     
#     if angle < 0:
#         angle += 2 * math.pi
#     
#     print(angle)

# In[103]:



#PlotPoints(Points, True)
#PlotPoints(PointsSorted, True)
#PlotPoints(PointsSortedClosed, True)


# In[95]:


# Repeat for all contours:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 
        
for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        Points = Points[0]
        
        # Get list of x- and y-coords of Points:
        Points_x = [Points[i][0] for i in range(len(Points))]
        Points_y = [Points[i][1] for i in range(len(Points))]

        Mean_x = sum(Points_x)/len(Points_x)
        Mean_y = sum(Points_y)/len(Points_y)

        OriginPt = [Mean_x, Mean_y]

        RefVector = [0, 1] # pointing up
        #RefVector = [1, 0] # pointing right
        #RefVector = [-1, 0] # pointing left

        PointsSorted = sorted(OtherPts, key=ClockwiseAngleAndDistance)
        #PointsSorted = sorted(Points, key=ClockwiseAngleAndDistance)
        #PointsSorted

        # Close the contour:
        PointsSortedClosed = cif.CloseContour(PointsSorted)
        #PointsSortedClosed
        
        PlotPoints(PointsSortedClosed, True)


# In[146]:


# Repeat for all contours:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 


AnnotatePts = True
ExportPlot = True

        
for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        Points = Points[0]
        
        # Sort Points by clockwise ordering:
        Sorted = cif.SortPointsClockwise(Points)

        # Close the contour:
        SortedClosed = cif.CloseContour(Sorted)
        
        PlotTitle = f'Sorted intersection points for slice {i} using clockwise ordering'
        
        cif.PlotPoints(SortedClosed, AnnotatePts, PlotTitle, ExportPlot)


# ### Clearly the sorting algorithm is not suitable.
# 
# ### Try something else.

# In[114]:


i = 17

AnnotatePtNums = True
#AnnotatePtNums = False

ExportPlot = True
ExportPlot = False

cif.PlotIntersectingPoints2D_OLD(MovContourData, i, dP, AnnotatePtNums, ExportPlot)


# In[112]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 
    
i = 17

Points = MovContourData['PointPCS'][i][0]

OriginPts = cif.GetCentroid(Points)
    
for i in range(1, len(Points)):
    print(cif.Point0IsClockwiseToPoint1(Points[i-1], Points[i], OriginPts))


# ### I'm not sure how to use the above function (modified from ciamej's contribution here:
# ### https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order)
# 
# ### Beside, I sense that it will give a similar result to the above results - i.e. sorting by angles.  So I'll try something else..

# In[149]:


"""
Idea for a new sorting algorithm:

1. Start with the first point in a list of points.  Append the first point to the a new list called Sorted.

2. Find the next closest point.  Append it to Sorted.

3. Repeat step 2 until there are no more points left.
"""

i = 12
i = 17

AnnotatePtNums = True
#AnnotatePtNums = False

ExportPlot = True
ExportPlot = False

cif.PlotIntersectingPoints2D_OLD(MovContourData, i, dP, AnnotatePtNums, ExportPlot)

Points = copy.deepcopy(MovContourData['PointPCS'][i][0])

# The first point:
Point0 = copy.deepcopy(Points[0])

# Remaining points:
RemPts = copy.deepcopy([Points[i] for i in range(1, len(Points))])

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 


# Initiate sorted list of Points:
Sorted = [copy.deepcopy(Point0)]

while RemPts:
    Lengths = cif.GetVectorLengths(Point0, RemPts)

    print(Lengths, '\n\n')

    # Get index of point in RemPts that is closest to Point0:
    ind = Lengths.index(min(Lengths))

    print(ind, '\n\n')

    # Re-define Point0:
    Point0 = RemPts[ind]
    
    #Sorted.append(RemPts[ind])
    Sorted.append(Point0)

    # Pop the ind^th point in RemPts:
    RemPts.pop(ind)
    
    


# In[150]:


AnnotatePts = True
ExportPlot = False
PlotTitle = f'Sorted intersection points for slice {i} using SortByClosest ordering'



cif.PlotPoints(Sorted, AnnotatePts, PlotTitle, ExportPlot)


# In[152]:


# Repeat for all contours:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 


AnnotatePts = True
ExportPlot = True

        
for i in range(len(MovContourData['PointPCS'])):
    Points = MovContourData['PointPCS'][i]
    
    if Points:
        Points = Points[0]
        
        # Sort Points by clockwise ordering:
        Sorted = cif.SortByClosest(Points)

        # Close the contour:
        SortedClosed = cif.CloseContour(Sorted)
        
        PlotTitle = f'Sorted intersection points for slice {i} using SortByClosest ordering'
        
        cif.PlotPoints(SortedClosed, AnnotatePts, PlotTitle, ExportPlot)


# In[ ]:





# def ReorderPtsClockwise():
#     
#     bool less(point a, point b)
# {
#     if (a.x - center.x >= 0 && b.x - center.x < 0)
#         return true;
#     if (a.x - center.x < 0 && b.x - center.x >= 0)
#         return false;
#     if (a.x - center.x == 0 && b.x - center.x == 0) {
#         if (a.y - center.y >= 0 || b.y - center.y >= 0)
#             return a.y > b.y;
#         return b.y > a.y;
#     }
# 
#     // compute the cross product of vectors (center -> a) x (center -> b)
#     int det = (a.x - center.x) * (b.y - center.y) - (b.x - center.x) * (a.y - center.y);
#     if (det < 0)
#         return true;
#     if (det > 0)
#         return false;
# 
#     // points a and b are on the same line from the center
#     // check which point is closer to the center
#     int d1 = (a.x - center.x) * (a.x - center.x) + (a.y - center.y) * (a.y - center.y);
#     int d2 = (b.x - center.x) * (b.x - center.x) + (b.y - center.y) * (b.y - center.y);
#     return d1 > d2;
# }
# 
