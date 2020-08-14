#!/usr/bin/env python
# coding: utf-8

# In[829]:


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



PointData0, FixContourData0,InterpData0, MovContourData0 = cif.CopyRois(FixDicomDir, MovDicomDir, FixRoiFpath, 
                                            dP, InterpolateAllPts, OnlyKeepClosestIntersPt,
                                            MaxDistToImagePlane, LogToConsole, PlotResults, 
                                            AnnotatePtNums, ExportResults)


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


# ### Aug 14:
# 
# ### Something is wrong when running GetClosestIntersectionPoints() including both intersection points...

# In[831]:


pd.DataFrame(InterpData0)


# In[834]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

interpdata = cif.GetIntersectingPointsInMovingPlanes_v2(InterpData0, MovDicomDir, 
                                           UseInterp, MustBeOnLineSeg,
                                           OnlyKeepClosestIntersPt)


# In[835]:


pd.DataFrame(interpdata)


# In[836]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

interpdata2 = cif.GetClosestIntersectionPoints(interpdata, MaxDistToImagePlane)


# In[ ]:





# In[811]:


pd.DataFrame(MovContourData0)


# In[644]:


pd.DataFrame(InterpData)


# In[649]:


UseInterp
MustBeOnLineSeg


# In[668]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

MovContourData0 = cif.GetIntersectingPointsInMovingPlanes_v1(InterpData=InterpData, 
                                                             MovingDicomDir=MovDicomDir,
                                                             UseInterp=True,
                                                             MustBeOnLineSeg=True)
pd.DataFrame(MovContourData0)


# In[671]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

InterpData0 = cif.GetIntersectingPointsInMovingPlanes_v2(InterpData=InterpData, 
                                                         MovingDicomDir=MovDicomDir,
                                                         UseInterp=True,
                                                         MustBeOnLineSeg=False)
pd.DataFrame(InterpData0)


# In[677]:


for i in range(len(InterpData0['AllIntPtsAllSlices'])):
    print(len(InterpData0['AllIntPtsAllSlices'][i]))


# In[682]:


for i in range(len(InterpData0['AllIntPtsAllSlices'][0])):
    print(len(InterpData0['AllIntPtsAllSlices'][0][i]))


# In[684]:


InterpData0['AllIntPtsAllSlices'][0][0]


# In[685]:


InterpData0['AllIntDtsAllSlices'][0][0]


# In[686]:


InterpData1 = cif.GetIntersectingPointsInMovingPlanes_v2(InterpData=InterpData, 
                                                         MovingDicomDir=MovDicomDir,
                                                         UseInterp=True,
                                                         MustBeOnLineSeg=True)
pd.DataFrame(InterpData1)


# In[695]:


InterpData1['AllIntPtsAllSlices'][0][0]


# In[696]:


InterpData1['AllIntDtsAllSlices'][0][0]


# In[698]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

InterpData1_2 = cif.GetClosestIntersectionPoints(InterpData=InterpData1, 
                                                 MaxDistToImagePlane=5)

pd.DataFrame(InterpData1_2)


# In[702]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

MovContourData1 = cif.GroupIntersPtsBySliceIndAndContourNum(InterpData=InterpData1_2, 
                                                            MovingDicomDir=MovDicomDir)

pd.DataFrame(MovContourData1)


# ### Which slice numbers in Fixed image contribute to the contours in MovContourData?

# In[704]:


for i in range(len(MovContourData1['InterpSliceInd'])):
    if MovContourData1['InterpSliceInd'][i]:
        MovSliceNo = MovContourData1['SliceNo'][i]
        
        InterpSliceInds = MovContourData1['InterpSliceInd'][i]
        
        print(f'Moving slice no {MovSliceNo} has points from interpolated contours at {InterpSliceInds}\n')


# In[725]:


MovContourData1['PointICS'][12]


# In[726]:


len(MovContourData1['PointICS'][12])


# In[730]:


MovContourData1['PointICS']
#len(MovContourData1['PointICS'])


# In[749]:


#MovContourData1['PointICS']
#len(MovContourData1['PointICS'])

for i in range(len(MovContourData1['PointICS'])):
    contours = MovContourData1['PointICS'][i]
    
    print(f'Slice no {i} has {len(contours)} contours')


# In[714]:


from RunSimpleElastixReg import RunSimpleElastixReg

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

RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm, LogToConsole=True)


# In[752]:


pd.DataFrame(FixContourData)


# ### FixContourData has interpolated contours.  I probably should instead keep FixContourData with only original contours and create a new dictionary called FixAndInterpContourData.
# 
# ### Create a new dictionary that only has the original contour data:

# In[754]:


FixContourData1 = {}

for Key, Values in FixContourData.items():
    OrigValues = [Values[i] for i in range(30)]
    
    FixContourData1.update({Key : OrigValues})

pd.DataFrame(FixContourData1)


# In[727]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
plot_slices = [12, 13, 14, 15]
plot_slices = [12, 13, 14, 15, 16, 17, 18]
plot_slices = [12]


# Choose the perspective:
perspective = 'axial'
#perspective = 'sagittal'
#perspective = 'coronal'

# Choose whether to export figure:
export_fig = True
export_fig = False

# Choose whether to plot contour points as lines or dots:
contours_as = 'lines'
#contours_as = 'dots'


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())             +'_Contour_interp_transformation_and_intersection_' + perspective + '.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                                 fix_pts=FixContourData1['PointICS'], 
                                                                 mov_pts=MovContourData1['PointICS'],
                                                                 export_fig=export_fig,
                                                                 export_fname=fig_fname,
                                                                 plot_slices=plot_slices,
                                                                 perspective=perspective,
                                                                 contours_as=contours_as,
                                                                 LogToConsole=True)


# In[779]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

FixPts = copy.deepcopy(FixContourData1['PointICS'])
MovPts = copy.deepcopy(MovContourData1['PointICS'])
Perspective = 'axial'
SlicesToPlot = [12]
ContoursAs = 'lines'
ExportFig = False
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())               +'_Contour_interp_transformation_and_intersection_' + perspective + '.png'
LogToConsole=True

Inds, Npas, Origins, Spacings, Orients, ARs, NumSlices, Txts, Contours, Nrows, Ncols = ruf.PrepDataForFixMovRegPlotWithContours(FixImage=FixIm, 
                                                        MovImage=MovIm, 
                                                        RegImage=RegIm, 
                                                        FixContours=FixContourData1['PointICS'], 
                                                        MovContours=MovContourData1['PointICS'],
                                                        SlicesToPlot=SlicesToPlot, 
                                                        Perspective=Perspective, 
                                                        LogToConsole=LogToConsole)


# In[805]:


#Inds
#Npas
#Origins
#Spacings
#Orients
#ARs
#NumSlices
#Txts
Contours[0]
Contours[1]
Contours[0][12]
#Contours[1][12]
#ind = 0
#ind = 1
#for i in range(len(Pts[ind])):
#    contours = Pts[ind][i] 
#    print(f'Slice no {i} has {len(contours)} contours')
#Nrows
#Ncols


# In[807]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

SlicesToPlot = [12]
SlicesToPlot = [12, 13, 14, 15, 16, 17, 18]

LogToConsole=True
LogToConsole=False

ExportFig = False
#ExportFig = True

ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())               + f'_Contour_interp_transformation_and_intersection_dP_{dP}_'               + perspective + '.png'

ruf.PlotFixMovRegImagesAndContours(FixIm, MovIm, RegIm, FixPts, MovPts,
                                   SlicesToPlot, Perspective, ContoursAs, 
                                   LogToConsole, ExportFig, ExportFname)


# In[789]:


FixPts[12]
len(FixPts[13])
np.array(FixPts[13]).shape

#MovPts[12]
#np.array(MovPts[12]).shape


# In[790]:


for i in range(len(FixPts)):
    print(np.array(FixPts[i]).shape)


# In[791]:


for i in range(len(MovPts)):
    print(np.array(MovPts[i]).shape)


# In[ ]:





# In[ ]:





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


# In[220]:


# Re-plot the interpolated contours with equal aspect ratio:

import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

FixedOrMoving = 'Fixed'
FixedOrMoving = 'Moving'

AnnotatePtNums = True
AnnotatePtNums = False

SubPlots = True
#SubPlots = False

ExportPlot = True
ExportPlot = False

if SubPlots:
    cif.PlotInterpContours2DSubPlots(InterpData=InterpData, 
                                     FixedOrMoving=FixedOrMoving,
                                     dP=dP, 
                                     AnnotatePtNums=AnnotatePtNums, 
                                     SubPlots=SubPlots,
                                     ExportPlot=ExportPlot)
else:
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

# In[221]:


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

if SubPlots:
    cif.PlotIntersectingPts2DSubPlots(MovContourData2, dP, AnnotatePtNums, SubPlots, ExportPlot)
else:
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


# ### Try the Travelling Salesperson Problem:
# 
# ### https://codereview.stackexchange.com/questions/81865/travelling-salesman-using-brute-force-and-heuristics

# In[158]:


len(Points)


# In[166]:


from itertools import permutations
import time


i = 12
#i = 17

Points = MovContourData['PointPCS'][i][0]

# Reduce Points to a shorter list points:
n = 11

points = [Points[i] for i in range(n)]

P = len(points)

print(f'There are {P} points in points.\n\n')

#print('points:\n\n', points, '\n\n')

# Start timing:
times = []
times.append(time.time())

Permutations = [perm for perm in permutations(points)]

R = len(Permutations)

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)

print(f'Took {Dtime} s to get {R} permuations of {P} points.')

#print('Permutations:\n\n')
#Permutations


"""
N = 3 took 0.0 s
N = 6 took 0.0 s
N = 10 took 1.5 s
N = 15 took .... too long!  I interrupted the kernel
"""


# ### It's clear that using bruteforce is not going to work.  It's taking impractically long to get the permutations for 15 points, so there's no chance of this approach being used for 100s let along 1000s of points.
# 
# ### The other algorithm listed on the stackexchange link called "optimised_travelling_salesman" is essentially the same as my SortByClosest function.
# 
# ### While waiting I've come up with a new idea - partially use SortByClosest.  But rather than sorting all points, only adjust the ordering if a line segment will exceed the min/mean segment length.
# 
# ### On the other hand, that will almost certainly fail for the case of slice 17 with the U-shaped points.
# 
# ### Need to come up with another strategy (05/08/2020).

# In[585]:


import pandas as pd

pd.DataFrame(InterpData)
#pd.DataFrame(FixContourData)
#pd.DataFrame(MovContourData)


# In[344]:


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
#cif.PlotInterpolatedContours3D(InterpData=InterpData,
cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
                                  FixContourData=FixContourData, 
                                  MovContourData=MovContourData,
                                  FixDicomDir=FixDicomDir,
                                  MovDicomDir=MovDicomDir,
                                  dP=dP, 
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[580]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

Convert2ICS = False
#Convert2ICS = True

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True


#cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
cif.PlotInterpolatedContours3D(InterpData=InterpData,
#cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
                                  FixContourData=FixContourData, 
                                  MovContourData=MovContourData,
                                  FixDicomDir=FixDicomDir,
                                  MovDicomDir=MovDicomDir,
                                  Convert2ICS=Convert2ICS,
                                  dP=dP, 
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[581]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

Convert2ICS = False
#Convert2ICS = True

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True


#cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
#cif.PlotInterpolatedContours3D(InterpData=InterpData,
cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
                                  FixContourData=FixContourData, 
                                  MovContourData=MovContourData,
                                  FixDicomDir=FixDicomDir,
                                  MovDicomDir=MovDicomDir,
                                  Convert2ICS=Convert2ICS,
                                  dP=dP, 
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[580]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

Convert2ICS = False
Convert2ICS = True

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True


#cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
cif.PlotInterpolatedContours3D(InterpData=InterpData,
#cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
                                  FixContourData=FixContourData, 
                                  MovContourData=MovContourData,
                                  FixDicomDir=FixDicomDir,
                                  MovDicomDir=MovDicomDir,
                                  Convert2ICS=Convert2ICS,
                                  dP=dP, 
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[590]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

Convert2ICS = False
Convert2ICS = True

PlotImagingPlanes = True
#PlotImagingPlanes = False

CombinePlots = True
CombinePlots = False

ExportPlot = False
#ExportPlot = True


#cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=False)
#cif.PlotInterpolatedContours3D(InterpData=InterpData,
cif.PlotInterpolatedContours3DNew(InterpData=InterpData, 
                                  FixContourData=FixContourData, 
                                  MovContourData=MovContourData,
                                  FixDicomDir=FixDicomDir,
                                  MovDicomDir=MovDicomDir,
                                  Convert2ICS=Convert2ICS,
                                  dP=dP, 
                                  PlotImagingPlanes=PlotImagingPlanes,
                                  CombinePlots=CombinePlots, ExportPlot=ExportPlot)


# In[338]:


#data = [(1,2,3), (10,20,30), (11, 22, 33), (110, 220, 330)]
#X,Y,Z = zip(*data)

X,Y,Z = zip(*FixContourData['PointPCS'][14][0])

X


# In[455]:


pd.DataFrame(InterpData)


# In[470]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

InterpData2 = copy.deepcopy(InterpData)

UseInterp = True
#UseInterp = False
    
#InterpData2, \
#MovContourData3, \
#AllIntersPoints, \
#AllIntersDistances, \
#IntersPoints, \
#IntersDistances = cif.GetIntersectingPointsInMovingPlanes(InterpData2, 
#                                                          MovDicomDir,
#                                                          UseInterp)

InterpData2 = cif.GetIntersectingPointsInMovingPlanes(InterpData2, 
                                                      MovDicomDir,
                                                      UseInterp)

AnnotatePtNums=True
AnnotatePtNums=False

SubPlots = True
#SubPlots = False

ExportPlot=True
ExportPlot=False

#cif.PlotIntersectingPts2D(MovContourData3, dP, AnnotatePtNums, SubPlots, ExportPlot)


# In[482]:


pd.DataFrame(InterpData2)


# In[479]:


#InterpData2['AllIntPtsAllSlices'][0]
print(len(InterpData2['MovInterpContourPts'][0]))
print(len(InterpData2['AllIntPtsAllSlices'][0]))
InterpData2['AllIntPtsAllSlices'][0][0]
print(len(InterpData2['AllIntPtsAllSlices'][0][0]))
#InterpData2['AllIntDtsAllSlices'][0]
#len(InterpData2['AllIntDtsAllSlices'][0])


# In[514]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

InterpData3 = copy.deepcopy(InterpData2)

InterpData3 = cif.GetClosestIntersectionPoints(InterpData3)


# In[515]:


pd.DataFrame(InterpData3)


# In[500]:


InterpData3['ClosestIntersDts']


# In[543]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif


#IntersPts, IntersDts = cif.GroupIntersPtsBySliceIndAndContourNum(InterpData3, MovDicomDir)
MovContourData3, CTBySliceByContour, PtsPCSBySliceByContour, PtsICSBySliceByContour, DtsBySliceByContour, ISIBySliceByContour = cif.GroupIntersPtsBySliceIndAndContourNum(InterpData3, MovDicomDir)


# In[521]:


#pd.DataFrame(IntersPts)
IntersPts
#len(IntersPts)


# In[542]:


PtsBySliceByContour
#len(PtsBySliceByContour)
PtsBySliceByContour[0]
ISIBySliceByContour

for i in range(len(ISIBySliceByContour)):
    #print(f'Slice {i}, ISI = {ISIBySliceByContour[i]}\n\n')
    #print(f'Slice {i}, Dts = {DtsBySliceByContour[i]}\n\n')
    #print(f'Slice {i}, Pts PCS = {PtsPCSBySliceByContour[i]}\n\n')
    #print(f'Slice {i}, Pts ICS = {PtsICSBySliceByContour[i]}\n\n')
    print(f'Slice {i}, CT = {CTBySliceByContour[i]}\n\n')


# In[535]:


pd.DataFrame(MovContourData)


# In[552]:


MovContourData3['PointPCS'][0][1]
#len(MovContourData3['PointPCS'][0])


# In[544]:


pd.DataFrame(MovContourData3)


# In[562]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

cif.Plot2DSubPlotsOfMovContourData(MovContourData3, dP, AnnotatePtNums=False, 
                                   SubPlots=True, ExportPlot=False)


# In[577]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

MustBeOnLineSeg = False
MustBeOnLineSeg = True

MaxDistToImagePlane = 5 # mm
#MaxDistToImagePlane = 2.5 # mm
#MaxDistToImagePlane = 0.5 # mm


InterpData2 = copy.deepcopy(InterpData)

UseInterp = True
#UseInterp = False

InterpData2 = cif.GetIntersectingPointsInMovingPlanes(InterpData2, 
                                                      MovDicomDir,
                                                      UseInterp,
                                                      MustBeOnLineSeg)

InterpData3 = copy.deepcopy(InterpData2)

InterpData3 = cif.GetClosestIntersectionPoints(InterpData3, MaxDistToImagePlane)

MovContourData3, CTBySliceByContour, PtsPCSBySliceByContour, PtsICSBySliceByContour, DtsBySliceByContour, ISIBySliceByContour = cif.GroupIntersPtsBySliceIndAndContourNum(InterpData3, MovDicomDir)

pd.DataFrame(InterpData2)
#pd.DataFrame(InterpData3)
#pd.DataFrame(MovContourData3)


# In[571]:


cif.Plot2DSubPlotsOfMovContourData(MovContourData3, dP, AnnotatePtNums=False, 
                                   SubPlots=True, ExportPlot=False)


# In[573]:


cif.GetLinePlaneIntersection(PlaneN=[0,0,1], PlaneP=[-10,0,0], 
                             LineP0=[4,3,2], LineP1=[4,3,-8], 
                             MustBeOnLineSeg=False)


# ### Aug 13 - Try something new:  Rather than finding where the transformed contours intercept the Moving image planes, create a list of a list of contour points, grouping points into different slices according to the slice that is nearest to the point.

# In[608]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

cif.GroupPointsIntoNearestPlanes(InterpData['MovOSContour1Pts'][0],
                                 DicomDir=MovDicomDir)


# ### Clearly this is taking too long.  It takes ~5 s to process one point, it will take 21 min to process the 253 points in that one contour.  Scaling this up to multiple contours let along multiple contours with a higher number of points is unrealistic.

# In[615]:


GridPtsPerPlane = cif.GetGridOfPointsInImagePlanes(DicomDir=MovDicomDir)


for i in range(len(GridPtsPerPlane)):
    Z = [GridPtsPerPlane[i][j][2] for j in range(len(GridPtsPerPlane[i]))]
    
    minZ = min(Z)
    maxZ = max(Z)
    
    print(f'minZ = {minZ}, maxZ = {maxZ}')


# ### I can't use the z-index alone to group points (as I suspected), since, for example, z = - 7 mm spans the first 5 slices.

# ### Are any points missing in the transformed columns of InterpData?

# In[618]:


for i in range(len(InterpData['InterpSliceInd'])):
    InterpSliceInd = InterpData['InterpSliceInd'][i]
    
    FixContour1Pts = InterpData['FixOSContour1Pts'][i]
    FixContour2Pts = InterpData['FixOSContour2Pts'][i]
    FixInterpContourPts = InterpData['FixInterpContourPts'][i]
    
    MovContour1Pts = InterpData['MovOSContour1Pts'][i]
    MovContour2Pts = InterpData['MovOSContour2Pts'][i]
    MovInterpContourPts = InterpData['MovInterpContourPts'][i]
    
    print(f'InterpSliceInd {InterpSliceInd} has {len(FixContour1Pts)} Fixed contour1 points')
    print(f'InterpSliceInd {InterpSliceInd} has {len(FixContour2Pts)} Fixed contour2 points')
    print(f'InterpSliceInd {InterpSliceInd} has {len(FixInterpContourPts)} Fixed interpolated contour points')
    print(f'InterpSliceInd {InterpSliceInd} has {len(MovContour1Pts)} Moving contour1 points')
    print(f'InterpSliceInd {InterpSliceInd} has {len(MovContour2Pts)} Moving contour2 points')
    print(f'InterpSliceInd {InterpSliceInd} has {len(MovInterpContourPts)} Moving interpolated contour points\n\n')


# ### So it looks like there are no missing points.
# 
# ### Are there any negative indices?

# In[639]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

for key, values in InterpData.items():
    if ( 'Contour' in key ) and ( 'Ind' in key ):
        #print(key + ':\n')
        #print(values, '\n')
        
        for i in range(len(values)):
            MinMaxInds = cif.GetMinMaxIndicesInList(values[i])
            
            InterpSliceInd = InterpData['InterpSliceInd'][i]
            
            #print(f'   InterpSliceInd {InterpSliceInd}:\n')
            print(key, f'at InterpSliceInd {InterpSliceInd}:\n')
            print(f'   Min, Max x-indices = {MinMaxInds[0:2]}')
            print(f'   Min, Max y-indices = {MinMaxInds[2:4]}')
            print(f'   Min, Max z-indices = {MinMaxInds[4:6]}')
            print(f'   Set of z-indices = {MinMaxInds[6::]}\n')
            
        


# ### There are no negative indices.

# In[591]:


pd.DataFrame(InterpData)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[403]:


pd.DataFrame(AllIntersPoints)


# In[404]:


pd.DataFrame(AllIntersDistances)


# In[389]:


[AllIntersDistances[s][0] for s in range(35)]


# In[447]:


# Initialise a list of points and distances for each imaging plane:
#IntersPoints = [[] for i in range(MovDims[2])]
#IntersDistances = [[] for i in range(MovDims[2])]
IntersPoints = [[] for i in range(35)]
IntersDistances = [[] for i in range(35)]
IntersIndices = [[] for i in range(35)]
#IntersPoints = []
#IntersDistances = []
#IntersIndices = []

S = len(AllIntersPoints) # num of slices
P = len(AllIntersPoints[0]) # num of points for any slice


for p in range(P):
    #Points = [AllIntersPoints[i][p] for i in range(S)]
    
    #Distances = [AllIntersDistances[i][p] for i in range(S)]
    
    Distances = []
    
    for s in range(S):
        if AllIntersDistances[s][p]:
            Distances.append(AllIntersDistances[s][p])
        else:
            # If the distance is [] arbitrarily append 100 to maintain the ordering:
            Distances.append(100)

    # The index of the closest point:
    ind = Distances.index(min(Distances))
    #IntersIndices.append(ind)
    #IntersIndices[ind].append(ind)

    MinD = Distances[ind]
    MinPt = AllIntersPoints[ind][p]
    
    print('ind =', ind, '\nMinD =', MinD, '\nMinPt =', MinPt)

    # Append closest point to IntersPoints and closest distance to 
    # IntersDistances:
    #if MinD >= 5:
    #    # Reject all points with a minimum distance >= 5 mm (the
    #    # sampling distance along z):
    #    IntersPoints.append([])
    #    IntersDistances.append([])
    #else:
    
    # Reject all points with a minimum distance >= 5 mm (the
    # sampling distance along z):
    if MinD < 5:
        ##IntersPoints[ind].append(AllIntersPoints[ind][p])
        ##IntersDistances[ind].append(AllIntersDistances[ind][p])
        ##IntersPoints.append(AllIntersPoints[ind][p])
        ##IntersDistances.append(AllIntersDistances[ind][p])
        #IntersPoints.append(MinPt)
        #IntersDistances.append(MinD)
        ##IntersPoints[ind][p].append(AllIntersPoints[ind][p])
        ##IntersDistances[ind][p].append(AllIntersDistances[ind][p])
        IntersPoints[ind].append(MinPt)
        IntersDistances[ind].append(MinD)
        IntersIndices[ind].append(ind)

        print(f'AllIntersPoints[{ind}][{p}] =', AllIntersPoints[ind][p])
        print(f'AllIntersDistances[{ind}][{p}] =', AllIntersDistances[ind][p], '\n\n')

        ##print(f'IntersPoints[{ind}] =', IntersPoints[ind])
        ##print(f'IntersDistances[{ind}] =', IntersDistances[ind], '\n\n')
        #print(f'IntersPoints[{p}] =', IntersPoints[p])
        #print(f'IntersDistances[{p}] =', IntersDistances[p], '\n\n')
        #print(f'IntersPoints[{ind}][{p}] =', IntersPoints[ind][p])
        #print(f'IntersDistances[{ind}][{p}] =', IntersDistances[ind][p], '\n\n')

        
print(f'IntersPoints =', IntersPoints)
print(f'IntersDistances =', IntersDistances)
print(f'IntersIndices =', IntersIndices)


# In[426]:


[AllIntersPoints[s][120] for s in range(35)]
[AllIntersDistances[s][120] for s in range(35)]


# In[449]:


IntersPoints
IntersDistances
#IntersIndices

#for p in range(len(IntersIndices)):
#    print(f'IntersIndices[{p}] = {IntersIndices[p]}, ',
#          f'\nIntersDistances[{p}] = {IntersDistances[p]}, ',
#          f'\nIntersPoints[{p}] = {IntersPoints[p]}\n\n')


# In[450]:


IntersDistances
IntersPoints


# In[451]:


pd.DataFrame(IntersPoints)


# In[452]:


pd.DataFrame(IntersDistances)


# In[453]:


pd.DataFrame(IntersIndices)


# In[ ]:




