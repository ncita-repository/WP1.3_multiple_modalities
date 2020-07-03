#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import packages and functions
import SimpleITK as sitk
print(sitk.Version())
import numpy as np
import time, os, copy
import importlib
import pydicom
#import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')
#from ipywidgets import interact, fixed
#from IPython.display import clear_output

#import GetImageAttributes
#importlib.reload(GetImageAttributes)
from GetImageAttributes import GetImageAttributes
from GetInputPoints import GetInputPoints
#import CreateInputFileForElastix
#importlib.reload(CreateInputFileForElastix)
from CreateInputFileForElastix import CreateInputFileForElastix
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
from GetOutputPoints import GetOutputPoints
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm' # tumour
#FixRoiFname = r'AIM_20200626_104631.dcm' # ventricles
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose which package to use to get the Image Plane Attributes:
#package = 'sitk'
package = 'pydicom'


# Get the Image Attributes for the images:
FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                              Package=package)

MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                              Package=package)


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

# Get some info on FixIm and MovIm:
#ruf.ShowImagesInfo(FixIm, MovIm)

# Register MovIm to FixIm:
#RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Get contour points into necessary arrays:
FixPtsPCS, FixPtsBySliceAndContourPCS,FixPtsICS, FixPtsBySliceAndContourICS,LUT = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                     Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings)


# # July 3:  Work with 2 contours - e.g. the 15th and 16th slices

# In[20]:


# Number of slices:
S = len(FixPtsBySliceAndContourPCS)

for s in range(S):
    # Number of lists of points for this slice:
    L = len(FixPtsBySliceAndContourPCS[s])
        
    if L:
        print(f'\nSlice {s} has {L} lists of points:')
        
        for l in range(L):
            P = len(FixPtsBySliceAndContourPCS[s][l])
            
            print(f'   List {l} has {P} points')


# In[12]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
plot_slices = [12, 13, 14, 15]
plot_slices = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
plot_slices = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
plot_slices = [11]
plot_slices = [12]
plot_slices = [13, 14, 15, 16, 17]


# Choose the perspective:
perspective = 'axial'
#perspective = 'sagittal'
#perspective = 'coronal'

# Choose whether to export figure:
#export_fig = True
export_fig = False

# Choose whether to plot contour points as lines or dots:
contours_as = 'lines'
#contours_as = 'dots'


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_' + perspective + '.png'


ruf.display_image_and_contours(image=FixIm, 
                               points=FixPtsBySliceAndContourICS, 
                               plot_slices=plot_slices,
                               perspective=perspective,
                               contours_as=contours_as,
                               export_fig=export_fig,
                               export_fname=fig_fname)


# ### Try interpolating between slice 14 (15/30) and slice 15 (16/30):

# In[39]:



ExportPlot = False
ExportPlot = True

Slices = [14, 15]

print(f'\nSlices = {Slices}')

# Initialise the list of contours
#Contours = []

# Note this assumes that each slice has only one contour - hence 
# indexing the 0^th contour in the s^th slice:
#Contours.append([FixPtsBySliceAndContourPCS[s][0] for s in Slices])
Contours = [FixPtsBySliceAndContourPCS[s][0] for s in Slices]

#Contours

Ncontours = [len(contour) for contour in Contours]

print('\nNcontours =', Ncontours)

#MaxNcontours = max(Ncontours)
#MaxNcontours

MaxInd = Ncontours.index(max(Ncontours))
MinInd = Ncontours.index(min(Ncontours))

print(f'\nContour {MaxInd} has the greatest number of contour points')


""" Find the nearest point pairs between the two contours """

# Initialise the indices of the points in contour MinInd that 
# are closest to each point in contour MaxInd:
NearestInds = []
    
# Loop through every point in the contour with the greatest 
# number of contour points:
for i in range(Ncontours[MaxInd]):
    point_i = Contours[MaxInd][i]
    
    # The distances between point i in contour MaxInd and every 
    # point in contour MinInd:
    distances = []
    
    
    
    # Loop through every point in the second contour:
    for j in range(Ncontours[MinInd]):
        point_j = Contours[MinInd][j]
        
        distance = ( (point_i[0] - point_j[0])**2 +  (point_i[1] - point_j[1])**2 + (point_i[2] - point_j[2])**2 )**(1/2)
        
        distances.append(distance)
        
    # The index of the nearest point in the second contour:
    NearestInds.append(distances.index(min(distances)))
    

print('\nNearestInds =', NearestInds)



""" Interpolate the points by taking the midway coordinates of each linked point """

InterpContour = []

# Loop through every point in the contour with the greatest 
# number of contour points:
for i in range(Ncontours[MaxInd]):
    point_i = Contours[MaxInd][i]
    
    point_j = Contours[MinInd][NearestInds[i]]
    
    interp_x = (point_i[0] + point_j[0])/2
    
    interp_y = (point_i[1] + point_j[1])/2
    
    interp_z = (point_i[2] + point_j[2])/2
    
    InterpContour.append([interp_x, interp_y, interp_z])
    

#print('\nInterpContour =', InterpContour)


MarkerSize = 5

# Create a figure with two subplots and the specified size:
plt.subplots(1, 1, figsize=(14, 14))

# Plot contour MaxInd:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in Contours[MaxInd]:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='blue');
plt.plot(X, Y, '.b', markersize=MarkerSize);


# Plot contour MinInd:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in Contours[MinInd]:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='green');
plt.plot(X, Y, '.g', markersize=MarkerSize);


# Plot the interpolated contour:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in InterpContour:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='red');
plt.plot(X, Y, '.r', markersize=MarkerSize);

plt.title(f'Interpolation (red) of contours from slice {Slices[MaxInd]} (blue) and {Slices[MinInd]} (green)')

# Create filename for exported figure:
FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())             + f'_Interpolation_of_contour_{Slices[MaxInd]}_and_{Slices[MinInd]}.png'

if ExportPlot:
    plt.savefig(FigFname, bbox_inches='tight')
    


# ### In order to fill missing points the interpolation has to be repeated - but instead of finding the closest points in the smaller set to each point in the larger set, find the closest points in the larger set for each point in the smaller set.
# 
# ### i.e. Two passes:  
# ### 1) Find closest points in contour 2 for every point in contour 1
# ### 2) Find closest points in contour 1 for every point in contour 2

# In[40]:


# Initialise the indices of the points in contour MaxInd that 
# are closest to each point in contour MinInd:
NearestInds = []
    
# Loop through every point in the contour with the greatest 
# number of contour points:
for i in range(Ncontours[MinInd]):
    point_i = Contours[MinInd][i]
    
    # The distances between point i in contour MinInd and every 
    # point in contour MaxInd:
    distances = []
    
    
    
    # Loop through every point in the second contour:
    for j in range(Ncontours[MaxInd]):
        point_j = Contours[MaxInd][j]
        
        distance = ( (point_i[0] - point_j[0])**2 +  (point_i[1] - point_j[1])**2 + (point_i[2] - point_j[2])**2 )**(1/2)
        
        distances.append(distance)
        
    # The index of the nearest point in the second contour:
    NearestInds.append(distances.index(min(distances)))
    

print('\nNearestInds =', NearestInds)



""" Interpolate the points by taking the midway coordinates of each linked point """

InterpContour = []

# Loop through every point in the contour with the smallest 
# number of contour points:
for i in range(Ncontours[MinInd]):
    point_i = Contours[MinInd][i]
    
    point_j = Contours[MaxInd][NearestInds[i]]
    
    interp_x = (point_i[0] + point_j[0])/2
    
    interp_y = (point_i[1] + point_j[1])/2
    
    interp_z = (point_i[2] + point_j[2])/2
    
    InterpContour.append([interp_x, interp_y, interp_z])
    

#print('\nInterpContour =', InterpContour)


MarkerSize = 5

# Create a figure with two subplots and the specified size:
plt.subplots(1, 1, figsize=(14, 14))

# Plot contour MaxInd:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in Contours[MaxInd]:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='blue');
plt.plot(X, Y, '.b', markersize=MarkerSize);


# Plot contour MinInd:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in Contours[MinInd]:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='green');
plt.plot(X, Y, '.g', markersize=MarkerSize);


# Plot the interpolated contour:

# Unpack tuple and store each x,y tuple in arrays 
# X and Y:
X = []
Y = []

for x, y, z in InterpContour:
    X.append(x)
    Y.append(y)

plt.plot(X, Y, linewidth=1, c='red');
plt.plot(X, Y, '.r', markersize=MarkerSize);

plt.title(f'Interpolation (red) of contours from slice {Slices[MinInd]} (blue) and {Slices[MaxInd]} (green)')

# Create filename for exported figure:
FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())             + f'_Interpolation_of_contour_{Slices[MinInd]}_and_{Slices[MaxInd]}.png'

if ExportPlot:
    plt.savefig(FigFname, bbox_inches='tight')


# ### Use functions:

# In[46]:


import importlib
import IndsOfNearestPoints
import InterpolatePoints
import PlotInterpolation
importlib.reload(IndsOfNearestPoints)
importlib.reload(InterpolatePoints)
importlib.reload(PlotInterpolation)

from IndsOfNearestPoints import IndsOfNearestPoints
from InterpolatePoints import InterpolatePoints
from PlotInterpolation import PlotInterpolation

ExportPlot = False
ExportPlot = True

Slices = [14, 15]

print(f'\nSlices = {Slices}')

# Initialise the list of contours
#Contours = []

# Note this assumes that each slice has only one contour - hence 
# indexing the 0^th contour in the s^th slice:
#Contours.append([FixPtsBySliceAndContourPCS[s][0] for s in Slices])
Contours = [FixPtsBySliceAndContourPCS[s][0] for s in Slices]

#Contours

Ncontours = [len(contour) for contour in Contours]

print('\nNcontours =', Ncontours)

#MaxNcontours = max(Ncontours)
#MaxNcontours

MaxInd = Ncontours.index(max(Ncontours))
MinInd = Ncontours.index(min(Ncontours))

print(f'\nContour {MaxInd} has the greatest number of contour points')


# Interpolate from contour 0 to 1:
NearestIndsContour1to2 = IndsOfNearestPoints(Points1=Contours[0], 
                                             Points2=Contours[1])

InterpContour1to2 = InterpolatePoints(Points1=Contours[0], 
                                      Points2=Contours[1], 
                                      Points2IndsNearestToPoint1=NearestIndsContour1to2)


# Interpolate from contour 1 to 0:
NearestIndsContour2to1 = IndsOfNearestPoints(Points1=Contours[1], 
                                             Points2=Contours[0])

InterpContour2to1 = InterpolatePoints(Points1=Contours[1], 
                                      Points2=Contours[0], 
                                      Points2IndsNearestToPoint1=NearestIndsContour2to1)

# Plot results:
PlotInterpolation(Points1=Contours[0], Points2=Contours[1], 
                  Interp1to2=InterpContour1to2, 
                  Interp2to1=InterpContour2to1, 
                  ExportPlot=ExportPlot)


# In[ ]:





# In[ ]:




