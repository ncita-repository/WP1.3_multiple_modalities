#!/usr/bin/env python
# coding: utf-8

# In[3]:


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
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)

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

# In[2]:


import importlib
import IndsOfNearestPoints
import InterpolateContours
import PlotInterpolation
importlib.reload(IndsOfNearestPoints)
importlib.reload(InterpolateContours)
importlib.reload(PlotInterpolation)

from IndsOfNearestPoints import IndsOfNearestPoints
from InterpolateContours import InterpolateContours
from PlotInterpolation import PlotInterpolation

ExportPlot = False
#ExportPlot = True

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
InterpContour1to2 = InterpolateContours(Contour1=Contours[0], 
                                        Contour2=Contours[1])


# Interpolate from contour 1 to 0:
InterpContour2to1 = InterpolateContours(Contour1=Contours[1], 
                                        Contour2=Contours[0])

# Plot results:
PlotInterpolation(Contour1=Contours[0], Contour2=Contours[1], 
                  Interp1to2=InterpContour1to2, 
                  Interp2to1=InterpContour2to1, 
                  ExportPlot=ExportPlot)


# ### July 6:  Up-sample the contours.
# 
# ### This is probably only necessary for contours with few points but since it's not computationally expensive, apply to all contours.  The reason is that the surface area of the catenoid (https://en.wikipedia.org/wiki/Minimal_surface_of_revolution) need to be minimised.  Taking the cummulative length of all lines that join the pair points between both contours will approximate the surface area - an approximation that will only work if the lines are very close.

# In[9]:


np.linspace(start=1, stop=2, num=6, endpoint=True)


# In[5]:


from ContourPerimeter import ContourPerimeter
from CummulativePerimeter import CummulativePerimeter

Slices = [14, 15]

print(f'\nSlices = {Slices}')

# Note this assumes that each slice has only one contour - hence 
# indexing the 0^th contour in the s^th slice:
#Contours.append([FixPtsBySliceAndContourPCS[s][0] for s in Slices])
Contours = [FixPtsBySliceAndContourPCS[s][0] for s in Slices]


N1 = len(Contours[0])
N2 = len(Contours[1])

#P1 = ContourPerimeter(Contour=Contours[0])
#P2 = ContourPerimeter(Contour=Contours[1])

CumP1 = CummulativePerimeter(Contour=Contours[0])
CumP2 = CummulativePerimeter(Contour=Contours[1])

# Average spacing between points:
d1 = CumP1[-1]/(N1 - 1)
d2 = CumP2[-1]/(N2 - 1)

print(f'Contour1 has {N1} points and perimeter {CumP1[-1]} with mean spacing {d1}')
print(f'Contour2 has {N2} points and perimeter {CumP2[-1]} with mean spacing {d2}')


# In[80]:


# Contour1 has the greater number of points, so re-sample Contour2 so that it has the same number
# of equidistant points:

""" Re-sample contour to N equidistant points """

def ResampleContour(Contour, N):
    # Import functions:
    from ContourPerimeter import ContourPerimeter
    
    # Compute the contour perimeter:
    P = ContourPerimeter(Contour=Contour)
    
    # The required point spacing:
    dP = P/N
    
    # Initialise NewContour and add the first point in Contour:
    #NewContour = []
    #NewContour.extend(Contour[0])
    
    # Initialise NewContour:
    NewContour = []
    
    # Loop from 1 to N for the N-1 remaining points in NewContour:
    #for i in range(1, N):
    #    # The previous (Point1) and current (Point2) points:
    #    Point1 = NewContour[i-1]
        
    
    # Loop through each segment in Contour:
    for s in range(len(Contour)):
        # The points that make up this segment:
        Point1 = Contour[s]
        Point2 = Contour[s+1]
        
        # The x1, y1 and z1, and x2, y2 and z2 for this segment:
        x1 = Point1[0]
        y1 = Point1[1]
        z1 = Point1[2]
        
        x2 = Point2[0]
        y2 = Point2[1]
        z2 = Point2[2]
        
        # The dx, dy and dz for this segment:
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        
        # The length for this segment:
        SegmentL = ( dx**2 + dy**2 + dz**2 )**(1/2)
        
        # Interpolate points along the segment by the incremental
        # factor "a" that achieves the required point spacing dP:
        a = dP/( (x1**2 + y1**2 + z1**2)**(1/2) ) - 1
        
        # Add a point to NewContour by "a" until the remaining 
        # segment length is less than dP:
        while SegmentL > dP:
            # The interpolated distances:
            x = x1 + a*dx
            y = y1 + a*dy
            z = z1 + a*dz
            
            NewContour.append([x, y, z])
            
            SegmentL = SegmentL - dP
            
        
        
        
        
        
    


# In[41]:


""" Up-sample two adjacent points """

def UpsamplePoints(Point1, Point2, N):
    # Import packages:
    import numpy as np
    
    # Create list of up-sampled arrays for x, y and z:
    X = np.linspace(start=Point1[0], stop=Point2[0], num=N)
    Y = np.linspace(start=Point1[1], stop=Point2[1], num=N)
    Z = np.linspace(start=Point1[2], stop=Point2[2], num=N)
    
    # Create a list of [x, y, z] points from each array of X, Y and Z coordinates:
    Points = []
    
    for i in range(N):
        Points.append([X[i], Y[i], Z[i]])
        
        
    #print('Original point pairs:\n\n', Point1, '\n', Point2)
    #print('\nUpsampled points:\n\n', Points)
        
    return Points
    

Points = UpsamplePoints(Point1=Contours[0][0], Point2=Contours[0][1], N=6)

#Points



""" Up-sample a contour """

def UpsampleContour(OrigContour, N):
    # Import function:
    #from UpsamplePoints import UpsamplePoints
    
    #print(f'There are {len(OrigContour)} points in the original contour')
          
    NewContour = []
    
    # Iterate over every two points in OrigContour:
    for i in range(len(OrigContour) - 1):
        #print(f'i = {i}')
        
        Point1 = OrigContour[i]
        
        Point2 = OrigContour[i+1]
        
        Points = UpsamplePoints(Point1=Point1, Point2=Point2, N=N)
        
        NewContour.extend(Points)
        
        
    #print(f'There are {len(NewContour)} points in the up-sampled contour, '
    #      f'i.e. {round(len(NewContour)/len(OrigContour), 1)} more points.')
    
    return NewContour
    
    
NewContour = UpsampleContour(OrigContour=Contours[0], N=10)



""" 
Estimate the area of the surface created by linking points from Contour1 to
points in Contour2, starting with the first point in Contour1 and a point in 
Contour2 given by the index StartInd2, and adding up all line lengths.

Note: There are N points in Contour1 and Contour2.
"""

def IntegrateLineLengths(Contour1, Contour2, StartInd2):
    
    # Initialise the cummulative line lengths:
    LineLength = 0
    
    for i in range(len(Contour1)):
        Point1 = Contour1[i]
        
        for j in range(len(Contour2)):
            # The index of the point for Contour2 starts
            # at StartInd2 and increments to N, then is 
            # "wrapped" back to 0 and so on until N 
            # iterations:
            if StartInd2 + j < len(Contour2):
                Ind2 = StartInd2 + j
            else:
                Ind2 = len(Contour2) - (StartInd2 + j)
                
            Point2 = Contour2[Ind2]
            
            # The vector length between Point1 and Point2:
            VectorL = ((Point1[0] - Point2[0])**2                        + (Point1[1] - Point2[1])**2                        + (Point1[2] - Point2[2])**2)**(1/2)
            
            # Add the line length between Point1 and Point2:
            LineLength = LineLength + VectorL
            
            
    return LineLength


#for i in range(5):
#    print(AddLineLengths(Contour1=Contours[0], Contour2=Contours[1], StartInd2=i))
    

""" 
Using the function AddLineLengths(), loop through each N possible 
combinations of starting index StartInd2, and work out which
starting index yields the minimum cummulative line length.
"""

def FindIndForMinCumLength(Contour1, Contour2):
    # Import function:
    from AddLineLengths import AddLineLengths
    
    # Create list to store N cummulative line lengths for all
    # N combinations of starting index for Contour2:
    CumLineLengths = []

    # Loop through all N possible starting indeces for Contour2:
    for i in range(len(Contour2)):
        CumLineLengths.append(AddLineLengths(Contour1=Contour1, Contour2=Contour2, StartInd2=i))

    
    # Return the index that yields the minimum cummulative line length:
    return CumLineLengths.index(min(CumLineLengths))


MinStartInd2 = FindIndForMinCumLength(Contour1=Contours[0], Contour2=Contours[1])

MinStartInd2



""" 
Defining Contour1 and OrigContour2 as the two contours, re-order the 
points in OrigContour that yielded the minimum cummulative line length.
Then the points in the Contour1 and NewContour2, when joined point-by-point, 
will create the minimal surface area (i.e. catenoid) of the two contours:

https://en.wikipedia.org/wiki/Minimal_surface_of_revolution
"""


def ReorderContour(OrigContour, StartInd):
    
    # Initialise NewContour:
    NewContour = []
    
    # Re-order OrigContour by starting at the index StartInd
    # then iterating through all N subsequent points, wrapping
    # the indeces past N:
    for i in range(len(OrigContour)):
        # The index of the point for OrigContour starts
        # at StartInd and increments to N, then is 
        # "wrapped" back to 0 and so on until N 
        # iterations:
        if StartInd + i < len(OrigContour):
            Ind = StartInd + i
        else:
            Ind = len(OrigContour) - (StartInd + i)

        NewContour.append(OrigContour[Ind])
        
    return NewContour


NewContour2 = ReorderContour(OrigContour=Contours[1], 
                             StartInd=FindIndForMinCumLength(Contour1=Contours[0], Contour2=Contours[1]))

NewContour2


# In[73]:


""" 
Create list of contour points containing N points from up-sampled Contour1
(USContour1) followed by N points from the modified up-sampled Contour2 
(NewUPContour2), where:

- USContour1 is up-sampled Contour1
- USContour2 is up-sampled Contour2
- NewUSContour2 is re-ordered points from USContour2 that minimises surface area 
formed by joining points in USContour1 and USContour2
"""

def ListOfPointsOfContours(Contour1, Contour2, N):
    # Import package:
    import time
    
    # Import functions:
    from UpsampleContour import UpsampleContour
    from FindIndForMinCumLength import FindIndForMinCumLength
    from ReorderContour import ReorderContour

    # Start timing:
    times = []
    times.append(time.time())
    
    if N > 1:
        print(f'\nUp-sampling contours by {N}x...')
        
        # Up-sample Contour1 and Contour2:
        Contour1 = UpsampleContour(Contour=Contour1, N=N)
        Contour2 = UpsampleContour(Contour=Contour2, N=N)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'Took {Dtime} s to up-sample contours.')
    
    
    print(f'\nMinimising area...')
    
    # Re-order Contour2:
    ROContour2 = ReorderContour(Contour=Contour2, 
                                StartInd=FindIndForMinCumLength(Contour1=Contour1,
                                                                Contour2=Contour2))
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'Took {Dtime} s to minimise area and re-order Contour2',
          f'({round(Dtime/len(Contour2), 3)} s per point).')
    
        
    # Initialise list of points, then extend USContour1 and ROUSContour2:
    Points = []
    Points.extend(Contour1)
    Points.extend(ROContour2)
    
    print(f'\nThere are {len(Contour1)} points in Contour1, {len(Contour2)} in Contour2',
          f'and {len(Points)} in the combined list of points.')
    
    return Points
    
    
N = 1

Points = ListOfPointsOfContours(Contour1=Contours[0], Contour2=Contours[1], N=N)

"""
Timings:

N    No of points   Total time [s]    Time per point [s]

1        265              2.4               0.022
2        525             32.7 (x13.6)       0.149 (x6.8)
3        789            125.6 (x52.3)       0.381 (x17.3)
4      1,060            291.5 (x121.5)      0.662 (x30.2)

"""


# ### Use new functions:

# In[74]:


from ParseTransformixOutput import ParseTransformixOutput
from ListOfPointsOfContours import ListOfPointsOfContours


""" Deform the points in Points """

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
CreateInputFileForElastix(Points=Points)

# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
PtNos, InInds, InPts_PCS, FixOutInds,OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

#PtNos, InInds, FixOutInds,\
#MovPtsPCS, MovPtsICS,\
#DefPCS, DefICS, MovOutInds,\
#MovPtsBySlicePCS, MovPtsBySliceICS = GetOutputPoints(InPts_ICS=FixPtsICS, Origin=MovOrigin, 
#                                                     Directions=MovDirs, Spacings=MovSpacings, 
#                                                     Dimensions=MovDims, LUT=LUT)


# In[75]:


len(OutPts_PCS)


# In[60]:


list(set(LUTDict['InSliceNo']))


# In[68]:


import math

# See what the contour points FixPtsBySliceAndContourPCS look like:

ExportPlot = False


# Number of points:
P = len(LUT)

SliceIndsWithContours = list(set(LUTDict['InSliceNo']))
NumSlicesWithContours = len(SliceIndsWithContours)

# Choose colours for different contours:
colours = ['blue', 'green', 'yellow', 'red']

# Slice no:
#s = 13





MarkerSize = 5

# Nrows and Ncols in figure:
Nrows = math.ceil(math.sqrt(NumSlicesWithContours))
Ncols = NumSlicesWithContours - Nrows

# Sub-plot number:
n = 0

# Initialise the figure:
plt.subplots(Nrows, Ncols, figsize=(Ncols*5, Nrows*5))

# Iterate through all slice indeces containing contours:
for i in SliceIndsWithContours:
    n += 1
    
    # Initialise this sub-plot:
    plt.subplot(Nrows, Ncols, n)
    plt.title(f'Slice no {i}')
    
    C = len(ContourData['InPointsPCS'][i])

    print(f'There are {C} contours in the {i}^th slice of ContourData.')
    
    # Iterate through all contours for this slice:
    for c in range(C):

        # Unpack tuple and store each x,y tuple in arrays 
        # X and Y:
        X = []
        Y = []

        for x, y, z in ContourData['InPointsPCS'][i][c]:
            X.append(x)
            Y.append(y)

        plt.plot(X, Y, linewidth=1, c=colours[c]);
        plt.plot(X, Y, 'b', markersize=MarkerSize, c=colours[c]);




# Create filename for exported figure:
FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())             + f'_.png'

if ExportPlot:
    plt.savefig(FigFname, bbox_inches='tight')


# ### July 8:  Work through Matt Orton's code and implement functions:

# In[18]:


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

import GetInputPoints
importlib.reload(GetInputPoints)

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
FixPtsPCS, FixPtsBySliceAndContourPCS,FixPtsICS, FixPtsBySliceAndContourICS,LUT, PointData, ContourData = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                                             Origin=FixOrigin, Directions=FixDirs,
                                             Spacings=FixSpacings)


# In[112]:


import pandas as pd

FixPtsPCS
FixPtsBySliceAndContourPCS
LUT
PointData
ContourData


# In[113]:


# Convert the dictionaries into data frames:
PointDataDF = pd.DataFrame(data=PointData)
ContourDataDF = pd.DataFrame(data=ContourData)

PointDataDF


# In[114]:


ContourDataDF


# In[279]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif    


# Define the slice to be interpolated:
InterpSliceInd = 14

print(f'Interpolating slice {InterpSliceInd}...')

# Get slice indices that contain contours:
SliceInds = PointData['InSliceNo']

# Get extent of region of contours by slice number 
# (equivalent to function _getExtentOfRegion in 
# generateInterpolationData.js):
#Extent = [SliceInds[0], SliceInds[-1]]
#
#print('Extent =', Extent)

# Get bounding slice indices:
BoundingSliceInds = cif.GetBoundingSliceInds(InterpInd=InterpSliceInd, 
                                             SliceInds=SliceInds, 
                                             ContourData=ContourData['InPointPCS'])

print('BoundingSliceInds =', BoundingSliceInds)

#print(f'Num of points in contour on slices {BoundingSliceInds[0]} and',
#      f'{BoundingSliceInds[1]} is {len(ContourData["InPointPCS"][BoundingSliceInds[0]][0])}',
#      f'and {len(ContourData["InPointPCS"][BoundingSliceInds[1]][0])}')

# Define the contours:
"""
Note: 

Need to add [0] to the end of ContourData['InPointPCS'][BoundingSliceInds[0]] since
the contour points are a list within the list for each slice.
"""
Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])

# Close the contours:
Contour1.append(Contour1[0])
Contour2.append(Contour2[0])

print(f'len(Contour1) (after closing) = {len(Contour1)}')
print(f'len(Contour2) (after closing) = {len(Contour2)}')

# Get the segment lengths (surplus to interpolation algorithm):
SegmentLengths1 = cif.GetSegmentLengths(Contour1)
SegmentLengths2 = cif.GetSegmentLengths(Contour2)

print(f'The min, mean and max segment lengths of Contour1 are {min(SegmentLengths1)},',
      f'{sum(SegmentLengths1)/len(SegmentLengths1)}, and {max(SegmentLengths1)}.')
print(f'The min, mean and max segment lengths of Contour2 are {min(SegmentLengths2)},',
      f'{sum(SegmentLengths2)/len(SegmentLengths2)}, and {max(SegmentLengths2)}.')

# Ensure that the points in both contours run clockwise:
Contour1 = cif.ReverseIfAntiClockwise(Contour=Contour1)
Contour2 = cif.ReverseIfAntiClockwise(Contour=Contour2)

# Get the cumulative perimeters of the two bounding contours:
CumPerim1 = cif.GetCumulativePerimeters(Contour=Contour1)
CumPerim2 = cif.GetCumulativePerimeters(Contour=Contour2)

print(f'len(CumPerim1) = {len(CumPerim1)}')
print(f'len(CumPerim2) = {len(CumPerim2)}')
#print('')
print(f'CumPerim1[-1] = {CumPerim1[-1]}, CumPerim2[-1] = {CumPerim2[-1]}')
#print('')
#print('CumPerim1 =\n\n', CumPerim1)
#print('')
#print('CumPerim2 =\n\n', CumPerim2)
#print('')

# Normalise the cumulative perimeters:
CumPerim1Norm = cif.NormaliseCumulativePerimeters(CumPerimeters=CumPerim1)
CumPerim2Norm = cif.NormaliseCumulativePerimeters(CumPerimeters=CumPerim2)

print(f'len(CumPerim1Norm) = {len(CumPerim1Norm)}')
print(f'len(CumPerim2Norm) = {len(CumPerim2Norm)}')
#print('')
#print('CumPerim1Norm =\n\n', CumPerim1Norm)
#print('')
#print('CumPerim2Norm =\n\n', CumPerim2Norm)
#print('')

# Select the desired minimum super-sampled inter-node spacing:
dP = 0.2
#dP = 0.5
#dP = 1
#dP = 1.5
#dP = 2
#dP = 50
#dP = 100
dP = 0.1

print(f'The chosen minimum inter-node spacing is {dP}.')

# The number of nodes required for Contour1 and Contour2 to
# achieve the minimum inter-node spacing of dP:
""" 
Note:
    a//b yields floor(a/b)
    - (-a//b) yields ceil(a/b)
"""
MinNumNodes1 = round(- (- CumPerim1[-1] // dP))
MinNumNodes2 = round(- (- CumPerim2[-1] // dP))

# The minimum number of nodes that will need to be added:
MinNumNodes = max(MinNumNodes1, MinNumNodes2)

print(f'MinNumNodes1 = {MinNumNodes1}')
print(f'MinNumNodes2 = {MinNumNodes2}')
print(f'MinNumNodes = {MinNumNodes}')

# The number of nodes that will have to be added to Contour1 and Contour2 
# are MinNumNodes plus the number of original nodes in Contour2 and Contour1
# respectively (i.e. the nodes not common to either)?:
NumNodesToAdd1 = MinNumNodes + len(Contour2)
NumNodesToAdd2 = MinNumNodes + len(Contour1)

print(f'Num of nodes to add to Contour1 is {NumNodesToAdd1}')
print(f'Num of nodes to add to Contour2 is {NumNodesToAdd2}')
print(f'The expected length of CumSSPerim1 is {NumNodesToAdd1 - 2} (NumNodesToAdd1 - 2)',
      f'+ {len(Contour1)} (len(Contour1)) = {NumNodesToAdd1 - 2 + len(Contour1)}')
print(f'The expected length of CumSSPerim2 is {NumNodesToAdd2 - 2} (NumNodesToAdd2 - 2)',
      f'+ {len(Contour2)} (len(Contour2)) = {NumNodesToAdd2 - 2 + len(Contour2)}')

# Get the cumulative super-sampled perimeters:
CumSSPerim1 = cif.GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm=CumPerim1Norm,
                                                          NumNodesToAdd=NumNodesToAdd1)
CumSSPerim2 = cif.GetNormCumulativeSuperSampledPerimeters(CumPerimetersNorm=CumPerim2Norm,
                                                          NumNodesToAdd=NumNodesToAdd2)

#print('')
#print(f'CumSSPerim1[-1] = {CumSSPerim1[-1]}, CumSSPerim2[-1] = {CumSSPerim2[-1]}')
print(f'len(CumSSPerim1) = {len(CumSSPerim1)}')
print(f'len(CumSSPerim2) = {len(CumSSPerim2)}')

# Get list of booleans for original (True) and added (False) nodes:
IsOrigNode1 = cif.GetIsOrigNode(Contour=Contour1, NumNodesToAdd=NumNodesToAdd1)
IsOrigNode2 = cif.GetIsOrigNode(Contour=Contour2, NumNodesToAdd=NumNodesToAdd2)

#print('')
print(f'len(IsOrigNode1) = {len(IsOrigNode1)}')
print(f'len(IsOrigNode2) = {len(IsOrigNode2)}')

# Get the nodes per segment:
#NodesToAddPerSegment1 = cif.GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim1, 
#                                                    IsOrigNode=IsOrigNode1)
#NodesToAddPerSegment2 = cif.GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim2,
#                                                    IsOrigNode=IsOrigNode2)

Inds1, IsOrigNodeSorted1, IndsOrigNodes1, NodesToAddPerSegment1= cif.GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim1, 
                              IsOrigNode=IsOrigNode1)
Inds2, IsOrigNodeSorted2, IndsOrigNodes2, NodesToAddPerSegment2= cif.GetNodesToAddPerSegment(CumSuperSampledPerim=CumSSPerim2, 
                              IsOrigNode=IsOrigNode2)

#print('')
print(f'len(NodesToAddPerSegment1) = {len(NodesToAddPerSegment1)}')
print(f'len(NodesToAddPerSegment2) = {len(NodesToAddPerSegment2)}')
#print('')
print(f'sum(NodesToAddPerSegment1) = {sum(NodesToAddPerSegment1)}')
print(f'sum(NodesToAddPerSegment2) = {sum(NodesToAddPerSegment2)}')
#print('')
#print(f'Inds1 =', Inds1)
#print('')
#print(f'IsOrigNodeSorted1 =', IsOrigNodeSorted1)
#print('')
#print(f'IndsOrigNodes1 =', IndsOrigNodes1)
#print('')
#print(f'NodesToAddPerSegment1 =', NodesToAddPerSegment1)
#print('')
#print('')
#print(f'Inds2 =', Inds2)
#print('')
#print(f'IsOrigNodeSorted2 =', IsOrigNodeSorted2)
#print('')
#print(f'IndsOrigNodes2 =', IndsOrigNodes2)
#print('')
#print(f'NodesToAddPerSegment2 =', NodesToAddPerSegment2)
#print('')

# Super-sample the contours:
SSContour1, SSIsOrigNode1 = cif.SuperSampleContour(Contour=Contour1, 
                                                   NodesToAddPerSegment=NodesToAddPerSegment1)
SSContour2, SSIsOrigNode2 = cif.SuperSampleContour(Contour=Contour2, 
                                                   NodesToAddPerSegment=NodesToAddPerSegment2)

print('')
print(f'len(SSContour1) = {len(SSContour1)}')
print(f'len(SSContour2) = {len(SSContour2)}')
print('')
print(f'len(SSIsOrigNode1) = {len(SSIsOrigNode1)}')
print(f'len(SSIsOrigNode2) = {len(SSIsOrigNode2)}')


# Start timing:
times = []
times.append(time.time())

# Find the starting index for SSContour2 that minimises the integrated line lengths
# (i.e. minimal area of the surfaces created by linking points in SSContour1 and
# SSContour2 N different ways for all N points in SSContour2):
print('\nSearching for the index for SSContour2 that minimises the surface area...')
IndMinCumL = cif.FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)

times.append(time.time())
Dtime = times[-1] - times[-2]
#Dtime = round(times[-1] - times[-2], 1)
print(f'Took {round(Dtime, 2)} s to find the index for SSContour2 that minimises the surface',
      f'area ({round(Dtime/len(SSContour1), 2)} s per iteration over {len(SSContour1)}',
      'possible shifts)') 
# 0.06 s for 230 possible shifts (0.09 s/iteration; dP = 50)
# 0.1 s for 290 possible shifts (0.11 s/iteration; dP = 2)
# 0.13 s for 352 possible shifts (0.18 s/iteration; dP = 1)
# 0.25 s for 476 possible shifts (0.18 s/iteration; dP = 0.5)
# 2.26 s for 849 possible shifts (0.18 s/iteration; dP = 0.2)
# 3.09 s for 1470 possible shifts (0.18 s/iteration; dP = 0.1)

print('')
print(f'IndMinCumL = {IndMinCumL}')

# Shift the points in SSContour2 (and also SSIsOrigNode2):
SSContour2Shifted = cif.ShiftContourPoints(Contour=SSContour2, Shift=IndMinCumL)
SSIsOrigNode2Shifted = cif.ShiftContourPoints(Contour=SSIsOrigNode2, Shift=IndMinCumL)

print('')
print(f'len(SSContour1)        = {len(SSContour1)}')
print(f'len(SSContour2)        = {len(SSContour2)}')
print(f'len(SSContour2Shifted) = {len(SSContour2Shifted)}')
print('')
print(f'len(SSIsOrigNode1)        = {len(SSIsOrigNode1)}')
print(f'len(SSIsOrigNode2)        = {len(SSIsOrigNode2)}')
print(f'len(SSIsOrigNode2Shifted) = {len(SSIsOrigNode2Shifted)}')

# Remove some added nodes:
OSContour1, OSContour2, OSIsOrigNode1, OSIsOrigNode2 = cif.ReduceNodesOfContours(SuperSampledContour1=SSContour1,
                                          SuperSampledContour2=SSContour2Shifted,
                                          SuperSampledIsOrigNode1=SSIsOrigNode1, 
                                          SuperSampledIsOrigNode2=SSIsOrigNode2Shifted)

print('')
print(f'len(OSContour1) = {len(OSContour1)}')
print(f'len(OSContour2) = {len(OSContour2)}')
print('')
print(f'len(OSIsOrigNode1)        = {len(OSIsOrigNode1)}')
print(f'len(OSIsOrigNode2)        = {len(OSIsOrigNode2)}')

# Define Contour1HasMorePts:
if len(Contour1) > len(Contour2):
    Contour1HasMorePts = True
else:
    Contour1HasMorePts = False


""" 
Note:  What if len(Contour1) = len(Contour2)?...
"""

# Define the fractional interpolated slice number:
FracInterpInd = (InterpSliceInd - BoundingSliceInds[0]) / (BoundingSliceInds[1] - BoundingSliceInds[0])

print(f'FracInterpInd = {FracInterpInd}')

# Linearly interpolate between the two over-sampled contours:
InterpContour = cif.InterpolateBetweenContours(OverSampledContour1=OSContour1,
                                               OverSampledContour2=OSContour2,
                                               OverSampledIsOrigNode1=OSIsOrigNode1,
                                               OverSampledIsOrigNode2=OSIsOrigNode2,
                                               FracInterpSliceInd=FracInterpInd,
                                               Contour1HasMorePts=Contour1HasMorePts)

print('')
print(f'len(Contour1)      = {len(Contour1)}')
print(f'len(Contour2)      = {len(Contour2)}')
print(f'len(InterpContour) = {len(InterpContour)}')


# Deviate from Matt's code - try interpolating using the super-sampled contours
# (i.e. before reducing the num of nodes):
InterpSSContour = cif.InterpolateBetweenContours(OverSampledContour1=SSContour1,
                                               OverSampledContour2=SSContour2,
                                               OverSampledIsOrigNode1=SSIsOrigNode1,
                                               OverSampledIsOrigNode2=SSIsOrigNode2,
                                               FracInterpSliceInd=FracInterpInd,
                                               Contour1HasMorePts=Contour1HasMorePts)


# In[229]:


#SSContour1
#print(NodesPerSegment1)
#print(NodesPerSegment2)

#print('CumSSPerim1 added:\n\n', [CumSSPerim1[i] for i in range(NumNodesToAdd1)])
#print('\n\nCumSSPerim1 original:\n\n', [CumSSPerim1[i] for i in range(NumNodesToAdd1, len(CumSSPerim1))])

#Contour1Shifted = cif.ShiftContourPoints(Contour=Contour1, Shift=3)

#print(len(Contour1))
#print(len(Contour1Shifted))

CumLineLengths = []

# Loop through all N possible starting indeces for Contour2:
for i in range(len(Contour2)):
    CumLineLengths.append(cif.IntegrateLineLengths(Contour1=Contour1, 
                                                   Contour2=Contour2, 
                                                   StartInd2=i))


# In[231]:


CumLineLengths.index(min(CumLineLengths))


# In[283]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])

# The actual contour that exists at slice InterpSliceInd:
ActualContour = copy.deepcopy(ContourData['InPointPCS'][InterpSliceInd][0])


ExportPlot = False
#ExportPlot = True

#Contour1Shifted = cif.ShiftContourPoints(Contour=Contour1, Shift=3)


#shift = cif.FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)
#print(f'shift = {shift}')
#SSContour2Shifted = cif.ShiftContourPoints(Contour=SSContour2, Shift=shift)


# Define contours to plot and their line colours:
#Contours = [Contour1, Contour1Shifted]
#Contours = [Contour1, SSContour1, OSContour1]
#Contours = [Contour1]
#Contours = [SSContour1]
#Contours = [OSContour1]
#Contours = [Contour2, SSContour2, OSContour2]
#Contours = [SSContour1, SSContour2, SSContour2Shifted]
#Contours = [Contour1, Contour2, InterpContour]

Contours = [Contour1, Contour2, InterpContour, ActualContour]
Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
          f'Contour at slice {BoundingSliceInds[1]}', 
          f'Interpolated contour at slice {InterpSliceInd}', 
          f'Actual contour at slice {InterpSliceInd}']

Contours = [Contour1, Contour2, InterpContour, InterpSSContour]
Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
          f'Contour at slice {BoundingSliceInds[1]}', 
          f'Interpolated contour at slice {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}']

Contours = [Contour1, Contour2, InterpContour, InterpSSContour, ActualContour]
Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
          f'Contour at slice {BoundingSliceInds[1]}', 
          f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Actual contour at {InterpSliceInd}']


Colours = ['b', 'g', 'r', 'm', 'y']
Shifts = [0, 2, 4, 6, 8] # to help visualise
Shifts = [0, 0, 0, 0, 0] # no shift

MarkerSize = 5

# Create a figure with two subplots and the specified size:
fig, ax = plt.subplots(1, 1, figsize=(14, 14))

for i in range(len(Contours)):
    # Unpack tuple and store each x,y tuple in arrays 
    # X and Y:
    X = []
    Y = []

    for x, y, z in Contours[i]:
        X.append(x)
        #X.append(x + Shifts[i])
        Y.append(y + Shifts[i])

    # Define linestyle:
    if 'nterp' in Labels[i]: # i.e. interpolated contour
        if 'uper' in Labels[i]: # super-sampled
            LineStyle='dashed'
        else: # reduced nodes
            LineStyle=(0, (5, 10)) # loosely dashed   
    else:
        LineStyle='solid'
        
    # Plot line:
    ax.plot(X, Y, linestyle=LineStyle, linewidth=1, c=Colours[i], label=Labels[i]);    
    legend = ax.legend(loc='upper left', fontsize='large')
    # Plot dots:
    if False:
        plt.plot(X, Y, '.', markersize=MarkerSize, c=Colours[i]);
    
    # Plot the first and last points with different markers to help
    # identify them:
    #plt.plot(X[0], Y[0], '>', markersize=MarkerSize + 5, c=Colours[i]);
    #plt.plot(X[-1], Y[-1], 's', markersize=MarkerSize + 5, c=Colours[i]);
    
    #P = list(range(len(X))) # list of point numbers
    #plt.plot(X, Y, [f'${p}$' for p in P], markersize=MarkerSize, c=Colours[i]);
    #plt.text(X, Y, [f'{p}' for p in P], fontsize=MarkerSize+5, c=Colours[i]);
    #plt.text(X, Y, [p for p in P], fontsize=MarkerSize+5, c=Colours[i]);
    
    # Annotate every 5th point with the point number:
    if False:
        for p in range(0, len(X), 5):
            plt.text(X[p], Y[p], p, fontsize=MarkerSize+5, c=Colours[i]);
        
    # Also annotate the last point:
    #plt.text(X[-1], Y[-1], len(X) - 1, fontsize=MarkerSize+5, c=Colours[i]);
    

plt.title(f'Interpolation of contour at slice {InterpSliceInd} between '           + f'slices {BoundingSliceInds[0]} and {BoundingSliceInds[1]} '           + f'(dP = {dP})')

# Create filename for exported figure:
FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())             + f'_Interpolation_at_slice_{InterpSliceInd}_of_contours_at_'             + f'{BoundingSliceInds[0]}_and_{BoundingSliceInds[1]}_dP_{dP}.png'

if ExportPlot:
    plt.savefig(FigFname, bbox_inches='tight')


# In[287]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])

# The actual contour that exists at slice InterpSliceInd:
ActualContour = copy.deepcopy(ContourData['InPointPCS'][InterpSliceInd][0])


ExportPlot = False
#ExportPlot = True

#Contour1Shifted = cif.ShiftContourPoints(Contour=Contour1, Shift=3)


#shift = cif.FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)
#print(f'shift = {shift}')
#SSContour2Shifted = cif.ShiftContourPoints(Contour=SSContour2, Shift=shift)


# Define contours to plot and their line colours:
#Contours = [Contour1, Contour1Shifted]
#Contours = [Contour1, SSContour1, OSContour1]
#Contours = [Contour1]
#Contours = [SSContour1]
#Contours = [OSContour1]
#Contours = [Contour2, SSContour2, OSContour2]
#Contours = [SSContour1, SSContour2, SSContour2Shifted]
#Contours = [Contour1, Contour2, InterpContour]

#Contours = [Contour1, Contour2, InterpContour, ActualContour]
#Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
#          f'Contour at slice {BoundingSliceInds[1]}', 
#          f'Interpolated contour at slice {InterpSliceInd}', 
#          f'Actual contour at slice {InterpSliceInd}']

#Contours = [Contour1, Contour2, InterpContour, InterpSSContour]
#Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
#          f'Contour at slice {BoundingSliceInds[1]}', 
#          f'Interpolated contour at slice {InterpSliceInd}', 
#          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}']

Contours = [Contour1, Contour2, InterpContour, InterpSSContour, ActualContour]
Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
          f'Contour at slice {BoundingSliceInds[1]}', 
          f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Actual contour at {InterpSliceInd}']
Colours = ['b', 'g', 'r', 'm', 'y']

Contours = [InterpContour, InterpSSContour, ActualContour]
Labels = [f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Actual contour at {InterpSliceInd}']
Colours = ['r', 'm', 'y']


#Shifts = [0, 2, 4, 6, 8] # to help visualise
#Shifts = [0, 0, 0, 0, 0] # no shift


cif.PlotInterpolationResults(Contours=Contours, Labels=Labels, 
                             Colours=Colours,
                             InterpSliceInd=InterpSliceInd,
                             BoundingSliceInds=BoundingSliceInds,
                             dP=dP, ExportPlot=ExportPlot)


# In[295]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

# Define the slice to be interpolated:
InterpSliceInd = 14

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2

InterpContour, InterpSSContour, ReducedInterpSSContour = cif.InterpolateContours(InterpSliceInd, PointData, ContourData, dP)

print(f'len(InterpContour)          = {len(InterpContour)}')
print(f'len(InterpSSContour)        = {len(InterpSSContour)}')
print(f'len(ReducedInterpSSContour) = {len(ReducedInterpSSContour)}')

Contour1 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[0]][0])
Contour2 = copy.deepcopy(ContourData['InPointPCS'][BoundingSliceInds[1]][0])

# The actual contour that exists at slice InterpSliceInd:
ActualContour = copy.deepcopy(ContourData['InPointPCS'][InterpSliceInd][0])


ExportPlot = False
#ExportPlot = True

#Contour1Shifted = cif.ShiftContourPoints(Contour=Contour1, Shift=3)


#shift = cif.FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)
#print(f'shift = {shift}')
#SSContour2Shifted = cif.ShiftContourPoints(Contour=SSContour2, Shift=shift)


# Define contours to plot, their labels and their line colours:
Contours = [Contour1, Contour2, InterpContour, InterpSSContour, ActualContour]
Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
          f'Contour at slice {BoundingSliceInds[1]}', 
          f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Actual contour at {InterpSliceInd}']
Colours = ['b', 'g', 'r', 'm', 'y']
#Shifts = [0, 2, 4, 6, 8] # to help visualise
Shifts = [0, 0, 0, 0, 0] # no shift

Contours = [InterpContour, InterpSSContour, ActualContour]
Labels = [f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Actual contour at {InterpSliceInd}']
Colours = ['r', 'm', 'y']
#Shifts = [0, 2, 4, 6, 8] # to help visualise
Shifts = [0, 0, 0, 0, 0] # no shift

Contours = [InterpContour, InterpSSContour, ReducedInterpSSContour, ActualContour]
Labels = [f'Interpolated contour at {InterpSliceInd}', 
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd} with \nreduced nodes',
          f'Actual contour at {InterpSliceInd}']
Colours = ['r', 'm', 'c', 'y']
#Shifts = [0, 2, 4, 6, 8] # to help visualise
Shifts = [0, 0, 2, 0, 0] # no shift


cif.PlotInterpolationResults(Contours=Contours, Labels=Labels, 
                             Colours=Colours, Shifts=Shifts,
                             InterpSliceInd=InterpSliceInd,
                             BoundingSliceInds=BoundingSliceInds,
                             dP=dP, ExportPlot=ExportPlot)


# In[ ]:





# In[ ]:





# In[215]:


test = [3, 1, 4, 2]

print(np.argsort(test))

print([i[0] for i in sorted(enumerate(test), key=lambda x:x[1])])


# In[ ]:





# In[17]:


Contour1.append(Contour1[0])
Contour1


# In[101]:


print(3/2)
print(3//2)
print(-3//2)
print(-(-3//2))


# In[ ]:




