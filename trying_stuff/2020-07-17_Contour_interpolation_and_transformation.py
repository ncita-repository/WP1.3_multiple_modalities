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
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('matplotlib', 'notebook')
#from ipywidgets import interact, fixed
#from IPython.display import clear_output

import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import GetInputPoints
importlib.reload(GetInputPoints)
import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)

#import GetImageAttributes
#importlib.reload(GetImageAttributes)
from GetImageAttributes import GetImageAttributes
from GetInputPoints import GetInputPoints
#import CreateInputFileForElastix
#importlib.reload(CreateInputFileForElastix)
from CreateInputFileForElastix import CreateInputFileForElastix
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
from ParseTransformixOutput import ParseTransformixOutput
from GetOutputPoints import GetOutputPoints
from PCStoICS import PCStoICS
import RegUtilityFuncs as ruf
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
FixPtsPCS, FixPtsBySliceAndContourPCS,FixPtsICS, FixPtsBySliceAndContourICS,LUT, PointData, FixContourData = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                                                Origin=FixOrigin, Directions=FixDirs,
                                                Spacings=FixSpacings)


# In[2]:


# Convert the dictionaries into data frames:
PointDataDF = pd.DataFrame(data=PointData)
FixContourDataDF = pd.DataFrame(data=FixContourData)

PointDataDF


# In[3]:


FixContourDataDF


# In[21]:


""" 
Note:

a//b yields floor(a/b)

- (-a//b) yields ceil(a/b)

"""

print(f'floor of {InterpSliceInd} is {InterpSliceInd//1}')
print(f'ceil of {InterpSliceInd} is {-(-InterpSliceInd//1)}')
print('')
print(f'int of floor of {InterpSliceInd} is {int(InterpSliceInd//1)}')
print(f'int of ceil of {InterpSliceInd} is {int(-(-InterpSliceInd//1))}')

int(InterpSliceInd//1)


# In[72]:


list(set(PointData['InSliceNo']))


# In[64]:


for i in range(len(FixContourData)):
    print(FixContourData['ContourType'][i])
    if FixContourData['ContourType'][i] == 1:
        print(True)
    else:
        print(False)


# In[4]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

# Define the slice(s) to be interpolated:
InterpSliceInd = 14
InterpSliceInd = 13.5

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Decide whether or not to interpolate between all contour point pairs
# or only those for which the node of the longest list of contour points
# was an original node (see Note 4 in the function 
# InterpolateBetweenContours):
InterpolateAllPts = True

# Get the indices of the slices that contain contours:
SlicesWithContours = cif.GetIndsOfSlicesWithContours(ContourData=FixContourData)

print('Contours exist on slices', SlicesWithContours)

#InterpContour = cif.InterpolateContours(FixContourData, PointData, InterpSliceInd, dP)

BoundingSliceInds,OSContour1,OSIsOrigNode1,OSContour2,OSIsOrigNode2,InterpContourPCS = cif.InterpolateContours(FixContourData, PointData, InterpSliceInd,
                                           dP, InterpolateAllPts)

# Continue only if None was not returned for InterpContourPCS:
if InterpContourPCS:
    #print(f'len(InterpContourPCS) = {len(InterpContourPCS)}')
    
    # The maximum contour number in the original contour data:
    LastCntNo = PointData['InContourNo'][-1]

    # Add interpolated contour to PointData:
    for i in range(len(InterpContourPCS)):
        PointData['PointNo'].append(PointData['PointNo'][-1] + 1)
        PointData['InSliceNo'].append(InterpSliceInd)
        PointData['InContourNo'].append(LastCntNo + 1)
        PointData['InContourType'].append(2)
        PointData['InPointIndex'].append([])
        PointData['InPointPCS'].append(InterpContourPCS[i])
        PointData['InPointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS[i], 
                                                Origin=FixOrigin, 
                                                Directions=FixDirs, 
                                                Spacings=FixSpacings))

    # Add interpolated contour to FixContourData:
    FixContourData['SliceNo'].append(InterpSliceInd)
    FixContourData['ContourType'].append(2)
    #FixContourData['PointPCS'].append(InterpContourPCS)
    #FixContourData['PointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS, 
    #                                           Origin=FixOrigin, 
    #                                           Directions=FixDirs, 
    #                                           Spacings=FixSpacings))
    """
    The list of contour points needs to be in a list structure to be consistent with 
    the other contour data in the dictionary.  Should verify this (July 23).
    """
    FixContourData['PointPCS'].append([InterpContourPCS])
    InterpContourICS = PCStoICS(Pts_PCS=InterpContourPCS, 
                                Origin=FixOrigin, 
                                Directions=FixDirs, 
                                Spacings=FixSpacings)
    FixContourData['PointICS'].append([InterpContourICS])

    """ GET ADDITIONAL DATA FOR PLOTTING """


    # Get the bounding slice indices for the slice at InterpSliceInd:
    BoundingSliceInds = cif.GetBoundingSliceInds(FixContourData, PointData, InterpSliceInd)

    Contour1 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[0]][0])
    Contour2 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[1]][0])

    if InterpSliceInd.is_integer():
        # The actual contour that exists at slice InterpSliceInd:
        ActualContour = copy.deepcopy(FixContourData['PointPCS'][InterpSliceInd][0])




    """ PLOT RESULTS """

    ExportPlot = False
    #ExportPlot = True

    #Contour1Shifted = cif.ShiftContourPoints(Contour=Contour1, Shift=3)


    #shift = cif.FindIndForMinCumLength(Contour1=SSContour1, Contour2=SSContour2)
    #print(f'shift = {shift}')
    #SSContour2Shifted = cif.ShiftContourPoints(Contour=SSContour2, Shift=shift)


    # Define contours to plot, their labels and their line colours:
    #Contours = [Contour1, Contour2, InterpContour, InterpSSContour, ActualContour]
    #Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
    #          f'Contour at slice {BoundingSliceInds[1]}', 
    #          f'Interpolated contour at {InterpSliceInd}', 
    #          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
    #          f'Actual contour at {InterpSliceInd}']
    #Colours = ['b', 'g', 'r', 'm', 'y']
    #Shifts = [0, 2, 4, 6, 8] # to help visualise
    #Shifts = [0, 0, 0, 0, 0] # no shift

    #Contours = [InterpContour, InterpSSContour, ActualContour]
    #Labels = [f'Interpolated contour at {InterpSliceInd}', 
    #          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
    #          f'Actual contour at {InterpSliceInd}']
    #Colours = ['r', 'm', 'y']
    #Shifts = [0, 2, 4, 6, 8] # to help visualise
    #Shifts = [0, 0, 0, 0, 0] # no shift

    #Contours = [InterpContour, InterpSSContour, ReducedInterpSSContour, ActualContour]
    #Labels = [f'Interpolated contour at {InterpSliceInd}', 
    #          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd}',
    #          f'Super-sampled interpolated \ncontour at slice {InterpSliceInd} with \nreduced nodes',
    #          f'Actual contour at {InterpSliceInd}']
    #Colours = ['r', 'm', 'c', 'y']
    #Shifts = [0, 2, 4, 6, 8] # to help visualise
    #Shifts = [0, 0, 2, 0, 0] # no shift

    #Contours = [Contour1, Contour2, InterpContour, ActualContour]
    #Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
    #          f'Contour at slice {BoundingSliceInds[1]}', 
    #          f'Interpolated contour at {InterpSliceInd}', 
    #          f'Actual contour at {InterpSliceInd}']
    #Colours = ['b', 'g', 'r', 'y']
    #Shifts = [0, 2, 4, 6] # to help visualise
    #Shifts = [0, 0, 0, 0] # no shift

    Contours = [Contour1, Contour2, InterpContour]
    Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
              f'Contour at slice {BoundingSliceInds[1]}', 
              f'Interpolated contour at {InterpSliceInd}']
    Colours = ['b', 'g', 'r', 'y']
    #Shifts = [0, 2, 4, 6] # to help visualise
    Shifts = [0, 0, 0, 0] # no shift
    if InterpSliceInd.is_integer():
        Contours.append(ActualContour)
        Labels.append(f'Actual contour at {InterpSliceInd}')

    cif.PlotInterpolationResults2D(Contours=Contours, Labels=Labels, 
                                   Colours=Colours, Shifts=Shifts,
                                   InterpSliceInd=InterpSliceInd,
                                   BoundingSliceInds=BoundingSliceInds,
                                   dP=dP, ExportPlot=ExportPlot)


# In[109]:


print(SlicesWithContours)

for i in range(len(SlicesWithContours) - 1):
    print(SlicesWithContours[i] + 0.5)


# ### Same as above but iterate through all slices with contours and interpolate at the midway slice:

# In[4]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif 

# Define the slice(s) to be interpolated:
#InterpSliceInd = 14
#InterpSliceInd = 13.5

# Chose the minimum inter-node distance for super-sampled contours:
dP = 0.2
#dP = 0.02

# Decide whether or not to interpolate between all contour point pairs
# or only those for which the node of the longest list of contour points
# was an original node (see Note 4 in the function 
# InterpolateBetweenContours):
InterpolateAllPts = True

# Get the indices of the slices that contain contours:
SlicesWithContours = cif.GetIndsOfSlicesWithContours(ContourData=FixContourData)

print('Contours exist on slices', SlicesWithContours)

# Initialise a dictionary to store the interpolation results:
InterpData = {'InterpSliceInd'       : [],
              'BoundingSliceInds'    : [],
              'FixOSIsOrigNode1'     : [],
              'FixOSIsOrigNode2'     : [],
              'FixOSContour1Inds'    : [],
              'FixOSContour2Inds'    : [],
              'FixInterpContourInds' : [],
              'FixOSContour1Pts'     : [],
              'FixOSContour2Pts'     : [],
              'FixInterpContourPts'  : []
              }

# Iterate through all SlicesWithContours and interpolate between each
# two slices:
for i in range(len(SlicesWithContours) - 1):
    InterpSliceInd = SlicesWithContours[i] + 0.5
    
    print(f'Attempting to interpolate at slice {InterpSliceInd}...')


    #InterpContour = cif.InterpolateContours(ContourData, PointData, InterpSliceInd, dP)
    
    BoundingSliceInds,    OSContour1,    OSIsOrigNode1,    OSContour2,    OSIsOrigNode2,    InterpContourPCS = cif.InterpolateContours(FixContourData, PointData, InterpSliceInd, 
                                               dP, InterpolateAllPts)

    # Continue only if None was not returned for InterpContourPCS:
    if InterpContourPCS:
        #print(f'len(InterpContourPCS) = {len(InterpContourPCS)}')

        # The maximum contour number in the original contour data:
        LastCntNo = PointData['InContourNo'][-1]

        # Add interpolated contour to PointData:
        for i in range(len(InterpContourPCS)):
            PointData['PointNo'].append(PointData['PointNo'][-1] + 1)
            PointData['InSliceNo'].append(InterpSliceInd)
            PointData['InContourNo'].append(LastCntNo + 1)
            PointData['InContourType'].append(2)
            PointData['InPointIndex'].append([])
            PointData['InPointPCS'].append(InterpContourPCS[i])
            PointData['InPointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS[i], 
                                                    Origin=FixOrigin, 
                                                    Directions=FixDirs, 
                                                    Spacings=FixSpacings))

        # Add interpolated contour to FixContourData:
        FixContourData['SliceNo'].append(InterpSliceInd)
        FixContourData['ContourType'].append(2)
        #FixContourData['PointPCS'].append(InterpContourPCS)
        #FixContourData['PointICS'].append(PCStoICS(Pts_PCS=InterpContourPCS, 
        #                                           Origin=FixOrigin, 
        #                                           Directions=FixDirs, 
        #                                           Spacings=FixSpacings))
        """
        The list of contour points needs to be in a list structure to be consistent with 
        the other contour data in the dictionary.  Should verify this (July 23).
        """
        FixContourData['PointPCS'].append([InterpContourPCS])
        InterpContourICS = PCStoICS(Pts_PCS=InterpContourPCS, 
                                    Origin=FixOrigin, 
                                    Directions=FixDirs, 
                                    Spacings=FixSpacings)
        FixContourData['PointICS'].append([InterpContourICS])
    
        
        # Store the interpolation output in InterpData:
        InterpData['InterpSliceInd'].append(InterpSliceInd)
        InterpData['BoundingSliceInds'].append(BoundingSliceInds)
        InterpData['FixOSIsOrigNode1'].append(OSIsOrigNode1)
        InterpData['FixOSIsOrigNode2'].append(OSIsOrigNode2)
        InterpData['FixOSContour1Pts'].append(OSContour1)
        InterpData['FixOSContour2Pts'].append(OSContour2)
        InterpData['FixInterpContourPts'].append(InterpContourPCS)


        """ GET ADDITIONAL DATA FOR PLOTTING """


        # Get the bounding slice indices for the slice at InterpSliceInd:
        BoundingSliceInds = cif.GetBoundingSliceInds(FixContourData, PointData, InterpSliceInd)

        Contour1 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[0]][0])
        Contour2 = copy.deepcopy(FixContourData['PointPCS'][BoundingSliceInds[1]][0])

        if InterpSliceInd.is_integer():
            # The actual contour that exists at slice InterpSliceInd:
            ActualContour = copy.deepcopy(FixContourData['PointPCS'][InterpSliceInd][0])




        """ PLOT RESULTS """

        ExportPlot = False
        #ExportPlot = True

        Contours = [Contour1, Contour2, InterpContourPCS]
        Labels = [f'Contour at slice {BoundingSliceInds[0]}', 
                  f'Contour at slice {BoundingSliceInds[1]}', 
                  f'Interpolated contour at {InterpSliceInd}']
        Colours = ['b', 'g', 'r', 'y']
        #Shifts = [0, 2, 4, 6] # to help visualise
        Shifts = [0, 0, 0, 0] # no shift
        if InterpSliceInd.is_integer():
            Contours.append(ActualContour)
            Labels.append(f'Actual contour at {InterpSliceInd}')

        cif.PlotInterpolationResults2D(Contours=Contours, Labels=Labels, 
                                       Colours=Colours, Shifts=Shifts,
                                       InterpSliceInd=InterpSliceInd,
                                       BoundingSliceInds=BoundingSliceInds,
                                       dP=dP, ExportPlot=ExportPlot)


# print(list(PointData.keys()), '\n')
# 
# for key, values in PointData.items():
#     print(key, '=', values)

# print(list(ContourData.keys()), '\n')
# 
# for key, values in ContourData.items():
#     print(key, '=', values)

# In[5]:


# Convert the dictionaries into data frames:
PointDataDF = pd.DataFrame(data=PointData)
FixContourDataDF = pd.DataFrame(data=FixContourData)
#InterpDataDF = pd.DataFrame(data=InterpData)

PointDataDF


# In[6]:


FixContourDataDF


# In[13]:


#InterpDataDF


# In[21]:


len(InterpData)


# ### Transform the contour points in InterpData (July 20)

# In[28]:


BaseKeys = ['OSContour1', 'OSContour2', 'InterpContour']

FixPtsKeys = ['Fix' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

FixIndsKeys = ['Fix' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]

MovPtsKeys = ['Mov' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

MovIndsKeys = ['Mov' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]


print(BaseKeys)
print(FixPtsKeys)
print(FixIndsKeys)
print(MovPtsKeys)
print(MovIndsKeys)


# In[7]:


# Create lists of keys that can be looped through:
BaseKeys = ['OSContour1', 'OSContour2', 'InterpContour']

FixIndsKeys = ['Fix' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
FixPtsKeys = ['Fix' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

MovIndsKeys = ['Mov' + BaseKeys[i] + 'Inds' for i in range(len(BaseKeys))]
MovPtsKeys = ['Mov' + BaseKeys[i] + 'Pts' for i in range(len(BaseKeys))]

# Loop through each row in InterpData and transform the points in FixPtsKeys:
for i in range(len(InterpData['InterpSliceInd'])):
    
    
    for k in range(len(FixPtsKeys)):
        Points = InterpData[FixPtsKeys[k]][i]
        
        # Create inputpoints.txt for Elastix, containing the contour points to be
        # transformed:
        CreateInputFileForElastix(Points=Points)

        # Transform MovingContourPts:
        TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

        # Parse outputpoints.txt:
        PtNos, FixInds, FixPts_PCS,        FixOutInds, MovPts_PCS,        Def_PCS, MovInds = ParseTransformixOutput()
        
        # Initialise the keys to be added to InterpData:
        #InterpData.update({MovIndsKeys[k] : [],
        #                   MovPtsKeys[k]  : []
        #                   })
        
        # Append the transformed results to InterpData:
        InterpData[FixIndsKeys[k]].append(FixInds) # the indices of the input (fixed) points
        #InterpData[MovIndsKeys[k]].append(MovInds) # the indices of the output (moving/transformed) points
        #InterpData[MovPtsKeys[k]].append(MovPts_PCS) # the output (moving/transformed) points
        
        # Get list of keys in InterpData:
        keys = list(InterpData.keys())
        
        if MovIndsKeys[k] in keys:
            InterpData[MovIndsKeys[k]].append(MovInds) # the indices of the output (moving/transformed) points
        else:
            InterpData.update({MovIndsKeys[k] : [MovInds]}) # the indices of the output (moving/transformed) points
            
        if MovPtsKeys[k] in keys:
            InterpData[MovPtsKeys[k]].append(MovPts_PCS) # the output (moving/transformed) points
        else:
            InterpData.update({MovPtsKeys[k] : [MovPts_PCS]}) # the output (moving/transformed) points


# In[56]:


for key, values in InterpData.items():
    print(key, f'has {len(values)} items in list')


# In[8]:


# Convert the InterpData into a data frame:
InterpDataDF = pd.DataFrame(data=InterpData)

InterpDataDF


# ### Find the intersection of the "lines" that would result from joining the points from one contour to another with the imaging planes in Moving Image.

# In[7]:


from GetVectorLength import GetVectorLength

print('MovDims =', MovDims)
print('')
print('MovOrigin =', MovOrigin)
print('')
print('MovSpacings =', MovSpacings)
print('')
print(f'MovDirs = [{MovDirs[0]}, {MovDirs[1]}, {MovDirs[2]}, \n', 
      f'         {MovDirs[3]}, {MovDirs[4]}, {MovDirs[5]}, \n', 
      f'         {MovDirs[6]}, {MovDirs[7]}, {MovDirs[8]}]')
print('')
print('Length of MovDirs_x =', GetVectorLength(MovDirs[0:3]))
print('')
print('Length of MovDirs_y =', GetVectorLength(MovDirs[3:6]))
print('')
print('Length of MovDirs_z =', GetVectorLength(MovDirs[6:9]))
print('')
print('Length of MovDirs =', ( MovDirs[0]**2 + MovDirs[1]**2 + MovDirs[2]**2 
       + MovDirs[3]**2 + MovDirs[4]**2 + MovDirs[5]**2 
       + MovDirs[6]**2 + MovDirs[7]**2 + MovDirs[8]**2)**(1/2))
print('')
print('Length of [MovDirs[0], MovDirs[4], MovDirs[8]] =', 
      GetVectorLength([MovDirs[0], MovDirs[4], MovDirs[8]]))


# In[15]:


MovIA_pydi = GetImageAttributes(MovDicomDir, 'pydicom')
MovIA_pydi


# In[16]:


MovIA_sitk = GetImageAttributes(MovDicomDir, 'sitk')
MovIA_sitk


# In[47]:


from NormaliseVector import NormaliseVector

print('Vector length of [MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]] \n=',
      GetVectorLength([MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]]))
print('')
print('Vector length of [MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]] \n=',
      GetVectorLength([MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]]))
print('')
print('Normalised [MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]] \n=',
      NormaliseVector([MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]]))
print('')
print('Normalised [MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]] \n=',
      NormaliseVector([MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]]))
print('')
print('Length of mormalised [MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]] \n=',
      GetVectorLength(NormaliseVector([MovIA_pydi[1][0], MovIA_pydi[1][4], MovIA_pydi[1][8]])))
print('')
print('Length of normalised [MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]] \n=',
      GetVectorLength(NormaliseVector([MovIA_sitk[1][0], MovIA_sitk[1][4], MovIA_sitk[1][8]])))


# In[52]:


"""
Each list in InterpData['MovOSContour1Pts'], InterpData['MovOSContour2Pts'] 
and InterpData['MovInterpContourPts'] has the list of contour points for each 
pair the contours used for each interpolated contour. 
"""

InterpData['MovOSContour1Pts'][0][0]


# In[111]:


import ContourInterpolatingFuncs
importlib.reload(ContourInterpolatingFuncs)
import ContourInterpolatingFuncs as cif

get_ipython().run_line_magic('matplotlib', 'notebook')

ExportPlot = False
#ExportPlot = True

cif.PlotInterpolationResults3D(InterpData=InterpData, dP=dP, ExportPlot=ExportPlot)


# In[9]:


InterpDataDF


# In[59]:


#print(j)

for i in range(len(InterpData['InterpSliceInd'])):
    print('len Contour1      =', len(InterpData['MovOSContour1Pts'][i]))
    print('len InterpContour =', len(InterpData['MovInterpContourPts'][i]))
    print('len Contour2      =', len(InterpData['MovOSContour2Pts'][i]), '\n\n')


# In[11]:


import GetLinePlaneIntersection
importlib.reload(GetLinePlaneIntersection)
from GetLinePlaneIntersection import GetLinePlaneIntersection


# The imaging plane normal:
PlaneNormal = MovDirs[6:9]

# Initialise a list of contours consisting of the intersecting points 
# for each imaging plane:
IntersContours = [[] for i in range(MovDims[2])]


# Loop through each set of interpolation data in InterpData:
for i in range(len(InterpData['InterpSliceInd'])):

    # Loop through each triplet of points in Contour1, InterpContour and Contour2:
    for j in range(len(InterpData['MovOSContour1Pts'][i])):

        Point1 = InterpData['MovOSContour1Pts'][i][j]
        Point2 = InterpData['MovInterpContourPts'][i][j]
        Point3 = InterpData['MovOSContour2Pts'][i][j]
        
        # Get the line directions for the line that joins Point1 and Point2, 
        # and the line that joins Point2 and Point3:
        LineDir1 = [Point2[n] - Point1[n] for n in range(len(Point1))]
        LineDir2 = [Point3[n] - Point2[n] for n in range(len(Point2))]
        
        # Store the intersection points for the Line1 and each image plane,
        # and for Line2 and each image plane:
        IntersPoints1 = []
        IntersPoints2 = []

        # Loop through all imaging planes in Moving Image and get the intersecting
        # points of Line1 and Line2 with each plane:
        for k in range(MovDims[2]):

            # Need a point lying on the plane - use the plane's origin:
            PlaneOrigin = [MovOrigin[0], MovOrigin[1], MovOrigin[2] + k*MovSpacings[2]]

            #print(f'Length of PlaneNormal is {GetVectorLength(PlaneNormal)}')

            IntersPoint1 = GetLinePlaneIntersection(PlaneNormal=PlaneNormal, 
                                                    PlanePoint=PlaneOrigin,
                                                    LineDirection=LineDir1,
                                                    LinePoint=Point1)
            
            IntersPoint2 = GetLinePlaneIntersection(PlaneNormal=PlaneNormal, 
                                                    PlanePoint=PlaneOrigin,
                                                    LineDirection=LineDir2,
                                                    LinePoint=Point2)
                

            IntersPoints1.append(IntersPoint1)
            IntersPoints2.append(IntersPoint2)
            
            if False:
                print('The intersection of Line1 joining point', Point1, 'in slice', 
                      InterpData['BoundingSliceInds'][i][0], 'and \npoint', Point2,
                      'in slice', InterpData['InterpSliceInd'][i], 'and the', 
                      f'{i}^th image plane with origin\n', PlaneOrigin, 'is', 
                      IntersPoint1, '.\n\n')

                print('The intersection of Line2 joining point', Point2, 'in slice', 
                      InterpData['InterpSliceInd'][i], 'and \npoint', Point3,
                      'in slice', InterpData['BoundingSliceInds'][i][1], 'and the',  
                      f'{i}^th image plane with origin\n', PlaneOrigin, 'is', 
                      IntersPoint2, '.\n\n')



        
        # Find if any of the intersection points lie between Point1 and Point2, and
        # between Point2 and Point3:
        for k in range(MovDims[2]):
            if ( IntersPoints1[k][2] > Point1[2] ) and ( IntersPoints1[k][2] < Point2[2] ):
                IntersContours[k].append(IntersPoints1[k])
                
            if ( IntersPoints2[k][2] > Point2[2] ) and ( IntersPoints2[k][2] < Point3[2] ):
                IntersContours[k].append(IntersPoints2[k])
                
                
# IntersContours has lists of points by contour. Modify it so that each slice has a list of
# a list of points (whereby the first list represents the list of points that belong to a
# contour, and the list of points follow), so that it's consistent with the structure of
# FixContourData:
IntersContoursPCS = []

for PtsThisSlice in IntersContours:    
    if PtsThisSlice: # i.e. if not empty
        IntersContoursPCS.append([PtsThisSlice])
    else:
        IntersContoursPCS.append([])

IntersContoursPCS


# In[ ]:





# In[12]:


FixContourDataDF


# In[15]:


# Create a dictionary to store contour points for Moving image (obtained by
# transforming the input/Fixed points, interpolating and finding intersection
# points with the image planes in Moving image) arranged by slices and contours:
"""
Note:
    Assign the number 3 for ContourType for all slices that have contours and 0
    for all slices that have no contours. 
    
    N.B. ContourType 0 -> No contours
                     1 -> Original contour
                     2 -> Interpolated contour
                     3 -> Contour obtained by transforming original or interpolated
                          points in the Fixed image domain and finding the intersection
                          with the Moving image planes
    
"""

# Create a list of slice numbers:
SliceNo = list(range(MovDims[2]))

# Initialise ContourType:
ContourType = []

for contour in IntersContours:
    if contour:
        ContourType.append(3)
        
    else:
        ContourType.append(0)

# Initialise IntersContoursICS:
IntersContoursICS = []

for i in range(len(IntersContoursPCS)):
    PtsThisSlicePCS = IntersContoursPCS[i]
    
    if PtsThisSlicePCS:
        # Convert points from PCS to ICS:
        """
        Note: There should only be one contour since the contour
        interpolation code can currently only deal with one contour
        per slice, but to future-proof this loop through all contours.
        """ 
        PtsThisSliceICS = []
        
        for c in range(len(PtsThisSlicePCS)):
            PtsThisSliceThisContourPCS = PtsThisSlicePCS[c]
            
            PtsThisSliceThisContourICS = PCStoICS(Pts_PCS=PtsThisSliceThisContourPCS, 
                                                  Origin=MovOrigin, 
                                                  Directions=MovDirs, 
                                                  Spacings=MovSpacings
                                                  )
                                                  
            PtsThisSliceICS.append(PtsThisSliceThisContourICS)
                                                  
            
        IntersContoursICS.append(PtsThisSliceICS)
                                 
        
    else:
        IntersContoursICS.append([])


MovContourData = {
                  'SliceNo'     : SliceNo, 
                  'ContourType' : ContourType,
                  'PointPCS'    : IntersContoursPCS,
                  'PointICS'    : IntersContoursICS
                   }


# In[16]:


MovContourDataDF = pd.DataFrame(data=MovContourData)

MovContourDataDF


# In[30]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf


# Plot the results:

# Chose which slices to plot:
PlotSlices = -1 # all
PlotSlices = [11]
PlotSlices = [12]
PlotSlices = [13]
PlotSlices = [12, 13, 14, 15, 16, 17]

# Choose the perspective:
Perspective = 'axial'
#Perspective = 'sagittal'
#Perspective = 'coronal'

# Choose whether to plot contour points as lines or dots:
ContoursAs = 'lines'
#ContoursAs = 'dots'

# Decide whether to export the figure or not:
ExportFig = False
#ExportFig = True

# Generate the filename for the exported figure:
FigFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime())            + f'_Result_of_registration_and_interpolation.png'


# Choose whether or not to log to the console:
LogToConsole = True
LogToConsole = False


ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, 
                                                                 mov_im=MovIm, 
                                                                 reg_im=RegIm,
                                                                 fix_pts=FixContourData['PointICS'], 
                                                                 mov_pts=MovContourData['PointICS'],
                                                                 export_fig=ExportFig,
                                                                 export_fname=FigFname,
                                                                 plot_slices=PlotSlices,
                                                                 perspective=Perspective,
                                                                 contours_as=ContoursAs,
                                                                 LogToConsole=LogToConsole)


# In[70]:


FixContourData['PointICS']


# In[35]:


i = 13

FixContourData['PointICS'][i]


# In[36]:



LenThisSlice = len(FixContourData['PointICS'][i])

print(f'LenThisSlice = {LenThisSlice}')

# Find out if this slice has LenThisSlice contours or if it
# has a single contour with LenThisSlice points:


# In[38]:


np.array(FixContourData['PointICS'][i]).shape


# In[39]:


np.array(FixContourData['PointICS']).shape


# In[66]:


print(np.array([]).shape)
print(len(np.array([]).shape))
print(np.array([]).shape[0])


# In[64]:


print(np.array([ [1,1,1], [1,1,1] ]).shape)
print(len(np.array([ [1,1,1], [1,1,1] ]).shape))
print(np.array([ [1,1,1], [1,1,1] ]).shape[0])
print(np.array([ [1,1,1], [1,1,1] ]).shape[1])


# In[68]:


print(np.array([ [ [1,1,1], [1,1,1], [1,1,1] ], [ [2,2,2], [2,2,2] ] ]).shape)
print(len(np.array([ [ [1,1,1], [1,1,1], [1,1,1] ], [ [2,2,2], [2,2,2] ] ]).shape))
print(np.array([ [ [1,1,1], [1,1,1], [1,1,1] ], [ [2,2,2], [2,2,2] ] ]).shape[0])


# In[17]:


for i in range(len(FixContourData['PointICS'])):
    print(np.array(FixContourData['PointICS'][i]).shape)


# In[18]:


for i in range(len(MovContourData['PointICS'])):
    print(np.array(MovContourData['PointICS'][i]).shape)


# In[ ]:




