#!/usr/bin/env python
# coding: utf-8

# # Guidance
# 
# 
# ## Data
# 
# This is a stripped down version of previous version "2020-06-24_Exercise_trying_3D_registration_using_SimpleElastix.ipynb" to make it easier to follow.  For details of what was tried, why, what worked/failed, etc. refer to the May 22 version.
# 
# The data used is from Subject 011 from the ACRIN-FMISO-Brain data set obtained here:
# 
# https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948305#a6b63241cb4843a291c45ef8a62d2073
# 
# You can selectively download data for that subject only.  Specifically two series from two sessions were used:
# 
# 8-T1 SE AXIAL POST FS FC-59362 (Series 8) within 04-10-1960-MRI Brain wwo Contrast-69626
# 
# and
# 
# 8-T1 SE AXIAL POST FC-81428 (Series 8) within 06-11-1961-MRI Brain wwo Contrast-79433
# 
# An ROI Collection was manually created within the OHIF-Viewer for the series from 1960.  The DICOM-RTSTRUCT file can be obtained here:
# 
# https://drive.google.com/drive/folders/1CDrwTYvlvBL6rpuyz29zxO6AU4Ancmd6?usp=sharing
# 
# Originally the 1960 series (with ROI Collection) was defined as the "fixed" image, and the 1961 series as "moving".  This is fine but requires the inverse transformation used to register moving to fixed to deform the fixed image's contour points to the moving image.
# 
# I'm sure this can be done but I don't yet know how, so in the meantime fixed will be defined as the 1961 series, and moving the 1960 one.
# 
# ## SimpleElastix
# 
# I recommend trying to install SimpleElastix using a pre-compiled Python module obtained here:
# 
# https://www.francescosantini.com/wp/2020/03/06/simpleelastix-precompiled-python-modules/
# 
# If it doesn't work you will have to compile it from source:
# 
# https://simpleelastix.readthedocs.io/GettingStarted.html
# 
# There's a "more official" wheel provided here:
# 
# https://sourceforge.net/projects/simpleelastix/
# 
# but it might be for an older version of SimpleElastix (and I think I had trouble with this wheel and ended up having to compile it - I hadn't yet discovered Francesco Santini's page).
# 
# ## Other packages
# 
# Most of the other packages can be pip installed straightforwardly.  
# 
# The other imports are functions that can be obtained here:
# 
# https://github.com/ncita-repository/WP1.3_multiple_modalities/tree/master/trying_stuff

# 

# In[2]:


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

#from GetDicomFpaths import GetDicomFpaths
#from GetDicoms import GetDicoms
#import PCStoICS
#importlib.reload(PCStoICS)
#from PCStoICS import PCStoICS
#import GetContourPtsForDicoms
#importlib.reload(GetContourPtsForDicoms)
#from GetContourPtsForDicoms import GetContourPtsForDicoms
#from ParseTransformixOutput import ParseTransformixOutput




# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm'
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose which package to use to get the Image Plane Attributes:
package = 'sitk'
package = 'pydicom'


# Get the Image Attributes for the images:
FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                              Package=package)

MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                              Package=package)

# Get contour points into necessary arrays and dictionaries:
#FixPts_PCS, FixPtsArr_PCS, FixPtsDict_PCS = GetContourPtsForDicoms(DicomDir=FixDicomDir,
#                                                                   RoiFpath=FixRoiFpath)

InPts_PCS, InPtsArr_PCS, InPtsDict_PCS,InPts_ICS, InPtsArr_ICS, InPtsDict_ICS = GetInputPoints(DicomDir=FixDicomDir, 
                                                        RoiFpath=FixRoiFpath,
                                                        Origin=FixOrigin,
                                                        Directions=FixDirs,
                                                        Spacings=FixSpacings)

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

CoordSys = 'PCS'
#CoordSys = 'ICS'

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
CreateInputFileForElastix(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath, 
                          Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings,
                          CoordSys=CoordSys)

# Register:
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)

# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
#PtNos, InInds, InPts_PCS, FixOutInds,\
#OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

PtNos, InInds, InPts_PCS, InPts_ICS,FixOutInds, OutPts_PCS, OutPts_ICS,Defs_PCS, Defs_ICS, MovOutInds,OutPtsArr_PCS, OutPtsArr_ICS = GetOutputPoints(Origin=MovOrigin, Directions=MovDirs,
                                               Spacings=MovSpacings, Dimensions=MovDims)

print('\nDone.')


# In[3]:


# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]

# Choose whether to export figure:
export_fig = True
export_fig = False


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_June17_Attempt2_eqs.png'

ruf.display_all_sitk_images_and_reg_results_with_all_contours(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                              fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                              export_fig=export_fig,
                                                              export_fname=fig_fname,
                                                              plot_slices=plot_slices)


# ## Success!
# 
# ## Next steps:
# 
# ## 1) Store the transformed points into an RTSTRUCT file for importing into the OHIF-Viewer to see how they appear there.
# 
# ## 2) Tidy up the contours:
# 
# ### i.   Some/most/all contours require "closing"
# ### ii.  Some require removal of un-closed segments (e.g. 18th slice)
# ### iii. Some require avoiding segments crossing through the outer boundary (e.g. 15-17th slice)
# ### iv. Some require some awkward closing (e.g. 13th slice)

# ## Step 1:  Create ROI for Moving Image
# 
# ### This is a bit trickier than simply modifying a copy of the fixed ROI, since there are more contour sequences within the transformed points than in the original ROI.

# # Use the fixed ROI as a template for the moving image:
# FixRoi = pydicom.dcmread(FixRoiFpath)
# 
# MovRoi = copy.deepcopy(FixRoi)
# 
# #MovRoi
# #MovRoi.PatientName

# # Read in the DICOMs:
# FixDicomFpaths, FixDicoms = GetDicoms(DirPath=FixDicomDir, SortMethod='slices', Debug=False)
# 
# MovDicomFpaths, MovDicoms = GetDicoms(DirPath=MovDicomDir, SortMethod='slices', Debug=False)
# 
# FixDicoms[0]

# MovDicoms[0]

# print('FixRoi:')
# FixRoi

# print('MovRoi:')
# MovRoi

# import CreateMovRoi
# importlib.reload(CreateMovRoi)
# from CreateMovRoi import CreateMovRoi
# 
# MovRoi = CreateMovRoi(FixRoiFpath=FixRoiFpath, MovDicomDir=MovDicomDir, 
#                       MovPtsArr_ICS=OutPtsArr_ICS, Debug=False)

# MovRoi

# # Sagittal views of Fixed and Moving Images

# FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, Package='pydicom')
# MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, Package='pydicom')
# 
# print('FixDims:', FixDims)
# print('\nMovDims:', MovDims)
# 
# print('\nFixDirs:\n', FixDirs[0:3], '\n', FixDirs[3:6], '\n', FixDirs[6:])
# print('\nMovDirs:\n', MovDirs[0:3], '\n', MovDirs[3:6], '\n', MovDirs[6:])
# 
# # Check the aspect ratios in different views:
# FixAxAR = FixSpacings[1]/FixSpacings[0]
# MovAxAR = MovSpacings[1]/MovSpacings[0]
# 
# FixSagAR = FixSpacings[1]/FixSpacings[2]
# MovSagAR = MovSpacings[1]/MovSpacings[2]
# #FixSagAR = FixSpacings[2]/FixSpacings[1]
# #MovSagAR = MovSpacings[2]/MovSpacings[1]
# 
# FixCorAR = FixSpacings[2]/FixSpacings[0]
# MovCorAR = MovSpacings[2]/MovSpacings[0]
# 
# print('\nFixAxAR =', FixAxAR)
# print('FixSagAR =', FixSagAR)
# print('FixCorAR =', FixCorAR)
# print('\nMovAxAR =', MovAxAR)
# print('MovSagAR =', MovSagAR)
# print('MovCorAR =', MovCorAR)
# 
# # Create 3D array of DICOM images:
# Fix3D = np.zeros(FixDims)
# Mov3D = np.zeros(MovDims)
# 
# for i in range(len(FixDicoms)):
#     Fix3D[:, :, i] = FixDicoms[i].pixel_array
# 
# for i in range(len(MovDicoms)):
#     Mov3D[:, :, i] = MovDicoms[i].pixel_array
#     
# 

# # plot 3 orthogonal slices
# plt.subplots(1, 3, figsize=(14, 5))
# 
# a1 = plt.subplot(1, 3, 1)
# plt.imshow(Fix3D[:, :, FixDims[2]//2])
# a1.set_aspect(FixAxAR)
# 
# a2 = plt.subplot(1, 3, 2)
# plt.imshow(Fix3D[:, FixDims[1]//2, :])
# a2.set_aspect(FixSagAR)
# 
# a3 = plt.subplot(1, 3, 3)
# plt.imshow(Fix3D[FixDims[0]//2, :, :].T)
# a3.set_aspect(FixCorAR)
# 
# plt.show()

# # Sagittal plots:
# 
# Nrows = FixDims[2] # same as the number of slices
# Nrows = FixDims[2]//2
# Ncols = 1
# 
# # Define the slice thickness along sagital cuts:
# #dSag = int(np.floor(FixDims[0]/FixDims[2])) # same as below
# dSag = FixDims[0]//FixDims[2] # same as above
# 
# # Create list of indices for the sagittal slices:
# #yInds = np.arange(0, FixDims[0]+1, dSag) # better to use linspace..
# yInds = np.linspace(0, FixDims[1]-1, num=Nrows, endpoint=True, retstep=False, dtype=int)
# 
# 
# Rotate = False
# Rotate = True
# 
# plt.subplots(Nrows, Ncols, figsize=(5, 5*Nrows))
# #plt.subplots(Nrows, Ncols, figsize=(5, 50))
# 
# i = 0
# 
# for s in range(Nrows):
#     i = i + 1 # increment sub-plot number
#     
#     #print(f'\nyInd = {yInds[s]}')
#     
#     if Rotate:
#         arr = np.rot90(Fix3D[:, yInds[s], :])
#     else:
#         arr = Fix3D[:, yInds[s], :]
#     
#     ax1 = plt.subplot(Nrows, Ncols, i)
#     plt.imshow(arr, cmap=plt.cm.Greys_r);
#     plt.title(f'Fixed slice {s+1}/{Nrows}')
#     plt.axis('off')
#     if Rotate:
#         ax1.set_aspect(1/FixSagAR)
#     else:
#         ax1.set_aspect(FixSagAR)
# 

# FixNda = sitk.GetArrayFromImage(FixIm)
# 
# print('FixDims         =', FixDims)
# print('FixNda.shape    =', FixNda.shape)
# print('FixIm.GetSize() =', FixIm.GetSize())
# print('')
# print('FixDirs              =\n', FixDirs)
# print('\nFixIm.GetDirection() =\n', FixIm.GetDirection())
# print('')
# print(FixDirs[0:2])
# print([FixDirs[0], FixDirs[2]])
# print(FixDirs[1:3])
# 
# print('\n', type(FixIm.GetDirection()))
# #print('\n', type(list(FixIm.GetDirection())))
# print('\n', type(FixDirs))
# 
# #FixDirs = FixIm.GetDirection()
# #FixDirs = list(FixIm.GetDirection())
# FixDirs = [FixIm.GetDirection()[i] for i in range(len(FixIm.GetDirection()))]
# 
# print('\nFixDirs =\n', FixDirs)
# 
# a = [1, 2, 3, 5, 6, 7]
# #a = list([1, 2, 3, 5, 6, 7])
# 
# for i in a:
#     FixDirs[i] = - FixDirs[i]
# 
# print('\nFixDirs =\n', FixDirs)

# In[43]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
#plot_slices = [12, 13, 14, 15]


# Choose the perspective:
perspective = 'axial'
#perspective = 'sagittal'
#perspective = 'coronal'

# Choose whether to export figure:
export_fig = True
export_fig = False

# Choose whether to plot contour points as lines or dots:
contours_as = 'lines'
contours_as = 'dots'


# Create filename for exported figure:
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_' + perspective + '.png'

#ruf.display_all_sitk_images_and_reg_results_with_all_contours_v2(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
#                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
#                                                                 export_fig=export_fig,
#                                                                 export_fname=fig_fname,
#                                                                 plot_slices=plot_slices,
#                                                                 perspective=perspective)

ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
                                                                 export_fig=export_fig,
                                                                 export_fname=fig_fname,
                                                                 plot_slices=plot_slices,
                                                                 perspective=perspective,
                                                                 contours_as=contours_as)


# ### June 24:  Done producing sagittal (and coronal).
# 
# ### Next work on tidying up the contours.

# ### June 26:  Try using different input contours

# In[17]:


import DownloadRoiFromXnat
importlib.reload(DownloadRoiFromXnat)
from DownloadRoiFromXnat import DownloadRoiFromXnat

# Define the REST variables that correspond to the ROI Collection
# to be downloaded:

#ProjLabel = 'SarcomaTCIA'
#SubjLabel = 'STS_001'
#ExpLabel = 'STS_001_MR1'
#DateTime = '20200311_211225'

ProjLabel = 'ACRIN'
SubjLabel = 'ACRIN-FMISO-Brain-011'
ExpLabel = 'ACRIN-FMISO-Brain-011_MR_4'
DateTime = '20200511_073405' # tumour
DateTime = '20200626_104631' # ventricles

# Define download directory:
#DownloadDir = 'cwd'
DownloadDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'

Roi = DownloadRoiFromXnat(ProjectLabel=ProjLabel, SubjectLabel=SubjLabel, 
                          ExperimentLabel=ExpLabel, DateTimeInLabel=DateTime,
                          DownloadDir=DownloadDir)

#Roi


# In[112]:


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
import GetInputPoints
importlib.reload(GetInputPoints)
import CreateInputFileForElastix
importlib.reload(CreateInputFileForElastix)
import TransformPoints
importlib.reload(TransformPoints)
import GetOutputPoints
importlib.reload(GetOutputPoints)
import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)

from GetImageAttributes import GetImageAttributes
from GetInputPoints import GetInputPoints
from CreateInputFileForElastix import CreateInputFileForElastix
from RunSimpleElastixReg import RunSimpleElastixReg
from TransformPoints import TransformPoints
from GetOutputPoints import GetOutputPoints
import RegUtilityFuncs as ruf



# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
#FixRoiFname = r'AIM_20200511_073405.dcm' # tumour
FixRoiFname = r'AIM_20200626_104631.dcm' # ventricles
FixRoiFpath = os.path.join(FixRoiDir, FixRoiFname)

# Chose which package to use to get the Image Plane Attributes:
#package = 'sitk'
package = 'pydicom'


# Get the Image Attributes for the images:
FixOrigin, FixDirs, FixSpacings, FixDims = GetImageAttributes(DicomDir=FixDicomDir, 
                                                              Package=package)

MovOrigin, MovDirs, MovSpacings, MovDims = GetImageAttributes(DicomDir=MovDicomDir, 
                                                              Package=package)

# Get contour points into necessary arrays and dictionaries:
#FixPts_PCS, FixPtsArr_PCS, FixPtsDict_PCS = GetContourPtsForDicoms(DicomDir=FixDicomDir,
#                                                                   RoiFpath=FixRoiFpath)

FixPtsPCS, FixPtsBySlicePCS,FixPtsICS, FixPtsBySliceICS,LUT = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                     Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings)

#print('\nLUT[0] =', LUT[0])

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

# Register:
RegIm, ElastixImFilt = RunSimpleElastixReg(FixIm, MovIm)


""" 
Should use the contour points in the Patient Coordinate System since it
seems that Elastix already performs a conversion to the Image Coordinate 
System.
"""
CoordSys = 'PCS'
#CoordSys = 'ICS'

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
CreateInputFileForElastix(Points=FixPtsPCS)

# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
#PtNos, InInds, InPts_PCS, FixOutInds,\
#OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

PtNos, InInds, FixOutInds,MovPtsPCS, MovPtsICS,DefPCS, DefICS, MovOutInds,MovPtsBySlicePCS, MovPtsBySliceICS = GetOutputPoints(InPts_ICS=FixPtsICS, Origin=MovOrigin, 
                                                     Directions=MovDirs, Spacings=MovSpacings, 
                                                     Dimensions=MovDims, LUT=LUT)

if False:
    """ 
    Rather than transforming all points at once, transform the points on a contour-by-contour
    basis.
    """

    # Store the outputs in a slice-by-slice array:
    #PtNosAll = []
    #InIndsAll = []
    #FixOutIndsAll = []
    #MovPtsAll_PCS = []
    #MovPtsAll_ICS = []
    #DefsAll_PCS = []
    #DefsAll_ICS = []
    #MovOutIndsAll = []

    PtNosBySlice = []
    InIndsBySlice = []
    FixOutIndsBySlice = []
    MovPtsBySlice_PCS = []
    MovPtsBySlice_ICS = []
    DefsBySlice_PCS = []
    DefsBySlice_ICS = []
    MovOutIndsBySlice = []


    for SliceNo in range(len(FixPtsBySlice_PCS)):
        print(f'\nSliceNo = {SliceNo}')

        PtsThisSlice_PCS = FixPtsBySlice_PCS[SliceNo]
        PtsThisSlice_ICS = FixPtsBySlice_ICS[SliceNo]

        # Proceed if PtsThisSlice_PCS is not empty:
        if PtsThisSlice_PCS:

            # Store output for each contour for this slice:
            PtNosByContour = []
            InIndsByContour = []
            FixOutIndsByContour = []
            MovPtsByContour_PCS = []
            MovPtsByContour_ICS = []
            DefsByContour_PCS = []
            DefsByContour_ICS = []
            MovOutIndsByContour = []


            # Loop through each array of points in PtsThisSlice_PCS 
            # (for each contour for this slice):
            for ContourNo in range(len(PtsThisSlice_PCS)):
                print(f'\nContourNo = {ContourNo}')

                PtsThisContour_PCS = PtsThisSlice_PCS[ContourNo]
                PtsThisContour_ICS = PtsThisSlice_ICS[ContourNo]

                # Create inputpoints.txt for Elastix, containing the contour points to be
                # transformed:
                #CreateInputFileForElastix(Points=PtsThisSlice_PCS)
                CreateInputFileForElastix(Points=PtsThisContour_PCS)

                # Transform MovingContourPts:
                TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

                #PtNos, InInds, FixOutInds,\
                #MovPts_PCS, MovPts_ICS,\
                #Defs_PCS, Defs_ICS, MovOutInds,\
                #MovPtsArr_PCS, MovPtsArr_ICS = GetOutputPoints(InPts_ICS=PtsThisSlice_ICS, Origin=MovOrigin, 
                #                                               Directions=MovDirs, Spacings=MovSpacings, 
                #                                               Dimensions=MovDims, LUT=LUT)
                PtNos, InInds, FixOutInds,                MovPts_PCS, MovPts_ICS,                Defs_PCS, Defs_ICS, MovOutInds,                MovPtsArr_PCS, MovPtsArr_ICS = GetOutputPoints(InPts_ICS=PtsThisContour_ICS, Origin=MovOrigin, 
                                                               Directions=MovDirs, Spacings=MovSpacings, 
                                                               Dimensions=MovDims, LUT=LUT)


                PtNosByContour.append(PtNos)
                InIndsByContour.append(InInds)
                FixOutIndsByContour.append(FixOutInds)
                MovPtsByContour_PCS.append(MovPts_PCS)
                MovPtsByContour_ICS.append(MovPts_ICS)
                DefsByContour_PCS.append(Defs_PCS)
                DefsByContour_ICS.append(Defs_ICS)
                MovOutIndsByContour.append(MovOutInds)



            #PtNosAll.extend(PtNosByContour)
            #InIndsAll.extend(InIndsByContour)
            #FixOutIndsAll.extend(FixOutIndsByContour)
            #MovPtsAll_PCS.extend(MovPtsByContour_PCS)
            #MovPtsAll_ICS.extend(MovPtsByContour_ICS)
            #DefsAll_PCS.extend(DefsByContour_PCS)
            #DefsAll_ICS.extend(DefsByContour_ICS)
            #MovOutIndsAll.extend(MovOutIndsByContour)

            PtNosBySlice.append(PtNosByContour)
            InIndsBySlice.append(InIndsByContour)
            FixOutIndsBySlice.append(FixOutIndsByContour)
            MovPtsBySlice_PCS.append(MovPtsByContour_PCS)
            MovPtsBySlice_ICS.append(MovPtsByContour_ICS)
            DefsBySlice_PCS.append(DefsByContour_PCS)
            DefsBySlice_ICS.append(DefsByContour_ICS)
            MovOutIndsBySlice.append(MovOutIndsByContour)

        else:
            print('There are no contour points for this slice.')

            PtNosBySlice.append([])
            InIndsBySlice.append([])
            FixOutIndsBySlice.append([])
            MovPtsBySlice_PCS.append([])
            MovPtsBySlice_ICS.append([])
            DefsBySlice_PCS.append([])
            DefsBySlice_ICS.append([])
            MovOutIndsBySlice.append([])
        

print('\nDone.')


# In[83]:


len(LUT)


# In[70]:


len(FixPtsBySlice_PCS[13])


# In[49]:


FixPtsBySlice_PCS[13][0]


# In[25]:


len(FixPtsBySlice_PCS)


# In[50]:


PtNosBySlice


# In[51]:


FixOutIndsBySlice


# In[52]:


MovPtsBySlice_PCS


# In[ ]:





# In[ ]:





# In[3]:


FixPtsArr_ICS


# In[4]:


for r in range(len(FixPtsArr_ICS)):
    print(f'r = {r}, len = {len(FixPtsArr_ICS[r])}')
    
#for r in range(len(MovPtsArr_ICS)):
#    print(f'r = {r}, len = {len(MovPtsArr_ICS[r])}')


# In[5]:


MovPtsArr_ICS


# In[6]:


#for r in range(len(FixPtsArr_ICS)):
#    print(f'r = {r}, len = {len(FixPtsArr_ICS[r])}')
    
for r in range(len(MovPtsArr_ICS)):
    print(f'r = {r}, len = {len(MovPtsArr_ICS[r])}')


# In[7]:


LUT


# In[11]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
#plot_slices = [12, 13, 14, 15]
#plot_slices = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]


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
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_' + perspective + '.png'

#ruf.display_all_sitk_images_and_reg_results_with_all_contours_v2(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
#                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
#                                                                 export_fig=export_fig,
#                                                                 export_fname=fig_fname,
#                                                                 plot_slices=plot_slices,
#                                                                 perspective=perspective)

ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
                                                                 fix_pts=FixPtsBySliceICS, mov_pts=MovPtsBySliceICS,
                                                                 export_fig=export_fig,
                                                                 export_fname=fig_fname,
                                                                 plot_slices=plot_slices,
                                                                 perspective=perspective,
                                                                 contours_as=contours_as)


# In[18]:


Roi


# In[19]:


ContourSequences = Roi.ROIContourSequence[0].ContourSequence

N = len(ContourSequences)

for i in range(N):
    
    print(f'i = {i}, UID = ', ContourSequences[i].ContourImageSequence[0].ReferencedSOPInstanceUID)


# In[21]:


import GetContourPtsForDicoms
importlib.reload(GetContourPtsForDicoms)
from GetContourPtsForDicoms import GetContourPtsForDicoms

Pts_PCS, PtsArr_PCS, PtsDict_PCS, PtsLUT = GetContourPtsForDicoms(DicomDir=FixDicomDir, 
                                                          RoiFpath=FixRoiFpath)

for i in range(len(PtsArr_PCS)):
    print(f'i = {i}, len of array = {len(PtsArr_PCS[i])}')


# In[36]:


PtsArr_PCS


# In[ ]:





# In[14]:


Defs_ICS


# In[89]:


import GetContourPtsForDicoms
importlib.reload(GetContourPtsForDicoms)
from GetContourPtsForDicoms import GetContourPtsForDicoms

PtsPCS, PtsBySlicePCS, LUT = GetContourPtsForDicoms(DicomDir=FixDicomDir, 
                                                    RoiFpath=FixRoiFpath)


# In[92]:


PtsPCS
len(PtsPCS)


# In[94]:


PtsBySlicePCS
len(PtsBySlicePCS)


# In[96]:


LUT
#len(LUT)


# In[200]:


import GetInputPoints
importlib.reload(GetInputPoints)
from GetInputPoints import GetInputPoints

FixPtsPCS, FixPtsBySlicePCS,FixPtsICS, FixPtsBySliceICS,LUT = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                     Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings)


# In[201]:


LUT


# In[100]:


FixPtsPCS
len(FixPtsPCS)


# In[102]:


FixPtsBySlicePCS
len(FixPtsBySlicePCS)


# In[115]:


FixPtsICS
len(FixPtsICS)


# In[117]:


FixPtsBySliceICS
len(FixPtsBySliceICS)


# In[118]:


LUT
#len(LUT)


# In[71]:


InInds


# In[122]:


import GetOutputPoints
importlib.reload(GetOutputPoints)
from GetOutputPoints import GetOutputPoints

PtNos, InInds, FixOutInds,MovPtsPCS, MovPtsICS,DefPCS, DefICS, MovOutInds,MovPtsBySlicePCS, MovPtsBySliceICS = GetOutputPoints(InPts_ICS=FixPtsICS, Origin=MovOrigin, 
                                                     Directions=MovDirs, Spacings=MovSpacings, 
                                                     Dimensions=MovDims, LUT=LUT)


# In[125]:


[MovOutInds[i][2] for i in range(len(MovOutInds))]


# In[126]:


[LUT[i][1] for i in range(len(LUT))]


# In[131]:


for i in range(len(LUT)):
    print(f'Point {LUT[i][0]} in fixed contour {LUT[i][1]} --> MovZind {MovOutInds[i][2]}')


# In[191]:


# Get list of all moving slice numbers (z-indeces in MovOutInds):
MovSliceNos = [MovOutInds[i][2] for i in range(len(MovOutInds))]

# Get list of unique slice numbers:
MovSliceNosSet = list(set(MovSliceNos))

#print('MovSliceNosSet =', MovSliceNosSet)
# [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

# Initialise indeces of the points that correspond to the slice 
# numbers in MovSliceNosSet, and the corresponding fixed contour 
# numbers organised by moving slices:
PtIndsByMovSlice = []
FixContourNosByMovSlice = []

# Get indeces of all slice nos in MovSliceNos that are equal to each 
# index in MovSliceNosSet:
for zind in MovSliceNosSet:
    inds = [i for i, e in enumerate(MovSliceNos) if e==zind]
    
    PtIndsByMovSlice.append(inds)
    
    FixContourNosByMovSlice.append([LUT[ind][1] for ind in inds])
    
#print('\nPtIndsByMovSlice =', PtIndsByMovSlice)
print('\nPtIndsByMovSlice[0:2] =', [PtIndsByMovSlice[i] for i in range(2)])

#print('\nFixContourNosByMovSlice =', FixContourNosByMovSlice)

# Initialise array of moving contour numbers organised by moving slices, 
# and a moving contour number counter:
MovContourNosByMovSlice = []
MovContourNo = 0

# Loop through each array in FixContourNoByMovSlice:
for FixContourNosThisSlice in FixContourNosByMovSlice:
    
    # Initialise array of moving contour numbers for this slice:
    MovContourNosThisSlice = []
    
    for i in range(len(FixContourNosThisSlice)):
        ThisFixContourNo = FixContourNosThisSlice[i]
        
        # If i=0 this is the start of a new contour.
        if i == 0:
            MovContourNosThisSlice.append(MovContourNo)
        
        else:
            PreviousFixContourNo = FixContourNosThisSlice[i-1]
            
            # If this FixContourNo is the same as the previous:
            if ThisFixContourNo == PreviousFixContourNo:
                # This point is part of the previous contour:
                MovContourNosThisSlice.append(MovContourNo)
                
            else:
                # This is the start of a different contour:
                MovContourNo = MovContourNo + 1
                
                MovContourNosThisSlice.append(MovContourNo)
                
    
    MovContourNosByMovSlice.append(MovContourNosThisSlice)
    
    # Next slice will be start of next contour no:
    MovContourNo = MovContourNo + 1
        

#print('\nMovContourNosByMovSlice =', MovContourNosByMovSlice)  
print('\nMovContourNosByMovSlice[0:2] =', [MovContourNosByMovSlice[i] for i in range(2)])  


# Initialise a list of moving contour nos in the order of the points 
# themselves, by using the indeces in PtIndsByMovSlice to index 
# MovContourNosBySlice:
MovContourNosByPoint = []
[MovContourNosByPoint.append([]) for i in range(len(MovOutInds))]

# Loop through each array in MovContourNosByMovSlice:
for i in range(len(MovContourNosByMovSlice)):
    MovContourNosThisSlice = MovContourNosByMovSlice[i]
    
    inds = PtIndsByMovSlice[i]
    
    #print('\nMovContourNosThisSlice =', MovContourNosThisSlice)
    
    #print('\ninds =', inds)
    
    # Loop through each array in MovContourNosThisSlice:
    for j in range(len(MovContourNosThisSlice)):
        ThisMovContourNo = MovContourNosThisSlice[j]
        
        # The index in the list of points for this contour:
        ind = inds[j]
        
        MovContourNosByPoint[ind] = ThisMovContourNo
        

#print('\nMovContourNosByPoint =', MovContourNosByPoint) 
#print('\nMovContourNosByPoint[0:2] =', [MovContourNosByPoint[i] for i in range(2)]) 
N = sum([len(MovContourNosByMovSlice[i]) for i in range(2)])
print(f'\nMovContourNosByPoint[0:{N}] =', [MovContourNosByPoint[i] for i in range(N)]) 


"""
Create a list of point indeces (for each moving slice) that contains 
a sub-list of point indeces for each contour using MovContourNosByPoint.
"""

# Number of slices:
S = len(PtIndsByMovSlice)

# Initialise array of point indeces and moving contour numbers 
# slice-by-slice along rows and contour by contour along columns:
PtIndsByMovSliceAndContour = []
MovContourNosByMovSliceAndContour = []

# Loop through each slice in PtIndsByMovSlice:
for s in range(S):
    #if s < 2:
    #    print('\ns =', s)

    # The point indices and moving contour numbers for this slice:
    PtIndsThisSlice = PtIndsByMovSlice[s]
    MovContourNosThisSlice = MovContourNosByMovSlice[s]
    
    #if s < 2:
    #    print('\nMovContourNosThisSlice =', MovContourNosThisSlice)
    #    print('\nPtIndsThisSlice        =', PtIndsThisSlice)
    
    # A list of unique contour numbers for this slice:
    UniqueContourNos = list(set(MovContourNosThisSlice))
    
    if len(UniqueContourNos) == 1:
        # All items in this slice belong to the same contour:
        PtIndsByMovSliceAndContour.append(PtIndsThisSlice)
        MovContourNosByMovSliceAndContour.append(MovContourNosThisSlice)
        
    else:
        # The point indices and contour numbers need to be split and grouped
        # by contour number.     
        # Get the indices of each distinct contour number in MovContourNosThisSlice:
        IndsByContour = [[i for i, e in enumerate(MovContourNosThisSlice) if e==c] for c in UniqueContourNos]
        
        # Use each set of indeces in IndsByContour to index the point indices and moving
        # contour numbers in PtIndsThisSlice and MovContourNosThisSlice:
        PtIndsThisSliceByContour = [[PtIndsThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]
        MovContourNosThisSliceByContour = [[MovContourNosThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]
        
        #if s < 2:
        #    print('\nMovContourNosThisSliceByContour =', MovContourNosThisSliceByContour)
        #    print('\nPtIndsThisSliceByContour =', PtIndsThisSliceByContour)
            
        # Append to PtIndsByMovSliceAndContour and MovContourNosByMovSliceAndContour:
        PtIndsByMovSliceAndContour.append(PtIndsThisSliceByContour)
        MovContourNosByMovSliceAndContour.append(MovContourNosThisSliceByContour)
        
    #if s < 2:
    #    print(f'\nMovContourNosByMovSliceAndContour[{s}] =', MovContourNosByMovSliceAndContour[s])
    #    print(f'\nPtIndsByMovSliceAndContour[{s}] =', PtIndsByMovSliceAndContour[s])     
    
print('\nMovContourNosByMovSliceAndContour =', MovContourNosByMovSliceAndContour) 
#print('\nMovContourNosByMovSliceAndContour[0:2] =', [MovContourNosByMovSliceAndContour[i] for i in range(2)]) 

print('\nPtIndsByMovSliceAndContour =', PtIndsByMovSliceAndContour) 
#print('\nPtIndsByMovSliceAndContour[0:2] =', [PtIndsByMovSliceAndContour[i] for i in range(2)])


# In[186]:


print(MovContourNosThisSlice)

UniqueContourNos = list(set(MovContourNosThisSlice))

print(UniqueContourNos)

IndsByContour = [[i for i, e in enumerate(MovContourNosThisSlice) if e==c] for c in UniqueContourNos]

print(IndsByContour)

PtIndsThisSliceByContour = [[PtIndsThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]
MovContourNosThisSliceByContour = [[MovContourNosThisSlice[i] for i in IndsByContour[j]] for j in range(len(IndsByContour))]

#print(PtIndsThisSliceByContour)
print(MovContourNosThisSliceByContour)


# In[171]:


[len(MovContourNosByMovSlice[i]) for i in range(2)]


# In[173]:


Ns = [len(MovContourNosByMovSlice[i]) for i in range(2)]

#N = 0
#N = [N + Ns[i] for i in range(len(Ns))]
#N

sum(Ns)


# In[149]:


import GetOutputPoints
importlib.reload(GetOutputPoints)
from GetOutputPoints import GetOutputPoints

PtNos, InInds, FixOutInds,MovPtsPCS, MovPtsICS,DefPCS, DefICS, MovOutInds,MovPtsBySlicePCS, MovPtsBySliceICS = GetOutputPoints(InPts_ICS=FixPtsICS, Origin=MovOrigin, 
                                                     Directions=MovDirs, Spacings=MovSpacings, 
                                                     Dimensions=MovDims, LUT=LUT)


# In[154]:


LUT


# In[151]:


list(np.random.choice(range(256), size=3))


# In[152]:


list(np.random.random(size=3) * 256)


# In[213]:


# Test latest refactored functions:

import GetInputPoints
importlib.reload(GetInputPoints)
from GetInputPoints import GetInputPoints

FixPtsPCS, FixPtsBySliceAndContourPCS,FixPtsICS, FixPtsBySliceAndContourICS,LUT = GetInputPoints(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath,
                     Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings)

import GetOutputPoints
importlib.reload(GetOutputPoints)
from GetOutputPoints import GetOutputPoints

PtNos, InInds, FixOutInds,MovPtsPCS, MovPtsICS,DefPCS, DefICS, MovOutInds,MovPtsBySliceAndContourPCS,MovPtsBySliceAndContourICS,LUT = GetOutputPoints(FixPts_ICS=FixPtsICS, MovOrigin=MovOrigin, 
                      MovDirs=MovDirs, MovSpacings=MovSpacings, 
                      MovDims=MovDims, LUT=LUT)


# In[217]:


LUT

# unique set of slice nos (z-indeces) for moving image:
set([LUT[i][7][2] for i in range(len(LUT))])


# In[214]:


MovPtsBySliceAndContourPCS


# ### July 2: Start again with the refactored code:

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

#from GetDicomFpaths import GetDicomFpaths
#from GetDicoms import GetDicoms
#import PCStoICS
#importlib.reload(PCStoICS)
#from PCStoICS import PCStoICS
#import GetContourPtsForDicoms
#importlib.reload(GetContourPtsForDicoms)
#from GetContourPtsForDicoms import GetContourPtsForDicoms
#from ParseTransformixOutput import ParseTransformixOutput




# Define FixedDicomDir and MovingDicomDir:
DicomDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011'
FixDicomDir = os.path.join(DicomDir, r'04-10-1960-MRI Brain wwo Contrast-69626\8-T1 SE AXIAL POST FS FC-59362')
MovDicomDir = os.path.join(DicomDir, r'06-11-1961-MRI Brain wwo Contrast-79433\8-T1 SE AXIAL POST FS FC-81428')

# Define the filepath to the ROI Collection (for the fixed) image:
FixRoiDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12\8 T1 SE AXIAL POST FS FC ROIs (source)'
FixRoiFname = r'AIM_20200511_073405.dcm'
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

# Create inputpoints.txt for Elastix, containing the contour points to be
# transformed:
#CreateInputFileForElastix(DicomDir=FixDicomDir, RoiFpath=FixRoiFpath, 
#                          Origin=FixOrigin, Directions=FixDirs, Spacings=FixSpacings,
#                          CoordSys=CoordSys)
CreateInputFileForElastix(Points=FixPtsPCS)


# Transform MovingContourPts:
TransformPoints(MovingImage=MovIm, TransformFilter=ElastixImFilt)

# Parse outputpoints.txt:
#PtNos, InInds, InPts_PCS, FixOutInds,\
#OutPts_PCS, Defs_PCS, MovOutInds = ParseTransformixOutput()

print('Parsing outputpoints.txt and grouping transformed points into slices and contours...')

PtNos, InInds, FixOutInds,MovPtsPCS, MovPtsICS,DefPCS, DefICS, MovOutInds,MovPtsBySliceAndContourPCS,MovPtsBySliceAndContourICS,LUT = GetOutputPoints(FixPts_ICS=FixPtsICS, MovOrigin=MovOrigin, 
                      MovDirs=MovDirs, MovSpacings=MovSpacings, 
                      MovDims=MovDims, LUT=LUT)

print('Done.')


# In[4]:


import RegUtilityFuncs
importlib.reload(RegUtilityFuncs)
import RegUtilityFuncs as ruf

# Chose which slices to plot:
plot_slices = -1 # plot all
plot_slices = [14]
plot_slices = [12, 13, 14, 15]
plot_slices = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]


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
fig_fname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_' + CoordSys             +'_reg_result_using_' + package + '_' + perspective + '.png'

#ruf.display_all_sitk_images_and_reg_results_with_all_contours_v2(fix_im=FixIm, mov_im=MovIm, reg_im=RegIm,
#                                                                 fix_pts=InPtsArr_ICS, mov_pts=OutPtsArr_ICS,
#                                                                 export_fig=export_fig,
#                                                                 export_fname=fig_fname,
#                                                                 plot_slices=plot_slices,
#                                                                 perspective=perspective)

ruf.display_all_sitk_images_and_reg_results_with_all_contours_v3(fix_im=FixIm, 
                                                                 mov_im=MovIm, 
                                                                 reg_im=RegIm,
                                                                 fix_pts=FixPtsBySliceAndContourICS, 
                                                                 mov_pts=MovPtsBySliceAndContourICS,
                                                                 export_fig=export_fig,
                                                                 export_fname=fig_fname,
                                                                 plot_slices=plot_slices,
                                                                 perspective=perspective,
                                                                 contours_as=contours_as)


# In[13]:


FixPtsBySliceAndContourICS


# In[18]:


print(FixPtsBySliceAndContourICS[0])

print(len(FixPtsBySliceAndContourICS[0]))


# In[16]:


fix_Ncontours = 0

for s in FixPtsBySliceAndContourICS:
    print(s)
    #print(len(FixPtsBySliceAndContourICS[s]))
    
    #fix_Ncontours = fix_Ncontours + len(FixPtsBySliceAndContourICS[s])


# In[10]:


list(np.random.choice(range(1), size=3))


# In[7]:


list(np.random.choice(range(2), size=3))


# In[11]:


[item/255 for item in list(np.random.choice(range(256), size=3))]


# In[ ]:




