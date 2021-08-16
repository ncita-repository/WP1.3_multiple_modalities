#!/usr/bin/env python
# coding: utf-8

# #### This notebook demonstrates the ability to copy a contour from one slice to another for the following 4 use cases:
# 
# #### 1. Copy a contour across slices in the same session/scan.
# 
# #### 2. Copy a contour across sessions/scans in the same study with the same/similar* image resolutions.
# 
# #### 3. Copy a contour across sessions/scans in the same study with different image resolutions and/or slice thickness.
# 
# #### 4. Copy a contour across studies (i.e. with different Frame of References).
# 
# ##### * Similar since for the ROI Collection I am working with there are no other scans with identical resolutions. 

# In[21]:


import importlib
import PlotDicomsAndRois_v2
importlib.reload(PlotDicomsAndRois_v2)
from PlotDicomsAndRois_v2 import PlotDicomsAndRois_v2 as PlotStacks
import os
import time
get_ipython().run_line_magic('matplotlib', 'inline')
#%matplotlib notebook


# ### Define the directories and filepaths:

# In[88]:


DicomRootDir = r'C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011\04-10-1960-MRI Brain wwo Contrast-69626'

S011_MR4_S8_DicomDir = os.path.join(DicomRootDir, r'8-T1 SE AXIAL POST FS FC-59362')
S011_MR4_S5_DicomDir = os.path.join(DicomRootDir, r'5-T1 SE AXIAL 3MM-81246')
S011_MR4_S4_DicomDir = os.path.join(DicomRootDir, r'4-T2 AXIAL FLAIR DARK FL-48830')

RoiRootDir = r'C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011'

S011_MR4_S8_RoiDir = os.path.join(RoiRootDir, 'Fixed_ROI_Collections')
S011_MR4_S8_RoiFpath = os.path.join(S011_MR4_S8_RoiDir, 'AIM_20200511_073405.dcm')

TargetRoiDir = os.path.join(RoiRootDir, 'New_ROI_Collections')


# In[91]:


#SourceRoi


# ### 1. Copy the contour on slice 15 to slice 19 within the same scan (series 8):

# In[85]:


import CopyRois_v3
importlib.reload(CopyRois_v3)
from CopyRois_v3 import CopyRois_v3

# Choose the slice numbers containing the contour to be copied from
# Source to Target:
"""
Note: CopySliceNums must be a list of integers - even if it contains
only a single integer.  For now only dealing with a list of length 1. 
"""
#CopyFromSliceNum = 0
CopyFromSliceNum = 15

CopyToSliceNum = 19

# Chose whether or not to export the Target ROI:
ExportRoi = False
#ExportRoi = True
 
LogToConsole = False
#LogToConsole = True

# Choose Structure Set Label (appears under Name in XNAT!):
TargetStructureSetLabel = f'Copy_from_MR4_S8_s{CopyFromSliceNum}_to_s{CopyToSliceNum}'

TargetRoi_no1 = CopyRois_v3(SourceDicomDir=S011_MR4_S8_DicomDir, 
                            SourceRoiFpath=S011_MR4_S8_RoiFpath, 
                            TargetDicomDir=S011_MR4_S8_DicomDir, 
                            TargetRoiDir=TargetRoiDir,
                            CopyFromSliceNum=CopyFromSliceNum, 
                            CopyToSliceNum=CopyToSliceNum,
                            TargetStructureSetLabel=TargetStructureSetLabel, 
                            ExportRoi=ExportRoi,
                            LogToConsole=LogToConsole)


# ### 2. Copy the contour on slice 15 from session 8 to slice 25 in session 5:

# In[87]:


import CopyRois_v3
importlib.reload(CopyRois_v3)
from CopyRois_v3 import CopyRois_v3

# Choose the slice numbers containing the contour to be copied from
# Source to Target:
"""
Note: CopySliceNums must be a list of integers - even if it contains
only a single integer.  For now only dealing with a list of length 1. 
"""
CopyFromSliceNum = 0
CopyFromSliceNum = 15

CopyToSliceNum = 25

# Chose whether or not to export the Target ROI:
ExportRoi = False
#ExportRoi = True
 
LogToConsole = False
#LogToConsole = True

# Choose Structure Set Label (appears under Name in XNAT!):
TargetStructureSetLabel = f'Copy_from_MR4_S8_s{CopyFromSliceNum}_to_S5_s{CopyToSliceNum}'


TargetRoi_no2 = CopyRois_v3(SourceDicomDir=S011_MR4_S8_DicomDir, 
                            SourceRoiFpath=S011_MR4_S8_RoiFpath, 
                            TargetDicomDir=S011_MR4_S5_DicomDir, 
                            TargetRoiDir=TargetRoiDir,
                            CopyFromSliceNum=CopyFromSliceNum, 
                            CopyToSliceNum=CopyToSliceNum,
                            TargetStructureSetLabel=TargetStructureSetLabel, 
                            ExportRoi=ExportRoi,
                            LogToConsole=LogToConsole)


# ### 3. Copy the contour on slice 15 from session 8 to slice 25 in session 4:

# In[90]:


import CopyRois_v3
importlib.reload(CopyRois_v3)
from CopyRois_v3 import CopyRois_v3

# Choose the slice numbers containing the contour to be copied from
# Source to Target:
"""
Note: CopySliceNums must be a list of integers - even if it contains
only a single integer.  For now only dealing with a list of length 1. 
"""
CopyFromSliceNum = 0
CopyFromSliceNum = 15

CopyToSliceNum = 10

# Chose whether or not to export the Target ROI:
ExportRoi = False
#ExportRoi = True
 
LogToConsole = False
#LogToConsole = True

# Choose Structure Set Label (appears under Name in XNAT!):
TargetStructureSetLabel = f'Copy_from_MR4_S8_s{CopyFromSliceNum}_to_S4_s{CopyToSliceNum}'


TargetRoi_no3 = CopyRois_v3(SourceDicomDir=S011_MR4_S8_DicomDir, 
                            SourceRoiFpath=S011_MR4_S8_RoiFpath, 
                            TargetDicomDir=S011_MR4_S4_DicomDir, 
                            TargetRoiDir=TargetRoiDir,
                            CopyFromSliceNum=CopyFromSliceNum, 
                            CopyToSliceNum=CopyToSliceNum,
                            TargetStructureSetLabel=TargetStructureSetLabel, 
                            ExportRoi=ExportRoi,
                            LogToConsole=LogToConsole)


# In[ ]:





# ### Import DICOMs and ROIs for plots:

# In[97]:


import pydicom
from GetDicoms import GetDicoms

# The DICOMs containing the original ROIs (Use Case #1):
S011_MR4_S8_DicomFpaths,S011_MR4_S8_Dicoms = GetDicoms(DirPath=S011_MR4_S8_DicomDir, 
                               SortMethod='slices', Debug=False)

# The original ROIs:
S011_MR4_S8_Roi = pydicom.dcmread(S011_MR4_S8_RoiFpath)

# The new ROI created for Use Case #1:
S011_MR4_S8_NewRoiFpath = os.path.join(TargetRoiDir, 
                                       'AIM_20200918_195841_copied_from_AIM_20200511_073405.dcm')

# The set of DICOMs for Use Case #2:
S011_MR4_S5_DicomFpaths,S011_MR4_S5_Dicoms = GetDicoms(DirPath=S011_MR4_S5_DicomDir, 
                               SortMethod='slices', Debug=False)

# The new ROI created for Use Case #2:
S011_MR4_S5_RoiFpath = os.path.join(TargetRoiDir, 
                                    'AIM_20200918_200347_copied_from_AIM_20200511_073405.dcm')

# The set of DICOMs for Use Case #3:
S011_MR4_S4_DicomFpaths,S011_MR4_S4_Dicoms = GetDicoms(DirPath=S011_MR4_S4_DicomDir, 
                               SortMethod='slices', Debug=False)

# The new ROI created for Use Case #3:
S011_MR4_S4_RoiFpath = os.path.join(TargetRoiDir, 
                                    'AIM_20200918_202506_copied_from_AIM_20200511_073405.dcm')


# Define location for plots to be exported to:
ExportDir = os.path.join(RoiRootDir, 'Results_of_Copy_ROI_use_cases')


# ### Plot the DICOMs for Subject 011 MR 4 Session 8 and the (original) ROIs:

# In[95]:


""" Choose plotting parameters: """

# Decide on which DICOMs and ROIs to plot:
""" Set RoiFpath = None if no ROIs to plot. """
DicomDir = S011_MR4_S8_DicomDir; DicomsToPlot = 'S011_MR4_S8'
RoiFpath = S011_MR4_S8_RoiFpath; RoisToPlot = 'S011_MR4_S8'

# Decide on SlicesToPlot:
"""
Note:  SlicesToPlot must either be -1, -2 or a list of integers.  

If SlicesToPlot is:
          -1 --> All slices will be plotted.
          -2 --> All slices containing contours will be plotted if 
                 RoiFpath is not None; All slices will be plotted if 
                 RoiFpath is None.
[13, 14, 15] --> Slices 13, 14 and 15 will be plotted.
        [13] --> Slice 13 will be plotted. Note that even a single 
                 integer must be enclosed in brackets.
"""
SlicesToPlot = -1 # all slices
SlicesToPlot = -2 # all slices containing contours
#SlicesToPlot = [13, 14, 15]
#SlicesToPlot = [13]

LogToConsole = False
#LogToConsole = True

# Decide whether to export plot:
ExportPlot = False
#ExportPlot = True

# Create filename for exported figure:
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_'       + DicomsToPlot + '_DICOMs_and_' + RoisToPlot + '_ROIs.png'

ExportFpath = os.path.join(ExportDir, ExportFname)


PlotStacks(DicomDir=DicomDir, RoiFpath=RoiFpath, 
           SlicesToPlot=SlicesToPlot, LogToConsole=LogToConsole, 
           ExportPlot=ExportPlot, ExportFname=ExportFpath)


# ### Plot the DICOMs for Subject 011 MR 4 Session 8 and the ROIs for Use Case #1:

# In[100]:


""" Choose plotting parameters: """

# Decide on which DICOMs and ROIs to plot:
""" Set RoiFpath = None if no ROIs to plot. """
DicomDir = S011_MR4_S8_DicomDir; DicomsToPlot = 'S011_MR4_S8'
RoiFpath = S011_MR4_S8_NewRoiFpath; RoisToPlot = 'S011_MR4_S8_New'

# Decide on SlicesToPlot:
"""
Note:  SlicesToPlot must either be -1, -2 or a list of integers.  

If SlicesToPlot is:
          -1 --> All slices will be plotted.
          -2 --> All slices containing contours will be plotted if 
                 RoiFpath is not None; All slices will be plotted if 
                 RoiFpath is None.
[13, 14, 15] --> Slices 13, 14 and 15 will be plotted.
        [13] --> Slice 13 will be plotted. Note that even a single 
                 integer must be enclosed in brackets.
"""
SlicesToPlot = -1 # all slices
SlicesToPlot = -2 # all slices containing contours
#SlicesToPlot = [13, 14, 15]
#SlicesToPlot = [13]

LogToConsole = False
#LogToConsole = True

# Decide whether to export plot:
ExportPlot = False
#ExportPlot = True

# Create filename for exported figure:
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_'       + DicomsToPlot + '_DICOMs_and_' + RoisToPlot + '_ROIs.png'

ExportFpath = os.path.join(ExportDir, ExportFname)


PlotStacks(DicomDir=DicomDir, RoiFpath=RoiFpath, 
           SlicesToPlot=SlicesToPlot, LogToConsole=LogToConsole, 
           ExportPlot=ExportPlot, ExportFname=ExportFpath)


# ### Plot the DICOMs for Subject 011 MR 4 Session 5 and the ROIs for Use Case #2:

# In[102]:


""" Choose plotting parameters: """

# Decide on which DICOMs and ROIs to plot:
""" Set RoiFpath = None if no ROIs to plot. """
DicomDir = S011_MR4_S5_DicomDir; DicomsToPlot = 'S011_MR4_S5'
RoiFpath = S011_MR4_S5_RoiFpath; RoisToPlot = 'S011_MR4_S5'

# Decide on SlicesToPlot:
"""
Note:  SlicesToPlot must either be -1, -2 or a list of integers.  

If SlicesToPlot is:
          -1 --> All slices will be plotted.
          -2 --> All slices containing contours will be plotted if 
                 RoiFpath is not None; All slices will be plotted if 
                 RoiFpath is None.
[13, 14, 15] --> Slices 13, 14 and 15 will be plotted.
        [13] --> Slice 13 will be plotted. Note that even a single 
                 integer must be enclosed in brackets.
"""
SlicesToPlot = -1 # all slices
SlicesToPlot = -2 # all slices containing contours
#SlicesToPlot = [13, 14, 15]
#SlicesToPlot = [13]

LogToConsole = False
#LogToConsole = True

# Decide whether to export plot:
ExportPlot = False
#ExportPlot = True

# Create filename for exported figure:
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_'       + DicomsToPlot + '_DICOMs_and_' + RoisToPlot + '_ROIs.png'

ExportFpath = os.path.join(ExportDir, ExportFname)


PlotStacks(DicomDir=DicomDir, RoiFpath=RoiFpath, 
           SlicesToPlot=SlicesToPlot, LogToConsole=LogToConsole, 
           ExportPlot=ExportPlot, ExportFname=ExportFpath)


# ### Plot the DICOMs for Subject 011 MR 4 Session 4 and the ROIs for Use Case #3:

# In[104]:


""" Choose plotting parameters: """

# Decide on which DICOMs and ROIs to plot:
""" Set RoiFpath = None if no ROIs to plot. """
DicomDir = S011_MR4_S4_DicomDir; DicomsToPlot = 'S011_MR4_S4'
RoiFpath = S011_MR4_S4_RoiFpath; RoisToPlot = 'S011_MR4_S4'

# Decide on SlicesToPlot:
"""
Note:  SlicesToPlot must either be -1, -2 or a list of integers.  

If SlicesToPlot is:
          -1 --> All slices will be plotted.
          -2 --> All slices containing contours will be plotted if 
                 RoiFpath is not None; All slices will be plotted if 
                 RoiFpath is None.
[13, 14, 15] --> Slices 13, 14 and 15 will be plotted.
        [13] --> Slice 13 will be plotted. Note that even a single 
                 integer must be enclosed in brackets.
"""
SlicesToPlot = -1 # all slices
SlicesToPlot = -2 # all slices containing contours
#SlicesToPlot = [13, 14, 15]
#SlicesToPlot = [13]

LogToConsole = False
#LogToConsole = True

# Decide whether to export plot:
ExportPlot = False
#ExportPlot = True

# Create filename for exported figure:
ExportFname = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + '_'       + DicomsToPlot + '_DICOMs_and_' + RoisToPlot + '_ROIs.png'

ExportFpath = os.path.join(ExportDir, ExportFname)


PlotStacks(DicomDir=DicomDir, RoiFpath=RoiFpath, 
           SlicesToPlot=SlicesToPlot, LogToConsole=LogToConsole, 
           ExportPlot=ExportPlot, ExportFname=ExportFpath)


# In[ ]:





# In[ ]:




