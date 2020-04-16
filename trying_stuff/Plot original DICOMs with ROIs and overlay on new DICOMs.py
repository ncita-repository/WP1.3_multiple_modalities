#!/usr/bin/env python
# coding: utf-8

# In[12]:


import importlib

import PlotDICOMsAndROIs
importlib.reload(PlotDICOMsAndROIs)
from PlotDICOMsAndROIs import PlotDICOMsAndROIs

import ModifyROIs
importlib.reload(ModifyROIs)
from ModifyROIs import ModifyROIs

import os


# # Define the directories and filepaths:

# In[15]:


# Define the directories and filepaths:

origDicomsDir = r'C:\Code\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\PGM-002_MR_1\T1post\PGM-002_MR_1\scans\11-T1post\resources\DICOM\files'

origRoiFpath = r'C:\Code\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\PGM-002_MR_1\T1post\AIM_20200123_145716\out\resources\RTSTRUCT\files\AIM_20200123_145716.dcm'

newDicomsDir = r'C:\Code\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\PGM-002_MR_1\T1pre_reg\PGM-002_MR_1\scans\19446-T1pre_reg\resources\DICOM\files'

newRoiDir = r'C:\Code\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\PGM-002_MR_1\T1post\AIM_20200123_145716\out\resources\RTSTRUCT\files'

# Define location for plots to be exported to:
exportDir = r'C:\Code\WP1.3_multiple_modalities\trying stuff\plots'

# Decide whether to export plot:
exportPlot = 0 # no
exportPlot = 1 # yes

# Print paths:
print('Directory of original series of DICOMs:\n\n', origDicomsDir, '\n')
print('Filepath of original DICOM-RTSTRUCT file (for the original series of DICOMs):\n\n', origRoiFpath, '\n')
print('Directory of new series of DICOMs:\n\n', newDicomsDir, '\n')
print('Directory where the new DICOM-RTSTRUCT file (for the new series of DICOMs) will be exported to:\n\n',       origRoiFpath, '\n')
if exportPlot:
    print('Plot will be exported to exportDir:\n\n', exportDir)
else:
    print('Plots will not be exported.')


# # First plot the original DICOM series and their corresponding ROIs:

# In[16]:


# First plot the original DICOM series and their corresponding ROIs:
PlotDICOMsAndROIs(origDicomsDir, origRoiFpath, exportPlot, exportDir)


# # Modify the original DICOM-RTSTUCT file for use with the new series of DICOMs:

# In[13]:


# Modify the original DICOM-RTSTUCT file for use with the new series of DICOMs:
newRoiFpath = ModifyROIs(origDicomsDir, origRoiFpath, newDicomsDir, newRoiDir)

print('Filepath of new DICOM-RTSTRUCT file (the ROIs from the original DICOMs adapted for the new DICOMs):\n\n',      newRoiFpath)


# # Plot the new DICOM series with the ROIs from the original DICOM series:

# In[17]:


# Plot the new DICOM series with the ROIs from the original DICOM series:
PlotDICOMsAndROIs(newDicomsDir, newRoiFpath, exportPlot, exportDir)


# In[ ]:




