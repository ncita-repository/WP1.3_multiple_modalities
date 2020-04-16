#!/usr/bin/env python
# coding: utf-8

# # Trying to remotely access DICOM tags, pixel arrays from DICOMs, and ROIs from DICOM-RTSTRUCT files (without downloading the files locally first)

# In[1]:


import xnat
#xnat.__version__ # 0.3.21
print(f'XNATpy version = {xnat.__version__}') # 0.3.21
import SimpleITK as sitk
import pydicom
import dicom2nifti
import numpy as np
#import math
import matplotlib.pyplot as plt
#from matplotlib import pyplot
import matplotlib.image as mpimg
import os, sys, time
import natsort
import dicom_functions
import explore_dicom_data
from myshow import myshow
import importlib


# # Connect to XNAT server:

# In[2]:


xnatAddress = 'http://10.1.1.17'

session = xnat.connect(xnatAddress, user='admin', password='admin')
#session = xnat.connect(xnatAddress, user='admin', password='admin', debug=True)


# # Define the REST path variables:

# In[3]:


# Define selected project and subject labels, and session number:
projectLabel = 'BrainTumorProg'
subjectLabel = 'PGM-002'
experimentNo = 0 # (= 1 in XNAT web app; e.g. ‘PGM-002_MR_1’)
scanNo = 0 # 0 -> T1post


# # Try to access DICOM tags and pixel arrays without downloading the files:

# In[4]:


scans = session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo]

scans


# In[5]:


# This reads only the header and not the pixel data:
# (see:
# https://xnat.readthedocs.io/en/latest/static/tutorial.html#accessing-xnat-files-as-local-files-partial-read)
# )
scans.read_dicom()


# In[6]:


# To read the pixel data:
scans.read_dicom(read_pixel_data=True)


# In[7]:


temp = scans.read_dicom(read_pixel_data=True).PixelData
len(temp)


# Looks like the pixel data is flat! 
# 
# Maybe pixel data is flattened like contour data is, so I'll try replacing "read_pixel_data=True" with
# "read_pixel_array=True":

# In[8]:


# Try this:
scans.read_dicom(read_pixel_array=True)


# Maybe the clue is here:
# 
# https://xnat.readthedocs.io/en/latest/static/tutorial.html#exploring-your-xnat-server
# 
# under "Accessing XNAT files as local files (partial read)"

# In[9]:


session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo]


# In[10]:


session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo].resources['DICOM']


# In[11]:


session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo].resources['DICOM'].files


# In[12]:


session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo].resources['DICOM'].files[0]


# In[13]:


fileObj = session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo].resources['DICOM'].files[0]

with fileObj.open() as fileObjOpen:
    fileObjOpen.read()
    
fileObjOpen


# In[14]:


len(fileObjOpen)


# In[15]:


fileObj = session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].scans[scanNo].resources['DICOM'].files[0]

with fileObj.open() as fileObjOpen:
    temp = fileObjOpen.read(3000)
    
temp


# In[16]:


len(temp)


# Ran out of ideas for now..
