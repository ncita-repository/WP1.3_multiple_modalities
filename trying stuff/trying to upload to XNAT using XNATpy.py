#!/usr/bin/env python
# coding: utf-8

# # Trying to upload DICOMs and ROIs (DICOM-RTSTRUCTs) to XNAT
# 
# ### https://xnat.readthedocs.io/en/latest/index.html

# In[2]:


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

# In[3]:


xnatAddress = 'http://10.1.1.17'
#xnatAddress = 'http://10.1.1.17/app/template/Index.vm'
#xnatAddress = 'http://10.1.1.17/app/template/Login.vm#!'
#xnatAddress = 'http://10.1.1.17:8080/xnat'
#xnatAddress = 'http://10.1.1.17:80/xnat'
#xnatAddress = 'http://localhost:8080/xnat'
#xnatAddress = 'http://10.1.1.17:8080/admin'
#xnatAddress = 'http://localhost:8080/admin'

#session = xnat.connect(xnatAddress, user='xnat', password='xnat')
#session = xnat.connect(xnatAddress, user='owner', password='owner')
session = xnat.connect(xnatAddress, user='admin', password='admin', debug=True)


# In[4]:


session


# # Define the REST path variables:

# In[5]:


# Define selected project and subject labels, and session number:
projectLabel = 'BrainTumorProg'
subjectLabel = 'PGM-002'
experimentNo = 0 # (= 1 in XNAT web app; e.g. ‘PGM-002_MR_1’)
scanNo = 0 # 0 -> T1post


# # Location of previously downloaded DICOM and DICOM-RTSTRUCT files:

# In[6]:


# The DICOM files were exported here:
dicomsPath = r'C:\Users\ctorti\Documents\GitHub\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\0\0\PGM-002_MR_1\scans\11-T1post\resources\DICOM\files'

# And the ROI DICOM-RTStruct file was exported here:
roiPath = r'C:\Users\ctorti\Documents\GitHub\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\0\0\AIM_20200123_145716\out\resources\RTSTRUCT\files'


# In[ ]:





# # Try uploading data to XNAT

# In[243]:


session


# In[244]:


dir(session.classes)


# Recall that the REST parameters were defined as such:
#     
# Define selected project and subject labels, and session number:
# projectLabel = 'BrainTumorProg'
# subjectLabel = 'PGM-002'
# experimentNo = 0 # (= 1 in XNAT web app; e.g. ‘PGM-002_MR_1’)
# scanNo = 0 # 0 -> T1post

# Recall that the DICOM-RTStruct file was downloaded to newPath using the command:
# ROIobject = session.projects[projectLabel].subjects[subjectLabel].experiments[experimentNo].assessors[scanNo].resources[1]
# ROIobject.download_dir(newPath)
# 
# So I'll have to (at least) create a new project, then new subject, then new experiment before uploading the ROI file.
# Not sure if I need to create an assessor and resource as well, or if it's automatically handled...

# In[240]:


# And newPath is:
newPath


# In[7]:


# The full filepath to the downloaded DICOM-RTStruct file is:
roiFpathDL = r'C:\Code\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\0\0\AIM_20200123_145716\out\resources\RTSTRUCT\files\AIM_20200123_145716.dcm' 


# In[8]:


# Define the new project name:
newProjectName = projectLabel + '_copy'

newProjectName


# The XNATpy guide's instructions on how to import DICOM files in an archived (zipped) file via prearchive:
# 
# prearchive_session = session.services.import_('/home/hachterberg/temp/ANONYMIZ.zip', \
#                                               project='brainimages', \
#                                               destination='/prearchive')
# 
# but I'm not sure if the project being imported to must already exist, or if one will be created if not...

# Following the guidance uner "Object creation" here:
# 
# https://xnat.readthedocs.io/en/latest/static/tutorial.html#exploring-your-xnat-server
# 
# A new subject is created as follows:
# 
# connection = xnat.connect('https://xnat.example.com')
# project = connection.projects['myproject']
# subject = connection.classes.SubjectData(parent=project, label='new_subject_label')
# 
# The parent of a subject would be a project, but I'm not sure what the parent of a project is...
# 
# Perhaps the session?
# 
# Yes found an example here:
# 
# https://groups.google.com/forum/#!msg/xnat_discussion/vYJuSZ6MbmI/dTnKfU5_BQAJ

# # Create a new project in XNAT:

# In[9]:


newProject = session.classes.ProjectData(parent=session, name=newProjectName)


# In[10]:


newProject


# It worked!  And I confirmed in XNAT web interface.

# # Create a new subject for the new project:

# In[11]:


newSubject = session.classes.SubjectData(label='PGM-002_copy', parent=newProject)


# In[12]:


newSubject


# It worked!  And I confirmed in XNAT web interface.

# # Create a new experiment for the new subject:

# In[313]:


newExperiment = session.classes.MrSessionData(label='MR Session', parent=newSubject)


# Not sure why this hasn't worked.  I was able to create a CR experiment and a MR experiment in XNAT web interface...

# In[17]:


# Charlie Moore has replied to my question here:

# https://groups.google.com/forum/#!msg/xnat_discussion/vYJuSZ6MbmI/dTnKfU5_BQAJ

# saying that "I believe XNAT prohibits use of the space character in experiment labels. 
# Could you try this again without the space in your “label” value?"

newExperiment = session.classes.MrSessionData(label='MR_Session', parent=newSubject)


# In[18]:


newExperiment


# # SUCCESS!  Seems I was confusing "Label" with "Experiment".  I thought I was supposed to assign the label "MR Session" in the same way that in the web interface I would need to select the type of experiment, e.g. "MR Session", "CT Session", etc.

# # Need to continue from here now that I've managed to create an experiment without using the WI.  Need to make sense of what I've done so far..

# In[ ]:





# In[ ]:





# # Using the web interface I created a new experiment called "PGM-002_copy_MR_2" selecting "MR Session" as type.  

# Then I found this thread:
# https://groups.google.com/forum/#!topic/xnat_discussion/dmM9_cwbn38
# 
# So will try using .classes.ResourceCatalog():

# In[283]:


newExperiment = session.classes.ResourceCatalog(parent=newSubject, label='PGM-002_copy_MR_1')


# In[284]:


newExperiment


# I didn't get an error, so I thought it worked.  But I didn't see any change in the web interface.  Then I realised that it was probably not correct since I need to create an MR Session.  Trying label='MR Session' this time:

# In[285]:


newExperiment = session.classes.ResourceCatalog(parent=newSubject, label='MR Session')


# In[314]:


newExperiment


# Although I didn't get an error when I re-defined newExperiment I did when I evaluated the variable itself.  And I also didn't see a new MR Session in the web interface.  So not sure what I've done.

# # Returning to the MR Session that I created using the web interface (WI), I'll see if I can upload an assessor to it.

# In[292]:


session.projects


# In[294]:


session.projects[newProjectName]


# In[13]:


newSubjectName = 'PGM-002_copy'
newExperimentName1 = 'PGM-002_copy_MR_1' # the name of the exp I tried to create in XNATpy
newExperimentName2 = 'PGM-002_copy_MR_2' # the name of the exp I create in the WI

session.projects[newProjectName].subjects[newSubjectName]


# In[14]:


newExperiment2 = session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName2]
newExperiment2


# In[304]:


newExperiment2.scans


# In[305]:


newExperiment2.accessors


# In[306]:


newExperiment2.resources


# Create an assessor following Hakim's guidance here:
# https://groups.google.com/forum/#!searchin/xnat_discussion/xnatpy$20status$20500|sort:date/xnat_discussion/N7i2tS5MleE/7Sxlr7JIBAAJ

# In[15]:


assessorLabel2 = 'ROI_tumour'
newExperiment2.create_assessor(assessorLabel2, type_='xnat:qcAssessmentData')  # Note the xsi_type here, this can be something else if you need it to be


# This has created an Assessment with the Experiment label "Auto QC" in the WI.

# In[16]:


# Now see if I can add files:

newExperiment2.create_assessor(assessorLabel2, type_='xnat:qcAssessmentData', files=roiFpathDL)  # Note the xsi_type here, this can be something else if you need it to be


# Clearly not like that.

# In[308]:


assessorLabel2 = 'ROI_tumour'
newExperiment2.create_assessor(assessorLabel2, type_='xnat:qcAssessmentData', label=assessorLabel2)


# So can't add a label using .create_assessor().  Next will try the "Alternative strategy (more generic)" provided by Hakim:

# In[309]:


assessor = session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName2]
assessor = session.classes.QcAssessmentData(parent=newExperiment2, label=assessorLabel2)
resource = session.classes.ResourceCatalog(parent=assessor, label='RESOURCE_LABEL') # was connection.classes... 
# but might have been a typo?


# In[310]:


assessor


# In[311]:


resource


# I'm not sure what this has done.  The WI looks the same other than perhaps a "+ Session QC Status" above Scans that may or may not have been present before.
# 
# So now need to figure out how to upload the DICOM-RTSTRUCT file to assessor.  Try importing via the prearchive:

# In[1]:


prearchive_session = session.services.import_(roiFpathDL,                                               project=newProjectName,                                               destination='/prearchive')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




