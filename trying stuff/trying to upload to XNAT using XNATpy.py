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


# #### Connect to XNAT server:

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


# #### Define the REST path variables:

# In[5]:


# Define selected project and subject labels, and session number:
projectLabel = 'BrainTumorProg'
subjectLabel = 'PGM-002'
experimentNo = 0 # (= 1 in XNAT web app; e.g. ‘PGM-002_MR_1’)
scanNo = 0 # 0 -> T1post


# #### Location of previously downloaded DICOM and DICOM-RTSTRUCT files:

# In[6]:


# The DICOM files were exported here:
dicomsPath = r'C:\Users\ctorti\Documents\GitHub\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\0\0\PGM-002_MR_1\scans\11-T1post\resources\DICOM\files'

# And the ROI DICOM-RTStruct file was exported here:
roiPath = r'C:\Users\ctorti\Documents\GitHub\WP1.3_multiple_modalities\trying stuff\XNAT downloads\BrainTumorProg\PGM-002\0\0\AIM_20200123_145716\out\resources\RTSTRUCT\files'


# In[ ]:





# #### Try uploading data to XNAT

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

# #### Create a new project in XNAT:

# In[9]:


newProject = session.classes.ProjectData(parent=session, name=newProjectName)


# In[10]:


newProject


# It worked!  And I confirmed in XNAT web interface.

# #### Create a new subject for the new project:

# In[11]:


newSubjectName = 'PGM-002_copy'

newSubject = session.classes.SubjectData(label='PGM-002_copy', parent=newProject)


# In[12]:


newSubject


# It worked!  And I confirmed in XNAT web interface.

# #### Create a new experiment for the new subject:
# 
# ##### Note:  Spaces are not permitted for labels!

# In[55]:


# Charlie Moore has replied to my question here:

# https://groups.google.com/forum/#!msg/xnat_discussion/vYJuSZ6MbmI/dTnKfU5_BQAJ

# saying that "I believe XNAT prohibits use of the space character in experiment labels. 
# Could you try this again without the space in your “label” value?"

newExperimentName = 'PGM-002_copy_MR_2'

newExperiment = session.classes.MrSessionData(label=newExperimentName, parent=newSubject)


# In[56]:


newExperiment


# #### SUCCESS!  Seems I was confusing "Label" with "Experiment".  I thought I was supposed to assign the label "MR Session" in the same way that in the web interface I would need to select the type of experiment, e.g. "MR Session", "CT Session", etc.

# #### Need to continue from here now that I've managed to create an experiment without using the WI.  Need to make sense of what I've done so far..

# In[ ]:





# In[ ]:





# Then I found this thread:
# https://groups.google.com/forum/#!topic/xnat_discussion/dmM9_cwbn38
# 
# So will try using .classes.ResourceCatalog():

# In[21]:


newExperiment2 = session.classes.ResourceCatalog(parent=newSubject, label='PGM-002_copy_MR_2_2')


# In[22]:


newExperiment2


# I didn't get an error, so I thought it worked.  But I didn't see any change in the web interface.  So not sure what it has done.

# #### Returning to the MR Session that I created above, I'll see if I can upload an assessor to it.

# In[23]:


session.projects


# In[24]:


session.projects[newProjectName]


# In[25]:


session.projects[newProjectName].subjects[newSubjectName]


# In[29]:


session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName]


# In[57]:


newExperiment.resources


# In[34]:


dir(session.classes)


# #### Following example of "Object creation" in XNATpy tutorial, I'll try to create a new RoiCollectionData:

# In[71]:


ROIcollection = session.classes.RoiCollectionData(parent=newExperiment, label='ROI_tumour')


# #### The above command didn't work previously (status 500 error).  So maybe that is a GET command rather than a PUT command, and worked this time because I've since created an assessor with the label ROI_tumour (see below). 

# In[72]:


ROIcollection


# In[78]:


# ROIcollection.items() # AttributeError: 'RoiCollectionData' object has no attribute 'items'
# ROIcollection.values() # AttributeError: 'RoiCollectionData' object has no attribute 'values'
# ROIcollection.elements() # AttributeError: 'RoiCollectionData' object has no attribute 'elements'


# In[ ]:





# #### The following cells were executed previously:

# Create an assessor following Hakim's guidance here:
# https://groups.google.com/forum/#!searchin/xnat_discussion/xnatpy$20status$20500|sort:date/xnat_discussion/N7i2tS5MleE/7Sxlr7JIBAAJ

# In[59]:


assessorLabel = 'ROI_tumour'
newExperiment.create_assessor(assessorLabel, type_='xnat:qcAssessmentData')  
# Note the xsi_type here, this can be something else if you need it to be


# This has created an Assessment with the Experiment label "Auto QC" in the WI.

# In[62]:


session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName].accessors[assessorLabel]


# In[63]:


session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName].accessors


# #### I don't understand why I got those outputs..

# In[60]:


# Now see if I can add files:

newExperiment.create_assessor(assessorLabel, type_='xnat:qcAssessmentData', files=roiFpathDL)  # Note the xsi_type here, this can be something else if you need it to be


# Clearly not like that.

# In[64]:


assessor = session.projects[newProjectName].subjects[newSubjectName].experiments[newExperimentName]
assessor


# #### Next will try the "Alternative strategy (more generic)" provided by Hakim:

# In[65]:


assessorLabel2 = 'ROI_tumour_2'

assessor = session.classes.QcAssessmentData(parent=newExperiment, label=assessorLabel2)


# This has created another Assessment with the Experiment label "Auto QC" in the WI.

# ## To summarise both of the following commands created an assessment with label "Auto QC":
# 
# ### newExperiment.create_assessor(assessorLabel, type_='xnat:qcAssessmentData')  
# 
# ### session.classes.QcAssessmentData(parent=newExperiment, label=assessorLabel2)

# In[81]:


# Try:
assessorLabel3 = 'ROI_tumour_3' # doesn't work (Received response with status code: 400)
#assessorLabel3 = 'ROI_tumour_2' # works (<RoiCollectionData ROI_tumour_2 (XNAT_E00026)>)


session.classes.RoiCollectionData(parent=newExperiment, label=assessorLabel3)


# In[82]:


# Try:
#xuassessorLabel3 = 'ROI_tumour_3' # doesn't work (Received response with status code: 400)
assessorLabel3 = 'ROI_tumour_2' # works (<RoiCollectionData ROI_tumour_2 (XNAT_E00026)>)


session.classes.RoiCollectionData(parent=newExperiment, label=assessorLabel3)


# #### see comment above (I think this is a GET rather than PUT command)

# In[66]:


# Continuing with the alternative strategy:
resource = session.classes.ResourceCatalog(parent=assessor, label='RESOURCE_LABEL') # was connection.classes... 
# but might have been a typo?


# #### I'm not sure what this has done.  The WI looks the same other than perhaps a "+ Session QC Status" above Scans that may or may not have been present before.

# In[67]:


assessor


# In[68]:


resource


# #### So now need to figure out how to upload the DICOM-RTSTRUCT file to assessor.  Try importing via the prearchive:

# #### Should project be the object (e.g. newProject) or the label (e.g. newProjectName)?

# In[83]:


# Setting project = newProjectName:

prearchive_session = session.services.import_(roiFpathDL,                                               project=newProjectName,                                               destination='/prearchive')


# #### Looking at the log files:
# 
# ##### access.log (tail /data/xnat/home/logs/access.log):
# 
# 2020-02-12 16:13:30,178 - admin 127.0.0.1 POST http://10.1.1.17/data/services/import?dest=%2Fprearchive&project=BrainTumorProg_copy
# 
# ##### dicom.log (tail /data/xnat/home/logs/dicom.log):
# 
# 2020-02-12 16:13:31,811 [pool-2-thread-1] ERROR org.nrg.dcm.xnat.DICOMSessionBuilder - Session builder not implemented for SOP class [1.2.840.10008.5.1.4.1.1.481.3] or modality [RTSTRUCT]
# 
# 

# In[84]:


# Setting project = newProject:

prearchive_session = session.services.import_(roiFpathDL,                                               project=newProject,                                               destination='/prearchive')


# #### Looking at the log files:
# 
# ##### access.log (tail /data/xnat/home/logs/access.log):
# 
# 2020-02-12 16:14:28,221 - admin 0:0:0:0:0:0:0:1 POST http://10.1.1.17/data/services/import?dest=%2Fprearchive&project=%3CProjectData+BrainTumorProg_copy+%28BrainTumorProg_copy%29%3E
# 
# This one looks more suspicious than the previous one (with project=newProjectName) because of "0:0:0:0:0:0:0:1" (as opposed to "127.0.0.1")
# 
# 
# ##### dicom.log (tail /data/xnat/home/logs/dicom.log):
# 
# 2020-02-12 16:14:28,565 [pool-2-thread-2] ERROR org.nrg.dcm.xnat.DICOMSessionBuilder - Session builder not implemented for SOP class [1.2.840.10008.5.1.4.1.1.481.3] or modality [RTSTRUCT]
# 
# 

# In[85]:


# This time trying without the pre-archive destination input:

prearchive_session = session.services.import_(roiFpathDL, project=newProjectName)


# #### Looking at the log files:
# 
# ##### access.log (tail /data/xnat/home/logs/access.log):
# 
# 2020-02-12 16:23:36,863 - admin 127.0.0.1 POST http://10.1.1.17/data/services/import?project=BrainTumorProg_copy
# 
# 
# ##### dicom.log (tail /data/xnat/home/logs/dicom.log):
# 
# 2020-02-12 16:23:37,206 [pool-2-thread-3] ERROR org.nrg.dcm.xnat.DICOMSessionBuilder - Session builder not implemented for SOP class [1.2.840.10008.5.1.4.1.1.481.3] or modality [RTSTRUCT]
# 
# 

# #### Trying to follow the advice by Hakim here:
# 
# ##### https://groups.google.com/forum/#!topic/xnat_discussion/dmM9_cwbn38

# In[44]:


assessorLabel = 'ROI_tumour'

assessor = session.classes.ResourceCatalog(parent=newExperiment, label=assessorLabel)


# In[42]:


accessorLabel = 'DICOM'

assessor = session.classes.ResourceCatalog(parent=newExperiment, label=assessorLabel)


# In[46]:


# This seems to have worked before:
#newExperiment.create_assessor(assessorLabel, type_='xnat:qcAssessmentData')  

# Tried others:
# newExperiment.create_assessor(assessorLabel) # <-- TypeError: create_assessor() missing 1 required positional argument: 'type_'
# newExperiment.create_assessor(assessorLabel, type_='xnat:RoiCollectionData')  # <-- XNATResponseError


# In[49]:


newExperiment.create_assessor(assessorLabel, type_='xnat:qcAssessmentData')


# In[51]:


resourceLabel = 'ROI_tumour'

#resource = session.classes.ResourceCatalog(parent=assessor, label='RESOURCE_LABEL') # was connection.classes... 
# but might have been a typo?
resource = session.classes.ResourceCatalog(parent=assessor, label=resourceLabel) 


# In[ ]:





# In[92]:


# Instead of newProject try newExperiment (probably won't work...): 

prearchive_session = session.services.import_(roiFpathDL,                                               project=newExperiment,                                               destination='/prearchive')


# #### Looking at the log files:
# 
# ##### access.log (tail /data/xnat/home/logs/access.log):
# 
# 2020-02-12 17:01:23,921 - admin 127.0.0.1 POST http://10.1.1.17/data/services/import?dest=%2Fprearchive&project=%3CMrSessionData+PGM-002_copy_MR_2+%28XNAT_E00015%29%3E
# 
# 
# 
# ##### dicom.log (tail /data/xnat/home/logs/dicom.log):
# 
# 2020-02-12 17:01:24,206 [pool-2-thread-4] ERROR org.nrg.dcm.xnat.DICOMSessionBuilder - Session builder not implemented for SOP class [1.2.840.10008.5.1.4.1.1.481.3] or modality [RTSTRUCT]
# 
# 

# In[93]:


# Instead of newProject try newExperiment (probably won't work...): 

prearchive_session = session.services.import_(roiFpathDL,                                               project=newExperimentName,                                               destination='/prearchive')


# #### Looking at the log files:
# 
# ##### access.log (tail /data/xnat/home/logs/access.log):
# 
# 2020-02-12 17:01:43,974 - admin 0:0:0:0:0:0:0:1 POST http://10.1.1.17/data/services/import?dest=%2Fprearchive&project=PGM-002_copy_MR_2
# 
# 
# 
# ##### dicom.log (tail /data/xnat/home/logs/dicom.log):
# 
# 2020-02-12 17:01:44,245 [pool-2-thread-1] ERROR org.nrg.dcm.xnat.DICOMSessionBuilder - Session builder not implemented for SOP class [1.2.840.10008.5.1.4.1.1.481.3] or modality [RTSTRUCT]
# 

# # 13/02/20:
# 
# #### While trying to upload a DICOM scan to the prearchive I've noticed that there are two entries in the prearchive showing errors on Feb 6 and 7.  
# 
# ##### Details:
# 
# ##### XNAT has encountered an error with your request:
# ##### Exception: org.xml.sax.SAXParseException; lineNumber: 1; columnNumber: 1; Premature end of file.
# ##### If this error continues to occur, please contact your system administrator with information about how to recreate the problem.

# # Giving up on trying to upload ROIs for now.  I'll try to upload DICOMs instead..

# In[86]:


dicomsPath


# In[88]:


os.listdir(dicomsPath)


# In[90]:


dicomsFpaths = dicom_functions.dc_filepaths(dicomsPath)
dicomsFpaths


# In[95]:


dicomsFpaths[0]


# In[91]:


newExperimentName


# In[96]:


# Try uploading the first DICOM file in dicomsFpaths to the prearchive:

prearchive_session = session.services.import_(dicomsFpaths[0],                                               project=newProjectName,                                               destination='/prearchive')


# #### There's now a "MR Session" under Experiment for subject PGM-002_copy.  I went to the prearchive, then clicked on "Review and Archieve" and noticed that the subject was not identified.  I was able to select PGM-002_copy from the pull-down list.  I wonder if I can add a subject name to the list of inputs to .services.import_()...

# #### Actually I can see that the archiving process including the assignment of a subject can be done as shown here:
# 
# ##### https://xnat.readthedocs.io/en/latest/static/tutorial.html#importing-data-into-xnat

# In[97]:


# Try uploading the second DICOM scan to the prearchive:

prearchive_session = session.services.import_(dicomsFpaths[1],                                               project=newProjectName,                                               destination='/prearchive')


# In[98]:


print(prearchive_session)


# In[99]:


# Get a list of sessions waiting for archiving:

session.prearchive.sessions()


# In[103]:


# The first two items are the errors mentioned above.  So archive the third entry only:

prearchive_session = session.prearchive.sessions()[2]

experiment = prearchive_session.archive(subject=newSubjectName, experiment=newExperimentName)


# #### The prearchive using the web interface is showing that there's a conflict.
# 
# ##### Details:
# 
# Current Warnings
# FAIL-4: Invalid modification of session subject via archive process.: PGM-002_MR_1 Already Exists for another Subject
# CONFLICT-14: Session already exists with matching files.
# CONFLICT-16: Session already contains a scan (11) with the same UID and number.
# 
# History
# 
# 02/12/2020 20:48:18
# Session already exists, retry with overwrite enabled

# In[104]:


session.prearchive.sessions()


# #### Seems that the scan I'm trying to archive is now listed at the top (the first entry) whereas previously it was the third item.  
# 
# #### This inconsistency is confusing..

# In[105]:


session.prearchive.sessions()[0]


# In[105]:


session.prearchive.sessions()[0]


# In[106]:


# That's the one I want to archive. So archive the first entry only:

prearchive_session = session.prearchive.sessions()[0]

experiment = prearchive_session.archive(subject=newSubjectName, experiment=newExperimentName)


# #### I got the same conflict mesesage

# In[107]:


prearchive_session


# In[108]:


# Try again but this time using the overwrite input

experiment = prearchive_session.archive(subject=newSubjectName, experiment=newExperimentName, overwrite='append')


# #### So that seems to have avoided the conflict error but I was expecting the second scan to be included in the same collection as the first one that I uploaded.
# 
# #### Instead it's been added to a new experiment with a different label (PGM-002_copy_MR_2).  
# 
# ##### Experiments 

# Date          Experiment        Project                Label 
# 
# 1996-08-13    MR Session        BrainTumorProg_copy    PGM-002_copy_MR_2
#               Auto QC           BrainTumorProg_copy    ROI_tumour
#               Auto QC           BrainTumorProg_copy    ROI_tumour_2
# 1996-08-13    MR Session        BrainTumorProg_copy    PGM-002_MR_1

# In[111]:


newProjectName


# In[109]:


newSubjectName


# In[110]:


newExperimentName


# #### So it seems that the scan that I moved to the archive using the web interface was assigned to a different experiment (PGM-002_MR_1) instead of the one defined by newExperimentName (PGM-002_copy_MR_2).

# In[112]:


# Try again with uploading the second DICOM scan to the prearchive:

prearchive_session = session.services.import_(dicomsFpaths[1],                                               project=newProjectName,                                               destination='/prearchive')


# In[113]:


# Try again but this time using the overwrite input but this time adding it to the experiment PGM-002_MR_1:

experiment = prearchive_session.archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='append')


# #### I received another conflict error:
# 
# FAIL-4: Invalid modification of session subject via archive process.: PGM-002_MR_1 Already Exists for another Subject
# CONFLICT-14: Session already exists with matching files.
# CONFLICT-16: Session already contains a scan (11) with the same UID and number.
# 
# History
# 
# 02/12/2020 21:35:16
# Session already exists with matching files.
# 
# #### Maybe because I forgot to select which file in the prearchive to archive!

# In[115]:


# Get a list of sessions waiting for archiving:

session.prearchive.sessions()


# In[116]:


session.prearchive.sessions()[2]


# In[117]:


# Try again this time selecting the correct file to archive:

session.prearchive.sessions()[2].archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='append')


# #### I got the same conflict error!

# In[118]:


# Get a list of sessions waiting for archiving:

session.prearchive.sessions()


# #### This time its the first entry!

# In[119]:


# Try again this time selecting the correct file to archive:

session.prearchive.sessions()[0].archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='append')


# In[120]:


# Get a list of sessions waiting for archiving:

session.prearchive.sessions()


# #### And now it's the 3rd entry again..

# In[121]:


# This time select the 3rd file and use 'delete' for overwrite:

session.prearchive.sessions()[2].archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='delete')


# In[122]:


# Try again this time using a noticeably different scan:

prearchive_session = session.services.import_(dicomsFpaths[20],                                               project=newProjectName,                                               destination='/prearchive')

prearchive_session


# In[123]:


session.prearchive.sessions()


# In[124]:


session.prearchive.sessions()[0].archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='append')


# In[125]:


session.prearchive.sessions()


# In[126]:


session.prearchive.sessions()[2].archive(subject=newSubjectName, experiment='PGM-002_MR_1', overwrite='delete')


# #### It seems it didn't work despite no errors.  The original (first) scan should have been replaced with the 21st scan but I see no evidence of this with the web interface!

# In[ ]:




