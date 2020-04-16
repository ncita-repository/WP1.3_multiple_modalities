#!/usr/bin/env python
# coding: utf-8

# # Overview of exercises:
# 
# ## Work with "Cristiano_test_data_2" provided by Simon from the XNAT_CRUK account.  In this data set the subject position is meant to be the same for all scans.  What varies is:
# 
# ### - Image contrast
# ### - Image FOV (resolution)
# ### - Different number of slices
# 
# #### Note:  All Dixon scans have the same FOV and number of planes.  Only the contrast varies.
# 
# ## Exercise 1 - Copy ROI Collection to scan with different resolution but same number of planes:
# 
# ### Copy ROI Collection from scan 3 ("t2_blade_tra", res = 320 x 320, 128 slices) to scan 4 ("t2_tra", res = 320 x 190, 128 slices). 
# 
# ## Exercise 2 - Copy ROI Collection to scan with different resolution and different number of planes:
# 
# ### Copy ROI Collection from scan 4 ("t2_tra", res = 320 x 190, 128 slices) to scan 7 ("Dixon_T1_opp", res = 384 x 240, 192 slices) 

# ### March 3:
# 
# ### Since I was unable to establish a connection to the XNAT CRUK site via XNATpy, I downloaded the entire patient collection using the web interface instead.
# 
# ### The root folder of the files are at:
# 
# ### C:\Data\Cristiano_test_data_2\20170302_112746_Aera
# 
# ### In my own XNAT instance I used the XNAT Compressed Uploader (Option 2 using the web interface) I selected the zip file that was downloaded in the XNAT_CRUK session.
# 
# ### I loaded the images in the OHIF viewer and selected Dixon_T1_opp (scan session #7), and drew some contours in the 11th frame.
# 
# ### March 4:
# 
# ### I later drew contours on slice 5 of scan 3 ("t2_blade_tra") and named it "ROIs_t2_blade_tra_slice_5". 
# 
# ### March 5:
# 
# ### I started wrapping up the individual sections of code into individual helper functions, correcting bugs along the way.  I succeeded in completing Exercise 1. 
# 
# ### March 6:
# 
# ### I continued wrapping up code into functions including one main function to call up the others.  I succeeded in completing Exercise 2. 

# In[ ]:


# Import packages and functions:

import DicomHelperFuncs
importlib.reload(DicomHelperFuncs)
from DicomHelperFuncs import CopyRoi


# ### Exercise 1:

# In[ ]:


# Run CopyROI() to do everything:
roi2, xnatDict = CopyRoi(fromRoiProjLab='TestData',                          fromRoiSubjLab='Volunteer_1',                          fromRoiExpLab='20170302_112746_Aera',                          fromRoiLabDateTime='20200304_110009',                          fromDicomScanLab='3',                          toDicomProjLab='TestData',                          toDicomSubjLab='Volunteer_1',                          toDicomExpLab='20170302_112746_Aera',                          toDicomScanLab='4',                          searchBy='z',                          del_dicoms=True,                          debug=False)


# ### Exercise 2:

# In[ ]:


# Run CopyROI() to do everything:
roi2, xnatDict = CopyRoi(fromRoiProjLab='TestData',                          fromRoiSubjLab='Volunteer_1',                          fromRoiExpLab='20170302_112746_Aera',                          fromRoiLabDateTime='20200306_114104',                          fromDicomScanLab='4',                          toDicomProjLab='TestData',                          toDicomSubjLab='Volunteer_1',                          toDicomExpLab='20170302_112746_Aera',                          toDicomScanLab='7',                          searchBy='z',                          del_dicoms=True,                          debug=False)

