#!/usr/bin/env python
# coding: utf-8

# ## Use case scenarios:
# 
# 1) Same DICOM Series/XNAT Scan, FoR, IPP, IOP, ST, PS (= Direct copy)
# 
# 2) Different DICOM Series/XNAT Scans; Same DICOM Study/XNAT Session, FoR, IPP, IOP, ST and PS
# 
#     a) Direct copy
#     
#     b) Relationship-preserving copy
# 
# 3) Different Series, PS and/or ST (if different ST they will likely have different IPP); Same Study, FoR, IOP 
# 
#     a) Direct copy
#     
#     b) Relationship-preserving copy
# 
# 4) Different Series, IOP, IPP, possibly different ST and PS; Same Study (= Relationship-preserving copy)
# 
# 5) Different or Same* Study, Different FoR (hence different IPP and IOP), and possibly different ST and/or PS (= Relationship-preserving copy)
# 

# # Possible classes within SimpleITK for image resampling:
# 
#     ResampleImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1ResampleImageFilter.html#details)
# 
#     GridImageSource (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1GridImageSource.html)
#     
#     
# # Other classes that may be useful for other needs:
# 
#     SimpleContourExtractor (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1SimpleContourExtractorImageFilter.html)
#     
#     BinaryProjectionImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1BinaryProjectionImageFilter.html)
#     
#     ContourDistanceImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1IsoContourDistanceImageFilter.html)
#     
#     LabelMapMaskImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1LabelMapMaskImageFilter.html)
#     
#     LabelMapToBinaryImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1LabelMapToBinaryImageFilter.html)
#     
#     PhysicalPointImageSource (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1PhysicalPointImageSource.html)
#     
#     RegionOfInterestImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1RegionOfInterestImageFilter.html#details)
#     
#     ExtractImageFilter (https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1ExtractImageFilter.html)
#     
#     TransformPhysicalPointToIndex() https://simpleitk.org/doxygen/v1_0/html/classitk_1_1simple_1_1Image.html#acfbba23767a42faa45d90d7548d2f53b
#     
#     TransformIndexToPhysicalPoint() https://simpleitk.org/doxygen/v1_0/html/classitk_1_1simple_1_1Image.html#aa48e3e8c7d03fdebf7a80e046617309a
#     
#     TransformContinuousIndexToPhysicalPoint() https://simpleitk.org/doxygen/v1_0/html/classitk_1_1simple_1_1Image.html#a6c89c8942bdf2c0e30ffc8a1e566fc0f
# 

# In[1]:


import os
import time
import copy
#import scipy
import numpy as np
import importlib
import pydicom
#import DicomTools
#importlib.reload(DicomTools)
#from DicomTools import GetDicoms
#from DicomTools import Get3DImage
#from DicomTools import GetImageAttributes
import SimpleITK as sitk
from RunSimpleElastixReg import RunSimpleElastixReg
import PlotDicomsAndRois_v2
importlib.reload(PlotDicomsAndRois_v2)
from PlotDicomsAndRois_v2 import PlotDicomsAndRois_v2 as PlotStacks
#%matplotlib inline
#% autoreload
#%DicomTools autoreload


# In[2]:


DicomRootDir = r"C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011"

MR4Dir = r"04-10-1960-MRI Brain wwo Contrast-69626"
MR12Dir = r"06-11-1961-MRI Brain wwo Contrast-79433"

MR4_S3_DicomDir = os.path.join(DicomRootDir, MR4Dir,
                               r"3-T2 TSE AXIAL-07507") # 

MR4_S4_DicomDir = os.path.join(DicomRootDir, MR4Dir,
                              r"4-T2 AXIAL FLAIR DARK FL-48830") # 416x512x21 .449x.449x5 mm

MR4_S5_DicomDir = os.path.join(DicomRootDir, MR4Dir,
                               r"5-T1 SE AXIAL 3MM-81246")        # 208x256x50 .898x.898x3 mm

MR4_S6_DicomDir = os.path.join(DicomRootDir, MR4Dir, 
                               r"6-ep2ddiff3scantracep2-55446") # 192x192x63 1.302x1.302x5 mm

MR4_S7_DicomDir = os.path.join(DicomRootDir, MR4Dir, 
                               r"7-ep2ddiff3scantracep2ADC-22197") # 

MR4_S8_DicomDir = os.path.join(DicomRootDir, MR4Dir, 
                               r"8-T1 SE AXIAL POST FS FC-59362") # 212x256x30 .898x.898x5 mm

MR4_S9_DicomDir = os.path.join(DicomRootDir, MR4Dir,
                               r"9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268") # 208x256x50 .898x.898x3 mm

MR12_S8_DicomDir = os.path.join(DicomRootDir, MR12Dir,
                                r"8-T1 SE AXIAL POST FS FC-81428") # 256x192x35 .938x.938x5 mm

RoiRootDir = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011"
RtsRootDir = os.path.join(RoiRootDir, "Fixed_RTS_Collections")
SegRootDir = os.path.join(RoiRootDir, "Fixed_SEG_Collections")

MR4_S8_RtsFpath = os.path.join(RtsRootDir, "MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") # tumour
MR4_S8_SegFpath = os.path.join(SegRootDir, "Seg_from_MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") # tumour

MR4_S9_RtsFpath = os.path.join(RtsRootDir, "MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") # tumour drawn on S9 stack
MR4_S9_SegFpath = os.path.join(SegRootDir, "Seg_from_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") # tumour drawn on S9 stack
MR4_S9_SegFpath = os.path.join(SegRootDir, "MR4_S9_tumour_SEG_20201105_115408.dcm")
MR4_S9_SegFpath = os.path.join(SegRootDir, "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")

MR4_S5_RtsFpath = os.path.join(RtsRootDir, "MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") # tumour drawn on S5 stack
MR4_S5_SegFpath = os.path.join(SegRootDir, "Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") # tumour drawn on S5 stack

MR4_S3_RtsFpath = os.path.join(RtsRootDir, "") 
MR4_S3_SegFpath = os.path.join(SegRootDir, "MR4_S3_brain_SEG_20201030_133110.dcm")
MR4_S3_SegFpath = os.path.join(SegRootDir, "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm")

RoiExportDir = os.path.join(RoiRootDir, "New_ROI_Collections")

PlotExportDir = r'C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011\Results_of_Copy_ROI_use_cases'


# ### Work on Case #1 - Use Series 9 (=Src) and 9 (=Trg):
# 
# #### 1) Create contours for Series 9 in OHIF-Viewer, download RTS from XNAT
# #### 2) Convert RTS to SEG
# #### 3) Copy RTS and SEG and modify tags to match the DICOMs
# #### 4) Export RTS and SEG, upload to XNAT using RoiUploadAssistant, import in OHIF-Viewer, verify that ROIs are rendered as expected

# SrcDicomDir = S011_MR4_S9_DicomDir
# SrcRtsFpath = S011_MR4_S9_RtsFpath
# SrcSegFpath = S011_MR4_S9_SegFpath
# 
# TrgDicomDir = S011_MR4_S5_DicomDir

# # Try this:
# # https://theaisummer.com/medical-image-coordinates/#the-coordinate-systems-in-medical-imaging
# # Stuff missing in code:
# # https://nipy.org/nibabel/coordinate_systems.html
# 
# aff_t1 = FixImNpa.affine
# aff_t2 = MovImNpa.affine
# inv_af_2 = np.linalg.inv(aff_t2)
# out_shape = FixImNpa.get_fdata().shape
# 
# # desired transformation
# T = inv_af_2.dot(aff_t1)
# 
# # apply transformation
# transformed_img = scipy.ndimage.affine_transform(MovImNpa.get_fdata(), T, output_shape=out_shape)

# ### Read in the DICOMs:

# In[15]:


S011_MR4_S9_Dicoms = GetDicoms(DicomDir=S011_MR4_S9_DicomDir, SortMethod='slices', Debug=False)
S011_MR4_S9_SitkIm, S011_MR4_S9_Npa = Get3DImage(DicomDir=S011_MR4_S9_DicomDir)


# ### Download the new contour ROI collection for Series 9:

# In[33]:


from DownloadRoiFromXnat import DownloadRoiFromXnat

ProjectLabel = 'ACRIN'

SubjectLabel = 'ACRIN-FMISO-Brain-011'

ExperimentLabel = 'ACRIN-FMISO-Brain-011_MR_4'

DateTimeInLabel = '20201014_163721'
DateTimeInLabel = '20201014'
DateTimeInLabel = '20200910'

FixRts = DownloadRoiFromXnat(ProjectLabel=ProjectLabel,
                             SubjectLabel=SubjectLabel,
                             ExperimentLabel=ExperimentLabel,
                             DateTimeInLabel=DateTimeInLabel,
                             #DownloadDir=RtsRootDir)
                             DownloadDir=r'C:\Temp')


# ### Don't know why I got the above error.  I even tried re-running the script 20200910_Download_ROI_collection_from_XNAT, which previously worked but for some reason didn't work today.
# 
# ### So I instead downloaded the RTS file from the XNAT UI.

# In[3]:


# Read in the RTS file:
S011_MR4_S9_Rts = pydicom.read_file(S011_MR4_S9_RtsFpath)

S011_MR4_S9_Rts


# ### Try running James' RTS-to-SEG command line conversion here:

# In[15]:


os.system('cmd /k "C:\Users\ctorti\etherj-cli-tools-1.0.0\etherj-cli-tools\bin\rtstoseg --help"')


# In[14]:


os.system('cmd /k r"C:\Users\ctorti\etherj-cli-tools-1.0.0\etherj-cli-tools\bin\rtstoseg --help"')


# In[16]:


os.system('cmd /k "C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin/rtstoseg --help"')


# In[17]:


print(os.system('cmd /k "C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin/rtstoseg --help"'))


# In[18]:


get_ipython().system('C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin/rtstoseg --help')


# ### --> Using ! to execute commands works! 

# In[27]:


BinDir = r'C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin'
#BinDir = r"C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin"

get_ipython().system('{BinDir}/rtstoseg --help')


# ### --> Can use variables using {}.

# In[28]:


get_ipython().system('{BinDir}/rtstoseg -s {S011_MR4_S9_DicomDir} {S011_MR4_S9_RoiFpath}')


# ### --> Although using !{BinDir}/rtstoseg worked when running the help command, but for some reason doesn't when trying to run a conversion...

# In[29]:


get_ipython().system('C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin/rtstoseg -s {S011_MR4_S9_DicomDir} {S011_MR4_S9_RoiFpath}')


# ### --> Seems that when defining paths with apostrophes, the variable paths are recognised using curly brackets but the presence of white spaces is a problem.  Defining the paths with quotations should get around the white space issue, but it seems that the use of curly brackets no longer works (seems it's not recognised as a variable any more).

# In[8]:


S011_MR4_S9_DicomDir


# In[10]:


S011_MR4_S9_RtsFpath


# ### Avoid use of variables:

# In[32]:


get_ipython().system('C:/Users/ctorti/etherj-cli-tools-1.0.0/etherj-cli-tools/bin/rtstoseg -s "C:\\\\Data\\\\Cancer imaging archive\\\\ACRIN-FMISO-Brain\\\\ACRIN-FMISO-Brain-011\\\\04-10-1960-MRI Brain wwo Contrast-69626\\\\9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268" "C:\\\\Data\\\\Data copied to Google Drive\\\\ACRIN-FMISO-Brain-011\\\\Fixed_RTS_Collections\\\\AIM_20201014_163721.dcm"')


# ### --> It looked like it should have worked since I get the same output when running it in a command shell, and despite the errors, a SEG file is generated.   But running the above cell doesn't seem to have generated any file (not in the bin folder, or the source DICOMs or the source RTS file).  
# 
# ### --> What works so far: Copying the RTS and DICOMs into the bin folder and running the script in a command shell!

# ### Read in the converted SEG file:

# In[4]:


# Read in the SEG file:
S011_MR4_S9_Seg = pydicom.read_file(S011_MR4_S9_SegFpath)


# In[5]:


S011_MR4_S9_Seg


# In[6]:


S011_MR4_S9_Seg.pixel_array


# In[7]:


print(S011_MR4_S9_Seg.Rows)
print(S011_MR4_S9_Seg.Columns)
print(np.shape(S011_MR4_S9_Seg.pixel_array))


# In[32]:


print(S011_MR4_S9_Rts.StructureSetLabel)
print(S011_MR4_S9_Rts.StructureSetROISequence[0].ROIName)
print(S011_MR4_S9_Rts.StructureSetDate)
print(S011_MR4_S9_Rts.StructureSetTime)


# In[20]:


S011_MR4_S9_Rts


# # USE CASE 1a:  Copy contour within same series (Series 9)

# In[158]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopyContourWithinSeries

NewRtsRoi = CopyContourWithinSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                    SrcRtsFpath=S011_MR4_S9_RtsFpath, 
                                    FromSliceNum=23, 
                                    ToSliceNum=22, 
                                    LogToConsole=False)


# In[153]:


NewRtsRoi


# In[159]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportRtsRoi

ExportRtsRoi(TrgRtsRoi=NewRtsRoi, 
             SrcRtsFpath=S011_MR4_S9_RtsFpath, 
             NamePrefix='s23_to_s22',
             ExportDir=RoiExportDir)


# In[136]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetDicomSOPuids
from RoiCopyTools import GetCIStoDcmInds

GetCIStoDcmInds(RtsRoi=S011_MR4_S9_Rts, SOPuids=GetDicomSOPuids(DicomDir=S011_MR4_S9_DicomDir))


# In[36]:


# The Source RTS ROI:
S011_MR4_S9_Rts


# In[56]:


# The new RTS ROI:
NewRtsRoi


# # Exported RTSTRUCT, uploaded to XNAT, imported in OHIF-Viewer and checked results.  All good.

# In[459]:


S011_MR4_S9_Dicoms[0]


# ### Plot the DICOMs with the original RTS contours:

# In[464]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndContours

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True


PlotDicomsAndContours(DicomDir=S011_MR4_S9_DicomDir, 
                      RtsFpath=S011_MR4_S9_RtsFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# ### Plot the DICOMs with the new RTS contours:

# In[466]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndContours

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True

# Filepath to the new RTS ROI:
NewRtsFpath = os.path.join(RoiExportDir, r"20201019_130006_Case_1a_s23_to_s22_from_AIM_20201014_163721.dcm")

PlotDicomsAndContours(DicomDir=S011_MR4_S9_DicomDir, 
                      RtsFpath=NewRtsFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# # USE CASE 1b:  Copy segmentation within same series (Series 9)

# In[44]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetDicomSOPuids
from RoiCopyTools import ContourImageSequencesToDicomSlice

# The SOP UIDs for S9:
S011_MR4_S9_SOPuids = GetDicomSOPuids(S011_MR4_S9_DicomDir)

#S011_MR4_S9_SOPuids


# In[43]:


# The ContourImageSequence-to-DICOMSliceInds:
S011_MR4_S9_CisToDcmInds = ContourImageSequencesToDicomSlice(S011_MR4_S9_Rts, S011_MR4_S9_SOPuids)

S011_MR4_S9_CisToDcmInds


# In[48]:


# The SOP UIDs of the contours in S9's RTS:
[S011_MR4_S9_SOPuids[i] for i in S011_MR4_S9_CisToDcmInds]


# In[49]:


S011_MR4_S9_Rts


# In[50]:


S011_MR4_S9_Seg


# In[65]:


S011_MR4_S9_Origin, S011_MR4_S9_Dirs,S011_MR4_S9_Spacings, S011_MR4_S9_Dims = GetImageAttributes(S011_MR4_S9_DicomDir, 'pydicom')

S011_MR4_S9_Dirs[:6]


# In[73]:


S011_MR4_S9_Seg.PixelData


# In[72]:


S011_MR4_S9_Seg.pixel_array


# In[74]:


np.shape(S011_MR4_S9_Seg.pixel_array)


# In[75]:


S011_MR4_S9_Seg.pixel_array[3]


# In[80]:


a = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 1]])

a


# In[81]:


np.insert(a, 1, [1, 1, 0], axis=0)


# In[84]:


a


# In[85]:


b = np.insert(a, 1, [1, 1, 0], axis=0)

b


# In[83]:


b = copy.deepcopy(a)

b


# orig = copy.deepcopy(S011_MR4_S9_Seg.pixel_array)
# 
# new = orig.append
# print(np.shape(orig))
# 
# 

# In[86]:


S011_MR4_S9_Dicoms = GetDicoms(S011_MR4_S9_DicomDir, 'slices', False)


# In[88]:


S011_MR4_S9_Dicoms[0]


# In[91]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopySegmentWithinScan

NewSegRoi = CopySegmentWithinScan(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                  SrcSegFpath=S011_MR4_S9_SegFpath, 
                                  FromSliceNum=23, 
                                  ToSliceNum=22, 
                                  LogToConsole=False)


# In[93]:


a = S011_MR4_S9_Seg.pixel_array

np.shape(a)


# In[97]:


b = np.append(a, a[-1], axis=0)


# In[98]:


for i in range(np.shape(a)[0]):
    print(np.shape(a[i]))


# In[99]:


np.shape(a[-1])


# In[102]:


np.dstack((a, a[-1]))


# In[103]:


b = a[-1]

print(a.shape)
print(b.shape)


# In[104]:


np.dstack((a, b))


# In[107]:


np.concatenate((a, b), axis=0)


# In[109]:


np.vstack((a,b))


# In[110]:


np.append(a, np.atleast_3d(b), axis=0)


# In[111]:


np.concatenate((a, b))


# In[113]:


c = np.reshape(a, (256, 208, 7))


# In[114]:


c.shape


# In[116]:


d = np.dstack((c, b))

d.shape


# In[117]:


e = np.reshape(d, (8, 256, 208))

e.shape


# In[119]:


c[:,:,0].shape


# ### Since the SEG file was converted from an RTS change the Series Description and Content Label to reflect this:

# In[169]:


S011_MR4_S9_Seg.SeriesDescription = 'SEG_from_' + S011_MR4_S9_Seg.SeriesDescription

S011_MR4_S9_Seg.SegmentSequence[0].SegmentLabel = 'SEG_from_' + S011_MR4_S9_Seg.SegmentSequence[0].SegmentLabel


# ### Export the original SEG conversion with modified Segment Label:

# In[170]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

ExportSegRoi(TrgSegRoi=S011_MR4_S9_Seg, 
             SrcSegFpath=S011_MR4_S9_SegFpath, 
             NamePrefix='Mod',
             ExportDir=SegRootDir)


# ### The RoiUploadAssistant tool uploaded the converted RTS-to-SEG file without issues and I was able to import the SEG in the OHIF-Viewer.

# ### Copy slice 23 to slice 22:

# In[171]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopySegmentWithinSeries

NewSegRoi = CopySegmentWithinSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                    SrcSegFpath=S011_MR4_S9_SegFpath, 
                                    FromSliceNum=23, 
                                    ToSliceNum=22, 
                                    LogToConsole=False)


# In[122]:


S011_MR4_S9_Seg


# In[172]:


NewSegRoi


# ### It seems that PixelData is wrong since the original and new arrays have the same number of elements (46592)...

# In[173]:


NewSegRoi


# In[174]:


# Since the SEG file was converted from an RTS change the Series Description and Content Label to reflect this:
NewSegRoi.SeriesDescription = 'SEG_from_' + NewSegRoi.SeriesDescription

NewSegRoi.SegmentSequence[0].SegmentLabel = 'SEG_from_' + NewSegRoi.SegmentSequence[0].SegmentLabel


# In[175]:


NewSegRoi


# In[177]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

ExportSegRoi(TrgSegRoi=NewSegRoi, 
             SrcSegFpath=S011_MR4_S9_SegFpath, 
             NamePrefix='s23_to_s22',
             ExportDir=RoiExportDir)


# # The RoiUploadAssistant tool uploaded the converted RTS-to-SEG file without issues and but I was unable to import the SEG in the OHIF-Viewer - get a revolving circle (as I did when I tried to export a SEG generated within the OHIF-Viewer).

# ### I tried to convert the modified SEG to RTS using James' tool but got error:
# 
# Exception in thread "main" java.lang.IllegalArgumentException: Requested frame outside input array: 8
#         at icr.etherj.dicom.SegUtils.unpackFrame(SegUtils.java:213)
#         at icr.etherj.dicom.SegUtils.unpackFrame(SegUtils.java:183)
#         at icr.etherj.dicom.impl.SegmentationToRtStructConverter.getFrameMask(SegmentationToRtStructConverter.java:671)
#         at icr.etherj.dicom.impl.SegmentationToRtStructConverter.processFrame(SegmentationToRtStructConverter.java:889)
#         at icr.etherj.dicom.impl.SegmentationToRtStructConverter.processFrames(SegmentationToRtStructConverter.java:937)
#         at icr.etherj.dicom.impl.SegmentationToRtStructConverter.convertImpl(SegmentationToRtStructConverter.java:457)
#         at icr.etherj.dicom.impl.SegmentationToRtStructConverter.convert(SegmentationToRtStructConverter.java:159)
#         at icr.etherj.dicom.impl.DefaultRoiConverter.toRtStruct(DefaultRoiConverter.java:84)
#         at icr.etherj.clients.SegToRts.run(SegToRts.java:215)
#         at icr.etherj.clients.SegToRts.main(SegToRts.java:85)
# 

# In[178]:


S011_MR4_S9_Seg


# In[179]:


NewSegRoi


# In[198]:


a = [1, 1, 1]

print(type(a))
print(sys.getsizeof(a))
print(sys.getsizeof(a[0]))
print('')

a = np.array(a, dtype=np.int32)

print(type(a))
print(a.dtype)
print(sys.getsizeof(a))
print(sys.getsizeof(a[0]))
print('')

a = list(a)

print(type(a))
#print(a[0].bit_length())
print(sys.getsizeof(a))
print(sys.getsizeof(a[0]))


# In[203]:


print(sys.getsizeof(1))
print(sys.getsizeof(10000000000))


# In[215]:


sys.getsizeof(NewSegRoi.PerFrameFunctionalGroupsSequence[0]              .FrameContentSequence[0]              .DimensionIndexValues[0])


# ### So it appears that tag values can misleadingly appear as strings (see https://github.com/pydicom/pydicom/issues/892).
# 
# ### Next I need to work out how to force the bit depth to be 32 rather than 28 bits.
# 
# ### The other issue is the number of elements in the Pixel Data array - they are the same in both original and modified ROIs...
# 
# ### Try using .tobytes() as described here:
# 
# ### https://pydicom.github.io/pydicom/dev/old/working_with_pixel_data.html

# In[258]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopySegmentWithinSeries

NewSegRoi = CopySegmentWithinSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                    SrcSegFpath=S011_MR4_S9_SegFpath, 
                                    FromSliceNum=23, 
                                    ToSliceNum=22, 
                                    LogToConsole=False)

NewSegRoi


# In[259]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

ExportSegRoi(TrgSegRoi=NewSegRoi, 
             SrcSegFpath=S011_MR4_S9_SegFpath, 
             NamePrefix='s23_to_s22',
             ExportDir=RoiExportDir)


# In[253]:


print(S011_MR4_S9_Seg.pixel_array.shape)
print(NewSegRoi.pixel_array.shape)


# In[254]:


print(7*256*208)
print(8*256*208)


# ### using NewSegRoi.PixelArray = New3dMask.tobytes() hasn't worked. <-- the issue was .PixelArray, which should have been .PixelData!!!
# 
# ### Prepare code to post to GitHub Issues:

# In[255]:


import pydicom
import numpy as np
import copy

# Filepath to SEG file:
SegFpath = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011\Fixed_SEG_Collections\Seg_AIM_20201014_163721.dcm"

# Import the Original SEG ROI:
OrigSegRoi = pydicom.dcmread(SegFpath)

# Use OrigSegRoi as a template for NewSegRoi: 
NewSegRoi = copy.deepcopy(OrigSegRoi)

# Next I increase the length of ReferencedInstanceSequence and 
# PerFrameFunctionalGroupsSequence by one.  But skipping that here
# since not required to demonstrate the issue.

# Modify the Number of Frames in NewSegRoi:
#NewSegRoi.NumberOfFrames = f"{int(OrigSegRoi.NumberOfFrames) + 1}"
NewSegRoi.NumberOfFrames = str(int(OrigSegRoi.NumberOfFrames) + 1)
    
# I also assign the appropriate SOP Instance UID for the
# DICOM associated with the new 2D mask by modifying the ReferencedSOPInstanceUID, 
# as well as other tags.  But skipping since not required to demonstrate this issue.

# The original 3D mask is:
Orig3dMask = OrigSegRoi.pixel_array

# The shape of Orig3dMask:
OrigShape = Orig3dMask.shape

# What follows next is a very unelegant way of copying one of the 2D masks and adding it
# to the 3D mask.  I originally thought this would do the trick: 
#
# New3dMask = np.append(Orig3dMask, Orig3dMask[-1], axis=0)
#
# but that didn't work.

# Re-shape New3DMask from (NumOfFrames, Rows, Columns) to 
# (Rows, Columns, NumOfFrames):
New3dMask = np.reshape(Orig3dMask, (OrigShape[1], OrigShape[2], OrigShape[0]))

# Stack the last 2D array to the 3D array:
New3dMask = np.dstack((New3dMask, New3dMask[:,:,-1]))

# The shape of New3dMask:
NewShape = New3dMask.shape

# Reshape New3dMask to the original form:
New3dMask = np.reshape(New3dMask, (NewShape[2], NewShape[0], NewShape[1]))

# The shape of New3dMask:
NewShape = New3dMask.shape

# Convert New3dMask to bytes:
NewSegRoi.PixelData = New3dMask.tobytes()

print('\nThe shape of the original 3D mask             =', OrigShape)
print('\nThe shape of the new 3D mask                  =', NewShape)
print('\nThe shape of the original PixelData tag value =', OrigSegRoi.pixel_array.shape)
print('\nThe shape of the new PixelData tag value      =', NewSegRoi.pixel_array.shape)


# In[256]:


NewSegRoi


# In[257]:


print(S011_MR4_S9_Seg.pixel_array.shape)
print(NewSegRoi.pixel_array.shape)


# In[234]:


print(np.__version__)
print(pydicom.__version__)
#print(copy.__version__)
import platform
print(platform.platform(),
      "\nPython", sys.version,
      "\npydicom", pydicom.__version__)


# ### This may be helpful:
# 
# ### https://github.com/pydicom/pydicom/issues/614

# In[243]:


FrameGenerator = pydicom.encaps.generate_pixel_data_frame(OrigSegRoi.PixelData)

FrameGenerator


# ### Just found this issue:
# 
# ### https://github.com/pydicom/pydicom/issues/1061

# In[240]:


print(OrigSegRoi.file_meta.TransferSyntaxUID)
print(NewSegRoi.file_meta.TransferSyntaxUID)


# In[242]:


print(OrigSegRoi[0x00280000:0x00290000])
print('')
print(NewSegRoi[0x00280000:0x00290000])


# ### Plot the DICOMs and segments of the unmodified SEG file (= conversion of RTS)

# In[441]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndSegments

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True


PlotDicomsAndSegments(DicomDir=S011_MR4_S9_DicomDir, 
                      SegFpath=S011_MR4_S9_SegFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# ### Plot the DICOMs and segments of the modified SEG file

# import RoiCopyTools
# importlib.reload(RoiCopyTools)
# from RoiCopyTools import PlotDicomsAndSegments
# 
# ExportPlot = False
# #ExportPlot = True
# 
# LogToConsole = False
# #LogToConsole = True
# 
# NewSegFpath = os.path.join(RoiExportDir, r'20201021_071652_Case_1b_s23_to_s22_from_Seg_AIM_20201014_163721.dcm')
# 
# PlotDicomsAndSegments(DicomDir=S011_MR4_S9_DicomDir, 
#                       SegFpath=NewSegFpath, 
#                       ExportPlot=ExportPlot, 
#                       ExportDir=PlotExportDir, 
#                       LogToConsole=LogToConsole)

# In[287]:


SegNpa = S011_MR4_S9_Seg.pixel_array

NewSegNpa = NewSegRoi.pixel_array

print(f'\nShape of SegNpa = {SegNpa.shape}')
print(f'\nShape of NewSegNpa = {NewSegNpa.shape}')


# In[296]:


#SegNpa[0]

#np.nonzero(SegNpa)

# The non-zero elements in the first frame of SegNpa:
non0 = np.nonzero(SegNpa[0])

#SegNpa[0, non0[0], non0[1]]

a = non0[0][0]
b = non0[0][-1]

c = non0[1][0]
d = non0[1][-1]

# The non-zero (cropped) array:
arr = SegNpa[0, a:b, c:d]

arr


# In[300]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[301]:


fig, ax = plt.subplots(1, 1)
ax = plt.subplot(1, 1, 1, aspect='equal')
ax.imshow(arr)


# In[303]:


# The non-zero elements in the first frame of NewSegNpa:
newnon0 = np.nonzero(NewSegNpa[0])

newnon0


# In[302]:


a = newnon0[0][0]
b = newnon0[0][-1]

c = newnon0[1][0]
d = newnon0[1][-1]

# The non-zero (cropped) array:
newarr = NewSegNpa[0, a:b, c:d]

newarr


# ### Has something gone wrong with the numpy reshaping stuff?...

# In[321]:


print('Orig3dMask.shape =', Orig3dMask.shape)

# The non-zero elements in the first frame of Orig3dMask:
non0 = np.nonzero(Orig3dMask[0])

print(f'\nlen(non0[0]) x len(non0[1]) = {len(non0[0])} x {len(non0[1])}')

a = non0[0][0]
b = non0[0][-1]

c = non0[1][0]
d = non0[1][-1]

# The non-zero (cropped) array:
Orig3dMaskCropped = Orig3dMask[0, a:b, c:d]

print('\nOrig3dMaskCropped.shape =', Orig3dMaskCropped.shape)

# Re-shape New3DMask from (NumOfFrames, Rows, Columns) to 
# (Rows, Columns, NumOfFrames):
New3dMask = np.reshape(Orig3dMask, (OrigShape[1], OrigShape[2], OrigShape[0]))

print('\n\nNew3dMask.shape =', New3dMask.shape)

# The non-zero elements in the first frame of New3dMask:
non0 = np.nonzero(New3dMask[:,:,0])

print(f'\nlen(non0[0]) x len(non0[1]) = {len(non0[0])} x {len(non0[1])}')

a = non0[0][0]
b = non0[0][-1]

c = non0[1][0]
d = non0[1][-1]

# The non-zero (cropped) array:
New3dMaskCropped = New3dMask[a:b, c:d, 0]

print('\nNew3dMaskCropped.shape =', New3dMaskCropped.shape)

fig, ax = plt.subplots(1, 2)

ax = plt.subplot(1, 2, 1, aspect='equal')
ax.imshow(Orig3dMaskCropped)

ax = plt.subplot(1, 2, 2, aspect='equal')
ax.imshow(New3dMaskCropped)


# In[372]:


print('Orig3dMask.shape =', Orig3dMask.shape, '\n')

#New3dMask = copy.deepcopy(Orig3dMask[0])
New3dMask = np.zeros((1, Orig3dMask.shape[1], Orig3dMask.shape[2]), dtype=int)

print('New3dMask.shape =', New3dMask.shape, '\n')

for i in range(Orig3dMask.shape[0]):
    New3dMask = np.append(New3dMask, Orig3dMask[i], axis=0)
    
    print('New3dMask.shape =', New3dMask.shape)


# In[331]:


print('Orig3dMask.shape =', Orig3dMask.shape, '\n')

#New3dMask = copy.deepcopy(Orig3dMask[0])
New3dMask = np.zeros((1, Orig3dMask.shape[1], Orig3dMask.shape[2]), dtype=int)

print('New3dMask.shape =', New3dMask.shape, '\n')

for i in range(Orig3dMask.shape[0]):
    New3dMask = np.concatenate((New3dMask, Orig3dMask[i]), axis=0)
    
    print('New3dMask.shape =', New3dMask.shape)


# # The problem might have been that the zeros array I initiated was int instead of boolean (binary)!!!

# In[413]:


print('Orig3dMask.shape =', Orig3dMask.shape, '\n')

#New3dMask = copy.deepcopy(Orig3dMask[0])
#New3dMask = np.zeros((Orig3dMask.shape[0] + 1, Orig3dMask.shape[1], Orig3dMask.shape[2]), dtype=int)
New3dMask = np.zeros((Orig3dMask.shape[0] + 1, Orig3dMask.shape[1], Orig3dMask.shape[2]), dtype=bool)

print('New3dMask.shape =', New3dMask.shape, '\n')

# Copy the first frame of Orig3dMask to the first frame in New3dMask:
New3dMask[0] = Orig3dMask[0]

print('New3dMask.shape =', New3dMask.shape, '\n')

# Copy all frames in Orig3dMask to subsequent frames:
for i in range(Orig3dMask.shape[0]):
    New3dMask[1+i] = Orig3dMask[i]
    
    print('New3dMask.shape =', New3dMask.shape)
    

print('\nNew3dMask.dtype =', New3dMask.dtype)


# In[414]:


f = 0
#f = 1

print('Orig3dMask.shape =', Orig3dMask.shape)

# The non-zero elements in frame f of Orig3dMask:
non0 = np.nonzero(Orig3dMask[f])

print(f'\nlen(non0[0]) x len(non0[1]) = {len(non0[0])} x {len(non0[1])}')

a = non0[0][0]
b = non0[0][-1]

c = non0[1][0]
d = non0[1][-1]

# The non-zero (cropped) array:
Orig2dMaskCropped = Orig3dMask[f, a:b, c:d]

print('\nOrig2dMaskCropped.shape =', Orig2dMaskCropped.shape)

# The non-zero elements in frame f of New3dMask:
non0 = np.nonzero(New3dMask[f])

print(f'\nlen(non0[0]) x len(non0[1]) = {len(non0[0])} x {len(non0[1])}')

a = non0[0][0]
b = non0[0][-1]

c = non0[1][0]
d = non0[1][-1]

# The non-zero (cropped) array:
New2dMaskCropped = New3dMask[f, a:b, c:d]

print('\nNew2dMaskCropped.shape =', New2dMaskCropped.shape)

fig, ax = plt.subplots(1, 2)

ax = plt.subplot(1, 2, 1, aspect='equal')
ax.imshow(Orig2dMaskCropped)

ax = plt.subplot(1, 2, 2, aspect='equal')
ax.imshow(New2dMaskCropped)


# In[412]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CropNonZerosIn2dMask

f = 0

print('Orig3dMask.shape =', Orig3dMask.shape)

print('\nCropping Orig2dMask to array of non-zero elements..')
    
Orig2dMaskCropped = CropNonZerosIn2dMask(Orig3dMask[f])

print('\nCropping New2dMask to array of non-zero elements..')
    
New2dMaskCropped = CropNonZerosIn2dMask(New3dMask[f])


fig, ax = plt.subplots(1, 2)

ax = plt.subplot(1, 2, 1, aspect='equal')
ax.imshow(Orig2dMaskCropped)

ax = plt.subplot(1, 2, 2, aspect='equal')
ax.imshow(New2dMaskCropped)


# ### Save the new 3D mask to NewSegRoi

# In[421]:


#import pydicom.pixel_data_handlers.numpy_handler.pack_bits as pack_bits
from pydicom.pixel_data_handlers.numpy_handler import pack_bits


# In[422]:


#NewSegRoi.PixelData = New3dMask.tobytes() # <-- still has problems using dtype = bool and dtype = int

#NewSegRoi.PixelData = pydicom.pixel_data_handlers.numpy_handler.pack_bits(New3dMask) # <-- New3dMask must be int or uint

packed = pack_bits(New3dMask.ravel())

NewSegRoi.PixelData = packed + b'\x00' if len(packed) % 2 else packed


# In[423]:


len(packed)


# In[424]:


NewSegRoi


# ### Repeat cropped plots to see if pixel data was encoded properly:

# In[425]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CropNonZerosIn2dMask

New3dMask = NewSegRoi.pixel_array

f = 0

print('Orig3dMask.shape =', Orig3dMask.shape)

print('\nCropping Orig2dMask to array of non-zero elements..')
    
Orig2dMaskCropped = CropNonZerosIn2dMask(Orig3dMask[f])

print('\nCropping New2dMask to array of non-zero elements..')
    
New2dMaskCropped = CropNonZerosIn2dMask(New3dMask[f])


fig, ax = plt.subplots(1, 2)

ax = plt.subplot(1, 2, 1, aspect='equal')
ax.imshow(Orig2dMaskCropped)

ax = plt.subplot(1, 2, 2, aspect='equal')
ax.imshow(New2dMaskCropped)


# In[426]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import Compare2dMasksFrom3dMasks

Compare2dMasksFrom3dMasks(OrigSegRoi, NewSegRoi, 0)


# ### Corrections made to ModifySegTagVals:

# In[427]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopySegmentWithinSeries

NewSegRoi = CopySegmentWithinSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                    SrcSegFpath=S011_MR4_S9_SegFpath, 
                                    FromSliceNum=23, 
                                    ToSliceNum=22, 
                                    LogToConsole=True)

NewSegRoi


# In[428]:


print('\nThe shape of the original PixelData tag value =', OrigSegRoi.pixel_array.shape)
print('\nThe shape of the new PixelData tag value      =', NewSegRoi.pixel_array.shape)


# In[429]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import Compare2dMasksFrom3dMasks

Compare2dMasksFrom3dMasks(OrigSegRoi, NewSegRoi, 0)


# In[430]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

ExportSegRoi(TrgSegRoi=NewSegRoi, 
             SrcSegFpath=S011_MR4_S9_SegFpath, 
             NamePrefix='Case_1b_s23_to_s22',
             ExportDir=RoiExportDir)


# In[438]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndSegments

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True

NewSegFpath = os.path.join(RoiExportDir, r'20201022_083705_Case_1b_s23_to_s22_from_Seg_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm')

PlotDicomsAndSegments(DicomDir=S011_MR4_S9_DicomDir, 
                      SegFpath=NewSegFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# # USE CASE 2a:  Direct copy contour across different series (but same FOR, IOP, IPP, ST and PS; Series 9 to Series 5)

# In[454]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import DirectCopyContourAcrossSeries

NewRtsRoi = DirectCopyContourAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                          SrcRtsFpath=S011_MR4_S9_RtsFpath, 
                                          FromSliceNum=27, 
                                          TrgDcmDir=S011_MR4_S5_DicomDir, 
                                          TrgRtsFpath=S011_MR4_S5_RtsFpath, 
                                          ToSliceNum=28, 
                                          LogToConsole=False)


# In[455]:


# Read in the original RTS for S5:
S011_MR4_S5_Rts = pydicom.dcmread(S011_MR4_S5_RtsFpath)

S011_MR4_S5_Rts


# In[456]:


NewRtsRoi


# In[471]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportRtsRoi

ExportRtsRoi(TrgRtsRoi=NewRtsRoi, 
             SrcRtsFpath=S011_MR4_S5_RtsFpath, 
             NamePrefix='Case_2a_DirectCopy_s23_to_s22',
             ExportDir=RoiExportDir)


# ### Plot the DICOMs with the original RTS contours:

# In[469]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndContours

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True


PlotDicomsAndContours(DicomDir=S011_MR4_S5_DicomDir, 
                      RtsFpath=S011_MR4_S5_RtsFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# ### Plot the DICOMs with the new RTS contours:

# In[472]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndContours

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True

# Filepath to the new RTS ROI:
NewRtsFname = r"20201022_112804_Case_2a_DirectCopy_s23_to_s22_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm"
NewRtsFpath = os.path.join(RoiExportDir, NewRtsFname)

PlotDicomsAndContours(DicomDir=S011_MR4_S5_DicomDir, 
                      RtsFpath=NewRtsFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# ### The new RTS was uploaded to XNAT using the RoiUploadAssistant tool, imported in the OHIF-Viewer and verified that they are as expected.

# # USE CASE 2b:  Relationship-preserving copy contour across different series (but same FOR, IOP, IPP, ST and PS; Series 9 to Series 5)

# In[523]:


S011_MR4_S9_Dicoms[0].FrameOfReferenceUID


# In[561]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetImageAttributes

MR4_S9_IPPs, MR4_S9_Dirs, MR4_S9_Spacings, MR4_S9_Dims = GetImageAttributes(S011_MR4_S9_DicomDir)
MR4_S5_IPPs, MR4_S5_Dirs, MR4_S5_Spacings, MR4_S5_Dims = GetImageAttributes(S011_MR4_S5_DicomDir)
MR4_S8_IPPs, MR4_S8_Dirs, MR4_S8_Spacings, MR4_S8_Dims = GetImageAttributes(S011_MR4_S8_DicomDir)
MR12_S8_IPPs, MR12_S8_Dirs, MR12_S8_Spacings, MR12_S8_Dims = GetImageAttributes(S011_MR12_S8_DicomDir)


# MR4_S9_IPPs

# MR4_S5_IPPs

# IPPs, Dirs, Spacings, Dims = GetImageAttributes(S011_MR4_S8_DicomDir)
# IPPs

# IPPs, Dirs, Spacings, Dims = GetImageAttributes(S011_MR4_S5_DicomDir)
# IPPs

# Comparisons = []
# 
# for i in range(len(MR4_S9_IPPs)):
#     if MR4_S9_IPPs[i] == MR4_S5_IPPs[i]:
#         Comparisons.append(True)
#         
#     else:
#         Comparisons.append(False)
#         
# Comparisons

# In[553]:


# Scan 9 and Scan 5 fulfil Case #2b:
if MR4_S9_Spacings == MR4_S5_Spacings:
    print(True)
else:
    print(False)


# In[555]:


# Scan 9 and Scan 8 fulfil Case #3:
if MR4_S9_Spacings == MR4_S8_Spacings:
    print(True)
else:
    print(False)


# In[ ]:


# No data available yet to fulfil Case #4...


# In[562]:


# MR4 Scan 9 and MR12 Scan 8 fulfil Case #5:
if MR4_S9_Spacings == MR12_S8_Spacings:
    print(True)
else:
    print(False)


# ### Work on Case 2b:

# In[596]:


MR4_S9_Seg = pydicom.dcmread(S011_MR4_S9_SegFpath)
MR4_S5_Seg = pydicom.dcmread(S011_MR4_S5_SegFpath)

MR4_S9_PixArr = MR4_S9_Seg.pixel_array
MR4_S5_PixArr = MR4_S5_Seg.pixel_array

MR4_S9_Dicoms = GetDicoms(DicomDir=S011_MR4_S9_DicomDir, SortMethod='slices', Debug=False)
MR4_S5_Dicoms = GetDicoms(DicomDir=S011_MR4_S5_DicomDir, SortMethod='slices', Debug=False)

MR4_S9_IPPs, MR4_S9_Dirs, MR4_S9_Spacings, MR4_S9_Dims = GetImageAttributes(S011_MR4_S9_DicomDir)
MR4_S5_IPPs, MR4_S5_Dirs, MR4_S5_Spacings, MR4_S5_Dims = GetImageAttributes(S011_MR4_S5_DicomDir)

#print(f'Shape of MR4 S9 DICOM PixelData = {MR4_S9_Dicoms[0].pixel_array.shape}')
#print(f'Shape of MR4 S5 DICOM PixelData = {MR4_S5_Dicoms[0].pixel_array.shape}')
print(f'Shape of MR4 S9 image      = {MR4_S9_Dims}')
print(f'Voxel dims of MR4 S9 image = {MR4_S9_Spacings}')
print(f'Origin of MR4 S9 image     = {MR4_S9_IPPs[0]}')
print(f'Shape of MR4 S9 SEG PixArr = {MR4_S9_PixArr.shape}')

print('')
print(f'Shape of MR4 S5 image      = {MR4_S5_Dims}')
print(f'Voxel dims of MR4 S5 image = {MR4_S5_Spacings}')
print(f'Origin of MR4 S5 image     = {MR4_S5_IPPs[0]}')
print(f'Shape of MR4 S5 SEG PixArr = {MR4_S5_PixArr.shape}')


# In[597]:


FrameNum = 0

non0 = np.nonzero(MR4_S9_PixArr[FrameNum])

print(f'First non-zero element in MR4_S9_PixArr is at [{FrameNum}, {non0[0][0]}, {non0[1][0]}]')


# In[598]:


for f in range(MR4_S9_PixArr.shape[0]):
    print(f'Frame {f} has shape {MR4_S9_PixArr[f].shape}')


# In[595]:


print(type(MR4_S9_PixArr))
print(MR4_S9_PixArr.dtype)
print(MR4_S9_PixArr.shape)
print(sys.getsizeof(MR4_S9_PixArr))
print('')
print(f'{sys.getsizeof(np.zeros((MR4_S9_Dims[2], MR4_S9_Dims[0], MR4_S9_Dims[1]), dtype=bool))} bytes with dtype=bool')
print(f'{sys.getsizeof(np.zeros((MR4_S9_Dims[2], MR4_S9_Dims[0], MR4_S9_Dims[1]), dtype=int))} bytes with dtype=int')


# ### Proposed workflow:
# 
# #### 1. If dealing with RTSTRUCT, first convert the RTSTRUCT to SEG
# 
# #### 2. Get the relationship between the segmentations and the DICOMs:
#     a. Get the SOPInstanceUID of the Source DICOM images
#     b. Get the ReferencedSOPInstanceUID of the Source segmentations
#     c. Determine the SegFrameNum-to-DicomSliceNum mapping from a and b
# 
# #### 3. Create a new Source 3D binary mask (= "labelmap") that consists of the original 2D masks re-mapped using SegFrameNum-to-DicomSliceNum. Now each "frame" in labelmap corresponds to its DICOM slice (the in-between frames are simply 2D masks of zeros).
# 
# #### 4. Re-sample the Source 3D labelmap  so that it has the voxel dimensions of the Target 3D image. This yields the Target 3D labelmap.
# 
# #### 5. Get a list of the SOPInstanceUID for the DICOM that corresponds to each non-zero 2D mask in Target 3D labelmap.
# 
# #### 6. If the original datatype was RTSTRUCT:
#     a. Convert each non-zero 2D mask in the 3D labelmap to a list of points in the PatientCoordinateSystem
#     b. Create a copy of the Source RTSTRUCT as a template for the Target RTSTRUCT
#     c. Modify the tags in Target RTSTRUCT accordingly (e.g. Assign the ReferencedSOPInstanceUIDs = the SOPInstanceUIDs from step 5) and export it
# 
# #### 7. If the original datatype was SEG:
#     a. Create a copy of the Source SEG as a template for the Target SEG
#     b. Modify the tags in Target SEG accordingly (e.g. Assign the ReferencedSOPInstanceUIDs = the SOPInstanceUID from step 5) and export it
# 
# #### If the FrameOfReferenceUID of the Source and Target scans are different then some additional steps will be needed:  namely registration and transforming of the Source labelmap to the Target domain.
# 

# In[605]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ConvertPixelArrayToLabelMap
from RoiCopyTools import GetPFFGStoDcmInds
from RoiCopyTools import GetDicomSOPuids

MR4_S9_LabMap = ConvertPixelArrayToLabelMap(PixelArray=MR4_S9_PixArr,
                                            NumOfSlices=MR4_S9_Dims[2],
                                            FrameToSliceInds=GetPFFGStoDcmInds(MR4_S9_Seg,
                                                                               GetDicomSOPuids(S011_MR4_S9_DicomDir)))

MR4_S5_LabMap = ConvertPixelArrayToLabelMap(PixelArray=MR4_S5_PixArr,
                                            NumOfSlices=MR4_S5_Dims[2],
                                            FrameToSliceInds=GetPFFGStoDcmInds(MR4_S5_Seg,
                                                                               GetDicomSOPuids(S011_MR4_S5_DicomDir)))

print(f'Shape of MR4 S9 image      = {MR4_S9_Dims}')
print(f'Voxel dims of MR4 S9 image = {MR4_S9_Spacings}')
print(f'Origin of MR4 S9 image     = {MR4_S9_IPPs[0]}')
print(f'Shape of MR4 S9 SEG PixArr = {MR4_S9_PixArr.shape}')
print(f'Shape of MR4 S9 SEG LabMap = {MR4_S9_LabMap.shape}')

print('')
print(f'Shape of MR4 S5 image      = {MR4_S5_Dims}')
print(f'Voxel dims of MR4 S5 image = {MR4_S5_Spacings}')
print(f'Origin of MR4 S5 image     = {MR4_S5_IPPs[0]}')
print(f'Shape of MR4 S5 SEG PixArr = {MR4_S5_PixArr.shape}')
print(f'Shape of MR4 S5 SEG LabMap = {MR4_S5_LabMap.shape}')


# In[619]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries

NewSeg = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                       SrcSegFpath=S011_MR4_S9_SegFpath, 
                                       FromSliceNum=28,
                                       TrgDcmDir=S011_MR4_S5_DicomDir, 
                                       TrgSegFpath=S011_MR4_S5_SegFpath, 
                                       LogToConsole=True)

NewSeg


# # 23/10:  I've made major changes to ModifySegTagVals but I haven't verified that the result is as expected.  I'll also need to ensure that the function still works in all the other functions that call it (or keep them as separate functions?)
# 
# # Also, the segments still need to be converted to contours..

# In[ ]:





# In[ ]:





# # USE CASE 2c:  Direct copy segmentation across different series (but same FOR, IOP, IPP, ST and PS; Series 9 to Series 5)

# In[449]:


# Read in the RTS-to-SEG conversion for S5:
S011_MR4_S5_Seg = pydicom.dcmread(S011_MR4_S5_SegFpath)

#S011_MR4_S5_Seg


# ### Since the SEG file was converted from an RTS change the Series Description and Content Label to reflect this:

# In[450]:


S011_MR4_S5_Seg.SeriesDescription = 'SEG_from_' + S011_MR4_S5_Seg.SeriesDescription

S011_MR4_S5_Seg.SegmentSequence[0].SegmentLabel = 'SEG_from_' + S011_MR4_S5_Seg.SegmentSequence[0].SegmentLabel


# ### Export the original SEG conversion with modified Segment Label:

# In[452]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

ExportSegRoi(TrgSegRoi=S011_MR4_S5_Seg, 
             SrcSegFpath=S011_MR4_S5_SegFpath, 
             NamePrefix='Mod',
             ExportDir=SegRootDir)


# ### Re-define the filepath for S5's SEG to use the modified SEG:

# In[453]:


S011_MR4_S5_SegFpath = os.path.join(SegRootDir, "Mod_Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm")


# ### The RoiUploadAssistant tool uploaded the converted RTS-to-SEG file without issues and I was able to import the SEG in the OHIF-Viewer.

# ### Plot the DICOMs and segments of the unmodified SEG file (= conversion of RTS)

# In[476]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndSegments

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True


PlotDicomsAndSegments(DicomDir=S011_MR4_S5_DicomDir, 
                      SegFpath=S011_MR4_S5_SegFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# In[478]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import DirectCopySegmentAcrossSeries

NewSegRoi = DirectCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                          SrcSegFpath=S011_MR4_S9_SegFpath, 
                                          FromSliceNum=27, 
                                          TrgDcmDir=S011_MR4_S5_DicomDir, 
                                          TrgSegFpath=S011_MR4_S5_SegFpath, 
                                          ToSliceNum=28, 
                                          LogToConsole=False)


# In[511]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import Compare2dMasksFrom3dMasks

OrigSegRoi = pydicom.dcmread(S011_MR4_S5_SegFpath)

Compare2dMasksFrom3dMasks(OrigSegRoi, NewSegRoi, S011_MR4_S5_DicomDir, S011_MR4_S5_DicomDir)


# ### Export the new SEG ROI:

# In[516]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ExportSegRoi

NewSegFpath = ExportSegRoi(TrgSegRoi=NewSegRoi, 
                           SrcSegFpath=S011_MR4_S9_SegFpath, 
                           NamePrefix='Case_2c_DirectCopy',
                           ExportDir=RoiExportDir)


# ### Plot the DICOMs and segments of the new SEG file

# In[520]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotDicomsAndSegments

ExportPlot = False
#ExportPlot = True

LogToConsole = False
#LogToConsole = True


PlotDicomsAndSegments(DicomDir=S011_MR4_S5_DicomDir, 
                      SegFpath=NewSegFpath, 
                      ExportPlot=ExportPlot, 
                      ExportDir=PlotExportDir, 
                      LogToConsole=LogToConsole)


# ### New SEG uploaded to XNAT using RoiUploadAssistant, imported in OHIF-Viewer and results look as expected.

# # USE CASE 2d:  Relationship-preserving copy segmentation across different series (but same FOR, IOP, IPP, ST and PS; Series 9 to Series 5)

# # 23/10:  I've made major changes to ModifySegTagVals but I haven't verified that the result is as expected.  I'll also need to ensure that the function still works in all the other functions that call it (or keep them as separate functions?)

# In[ ]:





# # USE CASE 3d:  Relationship-preserving copy segmentation across different series with different PS and/or ST (and hence IPPs), but same FOR and IOP (Series 9 to Series 6, and Series 6 to Series 9)

# ## Use Case 3d-i) Series 9 to Series 6

# In[650]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries

NewSeg = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                       SrcSegFpath=S011_MR4_S9_SegFpath, 
                                       FromSliceNum=28,
                                       TrgDcmDir=S011_MR4_S6_DicomDir, 
                                       TrgSegFpath='', 
                                       LogToConsole=True)

#NewSeg


# In[636]:


v0 = np.array([0.99726487384986, -0.0222629979573, -0.070477871046,               1.863426e-10, 0.95355613378866, -0.301215370979,               0.07391056342109235, 0.30039150892787814, 0.9509480374756568])

v1 = np.array([0.99726487428073, -0.0222629872367, -0.0704778683357,               -8.2676044e-09, 0.95355613776452, -0.3012153583927,               0.07391055760746114, 0.30039149710160407, 0.9509480416632908])

v0 - v1


# ### Possible classes within SimpleITK for image resampling - see list at top of notebook
#     
# ### Other classes that may be useful for other needs - see list at top of notebook

# In[647]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetSitkImageInfo

GetSitkImageInfo(S011_MR4_S9_DicomDir)


# In[648]:


GetSitkImageInfo(S011_MR4_S6_DicomDir)


# In[63]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries


Interpolation = 'Linear'
#Interpolation = 'BSpline'
#Interpolation = 'NearestNeighbour'

Method = 'ResampleImageFilter'
#Method = 'Resample1'
#Method = 'Resample2'
Method = 'Resample3'

LogToConsole = True
#LogToConsole = False

SrcIm, ResSrcIm,TrgIm = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                      SrcSegFpath=S011_MR4_S9_SegFpath, 
                                      FromSliceNum=28,
                                      TrgDcmDir=S011_MR4_S6_DicomDir, 
                                      TrgSegFpath='',
                                      Method=Method,
                                      Interpolation=Interpolation,
                                      LogToConsole=LogToConsole)

#ResSrcIm


# In[ ]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotSitkImages_Src_ResSrc_Trg

ExportPlot = False
#ExportPlot = True

LogToConsole = True
LogToConsole = False

dpi=80

PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcIm, 
                              ResSrcIm=ResSrcIm, 
                              TrgIm=TrgIm, 
                              Method=Method,
                              Interpolation=Interpolation,
                              LogToConsole=LogToConsole,
                              ExportPlot=ExportPlot,
                              ExportDir=PlotExportDir,
                              dpi=dpi)


# In[42]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetExtremePoints

GetExtremePoints(SrcIm)


# In[41]:


GetExtremePoints(TrgIm)


# In[3]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import GetImageAttributes
from RoiCopyTools import GetSitkImageInfo
from RoiCopyTools import ImportSitkImage

SrcIPPs, SrcDirs, SrcSpacings, SrcDims = GetImageAttributes(S011_MR4_S9_DicomDir)
TrgIPPs, TrgDirs, TrgSpacings, TrgDims = GetImageAttributes(S011_MR4_S6_DicomDir)

SrcSitkIm = ImportSitkImage(S011_MR4_S9_DicomDir)
TrgSitkIm = ImportSitkImage(S011_MR4_S6_DicomDir)


# In[97]:


import DicomTools
importlib.reload(DicomTools)
from DicomTools import GetImageAttributes
#from DicomTools import GetImageInfo
from DicomTools import ImportImage
from DicomTools import ImportDicoms
from DicomTools import GetDicomFpaths
from DicomTools import GetDicomSOPuids

SrcDicoms = ImportDicoms(S011_MR4_S9_DicomDir)
TrgDicoms = ImportDicoms(S011_MR4_S6_DicomDir)

print(SrcDims)
print(TrgDims)
print(SrcSitkIm.GetSize())
print(TrgSitkIm.GetSize())
print(SrcDicoms[0].Columns, SrcDicoms[0].Rows)
print(TrgDicoms[0].Columns, TrgDicoms[0].Rows)


# In[69]:


print(SrcSitkIm.GetSpacing())
print(TrgSitkIm.GetSpacing())
#print(SrcSitkIm.GetDimension())
#print(TrgSitkIm.GetDimension())
print(SrcIPPs[1][2] - SrcIPPs[0][2])
print(TrgIPPs[1][2] - TrgIPPs[0][2])
print(np.array(SrcIPPs[1]) - np.array(SrcIPPs[0]))
print(np.array(TrgIPPs[1]) - np.array(TrgIPPs[0]))


# In[70]:


SrcIPPs


# In[108]:


SrcSitkIPPs = []

for i in range(SrcDims[2]):
    SrcSitkIPPs.append(list(SrcSitkIm.TransformIndexToPhysicalPoint((0, 0, i))))
    
SrcSitkIPPs


# In[71]:


TrgIPPs


# In[107]:


TrgSitkIPPs = []

for i in range(TrgDims[2]):
    TrgSitkIPPs.append(list(TrgSitkIm.TransformIndexToPhysicalPoint((0, 0, i))))
    
TrgSitkIPPs


# In[75]:


np.diff(np.array(SrcIPPs), axis=0)[:,2]


# In[109]:


np.diff(np.array(SrcSitkIPPs), axis=0)[:,2]


# In[76]:


np.diff(np.array(TrgIPPs), axis=0)[:,2]


# In[110]:


np.diff(np.array(TrgSitkIPPs), axis=0)[:,2]


# In[92]:


SrcSOPs = GetDicomSOPuids(S011_MR4_S9_DicomDir)
TrgSOPs = GetDicomSOPuids(S011_MR4_S6_DicomDir)

SrcSOPs


# In[93]:


TrgSOPs


# In[95]:


print(f'Num of items in TrgSOPs = {len(TrgSOPs)}')
print(f'Num of unique items in TrgSOPs = {len(set(TrgSOPs))}')


# In[100]:


SrcDicomFpaths = GetDicomFpaths(S011_MR4_S9_DicomDir, SortMethod='slices', LogToConsole=True)
TrgDicomFpaths = GetDicomFpaths(S011_MR4_S6_DicomDir, SortMethod='slices', LogToConsole=True)


# In[99]:


SrcDicomFpaths


# In[98]:


TrgDicomFpaths


# In[103]:


for i in range(len(TrgDicomFpaths)):
    dicom = pydicom.read_file(TrgDicomFpaths[i])
    
    IPP = dicom.ImagePositionPatient
    
    SLoc = dicom.SliceLocation
    
    SThick = dicom.SliceThickness
    
    if i > 0:
        previousdicom = pydicom.read_file(TrgDicomFpaths[i - 1])
        
        previousIPP = previousdicom.ImagePositionPatient
        
        dIPPz = float(IPP[2]) - float(previousIPP[2])
    else:
        dIPPz = 0.0
        
        
    print(f'file {i}: IPPz = {IPP[2]}, dIPPz = {dIPPz}',
          f'SLoc = {SLoc}, SThick = {SThick}')


# In[114]:


TrgExtentZ = float(TrgIPPs[-1][2]) - float(TrgIPPs[0][2])

TrgSitkExtentZ = TrgSitkIPPs[-1][2] - TrgSitkIPPs[0][2]

print(TrgExtentZ)
print(TrgSitkExtentZ)
print('')
print(TrgExtentZ/len(TrgIPPs))
print(TrgSitkExtentZ/len(TrgSitkIPPs))
print('')
print(7.13211059570314/3)


# In[6]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import ResampleImage
from RoiCopyTools import PlotSitkImages_Src_ResSrc_Trg

Interpolation = 'linear'

ExportPlot = False

LogToConsole = False

ResImage = ResampleImage(Image=SrcSitkIm, RefImage=TrgSitkIm, Interpolation=Interpolation)

ResImage


# PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcSitkIm, ResSrcIm=ResImage, TrgIm=TrgSitkIm, 
#                               Interpolation=Interpolation, ExportPlot=ExportPlot, 
#                               ExportDir=PlotExportDir, LogToConsole=LogToConsole)

# In[116]:


SitkReader = sitk.ImageSeriesReader()
SitkNames = SitkReader.GetGDCMSeriesFileNames(S011_MR4_S9_DicomDir)

SitkReader.SetFileNames(SitkNames)
SrcSitkIm2 = SitkReader.Execute()
SrcSitkIm2 = sitk.Cast(SrcSitkIm2, sitk.sitkFloat32)

SrcSitk2IPPs = []

for i in range(SrcSitkIm2.GetSize()[2]):
    SrcSitk2IPPs.append(list(SrcSitkIm2.TransformIndexToPhysicalPoint((0, 0, i))))
    
SrcSitk2IPPs


# In[118]:


SrcIPPs


# In[115]:


SitkReader = sitk.ImageSeriesReader()
SitkNames = SitkReader.GetGDCMSeriesFileNames(S011_MR4_S6_DicomDir)

SitkReader.SetFileNames(SitkNames)
TrgSitkIm2 = SitkReader.Execute()
TrgSitkIm2 = sitk.Cast(TrgSitkIm2, sitk.sitkFloat32)

TrgSitk2IPPs = []

for i in range(TrgSitkIm2.GetSize()[2]):
    TrgSitk2IPPs.append(list(TrgSitkIm2.TransformIndexToPhysicalPoint((0, 0, i))))
    
TrgSitk2IPPs


# In[120]:


TrgSitkIPPs


# In[119]:


TrgIPPs


# In[17]:


S011_MR4_DicomDirs = [S011_MR4_S3_DicomDir, S011_MR4_S4_DicomDir, S011_MR4_S5_DicomDir,
                      S011_MR4_S6_DicomDir, S011_MR4_S7_DicomDir, S011_MR4_S8_DicomDir,
                      S011_MR4_S9_DicomDir]

SeriesNums = [3, 4, 5, 6, 7, 8, 9]


# # Check the Image attributes for all series for MR4

# In[19]:


import DicomTools
importlib.reload(DicomTools)
from DicomTools import GetImageAttributes

#Package = 'sitk'
LogToConsole = True
LogToConsole = False

for i in range(len(S011_MR4_DicomDirs)):
    print(f'\n\n*Series {SeriesNums[i]}*:\n')
    
    # Get the attributes using SimpleITK:
    Size, Spacings, SliceThick, Positions, Dirs = GetImageAttributes(S011_MR4_DicomDirs[i], 'sitk', LogToConsole)
    
    print(f'Size = {Size}\n')
    print(f'Spacings = {Spacings}\n')
    
    #IPPsZ = [Positions[i][2] for i in range(Size[2])]
    
    #print(f'\nIPPsZ = {IPPsZ}\n')
    
    # Get the attributes using Pydicom:
    Size, Spacings, SliceThick, Positions, Dirs = GetImageAttributes(S011_MR4_DicomDirs[i], 'pydicom', LogToConsole)
    
    print(f'SliceThickness = {SliceThick}\n')
    
    


# ### USE CASE 3d:  Relationship-preserving copy segmentation across different series with different PS and/or ST (and hence IPPs), but same FOR and IOP (Series 9 to Series 3, and Series 3 to Series 9)
# 
# ## Note change from Series 6 to Series 3 since Series 6 has repeated IPPs.

# ### Use Case 3d-i) Series 9 to Series 3

# In[253]:


import DicomTools
importlib.reload(DicomTools)
from DicomTools import GetImageAttributes

# Get the Image Attributes for Source and Target using SimpleITK:
SrcSitkSize, SrcSitkSpacing,SrcSitkIPP, SrcSitkDir = GetImageAttributes(DicomDir=S011_MR4_S9_DicomDir,
                                            Package='sitk')

TrgSitkSize, TrgSitkSpacing,TrgSitkIPP, TrgSitkDir = GetImageAttributes(DicomDir=S011_MR4_S3_DicomDir,
                                            Package='sitk')

print('Spacings from SimpleITK:')
print(f'   SrcSitkSpacing = {SrcSitkSpacing}')
print(f'   TrgSitkSpacing = {TrgSitkSpacing}\n\n')


# Get the image attributes using Pydicom:
SrcPydiSize, SrcPydiSpacing, SrcPydiIPP, SrcPydiDir = GetImageAttributes(S011_MR4_S9_DicomDir)

TrgPydiSize, TrgPydiSpacing, TrgPydiIPP, TrgPydiDir = GetImageAttributes(S011_MR4_S3_DicomDir)

print('\n\nSpacings from Pydicom:')
print(f'   SrcPydiSpacing = {SrcPydiSpacing}')
print(f'   TrgPydiSpacing = {TrgPydiSpacing}\n\n')


# In[255]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries


Interpolation = 'Linear'
#Interpolation = 'BSpline'
#Interpolation = 'NearestNeighbour'

Method = 'ResampleImageFilter_Identity' # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'ResampleImageFilter_Affine'   # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'Resample2_Origins'            # 2nd method, aligned by origins # <-- poor result
#Method = 'Resample2_Centres'            # 2nd method, aligned by centres # <-- poor result
#Method = 'Resample3_Origins'            # 3rd method, aligned by origins # <-- very poor result
#Method = 'Resample3_Centres'            # 3rd method, aligned by centres # <-- very poor result

LogToConsole = True
#LogToConsole = False

SrcIm, ResSrcIm,TrgIm = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                      SrcSegFpath=S011_MR4_S9_SegFpath, 
                                      FromSliceNum=28,
                                      TrgDcmDir=S011_MR4_S3_DicomDir, 
                                      TrgSegFpath=S011_MR4_S3_SegFpath,
                                      Method=Method,
                                      Interpolation=Interpolation,
                                      LogToConsole=LogToConsole)

#ResSrcIm


# In[256]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotSitkImages_Src_ResSrc_Trg

SrcImLabel = 'S9'
TrgImLabel = 'S3'

ExportPlot = False
ExportPlot = True

LogToConsole = True
LogToConsole = False

dpi=80

PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcIm, 
                              ResSrcIm=ResSrcIm, 
                              TrgIm=TrgIm, 
                              Method=Method,
                              Interpolation=Interpolation,
                              LogToConsole=LogToConsole,
                              ExportPlot=ExportPlot,
                              ExportDir=PlotExportDir,
                              SrcImLabel=SrcImLabel,
                              TrgImLabel=TrgImLabel,
                              dpi=dpi)


# ### Now resample labelmap

# In[5]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries


Interpolation = 'Linear'
#Interpolation = 'BSpline'
Interpolation = 'NearestNeighbour'

Method = 'ResampleImageFilter_Identity' # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'ResampleImageFilter_Affine'   # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'Resample2_Origins'            # 2nd method, aligned by origins # <-- poor result
#Method = 'Resample2_Centres'            # 2nd method, aligned by centres # <-- poor result
#Method = 'Resample3_Origins'            # 3rd method, aligned by origins # <-- very poor result
#Method = 'Resample3_Centres'            # 3rd method, aligned by centres # <-- very poor result

LogToConsole = True
LogToConsole = False

SrcLabMapIm, ResSrcLabMapIm,TrgLabMapIm, SrcLabMapToCopyIm,ResSrcLabMapToCopyIm = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S9_DicomDir, 
                                                     SrcSegFpath=S011_MR4_S9_SegFpath, 
                                                     FromSliceNum=26,
                                                     TrgDcmDir=S011_MR4_S3_DicomDir, 
                                                     TrgSegFpath=S011_MR4_S3_SegFpath,
                                                     Method=Method,
                                                     Interpolation=Interpolation,
                                                     LogToConsole=LogToConsole)

#ResSrcLabMapIm


# In[7]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import PlotSitkImages_Src_ResSrc_Trg

SrcImLabel = 'S9'
SrcImLabel = 'S9_OnlySegmentToCopy'
TrgImLabel = 'S3'


ExportPlot = False
ExportPlot = True

LogToConsole = True
LogToConsole = False

dpi=80

if False:
    PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcLabMapIm, 
                                  ResSrcIm=ResSrcLabMapIm, 
                                  TrgIm=TrgLabMapIm, 
                                  Method=Method,
                                  Interpolation=Interpolation,
                                  LogToConsole=LogToConsole,
                                  ExportPlot=ExportPlot,
                                  ExportDir=PlotExportDir,
                                  SrcImLabel=SrcImLabel,
                                  TrgImLabel=TrgImLabel,
                                  dpi=dpi)

PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcLabMapToCopyIm, 
                              ResSrcIm=ResSrcLabMapToCopyIm, 
                              TrgIm=TrgLabMapIm, 
                              Method=Method,
                              Interpolation=Interpolation,
                              LogToConsole=LogToConsole,
                              ExportPlot=ExportPlot,
                              ExportDir=PlotExportDir,
                              SrcImLabel=SrcImLabel,
                              TrgImLabel=TrgImLabel,
                              dpi=dpi)


# In[286]:


Interpolation


# In[41]:


print(sitk.GetArrayFromImage(SrcLabMapIm).max())
print(sitk.GetArrayFromImage(ResSrcLabMapIm).max())
print(sitk.GetArrayFromImage(TrgLabMapIm).max())


# In[42]:


nda = sitk.GetArrayFromImage(ResSrcLabMapIm)

for i in range(len(nda)):
    print(f'Slice no {i} has max = {nda[:,:,i].max()}')


# ### Should I be using BinaryImageToLabelMap filter?
# 
# https://github.com/SimpleITK/SimpleITK/issues/317
# 
# ### Look into LabelImageToLabelMapFilter, LabelMapContourOverlayImageFilter, LabelMapMaskImageFilter, LabelMapOverlayImageFilter, LabelMapToBinary, LabelMapToLabel, RegionOfInterest, RelabelLabelMap, Resample, Slice, Warp, DemonsRegistrationFilter, DICOMOrientImageFilter
# 
# https://simpleitk.org/doxygen/latest/html/namespaceitk_1_1simple.html#a5900764d7069ea24b3829b3cd23998e4

# In[13]:


#nda2list = sitk.GetArrayFromImage(SrcLabMapIm).tolist()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### Use Case 3d-ii) Series 3 to Series 9

# In[10]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries


#Interpolation = 'Linear'
#Interpolation = 'BSpline'
#Interpolation = 'NearestNeighbour'

Method = 'ResampleImageFilter_Identity' # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'ResampleImageFilter_Affine'   # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'Resample2_Origins'            # 2nd method, aligned by origins # <-- poor result
#Method = 'Resample2_Centres'            # 2nd method, aligned by centres # <-- poor result
#Method = 'Resample3_Origins'            # 3rd method, aligned by origins # <-- very poor result
#Method = 'Resample3_Centres'            # 3rd method, aligned by centres # <-- very poor result

LogToConsole = True
#LogToConsole = False

SrcIm, ResSrcIm, TrgIm,SrcLabMapIm, ResSrcLabMapIm,SrcLabMapToCopyIm, ResSrcLabMapToCopyIm,TrgLabMapIm = MappedCopySegmentAcrossSeries(SrcDcmDir=S011_MR4_S3_DicomDir, 
                                            SrcSegFpath=S011_MR4_S3_SegFpath, 
                                            FromSliceNum=11,
                                            TrgDcmDir=S011_MR4_S9_DicomDir, 
                                            TrgSegFpath=S011_MR4_S9_SegFpath,
                                            Method=Method,
                                            #Interpolation=Interpolation,
                                            LogToConsole=LogToConsole)

#ResSrcIm


# In[18]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import Plot_Src_ResSrc_Trg_Images

SrcImLabel = 'S3'
TrgImLabel = 'S9'

# Chose what to plot:
ImageType = 'Images'
ImageType = 'Labelmaps_w_all_Src_segs'
ImageType = 'Labelmaps_w_only_to_copy_Src_seg'

ExportPlot = False
ExportPlot = True

LogToConsole = True
LogToConsole = False

dpi=80

if ImageType == 'Images':
    SrcIm = SrcIm
    ResSrcIm = ResSrcIm
    TrgIm = TrgIm
    
elif ImageType == 'Labelmaps_w_all_Src_segs':
    SrcIm = SrcLabMapIm
    ResSrcIm = ResSrcLabMapIm
    TrgIm = TrgLabMapIm
    
else:
    SrcIm = SrcLabMapToCopyIm
    ResSrcIm = ResSrcLabMapToCopyIm
    TrgIm = TrgLabMapIm
    

Plot_Src_ResSrc_Trg_Images(SrcIm=SrcIm, 
                           ResSrcIm=ResSrcIm, 
                           TrgIm=TrgIm, 
                           SrcImLabel=SrcImLabel,
                           TrgImLabel=TrgImLabel,
                           ImageType=ImageType,
                           LogToConsole=LogToConsole,
                           ExportPlot=ExportPlot,
                           ExportDir=PlotExportDir,
                           dpi=dpi)


# ### Inspect MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm containing two ROIs:

# In[3]:


MR4_S3_B_SegFpath = os.path.join(SegRootDir, "MR4_S3_brain_SEG_20201030_133110.dcm")
MR4_S3_TandB_SegFpath = os.path.join(SegRootDir, "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm")

MR4_S9_T_SegFpath = os.path.join(SegRootDir, "MR4_S9_tumour_SEG_20201105_115408.dcm")
MR4_S9_TandB_SegFpath = os.path.join(SegRootDir, "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")

S3_B_Seg = pydicom.dcmread(MR4_S3_B_SegFpath)
S3_TandB_Seg = pydicom.dcmread(MR4_S3_TandB_SegFpath)

S9_T_Seg = pydicom.dcmread(MR4_S9_T_SegFpath)
S9_TandB_Seg = pydicom.dcmread(MR4_S9_TandB_SegFpath)

S3_B_Seg


# In[25]:


S3_TandB_Seg


# In[36]:


S9_T_Seg


# In[37]:


S9_TandB_Seg


# In[64]:


print(S3_B_Seg.pixel_array.shape)
print(S3_TandB_Seg.pixel_array.shape)
print(S9_T_Seg.pixel_array.shape)
print(S9_TandB_Seg.pixel_array.shape)


# In[83]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import MappedCopySegmentAcrossSeries


#Interpolation = 'Linear'
#Interpolation = 'BSpline'
#Interpolation = 'NearestNeighbour'

Method = 'ResampleImageFilter_Identity' # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'ResampleImageFilter_Affine'   # <-- reasonable results but could be better? e.g. slices 1 & 7
#Method = 'Resample2_Origins'            # 2nd method, aligned by origins # <-- poor result
#Method = 'Resample2_Centres'            # 2nd method, aligned by centres # <-- poor result
#Method = 'Resample3_Origins'            # 3rd method, aligned by origins # <-- very poor result
#Method = 'Resample3_Centres'            # 3rd method, aligned by centres # <-- very poor result

LogToConsole = True
LogToConsole = False

if False:
    SrcIm, ResSrcIm, TrgIm,    SrcLabMapIm, ResSrcLabMapIm,    SrcLabMapToCopyIm, ResSrcLabMapToCopyIm,    TrgLabMapIm = MappedCopySegmentAcrossSeries(SrcDcmDir=MR4_S3_DicomDir, 
                                                SrcSegFpath=MR4_S3_TandB_SegFpath, 
                                                FromSliceNum=11,
                                                FromSegLabel='tumour',
                                                TrgDcmDir=MR4_S9_DicomDir, 
                                                TrgSegFpath=MR4_S9_TandB_SegFpath,
                                                ToSegLabel='tumour',
                                                Method=Method,
                                                #Interpolation=Interpolation,
                                                LogToConsole=LogToConsole)

# Note change to ___LabMapIms:
SrcIm, ResSrcIm, TrgIm,SrcLabMapIms, ResSrcLabMapIms,SrcLabMapToCopyIm, ResSrcLabMapToCopyIm,TrgLabMapIms, NewTrgLabMapIms = MappedCopySegmentAcrossSeries(SrcDcmDir=MR4_S3_DicomDir, 
                                                              SrcSegFpath=MR4_S3_TandB_SegFpath, 
                                                              FromSliceNum=11,
                                                              FromSegLabel='tumour',
                                                              TrgDcmDir=MR4_S9_DicomDir, 
                                                              TrgSegFpath=MR4_S9_TandB_SegFpath,
                                                              ToSegLabel='tumour',
                                                              Method=Method,
                                                              #Interpolation=Interpolation,
                                                              LogToConsole=LogToConsole)


# In[58]:


sitk.GetArrayFromImage(SrcLabMapIms[0]).shape


# In[101]:


import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import Plot_Src_ResSrc_Trg_Images
from RoiCopyTools import Plot_Src_ResSrc_Trg_NewTrg_Images

SrcImLabel = 'S3'
TrgImLabel = 'S9'

# Chose what to plot:
ImageType = 'Images'
ImageType = 'Labelmaps_w_all_Src_segs'
ImageType = 'Labelmaps_w_only_to_copy_Src_seg'

ExportPlot = False
#ExportPlot = True

LogToConsole = True
LogToConsole = False

dpi=80

if ImageType == 'Images':
    Plot_Src_ResSrc_Trg_Images(SrcIm=SrcIm, 
                               ResSrcIm=ResSrcIm, 
                               TrgIm=TrgIm, 
                               SrcImLabel=SrcImLabel,
                               TrgImLabel=TrgImLabel,
                               ImageType=ImageType,
                               LogToConsole=LogToConsole,
                               ExportPlot=ExportPlot,
                               ExportDir=PlotExportDir,
                               dpi=dpi)
    
    
elif ImageType == 'Labelmaps_w_all_Src_segs':
    #Plot_Src_ResSrc_Trg_Images(SrcIm=SrcLabMapIms, 
    #                           ResSrcIm=ResSrcLabMapIms, 
    #                           TrgIm=TrgLabMapIms, 
    #                           SrcImLabel=SrcImLabel,
    #                           TrgImLabel=TrgImLabel,
    #                           ImageType=ImageType,
    #                           LogToConsole=LogToConsole,
    #                           ExportPlot=ExportPlot,
    #                           ExportDir=PlotExportDir,
    #                           dpi=dpi)
    
    Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm=SrcLabMapIms, 
                                      ResSrcIm=ResSrcLabMapIms, 
                                      TrgIm=TrgLabMapIms, 
                                      SrcImLabel=SrcImLabel,
                                      TrgImLabel=TrgImLabel,
                                      NewTrgIm=NewTrgLabMapIms,
                                      ImageType=ImageType,
                                      LogToConsole=LogToConsole,
                                      ExportPlot=ExportPlot,
                                      ExportDir=PlotExportDir,
                                      dpi=dpi)
    
else:
    #Plot_Src_ResSrc_Trg_Images(SrcIm=SrcLabMapToCopyIm, 
    #                           ResSrcIm=ResSrcLabMapToCopyIm, 
    #                           TrgIm=TrgLabMapIms, 
    #                           SrcImLabel=SrcImLabel,
    #                           TrgImLabel=TrgImLabel,
    #                           ImageType=ImageType,
    #                           LogToConsole=LogToConsole,
    #                           ExportPlot=ExportPlot,
    #                           ExportDir=PlotExportDir,
    #                           dpi=dpi)
    
    Plot_Src_ResSrc_Trg_NewTrg_Images(SrcIm=SrcLabMapToCopyIm, 
                                      ResSrcIm=ResSrcLabMapToCopyIm, 
                                      TrgIm=TrgLabMapIms, 
                                      NewTrgIm=NewTrgLabMapIms,
                                      SrcImLabel=SrcImLabel,
                                      TrgImLabel=TrgImLabel,
                                      ImageType=ImageType,
                                      LogToConsole=LogToConsole,
                                      ExportPlot=ExportPlot,
                                      ExportDir=PlotExportDir,
                                      dpi=dpi)


# ### I'm still not seeing what I'm expecting to see despite adding the ROIs as numpy arrays and as sitk images so it leads me to suspect that there's a problem with the inputs created in the main function.

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## Use Case 3d-ii) Series 6 to Series 9

# In[ ]:




