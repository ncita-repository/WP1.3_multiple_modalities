# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:12:23 2020

@author: ctorti
"""


"""
ROI copying tools



Make a Direct or Relationship-preserving copy of a single contour/segmentation 
in a given ROI/segment in a Source RTS/SEG for a given slice to one or more
contours/segmentations in a new (if making a relationship-preserving copy) or 
existing (if making a direct copy) ROI/segment in a new Target RTS/SEG either 
using the original Target RTS/SEG (if one has been provided) or using the 
Source RTS/SEG (if a Target RTS/SEG was not provided) as a template.  

For Direct copies the in-plane coordinates (currently the (x,y) coords) will be
preserved only, whilst for a Relationship-preserving copy, all coordinates are 
spatially preserved.

There are five possible Main Use Cases (indicated by the numbers 1 to 5)
and two Sub-cases (indicated by the letters a or b) that apply for a subset
of the Main Use Cases that are accommodated.  The characters i and ii will 
further split Main Use Case #3 into two: i for a sub-case that considers a
copy from a high to low voxel sampling image, and ii for the sub-case that
considers a copy from a low to a high voxel sampling image.  In all cases a 
single contour/segmentation will be copied from the Source RTS/SEG within the 
ROI/segment collection containing the ROIName/SegmentLabel that matches the
characters in a search string (FromSearchString):

1a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images are the same. 

1b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images are the same. 
    
2a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS.

2b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR, IOP, IPP, and VS. 

3a. Direct copy of the Source contour/segmentation on slice FromSliceNum to a 
    new/existing ROI/segment for slice ToSliceNum in Target. 
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).
    
3b. Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR and IOP but have 
    different VS (they will likely also have different IPP).

4.  Relationship-preserving copy of the Source contour/segmentation on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have the same FOR but have different IOP 
    and IPP, and possibly different VS.
   
5.  Relationship-preserving copy of the Source contour/contour on slice 
    FromSliceNum to a new ROI/segment for any number of slices in Target.
    The Source and Target images have different FOR, likely different IPP and 
    IOP, and possibly different VS.
    
FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
IOP = ImageOrientationPatient, VS = Voxel spacings
"""




"""
******************************************************************************
******************************************************************************
DEFINE FILE PATHS TO SOURCE AND TARGET (IF APPLICABLE) RTS/SEG FILES, 
TO SOURCE AND TARGET DICOM DIRECTORIES, AND OTHER INPUTS REQUIRED. 
******************************************************************************
******************************************************************************
"""


def CreateDictOfPaths():
    import os
    
    """ WARNING:
        Don't forget to add new directories/filepaths to the dictionary Dict
        below!
    """
    
    Dict = {}
    
    
    """ Subject TCGA-BB-A5HY in XNAT_TEST """
    
    Sub1_Dir = r"C:\Data\XNAT_TEST\TCGA-BB-A5HY"
    
    Dict['Sub1_Dir'] = Sub1_Dir
    
    """ Session 1 """
    
    Dict['Ses1_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\5_AX T2 BRAIN  PROPELLER\DICOM")
    
    Dict['Ses1_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\6_AX GRE BRAIN\DICOM")
    
    Dict['Ses1_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\7_AX T1 BRAIN\DICOM")
    
    Dict['Ses1_Ser8_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\8_Ax T2 FLAIR BRAIN POST\DICOM")
    
    Dict['Ses1_Ser9_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\9_AX T1 POST BRAIN\DICOM")
    
    
    """ Session 2 """
    
    Dict['Ses2_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\2_Head Routine  5.0  H30s\DICOM")
    
    Dict['Ses2_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\3_Head Routine  5.0  J30s  3\DICOM")
    
    Dict['Ses2_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\4_Head Routine  5.0  H60s\DICOM")
    
    Dict['Ses2_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\5_Head Routine  0.75  H20s\DICOM")

    Dict['Ses2_Ser2_LeftEyeCT_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hkhGnxQ_2kX4w5b39s0\RTSTRUCT",
                   r"AIM_20210201_000731.dcm")
    
    Dict['Ses2_Ser2_LeftEyeCT_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hlDxpyQ_7cJDZotLuP_icr_roiCollectionData\SEG",
                   r"SEG_20210213_140215.dcm") 
    
    Dict['Ses2_Ser2_RightEyeCT_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hkhGnxQ_2NZvk8B9pwA\RTSTRUCT",
                   r"AIM_20210201_000842.dcm") 
    
    Dict['Ses2_Ser2_RightEyeCT_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hlDxpyQ_AqEfVVmmysA_icr_roiCollectionData\SEG",
                   r"SEG_20210213_135857.dcm") 
    
    
    """ Session 3 """
    
    Dict['Ses3_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"6_T2 FLAIR Ax\DICOM")
    
    Dict['Ses3_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"7_Diff Ax\DICOM")
    
    Dict['Ses3_Ser13_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"13_T1 Ax Post fs IAC\DICOM")
    
    Dict['Ses3_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"14_T1 Ax Post Brain\DICOM")
    
    
    """ Session 4 """
    
    Dict['Ses4_Ser10_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans", r"10_T1  AXIAL 2MM POST repeat_S7_DIS3D\DICOM")
    
    Dict['Ses4_Ser11_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\11_T1 3D AX POST_S5_DIS3D\DICOM")
    
    Dict['Ses4_Ser12_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\12_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    
    Dict['Ses4_Ser13_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\13_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    
    Dict['Ses4_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\14_T2_FLAIR 3MM_S2_DIS3D\DICOM")
    
    Dict['Ses4_Ser11_Ventricles_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                    r"RoiCollection_hkhGnxQ_BNths5fNmA\RTSTRUCT",
                     r"AIM_20210131_230845.dcm") 
    
    Dict['Ses4_Ser11_Ventricles_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                    r"RoiCollection_hlDxpyQ_5Dz9OulRH6n_icr_roiCollectionData\SEG",
                     r"SEG_20210213_144855.dcm") 
    
    Dict['Ses4_Ser13_NasalCavity_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                   r"RoiCollection_hkhGnxQ_6daB7S8Eckp\RTSTRUCT",
                   r"AIM_20210131_214821.dcm") 
    
    Dict['Ses4_Ser13_NasalCavity_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                   r"RoiCollection_hlDxpyQ_AJ47Vd6ikeE_icr_roiCollectionData\SEG",
                   r"SEG_20210213_152014.dcm") 
    
    
    """ Session 5"""
    
    Dict['Ses5_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\6_T2 FLAIR Ax\DICOM")
    
    Dict['Ses5_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\7_CISS Ax\DICOM")
    
    Dict['Ses5_Ser12_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\12_T2 Ax fs Post Brain\DICOM")
    
    Dict['Ses5_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\14_T1 Ax Post fs IAC\DICOM")
    
    Dict['Ses5_Ser16_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\16_T1 Ax Post fs IAC repeat\DICOM")
    
    Dict['Ses5_Ser12_LeftEye_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hkhGnxQ_4rFYFKE1r6r\RTSTRUCT",
                   r"AIM_20210131_212509.dcm")
    
    Dict['Ses5_Ser12_LeftEye_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hlDxpyQ_8CVLpZmvRD8_icr_roiCollectionData\SEG",
                   r"SEG_20210213_153745.dcm")
    
    Dict['Ses5_Ser16_Tumour_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hkhGnxQ_79YdSO5eYPo\RTSTRUCT",
                   r"AIM_20210131_211622.dcm")
    
    Dict['Ses5_Ser16_Tumour_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hlDxpyQ_18XInHEwVjJ_icr_roiCollectionData\SEG",
                   r"SEG_20210213_155524.dcm")
    
    
    """ Session 6"""
    
    Dict['Ses6_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\3_Head Routine  5.0  H30f\DICOM")
    
    Dict['Ses6_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Dict['Ses6_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\5_Head Routine  0.75  H20f\DICOM")
    
    
    """ Session 7"""
    
    Dict['Ses7_Ser21_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session7\scans\21_T1  AXIAL 2MM POST_S6_DIS3D\DICOM")
    
    Dict['Ses7_Ser22_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\22_T1 3D AX POST_S5_DIS3D\DICOM")
    
    Dict['Ses7_Ser23_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\23_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    
    Dict['Ses7_Ser24_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\24_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    
    Dict['Ses7_Ser25_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\25_T2_FLAIR 3MM_S2_DIS3D\DICOM")

    
    """ Session 8"""
    
    Dict['Ses8_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\2_Head Routine  5.0  H30f\DICOM")
    
    Dict['Ses8_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Dict['Ses8_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\5_Head Routine  0.75  H20f\DICOM")
    
    
    """ Session 9"""
    
    Dict['Ses9_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\3_T2 FLAIR Ax\DICOM")
    
    Dict['Ses9_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\4_Diff Ax\DICOM")
    
    Dict['Ses9_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\5_Diff Ax_ADC\DICOM")
    
    Dict['Ses9_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\14_T2spc Ax\DICOM")
    
    Dict['Ses9_Ser15_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\15_T1 3D Ax Post\DICOM")
    
    
    """ Session 10"""
    
    Dict['Ses10_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\3\DICOM")
    
    Dict['Ses10_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\6\DICOM")
    
    Dict['Ses10_Ser19_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\19\DICOM")
    
    Dict['Ses10_Ser22_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\22\DICOM")
    
    Dict['Ses10_Ser3_RightEye_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session10", r"AIM_20210131_183753\RTSTRUCT", 
                   r"AIM_20210131_183753.dcm")
    
    Dict['Ses10_Ser3_LowerBrain_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session10", r"AIM_20210131_232127\RTSTRUCT",
                   r"AIM_20210131_232127.dcm")
    
    
    """ Session 11"""
    
    Dict['Ses11_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\2_Head Routine  0.75  H20s\DICOM")
    
    Dict['Ses11_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\3_Head Routine  5.0  H30s\DICOM")
    
    Dict['Ses11_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\5_Head Routine  5.0  J30s  3\DICOM")
    
    Dict['Ses11_Ser9_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\9_TEMPORAL BONES RIGHT  5.0  H30s\DICOM")
    
    Dict['Ses11_Ser10_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\10_TEMPORAL BONES RIGHT  0.6  H30s\DICOM") 
    
    
    
    """ Subject ACRIN-FMISO-Brain-011 in ACRIN-FMISO-Brain """
    
    Sub2_Dir = r"C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011"
    Sub2_RoiDir = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011"
    Sub2_RtsDir = os.path.join(Sub2_RoiDir, "Fixed_RTS_Collections")
    Sub2_SegDir = os.path.join(Sub2_RoiDir, "Fixed_SEG_Collections")
    
    Dict['Sub2_Dir'] = Sub2_Dir
    Dict['Sub2_RoiDir'] = Sub2_RoiDir
    Dict['Sub2_RtsDir'] = Sub2_RtsDir
    Dict['Sub2_SegDir'] = Sub2_SegDir

    MR4_Dir = r"04-10-1960-MRI Brain wwo Contrast-69626"
    # IPP =
    # [1.00, 0.00, 0.07,
    # -0.02, 0.95, 0.30,
    # -0.07, -0.30, 0.95]
    
    MR12_Dir = r"06-11-1961-MRI Brain wwo Contrast-79433"
    # IPP =
    # [0.99, 0.11, 0.07,
    # -0.10, 0.99, -0.03,
    # -0.07, 0.02, 1.00]
    
    Dict['MR4_Dir'] = MR4_Dir
    
    Dict['MR12_Dir'] = MR12_Dir
    
    
    """ DICOM directories for MR4 """
    
    Dict['MR4_Ser3_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"3-T2 TSE AXIAL-07507") 
    # 378, 448, 21
    # 0.51, 0.51, 7.5 mm
    # 5.0 mm
    
    Dict['MR4_Ser4_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"4-T2 AXIAL FLAIR DARK FL-48830") 
    # 416, 512, 21 pix
    # 0.45, 0.45, 5.0 mm
    # 5.0 mm
    
    Dict['MR4_Ser5_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"5-T1 SE AXIAL 3MM-81246")  
    # 208, 256, 50 
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    Dict['MR4_Ser7_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"7-ep2ddiff3scantracep2ADC-22197")
    # 192, 192, 21 pix
    # 1.3, 1.3, 7.5 mm
    # 5.0 mm
    
    Dict['MR4_Ser8_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"8-T1 SE AXIAL POST FS FC-59362") 
    # 212, 256, 30 pix
    # 0.9, 0.9, 5.0 mm
    # 5.0 mm
    
    Dict['MR4_Ser9_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268") 
    # 208, 256, 50 pix
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    
    """ RTS/SEG filepaths for MR4 """
    
    Dict['MR4_Ser3_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm")
    #= os.path.join(Sub2_RtsDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201112_162427.dcm") # tumour &
    # brain are not separate ROIs!
    
    Dict['MR4_Ser3_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    #= os.path.join(Sub2_SegDir, "MR4_S3_brain_SEG_20201030_133110.dcm")
    #= os.path.join(Sub2_SegDir,"MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") # one ROI!
    #= os.path.join(Sub2_SegDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") # two ROIs
    #= os.path.join(Sub2_SegDir, "MR4_S3_tumour_SEG_20201216_085421.dcm") # single frame on slice 10
    
    Dict['MR4_Ser5_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    
    Dict['MR4_Ser5_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S5_tumour_SEG_20201216_084441.dcm") 
    #= os.path.join(Sub2_SegDir, "Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    
    Dict['MR4_Ser8_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    
    Dict['MR4_Ser8_SegFpath']\
    = os.path.join(Sub2_SegDir, "Seg_from_MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    
    Dict['MR4_Ser9_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S9_tumour_and_brain_RTS_AIM_20201116_122858.dcm")
    #= os.path.join(Sub2_RtsDir, "MR4_S9_brain_RTS_AIM_20201112_142559.dcm")
    
    Dict['MR4_Ser9_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")
    #= os.path.join(Sub2_SegDir, "Seg_from_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    #= os.path.join(Sub2_SegDir, "MR4_S9_tumour_SEG_20201105_115408.dcm")
    
    
    
    """ DICOM directories for MR12 """
    
    Dict['MR12_Ser3_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"3-T2 TSE AX IPAT2-08659") 
    
    Dict['MR12_Ser4_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"4-T2 AXIAL FLAIR DARK FL-94212")
    
    Dict['MR12_Ser5_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"5-T1 SE AXIAL-43742")
    
    Dict['MR12_Ser8_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"8-T1 SE AXIAL POST FS FC-81428") 
    
    return Dict






def GetPathInputs(TestNum, UseOrigTrgRoi=False):
    """
    Convert a 3D SimpleITK image to a new pixel type.
    
    Inputs:
    ******                      
        
    TestNum : string
        The character string that denotes the test number to be run.
        
    UseOrigTrgRoi : boolean (optional; False by default)
        If True, if the value for the RTS/SEG filepath key is not None, and if
        a Direct copy is to be made, any existing contours/segmentations in the
        original Target RTS/SEG that match the desired RoiLabelToCopy input to
        CopyRoi() will be preserved.  If False, the value of the filepath key 
        to the original Target RTS/SEG will be over-written with None, so that 
        a new Target RTS/SEG file is created that does not contain any of the 
        original contours/segmentations even if they exist for the ROI of 
        interest.
        

    Outputs:
    *******
    
    Image : SimpleITK image
        The 3D image with modified pixel type.

    
    Note:
    ----
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """
    
    #import os
    #from copy import deepcopy
    
    Dict = CreateDictOfPaths()
    
    """ ******************************************************************* """
    """ TestNums for Testing                                                """
    """ ******************************************************************* """
    
    """ Contours """
    
    if TestNum == 'RR1':
        """ Rts Relationship-preserving copy test 1 """
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Alternative inputs to test other possible combinations:
        #FromRoiLabel = None # copy all ROIs
        #FromSliceNum = 10 # there are no contours on slice 10
        #FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1 # direct copy to next slice
        #ToSliceNum = 0
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'RR2':
        """ Rts Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR3':
        """ Rts Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RR4':
        """ Rts Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'RR5':
        """ Rts Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'RR6':
        """ Rts Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR7':
        """ Rts Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'RR8':
        """ Rts Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'RR9':
        """ Rts Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    #elif TestNum == 'RR10':
    #    """ Rts Relationship-preserving copy test 10 """
    #    
    #    SrcDcmKey = 'Ses2_Ser2_DcmDir'
    #    SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
    #    
    #    TrgDcmKey = 'Ses3_Ser6_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Left eye'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    SrcLabel = 'Ses2_Ser2'
    #    TrgLabel = 'Ses3_Ser6'
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' 3_6'
    #
    #elif TestNum == 'RR11':
    #    """ Rts Relationship-preserving copy test 11 """
    #    
    #    SrcDcmKey = 'Ses2_Ser2_DcmDir'
    #    SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
    #    
    #    TrgDcmKey = 'Ses3_Ser6_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Right eye'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    SrcLabel = 'Ses2_Ser2'
    #    TrgLabel = 'Ses3_Ser6'
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'RR10':
        """ Rts Relationship-preserving copy test 10 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR11':
        """ Rts Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'RR12':
        """ Rts Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR13':
        """ Rts Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR14':
        """ Rts Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR15':
        """ Rts Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
        
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmKey = 'Ses4_Ser13_DcmDir'
    #    SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
    #    
    #    TrgDcmKey = 'MR12_Ser3_DcmDir'
    #    TrgDcmKey = 'MR12_Ser4_DcmDir'
    #    TrgDcmKey = 'MR12_Ser8_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Nasal cavity'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
    #    TxtToAddToRoiLabel = ' MR12 S8'
    
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmKey = 'Ses4_Ser11_DcmDir'
    #    SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
    #    
    #    TrgDcmKey = 'MR12_Ser3_DcmDir'
    #    TrgDcmKey = 'MR12_Ser4_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    FromRoiLabel = 'Ventricles'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
     
    if TestNum == 'RD1':
        """ Rts Direct copy test 1 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser2_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        #FromSliceNum = 6
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'RD2':
        """ Rts Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'RD3':
        """ Rts Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RD4':
        """ Rts Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'RD5':
        """ Rts Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
        
    
    """ Segmentations """
    
    
    if TestNum == 'SR1':
        """ Seg Relationship-preserving copy test 1 """
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy segmentations on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'SR2':
        """ Seg Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR3':
        """ Seg Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SR4':
        """ Seg Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'SR5':
        """ Seg Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'SR6':
        """ Seg Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR7':
        """ Seg Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'SR8':
        """ Seg Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'SR9':
        """ Seg Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    
    elif TestNum == 'SR10':
        """ Seg Relationship-preserving copy test 10 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR11':
        """ Seg Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'SR12':
        """ Seg Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR13':
        """ Seg Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR14':
        """ Seg Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR15':
        """ Seg Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
        
     
    if TestNum == 'SD1':
        """ Seg Direct copy test 1 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser2_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        #FromSliceNum = 6 # <-- no segmentation on slice 6?
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'SD2':
        """ Seg Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'SD3':
        """ Seg Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SD4':
        """ Seg Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'SD5':
        """ Seg Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
    
    
    
    
    
    
    
        """ ************************************************************** """
        """ ************************************************************** """
        """ TestNums for Development                                       """
        """ ************************************************************** """
        """ ************************************************************** """
    
    elif TestNum in ['1R', '1S']:
        """ UseCase 1 (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 29
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
    elif TestNum in ['2aR', '2aS']:
        """ UseCase 2a (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser5_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser5_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser5_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 22
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
    elif TestNum in ['2bR', '2bS']:
        """ UseCase 2b (relationship-preserving copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser5_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser5_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser5_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 27
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['3aiR', '3aiS']:
        """ UseCase 3ai (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser3_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser3_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser3_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 25
        ToSliceNum = 20
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
    elif TestNum in ['3aiiR', '3aiiS']:
        """ UseCase 3aii (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser3_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser3_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser3_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = 26
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
    elif TestNum in ['3biR', '3biS']:
        """ UseCase 3bi (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser3_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser3_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser3_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 26
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['3biiR', '3biiS']:
        """ UseCase 3bii (relationship-preserving copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser3_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser3_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser3_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    elif TestNum in ['5bR', '5bS']:
        """ UseCase 5b (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR12_Ser8_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                raise Exception('A RTS does not yet exist for MR12 Ser8.')
            else:
                raise Exception('A SEG does not yet exist for MR12 Ser8.')
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    
    
    """ Use the keys to select the directories and file paths from Dict. """
    SrcDcmDir = Dict[SrcDcmKey]
    SrcRoiFpath = Dict[SrcRoiKey]
    
    TrgDcmDir = Dict[TrgDcmKey]
    if TrgRoiKey:
        TrgRoiFpath = Dict[TrgRoiKey]
    else:
        TrgRoiFpath = None
    
    
    
    return SrcRoiFpath, FromSliceNum, FromRoiLabel, ToSliceNum, SrcDcmDir,\
           TrgDcmDir, TrgRoiFpath, TxtToAddToRoiLabel, SrcLabel, TrgLabel
        



"""
******************************************************************************
******************************************************************************
CHECK VALIDITY OF INPUTS TO COPYROI
******************************************************************************
******************************************************************************
"""

def CheckValidityOfInputs(SrcRoi, FromRoiLabel, FromSliceNum, ToSliceNum, 
                          TrgRoi=None, LogToConsole=False):
    """
    Establish whether the inputs FromRoiLabel, FromSliceNum, ToSliceNum,
    SrcRoi and TrgRoi are valid combinations. If not raise exception.
    
    Inputs:
    ******
    
    SrcRoi : Pydicom object
        Source RTS/SEG object.
        
    FromRoiLabel : string
        All or part of the Source ROI Name of the ROI containing the contour(s)
        to be copied.
    
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
    
    ToSliceNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If ToSliceNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    
    TrgRoi : Pydicom object (optional; None by default)
        Target RTS/SEG object that the contour(s)/segmentation(s) is to be  
        copied to. TrgRoi is only non-None if a Direct copy of the contour/
        segmentation is to be made and added to an existing RTS/SEG.  
                        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Outputs:
    ******
        
    IsValid : boolean
        True if the combinations of inputs are valid.
    """
    
    if TrgRoi:
        from DicomTools import IsSameModalities
        
        SrcModality = SrcRoi.Modality
        TrgModality = TrgRoi.Modality
        
        if not IsSameModalities(SrcRoi, TrgRoi):
            msg = f"The Source ({SrcModality}) and Target ({TrgModality}) " \
                  + "modalities are different. This is not valid."
            
            raise Exception(msg)
            
            
    if FromRoiLabel == None:
        msg = f'WARNING: FromRoiLabel = {FromRoiLabel}, so all contours/'\
              + 'segmentations in all ROIs/segments in the RTS/SEG are to'\
              + f' be copied. But '
                  
        if FromSliceNum != None and ToSliceNum != None:
            msg += f'FromSliceNum (= {FromSliceNum} != None) and ToSliceNum '\
                   + f'(= {ToSliceNum} != None) are defined. '
            msg += 'This is not a valid combination of inputs.'
            
            raise Exception(msg)
            
        elif FromSliceNum != None:
            msg += f'FromSliceNum (= {FromSliceNum} != None) is defined.'
            msg += 'This is not a valid combination of inputs.'
            
            raise Exception(msg)
            
        elif ToSliceNum != None:
            msg += f'ToSliceNum (= {ToSliceNum} != None) is defined.'
            msg += 'This is not a valid combination of inputs.'
        
            raise Exception(msg)
    
    if FromSliceNum == None and ToSliceNum != None:
        msg = f'WARNING: FromSliceNum = {FromSliceNum}, so all contours/'\
              + 'segmentations in the ROI/segment whose label matches '\
              + f'{FromRoiLabel} are to be copied. But ToSliceNum (= '\
              + f'{ToSliceNum} != None) is defined. This is not a valid '\
              + 'combination of inputs.'
              
        raise Exception(msg)
    
    return True




        
"""
******************************************************************************
******************************************************************************
WHICH USE CASE APPLIES?
******************************************************************************
******************************************************************************
"""

def WhichUseCase(FromSliceNum, FromRoiLabel, ToSliceNum, SrcDcmDir, TrgDcmDir, 
                 ForceRegistration, LogToConsole=False):
    """
    Determine which of 5 Use Cases applies.
    
    Inputs:
    ******
        
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
        
    FromRoiLabel : string
        All or part of the Source ROI Name of the ROI containing the contour(s)
        to be copied.
    
    ToSliceNum : integer or None
        The slice index within the Target DICOM stack corresponding to the 
        contour/segmentation to be copied to (applies only for the case of 
        direct copies of single entities). The index is zero-indexed.
        If ToSliceNum = None, a relationship-preserving copy is to be made, 
        since the slice location(s) where the contour(s)/segmentation(s) will 
        be copied to will not depend on user input.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
                       
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    ForceRegistration : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.

              
    Outputs:
    ******
        
    UseCaseThatApplies : string
        String that denotes the Use Case that applies.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
    """
    
    from DicomTools import GetDicomUids
    from ImageTools import GetImageAttributes
    from GeneralTools import AreListsEqualToWithinEpsilon#, PrintTitle
    
    
    """ Get the DICOM UIDs. """
    SrcStudyuid, SrcSeriesuid, SrcFORuid,\
    SrcSOPuids = GetDicomUids(DicomDir=SrcDcmDir)
    
    TrgStudyuid, TrgSeriesuid, TrgFORuid,\
    TrgSOPuids = GetDicomUids(DicomDir=TrgDcmDir)
        
    if ToSliceNum:
        if ToSliceNum > len(TrgSOPuids) - 1:
            msg = f'The input argument "ToSliceNum" = {ToSliceNum} exceeds '\
                  + f'the number of Target DICOMs ({len(TrgSOPuids)}).'
            
            raise Exception(msg)
    
    
    """ Get the Image Attributes for Source and Target using SimpleITK. """
    SrcSize, SrcSpacing, SrcST, SrcIPPs, SrcDirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=SrcDcmDir, Package='sitk')
    
    TrgSize, TrgSpacing, TrgST, TrgIPPs, TrgDirs,\
    ListOfWarnings = GetImageAttributes(DicomDir=TrgDcmDir, Package='sitk')
    
    
    """ Are the Source and Target part of the same study, series and do they 
    have the same SOP UIDs? """
    if SrcStudyuid == TrgStudyuid:
        SameStudy = True
    else:
        SameStudy = False
        
    if SrcSeriesuid == TrgSeriesuid:
        SameSeries = True
    else:
        SameSeries = False
        
    if SrcFORuid == TrgFORuid:
        SameFOR = True
    else:
        SameFOR = False
        
    if SrcSOPuids == TrgSOPuids:
        SameSOPs = True
    else:
        SameSOPs = False
    
    """ Are the Source and Target spacings, directions and origins equal to  
    within 1e-06? """
    SameSpacings = AreListsEqualToWithinEpsilon(SrcSpacing, TrgSpacing)
    SameDirs = AreListsEqualToWithinEpsilon(SrcDirs, TrgDirs)
    SameOrigins = AreListsEqualToWithinEpsilon(SrcIPPs[0], TrgIPPs[0])
    SameST = AreListsEqualToWithinEpsilon([SrcST], [TrgST])
    
    if SrcSize == TrgSize:
        SameSize = True
    else:
        SameSize = False
    
        
    """ Check which case follows. """
    if SameFOR:
        if SameDirs:
            """ 
            Comment 11/02:
                
            It's not entirely necessary to resample the Source image to the 
            Target domain just because their origins differ. It depends on
            whether the slice in Source coincide with slices in Target.
            
            Even if the slices don't coincide, an approximation can be made as
            follows.  When dealing contours, the in-plane coordinates don't 
            need to be changed, and the destination slice can be determined by 
            matching on the out-of-plane (z) coordinates as best as possible. 
            When dealing with segmentations, the destination slice can be 
            determined by matching on the z components of the IPPs as best as
            possible.  Once the destination slice number is known, the 
            necessary shift in pixels for the in-plane components can be 
            applied accordingly. 
            
            Perhaps the approximation can be justified if the difference 
            between the planes in Source and Target are close to an integer.
            If the difference is large, the Source image should be resampled to
            the Target image.  For now I will treat UseCase 2 data sets as 
            UseCase 3 if their origins are not the same. Further complexity
            (taking into account the error in the above approximation) can be
            considered for future work. """
            if SameSpacings and SameOrigins:
            #if SameSpacings:
                if SameSeries:
                    UseCaseThatApplies = '1' # Direct copy
                else:
                    if isinstance(ToSliceNum, int):
                        UseCaseThatApplies = '2a' # Direct copy
                    else:
                        UseCaseThatApplies = '2b' # Relationship-preserving copy
            else:
                """ Comment 11/02:
                    
                Arguably resampling is only required if the spacings along any
                direction are different when dealing with segmentations.  When
                dealing with contours, resampling is probably unnecessary if
                the spacings are different only in the out-of-plane (z)
                direction. """
                if isinstance(ToSliceNum, int):
                    UseCaseThatApplies = '3a' # Direct copy
                else:
                    UseCaseThatApplies = '3b' # Relationship-preserving copy
        else:
            if isinstance(ToSliceNum, int):
                UseCaseThatApplies = '4a' # Direct copy
            else:
                UseCaseThatApplies = '4b' # Relationship-preserving copy
    else:
        if isinstance(ToSliceNum, int):
            UseCaseThatApplies = '5a' # Direct copy
        else:
            UseCaseThatApplies = '5b' # Relationship-preserving copy
    
        
    """ Print comparisons to the console. """
    if LogToConsole:
        from ImageTools import ImportImage, GetImageExtent
        
        SrcIm = ImportImage(SrcDcmDir)
        TrgIm = ImportImage(TrgDcmDir)
        
        """ Get the image extents. """
        SrcImExtent = GetImageExtent(SrcIm)
        TrgImExtent = GetImageExtent(TrgIm)
        
        """ Are the image extents equal to within 1e-06? """
        SameExtents = AreListsEqualToWithinEpsilon(SrcImExtent, TrgImExtent)
    
        #CompareImageAttributes(SrcDcmDir, TrgDcmDir)
        
        print('\n\n', '-'*120)
        print('Results of running WhichUseCase():')
        print('The Source and Target are in the/have:')
        
        if SameStudy:
            print('   * same study')
        else:
            print('   * different studies')
        
        if SameFOR:
            print('   * same Frame of Reference')
        else:
            print('   * different Frames of Reference')
            
        if SameSeries:
            print('   * same series')
        else:
            print('   * different series')
            
        if SameSOPs:
            print('   * same DICOMs')
        else:
            print('   * different DICOMs')
        
        if SameSize:
            print(f'   * same image size\n      Source/Target: {SrcSize}')
        else:
            print(f'   * different image sizes:\n      Source: {SrcSize}',
                  f'\n      Target: {TrgSize}')
        
        if SameSpacings:
            print(f'   * same voxel sizes\n      Source/Target: {SrcSpacing}')
        else:
            print(f'   * different voxel sizes:\n      Source: {SrcSpacing}',
                  f'\n      Target: {TrgSpacing}')
        
        if SameST:
            print(f'   * same slice thickness\n      Source/Target: {SrcST}')
        else:
            print(f'   * different slice thickness:\n      Source: {SrcST}',
                  f'\n      Target: {TrgST}')
        
        if SameExtents:
            print(f'   * same image extents\n      Source/Target: {SrcImExtent}')
        else:
            print(f'   * different image extents:\n      Source: {SrcImExtent}',
                  f'\n      Target: {TrgImExtent}')
        
        if SameOrigins:
            print(f'   * same origin\n      Source/Target: {SrcIPPs[0]}')
        else:
            print(f'   * different origins:\n      Source: {SrcIPPs[0]}',
                  f'\n      Target: {TrgIPPs[0]}')
            
        if SameDirs:
            print('   * same patient orientation\n',
                  f'      Source/Target: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'                      {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'                      {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]')
        else:
            print('   * different patient orientations:\n',
                  f'      Source: [{SrcDirs[0]}, {SrcDirs[1]}, {SrcDirs[2]},\n',
                  f'               {SrcDirs[3]}, {SrcDirs[4]}, {SrcDirs[5]},\n',
                  f'               {SrcDirs[6]}, {SrcDirs[7]}, {SrcDirs[8]}]\n',
                  f'      Target: [{TrgDirs[0]}, {TrgDirs[1]}, {TrgDirs[2]},\n',
                  f'               {TrgDirs[3]}, {TrgDirs[4]}, {TrgDirs[5]},\n',
                  f'               {TrgDirs[6]}, {TrgDirs[7]}, {TrgDirs[8]}]')
        
    
    if True:#LogToConsole:
        #PrintTitle(f'Case {UseCase} applies.')
        if ForceRegistration and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
            UseCaseToApply = UseCaseThatApplies.replace('3', '5').replace('4', '5')
            
            msg = f'\n*UseCase {UseCaseThatApplies} applies but since '\
                  + f'ForceRegistration = {ForceRegistration}, UseCase '\
                  + f'{UseCaseToApply} will be applied.'
            
            print(msg)
        else:
            UseCaseToApply = UseCaseThatApplies
            
            print(f'\n*UseCase {UseCaseThatApplies} applies.')
            
        print('-'*120)
        
    return UseCaseThatApplies, UseCaseToApply






"""
******************************************************************************
******************************************************************************
COPY A CONTOUR
******************************************************************************
******************************************************************************
"""

def CopyRts(SrcRtsFpath, FromSliceNum, FromRoiLabel, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgRtsFpath=None, ToSliceNum=None, 
            Interpolation='BlurThenLinear', Variance=(1,1,1), 
            ThreshPreRes=0.75, ThreshPostTx=0.05, Tx='affine',
            ForceRegistration=False, AddTxtToRoiLabel='', LogToConsole=False,
            ListOfTimings=[]):
    """
    Note 01/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific contour
            2. All contours in a specific ROI
            3. All contours in all ROIs
        
     
    Inputs:
    ******
    
    SrcRtsFpath : string
        Filepath of the Source RTS file.
    
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        contour to be copied (applies for the case of direct copies of a single 
        contour). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
        
    FromRoiLabel : string
        All or part of the Source ROI Name of the ROI containing the contour(s)
        to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target RTS file that the contour(s) is/are to be copied 
        to. TrgRtsFpath != None only if an existing Target RTS exists and is to
        be added to. An existing RTS will be added to only for Direct copy 
        operations.     
        
    ToSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the contour will be
        copied to (applies only for the case of direct copies of single 
        contours). The index is zero-indexed.
        If ToSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the contour(s) will be copied to will
        not depend on user input.
    
    Interpolation : string
        The type of interpolation to be used if the Source labelmap image(s)
        is/are to be resampled to the Target grid.  Acceptable inputs are:
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    Variance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if a 
        Gaussian blurring + linearly resampling approach is taken for 
        resampling of the Source labelmap image(s).
        
    ThreshPreRes : float (optional; 0.75 by default)
        The threshold level used to perform binary thresholding prior to 
        resampling. Example use is on a labelmap image if a Gaussian blurring 
        + linearly resampling approach is taken to combat aliasing effects.
    
    ThreshPostTx : float (optional; 0.05 by default)
        The threshold level used to perform binary thresholding after 
        registration transformation.  
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    ForceRegistration : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    AddTxtToRoiLabel : string (optional, '' by default)
        String of text to add to the ROI Name of the new Target RTS.
                           
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    LogOfTimings : list of strings (optional; empty by default)
        A list of the time to execute certain tasks during the calling of 
        CopyRts().
          
                  
    Outputs:
    *******
        
    TrgRts : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target RTS object.
        
    ListOfTimings : list of strings
        A list (for each message) of timings for individual operations.
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours will appear displaced w.r.t. 
    the anatomical features highlighted in the Source image.
    """
    
    import importlib
    import RtsTools
    importlib.reload(RtsTools)
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    #from RtsTools import GetPtsByCntByRoi
    from RtsTools import GetRtsDataOfInterest, ProportionOfRoisInExtent
    from RtsTools import CreateRts
    from GeneralTools import UniqueItems, PrintIndsByRoi#, PrintTitle
    from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    from GeneralTools import ShiftPtsByCntByRoi
    #from GeneralTools import ShiftFrame
    #from DicomTools import GetRoiNum
    #from ConversionTools import Points2ContourData
    from ConversionTools import PtsByCntByRoi2CntDataByCntByRoi
    from ImageTools import ImportImage
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print(f'Running CopyRts():')
        
        SrcRts = dcmread(SrcRtsFpath)
        
        SrcRoiLabels = GetRoiLabels(SrcRts)
        
        print(f'\nSrcRoiLabels = {SrcRoiLabels}')
        
        if TrgRtsFpath:
            TrgRts = dcmread(TrgRtsFpath)
            
            TrgRoiLabels = GetRoiLabels(TrgRts)
        
            print(f'\nTrgRoiLabels = {TrgRoiLabels}')
    
    
    #""" Log timing messages. """
    #TimingMsgs = []
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    #""" Check the validity of the combination of inputs. """
    #CheckValidityOfInputs(FromRoiLabel, FromSliceNum, ToSliceNum, LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(FromSliceNum, FromRoiLabel, ToSliceNum, 
    #                              SrcDcmDir, TrgDcmDir, ForceRegistration,
    #                              LogToConsole=True)
    
    #""" Manipulate UseCase if ForceRegistration = True and UseCase was 3 or 4. """
    #if ForceRegistration and not '5' in UseCase:
    #    msg = f'\nUseCase was {UseCase} but will be treated as UseCase '
    #    
    #    UseCase.replace('3', '5').replace('4', '5')
    #    
    #    msg = f'{UseCase}.'
    #    
    #    print(msg)
    
    """ Get the data of interest from the Source RTS. """
    SrcPtsByCntByRoi,\
    SrcC2SindsByRoi = GetRtsDataOfInterest(SrcRtsFpath, FromSliceNum, 
                                           FromRoiLabel, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source RTS file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    """ Raise exception if there is no data of interest from the Source RTS."""
    if not UniqueItems(SrcC2SindsByRoi):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        from RtsTools import GetNumOfContoursByRoi, GetCStoSliceIndsByRoi
        
        SrcRts = dcmread(SrcRtsFpath)
        
        Names = GetRoiLabels(SrcRts)
        N = len(Names)
        
        NumbyRoi = GetNumOfContoursByRoi(SrcRts)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        C2SindsByRoi = GetCStoSliceIndsByRoi(SrcRts, SOPuids)
        
        msg = f"There are no contours on slice {FromSliceNum} in any ROI in "\
              + f"the Source RTS. There are {N} ROIs in the Source RTS:"
        
        for i in range(N):
            msg += f"\nROI {i+1} with name '{Names[i]}' has {NumbyRoi[i]} "\
                   + f"contours that correspond to slices {C2SindsByRoi[i]}" 
        
        
        raise Exception(msg)
    
    """
    Determine whether the contours that make up the ROI(s) intersect with the
    Target image extent (14/01/2021). 
    """
    FracProp = ProportionOfRoisInExtent(PtsByCntByRoi=SrcPtsByCntByRoi,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source ROI(s) intersect the '\
          + 'Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
        
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the points in the ROI lie ',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the points in the ROI lie within the physical extent '\
              + 'of the Target image.'
        raise Exception(msg)
        return None
    
    
    """ Import the 3D images. """
    SrcIm = ImportImage(SrcDcmDir)
    TrgIm = ImportImage(TrgDcmDir)
    
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if UseCase in ['1', '2a', '3a']: # 08/02
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a contour without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        
        """ The contour on FromSliceNum is to be copied to ToSliceNum. Modify 
        the contour-to-slice index (from FromSliceNum) to ToSliceNum. """
        SrcC2SindsByRoi = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcC2SindsByRoi,
                                                   IndToReplace=FromSliceNum, 
                                                   ReplacementInd=ToSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
                
        if TrgRtsFpath:
            """ Direct copy of a contour to an existing RTS. """
            
            """ An existing Target RTS is to be modified.  Since this is a  
            Direct copy, any existing contours in Target with ROI label  
            matching FromRoiLabel are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPtsByCntByRoi,\
            TrgC2SindsByRoi = GetRtsDataOfInterest(TrgRtsFpath, FromSliceNum,
                                                   FromRoiLabel, TrgDcmDir,
                                                   LogToConsole)
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target RTS file for '\
                  + 'existing contours within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            
            if LogToConsole:
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}.')
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The contour will be copied to a different slice location, so the
        z-components of the coordinates will need to be modified. """
        
        """ Shift the out-of-plane (z) components of the points in
        SrcPtsByCntByRoi. """
        SrcPtsByCntByRoi = ZshiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi,
                                               NewZind=ToSliceNum,
                                               DicomDir=SrcDcmDir)
        
        """ Shift the points to account for any difference in the origin of
        SrcIm and TrgIm, and to apply the necessary shift from FromSliceNum to
        ToSliceNum. """
        SrcPtsByCntByRoi,\
        SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                             C2SindsByRoi=SrcC2SindsByRoi, 
                                             SrcImage=SrcIm, 
                                             SrcSliceNum=FromSliceNum, 
                                             TrgImage=TrgIm, 
                                             TrgSliceNum=ToSliceNum, 
                                             RefImage=TrgIm, ShiftInX=False,  
                                             ShiftInY=False, ShiftInZ=True, 
                                             Fractional=False, 
                                             LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcC2SindsByRoi = {SrcC2SindsByRoi}.')
        
        
        if TrgRtsFpath:
            """ An existing Target RTS is to be modified. 
            Concatenate TrgC2SindsByRoi and SrcC2SindsByRoi, and 
            TrgPtsByCntByRoi and SrcPtsByCntByRoi. """
            
            TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + SrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
            
            TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + SrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
        
        else:
            """ A Target RTS was not provided. """
            
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
            if LogToConsole:
                print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
                print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
           
        
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi, and 
            TrgC2SindsByRoi is SrcC2SindsByRoi. """
            
            TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            TrgC2SindsByRoi = deepcopy(SrcC2SindsByRoi)
        
        
        
        
        if False:
            """ Shift the points to account for any difference in the origin (i.e.
            slice number 0) of SrcIm and TrgIm, and to apply the necessary shift in 
            the z components of the points, and to the indices in SrcC2SindsByRoi. """
            SrcPtsByCntByRoi,\
            SrcC2SindsByRoi = ShiftPtsByCntByRoi(PtsByCntByRoi=SrcPtsByCntByRoi, 
                                                 C2SindsByRoi=SrcC2SindsByRoi, 
                                                 SrcImage=SrcIm, 
                                                 SrcSliceNum=0, 
                                                 TrgImage=TrgIm, 
                                                 TrgSliceNum=0, 
                                                 RefImage=TrgIm, ShiftInX=False,  
                                                 ShiftInY=False, ShiftInZ=True, 
                                                 Fractional=False, 
                                                 LogToConsole=LogToConsole)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        #TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1(SrcC2SindsByRoi, 
        #                                                        SrcIm, TrgIm)
        
        TrgC2SindsByRoi = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcC2SindsByRoi, 
                                                                SrcIm, TrgIm)
        
        """ TrgPtsByCntByRoi is simply SrcPtsByCntByRoi. """
        TrgPtsByCntByRoi = deepcopy(SrcPtsByCntByRoi)
            
        """ Convert points to contour data format. """
        TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        
        if LogToConsole:
            print(f'\nSrcC2SindsByRoi = {SrcC2SindsByRoi}')
            print(f'TrgC2SindsByRoi (= shifted SrcC2SindsByRoi) =',
                  f'{SrcC2SindsByRoi}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        from ImageTools import GetImageInfo, ResampleLabmapImByRoi
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PtsByCntByRoi2PixArrByRoi
        from ConversionTools import PixArrByRoi2LabmapImByRoi
        from ConversionTools import PixArrByRoi2PtsByCntByRoi
        #from RtsTools import GetMaskFromContoursForSliceNum
        #from SegTools import ProportionOfSegsInExtent
        #from GeneralTools import NumOfListsAtDepthTwo
        
        
        #""" Import the 3D images. """
        #SrcIm = ImportImage(SrcDcmDir)
        #TrgIm = ImportImage(TrgDcmDir)
        
        """ Convert the RTS data of interest to a list of 3D pixel arrays for 
        each ROI in SrcPtsByCntByRoi. """
        SrcPixArrByRoi = PtsByCntByRoi2PixArrByRoi(SrcPtsByCntByRoi, SrcIm, 
                                                   LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source points to pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Determine whether the pixels that make up SrcPixArrByRoi intersect 
        with the Target image extent. 
        
        Note:
            This is redundant since this was already checked for the points.  
            Only useful as a confirmation that the conversion from points to 
            pixel arrays was correct.
        
        FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrByRoi, 
                                            F2SindsBySeg=SrcC2SindsByRoi,
                                            SrcDicomDir=SrcDcmDir,
                                            TrgDicomDir=TrgDcmDir,
                                            LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'\n{round(FracProp*100, 1)}% of the voxels in the',
                  'segmentation to be copied lie within the physical extent',
                  'of the Target image.')
        
        if FracProp == 0:
            msg = 'None of the voxels in the segmentation to be copied lie '\
                  + 'within the physical extent of the Target image.'
            raise Exception(msg)
            return None
        """
        
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabmapImByRoi = PixArrByRoi2LabmapImByRoi(PixArrByRoi=SrcPixArrByRoi, 
                                                     F2SindsByRoi=SrcC2SindsByRoi, 
                                                     RefIm=SrcIm,
                                                     LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              + 'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
        """ Get the frame-to-slice indices in the order associated with
        SrcLabmapImByRoi.
        Comment:
            Not sure the potential ordering difference between SrcC2SindsByRoi
            and SrcF2SindsByRoi is important for further steps.. """
        SrcF2SindsByRoi = []

        for SrcLabmapIm in SrcLabmapImByRoi:
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(SrcLabmapIm, LogToConsole=False)
            
            SrcF2SindsByRoi.append(F2Sinds)
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        
        """ Resample the source labelmaps.
        Comment:
            Should use SrcF2SindsByRoi rather than SrcC2SindsByRoi?..."""
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabmapImByRoi, ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = ResampleLabmapImByRoi(LabmapImByRoi=SrcLabmapImByRoi, 
                                                   #F2SindsByRoi=SrcC2SindsByRoi,
                                                   F2SindsByRoi=SrcF2SindsByRoi,
                                                   SrcImage=SrcIm, 
                                                   TrgImage=TrgIm,
                                                   Interpolation=Interpolation,
                                                   Variance=Variance,
                                                   ThreshLevel=ThreshPreRes,
                                                   LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
    
        
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from ImageTools import TransformLabmapImByRoi
    
        """ Transform SrcLabmapImByRoi using the transformation that registers
        SrcIm to TrgIm. 
        
        Note:
            Even though the data is being transformed the prefix 'Res' will be
            used since further operations below will need to be applied to
            either resampled or transformed data. """
        
        ResSrcLabmapImByRoi, ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = TransformLabmapImByRoi(LabmapImByRoi=SrcLabmapImByRoi,
                                                    F2SindsByRoi=SrcF2SindsByRoi,
                                                    FixImage=TrgIm, MovImage=SrcIm, 
                                                    Tx=Tx, ThreshLevel=ThreshPostTx,
                                                    LogToConsole=LogToConsole)
        
        
        
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
    
    if UseCaseToApply in ['3a', '4a', '5a']:
        """ Direct copy with resampling/registration and averaging of multi-
        framed pixel arrays.
        
        Note:
        
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled/registered labelmap image).  So if after resampling there are 
        more than one frame, the frames will be averaged or a logical OR 
        operation will be performed.
        """
        
        from GeneralTools import ShiftFramesInPixArrBySeg
        from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrByRoi) # number of ROIs
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            """ Perform averaging or OR operation? """
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                f'Prior to {ReduceUsing} operation')
            
            
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5
    
                ResSrcPixArrByRoi,\
                ResF2SindsByRoi = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                                         F2SindsBySeg=ResSrcF2SindsByRoi,
                                                         MakeBinary=True, 
                                                         BinaryThresh=BinaryThresh,
                                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                #print(f'\n\n\nBefore running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
                    
                ResSrcPixArrByRoi,\
                ResF2SindsByRoi = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi,
                                                       F2SindsBySeg=ResSrcF2SindsByRoi,
                                                       LogToConsole=LogToConsole)
                
                #print(f'\nAfter running OrPixArrByRoi():')
                #print(f'ResSrcPixArrByRoi = {ResSrcPixArrByRoi}')
                #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
                #for r in range(len(ResSrcPixArrByRoi)):
                #    print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
            
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrByRoi[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment.')
                
                PlotPixArrBySeg(ResSrcPixArrByRoi, ResSrcF2SindsByRoi, 
                                'After {ReduceUsing} operation')
        
        
        
        """ Shift the out-of-plane elements in ResSrcPixArrByRoi to account for
        the shift from FromSliceNumber to ToSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be FromSliceNum but instead
        ResSrcF2SindsByRoi[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ FromSliceNum --> ResFromSliceNum following resampling or 
        transformation. """
        ResFromSliceNum = ResSrcF2SindsByRoi[0][0]
        
        if True:#LogToConsole:
            #print(f'\nFromSliceNum = {FromSliceNum} --> ResFromSliceNum =',
            #      f'{ResFromSliceNum}')
            #print(f'ToSliceNum = {ToSliceNum}')
            print(f'\nFromSliceNum = {FromSliceNum} --> ResFromSliceNum =',
                  f'{ResFromSliceNum} following resampling/transformation -->',
                  f'ToSliceNum = {ToSliceNum} for direct copy')
            print('\nResSrcF2SindsByRoi prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsByRoi}')
        
        ResSrcPixArrByRoi,\
        ResSrcF2SindsByRoi = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrByRoi, 
                                                      F2SindsBySeg=ResSrcF2SindsByRoi, 
                                                      ##SrcImage=SrcIm, 
                                                      ##SrcSliceNum=FromSliceNum,
                                                      ##TrgImage=TrgIm,
                                                      #SrcImage=ResSrcLabmapImByRoi[0], 
                                                      #SrcSliceNum=ResFromSliceNum,
                                                      #TrgImage=ResSrcLabmapImByRoi[0],
                                                      SrcImage=TrgIm, 
                                                      SrcSliceNum=ResFromSliceNum,
                                                      TrgImage=TrgIm,
                                                      TrgSliceNum=ToSliceNum, 
                                                      RefImage=TrgIm,
                                                      ShiftInX=False,
                                                      ShiftInY=False,
                                                      ShiftInZ=True,
                                                      Fractional=False,
                                                      LogToConsole=LogToConsole)
        
        if True:#LogToConsole:
            print(f'\nAfter running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
            print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
            for r in range(len(ResSrcPixArrByRoi)):
                print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
                print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
        
        """ Convert ResSrcPixArrByRoi to a list (for each ROI) of a list (for
        each contour) of ContourData, and a list (for each ROI) of a list (for
        each contour) of a list (for each point) of coordinates. """
        
        #print(f'\nBefore running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcF2SindsByRoi = {ResSrcF2SindsByRoi}')
        #print(f'len(ResSrcPixArrByRoi) = {len(ResSrcPixArrByRoi)}')
        #for r in range(len(ResSrcPixArrByRoi)):
        #    #print(f'type(ResSrcPixArrByRoi[{r}]) = {type(ResSrcPixArrByRoi[r])}')
        #    print(f'ResSrcPixArrByRoi[{r}].shape = {ResSrcPixArrByRoi[r].shape}')
        
        """ Is Thresh = 0.5 below too high? """
            
        ResSrcPtsByCntByRoi, ResSrcCntDataByCntByRoi,\
        ResSrcC2SindsByRoi = PixArrByRoi2PtsByCntByRoi(PixArrByRoi=ResSrcPixArrByRoi,
                                                       F2SindsByRoi=ResSrcF2SindsByRoi,
                                                       DicomDir=TrgDcmDir,
                                                       Thresh=0.5,
                                                       LogToConsole=LogToConsole)
        
        #print(f'\nAfter running PixArrByRoi2PtsByCntByRoi():')
        #print(f'ResSrcC2SindsByRoi = {ResSrcC2SindsByRoi}')
        #print(f'len(ResSrcPtsByCntByRoi) = {len(ResSrcPtsByCntByRoi)}')
        #for r in range(len(ResSrcPtsByCntByRoi)):
        #    C = len(ResSrcPtsByCntByRoi[r])
        #    print(f'   There are {C} contours in the {r}^th ROI')
        #    for c in range(C):
        #        P = len(ResSrcPtsByCntByRoi[r][c])
        #        print(f'      There are {P} points in the {c}^th contour')
            
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to convert the Source pixel arrays to points.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if LogToConsole:
            #print('\n\n', '-'*120)
            print(f'\n\nAfter converting ResPixArrToCopy to ContourData and',
                  f'points, ResSrcC2SindsByRoi =')
            PrintIndsByRoi(ResSrcC2SindsByRoi)
            print('')
        
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgRtsFpath:
                """ An existing Target RTS is to be modified. 
                Concatenate TrgC2SindsByRoi and ResSrcC2SindsByRoi, and
                TrgPtsByCntByRoi and ResSrcPtsByCntByRoi. """
                
                TrgC2SindsByRoi = [TrgC2SindsByRoi[r] + ResSrcC2SindsByRoi[r] for r in range(len(TrgC2SindsByRoi))]
                
                TrgPtsByCntByRoi = [TrgPtsByCntByRoi[r] + ResSrcPtsByCntByRoi[r] for r in range(len(TrgPtsByCntByRoi))]
                
            else:
                """ A Target RTS was not provided. """
                
                TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
                
                TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
                
            
            """ Convert points to contour data format. """
            TrgCntDataByCntByRoi = PtsByCntByRoi2CntDataByCntByRoi(TrgPtsByCntByRoi)
        
        else:
            """ UseCase is '3b', '4b' or '5b'. """
            
            TrgPtsByCntByRoi = deepcopy(ResSrcPtsByCntByRoi)
            
            TrgCntDataByCntByRoi = deepcopy(ResSrcCntDataByCntByRoi)
            
            TrgC2SindsByRoi = deepcopy(ResSrcC2SindsByRoi)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'\n\n\nPrior to running CreateRts():')
        print(f'TrgC2SindsByRoi = {TrgC2SindsByRoi}')
        for r in range(len(TrgPtsByCntByRoi)):
            C = len(TrgPtsByCntByRoi[r])
            print(f'   There are {C} in the {r}^th ROI')
            for c in range(len(TrgPtsByCntByRoi[r])):
                P = len(TrgPtsByCntByRoi[r][c])
                print(f'      There are {P} pts in the {c}^th contour')
    
    
    """ Create the Target RTS. """
    TrgRts = CreateRts(SrcRtsFpath, TrgRtsFpath, TrgCntDataByCntByRoi,
                       TrgPtsByCntByRoi, TrgC2SindsByRoi, TrgDcmDir,
                       AddTxtToRoiLabel, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the RTS object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print('\nEnd of CopyRts().')
        print('-'*120)
        
    return TrgRts, ListOfTimings
    #return TrgRts, TrgPtsByCntByRoi, TrgC2SindsByRoi, ListOfTimings






"""
******************************************************************************
******************************************************************************
COPY A SEGMENTATION
******************************************************************************
******************************************************************************
"""

def CopySeg(SrcSegFpath, FromSliceNum, FromSegLabel, SrcDcmDir, TrgDcmDir,
            UseCaseToApply, TrgSegFpath=None, ToSliceNum=None, 
            Interpolation='BlurThenLinear', Variance=(1,1,1), 
            ThreshPreRes=0.75, ThreshPostTx=0.05, Tx='affine',
            ForceRegistration=False, AddTxtToSegLabel='', LogToConsole=False,
            ListOfTimings=[]):
    """
    Note 09/02/2021:
        This function will, depending on the inputs, copy either:
            1. A specific segmentation
            2. All segmentations in a specific segment
            3. All segmentations in all segments
        
     
    Inputs:
    ******
    
    SrcSegFpath : string
        Filepath of the Source SEG file.
    
    FromSliceNum : integer or None
        The slice indeces within the Source DICOM stack corresponding to the
        segmentation to be copied (applies for the case of direct copies of a  
        single segmentation). The index is zero-indexed.
        If FromSliceNum = None, a relationship-preserving copy will be made.
        
    FromSegLabel : string
        All or part of the Source segment Description of the segment containing 
        the segmentation(s) to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
    
    UseCaseToApply : string
        String that denotes the Use Case that will be applied.  This will be
        the same as UseCaseThatApplies unless ForceRegistration is True and
        UseCaseThatApplies is '3a', '3b', '4a' or '4b'.
        
    TrgRtsFpath : string or None (optional; None by default)
        Filepath to the Target SEG file that the segmentation(s) is/are to be 
        copied to. TrgSegFpath != None only if an existing Target SEG exists 
        and is to added to. An existing SEG will be added to only for Direct 
        copy operations.     
        
    ToSliceNum : integer (optional; None by default)
        The slice index within the Target DICOM stack where the segmentation 
        will be copied to (applies only for the case of direct copies of single 
        segmentations). The index is zero-indexed.
        If ToSliceNum = None, a relationship-preserving copy will be made, 
        since the slice location(s) where the segmentation(s) will be copied to 
        will not depend on user input.
        
    Interpolation : string
        The type of interpolation to be used if the Source labelmap image(s)
        is/are to be resampled to the Target grid.  Acceptable inputs are:
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    Variance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if a 
        Gaussian blurring + linearly resampling approach is taken for 
        resampling of the Source labelmap image(s).
        
    ThreshPreRes : float (optional; 0.75 by default)
        The threshold level used to perform binary thresholding prior to 
        resampling. Example use is on a labelmap image if a Gaussian blurring 
        + linearly resampling approach is taken to combat aliasing effects.
    
    ThreshPostTx : float (optional; 0.05 by default)
        The threshold level used to perform binary thresholding after 
        registration transformation.
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    ForceRegistration : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
        
    AddTxtToRoiLabel : string (optional, '' by default)
        String of text to add to the ROI Name of the new Target SEG.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    ListOfTimings : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopySeg().
          
                  
    Outputs:
    *******
        
    TrgSeg : Pydicom object
        The new (if making a Relationship-preserving copy) or modified (if
        making a Direct copy) Target SEG object.
      
    TimingMsgs : list of strings
        A list (for each message) of timings for individual operations.
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target segmentations will appear displaced  
    w.r.t. the anatomical features highlighted in the Source image.
    """
    
    import importlib
    import ConversionTools
    importlib.reload(ConversionTools)
    import GeneralTools
    importlib.reload(GeneralTools)
    import ImageTools
    importlib.reload(ImageTools)
    import SegTools
    importlib.reload(SegTools)
    
    import time
    from copy import deepcopy
    from pydicom import dcmread
    from ImageTools import ImportImage#, GetImageInfo
    from SegTools import GetSegDataOfInterest, ProportionOfSegsInExtent
    from SegTools import CreateSeg
    from GeneralTools import UniqueItems#, PrintIndsByRoi#, PrintTitle
    #from GeneralTools import ZshiftPtsByCntByRoi#, GetPixelShiftBetweenSlices
    #from GeneralTools import ShiftFrame
    
    if LogToConsole:
        from DicomTools import GetRoiLabels
        
        print('\n\n', '-'*120)
        print(f'Running CopySeg():')
        
        SrcSeg = dcmread(SrcSegFpath)
        
        SrcSegLabels = GetRoiLabels(SrcSeg)
        
        print(f'\nSrcSegLabels = {SrcSegLabels}')
        
        if TrgSegFpath:
            TrgSeg = dcmread(TrgSegFpath)
            
            TrgSegLabels = GetRoiLabels(TrgSeg)
        
            print(f'\nTrgSegLabels = {TrgSegLabels}')
    
    
    #""" Log timing messages. """
    #TimingMsgs = []
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    #""" Check the validity of the combination of inputs. """
    #CheckValidityOfInputs(FromRoiLabel, FromSliceNum, ToSliceNum, LogToConsole)
    
    #""" Determine which Use Case to apply. """
    #UseCaseThatApplies,\
    #UseCaseToApply = WhichUseCase(FromSliceNum, FromSegLabel, ToSliceNum, 
    #                              SrcDcmDir, TrgDcmDir, ForceRegistration,
    #                              LogToConsole=True)
    
    #""" Manipulate UseCase if ForceRegistration = True and UseCase was 3 or 4. """
    #if ForceRegistration and not '5' in UseCase:
    #    msg = f'\nUseCase was {UseCase} but will be treated as UseCase '
    #    
    #    UseCase.replace('3', '5').replace('4', '5')
    #    
    #    msg = f'{UseCase}.'
    #    
    #    print(msg)
        
    """ Get the data of interest from the Source SEG. """
    SrcPixArrBySeg,\
    SrcF2SindsBySeg = GetSegDataOfInterest(SrcSegFpath, FromSliceNum, 
                                           FromSegLabel, SrcDcmDir, 
                                           LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to parse the Source SEG file for the ROI.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    """ Raise exception if there is no data of interest from the Source SEG."""
    if not UniqueItems(SrcF2SindsBySeg):
        from DicomTools import GetRoiLabels, GetDicomSOPuids
        from SegTools import GetPFFGStoSliceIndsBySeg
        
        SrcSeg = dcmread(SrcSegFpath)
        
        Labels = GetRoiLabels(SrcSeg)
        L = len(Labels)
        
        SOPuids = GetDicomSOPuids(SrcDcmDir)
        
        F2SindsBySeg = GetPFFGStoSliceIndsBySeg(SrcSeg, SOPuids)
        
        msg = f"There are no segmentations on slice {FromSliceNum} in any "\
              + f"segment in the Source SEG. There are {L} segments in the "\
              + f"Source SEG:"
        
        for i in range(L):
            msg += f"\nSegment {i+1} with label '{Labels[i]}' has "\
                   + f"{len(F2SindsBySeg[i])} segmentations that correspond "\
                   + f"to slices {F2SindsBySeg[i]}" 
                   
        raise Exception(msg)
    
    """
    Determine whether the segmentations that make up the segment(s) intersect 
    with the Target image extent. 
    """
    FracProp = ProportionOfSegsInExtent(PixArrBySeg=SrcPixArrBySeg,
                                        F2SindsBySeg=SrcF2SindsBySeg, 
                                        SrcDicomDir=SrcDcmDir,
                                        TrgDicomDir=TrgDcmDir,
                                        LogToConsole=LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to check whether the Source segment(s) intersect '\
          + 'the Target grid space.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print(f'\n{round(FracProp*100, 1)}% of the indices in the SEG lie',
              'within the physical extent of the Target image.')
    
    if FracProp == 0:
        msg = 'None of the indices in the SEG lie within the physical extent'\
              + 'of the Target image.'
        raise Exception(msg)
        return None
    
    """ Import the 3D images. """
    SrcIm = ImportImage(SrcDcmDir)
    TrgIm = ImportImage(TrgDcmDir)
    
    
    """
    **************************************************************************
    What follows is UseCase-specific.
    **************************************************************************
    """
    
    #if UseCase in ['1', '2a', '3a']: # 08/02
    if UseCaseToApply in ['1', '2a']: # 08/02
        """ Direct copy of a segmentation without resampling or registration
        transformation. """
        
        from GeneralTools import ReplaceIndInC2SindsByRoi
        from GeneralTools import ShiftFramesInPixArrBySeg
        
        """ The segmentation on FromSliceNum is to be copied to ToSliceNum. 
        Modify the frame-to-slice index (from FromSliceNum) to ToSliceNum. """
        SrcF2SindsBySeg = ReplaceIndInC2SindsByRoi(C2SindsByRoi=SrcF2SindsBySeg,
                                                   IndToReplace=FromSliceNum, 
                                                   ReplacementInd=ToSliceNum)
        
        if LogToConsole:
                print(f'\nModfied SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        if TrgSegFpath:
            """ Direct copy of a segmentation to an existing SEG. """
            
            """ An existing Target SEG is to be modified.  Since this is a  
            Direct copy, any existing segmentations in Target with segment 
            label matching FromSegLabel are to be preserved. """
            
            """ Get the data of interest from the Target RTS. """
            TrgPixArrBySeg,\
            TrgF2SindsBySeg = GetSegDataOfInterest(TrgSegFpath, FromSliceNum,
                                                   FromSegLabel, TrgDcmDir,
                                                   LogToConsole)
            
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            msg = f'Took {Dtime} s to parse the Target SEG file for '\
                  + 'existing segmentations within the ROI.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
            
            if LogToConsole:
                print(f'TrgF2SindsByRoi = {TrgF2SindsBySeg}.')
        

    #if UseCase in ['1', '2a']:
        #""" Direct copy without resampling or registration transformation.
        """
        The segmentation will be copied to a different slice location, so the
        z-components of the indices will need to be modified. """
        
        """ Shift the in-plane (x & y) and out-of-plane (z) components of the 
        indices in SrcPixArrBySeg to account for the shift from FromSliceNum to
        ToSliceNum. """
        SrcPixArrBySeg,\
        SrcF2SindsBySeg = ShiftFramesInPixArrBySeg(PixArrBySeg=SrcPixArrBySeg, 
                                                   F2SindsBySeg=SrcF2SindsBySeg, 
                                                   SrcImage=SrcIm, 
                                                   SrcSliceNum=FromSliceNum,
                                                   TrgImage=TrgIm, 
                                                   TrgSliceNum=ToSliceNum, 
                                                   RefImage=TrgIm,
                                                   ShiftInX=True,
                                                   ShiftInY=True,
                                                   ShiftInZ=True,
                                                   Fractional=False,
                                                   LogToConsole=LogToConsole)
        
        if LogToConsole:
                print(f'\nShifted SrcF2SindsBySeg = {SrcF2SindsBySeg}.')
                
        
        if TrgSegFpath:
            """ An existing Target SEG is to be modified. 
            Concatenate TrgF2SindsBySeg and SrcF2SindsBySeg, and 
            TrgPixArrBySeg and SrcPixArrBySeg. """
            
            TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + SrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
            
            TrgPixArrBySeg = [TrgPixArrBySeg[s] + SrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
        
        else:
            """ A Target SEG was not provided. """
            
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            
            if LogToConsole:
                print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
                print(f'TrgF2SindsBySeg = {TrgF2SindsBySeg}')
    
    
    
    if UseCaseToApply == '2b':
        """ Relationship-preserving copy without resampling or registration
        transformation. """
        
        if False:
            """ TrgPixArrBySeg is simply SrcPixArrBySeg, and 
            TrgF2SindsBySeg is SrcF2SindsBySeg. """
            
            TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
            TrgF2SindsBySeg = deepcopy(SrcF2SindsBySeg)
            
        
        #from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v1
        from GeneralTools import GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2
        
        """ Shift the indices in SrcC2SindsByRoi to account for any differences
        in the origin of SrcIm and TrgIm. (new 13/02/21) """
        TrgF2SindsBySeg = GetSrcC2SindsByRoi2TrgC2SindsByRoi_v2(SrcF2SindsBySeg, 
                                                                SrcIm, TrgIm)
        
        
        """ TrgPixArrBySeg is simply SrcPixArrBySeg. """
        TrgPixArrBySeg = deepcopy(SrcPixArrBySeg)
        
        if LogToConsole:
            print(f'\nSrcF2SindsBySeg = {SrcF2SindsBySeg}')
            print(f'TrgF2SindsBySeg (= shifted SrcF2SindsBySeg) =',
                  f'{SrcF2SindsBySeg}.')
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copy with resampling or
        registration transformation. """
        
        #import numpy as np
        #from ImageTools import GaussianBlurImage, BinaryThresholdImage 
        from ImageTools import ResampleLabmapImByRoi#, GetImageInfo
        #from ConversionTools import ConvertImagePixelType
        #from ConversionTools import PixArr2Image, Image2PixArr
        from ConversionTools import PixArrByRoi2LabmapImByRoi
        #from GeneralTools import NumOfListsAtDepthTwo
        
        """ Convert pixel arrays to binary labelmap images. """
        SrcLabmapImBySeg = PixArrByRoi2LabmapImByRoi(PixArrByRoi=SrcPixArrBySeg, 
                                                     F2SindsByRoi=SrcF2SindsBySeg, 
                                                     RefIm=SrcIm,
                                                     LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to generate labelmap images from the Source '\
              'pixel arrays.\n'
        ListOfTimings.append(msg)
        print(f'\n*{msg}')
    
    
    
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b']:
        """ Direct or relationship-preserving copy with resampling. """
        
        
        """ Resample the source labelmaps. """
        
        #Interp = 'NearestNeighbor'
        #Interp = 'LabelGaussian'
        
        ResSrcLabmapImBySeg, ResSrcPixArrBySeg,\
        ResSrcF2SindsBySeg = ResampleLabmapImByRoi(LabmapImByRoi=SrcLabmapImBySeg, 
                                                   #F2SindsByRoi=SrcC2SindsByRoi,
                                                   F2SindsByRoi=SrcF2SindsBySeg,
                                                   SrcImage=SrcIm, 
                                                   TrgImage=TrgIm,
                                                   Interpolation=Interpolation,
                                                   Variance=Variance,
                                                   ThreshLevel=ThreshPreRes,
                                                   LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        msg = f'Took {Dtime} s to resample the Source labelmap images.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
    if UseCaseToApply in ['5a', '5b']:
        """ Direct or relationship-preserving copying with registration
        transformation. """
        
        from ImageTools import TransformLabmapImByRoi
    
        """ Transform SrcLabmapImBySeg using the transformation that registers
        SrcIm to TrgIm. 
        
        Note:
            Even though the data is being transformed the prefix 'Res' will be
            used since further operations below will need to be applied to
            either resampled or transformed data. """
        
        ResSrcLabmapImBySeg, ResSrcPixArrBySeg,\
        ResSrcF2SindsBySeg = TransformLabmapImByRoi(LabmapImByRoi=SrcLabmapImBySeg,
                                                    F2SindsByRoi=SrcF2SindsBySeg,
                                                    FixImage=TrgIm, MovImage=SrcIm, 
                                                    Tx=Tx, ThreshLevel=ThreshPostTx, 
                                                    LogToConsole=LogToConsole)
        
        #times.append(time.time())
        #Dtime = round(times[-1] - times[-2], 1)
        #print(f'\n*Took {Dtime} s to transform the Source labelmap images',
        #      'from the registration transformation.')
        
    
    if UseCaseToApply in ['3a', '4a', '5a']:
        """ Direct copy with resampling/registration and averaging of multi-
        framed pixel arrays.
        
        Note:
        
        For a Direct copy, one contour/segment is always copied to a single
        contour/segmentation (irrespective of how many frames are in the 
        resampled/registered labelmap image).  So if after resampling there are 
        more than one frame, the frames will be averaged or a logical OR 
        operation will be performed.
        """
        
        from GeneralTools import ShiftFramesInPixArrBySeg
        from PlottingTools import PlotPixArrBySeg
        
        R = len(ResSrcPixArrBySeg) # number of segments
        
        """ Get the number of frames in each pixel array. """
        NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
        MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
        
        if MaxNumOfFramesBySeg > 1:
            """ Reduce the number of frames in each pixel array to 1. """
            
            """ Perform averaging or OR operation? """
            #ReduceUsing = 'mean'
            ReduceUsing = 'OR'
            
            if LogToConsole:
                print(f'\nThere are {NumOfFramesBySeg} frames in each',
                      f'ROI/segment. Using {ReduceUsing} operation to reduce',
                      'to single-framed pixel array(s)')
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                f'Prior to {ReduceUsing} operation')
                
            if ReduceUsing == 'mean':
                from GeneralTools import MeanFrameInPixArrBySeg
                
                """ 14/02:  There sems to be a problem with the averaging 
                operation that results in many more than 1 F2Sind and many more
                than 1 contour.. """
                
                """ The threshold used when converting the averaged (non-binary)  
                pixel array to a binary pixel array. """
                BinaryThresh = 0.5 # <-- is this too high?
    
                ResSrcPixArrBySeg,\
                ResF2SindsBySeg = MeanFrameInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                                         F2SindsBySeg=ResSrcF2SindsBySeg,
                                                         MakeBinary=True, 
                                                         BinaryThresh=BinaryThresh,
                                                         LogToConsole=LogToConsole)
            
            else:
                from GeneralTools import OrFrameOfPixArrBySeg
                
                ResSrcPixArrBySeg,\
                ResSrcF2SindsBySeg = OrFrameOfPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg,
                                                          F2SindsBySeg=ResSrcF2SindsBySeg,
                                                          LogToConsole=LogToConsole)
        
            
            if LogToConsole:
                NumOfFramesBySeg = [ResSrcPixArrBySeg[r].shape[0] for r in range(R)]
        
                MaxNumOfFramesBySeg = max(NumOfFramesBySeg)
                
                print(f'\nFollowing {ReduceUsing} operation there are',
                      f'{NumOfFramesBySeg} frames in each ROI/segment, and',
                      f'the new ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}.')
                
                PlotPixArrBySeg(ResSrcPixArrBySeg, ResSrcF2SindsBySeg, 
                                'After {ReduceUsing} operation')
                
                
        
        """ Shift the out-of-plane elements in ResSrcPixArrBySeg to account for
        the shift from FromSliceNumber to ToSliceNum. 
        Note:
            The function below was incorrectly applying in-plane shifts as well
            up until the added inputs ShiftInX / Y / Z was added on 10/02. 
        
        14/02: I think SrcSliceNum shouldn't be FromSliceNum but instead
        ResSrcF2SindsBySeg[0][0], since this is the frame-to-slice number of the
        resampled/transformed labelmap.  The slice number of the original
        (Source) labelmap is not relevant.
        
        Likewise, the inputs SrcImage and TrgImage should both be TrgImage 
        since the resampled labelmap shares the TrgImage grid, and it will be
        shifted to another location, still within the TrgImage grid.
        """
        
        """ FromSliceNum --> ResFromSliceNum following resampling or 
        transformation. """
        ResFromSliceNum = ResSrcF2SindsBySeg[0][0]
        
        if LogToConsole:
            #print(f'\nFromSliceNum = {FromSliceNum} --> ResFromSliceNum =',
            #      f'{ResFromSliceNum}')
            #print(f'ToSliceNum = {ToSliceNum}')
            print(f'\nFromSliceNum = {FromSliceNum} --> ResFromSliceNum =',
                  f'{ResFromSliceNum} following resampling/transformation -->',
                  f'ToSliceNum = {ToSliceNum} for direct copy')
            print('\nResSrcF2SindsBySeg prior to shifting of out-of-plane',
                  f'elements: {ResSrcF2SindsBySeg}')
            
        ResSrcPixArrBySeg,\
        ResSrcF2SindsBySeg = ShiftFramesInPixArrBySeg(PixArrBySeg=ResSrcPixArrBySeg, 
                                                      F2SindsBySeg=ResSrcF2SindsBySeg, 
                                                      SrcImage=TrgIm, 
                                                      SrcSliceNum=ResFromSliceNum,
                                                      TrgImage=TrgIm, 
                                                      TrgSliceNum=ToSliceNum, 
                                                      RefImage=TrgIm,
                                                      ShiftInX=False,
                                                      ShiftInY=False,
                                                      ShiftInZ=True,
                                                      Fractional=False,
                                                      LogToConsole=LogToConsole)
        
        if LogToConsole:
            print(f'\nAfter running ShiftFramesInPixArrBySeg():')
            print(f'ResSrcF2SindsBySeg = {ResSrcF2SindsBySeg}')
            print(f'len(ResSrcPixArrBySeg) = {len(ResSrcPixArrBySeg)}')
            for r in range(len(ResSrcPixArrBySeg)):
                print(f'type(ResSrcPixArrBySeg[{r}]) = {type(ResSrcPixArrBySeg[r])}')
                print(f'ResSrcPixArrBySeg[{r}].shape = {ResSrcPixArrBySeg[r].shape}')
                
    
    if UseCaseToApply in ['3a', '3b', '4a', '4b', '5a', '5b']:
        """ Direct or relationship-preserving copying with resampling or
        registration-transformation. """
    
        if UseCaseToApply in ['3a', '4a', '5a']:
            """ Direct copying with resampling/transformation. """
            
            if TrgSegFpath:
                """ An existing Target SEG is to be modified. 
                Concatenate TrgF2SindsBySeg and ResSrcF2SindsBySeg, and
                TrgPixArrBySeg and ResSrcPixArrBySeg. """
                
                TrgF2SindsBySeg = [TrgF2SindsBySeg[s] + ResSrcF2SindsBySeg[s] for s in range(len(TrgF2SindsBySeg))]
                
                TrgPixArrBySeg = [TrgPixArrBySeg[s] + ResSrcPixArrBySeg[s] for s in range(len(TrgPixArrBySeg))]
                
            else:
                """ A Target SEG was not provided. """
                
                TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
                
                TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
        
        
        else:
            """ UseCaseToApply is '3b', '4b' or '5b'. """
            
            TrgPixArrBySeg = deepcopy(ResSrcPixArrBySeg)
            
            TrgF2SindsBySeg = deepcopy(ResSrcF2SindsBySeg)
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[0], 1)
    msg = f'Took {Dtime} s to copy the ROI(s) from Source to Target.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    #print(f'\n\n\n\n\nTrgPixArrBySeg = {TrgPixArrBySeg}')
    
    
    """ Create/overwrite the Target SEG. """
    
    TrgSeg = CreateSeg(SrcSegFpath, TrgSegFpath, TrgPixArrBySeg, 
                       TrgF2SindsBySeg, TrgDcmDir, FromSegLabel,
                       AddTxtToSegLabel, LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    msg = f'Took {Dtime} s to create the SEG object.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    if LogToConsole:
        print('\nEnd of CopySeg().')
        print('-'*120)
        
    return TrgSeg, ListOfTimings
    #return TrgSeg, TrgPixArrBySeg, TrgF2SindsBySeg, ListOfTimings





"""
******************************************************************************
******************************************************************************
COPY A CONTOUR / SEGMENTATION
******************************************************************************
******************************************************************************
"""


def CopyRoi(SrcRoiFpath, FromSliceNum, FromRoiLabel, SrcDcmDir, TrgDcmDir,
            TrgRoiFpath=None, ToSliceNum=None, Interpolation='BlurThenLinear',
            Variance=(1,1,1), ThreshPreRes=0.75, ThreshPostTx=0.05, Tx='affine',
            ForceRegistration=False, TxtToAddToRoiLabel='', LogToConsole=False,
            ExportLogFiles=False):
    """
    
    Inputs:
    ******
    
    SrcRoiFpath : string
        Filepath of Source RTS/SEG file.
    
    FromSliceNum : integer
        Slice index of the Source DICOM stack corresponding to the contour/
        segmentation to be copied (counting from 0).
        
    FromRoiLabel : string
        All or part of the Source ROIName/SegmentLabel of the ROI/segment 
        containing the contour/segmentation to be copied.
    
    SrcDcmDir : string
        Directory containing the Source DICOMs.
    
    TrgDcmDir : string
        Directory containing the Target DICOMs.
        
    TrgRoiFpath : string or None (optional; None by default)
        Filepath to the Target RTS/SEG file that the contour(s)/segmentation(s) 
        is/are to be copied to. TrgRoiFpath != None only if an existing Target 
        RTS/SEG exists and is to be added to. An existing RTS/SEG will be added 
        to only for Direct copy operations.        
        
    ToSliceNum : integer (optional; None by default)
        Slice index within the Target DICOM stack where the contour/
        segmentation is to be copied to (counting from 0).  This only applies 
        for Direct copies, hence the default value None.
        
    Interpolation : string
        The type of interpolation to be used if the Source labelmap image(s)
        is/are to be resampled to the Target grid.  Acceptable inputs are:
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    Variance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if a 
        Gaussian blurring + linearly resampling approach is taken for 
        resampling of the Source labelmap image(s).
        
    ThreshPreRes : float (optional; 0.75 by default)
        The threshold level used to perform binary thresholding prior to 
        resampling. Example use is on a labelmap image if a Gaussian blurring 
        + linearly resampling approach is taken to combat aliasing effects.
    
    ThreshPostTx : float (optional; 0.05 by default)
        The threshold level used to perform binary thresholding after 
        registration transformation.
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    ForceRegistration : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
        
    TxtToAddToRoiLabel : string (optional, '' by default)
        String of text to pass to CreateRts()/CreateSeg(). The string will be
        appended to the RTS StructureSetLabel or SEG SeriesDescription, and to 
        the filename of the exported (new) Target RTS/SEG.
            
    LogToConsole : boolean (default False)
        If True, some results will be logged to the console.
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
          
                  
    Outputs:
    *******
        
    TrgRoi : Pydicom object
        The new (if making a Relationship-preserving copy) or the modified (if
        making a Direct copy of a contour/segmentation to an existing RTS/SEG) 
        Target RTS/SEG object.
        
    DictOfInputs : dictionary
        A dictionary containing the inputs that were called to CopyRoi().
        
    ListOfInputs : list of strings
        A list containing the inputs that were called to CopyRoi().
        
    ListOfTimings : list of strings
        A list of the time to execute certain tasks during the calling of 
        CopyRoi() and subsequent tasks called within CopyRts() or CopySeg().
        
        
    Notes:
    *****
    
    An example use of the ForceRegistration option is to overcome translations 
    in the patient that are not reflected in the ImagePositionPatient tag. If 
    using resampling techniques the Target contours/segmentations will appear 
    displaced w.r.t. the anatomical features highlighted in the Source image.
    """

    import time
    import os
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    from GeneralTools import ExportDictionaryToJson, ExportListToTxt
    
    RunDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    ListOfInputs = [f'RunDateTime: {RunDateTime}\n',
                    f'SrcRoiFpath = {SrcRoiFpath}\n',
                    f'FromSliceNum = {FromSliceNum}\n',
                    f'FromRoiLabel = {FromRoiLabel}\n',
                    f'SrcDcmDir = {SrcDcmDir}\n',
                    f'TrgDcmDir = {TrgDcmDir}\n',
                    f'TrgRoiFpath = {TrgRoiFpath}\n',
                    f'ToSliceNum = {ToSliceNum}\n',
                    f'Interpolation = {Interpolation}\n',
                    f'Variance = {Variance}\n',
                    f'ThreshPreRes = {ThreshPreRes}\n',
                    f'ThreshPostTx = {ThreshPostTx}\n',
                    f'Tx = {Tx}\n',
                    f'ForceRegistration = {ForceRegistration}\n',
                    f'TxtToAddToRoiLabel = {TxtToAddToRoiLabel}\n',
                    f'LogToConsole = {LogToConsole}']
    
    DictOfInputs = {'RunDateTime' : RunDateTime,
                    'SrcRoiFpath' : SrcRoiFpath,
                    'FromSliceNum' : FromSliceNum,
                    'FromRoiLabel' : FromRoiLabel,
                    'SrcDcmDir' : SrcDcmDir,
                    'TrgDcmDir' : TrgDcmDir,
                    'TrgRoiFpath' : TrgRoiFpath,
                    'ToSliceNum' : ToSliceNum,
                    'Interpolation' : Interpolation,
                    'Variance' : Variance,
                    'ThreshPreRes' : ThreshPreRes,
                    'ThreshPostTx' : ThreshPostTx,
                    'Tx' : Tx,
                    'ForceRegistration' : ForceRegistration,
                    'TxtToAddToRoiLabel' : TxtToAddToRoiLabel,
                    'LogToConsole' : LogToConsole}
    
    
    ListOfTimings = [f'RunDateTime: {RunDateTime}\n']
    
    times = []
    
    
    """ Start timing. """
    times.append(time.time())
    
    """ Establish whether the inputs are valid. """
    CheckValidityOfInputs(SrcRoiFpath, FromRoiLabel, FromSliceNum, ToSliceNum, 
                          TrgRoiFpath, LogToConsole)
    
    times.append(time.time())
    Dtime = round(1000*(times[-1] - times[-2]), 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} ms to check the validity of inputs to CopyRoi().\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        
    """ Determine which Use Case to apply. """
    UseCaseThatApplies,\
    UseCaseToApply = WhichUseCase(FromSliceNum, FromRoiLabel, ToSliceNum, 
                                  SrcDcmDir, TrgDcmDir, ForceRegistration,
                                  LogToConsole=True)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        #print(f'\n\nDone.  Took {Dtime} s to run.')
        msg = f'Took {Dtime} s to determine which UseCase applies.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
    
    DictOfInputs['UseCaseThatApplies'] = UseCaseThatApplies
    DictOfInputs['UseCaseToApply'] = UseCaseToApply
    
    #print(f'\n\n\nSrcRoiFpath = {SrcRoiFpath}')
    
    Modality = dcmread(SrcRoiFpath).Modality
    
    if Modality == 'RTSTRUCT':        
        TrgRoi, ListOfTimings\
        = CopyRts(SrcRoiFpath, FromSliceNum, FromRoiLabel, SrcDcmDir, 
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, ToSliceNum, 
                  Interpolation, Variance, ThreshPreRes, ThreshPostTx, Tx,
                  ForceRegistration, TxtToAddToRoiLabel, LogToConsole,
                  ListOfTimings)
    
    elif Modality == 'SEG':
        TrgRoi, ListOfTimings\
        = CopySeg(SrcRoiFpath, FromSliceNum, FromRoiLabel, SrcDcmDir,
                  TrgDcmDir, UseCaseToApply, TrgRoiFpath, ToSliceNum, 
                  Interpolation, Variance, ThreshPreRes, ThreshPostTx, Tx,
                  ForceRegistration, TxtToAddToRoiLabel, LogToConsole, 
                  ListOfTimings)
        
    else:
        msg = f'The Source modality ({Modality}) must be "RTS" or "SEG".'
        
        raise Exception(msg)

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        msg = f'Took {Dtime} s to copy the ROI(s) and create the {Modality}'\
              + ' object.\n'  
        ListOfTimings.append(msg)
        print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a sub-
        directory "logs" in the current working directory. """
        JsonFname = RunDateTime + '_DictOfInputs.json'
        TxtFname = RunDateTime + '_ListOfInputs.txt'
        
        CWD = os.getcwd()
        
        ExportDir = os.path.join(CWD, 'logs')
        
        JsonFpath = os.path.join(ExportDir, JsonFname)
        TxtFpath = os.path.join(ExportDir, TxtFname)
        
        ExportDictionaryToJson(DictOfInputs, JsonFname, ExportDir)
        ExportListToTxt(ListOfInputs, TxtFname, ExportDir)
            
        print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
              '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
    return TrgRoi, DictOfInputs, ListOfInputs, ListOfTimings







"""
******************************************************************************
******************************************************************************
CHECK FOR ERRORS IN NEW TARGET RTS / SEG
******************************************************************************
******************************************************************************
"""

def ErrorCheckRoi(Roi, DicomDir, LogToConsole=False, DictOfInputs=None,
                  ExportLogFiles=False):
    """
    Check a RTS/SEG for errors in dependencies based on provided directory of
    the DICOMs that relate to the RTS/SEG.  
    
    Inputs:
    ******
    
    Roi : Pydicom object
        RTS/SEG ROI object to be error checked.
        
    DicomDir : string
        Directory containing the DICOMs that relate to Roi.
        
    LogToConsole : boolean (default False)
        Denotes whether some results will be logged to the console.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If !=None, it is a dictionary containing the inputs that were called to 
        CopyRoi().
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
                           
    
    Outputs:
    *******
    
    LogList : list of strings
        A list of strings that describe info or any errors that are found.
    
    Nerrors : integer
        The number of errors found.
    """
    
    import time
    import os
    from GeneralTools import ExportListToTxt
    
    Modality = Roi.Modality
    
    """ Start timing. """
    times = []
    times.append(time.time())
    
    if Modality == 'RTSTRUCT':
        #import RtsTools
        #import importlib
        #importlib.reload(RtsTools)
        from RtsTools import ErrorCheckRts
        
        LogList, Nerrors = ErrorCheckRts(Roi, DicomDir, LogToConsole)
        
    elif Modality == 'SEG':
        #import importlib
        #import SegTools
        #importlib.reload(SegTools)
        from SegTools import ErrorCheckSeg
        
        LogList, Nerrors = ErrorCheckSeg(Roi, DicomDir, LogToConsole)
        
    else:
        msg = f'The modality ({Modality}) must be either "RTS" or "SEG".'
        
        raise Exception(msg)
        
        LogList = None
        Nerrors = None
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to error check the {Modality} object.\n')
    
    
    if ExportLogFiles:
        """ Export log of error check results. """
        if not DictOfInputs == None:
            DateTime = DictOfInputs['RunDateTime']
        else:
            DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        Fname = DateTime + '_ErrorCheckLog.txt'
        
        CWD = os.getcwd()
        ExportDir = os.path.join(CWD, 'logs')
        
        Fpath = os.path.join(ExportDir, Fname)
        
        ExportListToTxt(LogList, Fname, ExportDir)
            
        print('Log of error checks saved to:\n', Fpath, '\n')
    
    return LogList, Nerrors






"""
******************************************************************************
******************************************************************************
EXPORT NEW TARGET RTS / SEG TO DISK
******************************************************************************
******************************************************************************
"""

def ExportTrgRoi(TrgRoi, SrcRoiFpath, ExportDir, TxtToAddToFname='',
                 DictOfInputs=None):
    """
    Export RTS/SEG to disk.  
    
    Inputs:
    ******
    
    TrgRoi : Pydicom object
        Target RTS/SEG ROI object to be exported.
        
    SrcRoiFpath : string
        Full path of the Source RTS/SEG file (used to generate the filename of
        the new RTS/SEG file).
                              
    ExportDir : string
        Directory where the new RTS/SEG is to be exported.
    
    NamePrefix : string (optional; '' by default)
        Prefix to be added to the assigned filename (after the DateTime stamp), 
        e.g. 'Case3b-i'.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If not None, a dictionary containing the inputs that were called to 
        CopyRoi().
                           
    
    Outputs:
    *******
    
    TrgRoiFpath : string
        Full path of the exported Target RTS/SEG file.
    """
    
    import os
    import time
    
    # Get the filename of the original RTS/SEG file:
    #SrcRoiFname = os.path.split(SrcRoiFpath)[1]
    
    if not os.path.isdir(ExportDir):
        os.mkdir(ExportDir)
    
    if not DictOfInputs == None:
        DateTime = DictOfInputs['RunDateTime']
    else:
        DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    #FnamePrefix = DateTime + '_' + NamePrefix + '_from_' 
    FnamePrefix = DateTime + '_' + TxtToAddToFname.replace(' ', '_')
    
    # Create a new filename (this will appear under Label in XNAT):
    #TrgRoiFname = FnamePrefix + SrcRoiFname
    #TrgRoiFname = FnamePrefix
    TrgRoiFname = FnamePrefix + '.dcm'
    
    TrgRoiFpath = os.path.join(ExportDir, TrgRoiFname)
    
    TrgRoi.save_as(TrgRoiFpath)
        
    print(f'New Target {TrgRoi.Modality} exported to:\n', TrgRoiFpath, '\n')
    
    return TrgRoiFpath







"""
******************************************************************************
******************************************************************************
RUN COPYROI()
******************************************************************************
******************************************************************************
"""


def RunCopyRoi(TestNums, LogToConsole=False, TxtToAddToRoiLabel='', 
               Interpolation='BlurThenLinear', Variance=(1,1,1), 
               ThreshPreRes=0.75, ThreshPostTx=0.05, Tx='affine', 
               ForceRegistration=False, ExportRoi=True, PlotResults=False, 
               PlotAllSlices=False, ExportPlot=False, ExportLogFiles=True):
    """
    
    Inputs:
    ******
    
    TestNums : list of strings
        A list of strings denoting the tests to run (as defined in 
        CreateDictOfPaths()). Several tests can be run in series by providing
        a list of comma-separated strings (e.g. ['SR1', 'SR2', 'SR3']). To 
        run a single test a single-itemed list must be provided (e.g. ['RD4']).
    
    TxtToAddToRoiLabel : string (optional, '' by default)
        String of text to pass to CreateRts()/CreateSeg(). The string will be
        appended to the RTS StructureSetLabel or SEG SeriesDescription, and to 
        the filename of the exported (new) Target RTS/SEG.
    
    LogToConsole : boolean (optional; False by default)
        If True, intermediate results will be logged to the console during the
        running of CopyRoi().
    
    Interpolation : string (optional; 'BlurThenLinear' by default)
        The type of interpolation to be used if the Source labelmap image(s)
        is/are to be resampled to the Target grid.  Acceptable inputs are:
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'LabelGaussian' (or 'Gaussian' or 'gaussian')
        - 'BlurThenLinear' (or 'Linear' or 'linear') after Gaussian blurring 
        (followed by binary thresholding) (Default value)
    
    Variance : tuple of floats (optional; (1, 1, 1) by default)
        A tuple (for each dimension) of the variance to be applied if a 
        Gaussian blurring + linearly resampling approach is taken for 
        resampling of the Source labelmap image(s).
        
    ThreshPreRes : float (optional; 0.75 by default)
        The threshold level used to perform binary thresholding prior to 
        resampling. Example use is on a labelmap image if a Gaussian blurring 
        + linearly resampling approach is taken to combat aliasing effects.
    
    ThreshPostTx : float (optional; 0.05 by default)
        The threshold level used to perform binary thresholding after 
        registration transformation.
    
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    ForceRegistration : boolean (optional; False by default)
        If True the Source image will be registered to the Target image, and 
        the Source labelmap will be transformed to the Target image grid
        accordingly.  
    
    ExportRoi : boolean (optional; True by default)
        If True, the new Target RTS/SEG will be exported to a directory called
        "new_RTS"/"new_SEG" in the current working directory.
        
    PlotResults : boolean (optional; False by default)
        If True, the new RTS/SEG will be parsed and all contours/segmentations
        will be plotted.
    
    PlotAllSlices : boolean (optional; False by default)
        If True, and if PlotResults is True, all slices in Source and Target
        will be plotted. If False, only slices containing contours/segmentations
        will be plotted.
    
    ExportPlot : boolean (optional; False by default)
        If True, and if PlotResults is True, the plot will be exported to a .jpg
        in a directory called "plots_RTS"/"plots_SEG" in the current working 
        directory.
    
    ExportLogFiles : boolean (default False)
        If True, log files will be exported.
    """

    import time
    import os
    from pydicom import dcmread
    #from DicomTools import IsSameModalities
    from GeneralTools import PrintTitle, ExportListToTxt, ExportDictionaryToJson
    
    if not isinstance(TestNums, list):
        msg = 'The input argument "TestNums" must be a list of character '\
              + 'strings.'
        print(msg)
    
    CWD = os.getcwd()
    RtsExportDir = os.path.join(CWD, 'new_RTS')
    SegExportDir = os.path.join(CWD, 'new_SEG')
    RtsPlotExportDir = os.path.join(CWD, 'plots_RTS')
    SegPlotExportDir = os.path.join(CWD, 'plots_SEG')
    LogExportDir = os.path.join(CWD, 'logs')
    
    ListOfTimings = [] # list of strings (messages outputted to console)
    
    """ Start timing. """
    Times = []
    Times.append(time.time())
    
    
    for TestNum in TestNums:
        TestRunTimes = [] # this this Test run only
        TestRunTimes.append(time.time())
        
        PrintTitle(f'Running Test {TestNum}:')
    
        SrcRoiFpath, FromSliceNum, FromRoiLabel, ToSliceNum, SrcDcmDir,\
        TrgDcmDir, TrgRoiFpath, TxtToAddToRoiLabel, SrcLabel,\
        TrgLabel = GetPathInputs(TestNum)
        
        Times.append(time.time())
        Dtime = round(1000*(Times[-1] - Times[-2]), 1)
        msg = f'Took {Dtime} ms to get the inputs to CopyRoi().\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if LogToConsole:
            print(f'Inputs to CopyRoi():\n FromRoiLabel = {FromRoiLabel}\n',
                      f'FromSliceNum = {FromSliceNum}\n',
                      f'ToSliceNum = {ToSliceNum}\n')
        
    
        """ Import the RTS or SEGs: """
        SrcRoi = dcmread(SrcRoiFpath)
        SrcModality = SrcRoi.Modality
        
        #print(f'\n\nSrcRoiFpath = {SrcRoiFpath}')
        #print(f'SrcModality = {SrcModality}')
    
        if TrgRoiFpath:
            TrgRoi = dcmread(TrgRoiFpath)
        else:
            TrgRoi = None
        
        """ New ROI label for RTS StructureSetLabel / SEG Series Description 
        and RTS/SEG filename: """
        if 'RR' in TestNum or 'RD' in TestNum:
            NewRoiLabel = SrcRoi.StructureSetLabel + TxtToAddToRoiLabel
        else:
            """ 'SR' in TestNum or 'SD' in TestNum: """
            NewRoiLabel = SrcRoi.SeriesDescription + TxtToAddToRoiLabel
            
        
        """ Text to add to file names: """
        TxtToAddToFname = f'Test_{TestNum}_{NewRoiLabel}'
        
        
        """ Copy Contours/Segmentations """
        
        if SrcModality == 'RTSTRUCT':
            NewTrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
            = CopyRoi(SrcRoiFpath, FromSliceNum, FromRoiLabel, SrcDcmDir, 
                      TrgDcmDir, TrgRoiFpath, ToSliceNum, Interpolation, 
                      Variance, ThreshPreRes, ThreshPostTx, Tx,
                      ForceRegistration, TxtToAddToRoiLabel, LogToConsole, 
                      False)
        else:
            NewTrgRoi, DictOfInputs, ListOfInputs, ListOfTimings\
            = CopyRoi(SrcRoiFpath, FromSliceNum, FromRoiLabel, SrcDcmDir,  
                      TrgDcmDir, TrgRoiFpath, ToSliceNum, Interpolation,  
                      Variance, ThreshPreRes, ThreshPostTx, Tx, 
                      ForceRegistration, TxtToAddToRoiLabel, LogToConsole, 
                      False)
        
        
        if ExportLogFiles:
            """ Export DictOfInputs as a json, and ListOfInputs as a txt, to a 
            sub-directory "logs" in the current working directory. """
            
            if not DictOfInputs == None:
                DateTime = DictOfInputs['RunDateTime']
            else:
                DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
                
            JsonFname = DateTime + f'_TestRun_{TestNum}_DictOfInputs.json'
            TxtFname = DateTime + f'_TestRun_{TestNum}_ListOfInputs.txt'
            
            CWD = os.getcwd()
            
            ExportDir = os.path.join(CWD, 'logs')
            
            JsonFpath = os.path.join(ExportDir, JsonFname)
            TxtFpath = os.path.join(ExportDir, TxtFname)
            
            ExportDictionaryToJson(DictOfInputs, JsonFname, ExportDir)
            ExportListToTxt(ListOfInputs, TxtFname, ExportDir)
                
            print('Dictionary of inputs to CopyRoi() saved to:\n', JsonFpath, 
                  '\nand list of inputs saved to:\n', TxtFpath, '\n')
        
        
        
        """ Error check RTS/SEG """
        Times.append(time.time())
        
        ErrorList, Nerrors = ErrorCheckRoi(NewTrgRoi, TrgDcmDir, LogToConsole,
                                           DictOfInputs, False)
        
        Times.append(time.time())
        Dtime = round(Times[-1] - Times[-2], 1)
        msg = f'Took {Dtime} s to error check the {SrcModality}.  There were '\
              + f'{Nerrors} errors found.\n'
        ListOfTimings.append(msg)
        print(f'*{msg}')
        
        if ExportLogFiles:
            """ Export log of error check results. """
            
            Fname = DateTime + f'_TestRun_{TestNum}_ErrorCheckLog.txt'
            
            CWD = os.getcwd()
            ExportDir = os.path.join(CWD, 'logs')
            
            Fpath = os.path.join(ExportDir, Fname)
            
            ExportListToTxt(ErrorList, Fname, ExportDir)
                
            print('Log of error checks saved to:\n', Fpath, '\n')
        
        
    
        """ Export RTS/SEG """
    
        if ExportRoi:# and not Nerrors:
            UseCaseThatApplies = DictOfInputs['UseCaseThatApplies']
            UseCaseToApply = DictOfInputs['UseCaseToApply']
            
            if ForceRegistration and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
                TxtToAddToFname += f'_ForcedRegistration_{Tx}'
                
            if not ForceRegistration and UseCaseThatApplies in ['3a', '3b', '4a', '4b']:
                TxtToAddToFname += f'_{Interpolation}'
                
            if UseCaseToApply in ['5a', '5b']:
                TxtToAddToFname += f'_{Tx}'
            
            if NewTrgRoi.Modality == 'RTSTRUCT':
                NewTrgRoiFpath = ExportTrgRoi(NewTrgRoi, SrcRoiFpath, RtsExportDir, 
                                              TxtToAddToFname, DictOfInputs)
            else:
                NewTrgRoiFpath = ExportTrgRoi(NewTrgRoi, SrcRoiFpath, SegExportDir, 
                                              TxtToAddToFname, DictOfInputs)
        
        
        """ Plot Contours """
    
        if PlotResults:
            if ExportPlot:
                dpi = 120
            else:
                dpi = 80
    
    
    
            if TrgRoi:
                ListOfRois = [SrcRoi, TrgRoi, NewTrgRoi]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'Original Target', 'New Target']
            else:
                ListOfRois = [SrcRoi, NewTrgRoi]
    
                ListOfDicomDirs = [SrcDcmDir, TrgDcmDir]
    
                ListOfPlotTitles = ['Source', 'New Target']
    
    
            #TxtToAdd = f'(RunCase = {RunCase}, \nFromSliceNum = {FromSliceNum}, '\
            #           + f'\nToSliceNum = {ToSliceNum})'
            
            #print(f'\nListOfDicomDirs = {ListOfDicomDirs}')
            
            Times.append(time.time())
            
            if NewTrgRoi.Modality == 'RTSTRUCT':
                from PlottingTools import PlotContoursFromListOfRtss_v1
                
                PlotContoursFromListOfRtss_v1(ListOfRois, ListOfDicomDirs, 
                                              ListOfPlotTitles,
                                              PlotAllSlices=PlotAllSlices,
                                              AddTxt=TxtToAddToFname, 
                                              ExportPlot=ExportPlot, 
                                              ExportDir=RtsPlotExportDir, 
                                              dpi=dpi, LogToConsole=False)
    
            else:
                from PlottingTools import PlotPixArrsFromListOfSegs_v1
                
                PlotPixArrsFromListOfSegs_v1(ListOfRois, ListOfDicomDirs, 
                                             ListOfPlotTitles,
                                             PlotAllSlices=PlotAllSlices, 
                                             AddTxt=TxtToAddToFname, 
                                             ExportPlot=ExportPlot, 
                                             ExportDir=SegPlotExportDir,  
                                             dpi=dpi, LogToConsole=False)
                
            Times.append(time.time())
            Dtime = round(Times[-1] - Times[-2], 1)
            msg = f'Took {Dtime} s to plot the results.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
        
            TestRunTimes.append(time.time())
            Dtime = round(TestRunTimes[-1] - TestRunTimes[-2], 1)
            msg = f'Took {Dtime} s to run test {TestNum}.\n'
            ListOfTimings.append(msg)
            print(f'*{msg}')
    
        
    print('All Test run iterations complete.\n')
    
    Times.append(time.time())
    Dtime = round(Times[-1] - Times[0], 1)
    msg = f'Took {Dtime} s to run all tests.\n'
    ListOfTimings.append(msg)
    print(f'*{msg}')
    
    
    if ExportLogFiles:
        """ Export the ListOfTimings """
        
        Fname = DateTime + f'_TestRun_{TestNum}_ListOfTimings.txt'
        
        Fpath = os.path.join(LogExportDir, Fname)
        
        ExportListToTxt(ListOfTimings, Fname, LogExportDir)
        
        print('Log file saved to:\n', Fpath, '\n')
    
    
    return 
    #return NewTrgRoi
    








"""
******************************************************************************
******************************************************************************
EXPORT A DICOM SPATIAL OR DEFORMABLE SPATIAL REGISTRATION OBJECT (SRO) TO DISK
******************************************************************************
******************************************************************************
"""


def ExportDro_OLD(SrcRoi, TrgRoi, ExportDir, Description=''):
    """
    Export a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO) to disk.  
    
    Inputs:
    ******
    
    SrcRoi : Pydicom object
        Source RTS/SEG ROI object. Any tags associated with the "moving" image
        will be assigned to the corresponding tag values in SrcRoi. 
        
    TrgRoi : Pydicom object
        The Target RTS/SEG ROI object that was created from the copying 
        process.  Any tags associated with the "fixed" image will be assigned
        to the corresponding tag values in TrgRoi.
            
    ExportDir : string
        Directory where the new DRO is to be exported.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
                           
    
    Outputs:
    *******
    
    DroFpath : string
        Full path of the exported DRO.
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #from copy import deepcopy
    from GeneralTools import ParseTransformParameterFile
    from DicomTools import ImportDicoms
    
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0']
    
    
    
    # Sample DICOM Registration objects:
    DroRootDir = r'C:\Users\ctorti\Documents\DICOM Registration Object\IJ_923_dicom-sro.tar\IJ_923_dicom-sro\dicom-sro'
    
    DroDir1 = r"rect-offset-sro-rigid" # Spatial Registration object
    DroDir2 = r"sphere-centered-sro-deformable" # Deformable Spatial 
    # Registration object
    
    DroFname1 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm"
    DroFname2 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm"
    
    DroFpath1 = os.path.join(DroRootDir, DroDir1, DroFname1)
    DroFpath2 = os.path.join(DroRootDir, DroDir2, DroFname2)
    
    #Dro1 = dcmread(DroFpath1) # Spatial Registration Storage
    #Dro2 = dcmread(DroFpath2) # Deformable Spatial Registration Storage
    
    """
    If the registration was non-elastic use Dro1 as a template.
    If the registration was elastic use Dro2.
    """
    
    if True:
        Dro = dcmread(DroFpath1)
    else:
        Dro = dcmread(DroFpath2)
    
    
    CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    CurrentTime = time.strftime("%H%M%S", time.gmtime())
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Dro.InstanceCreationData = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    if TrgRoi.Modality == 'RTSTRUCT':
        ContentDate = TrgRoi.StructureSetDate
        ContentTime = TrgRoi.StructureSetTime
        SeriesDescription = TrgRoi.StructureSetLabel
        
        TrgRefSOPClassUID = TrgRoi.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPClassUID
        
        SrcRefSOPClassUID = SrcRoi.ReferencedFrameOfReferenceSequence[0]\
                                  .RTReferencedStudySequence[0]\
                                  .RTReferencedSeriesSequence[0]\
                                  .ContourImageSequence[0]\
                                  .ReferencedSOPClassUID
        
        
                                          
    else:
        ContentDate = TrgRoi.ContentDate
        ContentTime = TrgRoi.ContentTime
        SeriesDescription = TrgRoi.SeriesDescription
        
        TrgRefSOPClassUID = TrgRoi.ReferencedSeriesSequence[0]\
                                  .ReferencedInstanceSequence[0]\
                                  .ReferencedSOPClassUID
        
        SrcRefSOPClassUID = SrcRoi.ReferencedSeriesSequence[0]\
                                  .ReferencedInstanceSequence[0]\
                                  .ReferencedSOPClassUID
        
        
    
    
    Dro.StudyDate = TrgRoi.StudyDate
    Dro.SeriesDate = TrgRoi.SeriesDate
    Dro.ContentDate = ContentDate
    Dro.StudyTime = TrgRoi.StudyTime
    Dro.SeriesTime = TrgRoi.SeriesTime
    Dro.ContentTime = ContentTime
    Dro.Manufacturer = TrgRoi.Manufacturer
    """Consider modifying StudyDescription (= '' in the template DRO."""
    Dro.SeriesDescription = SeriesDescription
    Dro.ManufacturerModelName = TrgRoi.ManufacturerModelName
    Dro.PatientName = TrgRoi.PatientName
    Dro.PatientID = TrgRoi.PatientID
    Dro.PatientBirthDate = TrgRoi.PatientBirthDate
    Dro.PatientSex = TrgRoi.PatientSex
    Dro.PatientAge = TrgRoi.PatientAge
    Dro.StudyInstanceUID = TrgRoi.StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    Dro.SeriesInstanceUID = generate_uid()
    Dro.StudyID = TrgRoi.StudyID
    Dro.SeriesNumber = TrgRoi.SeriesNumber
    Dro.InstanceNumber = TrgRoi.InstanceNumber
    Dro.FrameOfReferenceUID = TrgRoi.FrameOfReferenceUID
    Dro.PositionReferenceIndicator = TrgRoi.PositionReferenceIndicator
    """Keep ContentLabel as 'REGISTRATION'."""
    Dro.ContentDescription = Description
    if TrgRoi.Modality == 'SEG':
        Dro.ContentCreatorName = TrgRoi.ContentCreatorName
    
    
    
    """ Modify the RegistrationSequence for the fixed domain. """
    
    Dro.RegistrationSequence[0]\
       .FrameOfReferenceUID = TrgRoi.FrameOfReferenceUID
    
    
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    
    Need to cover 'RIGID_SCALE' and 'AFFINE' cases...
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = 
    """
    
    """
    Modify the FrameOfReferenceTransformationMatrix.  No change required to
     that of the template DRO is the identity transformation is desired, i.e.
    ['1.0', '0.0', '0.0', '0.0', 
     '0.0', '1.0', '0.0', '0.0', 
     '0.0', '0.0', '1.0', '0.0', 
     '0.0', '0.0', '0.0', '1.0']
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = 
    """
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the RegistrationSequence for the moving domain. """
    
    Dro.RegistrationSequence[1]\
       .FrameOfReferenceUID = SrcRoi.FrameOfReferenceUID
    
    
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    
    Need to cover 'RIGID_SCALE' and 'AFFINE' cases...
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = 
    """
    
    """ Modify the FrameOfReferenceTransformationMatrix. """
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxParams
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ continue with Series and Instance Reference Macro attributes, as 
    covered in 10.3 (page 15) of the DICOM Supplement 73. """
    
    
    """ Modify the ReferencedSeriesSequence for the fixed domain. """
    
    Dro.ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = TrgRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the first 
    ReferencedSeriesSequence referenced the StudyInstanceUID of the fixed 
    series."""
    Dro.ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = TrgRoi.StudyInstanceUID
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgRoi.SeriesInstanceUID
    
    
    
    
    
    """ Modify the ReferencedSeriesSequence for the moving domain. """
    
    Dro.ReferencedSeriesSequence[1]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = SrcRefSOPClassUID
    
    Dro.ReferencedSeriesSequence[1]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = SrcRoi.StudyInstanceUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the second 
    ReferencedSeriesSequence referenced the SOPInstanceUID of 31st DICOM (in a
    collection of 100 DICOMs) of the moving series.
    Dro.ReferencedSeriesSequence[1]\
       .SeriesInstanceUID = SrcRoi.?????????
    """
    
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    fixed domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = TrgRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the first 
    StudiesContainingOtherReferencedInstancesSequence referenced the 
    SOPInstanceUID of the last DICOM in the fixed series.
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = ?????
    """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgRoi.SeriesInstanceUID
    
    
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    moving domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPClassUID = SrcRefSOPClassUID
    
    """In the template DRO the ReferencedSOPInstanceUID in the second 
    StudiesContainingOtherReferencedInstancesSequence referenced the 
    SOPInstanceUID of the last DICOM in the moving series.
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .ReferencedInstanceSequence[0]\
       .ReferencedSOPInstanceUID = ?????
    """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcRoi.SeriesInstanceUID
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    """
    Delete all GroupLength tags?
    """
    
    
    DroFname = CurrentDateTime + '_' + Description.replace(' ', '_') + '.dcm'
    
    DroFpath = os.path.join(ExportDir, DroFname)
    
    Dro.save_as(DroFpath)
        
    print('\nNew DICOM Registration Object exported to:\n', DroFpath)
    
    return DroFpath





def CreateDro(SrcDicomDir, TrgDicomDir, Description='', LogToConsole=False):
    """
    Create a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO).  
    
    *COMMENT*:  At the moment only the Spatial Registration Object Storage IOD
    is covered.
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    from copy import deepcopy
    from DicomTools import ImportDicoms
    from GeneralTools import ParseTransformParameterFile
    from GeneralTools import GetTxMatrixType
    
    # Start timing:
    times = []
    times.append(time.time())
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
        
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to parse the transform parameter file.')

    TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0']
    
    
    
    # Sample DICOM Registration objects:
    DroRootDir = r'C:\Users\ctorti\Documents\DICOM Registration Object\IJ_923_dicom-sro.tar\IJ_923_dicom-sro\dicom-sro'
    
    DroDir1 = r"rect-offset-sro-rigid" # Spatial Registration object
    DroDir2 = r"sphere-centered-sro-deformable" # Deformable Spatial 
    # Registration object
    
    DroFname1 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm"
    DroFname2 = "2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm"
    
    DroFpath1 = os.path.join(DroRootDir, DroDir1, DroFname1)
    DroFpath2 = os.path.join(DroRootDir, DroDir2, DroFname2)
    
    #Dro1 = dcmread(DroFpath1) # Spatial Registration Storage
    #Dro2 = dcmread(DroFpath2) # Deformable Spatial Registration Storage
    
    """
    If the registration was non-elastic use Dro1 as a template.
    If the registration was elastic use Dro2.
    """
    
    if True:
        Dro = dcmread(DroFpath1)
    else:
        Dro = dcmread(DroFpath2)
    
    
    CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    CurrentTime = time.strftime("%H%M%S", time.gmtime())
    #CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Dro.InstanceCreationData = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.StudyDate = TrgDicoms[0].StudyDate
    Dro.SeriesDate = TrgDicoms[0].SeriesDate
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.StudyTime = TrgDicoms[0].StudyTime
    Dro.SeriesTime = TrgDicoms[0].SeriesTime
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.Manufacturer = TrgDicoms[0].Manufacturer
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    Dro.SeriesDescription = TrgDicoms[0].SeriesDescription
    Dro.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
    Dro.PatientName = TrgDicoms[0].PatientName
    Dro.PatientID = TrgDicoms[0].PatientID
    Dro.PatientBirthDate = TrgDicoms[0].PatientBirthDate
    Dro.PatientSex = TrgDicoms[0].PatientSex
    Dro.PatientAge = TrgDicoms[0].PatientAge
    Dro.StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    Dro.SeriesInstanceUID = generate_uid()
    Dro.StudyID = TrgDicoms[0].StudyID
    Dro.SeriesNumber = TrgDicoms[0].SeriesNumber
    Dro.InstanceNumber = TrgDicoms[0].InstanceNumber
    Dro.FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    Dro.PositionReferenceIndicator = TrgDicoms[0].PositionReferenceIndicator
    """ Keep ContentLabel as 'REGISTRATION'. """
    """ Keep ContentCreatorName as ''. """
    
    
    
    """ Modify the RegistrationSequence for the fixed domain. """
    
    Dro.RegistrationSequence[0]\
       .FrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    
    """
    Since the first sequence in RegistrationSequence will be the fixed image
    domain, FrameOfReferenceTransformationMatrix will be the identity matrix
    ['1.0', '0.0', '0.0', '0.0', 
     '0.0', '1.0', '0.0', '0.0', 
     '0.0', '0.0', '1.0', '0.0', 
     '0.0', '0.0', '0.0', '1.0'].
    
    Hence no change is required FrameOfReferenceTransformationMatrix for this
    sequence.
    """
    
    
    """
    Since the FrameOfReferenceTransformationMatrix is the identity matrix, the 
    FrameOfReferenceTransformationMatrixType will be 'RIGID'.  Hence no change
    is required to FrameOfReferenceTransformationMatrixType for this sequence.
    """
    
    
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the RegistrationSequence for the moving domain. """
    
    Dro.RegistrationSequence[1]\
       .FrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    
    """ Modify the FrameOfReferenceTransformationMatrix. """
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxParams
       
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' = value in template DRO
    'RIGID_SCALE'
    'AFFINE'    
    """
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
        
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the ReferencedSeriesSequence for the fixed domain. """
    
    print('\nModify the ReferencedSeriesSequence for the fixed domain.')
    
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(TrgDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in ReferencedSeriesSequence for',
              'the fixed series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = TrgDicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = TrgDicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the fixed series.')
    
    
    
    
    """ Modify the ReferencedSeriesSequence for the moving domain. """
    
    print('\nModify the ReferencedSeriesSequence for the moving domain.')
    
    Dro.ReferencedSeriesSequence[1]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[1]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(SrcDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[1]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[1]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[1]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in ReferencedSeriesSequence for',
              'the moving series.')
        
    # Modify the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[1]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = SrcDicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[1]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = SrcDicoms[i].SOPInstanceUID

    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the moving series.')
        
    
    print('\nModify the StudiesContainingOtherReferencedInstancesSequence',
          'for the fixed domain.')
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    fixed domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
       
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(TrgDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'fixed series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = TrgDicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = TrgDicoms[i].SOPInstanceUID
    
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to modify the',
              'ReferencedInstanceSequences in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'fixed series.')
        
        
    print('\nModify the StudiesContainingOtherReferencedInstancesSequence for',
          'the moving domain.')
    
    """ Modify the StudiesContainingOtherReferencedInstancesSequence for the 
    moving domain. """
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .StudyInstanceUID = SrcDicoms[0].StudyInstanceUID
       
    Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
    
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    ReqRIS = len(SrcDicoms)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              'ReferencedInstanceSequence in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'moving series.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = SrcDicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[1]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = SrcDicoms[i].SOPInstanceUID
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to modify the',
              'ReferencedInstanceSequences in',
              'StudiesContainingOtherReferencedInstancesSequence for the',
              'moving series.')
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    print('Deleting GroupLength tags.')
    
    """
    The SRO template contains 4 retired tags with name "GroupLength". 
    Delete them.
    """
    
    del Dro[0x08, 0x00]
    del Dro[0x10, 0x00]
    del Dro[0x20, 0x00]
    del Dro[0x70, 0x00]
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if True:#LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
        
        
    return Dro





def ExportDro(Dro, ExportDir, Description='', LogToConsole=False):
    """
    Export a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO) to disk.  
    
    *COMMENT*:  At the moment only the Spatial Registration Object Storage IOD
    is covered.
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Registration Object to export.
        
    ExportDir : string
        Directory where the new DRO is to be exported.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    DroFpath : string
        Full path of the exported DRO.
    """
    
    import os
    import time
    
    CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    DroFname = CurrentDateTime + '_' + Description.replace(' ', '_') + '.dcm'
    
    DroFpath = os.path.join(ExportDir, DroFname)
    
    print('Exporting DRO.')
    
    Dro.save_as(DroFpath)
        
    print('\nNew DICOM Registration Object exported to:\n', DroFpath)
    
    return DroFpath





