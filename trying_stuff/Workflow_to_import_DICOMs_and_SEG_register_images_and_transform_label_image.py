# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 16:02:45 2021

@author: ctorti
"""

""" Code posted to SimpleITK GitHub issue page:
    https://github.com/SimpleITK/SimpleITK/issues/1358
"""



import SimpleITK as sitk
import numpy as np
from pydicom import dcmread
import os
import time
import matplotlib.pyplot as plt
from ipywidgets import interact, fixed
%matplotlib inline # for Anaconda Jupyter

def ImportDicomSeriesAsSitkImage(DicomDir, PixType='sitkFloat32'):
    """ Import a DICOM series as a SimpleITK image. """
    
    # Only a limited number of pixel types are catered for:
    PixTypesDict = {'sitkUInt8': sitk.sitkUInt8,
                    'sitkUInt16' : sitk.sitkUInt16,
                    'sitkUInt32' : sitk.sitkUInt32,
                    'sitkFloat32' : sitk.sitkFloat32,
                    'sitkFloat64' : sitk.sitkFloat64,
                    'sitkVectorUInt8' : sitk.sitkVectorUInt8,
                    'sitkVectorUInt16' : sitk.sitkVectorUInt16,
                    'sitkVectorFloat64' : sitk.sitkVectorFloat64}
    
    PixTypes = list(PixTypesDict.keys())
    
    if PixType in PixTypes:
        PixType = PixTypesDict[PixType]
    else:
        raise Exception("'PixType' '{PixType}' is not allowed.")
    
    
    Reader = sitk.ImageSeriesReader()
    Fnames = Reader.GetGDCMSeriesFileNames(DicomDir)
    Reader.SetFileNames(Fnames)
    
    Image = Reader.Execute()
    Image = sitk.Cast(Image, PixType)
    
    return Image


def GetDicomSOPuids(DicomDir):
    """ Get a list of SOP Instance UIDs for a DICOM series using SimpleITK's 
    ImageSeriesReader. """
    
    SeriesReader = sitk.ImageSeriesReader()
    
    FilePaths = SeriesReader.GetGDCMSeriesFileNames(DicomDir)
    
    FileReader = sitk.ImageFileReader()
    
    SOPuids = [] 
    
    for fpath in FilePaths:
        FileReader.SetFileName(fpath)

        FileReader.ReadImageInformation()
    
        SOPuids.append(FileReader.GetMetaData('0008|0018'))
        
    return SOPuids


def GetRSOPuidsInPFFGS(Seg):
    """ Get the list of ReferencedSOPInstanceUIDs in the 
    Per-FrameFunctionalGroupsSequence of a SEG. """
    
    RSOPuids = [] 
    
    sequences = Seg.PerFrameFunctionalGroupsSequence
    
    for sequence in sequences:
        uid = sequence.DerivationImageSequence[0]\
                      .SourceImageSequence[0]\
                      .ReferencedSOPInstanceUID
        
        RSOPuids.append(uid)
        
    return RSOPuids


def GetFrame2SliceInds(Seg, SOPuids):
    """ Get a list of the slice numbers that correspond to each  
    Per-FrameFunctionalGroupsSequence in a SEG. """
    
    # Get the list of ReferencedSOPInstanceUIDs:
    RSOPuids = GetRSOPuidsInPFFGS(Seg)
    
    F2Sinds = [] 
    
    for i in range(len(RSOPuids)):
        F2Sinds.append(SOPuids.index(RSOPuids[i]))
        
    return F2Sinds


def ImportLabIm(SegFpath, DcmDir):
    DcmIm = ImportDicomSeriesAsSitkImage(DcmDir)
    
    Seg = dcmread(SegFpath)
    
    PixArr = Seg.pixel_array
    
    SOPuids = GetDicomSOPuids(DcmDir)
    
    F2Sinds = GetFrame2SliceInds(Seg, SOPuids)
    
    # Number of frames, rows and columns:
    F, R, C = PixArr.shape
    
    # Number of slices:
    S = len(SOPuids)
    
    LabArr = np.zeros((S, R, C), dtype='uint')
    
    for i in range(F):
        s = F2Sinds[i]

        LabArr[s] += PixArr[i]
        
    LabIm = sitk.GetImageFromArray(LabArr)
    LabIm.SetOrigin(DcmIm.GetOrigin())
    LabIm.SetSpacing(DcmIm.GetSpacing())
    LabIm.SetDirection(DcmIm.GetDirection())
    
    return LabIm


def RegisterImages(FixIm, MovIm, FinalInterp):
    """ FinalInterp currently only accepts 'Linear' or 'NearestNeighbor'. """
    
    if FinalInterp == 'Linear':
        Interpolator = sitk.sitkLinear
    elif FinalInterp == 'NearestNeighbor':
        Interpolator = sitk.sitkNearestNeighbor
    else:
        raise Exception("'FinalInterp' must either be 'Linear' or 'NearestNeighbor'.")
        
    
    InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm, 
                                                  sitk.AffineTransform(FixIm.GetDimension()), 
                                                  sitk.CenteredTransformInitializerFilter.MOMENTS)
    
    RegMethod = sitk.ImageRegistrationMethod()

    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(0.05)

    RegMethod.SetInterpolator(sitk.sitkLinear)

    """ Optimizer settings. """
    RegMethod.SetOptimizerAsGradientDescent(learningRate=1.0, 
                                            numberOfIterations=500, 
                                            estimateLearningRate=RegMethod.Once)
    RegMethod.SetOptimizerScalesFromPhysicalShift()

    """ Setup for the multi-resolution framework. """     
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    #RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    RegMethod.SetInitialTransform(InitialTx)

    FinalTx = RegMethod.Execute(FixIm, MovIm)

    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, Interpolator, 0.0, MovIm.GetPixelID())
    
    return RegIm, FinalTx


def GetIndOfNearestSlice(Im0, Ind0, Im1, UseCentre=True):
    """ Get the index of the slice in Im1 that is nearest in z-position to the Ind0^th slice in Im0. 
    
    Returns:
    *******
    
    Ind1 : integer
        The slice index in Im1 whose z-position is nearest to Im0[Ind0].
    
    MinDiff : float
        The difference in z-position between Im0[Ind0] and Im1[Ind1] in mm.
    """
    
    # The z-position of Im0 at Ind0:
    z0 = Im0.TransformIndexToPhysicalPoint([0,0,Ind0])[2]
    
    R1, C1, S1 = Im1.GetSize()
    
    if UseCentre:
        i = R1//2
        j = C1//2
    else:
        i = 0
        j = 0
    
    # All z-positions in Im1: """
    Z1 = [Im1.TransformIndexToPhysicalPoint([i,j,k])[2] for k in range(S1)]
    
    Diffs = [abs(z0 - Z1[k]) for k in range(S1)]
    
    MinDiff = min(Diffs)
    
    Ind1 = Diffs.index(MinDiff)
    
    return Ind1, MinDiff


def Image2PixArr(Image):
    """ Convert a SimpleITK image to a Numpy pixel array and get the non-empty
    frame-to-slice indices. """
    
    Nda = sitk.GetArrayFromImage(Image)
    
    F2Sinds = []

    for i in range(Nda.shape[0]):
        Sum = np.amax(Nda[i])

        if Sum:
            F2Sinds.append(i)
    
    PixArr = np.zeros((len(F2Sinds), Nda.shape[1], Nda.shape[2]))
    
    for i in range(len(F2Sinds)):
        PixArr[i] = Nda[F2Sinds[i]]    
    
    return PixArr, F2Sinds


def PlotTwoImages(Title0, Im0, Ind0, Title1, Im1, Ind1):
    PixArr0 = sitk.GetArrayViewFromImage(Im0)[Ind0,:,:]
    PixArr1 = sitk.GetArrayViewFromImage(Im1)[Ind1,:,:]
    
    z0 = Im0.TransformIndexToPhysicalPoint([0,0,Ind0])[2]
    z1 = Im1.TransformIndexToPhysicalPoint([0,0,Ind1])[2]
    
    Title0 = f'{Title0}\nz = {round(z0, 2)} mm'
    Title1 = f'{Title1}\nz = {round(z1, 2)} mm'
    
    plt.subplots(1, 2, figsize=(15,8))
    
    plt.subplot(1, 2, 1)
    plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
    plt.title(Title0)
    
    plt.subplot(1, 2, 2)
    plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
    plt.title(Title1)
    
    plt.show()
    
    return


def PlotThreeImages(Title0, Im0, Ind0, Title1, Im1, Ind1, Title2, Im2, Ind2):
    PixArr0 = sitk.GetArrayViewFromImage(Im0)[Ind0,:,:]
    PixArr1 = sitk.GetArrayViewFromImage(Im1)[Ind1,:,:]
    PixArr2 = sitk.GetArrayViewFromImage(Im2)[Ind2,:,:]
    
    z0 = Im0.TransformIndexToPhysicalPoint([0,0,Ind0])[2]
    z1 = Im1.TransformIndexToPhysicalPoint([0,0,Ind1])[2]
    z2 = Im2.TransformIndexToPhysicalPoint([0,0,Ind2])[2]
    
    Title0 = f'{Title0}\nz = {round(z0, 2)} mm'
    Title1 = f'{Title1}\nz = {round(z1, 2)} mm'
    Title2 = f'{Title2}\nz = {round(z2, 2)} mm'
    
    plt.subplots(1, 3, figsize=(15,8))
    
    plt.subplot(1, 3, 1)
    plt.imshow(PixArr0, cmap=plt.cm.Greys_r);
    plt.title(Title0)
    
    plt.subplot(1, 3, 2)
    plt.imshow(PixArr1, cmap=plt.cm.Greys_r);
    plt.title(Title1)
    
    plt.subplot(1, 3, 3)
    plt.imshow(PixArr2, cmap=plt.cm.Greys_r);
    plt.title(Title2)
    
    plt.show()
    
    return


def PlotRegResults(FixIm, MovIm, RegIm, FixInd, FixTitle='Fixed image', 
                   MovTitle='Moving image', RegTitle='Registered image'):
    
    MovInd, zDiff = GetIndOfNearestSlice(FixIm, FixInd, MovIm, True)
    
    FixPixArr = sitk.GetArrayViewFromImage(FixIm)[FixInd,:,:]
    MovPixArr = sitk.GetArrayViewFromImage(MovIm)[MovInd,:,:]
    RegPixArr = sitk.GetArrayViewFromImage(RegIm)[FixInd,:,:]
    
    FixZ = FixIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    MovZ = MovIm.TransformIndexToPhysicalPoint([0,0,MovInd])[2]
    RegZ = RegIm.TransformIndexToPhysicalPoint([0,0,FixInd])[2]
    
    FixTitle = f'{FixTitle}\nz = {round(FixZ, 2)} mm'
    MovTitle = f'{MovTitle}\nz = {round(MovZ, 2)} mm'
    RegTitle = f'{RegTitle}\nz = {round(RegZ, 2)} mm'
        
    plt.subplots(2, 3, figsize=(18,12), dpi=80)
    
    alpha = 0.5
    
    BlendedIm = (1.0 - alpha)*FixIm[:,:,FixInd] + alpha*RegIm[:,:,FixInd]
    
    DiffIm = FixIm[:,:,FixInd] - RegIm[:,:,FixInd]
    
    plt.subplot(2, 3, 1)
    plt.imshow(MovPixArr, cmap=plt.cm.Greys_r);
    plt.title(MovTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 2)
    plt.imshow(FixPixArr, cmap=plt.cm.Greys_r);
    plt.title(FixTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 3)
    plt.imshow(RegPixArr, cmap=plt.cm.Greys_r);
    plt.title(RegTitle)
    plt.axis('off')
    
    plt.subplot(2, 3, 4)
    plt.axis('off')
    
    plt.subplot(2, 3, 5)
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.title('Blended image')
    plt.axis('off')
    
    plt.subplot(2, 3, 6)
    plt.imshow(sitk.GetArrayViewFromImage(DiffIm), cmap=plt.cm.Greys_r);
    plt.title('Difference image')
    plt.axis('off')
    
    plt.show()
    
    return




"""
Modify the directory and file paths accordingly:

SrcDcmDir = r'C:\.....\2021-03-25_Cropped_segments\Ses4_Ser11_DICOMs'
TrgDcmDir = r'C:\.....\2021-03-25_Cropped_segments\Ses3_Ser6_DICOMs'

SrcSegFpath = r'C:\.....\2021-03-25_Cropped_segments\Ses4_Ser11_Ventricles_SEG_20210213_144855.dcm'

SrcDcmIm = ImportDicomSeriesAsSitkImage(SrcDcmDir)
TrgDcmIm = ImportDicomSeriesAsSitkImage(TrgDcmDir)

SrcLabIm = ImportLabIm(SrcSegFpath, SrcDcmDir)

The Source and Target DICOM images can be viewed as follows:

interact(PlotTwoImages,
         Title0='Source DICOM image',
         Im0=fixed(SrcDcmIm),
         Ind0=(0, SrcDcmIm.GetSize()[2]-1),
         Title1='Target DICOM image',
         Im1=fixed(TrgDcmIm),
         Ind1=(0, TrgDcmIm.GetSize()[2]-1));

         
Next register SrcDcmIm to TrgDcmIm:

RegDcmIm, FinalTx = RegisterImages(FixIm=TrgDcmIm, MovIm=SrcDcmIm, FinalInterp='NearestNeighbor')

Display the result:

PlotRegResults(TrgDcmIm, SrcDcmIm, RegDcmIm, 19)

interact(PlotThreeImages,
         Title0='Source DICOM image',
         Im0=fixed(SrcDcmIm),
         Ind0=(0, SrcDcmIm.GetSize()[2]-1),
         Title1='Target DICOM image',
         Im1=fixed(TrgDcmIm),
         Ind1=(0, TrgDcmIm.GetSize()[2]-1),
         Title2='Registered Source DICOM image',
         Im2=fixed(RegDcmIm),
         Ind2=(0, RegDcmIm.GetSize()[2]-1));


Finally, resample SrcLabIm using FinalTx:

TxSrcLabIm = sitk.Resample(SrcLabIm, TrgDcmIm, FinalTx, sitk.sitkNearestNeighbor, 0.0, SrcLabIm.GetPixelID())

And display the result:

SrcPixArr, SrcLabF2Sinds = Image2PixArr(SrcLabIm)
TxSrcPixArr, TxSrcLabF2Sinds = Image2PixArr(TxSrcLabIm)

a = SrcLabF2Sinds[0]
b = SrcLabF2Sinds[-1]

c = TxSrcLabF2Sinds[0]
d = TxSrcLabF2Sinds[-1]

interact(PlotTwoImages,
         Title0='Source label image',
         Im0=fixed(SrcLabIm),
         Ind0=(a, b),
         Title1='Transformed label image',
         Im1=fixed(TxSrcLabIm),
         Ind1=(c, d));
"""
