# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:35:37 2021

@author: ctorti
"""


""" Functions that create DROs from SimpleITK Transforms. """

import SimpleITK as sitk
import numpy as np
import os
import time
from datetime import datetime
from pydicom import dcmread
from pydicom.uid import generate_uid
#from copy import deepcopy

from dro_tools.general import modify_ref_im_seq#, modify_ref_ser_seq
#from dro_tools.general import modify_stu_con_oth_ref_ins_seq
#from dro_tools.matrices import is_matrix_orthonormal, is_matrix_orthogonal
from dro_tools.matrices import get_tx_matrix_type
from dicom_tools.imports import import_dicoms


from general_tools.general import reduce_list_of_str_floats_to_16

def create_def_dro(SrcDicomDir, TrgDicomDir, 
                   GridOrig, GridDir, GridDims, GridRes, VectGridData,
                   TxMatrix=None, Description='', LogToConsole=False):
    """
    Create a DICOM Deformable Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    SrcDicomDir : str
        Directory containing Source DICOMs.
    TrgDicomDir : str
        Directory containing Target DICOMs.
    GridOrig : list of float strs
        List of float strings representing the origin of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    GridDir : list of float strs
        List of float strings representing the direction of the BSpline grid 
        used to deform the moving (Source) image to the fixed (Target) image.
    GridDims : list of ints
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    GridRes : list of floats
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image.
    VectGridData : list of float strs
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image.
    TxMatrix : list of strs or None, optional (None by default)
        List of strings representing the non-deformable transformation 
        parameters to be applied prior to the bspline registration. TxMatrix
        will be stored in PreDeformationMatrixRegistrationSequence.  If None, 
        the identity matrix will be used.
    Description : str, optional ('' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    Dro : Pydicom Object
        The DICOM Registration Object.
    
    Note
    ----
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    """
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DROs')
    
    DroFname = 'sample_deformable_dro.dcm'
    DroFpath = os.path.join(DroDir, DroFname)
    
    #""" Check if the sample DRO has already been downloaded: """   
    #if os.path.isfile(DroFpath):
    #    print(f'Sample DRO already downloaded to:\n {DroDir}\n')
    #else:
    #    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    #    print(f'Downloading DROs from {url}...\n')
    #    DownloadSampleDros()
    
    """ Decide whether to reference all SOPInstanceUIDs in SrcDicoms and 
    TrgDicoms within ReferencedImageSequence (not required), or just one
    (for compactness). """
    #RefAllSOPs = True
    RefAllSOPs = False
    
    if LogToConsole:
        if RefAllSOPs:
            print('Note: All SOPInstanceUIDs of both series will be referenced',
                  'in DeformableRegistrationSequence[0].ReferencedImageSequence.')
        else:
            print('Note: Only the first SOPInstanceUID of both series will be',
                  'referenced in',
                  'DeformableRegistrationSequence[0].ReferencedImageSequence.')
    
    # Start timing:
    times = []
    times.append(time.time())
    
    SrcDicoms = import_dicoms(SrcDicomDir)
    TrgDicoms = import_dicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
    
    # Read in the DRO template:
    Dro = dcmread(DroFpath)
    
    #CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    #CurrentTime = time.strftime("%H%M%S", time.gmtime())
    ##CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    TimeNow = datetime.datetime.now()
    CurrentDate = TimeNow.strftime('%Y%m%d')
    CurrentTime = TimeNow.strftime('%H%M%S.%f')
    
    Dro.InstanceCreationDate = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    Dro.ContentDate = SrcDicoms[0].ContentDate
    Dro.ContentTime = SrcDicoms[0].ContentTime
    Dro.StudyDate = SrcDicoms[0].StudyDate
    try:
        SeriesDate = SrcDicoms[0].SeriesDate
        try:
            Dro.SeriesDate = SeriesDate
        except AttributeError:
            pass
    except AttributeError:
        pass
    #Dro.ContentDate = SrcDicoms[0].ContentDate
    Dro.ContentDate = CurrentDate
    Dro.StudyTime = SrcDicoms[0].StudyTime
    try:
        SeriesTime = SrcDicoms[0].SeriesTime
        try:
            Dro.SeriesTime = SeriesTime
        except AttributeError:
            pass
    except AttributeError:
        pass
    #Dro.ContentTime = SrcDicoms[0].ContentTime
    Dro.ContentTime = CurrentTime
    Dro.AccessionNumber = ''
    #Dro.Manufacturer = SrcDicoms[0].Manufacturer
    Dro.Manufacturer = 'NCITA'
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    try:
        StudyDesc = SrcDicoms[0].StudyDescription
        try:
            Dro.StudyDescription = StudyDesc
        except AttributeError:
            pass
    except AttributeError:
        pass
    try:
        SeriesDesc = SrcDicoms[0].SeriesDescription
        try:
            Dro.SeriesDescription = SeriesDesc
        except AttributeError:
            pass
    except AttributeError:
        pass
    #Dro.ManufacturerModelName = SrcDicoms[0].ManufacturerModelName
    Dro.ManufacturerModelName = ''
    Dro.PatientName = SrcDicoms[0].PatientName
    Dro.PatientID = SrcDicoms[0].PatientID
    Dro.PatientBirthDate = SrcDicoms[0].PatientBirthDate
    Dro.PatientSex = SrcDicoms[0].PatientSex
    try:
        PatientAge = SrcDicoms[0].PatientAge
        try:
            Dro.PatientAge = PatientAge
        except AttributeError:
            pass
    except AttributeError:
        pass
    """ Keep De-identificationMethod as it is, but consider modifying to 
    represent the actual methods used. """
    """ Delete PrivateCreator and PrivateTagData tags (not in the DICOM 
    standard): """
    del Dro[0x0013, 0x0010]
    del Dro[0x0013, 0x1010]
    del Dro[0x0013, 0x1013]
    Dro.SoftwareVersions = ''
    Dro.StudyInstanceUID = SrcDicoms[0].StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    Dro.SeriesInstanceUID = generate_uid()
    Dro.StudyID = SrcDicoms[0].StudyID
    Dro.SeriesNumber = SrcDicoms[0].SeriesNumber
    Dro.InstanceNumber = SrcDicoms[0].InstanceNumber
    Dro.FrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    Dro.PositionReferenceIndicator = SrcDicoms[0].PositionReferenceIndicator
    """ Keep LongitudinalTemporalInformationModified as 'MODIFIED'. """
    
    # Modify the ReferencedSeriesSequence:
    """
    Note:
        ReferencedSeriesSequence will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=SrcDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    Seq = modify_ref_im_seq(Dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence,
                         Dicoms=SrcDicoms, RefAllSOPs=RefAllSOPs, 
                         LogToConsole=LogToConsole)
    
    Dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence = Seq
    
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
    
    # Modify the StudiesContainingOtherReferencedInstancesSequence. 
    """
    Note:
        StudiesContainingOtherReferencedInstancesSequence will reference the
        fixed (target) series.
    """
    #Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
    #                               RefAllSOPs=RefAllSOPs, 
    #                               LogToConsole=LogToConsole)
    Seq = modify_ref_im_seq(Dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value,
                         Dicoms=TrgDicoms, RefAllSOPs=RefAllSOPs, 
                         LogToConsole=LogToConsole)
    
    Dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value = Seq
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
    
    # Modify the DeformableRegistrationSequence:
    """
    Note:
        The sample DRO contains two sequences in DeformableRegistrationSequence
        - the first relates to the fixed image series, the second to the moving
        image series. 
        
        The DICOM standard seems to suggest that only the series to be
        registered is to be referenced in DeformableRegistrationSequence:
        
        https://dicom.innolitics.com/ciods/deformable-spatial-registration/deformable-spatial-registration/00640002
        
        The first sequence even contains MatrixRegistrationSequence, which is
        not a child of DeformableRegistrationSequence.  It seems that they 
        wanted to reference both moving and fixed series, but it's not clear 
        why since the moving series was already fully referenced in 
        ReferencedSeriesSequence, and the fixed in 
        StudiesContainingOtherReferencedInstanceSequence.
        
        But since I've only found examples of DROs (spatial and deformable)
        that contain two sequences in DeformableRegistrationSequence (or
        RegistrationSequence for spatial DROs) I will maintain this structure
        in case there are best-practices that are not necessarily to the 
        standard.
    """
    if LogToConsole:
        print('\nModify the DeformableRegistrationSequence.')
    
    #""" Remove the first sequence in DeformableRegistrationSequence (which
    #references the fixed image series). """
    #Dro.DeformableRegistrationSequence.pop(0)
    
    
    # Modify the first ReferencedImageSequence:
    """
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    Seq = modify_ref_im_seq(Dro.DeformableRegistrationSequence[0].ReferencedImageSequence,
                         Dicoms=TrgDicoms, RefAllSOPs=True, 
                         LogToConsole=LogToConsole)
    
    Dro.DeformableRegistrationSequence[0].ReferencedImageSequence = Seq
    
    # > Modify the SourceFrameOfReferenceUID:
    """
    Note:
    The DICOM standard defines registration from Source to Registered RCS 
    (Referenced Coordinate System).
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_A.39.html
    """
    Dro.DeformableRegistrationSequence[0]\
       .SourceFrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    
    # Delete MatrixRegistrationSequence:
    del Dro[0x64, 0x02][0][0x70, 0x0309]
    
    # Modify the second ReferencedImageSequence. 
    """
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    Seq = modify_ref_im_seq(Dro.DeformableRegistrationSequence[1].ReferencedImageSequence,
                         Dicoms=SrcDicoms, RefAllSOPs=True, 
                         LogToConsole=LogToConsole)
    
    Dro.DeformableRegistrationSequence[1].ReferencedImageSequence = Seq
    
    # > Modify the SourceFrameOfReferenceUID:
    Dro.DeformableRegistrationSequence[1]\
       .SourceFrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    
    # > Modify the DeformableRegistrationGridSequence:
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImagePositionPatient = reduce_list_of_str_floats_to_16(GridOrig)
    
    if len(GridDir) > 6:
        # Ignore the direction cosines along the z-direction:
        GridDir = GridDir[0:6]
        
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImageOrientationPatient = reduce_list_of_str_floats_to_16(GridDir)
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridDimensions = GridDims
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridResolution = GridRes
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .VectorGridData = VectGridData
    
    # > Modify PreDeformationMatrixRegistrationSequence:
    Dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxMatrix
    
    # > Modify the FrameOfReferenceTransformationMatrixType:
    """ 
    Acceptable values:
        'RIGID' (= value in template DRO)
        'RIGID_SCALE' (Similarity transform)
        'AFFINE' 
    Note:  Need to verify that get_tx_matrix_type covers all 3 options.
    """
    
    #if TxParams != None: # 04/06/21
    #    print(f'TxParams = {TxParams}')
    #    
    #    Dro.DeformableRegistrationSequence[1]\
    #       .PreDeformationMatrixRegistrationSequence[0]\
    #       .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(TxParams, 
    #                                                                   LogToConsole)
    #
    #    times.append(time.time())
    #    Dtime = round(times[-1] - times[-2], 3)
    #    if LogToConsole:
    #        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
    
    # 04/06/21:
    Dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrixType\
           = get_tx_matrix_type(TxMatrix, LogToConsole)
           
    """ > Consider adding details for the FrameOfReferenceTransformationComment
    and RegistrationTypeCodeSequence (optional). 
    Dro.DeformableRegistrationSequence[1]\
       .FrameOfReferenceTransformationComment = ''
    
    Dro.DeformableRegistrationSequence[1]\
       .RegistrationTypeCodeSequence = ''
    """
    
    """ > Leave PostDeformationMatrixRegistrationSequence unchanged. """
    
    """ Keep ContentLabel as 'REGISTRATION'. """
    Dro.ContentCreatorName = 'NCITA'
    """ Delete PrivateCreator tag (not in the DICOM standard). """
    del Dro[0x3773, 0x01]
    
    """ > Consider adding UsedFiducialsSequence (optional). 
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
    
    Spatial Fiducials Series Module:
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
    """
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    Dtime = round(times[-1] - times[0], 1)
    if LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
    
    return Dro

def create_def_dro_from_tx(SrcDicomDir, TrgDicomDir, BsplineTx, PreRegTx=None,
                           Description='', LogToConsole=False):
    """
    Create a DICOM Deformable Spatial Registration Object (DRO) from a 
    SimpleITK BSpline transformation.  
    
    Parameters
    ----------
    SrcDicomDir : str
        Directory containing Source DICOMs.
    TrgDicomDir : str
        Directory containing Target DICOMs.
    BsplineTx : SimpleITK Transform
        The final transform following a deformable image registration with a
        BSpline Transform.
    PreRegTx : SimpleITK Transform or None, optional (None by default)
        The transform used prior to registration, and to be stored in 
        PreDeformationMatrixRegistrationSequence.  If None, the identity 
        matrix will be stored in PreDeformationMatrixRegistrationSequence.
    Description : str, optional ('' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    Dro : Pydicom Object
        The DICOM Registration Object.
    
    Note
    ----
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    VectorGridData has VR = OF:
    
    A stream of 32-bit IEEE 754:1985 floating point words. OF is a VR that 
    requires byte swapping within each 32-bit word when changing byte ordering 
    (see Section 7.3).
    
    http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html
    """
    
    """ The final transform after executing the registration is of class
    SimpleITK.SimpleITK.Transform but needs to be of class
    SimpleITK.SimpleITK.BSplineTransform. Check and change if needed: """
    if type(BsplineTx) == sitk.SimpleITK.Transform:
        BsplineTx = sitk.BSplineTransform(BsplineTx)
    
    CoeffImX, CoeffImY, CoeffImZ = BsplineTx.GetCoefficientImages()
    
    # Zip the vectors to a list of x, y, z tuples:
    ListOfTuples = list(zip(list(sitk.GetArrayViewFromImage(CoeffImX).flatten()),
                            list(sitk.GetArrayViewFromImage(CoeffImY).flatten()),
                            list(sitk.GetArrayViewFromImage(CoeffImZ).flatten())))
    
    FlatListOfVects = [item for tup in ListOfTuples for item in tup]
    
    VectGridData = np.array(FlatListOfVects).tobytes()
    
    GridOrig = [str(item) for item in CoeffImX.GetOrigin()]
    GridDir = [str(item) for item in CoeffImX.GetDirection()]
    GridDims = list(CoeffImX.GetSize())
    GridRes = list(CoeffImX.GetSpacing())
    
    if len(GridDir) > 6:
        # Ignore the direction cosines along the z-direction:
        GridDir = GridDir[0:6]
                                    
    print(f'type(GridOrig) = {type(GridOrig)}')
    print(f'GridOrig = {GridOrig}\n')
    
    print(f'type(GridDir) = {type(GridDir)}')
    print(f'GridDir = {GridDir}\n')
    
    print(f'type(GridDims) = {type(GridDims)}')
    print(f'GridDims = {GridDims}\n')
    
    print(f'type(GridRes) = {type(GridRes)}')
    print(f'GridRes = {GridRes}\n')
    
    print(f'type(VectGridData) = {type(VectGridData)}')
    #print(f'VectGridData = {VectGridData}\n')
    
    # Get the transformation matrix (TxMatrix) from the pre-registration
    # transform parameters: """
    if PreRegTx == None:
        #PreRegTxParams = None
        
        TxMatrix = ['1.0', '0.0', '0.0', '0.0',
                    '0.0', '1.0', '0.0', '0.0',
                    '0.0', '0.0', '1.0', '0.0',
                    '0.0', '0.0', '0.0', '1.0']
    else:
        TxParams = PreRegTx.GetParameters()
        
        # Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
        # Reference Transformation Matrix 
        # (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
        #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
        TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
                    str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
                    str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
                    str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    Dro = create_def_dro(SrcDicomDir, TrgDicomDir, GridOrig, GridDir, GridDims, 
                         GridRes, VectGridData, TxMatrix, Description, 
                         LogToConsole)
    
    return Dro