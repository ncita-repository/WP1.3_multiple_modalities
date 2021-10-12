# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:35:37 2021

@author: ctorti
"""


import os
import time
from datetime import datetime
from pydicom import dcmread
from pydicom.uid import generate_uid
import SimpleITK as sitk
import numpy as np

from dro_tools.general import modify_ref_im_seq#, modify_ref_ser_seq
#from dro_tools.general import modify_stu_con_oth_ref_ins_seq
#from dro_tools.matrices import is_matrix_orthonormal, is_matrix_orthogonal
from dro_tools.matrices import get_tx_matrix_type
#from dicom_tools.imports import import_dicoms
from general_tools.general import reduce_list_of_str_floats_to_16

def create_def_dro_OBSOLETE(
        srcDicoms, trgDicoms, gridOrig, gridDir, gridDims, gridRes, 
        vectGridData, sampleDroDir, txMatrix=None, description='', p2c=False
        ):
    """
    SEE: create_def_dro() in dro_tools.create_dro.py
    
    Create a DICOM Deformable Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    srcDicoms : list of Pydicom Objects
        List of Source DICOMs.
    trgDicoms : list of Pydicom Objects
        List of Target DICOMs.
    gridOrig : list of strs
        List of float strings representing the origin of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    gridDir : list of strs
        List of float strings representing the direction of the BSpline grid 
        used to deform the moving (Source) image to the fixed (Target) image.
    gridDims : list of ints
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    gridRes : list of floats
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image.
    vectGridData : list of strs
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image.
    sampleDroDir : str
        Directory path to the sample DROs.
    txMatrix : list of strs or None, optional
        List of strings representing the non-deformable transformation 
        parameters to be applied prior to the bspline registration. txMatrix
        will be stored in PreDeformationMatrixRegistrationSequence.  If None, 
        the identity matrix will be used. Default value is None.
    description : str, optional
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO. Default value is ''.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        Default value is False.
    
    Returns
    -------
    dro : Pydicom object
        The DICOM Registration Object.
    
    Notes
    -----
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    """
    
    #cwd = os.getcwd()

    #droDir = os.path.join(cwd, 'sample_DROs')
    
    droFname = 'sample_deformable_dro.dcm'
    #droFpath = os.path.join(droDir, droFname)
    droFpath = os.path.join(sampleDroDir, droFname)
        
    #""" Check if the sample DRO has already been downloaded: """   
    #if os.path.isfile(droFpath):
    #    print(f'Sample DRO already downloaded to:\n {droDir}\n')
    #else:
    #    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    #    print(f'Downloading DROs from {url}...\n')
    #    DownloadSampleDros()
    
    """
    Decide whether to reference all SOPInstanceUIDs in srcDicoms and trgDicoms
    within ReferencedImageSequence (not required), or just one (for compactness).
    """
    #refAllSOPs = True
    refAllSOPs = False
    
    if p2c:
        if refAllSOPs:
            print('Note: All SOPInstanceUIDs of both series will be referenced',
                  'in DeformableRegistrationSequence[0].ReferencedImageSequence.')
        else:
            print('Note: Only the first SOPInstanceUID of both series will be',
                  'referenced in',
                  'DeformableRegistrationSequence[0].ReferencedImageSequence.')
    
    # Start timing:
    times = []
    times.append(time.time())
    
    """
    srcDicoms = import_dcms(srcDicomDir)
    trgDicoms = import_dcms(trgDicomDir)
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n*Took {dtime} s to import the source and target DICOMs.')
    """
    
    # Read in the DRO template.
    dro = dcmread(droFpath)
    
    #currentDate = time.strftime("%Y%m%d", time.gmtime())
    #currentTime = time.strftime("%H%M%S", time.gmtime())
    ##currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    timeNow = datetime.datetime.now()
    currentDate = timeNow.strftime('%Y%m%d')
    currentTime = timeNow.strftime('%H%M%S.%f')
    
    dro.InstanceCreationDate = currentDate
    dro.InstanceCreationTime = currentTime
    
    # Generate a new SOPInstanceUID:
    dro.SOPInstanceUID = generate_uid()
    dro.file_meta.MediaStorageSOPInstanceUID = dro.SOPInstanceUID
    
    dro.ContentDate = srcDicoms[0].ContentDate
    dro.ContentTime = srcDicoms[0].ContentTime
    dro.StudyDate = srcDicoms[0].StudyDate
    try:
        seriesDate = srcDicoms[0].SeriesDate
        try:
            dro.SeriesDate = seriesDate
        except AttributeError:
            pass
    except AttributeError:
        pass
    #dro.ContentDate = srcDicoms[0].ContentDate
    dro.ContentDate = currentDate
    dro.StudyTime = srcDicoms[0].StudyTime
    try:
        seriesTime = srcDicoms[0].SeriesTime
        try:
            dro.SeriesTime = seriesTime
        except AttributeError:
            pass
    except AttributeError:
        pass
    #dro.ContentTime = srcDicoms[0].ContentTime
    dro.ContentTime = currentTime
    dro.AccessionNumber = ''
    #dro.Manufacturer = srcDicoms[0].Manufacturer
    dro.Manufacturer = 'NCITA'
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    try:
        studyDesc = srcDicoms[0].StudyDescription
        try:
            dro.StudyDescription = studyDesc
        except AttributeError:
            pass
    except AttributeError:
        pass
    try:
        seriesDesc = srcDicoms[0].SeriesDescription
        try:
            dro.SeriesDescription = seriesDesc
        except AttributeError:
            pass
    except AttributeError:
        pass
    #dro.ManufacturerModelName = srcDicoms[0].ManufacturerModelName
    dro.ManufacturerModelName = ''
    dro.PatientName = srcDicoms[0].PatientName
    dro.PatientID = srcDicoms[0].PatientID
    dro.PatientBirthDate = srcDicoms[0].PatientBirthDate
    dro.PatientSex = srcDicoms[0].PatientSex
    try:
        patientAge = srcDicoms[0].PatientAge
        try:
            dro.PatientAge = patientAge
        except AttributeError:
            pass
    except AttributeError:
        pass
    """ 
    Keep De-identificationMethod as it is, but consider modifying to 
    represent the actual methods used. 

    Delete PrivateCreator and PrivateTagData tags (not in the DICOM 
    standard): 
    """
    del dro[0x0013, 0x0010]
    del dro[0x0013, 0x1010]
    del dro[0x0013, 0x1013]
    dro.SoftwareVersions = ''
    dro.StudyInstanceUID = srcDicoms[0].StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    dro.SeriesInstanceUID = generate_uid()
    dro.StudyID = srcDicoms[0].StudyID
    dro.SeriesNumber = srcDicoms[0].SeriesNumber
    dro.InstanceNumber = srcDicoms[0].InstanceNumber
    dro.FrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
    dro.PositionReferenceIndicator = srcDicoms[0].PositionReferenceIndicator
    """ Keep LongitudinalTemporalInformationModified as 'MODIFIED'. """
    
    
    
    """ 
    Modify the ReferencedSeriesSequence. 
    
    Note:
        ReferencedSeriesSequence will reference the moving (source) series.
    """
    #dro = ModifyRefSerSeq(dro, SeqNum=0, Dicoms=srcDicoms, 
    #                      refAllSOPs=refAllSOPs, p2c=p2c)
    
    seq = modify_ref_im_seq(
        seq=dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence,
        dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
    
    dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence = seq
    
    dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = srcDicoms[0].SeriesInstanceUID
    
    
    """ 
    Modify the StudiesContainingOtherReferencedInstancesSequence. 
    
    Note:
        StudiesContainingOtherReferencedInstancesSequence will reference the
        fixed (target) series.
    """
    #dro = modify_stu_con_oth_ref_ins_seq(dro, seqNum=0, dicoms=trgDicoms, 
    #                               refAllSOPs=refAllSOPs, 
    #                               p2c=p2c)
    
    seq = modify_ref_im_seq(
        seq=dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value,
        dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
    
    dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value = seq
    
    dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = trgDicoms[0].SeriesInstanceUID
    
    dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .StudyInstanceUID = trgDicoms[0].StudyInstanceUID
    
    
    
    """ Modify the DeformableRegistrationSequence. 
    
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
    if p2c:
        print('\nModify the DeformableRegistrationSequence.')
    
    #""" Remove the first sequence in DeformableRegistrationSequence (which
    #references the fixed image series). """
    #dro.DeformableRegistrationSequence.pop(0)
    
    
    """ 
    Modify the first ReferencedImageSequence. 
    
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #dro = modify_ref_ser_seq(dro, seqNum=0, Dicoms=trgDicoms, 
    #                      refAllSOPs=refAllSOPs, p2c=p2c)
    seq = modify_ref_im_seq(
        seq=dro.DeformableRegistrationSequence[0].ReferencedImageSequence,
        dicoms=trgDicoms, refAllSOPs=True, p2c=p2c
        )
    
    dro.DeformableRegistrationSequence[0].ReferencedImageSequence = seq
    
    
    """ > Modify the SourceFrameOfReferenceUID.
    Note:
    
    The DICOM standard defines registration from Source to Registered RCS 
    (Referenced Coordinate System).
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_A.39.html
    """
    dro.DeformableRegistrationSequence[0]\
       .SourceFrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
    
    
    # Delete MatrixRegistrationSequence.
    del dro[0x64, 0x02][0][0x70, 0x0309]
    
    
    
    
    
    
    """ 
    Modify the second ReferencedImageSequence. 
    
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #dro = modify_ref_ser_seq(dro, seqNum=1, Dicoms=srcDicoms, 
    #                      refAllSOPs=refAllSOPs, p2c=p2c)
    seq = modify_ref_im_seq(
        seq=dro.DeformableRegistrationSequence[1].ReferencedImageSequence,
        dicoms=srcDicoms, refAllSOPs=True, p2c=p2c
        )
    
    dro.DeformableRegistrationSequence[1].ReferencedImageSequence = seq
    
    
    """ > Modify the SourceFrameOfReferenceUID. """
    dro.DeformableRegistrationSequence[1]\
       .SourceFrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
    
    
    """ > Modify the DeformableRegistrationGridSequence. """
    dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImagePositionPatient = reduce_list_of_str_floats_to_16(gridOrig)
    
    if len(gridDir) > 6:
        # Ignore the direction cosines along the z-direction:
        gridDir = gridDir[0:6]
        
    dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImageOrientationPatient = reduce_list_of_str_floats_to_16(gridDir)
    
    dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridDimensions = gridDims
    
    dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridResolution = gridRes
    
    dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .VectorGridData = vectGridData
    
    
    
    """ > Modify PreDeformationMatrixRegistrationSequence. """
    dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrix = txMatrix
    
    """ > Modify the FrameOfReferenceTransformationMatrixType.  Acceptable 
    values:
    'RIGID' (= value in template DRO)
    'RIGID_SCALE' (Similarity transform)
    'AFFINE' 
    
    Note:  Need to verify that get_tx_matrix_type covers all 3 options.
    """
    
    #if txParams != None: # 04/06/21
    #    print(f'txParams = {txParams}')
    #    
    #    dro.DeformableRegistrationSequence[1]\
    #       .PreDeformationMatrixRegistrationSequence[0]\
    #       .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(txParams, 
    #                                                                   p2c)
    #
    #    times.append(time.time())
    #    dtime = round(times[-1] - times[-2], 3)
    #    if p2c:
    #        print(f'\n   *Took {dtime} s to get the transform matrix type.')
    
    # 04/06/21:
    dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(
           txMatrix, p2c
           )
           
    """ > Consider adding details for the FrameOfReferenceTransformationComment
    and RegistrationTypeCodeSequence (optional). 
    dro.DeformableRegistrationSequence[1]\
       .FrameOfReferenceTransformationComment = ''
    
    dro.DeformableRegistrationSequence[1]\
       .RegistrationTypeCodeSequence = ''
    """
    
    
    """ > Leave PostDeformationMatrixRegistrationSequence unchanged. """
    
    
    """ Keep ContentLabel as 'REGISTRATION'. """
    dro.ContentCreatorName = 'NCITA'
    # Delete PrivateCreator tag (not in the DICOM standard).
    del dro[0x3773, 0x01]
    
    
    """ > Consider adding UsedFiducialsSequence (optional). 
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
    
    Spatial Fiducials Series Module:
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
    """
    
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
        
    
    dtime = round(times[-1] - times[0], 1)
    if p2c:
        print(f'\n   *Took a total of {dtime} s to create the DRO.')
        
        
    return dro

def create_def_dro_from_tx_OBSOLETE(
        srcDicoms, trgDicoms, regTx, sampleDroDir, preRegTx=None, 
        description='', p2c=False
        ):
    """
    SEE: create_def_dro() in dro_tools.create_dro.py
    
    Create a DICOM Deformable Spatial Registration Object (DRO) from a 
    SimpleITK BSpline transformation.  
    
    Parameters
    ----------
    srcDicoms : list of Pydicom Objects
        List of Source DICOMs.
    trgDicoms : list of Pydicom Objects
        List of Target DICOMs.
    regTx : SimpleITK Transform
        The registration transform.
    sampleDroDir : str
        Directory path to the sample DROs.
    preRegTx : SimpleITK Transform or None, optional
        The transform used prior to registration, and to be stored in 
        PreDeformationMatrixRegistrationSequence.  If None, the identity 
        matrix will be stored in PreDeformationMatrixRegistrationSequence.
        The default value is None.
    description : str, optional ('' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    dro : Pydicom Object
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
    if type(regTx) == sitk.SimpleITK.Transform:
        regTx = sitk.BSplineTransform(regTx)
    
    coeffImX, coeffImY, coeffImZ = regTx.GetCoefficientImages()
    
    # Zip the vectors to a list of x, y, z tuples:
    listOfTuples = list(zip(
        list(sitk.GetArrayViewFromImage(coeffImX).flatten()),
        list(sitk.GetArrayViewFromImage(coeffImY).flatten()),
        list(sitk.GetArrayViewFromImage(coeffImZ).flatten())
        ))
    
    flatListOfVects = [item for tup in listOfTuples for item in tup]
    
    vectGridData = np.array(flatListOfVects).tobytes()
    
    gridOrig = [str(item) for item in coeffImX.GetOrigin()]
    gridDir = [str(item) for item in coeffImX.GetDirection()]
    gridDims = list(coeffImX.GetSize())
    gridRes = list(coeffImX.GetSpacing())
    
    if len(gridDir) > 6:
        # Ignore the direction cosines along the z-direction:
        gridDir = gridDir[0:6]
                                    
    print(f'type(gridOrig) = {type(gridOrig)}')
    print(f'gridOrig = {gridOrig}\n')
    
    print(f'type(gridDir) = {type(gridDir)}')
    print(f'gridDir = {gridDir}\n')
    
    print(f'type(gridDims) = {type(gridDims)}')
    print(f'gridDims = {gridDims}\n')
    
    print(f'type(gridRes) = {type(gridRes)}')
    print(f'gridRes = {gridRes}\n')
    
    print(f'type(vectGridData) = {type(vectGridData)}')
    #print(f'vectGridData = {vectGridData}\n')
    
    # Get the transformation matrix (txMatrix) from the pre-registration
    # transform parameters: """
    if preRegTx == None:
        #preRegTxParams = None
        
        txMatrix = ['1.0', '0.0', '0.0', '0.0',
                    '0.0', '1.0', '0.0', '0.0',
                    '0.0', '0.0', '1.0', '0.0',
                    '0.0', '0.0', '0.0', '1.0']
    else:
        txParams = preRegTx.GetParameters()
        
        # Add the row vector [0, 0, 0, 1] to txParams to complete the Frame of
        # Reference Transformation Matrix 
        # (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
        #txParams = [str(item) for item in txParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
        txMatrix = [
            str(txParams[0]), str(txParams[1]), str(txParams[2]), '0.0',
            str(txParams[3]), str(txParams[4]), str(txParams[5]), '0.0',
            str(txParams[6]), str(txParams[7]), str(txParams[8]), '0.0',
            str(txParams[9]), str(txParams[10]), str(txParams[11]), '1.0'
            ]
    
    dro = create_def_dro(
        srcDicoms, trgDicoms, gridOrig, gridDir, gridDims, gridRes, 
        vectGridData, sampleDroDir, txMatrix, description, p2c
        )
    
    return dro