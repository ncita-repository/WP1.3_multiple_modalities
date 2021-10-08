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

#from dro_tools.general import modify_ref_im_seq
from dro_tools.general import modify_ref_ser_seq
from dro_tools.general import modify_stu_con_oth_ref_ins_seq
#from dro_tools.matrices import is_matrix_orthonormal, is_matrix_orthogonal
from dro_tools.matrices import get_tx_matrix_type

def create_spa_dro_OBSOLETE(srcDataset, trgDataset, newDataset, params, description=''):
    """
    SEE: create_spa_dro() in dro_tools.create_dro.py
    
    Create a DICOM Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    description : str, optional ('' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    Returns
    -------
    dro : Pydicom Object
        The DICOM Registration Object.
    
    Note
    ----
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    """
    Decide whether to reference all SOPInstanceUIDs in srcDicoms and trgDicoms
    within ReferencedInstanceSequence (not required), or just one (for 
    compactness) to record the SeriesInstanceUIDs of the two series.
    """
    #refAllSOPs = True
    refAllSOPs = False
    
    sampleDroDir = params.cfgDict['sampleDroDir']
    p2c = params.cfgDict['p2c']
    
    # The tuple of the registration transform parameters:
    txParams = newDataset.txParams
    
    srcDicoms = srcDataset.dicoms
    trgDicoms = trgDataset.dicoms
    
    #cwd = os.getcwd()

    #droDir = os.path.join(cwd, 'sample_DROs')
    
    droFname = 'sample_spatial_dro.dcm'
    #droFpath = os.path.join(droDir, droFname)
    droFpath = os.path.join(sampleDroDir, droFname)
        
    #""" Check if the sample DRO has already been downloaded: """   
    #if os.path.isfile(NewSpaDroFpath):
    #    print(f'Sample DRO already downloaded to:\n {droDir}\n')
    #else:
    #    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    #    print(f'Downloading DROs from {url}...\n')
    #    DownloadSampleDros()
    
    if p2c:
        if refAllSOPs:
            print('Note: All SOPInstanceUIDs of both series will be referenced',
                  'in ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
        else:
            print('Note: Only the first SOPInstanceUID of both series will be',
                  'referenced in',
                  'ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
    
    # Start timing:
    times = []
    times.append(time.time())
    
    """
    srcDicoms = import_dicoms(srcDicomDir)
    trgDicoms = import_dicoms(trgDicomDir)
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n*Took {dtime} s to import the source and target DICOMs.')
    """
        
    # Add the row vector [0, 0, 0, 1] to txParams to complete the Frame of
    # Reference Transformation Matrix:
    # (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
    #txParams = [str(item) for item in txParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
    txMatrix = [str(txParams[0]), str(txParams[1]), str(txParams[2]), '0.0',
                str(txParams[3]), str(txParams[4]), str(txParams[5]), '0.0',
                str(txParams[6]), str(txParams[7]), str(txParams[8]), '0.0',
                str(txParams[9]), str(txParams[10]), str(txParams[11]), '1.0']
    
    # Read in the DRO template:
    dro = dcmread(droFpath)
    
    #currentDate = time.strftime("%Y%m%d", time.gmtime())
    #currentTime = time.strftime("%H%M%S", time.gmtime())
    ##currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    timeNow = datetime.now()
    currentDate = timeNow.strftime('%Y%m%d')
    currentTime = timeNow.strftime('%H%M%S.%f')
    
    dro.InstanceCreationDate = currentDate
    dro.InstanceCreationTime = currentTime
    
    # Generate a new SOPInstanceUID:
    dro.SOPInstanceUID = generate_uid()
    dro.file_meta.MediaStorageSOPInstanceUID = dro.SOPInstanceUID
    
    
    dro.ContentDate = trgDicoms[0].ContentDate
    dro.ContentTime = trgDicoms[0].ContentTime
    dro.StudyDate = trgDicoms[0].StudyDate
    try:
        dro.SeriesDate = trgDicoms[0].SeriesDate
    except AttributeError:
        pass
    #dro.ContentDate = trgDicoms[0].ContentDate
    dro.ContentDate = currentDate
    dro.StudyTime = trgDicoms[0].StudyTime
    try:
        dro.SeriesTime = trgDicoms[0].SeriesTime
    except AttributeError:
        pass
    #dro.ContentTime = trgDicoms[0].ContentTime
    dro.ContentTime = currentTime
    #dro.Manufacturer = trgDicoms[0].Manufacturer
    dro.Manufacturer = 'NCITA'
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    try:
        dro.SeriesDescription = trgDicoms[0].SeriesDescription
    except AttributeError:
        pass
    #dro.ManufacturerModelName = trgDicoms[0].ManufacturerModelName
    dro.ManufacturerModelName = ''
    dro.PatientName = trgDicoms[0].PatientName
    dro.PatientID = trgDicoms[0].PatientID
    dro.PatientBirthDate = trgDicoms[0].PatientBirthDate
    dro.PatientSex = trgDicoms[0].PatientSex
    try:
        dro.PatientAge = trgDicoms[0].PatientAge
    except AttributeError:
        pass
    dro.SoftwareVersions = ''
    dro.StudyInstanceUID = trgDicoms[0].StudyInstanceUID
    # Generate a new UID for SeriesInstanceUID:
    dro.SeriesInstanceUID = generate_uid()
    dro.StudyID = trgDicoms[0].StudyID
    dro.SeriesNumber = trgDicoms[0].SeriesNumber
    dro.InstanceNumber = trgDicoms[0].InstanceNumber
    dro.FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
    dro.PositionReferenceIndicator = trgDicoms[0].PositionReferenceIndicator
    """ Keep ContentLabel as 'REGISTRATION'. """
    """ Keep ContentCreatorName as ''. """
    
    # Modify the RegistrationSequence for the fixed domain:
    dro.RegistrationSequence[0]\
       .FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
    
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
    Consider adding FrameOfReferenceTransformationComment:
    dro.RegistrationSequence[0]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = 
    """
    
    # Modify the RegistrationSequence for the moving domain:
    dro.RegistrationSequence[1]\
       .FrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
    
    # Modify the FrameOfReferenceTransformationMatrix:
    dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrix = txMatrix # was txParams (06/05/21)
       
    # Modify the FrameOfReferenceTransformationMatrixType.  
    """
    Acceptable values:
     'RIGID' (= value in template DRO)
     'RIGID_SCALE' (Similarity transform)
     'AFFINE' 
     Note:  Need to verify that get_tx_matrix_type covers all 3 options.
    """
    dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(
           txParams, p2c
           )
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n   *Took {dtime} s to get the transform matrix type.')
        
    """ 
    Consider adding FrameOfReferenceTransformationComment:
    dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = ''
    """
     
    # Modify the ReferencedSeriesSequence:
    """
    Note:
        ReferencedSeriesSequence will have two items: The first for the fixed 
        series, and second for the moving series.
    """
    # Modify ReferencedSeriesSequence for the fixed series:
    dro = modify_ref_ser_seq(
        dro, seqNum=0, dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
    
    # Modify ReferencedSeriesSequence for the moving series:
    dro = modify_ref_ser_seq(
        dro, seqNum=1, dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
     
    # Modify the StudiesContainingOtherReferencedInstancesSequence:
    """
    Note:
        StudiesContainingOtherReferencedInstancesSequence will have two items: 
        The first for the fixed series, and second for the moving series.
    """
    # Modify StudiesContainingOtherReferencedInstancesSequence for the fixed 
    # series:
    dro = modify_stu_con_oth_ref_ins_seq(
        dro, seqNum=0, dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
    
    # Modify StudiesContainingOtherReferencedInstancesSequence for the moving 
    # series:
    dro = modify_stu_con_oth_ref_ins_seq(
        dro, seqNum=1, dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
        )
    
    """ > Consider adding UsedFiducialsSequence (optional). 
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
    
    Spatial Fiducials Series Module:
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
    """
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    if p2c:
        print('Deleting GroupLength tags.')
    
    """
    The SRO template contains 4 retired tags with name "GroupLength". 
    Delete them.
    """
    del dro[0x08, 0x00]
    del dro[0x10, 0x00]
    del dro[0x20, 0x00]
    del dro[0x70, 0x00]
    
    #times.append(time.time())
    #dtime = round(times[-1] - times[-2], 1)
    #if p2c:
    #    print(f'\n   *Took {dtime} s to delete all GroupLength tags.')
        
    dtime = round(times[-1] - times[0], 1)
    if p2c:
        print(f'\n   *Took a total of {dtime} s to create the DRO.')
    
    return dro

def create_spa_dro_from_tx_OBSOLETE(
        srcDataset, trgDataset, newDataset, params, description=''
        ):
    """
    SEE: create_spa_dro() in dro_tools.create_dro.py
    
    Create a DICOM Spatial Registration Object (DRO) from a SimpleITK 
    transform.  
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    description : str, optional ('' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    Returns
    -------
    dro : Pydicom Object
        The DICOM Registration Object.
    
    Note
    ----
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    # The registration and pre-registration transforms:
    regTx = newDataset.finalTx
    
    p2c = params.cfgDict['p2c']
    
    # The transform parameters:
    txParams = regTx.GetParameters()
    
    # Add the row vector [0, 0, 0, 1] to txParams to complete the Frame of
    # Reference Transformation Matrix 
    # (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
    #txParams = [str(item) for item in txParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
    txMatrix = [str(txParams[0]), str(txParams[1]), str(txParams[2]), '0.0',
                str(txParams[3]), str(txParams[4]), str(txParams[5]), '0.0',
                str(txParams[6]), str(txParams[7]), str(txParams[8]), '0.0',
                str(txParams[9]), str(txParams[10]), str(txParams[11]), '1.0']
    
    if p2c:
        print(f'txParams = {txParams}\n')
        print(f'txMatrix = {txMatrix}\n')
    
    dro = create_spa_dro(
        srcDataset, trgDataset, newDataset, params, description=''
        )
        
    return dro
