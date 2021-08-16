# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 12:46:13 2021

@author: ctorti
"""




"""
******************************************************************************
******************************************************************************
DICOM SPATIAL OR DEFORMABLE SPATIAL REGISTRATION OBJECT (SRO) RELATED FUNCTIONS
******************************************************************************
******************************************************************************
"""


def DownloadSampleDros():
    """
    08/07/21:
        See download_samples in dro_tools.download.py
        
    Download sample DICOM Registration Objects from:
    https://www.insight-journal.org/browse/publication/923
    """
    
    import os
    from pathlib import Path
    import time
    import requests
    import zipfile
    import tarfile
    import shutil
    from GeneralTools import DeleteFile, DeleteDirectory
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DRO')
    
    # Create the directory if it doesn't exist
    if not os.path.isdir(DroDir):
        #os.mkdir(ExportDir)
        Path(DroDir).mkdir(parents=True)
    
    # Start timing:
    times = []
    times.append(time.time())
    
    """ Download the zip file:"""
    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    filerequest = requests.get(url)
    
    ZipFname = r'Insight Journal_923_1.zip'
    ZipFpath = os.path.join(DroDir, ZipFname)

    with open(ZipFpath, 'wb') as file:
        file.write(filerequest.content)
        print(f"Downloaded '{ZipFname}' from {url}\n")

    if not file.closed:
        file.close()
        
    """ Extract the zip file. """
    ZipFile = zipfile.ZipFile(ZipFpath, mode='r')
    ZipFile.extractall(DroDir)
    print(f"Extracted '{ZipFname}' to {DroDir}\n")
    
    """ Uncompress the tar file: """
    TarFname = 'dicom-sro.tar.gz'
    TarFpath = os.path.join(DroDir, TarFname)

    with tarfile.open(TarFpath) as file:
        file.extractall(DroDir)
        print(f"Extracted '{TarFname}' to {DroDir}\n")

    if not file.closed:
        file.close()
    
    
    # The Spatial Registration Storage DRO:
    SpaDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337388794.677.609_0001_000001_137960834102c0.dcm'
    SpaDroFpath = os.path.join(DroDir, 'dicom-sro', 'rect-offset-sro-rigid', 
                               SpaDroFname)

    # The Deformable Spatial Registration Storage DRO:
    DefDroFname = '2.16.840.1.114362.1.6.0.6.13712.7616481285.337389158.629.1276_0001_000001_137960823902bf.dcm'
    DefDroFpath = os.path.join(DroDir, 'dicom-sro', 
                               'sphere-centered-sro-deformable', DefDroFname)

    NewSpaDroFname = 'sample_spatial_dro.dcm'
    NewSpaDroFpath = os.path.join(DroDir, NewSpaDroFname)
    
    NewDefDroFname = 'sample_deformable_dro.dcm'
    NewDefDroFpath = os.path.join(DroDir, NewDefDroFname)
    
    """ Copy the sample DROs to DroDir: """
    shutil.copy(SpaDroFpath, NewSpaDroFpath)
    shutil.copy(DefDroFpath, NewDefDroFpath)

    print(f"Copied sample DROs to {DroDir}\n")
    
    """ Delete unnecessary files and directories: """
    print('Tidying up..\n')
    
    fpath = os.path.join(DroDir, 'dicom-sro.tar.gz')
    DeleteFile(fpath)
    
    fpath = os.path.join(DroDir, 'Insight Journal_923_1.zip')
    DeleteFile(fpath)
    
    fpath = os.path.join(DroDir, 'Git_8118756.zip')
    DeleteFile(fpath)
    
    fpath = os.path.join(DroDir, 'ITKIOTransformDCMTK.pdf')
    DeleteFile(fpath)
    
    dirpath = os.path.join(DroDir, 'dicom-sro')       
    DeleteDirectory(dirpath)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    print(f'\n*Took {Dtime} s to download the sample DROs.')


def ModifyRefImSeq(Sequence, Dicoms, RefAllSOPs=False, LogToConsole=False):
    """
    08/07/21:
        See modify_ref_im_seq in dro_tools.general.py
      
    Modify a ReferencedImageSequence (or ReferencedInstanceSequence).  
    
    Inputs:
    ******
    
    Sequence : Sequence of a Pydicom object
        A ReferencedImageSequence to be modified.
    
    SeqNum : integer
        The (zero-indexed) integer of the sequence to modify.
    
    Dicoms : list of Pydicom objects
        The DICOMs that related to the sequence.
    
    RefAllSOPs : boolean (optional; False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Sequence : Sequence of a Pydicom object
        The modified ReferencedImageSequence.
    """
    
    import time
    #from copy import deepcopy
    
    # Start timing:
    times = []
    times.append(time.time())
       
    # The number of sequences:
    N0 = len(Sequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        N1 = len(Dicoms)
        
        if LogToConsole:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence,',
                  f'and {N1} (= no. of DICOMs) required sequences.')
    else:
        N1 = 1
        
        if LogToConsole:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence',
                  f'and {N1} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if N1 > N0:
        for i in range(N1 - N0):           
            Sequence.append(Sequence[-1])          
    else:
        for i in range(N0 - N1):
            Sequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to ReferencedImageSequence.')
        
    # Update the sequences:
    for i in range(N1):
        Sequence[i].ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Sequence[i].ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update ReferencedImageSequence.')
    
    return Sequence


def ModifyRefSerSeq(Dro, SeqNum, Dicoms, RefAllSOPs=False, LogToConsole=False):
    """
    08/07/21:
        See modify_ref_ser_seq in dro_tools.general.py
      
    Modify the ReferencedSeriesSequence in a DICOM Spatial or Deformable
    Spatial Registration Object (DRO) for a chosen sequence number.  
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Registration Object to be modified.
    
    SeqNum : integer
        The (zero-indexed) integer of the sequence to modify.
    
    Dicoms : list of Pydicom objects
        The DICOMs that related to the sequence.
    
    RefAllSOPs : boolean (optional; False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The modified DICOM Registration Object.
    """
    
    import time
    from copy import deepcopy
    
    # Start timing:
    times = []
    times.append(time.time())
    
    Dro.ReferencedSeriesSequence[SeqNum]\
       .SeriesInstanceUID = Dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[SeqNum]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(Dicoms)
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (= no. of DICOMs) required sequences.')
    else:
        ReqRIS = 1
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[SeqNum]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[SeqNum]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[SeqNum]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[SeqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[SeqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update',
              f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence.')
    
    return Dro


def ModifyStuConOthRefInsSeq(Dro, SeqNum, Dicoms, RefAllSOPs=False,
                             LogToConsole=False):
    """
    08/07/21:
        See modify_stu_con_oth_ref_ins_seq in dro_tools.general.py
      
    Modify the StudiesContainingOtherReferencedInstancesSequence in a DICOM 
    Spatial or Deformable Spatial Registration Object (DRO) for a chosen 
    sequence number.  
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Registration Object to be modified.
    
    SeqNum : integer
        The (zero-indexed) integer of the sequence to modify.
    
    Dicoms : list of Pydicom objects
        The DICOMs that related to the sequence.
    
    RefAllSOPs : boolean (optional; False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The modified DICOM Registration Object.
    """
    
    import time
    from copy import deepcopy
    
    # Start timing:
    times = []
    times.append(time.time())
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
       .StudyInstanceUID = Dicoms[0].StudyInstanceUID
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = Dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(Dicoms)
        
        if LogToConsole:
            msg = f'\nThere are {OrigRIS} sequences in '\
                  + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
                  + f'.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                  + f' and {ReqRIS} (= no. of DICOMs) required sequences.'
            print(msg)
    else:
        ReqRIS = 1
        
        if LogToConsole:
            msg = f'\nThere are {OrigRIS} sequences in '\
                  + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
                  + f'.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                  + f' and {ReqRIS} (the first DICOM) required sequence.'
            print(msg)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        msg = f'\n   *Took {Dtime} s to add sequences to '\
              + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
              + f'.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        msg = f'\n   *Took {Dtime} s to update '\
              + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
              + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
    
    return Dro


def CreateSpaDro(SrcDicomDir, TrgDicomDir, TxParams, Description='', 
                 LogToConsole=False):
    """
    08/07/21:
        See create_spa_dro in dro_tools.tx_to_spa_dro.py
     
    Create a DICOM Spatial Registration Object (DRO).  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    TxParams : list of floats
        List of floats representing the non-deformable transformation 
        parameters from the SimpleElastix image filter or SimpleITK transform.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    import os
    import time
    from datetime import datetime
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #from copy import deepcopy
    from DicomTools import ImportDicoms
    from GeneralTools import GetTxMatrixType
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DROs')
    
    DroFname = 'sample_spatial_dro.dcm'
    DroFpath = os.path.join(DroDir, DroFname)
        
    
    #""" Check if the sample DRO has already been downloaded: """   
    #if os.path.isfile(NewSpaDroFpath):
    #    print(f'Sample DRO already downloaded to:\n {DroDir}\n')
    #else:
    #    url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
    #    print(f'Downloading DROs from {url}...\n')
    #    DownloadSampleDros()
    
    
    """ Decide whether to reference all SOPInstanceUIDs in SrcDicoms and 
    TrgDicoms within ReferencedInstanceSequence (not required), or just one
    (for compactness), for the sake of maintaining the SeriesInstanceUIDs of 
    the two series. """
    #RefAllSOPs = True
    RefAllSOPs = False
    
    if LogToConsole:
        if RefAllSOPs:
            print('Note: All SOPInstanceUIDs of both series will be referenced',
                  'in ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
        else:
            print('Note: Only the first SOPInstanceUID of both series will be',
                  'referenced in',
                  'ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
    
    # Start timing:
    times = []
    times.append(time.time())
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
        
    
    """ Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
    Reference Transformation Matrix 
    (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
    """
    #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
    TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
                str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
                str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
                str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    
    """ Read in the DRO template. """
    Dro = dcmread(DroFpath)
    
    #CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    #CurrentTime = time.strftime("%H%M%S", time.gmtime())
    ##CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    TimeNow = datetime.now()
    CurrentDate = TimeNow.strftime('%Y%m%d')
    CurrentTime = TimeNow.strftime('%H%M%S.%f')
    
    Dro.InstanceCreationDate = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.StudyDate = TrgDicoms[0].StudyDate
    try:
        Dro.SeriesDate = TrgDicoms[0].SeriesDate
    except AttributeError:
        pass
    #Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentDate = CurrentDate
    Dro.StudyTime = TrgDicoms[0].StudyTime
    try:
        Dro.SeriesTime = TrgDicoms[0].SeriesTime
    except AttributeError:
        pass
    #Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.ContentTime = CurrentTime
    #Dro.Manufacturer = TrgDicoms[0].Manufacturer
    Dro.Manufacturer = 'NCITA'
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    try:
        Dro.SeriesDescription = TrgDicoms[0].SeriesDescription
    except AttributeError:
        pass
    #Dro.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
    Dro.ManufacturerModelName = ''
    Dro.PatientName = TrgDicoms[0].PatientName
    Dro.PatientID = TrgDicoms[0].PatientID
    Dro.PatientBirthDate = TrgDicoms[0].PatientBirthDate
    Dro.PatientSex = TrgDicoms[0].PatientSex
    try:
        Dro.PatientAge = TrgDicoms[0].PatientAge
    except AttributeError:
        pass
    Dro.SoftwareVersions = ''
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
       .FrameOfReferenceTransformationMatrix = TxMatrix # was TxParams (06/05/21)
       
    """
    Modify the FrameOfReferenceTransformationMatrixType.  Acceptable values:
    'RIGID' (= value in template DRO)
    'RIGID_SCALE' (Similarity transform)
    'AFFINE' 
    
    Note:  Need to verify that GetTxMatrixType covers all 3 options.
    """
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .MatrixSequence[0]\
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(TxParams, 
                                                                   LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
        
    """ 
    Consider adding FrameOfReferenceTransformationComment. 
    
    Dro.RegistrationSequence[1]\
       .MatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = ''
    """
    
    
    
    """ 
    Modify the ReferencedSeriesSequence. 
    
    Note:
        ReferencedSeriesSequence will have two items: The first for the fixed 
        series, and second for the moving series.
    """
    # Modify ReferencedSeriesSequence for the fixed series:
    Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
                          RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    
    # Modify ReferencedSeriesSequence for the moving series:
    Dro = ModifyRefSerSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
                          RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    
    
    """ 
    Modify the StudiesContainingOtherReferencedInstancesSequence. 
    
    Note:
        StudiesContainingOtherReferencedInstancesSequence will have two items: 
        The first for the fixed series, and second for the moving series.
    """
    # Modify StudiesContainingOtherReferencedInstancesSequence for the fixed 
    # series:
    Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
                                   RefAllSOPs=RefAllSOPs, 
                                   LogToConsole=LogToConsole)
    
    # Modify StudiesContainingOtherReferencedInstancesSequence for the moving 
    # series:
    Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
                                   RefAllSOPs=RefAllSOPs, 
                                   LogToConsole=LogToConsole)
    
    
    """ > Consider adding UsedFiducialsSequence (optional). 
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
    
    Spatial Fiducials Series Module:
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
    """
    
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    if LogToConsole:
        print('Deleting GroupLength tags.')
    
    """
    The SRO template contains 4 retired tags with name "GroupLength". 
    Delete them.
    """
    
    del Dro[0x08, 0x00]
    del Dro[0x10, 0x00]
    del Dro[0x20, 0x00]
    del Dro[0x70, 0x00]
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if LogToConsole:
    #    print(f'\n   *Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
        
        
    return Dro


def CreateSpaDroFromSitkTx(SrcDicomDir, TrgDicomDir, SitkTx, Description='', 
                           LogToConsole=False):
    """
    08/07/21:
        See create_spa_dro_from_tx in dro_tools.tx_to_spa_dro.py
     
    Create a DICOM Spatial Registration Object (DRO) from a SimpleITK 
    transform.  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    SitkTx : SimpleITK Transform
        The SimpleITK transform containing the transform parameters.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    # The transform parameters:
    TxParams = SitkTx.GetParameters()
    
    """ Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
    Reference Transformation Matrix 
    (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
    """
    #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
    TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
                str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
                str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
                str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    if LogToConsole:
        print(f'TxParams = {TxParams}\n')
        print(f'TxMatrix = {TxMatrix}\n')
    
    Dro = CreateSpaDro(SrcDicomDir, TrgDicomDir, TxParams, Description, 
                       LogToConsole)
        
    return Dro


def CreateDefDro_OLD(SrcDicomDir, TrgDicomDir, GridOrig, GridDir, GridDims, GridRes,
                 VectGridData, TxParams=None, Description='', LogToConsole=False):
    """
    Create a DICOM Deformable Spatial Registration Object (DRO).  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    GridOrig : list of float strings
        List of float strings representing the origin of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    
    GridDir : list of float strings
        List of float strings representing the direction of the BSpline grid 
        used to deform the moving (Source) image to the fixed (Target) image.
    
    GridDims : list of integers
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    
    GridRes : list of floats
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image.
    
    VectGridData : list of float strings
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image.
    
    TxParams : list of floats or None (optional; None by default)
        List of floats representing the non-deformable transformation 
        parameters from the SimpleElastix image filter or SimpleITK transform
        used prior to the bspline transformation used to register the moving 
        (Source) image to the fixed (Target) image; and will be stored in 
        PreDeformationMatrixRegistrationSequence.  If None, the identity 
        matrix will be used.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    from copy import deepcopy
    from DicomTools import ImportDicoms
    from GeneralTools import GetTxMatrixType
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DRO')
    
    #NewSpaDroFname = 'sample_spatial_dro.dcm'
    #NewSpaDroFpath = os.path.join(DroDir, NewSpaDroFname)
    
    NewDefDroFname = 'sample_deformable_dro.dcm'
    NewDefDroFpath = os.path.join(DroDir, NewDefDroFname)
        
    
    """ Check if the sample DRO has already been downloaded: """   
    if os.path.isfile(NewDefDroFpath):
        print(f'Sample DRO already downloaded to:\n {DroDir}\n')
    else:
        url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
        print(f'Downloading DROs from {url}...\n')
        DownloadSampleDros()
    
    
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
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
    
    
    if TxParams == None:
        TxMatrix = ['1.0', '0.0', '0.0', '0.0',
                    '0.0', '1.0', '0.0', '0.0',
                    '0.0', '0.0', '1.0', '0.0',
                    '0.0', '0.0', '0.0', '1.0']
    else:
        """ Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
        Reference Transformation Matrix 
        (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
        """
        #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
        TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
                    str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
                    str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
                    str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    
    """ Select an appropriate DRO template. """
    #if VectGridData == None:
    #    """ The registration was non-elastic. Use a DRO with Spatial 
    #    Registration Storage SOP Class UID. """
    #    Dro = dcmread(NewSpaDroFpath)
    #else:
    #    """ The registration was elastic. Use a DRO with Deformable Spatial 
    #    Registration Storage SOP Class UID. """
    #    Dro = dcmread(NewDefDroFpath)
    Dro = dcmread(NewDefDroFpath)
    
    
    CurrentDate = time.strftime("%Y%m%d", time.gmtime())
    CurrentTime = time.strftime("%H%M%S", time.gmtime())
    #CurrentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    
    Dro.InstanceCreationDate = CurrentDate
    Dro.InstanceCreationTime = CurrentTime
    
    # Generate a new SOPInstanceUID:
    Dro.SOPInstanceUID = generate_uid()
    Dro.file_meta.MediaStorageSOPInstanceUID = Dro.SOPInstanceUID
    
    
    Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.StudyDate = TrgDicoms[0].StudyDate
    try:
        Dro.SeriesDate = TrgDicoms[0].SeriesDate
    except AttributeError:
        pass
    #Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentDate = CurrentDate
    Dro.StudyTime = TrgDicoms[0].StudyTime
    try:
        Dro.SeriesTime = TrgDicoms[0].SeriesTime
    except AttributeError:
        pass
    #Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.ContentTime = CurrentTime
    #Dro.Manufacturer = TrgDicoms[0].Manufacturer
    Dro.Manufacturer = 'NCITA'
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    try:
        Dro.SeriesDescription = TrgDicoms[0].SeriesDescription
    except AttributeError:
        pass
    #Dro.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
    Dro.ManufacturerModelName = ''
    Dro.PatientName = TrgDicoms[0].PatientName
    Dro.PatientID = TrgDicoms[0].PatientID
    Dro.PatientBirthDate = TrgDicoms[0].PatientBirthDate
    Dro.PatientSex = TrgDicoms[0].PatientSex
    try:
        Dro.PatientAge = TrgDicoms[0].PatientAge
    except AttributeError:
        pass
    Dro.SoftwareVersions = ''
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
    
    
    
    
    
    
    """ Modify the DeformableRegistrationSequence. """
    if LogToConsole:
        print('\nModify the DeformableRegistrationSequence.')
    
    """ > Modify the SourceFrameOfReferenceUID.
    Note:
    
    The DICOM standard defines registration from Source to Registered RCS 
    (Referenced Coordinate System).
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_A.39.html
    """
    Dro.DeformableRegistrationSequence[0]\
       .SourceFrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    
    
    
    """ > Modify PreDeformationMatrixRegistrationSequence. """
    Dro.DeformableRegistrationSequence[0]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxMatrix
    
    """ > Modify the FrameOfReferenceTransformationMatrixType.  Acceptable 
    values:
    'RIGID' (= value in template DRO)
    'RIGID_SCALE' (Similarity transform)
    'AFFINE' 
    
    Note:  Need to verify that GetTxMatrixType covers all 3 options.
    """
    
    print(f'TxParams = {TxParams}')
    
    Dro.DeformableRegistrationSequence[0]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(TxParams, 
                                                                   LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
    
    
    """ > Consider adding details for the FrameOfReferenceTransformationComment
    and RegistrationTypeCodeSequence (optional). 
    Dro.DeformableRegistrationSequence[0]\
       .FrameOfReferenceTransformationComment = ''
    
    Dro.DeformableRegistrationSequence[0]\
       .RegistrationTypeCodeSequence = ''
    """
    
    
    """ > Leave PostDeformationMatrixRegistrationSequence unchanged. """
    
    
    
    """ > Modify the DeformableRegistrationGridSequence. """
    Dro.DeformableRegistrationSequence[0]\
       .DeformableRegistrationGridSequence[0]\
       .ImagePositionPatient = GridOrig
    
    Dro.DeformableRegistrationSequence[0]\
       .DeformableRegistrationGridSequence[0]\
       .ImageOrientationPatient = GridDir
    
    Dro.DeformableRegistrationSequence[0]\
       .DeformableRegistrationGridSequence[0]\
       .GridDimensions = GridDims
    
    Dro.DeformableRegistrationSequence[0]\
       .DeformableRegistrationGridSequence[0]\
       .GridResolution = GridRes
    
    Dro.DeformableRegistrationSequence[0]\
       .DeformableRegistrationGridSequence[0]\
       .VectorGridData = VectGridData
    
    
    
    """ 
    Modify the ReferencedSeriesSequence. 
    
    Note:
        ReferencedSeriesSequence will have two items: The first for the fixed 
        series, and second for the moving series.
        
        Each sequence in ReferencedSeriesSequence should contain a
        ReferencedInstanceSequence (containing a ReferencedSOPClassUID and a
        ReferencedSOPInstanceUID) and a SeriesInstanceUID, e.g.:
        
        (0008, 1115)  Referenced Series Sequence   2 item(s) ---- 
        (0008, 114a)  Referenced Instance Sequence   1 item(s) ---- 
           (0008, 1150) Referenced SOP Class UID            UI: CT Image Storage
           (0008, 1155) Referenced SOP Instance UID         UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.19039.1372783046.430261
           ---------
        (0020, 000e) Series Instance UID                 UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.19039.1372783046.430164
        ---------
        (0008, 114a)  Referenced Instance Sequence   1 item(s) ---- 
           (0008, 1150) Referenced SOP Class UID            UI: CT Image Storage
           (0008, 1155) Referenced SOP Instance UID         UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.22085.1372783046.709969
           ---------
        (0020, 000e) Series Instance UID                 UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.22085.1372783046.709933
        ---------
        
        https://dicom.innolitics.com/ciods/deformable-spatial-registration/common-instance-reference/00081115
        
        But the sample DRO doesn't have a ReferencedInstanceSequence within
        each ReferencedSeriesSequence, but instead contains it's children
        (ReferencedSOPClassUID and ReferencedSOPInstanceUID), followed by a
        SeriesInstanceUID, followed by another pair of children of
        ReferencedInstanceSequence, e.g.:
        
        (0008, 1115)  Referenced Series Sequence   3 item(s) ---- 
        (0008, 1150) Referenced SOP Class UID            UI: CT Image Storage
        (0008, 1155) Referenced SOP Instance UID         UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.25113.1372783047.992809
        ---------
        (0020, 000e) Series Instance UID                 UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.25113.1372783047.992704
        ---------
        (0008, 1150) Referenced SOP Class UID            UI: CT Image Storage
        (0008, 1155) Referenced SOP Instance UID         UI: 1.2.826.0.1.3680043.8.274.1.1.8323328.19039.1372783046.430269
        ---------
        
        So copy StudiesContainingOtherReferencedInstancesSequence (which does
        follow the DICOM standard) to ReferencedSeriesSequence, then modify
        the sequences.
    """
    
    print('\nBefore:\n')
    for seq in Dro.ReferencedSeriesSequence:
        print(seq, '\n')
    
    #Dro.ReferencedSeriesSequence = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence)
    Dro.ReferencedSeriesSequence[0] = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[0])
    Dro.ReferencedSeriesSequence[1] = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[1])
    
    
    print('\nAfter copying StudiesContainingOtherReferencedInstancesSequence:\n')
    for seq in Dro.ReferencedSeriesSequence:
        print(seq, '\n')
    
    """ Delete StudyInstanceUID from each sequence in ReferencedSeriesSequence
    """
    for i in range(len(Dro.ReferencedSeriesSequence)):
        del Dro.ReferencedSeriesSequence[i][0x20, 0x0d]
    
    print('\nAfter deleting StudyInstanceUID:\n')
    for seq in Dro.ReferencedSeriesSequence:
        print(seq, '\n')
        
    # Modify ReferencedSeriesSequence for the fixed series:
    Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
                          RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    
    # Modify ReferencedSeriesSequence for the moving series:
    Dro = ModifyRefSerSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
                          RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    
    
    """ 
    Modify the StudiesContainingOtherReferencedInstancesSequence. 
    
    Note:
        StudiesContainingOtherReferencedInstancesSequence will have two items: 
        The first for the fixed series, and second for the moving series.
    """
    # Modify StudiesContainingOtherReferencedInstancesSequence for the fixed 
    # series:
    Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
                                   RefAllSOPs=RefAllSOPs, 
                                   LogToConsole=LogToConsole)
    
    # Modify StudiesContainingOtherReferencedInstancesSequence for the moving 
    # series:
    Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
                                   RefAllSOPs=RefAllSOPs, 
                                   LogToConsole=LogToConsole)
    
    
    
    
    """ > Consider adding UsedFiducialsSequence (optional). 
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
    
    Spatial Fiducials Series Module:
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
    """
    
    
    
    """
    https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
    """
    
    if LogToConsole:
        print('Deleting GroupLength tags.')
    
    """
    The SRO template contains 5 retired tags with name "GroupLength". 
    Delete them.
    """
    
    del Dro[0x08, 0x00]
    del Dro[0x10, 0x00]
    del Dro[0x20, 0x00]
    del Dro[0x64, 0x00]
    del Dro[0x70, 0x00]
    
    #times.append(time.time())
    #Dtime = round(times[-1] - times[-2], 1)
    #if LogToConsole:
    #    print(f'\n   *Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
        
        
    return Dro



def CreateDefDro(SrcDicomDir, TrgDicomDir, 
                 GridOrig, GridDir, GridDims, GridRes, VectGridData,
                 TxMatrix=None, Description='', LogToConsole=False):
    """
    08/07/21:
        See create_def_dro in dro_tools.tx_to_def_dro.py
     
    Create a DICOM Deformable Spatial Registration Object (DRO).  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    GridOrig : list of float strings
        List of float strings representing the origin of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    
    GridDir : list of float strings
        List of float strings representing the direction of the BSpline grid 
        used to deform the moving (Source) image to the fixed (Target) image.
    
    GridDims : list of integers
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image.
    
    GridRes : list of floats
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image.
    
    VectGridData : list of float strings
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image.
    
    TxMatrix : list of strings or None (optional; None by default)
        List of strings representing the non-deformable transformation 
        parameters to be applied prior to the bspline registration. TxMatrix
        will be stored in PreDeformationMatrixRegistrationSequence.  If None, 
        the identity matrix will be used.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    """
    
    import os
    import time
    import datetime
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    #from copy import deepcopy
    from DicomTools import ImportDicoms
    #import importlib
    #import GeneralTools
    #importlib.reload(GeneralTools)
    from GeneralTools import GetTxMatrixType, ReduceListOfStringFloatsTo16
    
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
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
    
    
    # The following was moved to CreateDefDroFromBsplineTx(): (04/06/21) 
    #if TxParams == None:
    #    TxMatrix = ['1.0', '0.0', '0.0', '0.0',
    #                '0.0', '1.0', '0.0', '0.0',
    #                '0.0', '0.0', '1.0', '0.0',
    #                '0.0', '0.0', '0.0', '1.0']
    #else:
    #    """ Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
    #    Reference Transformation Matrix 
    #    (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
    #    """
    #    #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
    #    TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
    #                str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
    #                str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
    #                str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    
    """ Read in the DRO template. """
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
    
    
    
    """ 
    Modify the ReferencedSeriesSequence. 
    
    Note:
        ReferencedSeriesSequence will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=SrcDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    
    Seq = ModifyRefImSeq(Dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence,
                         Dicoms=SrcDicoms, RefAllSOPs=RefAllSOPs, 
                         LogToConsole=LogToConsole)
    
    Dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence = Seq
    
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
    
    
    """ 
    Modify the StudiesContainingOtherReferencedInstancesSequence. 
    
    Note:
        StudiesContainingOtherReferencedInstancesSequence will reference the
        fixed (target) series.
    """
    #Dro = ModifyStuConOthRefInsSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
    #                               RefAllSOPs=RefAllSOPs, 
    #                               LogToConsole=LogToConsole)
    
    Seq = ModifyRefImSeq(Dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value,
                         Dicoms=TrgDicoms, RefAllSOPs=RefAllSOPs, 
                         LogToConsole=LogToConsole)
    
    Dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value = Seq
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[0]\
       .StudyInstanceUID = TrgDicoms[0].StudyInstanceUID
    
    
    
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
    if LogToConsole:
        print('\nModify the DeformableRegistrationSequence.')
    
    #""" Remove the first sequence in DeformableRegistrationSequence (which
    #references the fixed image series). """
    #Dro.DeformableRegistrationSequence.pop(0)
    
    
    """ 
    Modify the first ReferencedImageSequence. 
    
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=0, Dicoms=TrgDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    Seq = ModifyRefImSeq(Dro.DeformableRegistrationSequence[0].ReferencedImageSequence,
                         Dicoms=TrgDicoms, RefAllSOPs=True, 
                         LogToConsole=LogToConsole)
    
    Dro.DeformableRegistrationSequence[0].ReferencedImageSequence = Seq
    
    
    """ > Modify the SourceFrameOfReferenceUID.
    Note:
    
    The DICOM standard defines registration from Source to Registered RCS 
    (Referenced Coordinate System).
    
    http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_A.39.html
    """
    Dro.DeformableRegistrationSequence[0]\
       .SourceFrameOfReferenceUID = TrgDicoms[0].FrameOfReferenceUID
    
    
    """ Delete MatrixRegistrationSequence. """
    del Dro[0x64, 0x02][0][0x70, 0x0309]
    
    
    
    
    
    
    """ 
    Modify the second ReferencedImageSequence. 
    
    Note:
        ReferencedImageSequence[0] will reference the fixed (target) series and
        ReferencedImageSequence[1] will reference the moving (source) series.
    """
    #Dro = ModifyRefSerSeq(Dro, SeqNum=1, Dicoms=SrcDicoms, 
    #                      RefAllSOPs=RefAllSOPs, LogToConsole=LogToConsole)
    Seq = ModifyRefImSeq(Dro.DeformableRegistrationSequence[1].ReferencedImageSequence,
                         Dicoms=SrcDicoms, RefAllSOPs=True, 
                         LogToConsole=LogToConsole)
    
    Dro.DeformableRegistrationSequence[1].ReferencedImageSequence = Seq
    
    
    """ > Modify the SourceFrameOfReferenceUID. """
    Dro.DeformableRegistrationSequence[1]\
       .SourceFrameOfReferenceUID = SrcDicoms[0].FrameOfReferenceUID
    
    
    """ > Modify the DeformableRegistrationGridSequence. """
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImagePositionPatient = ReduceListOfStringFloatsTo16(GridOrig)
    
    if len(GridDir) > 6:
        # Ignore the direction cosines along the z-direction:
        GridDir = GridDir[0:6]
        
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .ImageOrientationPatient = ReduceListOfStringFloatsTo16(GridDir)
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridDimensions = GridDims
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .GridResolution = GridRes
    
    Dro.DeformableRegistrationSequence[1]\
       .DeformableRegistrationGridSequence[0]\
       .VectorGridData = VectGridData
    
    
    
    """ > Modify PreDeformationMatrixRegistrationSequence. """
    Dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrix = TxMatrix
    
    """ > Modify the FrameOfReferenceTransformationMatrixType.  Acceptable 
    values:
    'RIGID' (= value in template DRO)
    'RIGID_SCALE' (Similarity transform)
    'AFFINE' 
    
    Note:  Need to verify that GetTxMatrixType covers all 3 options.
    """
    
    #if TxParams != None: # 04/06/21
    #    print(f'TxParams = {TxParams}')
    #    
    #    Dro.DeformableRegistrationSequence[1]\
    #       .PreDeformationMatrixRegistrationSequence[0]\
    #       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(TxParams, 
    #                                                                   LogToConsole)
    #
    #    times.append(time.time())
    #    Dtime = round(times[-1] - times[-2], 3)
    #    if LogToConsole:
    #        print(f'\n   *Took {Dtime} s to get the transform matrix type.')
    
    # 04/06/21:
    Dro.DeformableRegistrationSequence[1]\
       .PreDeformationMatrixRegistrationSequence[0]\
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(TxMatrix, 
                                                                   LogToConsole)
           
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


def CreateDefDroFromBsplineTx_OLD(SrcDicomDir, TrgDicomDir, BsplineTx, 
                              TxParams=None, Description='', 
                              LogToConsole=False):
    """
    Create a DICOM Deformable Spatial Registration Object (DRO) from a 
    SimpleITK BSpline transformation.  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    BsplineTx : SimpleITK Transform
        The final transform following a deformable image registration with a
        BSpline Transform.
    
    TxParams : list of floats or None (optional; None by default)
        List of floats representing the non-deformable transformation 
        parameters from the SimpleElastix image filter or SimpleITK transform
        used prior to the BSpline transformation used to register the moving 
        (Source) image to the fixed (Target) image; and will be stored in 
        PreDeformationMatrixRegistrationSequence.  If None, the identity 
        matrix will be used.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    VectorGridDaa has VR = OF:
    
    A stream of 32-bit IEEE 754:1985 floating point words. OF is a VR that 
    requires byte swapping within each 32-bit word when changing byte ordering 
    (see Section 7.3).
    
    http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html
    """
    
    import SimpleITK as sitk
    import numpy as np
    
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
    
    Dro = CreateDefDro(SrcDicomDir, TrgDicomDir, GridOrig, GridDir, GridDims, 
                       GridRes, VectGridData, TxParams, Description, 
                       LogToConsole)
    
    return Dro


def CreateDefDroFromBsplineTx(SrcDicomDir, TrgDicomDir, BsplineTx, 
                              PreRegTx=None, Description='', 
                              LogToConsole=False):
    """
    08/07/21:
        See create_def_dro_from_tx in dro_tools.tx_to_def_dro.py
     
    Create a DICOM Deformable Spatial Registration Object (DRO) from a 
    SimpleITK BSpline transformation.  
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    BsplineTx : SimpleITK Transform
        The final transform following a deformable image registration with a
        BSpline Transform.
    
    PreRegTx : SimpleITK Transform or None (optional; None by default)
        The transform used prior to registration, and to be stored in 
        PreDeformationMatrixRegistrationSequence.  If None, the identity 
        matrix will be stored in PreDeformationMatrixRegistrationSequence.
    
    Description : string (optional; '' by default)
        Description to be added to ContentDescription, and used in the file
        name of the exported DRO.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    Dro : Pydicom object
        The DICOM Registration Object.
    
    
    Notes:
    *****
    
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
    
    import SimpleITK as sitk
    import numpy as np
    
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
    
    """ Get the transformation matrix (TxMatrix) from the pre-registration
    transform parameters: """
    if PreRegTx == None:
        #PreRegTxParams = None
        
        TxMatrix = ['1.0', '0.0', '0.0', '0.0',
                    '0.0', '1.0', '0.0', '0.0',
                    '0.0', '0.0', '1.0', '0.0',
                    '0.0', '0.0', '0.0', '1.0']
    else:
        TxParams = PreRegTx.GetParameters()
        
        """ Add the row vector [0, 0, 0, 1] to TxParams to complete the Frame of
        Reference Transformation Matrix 
        (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
        """
        #TxParams = [str(item) for item in TxParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
        TxMatrix = [str(TxParams[0]), str(TxParams[1]), str(TxParams[2]), '0.0',
                    str(TxParams[3]), str(TxParams[4]), str(TxParams[5]), '0.0',
                    str(TxParams[6]), str(TxParams[7]), str(TxParams[8]), '0.0',
                    str(TxParams[9]), str(TxParams[10]), str(TxParams[11]), '1.0']
    
    
    Dro = CreateDefDro(SrcDicomDir, TrgDicomDir, GridOrig, GridDir, GridDims, 
                       GridRes, VectGridData, TxMatrix, Description, 
                       LogToConsole)
    
    return Dro




"""
******************************************************************************
******************************************************************************
PARSE A DICOM SPATIAL OR DEFORMABLE SPATIAL REGISTRATION OBJECT (SRO) FROM 
DISK TO OBTAIN THE REQUIRED PARAMETERS TO CREATE A SimpleITK TRANSFORM
******************************************************************************
******************************************************************************
"""

def CreateSitkTxFromTxMatrix_INCOMPLETE(TxMatrix, Transform):
    """
    Create a SimpleITK transform from a list of elements that make up the
    transform matrix.
    
    Inputs:
    ******
    
    TxMatrix : list of float strings
        List of float strings representing the a transformation (e.g. as parsed 
        from FrameOfReferenceTransformationMatrix of a DRO).
    
    Transform : string
        The transform to create. Acceptable values include:
        - 'rigid' or 'RIGID'
        - 'rigid_scale' or 'RIGID_SCALE'
        - 'affine' or 'AFFINE'
    
        
    Outputs:
    *******
    
    SitkTx : SimpleITK Transform
        The SimpleITK Transform.
    
    
    Note:
    ****
    
    TxMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. Hence
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    """
    
    import SimpleITK as sitk
            
    if Transform.lower() == 'rigid':
        SitkTx = sitk.Euler3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:6]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
        
    
    elif Transform.lower() == 'rigid_scale':
        SitkTx = sitk.Similarity3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
        
    elif Transform.lower() == 'affine':
        SitkTx = sitk.AffineTransform(3)
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        SitkTx.SetParameters(TxParams) # 07/05/21
        
    else:
        msg = f"Transform = {Transform}. Only 'rigid', 'rigid_scale' and "\
              + "'affine' transforms are supported at present."
        
        raise Exception(msg)
    
    return SitkTx


def CreateSitkTxFromSpaDro(Dro, LogToConsole=False):
    """
    08/07/21:
        See create_tx_from_spa_dro in dro_tools.dro_to_tx.py
    
    Create a SimpleITK transformation based on parameters parsed from
    a DICOM Spatial Registration Object (DRO).  
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Spatial Registration Object.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    SitkTx : SimpleITK Transform
        A transform based on the parameters stored in Dro.
    
    
    Notes:
    *****
    
    It will be assumed that RegistrationSequence has two items: the first for
    the Fixed scan and thev second for the Moving scan.  While the 
    FrameOfReferenceMatrix in the second item is necessary, the first provides
    a convenient means of arciving the FrameOfReferenceUID of the Fixed scan.
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    import SimpleITK as sitk
    #import numpy as np
    
    MatrixType = Dro.RegistrationSequence[1]\
                    .MatrixRegistrationSequence[0]\
                    .MatrixSequence[0]\
                    .FrameOfReferenceTransformationMatrixType
    
    Matrix = [float(s) for s in Dro.RegistrationSequence[1]\
                                   .MatrixRegistrationSequence[0]\
                                   .MatrixSequence[0]\
                                   .FrameOfReferenceTransformationMatrix]
    
    if LogToConsole:
        print(f'MatrixType from DRO = {MatrixType}')
        print(f'Matrix from DRO = {Matrix}')
    
    if MatrixType == 'RIGID':
        SitkTx = sitk.Euler3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:6]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
        
    
    elif MatrixType == 'RIGID_SCALE':
        SitkTx = sitk.Similarity3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
        
    elif MatrixType == 'AFFINE':
        SitkTx = sitk.AffineTransform(3)
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        Matrix = [Matrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        SitkTx.SetParameters(Matrix) # 07/05/21
        
    else:
        msg = f"MatrixType = {MatrixType} is not a recognised value. Only "\
              + "'RIGID', 'RIGID_SCALE' and 'AFFINE' are acceptable values."
        
        raise Exception(msg)
    
    if LogToConsole:
        print(f'Parameters in new SitkTx = {SitkTx.GetParameters()}\n')
    
    return SitkTx



def CreatePreDefSitkTxFromDefDro(Dro, LogToConsole=False):
    """
    08/07/21:
        See create_pre_tx_from_def_dro in dro_tools.dro_to_tx.py
    
    Create a SimpleITK transformation based on parameters parsed from
    PreDeformationMatrixRegistrationSequence in a DICOM Deformable Spatial 
    Registration Object (DRO).  
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Deformable Registration Object.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    PreDefTx : SimpleITK Transform
        A transform based on the parameters stored in Dro.
    
    
    Notes:
    *****
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    Even thought the FrameOfReferenceTransformationMatrixType may be 'RIGID',
    'RIGID_SCALE' or 'AFFINE', create an affine transform for simplicity.
    
    TxMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. Hence
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    """
    
    import SimpleITK as sitk
    
    TxType = Dro.DeformableRegistrationSequence[1]\
                .PreDeformationMatrixRegistrationSequence[0]\
                .FrameOfReferenceTransformationMatrixType
    
    TxMatrix = Dro.DeformableRegistrationSequence[1]\
                  .PreDeformationMatrixRegistrationSequence[0]\
                  .FrameOfReferenceTransformationMatrix
    
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]]
    
    if LogToConsole:
        print(f'TxType = {TxType}\n')
        print(f'TxMatrix = {TxMatrix}\n')
        print(f'TxParams = {TxParams}\n')
    
    #if TxType == 'RIGID':
    #    PreDefTx = sitk.Euler3DTransform()
    #elif TxType == 'RIGID_SCALE':
    #    PreDefTx = sitk.Similarity3DTransform()
    #elif TxType == 'AFFINE':
    #    PreDefTx = sitk.AffineTransform(3)
    #else:
    #    msg = "FrameOfReferenceTransformationMatrixType is not one of the "\
    #          + f"expected values: 'RIGID', 'RIGID_SCALE' or 'AFFINE'."
    #    raise Exception(msg)
    
    PreDefTx = sitk.AffineTransform(3)
    PreDefTx.SetParameters(TxParams)
    
    return PreDefTx


def CreateSitkTxFromDefDro(Dro, RefIm=None, LogToConsole=False):
    """
    08/07/21:
        See create_tx_from_def_dro in dro_tools.dro_to_tx.py
    
    Create a SimpleITK BSpline transformation based on parameters parsed from
    a DICOM Deformable Spatial Registration Object (DRO).  
    
    Inputs:
    ******
    
    Dro : Pydicom object
        The DICOM Deformable Registration Object.
    
    RefIm : SimpleITK Image (optional; None by default)
        A reference image (e.g. fixed image) that shares the same direction 
        cosines as the BSpline grid.  If present, the direction of RefIm will
        be used to set that of the grid.  If None, the direction cosines along
        x and y, as stored in 
        
        DeformableRegistrationSequence[1]\
        .DeformableRegistrationGridSequence[0]\
        .ImageOrientationPatient
        
        will be used to determine the direction cosine along the z direction.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    SitkTx : SimpleITK Transform
        A transform based on the parameters stored in Dro.
    
    
    Notes:
    *****
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    In order to set the direction of the coefficient image below, the direction
    cosine along the z-direction is obtained by taking the cross product of the
    direction cosines along the x and y directions.  This is the simplest 
    solution but due to the 16 character limit when storing ImagePositionPatient
    and ImageOrientationPatient there are small differences from the values
    that would be obtained from FixIm.GetDirection().
    
    e.g. the difference between the direction cosines along x, y and z obtained
    from IOP and FixIm.GetDirection() for a sample data was:
        
    [-5.83795068e-09,  2.60326283e-08, -2.79367864e-08, 
    -5.55216258e-09, -1.82105636e-08,  2.78585151e-08,  
    2.62768308e-08, -2.83423783e-08, -2.48949777e-08]
    """
    
    import SimpleITK as sitk
    import numpy as np
    
    GridOrig = [float(item) for item in Dro.DeformableRegistrationSequence[1]\
                                           .DeformableRegistrationGridSequence[0]\
                                           .ImagePositionPatient]
    
    if LogToConsole:
        print(f'type(GridOrig) = {type(GridOrig)}')
        print(f'GridOrig = {GridOrig}\n')
    
    if RefIm == None:
        GridDir = [float(item) for item in Dro.DeformableRegistrationSequence[1]\
                                              .DeformableRegistrationGridSequence[0]\
                                              .ImageOrientationPatient]
        
        # Add the direction cosine along the z-direction:
        GridDir.extend(list(np.cross(np.array(GridDir[0:3]), np.array(GridDir[3:]))))
    else:
        GridDir = RefIm.GetDirection()
    
    if LogToConsole:
        print(f'type(GridDir) = {type(GridDir)}')
        print(f'GridDir = {GridDir}\n')
    
    GridSize = [item for item in Dro.DeformableRegistrationSequence[1]\
                                    .DeformableRegistrationGridSequence[0]\
                                    .GridDimensions]
    
    if LogToConsole:
        print(f'type(GridSize) = {type(GridSize)}')
        print(f'GridSize = {GridSize}\n')
    
    GridSpacing = [item for item in Dro.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .GridResolution]
    
    if LogToConsole:
        print(f'type(GridSpacing) = {type(GridSpacing)}')
        print(f'GridSpacing = {GridSpacing}\n')
    
    #dim = len(GridSize)
    #
    #MeshSize = [item - dim for item in GridSize]
    
    #VectGridData = [item for item in Dro.DeformableRegistrationSequence[1]\
    #                                    .DeformableRegistrationGridSequence[0]\
    #                                    .VectorGridData]
    
    #print(f'type(VectGridData) = {type(VectGridData)}')
    #print(f'VectGridData = {VectGridData}\n')
    
    # Convert from bytes to list of floats:
    FlatList = list(np.frombuffer(Dro.DeformableRegistrationSequence[1]\
                                     .DeformableRegistrationGridSequence[0]\
                                     .VectorGridData, dtype=np.float64))
    
    # Unflatten list to (x, y, z) tuples:
    ListOfTuples = zip(*[iter(FlatList)]*3)
    
    # Unzip from x, y, z to lists for x, y and z:
    CoeffX, CoeffY, CoeffZ = zip(*ListOfTuples)
    
    CoeffArrX = np.reshape(np.array(CoeffX), GridSize)
    CoeffArrY = np.reshape(np.array(CoeffY), GridSize)
    CoeffArrZ = np.reshape(np.array(CoeffZ), GridSize)
    
    if LogToConsole:
        print(f'CoeffArrX.shape = {CoeffArrX.shape}')
        print(f'CoeffArrY.shape = {CoeffArrY.shape}')
        print(f'CoeffArrZ.shape = {CoeffArrZ.shape}\n')
    
    CoeffImX = sitk.GetImageFromArray(CoeffArrX)
    CoeffImX.SetOrigin(GridOrig)
    CoeffImX.SetDirection(GridDir)
    CoeffImX.SetSpacing(GridSpacing)
    
    CoeffImY = sitk.GetImageFromArray(CoeffArrY)
    CoeffImY.SetOrigin(GridOrig)
    CoeffImY.SetDirection(GridDir)
    CoeffImY.SetSpacing(GridSpacing)
    
    CoeffImZ = sitk.GetImageFromArray(CoeffArrZ)
    CoeffImZ.SetOrigin(GridOrig)
    CoeffImZ.SetDirection(GridDir)
    CoeffImZ.SetSpacing(GridSpacing)
    
    SitkTx = sitk.BSplineTransform([CoeffImX, CoeffImY, CoeffImZ])
    
    return SitkTx




"""
******************************************************************************
******************************************************************************
EXPORT A DICOM SPATIAL OR DEFORMABLE SPATIAL REGISTRATION OBJECT (SRO) TO DISK
******************************************************************************
******************************************************************************
"""

def ExportDro(Dro, ExportDir, Fname='', DictOfInputs=None):
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
    
    Fname : string (optional; '' by default)
        Filename of exported DRO file.
    
    DictOfInputs : dictionary or None (optional; None by default)
        If not None, a dictionary containing the inputs that were called to 
        CopyRoi().
    
    
    Outputs:
    *******
    
    DroFpath : string
        Full path of the exported DRO.
    """
    
    import os
    import time
    from datetime import datetime
    from pathlib import Path
    
    if not os.path.isdir(ExportDir):
        #os.mkdir(ExportDir)
        Path(ExportDir).mkdir(parents=True)
        
    if Fname == '':
        if DictOfInputs == None:
            DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        else:
            DateTime = DictOfInputs['RunDateTime']
            
            """ Remove the decimal from the seconds (e.g. 20210611_145930.766799
            to 20210611_145930): """
            # Convert from string to datetime:
            dateTime = datetime.strptime(DateTime, "%Y%m%d_%H%M%S.%f")
            
            # Convert back to string without fractional seconds:
            DateTime = datetime.strftime(dateTime, "%Y%m%d_%H%M%S")
        
        Fname = DateTime + '_' + '.dcm'
    else:
        Fname = Fname.replace(' ', '_')
    
    if not '.dcm' in Fname:
        Fname += '.dcm'
    
    DroFpath = os.path.join(ExportDir, Fname)
    
    #print('Exporting DRO.')
    
    Dro.save_as(DroFpath)
        
    print(f'\nNew DICOM Registration Object exported to:\n {DroFpath}\n')
    
    return DroFpath