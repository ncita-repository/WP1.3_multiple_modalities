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





def CreateDro(SrcDicomDir, TrgDicomDir, TxParams, Description='', 
              LogToConsole=False):
    """
    Create a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO).  
    
    *NOTE*:  At the moment only the Spatial Registration Object Storage IOD is
             catered for. 
        
    
    Inputs:
    ******
    
    SrcDicomDir : string
        Directory containing Source DICOMs.
    
    TrgDicomDir : string
        Directory containing Target DICOMs.
    
    TxParams : list of strings
        The registration transform parameters.
    
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
    
    """
    
    import os
    import time
    from pydicom import dcmread
    from pydicom.uid import generate_uid
    from copy import deepcopy
    from DicomTools import ImportDicoms
    from ImageTools import GetTxMatrixType
    
    cwd = os.getcwd()

    DroDir = os.path.join(cwd, 'sample_DRO')
    
    NewSpaDroFname = 'sample_spatial_dro.dcm'
    NewSpaDroFpath = os.path.join(DroDir, NewSpaDroFname)
    
    NewDefDroFname = 'sample_deformable_dro.dcm'
    NewDefDroFpath = os.path.join(DroDir, NewDefDroFname)
        
    
    """ Check if the sample DROs have already been downloaded: """   
    if os.path.isfile(NewSpaDroFpath) and os.path.isfile(NewDefDroFpath):
        print(f'Sample DROs already downloaded to:\n {DroDir}\n')
    else:
        url = 'https://www.insight-journal.org/download/fulldownload/zip/923/1'
        print(f'Downloading DROs from {url}...\n')
        DownloadSampleDros()
    
    
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
    
    
    """
    If the registration was non-elastic use Dro1 as a template.
    If the registration was elastic use Dro2.
    
    For now use if True statement, to be replaced with more intelligent logic.
    """
    if True:
        Dro = dcmread(NewSpaDroFpath)
    else:
        Dro = dcmread(NewDefDroFpath)
    
    
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
    #Dro.ContentDate = TrgDicoms[0].ContentDate
    Dro.ContentDate = CurrentDate
    Dro.StudyTime = TrgDicoms[0].StudyTime
    Dro.SeriesTime = TrgDicoms[0].SeriesTime
    #Dro.ContentTime = TrgDicoms[0].ContentTime
    Dro.ContentTime = CurrentTime
    Dro.Manufacturer = TrgDicoms[0].Manufacturer
    """ Consider modifying StudyDescription (= '' in the template DRO. """
    Dro.SeriesDescription = TrgDicoms[0].SeriesDescription
    Dro.ManufacturerModelName = TrgDicoms[0].ManufacturerModelName
    Dro.PatientName = TrgDicoms[0].PatientName
    Dro.PatientID = TrgDicoms[0].PatientID
    Dro.PatientBirthDate = TrgDicoms[0].PatientBirthDate
    Dro.PatientSex = TrgDicoms[0].PatientSex
    Dro.PatientAge = TrgDicoms[0].PatientAge
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
       .FrameOfReferenceTransformationComment = 
    """
    
    
    
    """ Modify the ReferencedSeriesSequence for the fixed domain. """
    if LogToConsole:
        print('\nModify the ReferencedSeriesSequence for the fixed domain.')
    
    Dro.ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = TrgDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(TrgDicoms)
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  'ReferencedSeriesSequence[0].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (= no. of Target DICOMs) required sequences.')
    else:
        ReqRIS = 1
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  'ReferencedSeriesSequence[0].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (the first Target DICOM) required sequence.')
    
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
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
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
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the fixed series.')
    
    
    
    
    """ Modify the ReferencedSeriesSequence for the moving domain. """
    if LogToConsole:
        print('\nModify the ReferencedSeriesSequence for the moving domain.')
    
    Dro.ReferencedSeriesSequence[1]\
       .SeriesInstanceUID = SrcDicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[1]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(SrcDicoms)
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  'ReferencedSeriesSequence[1].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (= no. of Moving DICOMs) required sequences.')
    else:
        ReqRIS = 1
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  'ReferencedSeriesSequence[1].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (the first Moving DICOM) required sequence.')
    
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
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
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
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the moving series.')
        
    
    
    """ Delete the StudiesContainingOtherReferencedInstancesSequence. """
    if LogToConsole:
        print('\nDeleting the StudiesContainingOtherReferencedInstancesSequence.')
    
    del Dro[0x08, 0x1200]
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to delete',
              'StudiesContainingOtherReferencedInstancesSequence.')
    
    
    
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
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to delete all GroupLength tags.')
        
    
    Dtime = round(times[-1] - times[0], 1)
    if LogToConsole:
        print(f'\n   *Took a total of {Dtime} s to create the DRO.')
        
        
    return Dro





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
    from pathlib import Path
    
    if not os.path.isdir(ExportDir):
        #os.mkdir(ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    if not DictOfInputs == None:
        DateTime = DictOfInputs['RunDateTime']
    else:
        DateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
    if Fname == '':
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


