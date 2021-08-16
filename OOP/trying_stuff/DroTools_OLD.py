# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 12:51:41 2021

@author: ctorti
"""




def CreateDro_OLD(SrcDicomDir, TrgDicomDir, Description='', LogToConsole=False):
    """
    Create a DICOM Spatial Registration Object or Deformable Spatial
    Registration Object (DRO).  
    
    *NOTES*:  
        1. At the moment only the Spatial Registration Object Storage IOD is
        is catered for. 
        
        2. Parsing of the transform parameter file (TransformParameters.0.txt) 
        is done within this function since it doesn't seem possible to access 
        the parameters from the Elastix filter itself in SimpleElastix. The
        function RegisterImages() should be modified to parse the file and 
        return the transform parameters as an output.
        
        If SimpleITK is used to perform registration using RegisterImagesSitk()
        the transform paramters are returned as an output. Once the Elastix-
        based registration function has been modified to output the transform
        parameters, they should be returned by CopyRts() and CopySeg(),
        irrespective of whether SE or SITK is used to perform registration, so
        that they can be provided as input to CreateDro().
        
        
    
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
    
    """ Decide whether to reference all SOPInstanceUIDs in SrcDicoms and 
    TrgDicoms within ReferencedInstanceSequence (not required), or just one
    (for compactness), for the sake of maintaining the SeriesInstanceUIDs of 
    the two series. """
    #RefAllSOPs = True
    RefAllSOPs = False
    
    if RefAllSOPs:
        print('Note: All SOPInstanceUIDs of both series will be referenced',
              'in ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
    else:
        print('Note: Only the first SOPInstanceUID of both series will be',
              'referenced in ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
    
    # Start timing:
    times = []
    times.append(time.time())
    
    SrcDicoms = ImportDicoms(SrcDicomDir)
    TrgDicoms = ImportDicoms(TrgDicomDir)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if True:#LogToConsole:
        print(f'\n*Took {Dtime} s to import the 3D images.')
        
    # The Transform parameter file:
    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    # Parse the Transform parameter file:
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
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
       .FrameOfReferenceTransformationMatrix = TxParams
       
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
       .FrameOfReferenceTransformationMatrixType = GetTxMatrixType(LogToConsole)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
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
    if RefAllSOPs:
        ReqRIS = len(TrgDicoms)
        
        print(f'\nThere are {OrigRIS} sequences in',
              'ReferencedSeriesSequence[0].ReferencedInstanceSequence, and',
              f'{ReqRIS} (= no. of Target DICOMs) required sequences.')
    else:
        ReqRIS = 1
    
        print(f'\nThere are {OrigRIS} sequences in',
              'ReferencedSeriesSequence[0].ReferencedInstanceSequence, and',
              f'{ReqRIS} (the first Target DICOM) required sequence.')
    
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
    Dtime = round(times[-1] - times[-2], 3)
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
    if RefAllSOPs:
        ReqRIS = len(SrcDicoms)
        
        print(f'\nThere are {OrigRIS} sequences in',
              'ReferencedSeriesSequence[1].ReferencedInstanceSequence, and',
              f'{ReqRIS} (= no. of Moving DICOMs) required sequences.')
    else:
        ReqRIS = 1
    
        print(f'\nThere are {OrigRIS} sequences in',
              'ReferencedSeriesSequence[1].ReferencedInstanceSequence, and',
              f'{ReqRIS} (the first Moving DICOM) required sequence.')
    
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
    Dtime = round(times[-1] - times[-2], 3)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to update the ReferencedInstanceSequences',
              'for the moving series.')
        
    
    
    """ Delete the StudiesContainingOtherReferencedInstancesSequence. """
    
    print('\nDeleting the StudiesContainingOtherReferencedInstancesSequence.')
    
    del Dro[0x08, 0x1200]
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if True:#LogToConsole:
        print(f'\n   *Took {Dtime} s to delete',
              'StudiesContainingOtherReferencedInstancesSequence.')
    
    
    
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