# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:44:06 2020

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
GENERAL DICOM FUNCTIONS
******************************************************************************
******************************************************************************
"""


def GetDicomFpaths_OLD(DicomDir, SortMethod='slices', WalkThroughSubDirs=False,
                   LogToConsole=False):
    """
    NOTE: This function is less efficient than the new one, which uses 
    SimpleITK's ImageSeriesReader.
    
    Get a sorted list of full filepaths for all DICOM files in a directory
    
    Inputs:
    ******
    
    DicomDir : string
        Path to directory containing DICOMs.
        
    SortMethod : string (optional, 'slices' by default) 
        Denotes the method used for sorting of files fetched using 
        GetDicomFpaths() with acceptable values:
            -- 'none' (or 'None' or any other strings) = no sorting
            -- 'natural' (or 'Natural') = using natsort to sort the filenames
            -- 'slices' (or 'Slices') = sorted along the slice direction 
    
    WalkThroughSubDirs : boolean (optional; False by default)
        If True DICOMs will be searched by walking through sub-directories
        within DicomDir.  If False, the search will only happen at the depth
        of DicomDir.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    FilePaths : list of strings 
        List of the full filepaths of all DICOM files in DicomDir either:
            -- without sorting
            -- sorted in a natural way, 
            e.g. 3-1, 3-2, ..., 3-9, 3-10, 3-11, ...
            rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...
            -- sorted along the scan direction in ascending order
    """
    
    import os
    from natsort import natsorted
    from pydicom import read_file
    import numpy as np
    
    #print('\nSortMethod (input for GetDicomFpaths()) =', SortMethod)
    
    # Create an empty list of all Dicom file paths:
    FilePaths = []
    
    if WalkThroughSubDirs:
        """ Use os.walk to get list of file paths of .dcm files at any depth. """
        for DirName, SubDirList, FileList in os.walk(DicomDir):
            for FileName in FileList:
                if '.dcm' in FileName.lower():  # check for Dicom files only
                    FilePaths.append(os.path.join(DirName, FileName))
    else:
        FileList = os.listdir(DicomDir)
        
        for FileName in FileList:
            if '.dcm' in FileName.lower():  # check for Dicom files only
                FilePaths.append(os.path.join(DicomDir, FileName)) 
        
    #print('FilePaths =', FilePaths)
    
    if FilePaths==[]:
        print(f'\nNo DICOMs were found in {DicomDir}.')
    else:
        if SortMethod in ['natural', 'Natural']:
            # Sort files in natural way:
            FilePaths = natsorted(FilePaths)
            
        if SortMethod in ['slices', 'Slices']:
            # Create an empty array to store the Positions:
            Positions = []
            
            for FilePath in FilePaths:
                #print('\nFilePath =', FilePath)
                
                # Read in Dicom:
                Dicom = read_file(FilePath)
                
                # Parse the Image Position (Patient):
                Position = [float(item) for item in Dicom.ImagePositionPatient]
                
                # Add Position for this Dicom to Positions:
                Positions.append(Position)
                
            # Convert to numpy array:
            Positions = np.array(Positions)
            
            #print('Positions =', Positions)
    
            # Get the maximum change of Positions along x, y and z:
            dPositions = np.max(Positions, axis=0) - np.min(Positions, axis=0)
            
            if LogToConsole:
                print('\n   The maximum change of Positions along x, y and z',
                      f'= {dPositions}')
    
            # Get the index of the axis with the greatest change of Position:
            ind = np.where(dPositions==np.max(dPositions))[0][0]
            
            if LogToConsole:
                print('\n   The index of the axis with the greatest change of',
                      f'Position = {ind}')
            
            # Convert ind to 'x', 'y' or 'z':
            if ind == 0:
                SortAxis = 'x'
            if ind == 1:
                SortAxis = 'y'
            if ind == 2:
                SortAxis = 'z'
                
            if LogToConsole:
                print('\n   The DICOM filepaths will be sorted along the ' \
                      + SortAxis + '-axis.')
            
                print('\n   The pre-sorted ' + SortAxis + '-positions are:\n')
                print([round(Positions[i][ind], 2) for i in range(len(Positions))])
            
            """
            Solution on how to sort along a specific column:
            
            https://stackoverflow.com/questions/22698687/how-to-sort-2d-array-numpy-ndarray-based-to-the-second-column-in-python
            
            The following are equivalent:
            
            Positions[Positions[:, 2].argsort()]
            Positions[np.argsort(Positions[:, 1])]
            
            """
            
            SortInds = Positions[:, ind].argsort()
            
            #PositionsSorted = Positions[SortInds]
            
            if LogToConsole:
                print('\n   The sorted ' + SortAxis + '-positions are:\n')
                print([round(Positions[i][ind], 2) for i in SortInds])
            
            # Sort FilePaths with SortInds:
            #FilePaths = FilePaths[SortInds] # <-- TypeError: list indices must be 
            # integers or slices, not list
            FilePaths = [FilePaths[i] for i in SortInds]
        
    
    return FilePaths








def GetDicomFpaths(DicomDir):
    """
    Get a list of full filepaths for all DICOM files in a directory using
    SimpleITK's ImageSeriesReader.
    
    Inputs:
    ******
    
    DicomDir : string
        Path to directory containing DICOMs.
    
    
    Outputs:
    *******
    
    FilePaths : list of strings 
        List of the full filepaths of all DICOM files in DicomDir.
    """
    
    import SimpleITK as sitk
    
    SeriesReader = sitk.ImageSeriesReader()
    
    FilePaths = SeriesReader.GetGDCMSeriesFileNames(DicomDir)
    
    return FilePaths


    
    


def ImportDicoms(DicomDir, SortMethod='slices', LogToConsole=False):
    """
    Import DICOM objects from a directory containing DICOM files.
    
    Inputs:
    *******
    
    DicomDir : string
        Path to directory containing DICOM files.
        
    SortMethod : string (optional, 'slices' by default) 
        Denotes the method used for sorting of files fetched using 
        GetDicomFpaths() with acceptable values:
            -- 'none' (or 'None' or any other strings) = no sorting
            -- 'natural' (or 'Natural') = using natsort to sort the filenames
            -- 'slices' (or 'Slices') = sorted along the slice direction 
        
    Outputs:
    *******
    
    Dicoms : list of Pydicom objects
        List of DICOM objects.
    """
    
    from pydicom import read_file
    
    #print('\nSortMethod (input for GetDicoms()) =', SortMethod)
    
    #fpaths = GetDicomFpaths(DicomDir, SortMethod, LogToConsole)
    fpaths = GetDicomFpaths(DicomDir)
    
    Dicoms = []
    
    for f in range(len(fpaths)):
        fpath = fpaths[f]
        
        dicom = read_file(fpath)
        
        Dicoms.append(dicom)
        
    return Dicoms






def ImportDicom(DicomDir, SortMethod='slices', InStackNum=0, 
                LogToConsole=False):
    """
    Import a single DICOM object from a directory containing DICOM file(s).
    
    Inputs:
    *******
    
    DicomDir : string
        Path to directory containing DICOM files.
        
    SortMethod : string (optional, 'slices' by default) 
        Denotes the method used for sorting of files fetched using 
        GetDicomFpaths() with acceptable values:
            -- 'none' (or 'None' or any other strings) = no sorting
            -- 'natural' (or 'Natural') = using natsort to sort the filenames
            -- 'slices' (or 'Slices') = sorted along the slice direction
            
    InStackNum : integer (optional; 0 by default)
        The in-stack index of the file within the sorted list of DICOM 
        filepaths to import.  The first file is imported by default.
    
        
    Outputs:
    *******
    
    Dicoms : list of Pydicom objects
        List of DICOM objects.
    """
    
    from pydicom import read_file
    
    #fpaths = GetDicomFpaths(DicomDir, SortMethod, LogToConsole)
    fpaths = GetDicomFpaths(DicomDir)
    
    if InStackNum < len(fpaths):
        Dicom = read_file(fpaths[InStackNum])
        
        return Dicom
    
    else:
        msg = 'It is not possible to import a DICOM at InStackNum = '\
              + f'{InStackNum} since there are only {len(fpaths)} DICOMs in '\
              + 'the specified directory.'
              
        raise Exception(msg)
        
        return None






def GetDicomSOPuids_OLD(DicomDir):
    """
    NOTE: This function is less efficient than the new one, which uses 
    SimpleITK's ImageSeriesReader.
    
    Get a list of SOP Instance UIDs for a DICOM series.
    
    Inputs:
    ******
    
    DicomDir : string
        Directory containing the DICOMs.
        
    
    Outputs:
    *******
    
    SOPuids : list of strings
        List of the SOP UIDs of the DICOMs.
    """
    
    Dicoms = ImportDicoms(DicomDir=DicomDir, SortMethod='slices', 
                          LogToConsole=False)
    
    SOPuids = [] 
    
    for dicom in Dicoms:
        SOPuids.append(dicom.SOPInstanceUID)
        
    return SOPuids






def GetDicomSOPuids(DicomDir):
    """
    Get a list of SOP Instance UIDs for a DICOM series using SimpleITK's 
    ImageSeriesReader.
    
    Inputs:
    ******
    
    DicomDir : string
        Directory containing the DICOMs.
        
    
    Outputs:
    *******
    
    SOPuids : list of strings
        List of the SOP UIDs of the DICOMs.
    """
    
    import SimpleITK as sitk
    
    SeriesReader = sitk.ImageSeriesReader()
    
    FilePaths = SeriesReader.GetGDCMSeriesFileNames(DicomDir)
    
    FileReader = sitk.ImageFileReader()
    
    SOPuids = [] 
    
    for fpath in FilePaths:
        FileReader.SetFileName(fpath)

        FileReader.ReadImageInformation()
    
        SOPuids.append(FileReader.GetMetaData('0008|0018'))
        
    return SOPuids








def GetDicomUids_OLD(DicomDir):
    """
    NOTE: This function is less efficient than the new one, which uses 
    SimpleITK's ImageSeriesReader.
    
    Get the Study, Series, Frame of Reference and SOP Instance UIDs for a DICOM 
    series.
    
    Inputs:
    ******
    
    DicomDir : string
        Directory containing the DICOMs.
        
    
    Outputs:
    *******
    
    Studyuid : string
        Study Instance UID of the DICOM series.
        
    Seriesuid : string
        Series Instance UID of the DICOM series.
        
    FORuid : string
        Frame of Reference UID of the DICOM series.
        
    SOPuids : list of strings
        List of the SOP UIDs of the DICOMs.
    """
    
    Dicoms = ImportDicoms(DicomDir=DicomDir, SortMethod='slices', 
                          LogToConsole=False)
    
    SOPuids = [] 
    
    for dicom in Dicoms:
        SOPuids.append(dicom.SOPInstanceUID)
        
    Studyuid = Dicoms[0].StudyInstanceUID
    Seriesuid = Dicoms[0].SeriesInstanceUID
    FORuid = Dicoms[0].FrameOfReferenceUID
        
    return Studyuid, Seriesuid, FORuid, SOPuids






def GetDicomUids(DicomDir):
    """
    Get the Study, Series, Frame of Reference and SOP Instance UIDs for a DICOM 
    series.
    
    Inputs:
    ******
    
    DicomDir : string
        Directory containing the DICOMs.
        
    
    Outputs:
    *******
    
    Studyuid : string
        Study Instance UID of the DICOM series.
        
    Seriesuid : string
        Series Instance UID of the DICOM series.
        
    FORuid : string
        Frame of Reference UID of the DICOM series.
        
    SOPuids : list of strings
        List of the SOP UIDs of the DICOMs.
    """
    
    import SimpleITK as sitk
    
    SeriesReader = sitk.ImageSeriesReader()
    
    FilePaths = SeriesReader.GetGDCMSeriesFileNames(DicomDir)
    
    FileReader = sitk.ImageFileReader()
    
    SOPuids = [] 
    
    for fpath in FilePaths:
        FileReader.SetFileName(fpath)

        FileReader.ReadImageInformation()
    
        SOPuids.append(FileReader.GetMetaData('0008|0018'))
        
    Studyuid = FileReader.GetMetaData('0020|000d')
    Seriesuid = FileReader.GetMetaData('0020|000e')
    FORuid = FileReader.GetMetaData('0020|0052')
        
    return Studyuid, Seriesuid, FORuid, SOPuids







def GetRoiLabels(Roi):
    """
    Get the list of contour/segment labels in a RTS/SEG.
    
    Inputs:
    ******
    
    Roi : Pydicom object
        Object from a RTS/SEG file.
    
        
    Outputs:
    *******
    
    Labels : list of strings
        List (for each ROI/segment) of ROIName/SegmentLabel tag values.
    """
    
    Labels = [] 
    
    if 'RTSTRUCT' in Roi.Modality:
        sequences = Roi.StructureSetROISequence
        
        for sequence in sequences:
            label = sequence.ROIName
            
            Labels.append(label)
            
    elif 'SEG' in Roi.Modality:
    
        sequences = Roi.SegmentSequence
        
        for sequence in sequences:
            label = sequence.SegmentLabel
            
            Labels.append(label) 
            
    else:
        msg = "The ROI object is not recognised as being RTS or SEG."
        
        raise Exception(msg)
        
    return Labels








def GetRoiNums(Roi, SearchString):
    """
    Get the ROI/segment number(s) that matches SearchString.
    
    Inputs:
    ******
    
    Roi : Pydicom object
        RTS or SEG object.
        
    SearchString : string
        All or part of the ROIName/SegmentLabel of the ROI/segment of interest.
        If SearchString is an empty string return the first ROI/segment number. 
                       
                            
    Outputs:
    *******
    
    RoiNums : list of integers
        A list of indices corresponding to the ROI(s) that matche SearchString.
    """
    
    RoiLabels = GetRoiLabels(Roi)
    
    #print('\nRoiLabels =', RoiLabels)
    
    #print('\nRoiLabel to find =', SearchString)
    
    if not RoiLabels:
        msg = 'There are no ROIs in the RTS/SEG.'
            
        raise Exception(msg)
    
    
    if SearchString == "":
        RoiNum = 0
        
    else:
        #RoiNum = RoiLabels.index(SearchString)
        #RoiNum = [i for i, RoiLabel in enumerate(RoiLabels) if SearchString in RoiLabel][0] # 13/01/2021
        RoiNums = [i for i, RoiLabel in enumerate(RoiLabels) if SearchString in RoiLabel] # 13/01/2021
        
        #print(f'\nRoiNum = {RoiNum}')
        
        #if RoiNum:
        #    RoiNum = RoiNum[0] # 13/01/2021
        #    
        #else:
        #    msg = "There is no ROI/segment in the RTS/SEG containing "\
        #          + f"names/labels {RoiLabels} \nthat matches 'SearchString' "\
        #          + f"= '{SearchString}'."
        #    
        #    raise Exception(msg)
            
        if not RoiNums:
            msg = "There is no ROI/segment in the RTS/SEG containing "\
                  + f"names/labels {RoiLabels} \nthat matches 'SearchString' "\
                  + f"= '{SearchString}'."
            
            raise Exception(msg)
            
    return RoiNums









def IsSameModalities(Roi0, Roi1):
    """
    Compare Modality tag of two RTS or two SEG objects. 
    """

    if Roi0.Modality == Roi1.Modality:
        return True
    
    else:
        return False










def InspectDicomStudyDir(StudyDir):
    import os
    from natsort import natsorted
    from ImageTools import GetImageAttributes
    
    # First get a list of DICOM scan directories
    ScanDirs = []
    ScanDirNames = []

    # Use os to get list of file paths of Dicom-only files:
    for DirName, SubDirList, FileList in os.walk(StudyDir):
        #print('DirName =', DirName)
        #print('SubDirList =', SubDirList)
        #print('FileList =', FileList)
        #print('\n\n')
        
        DicomsFound = 0
        
        for file in FileList:
            if '.dcm' in file:
                #print(FileList, 'contains at least 1 DICOM')
                DicomsFound += 1
                
        
        if DicomsFound:
            ScanDirs.append(DirName)
            
            # DirName is a directory containing DICOMs
            RootDir, ScanDir = os.path.split(DirName)
    
            ScanDirNames.append(ScanDir)
    
    
    # Sort the list of directory names using natsort:
    ScanDirNames = natsorted(ScanDirNames)
    
    """ 
    The DICOM scan directory names have the form:
        'N-Label'
        
    e.g.:
        '2-T1 SE SAG IPATX2-05296'
        
    where scan number = 2
    and scan label = 'T1 SE SAG IPATX2-05296'
    
    Get the scan numbers and labels separately.
    """
    ScanNums = []
    ScanLabels = []
    
    for ScanDirName in ScanDirNames:
        # Index of first '-':
        ind = ScanDirName.index('-')
        
        ScanNums.append(ScanDirName[0:ind])
        
        ScanLabels.append(ScanDirName[ind+1::])
    
    
    # Initialise the Origins, Directions, Spacings, Dims, and FORuids:
    Origins = []
    Dirs = []
    Spacings = []
    Size = []
    FORuids = []
    
    # Loop through each Scan directory:
    for ScanDir in ScanDirs:
    
        # Get the Image Attributes:
        Studyuid, Seriesuid, FORuid, size, spacings, ST,\
        IPPs, dirs, ListOfWarnings = GetImageAttributes(ScanDir)
        
        # Get the DICOMs:
        DicomFpaths, Dicoms = ImportDicoms(DirPath=ScanDir, 
                                           SortMethod='slices', 
                                           Debug=False)
        
        Origins.append(IPPs[0])
        Dirs.append(dirs)
        Spacings.append(spacings)
        Size.append(size)
        FORuids.append(FORuid)
        
    
    # Create a dictionary:
    StudyDict = {'Scan num':ScanNums,
                 'Scan label':ScanLabels,
                 #'FOR UID':FORuids,
                 'Size':Size,
                 'Spacings':Spacings,
                 'Origin':Origins,
                 'Directions':Dirs
                 }
    
    #return StudyDict, ScanNums, ScanLabels, FORuids, Origins, Dirs, Dims, Spacings
    return StudyDict
    
    



def InspectDicomSubjectDir(SubjectDir): # THIS IS NOT COMPLETE!
    """
    This is not complete.
    """
    
    import os
    from natsort import natsorted
    
    # Get a list of study directories:
    StudyDirs = os.list(SubjectDir)
    
    # Sort using natsort:
    StudyDirs = natsorted(StudyDirs)
    
    
    """
    The directory names of the study directories have the form:
        'MM-DD-YYYY-MOD StudyLabel'
        
    where MOD is the modality (e.g. 'MRI', 'PET', 'CT')
    
    Get the date, modality and study label separately.
    """
    
    StudyDates = []
    Modalities = []
    StudyLabels = []
    
    for StudyDir in StudyDirs:
        # Get the study directory names:
        StudyDirRoot, StudyDirName = os.path.split(StudyDir)

        #StudyDates.append()
    

    return





#def DoesDirContainDicoms(Directory):
#    """ Check if a directory contains any DICOMs. Return True if there is one
#    or more .dcm files, False otherwise. Also return the number of files. """
#    
#    import os
#    
#    FileList = os.listdir(Directory)  
    




def ConvertNumpyDtype(ArrayIn, ConvertTo, LogToConsole=False):
    """
    Convert Numpy array from one data type to another, preserving the
    dynamic range.
    
    Inputs:
    ******
    
    ArrayIn : Numpy array 
        Numpy array whose data type is to be converted.  Acceptable data types
        include:
          - 'uint8'
          - 'uint16'
          - 'int8'
          - 'int16'
          - 'float32'
          - 'float64'
                        
    ConvertTo : string
        String denoting the data type to convert ArrayIn to.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
        
    Outputs:
    ********
    
    ArrayOut : Numpy array
        Numpy array with converted data type.
    """
    import numpy as np
    
    # Get info on input array, first assuming that it is an int dtype:
    try:
        info = np.iinfo(ArrayIn.dtype)
        
    # If ValueError arises, assume it is a float:
    except ValueError:
        info = np.finfo(ArrayIn.dtype)
      
    
    # Determine value to scale array up to and desired Numpy data type:
    if ConvertTo == 'uint8':
        MinValue = 0
        MaxValue = 255
        NpDtype = np.uint8
    
    elif ConvertTo == 'uint16':
        MinValue = 0
        MaxValue = 65535
        NpDtype = np.uint16
    
    elif ConvertTo == 'int8':
        MinValue = -128
        MaxValue = 127
        NpDtype = np.int8
    
    elif ConvertTo == 'int16':
        MinValue = -32768
        MaxValue = 32767
        NpDtype = np.int16
        
    elif ConvertTo == 'float32':
        MinValue = -3.4028235e+38
        MaxValue = -MinValue
        NpDtype = np.float32
    
    elif ConvertTo == 'float64':
        MinValue = -1.7976931348623157e+308
        MaxValue = -MinValue
        NpDtype = np.float64
    
    else:
        msg = f"ConvertTo (= {ConvertTo}) must be 'uint8', 'uint16', 'int8', "\
              + "'int16', 'float32' or 'float64'."
        raise Exception(msg)
        
    
    
    if LogToConsole:    
        print('Before dtype:', ArrayIn.dtype, '\n')
        print(info, '\n')
        print(f'ConvertTo = {ConvertTo}\n')
        print(np.iinfo(NpDtype), '\n')
    
    
    #if LogToConsole:    
    #    print(f'Max value possible for {ConvertTo} = {MaxValue}\n')
    
    # Get min/max values:
    ArrayMin = np.min(ArrayIn)
    ArrayMax = np.max(ArrayIn)
    
    if LogToConsole:
        print(f'Min/Max of input array = [{ArrayMin}, {ArrayMax}]\n')
    
    
    if ArrayMin < MinValue or ArrayMax > MaxValue:
        # Deal with negative values if ConvertTo is an unsigned integer
        # (by subtracting the min value):
        if 'u' in ConvertTo:
            Array = ArrayIn - ArrayMin
    
            # Get min/max values:
            ArrayMin = np.min(Array)
            ArrayMax = np.max(Array)
            
            if LogToConsole:
                print('Min/Max of array after subtracting Min =',
                      f'[{ArrayMin}, {ArrayMax}]\n')
            
        else:
            Array = ArrayIn - 0
        
        # Normalise Array to [0, 1] or [-1, 1] (depends on whether the input array
        # is signed or unsigned):
        #Array = ArrayIn.astype(np.float64) / info.max # is conversion to float useful?
        #Array = ArrayIn / info.max 
        Array = Array / MaxValue
        
        # Get min/max values:
        ArrayMin = np.min(Array)
        ArrayMax = np.max(Array)
        
        if LogToConsole:
            print(f'Min/Max of Array after normalising = [{ArrayMin}, {ArrayMax}]\n')
            
        
        # Scale Array by MaxValue:
        Array = MaxValue * Array
        
        # Get min/max values:
        ArrayMin = np.min(Array)
        ArrayMax = np.max(Array)
        
        if LogToConsole:
            print(f'Min/Max of array after scaling by {MaxValue} =',
                  f'[{ArrayMin}, {ArrayMax}]\n')
    
        # Convert Array to desired dtype:
        ArrayOut = Array.astype(NpDtype)
    else:
        # Convert ArrayIn to desired dtype:
        ArrayOut = ArrayIn.astype(NpDtype)
    
    
    if LogToConsole:
        # Get min/max values:
        ArrayMin = np.min(ArrayOut)
        ArrayMax = np.max(ArrayOut)
    
        print(f'Min/Max of array after converting to {ConvertTo}',
              f'= [{ArrayMin}, {ArrayMax}]\n')
    
        # Re-check info:
        try:
            info = np.iinfo(ArrayOut.dtype)
            
        # If ValueError arises, assume it is a float:
        except ValueError:
            info = np.finfo(ArrayOut.dtype)
        
        print('After dtype:', ArrayOut.dtype, '\n')
        #print(f'After: {info}\n')
    
    return ArrayOut





def CreateDicomsFromCroppedSitkIm(CroppedIm, FirstZind, OrigDicomDir, 
                                  TxtToAddToSeriesDesc, NewSeriesNum,
                                  SortMethod='slices', LogToConsole=False):
    """
    27/05/21: The DICOMs created are not displaying properly, suggesting an 
    issue with the conversion to bytes. I've not figured out what the problem
    is.
    
    Create new DICOM objects from a SimpleITK image using as templates, the 
    original DICOMs (that were imported as a SimpleITK image and modified in
    some way).
    
    Inputs:
    *******
    
    CroppedIm : SimpleITK image
        The cropped image.
    
    FirstZind : integer
        The slice number within the original DICOM series that corresponds to
        the first frame in CroppedIm.
    
    OrigDicomDir : string
        Path to directory containing the original DICOM files.
    
    TxtToAddToSeriesDesc : string
        Text that will be added to the original SeriesDescription for all
        new DICOMs (e.g. 'cropped').
    
    NewSeriesNum : string
        The new series number to assign to the new DICOM series.
    
        
    Outputs:
    *******
    
    Dicoms : list of Pydicom objects
        List of DICOM objects related to the frames in SitkIm.
    
    
    Note:
    ****
    
    Only the following tags are modified (all else are either the same or not
    significant enough to change):
        SOPInstanceUID
        MediaStorageSOPInstanceUID
        ContentDate
        ContentTime
        SeriesDescription
        SeriesInstanceUID
        SeriesNumber
        ImagePositionPatient
        SliceLocation
        Rows
        Columns
        PixelData
    """
    
    import datetime
    #import numpy as np
    #from pydicom import dcmread
    from pydicom.uid import generate_uid
    #from pydicom.pixel_data_handlers.numpy_handler import pack_bits
    from pydicom.dataset import validate_file_meta
    import SimpleITK as sitk
    
    # Convert CroppedIm to int16:
    #CroppedIm = sitk.Cast(CroppedIm, sitk.sitkInt16)
    
    OrigDicoms = ImportDicoms(OrigDicomDir, SortMethod, LogToConsole)
    
    if False:
        if len(OrigDicoms) != CroppedIm.GetSize()[2]:
            msg = f"The number of original DICOMs, {len(OrigDicoms)}, does not "\
                  + f"match the number of frames in CroppedIm, {CroppedIm.GetSize()[2]}."
            raise Exception(msg)
            
            """ Get the list of all z-positions for the original series: """
            OrigZs = []
            
            for dicom in OrigDicoms:
                OrigZs.append(float(dicom.ImagePositionPatient[2]))
    
    
    TimeNow = datetime.datetime.now()
    CurrentDate = TimeNow.strftime('%Y%m%d')
    CurrentTime = TimeNow.strftime('%H%M%S.%f')
    
    # Generate a new SeriesInstanceUID:
    NewSeriesInstanceUID = generate_uid()
    
    R, C, F = CroppedIm.GetSize()
    
    if type(NewSeriesNum) != str:
        NewSeriesNum = str(NewSeriesNum)
        
    # Get the data type of the original DICOMs:
    OrigDtype = str(OrigDicoms[0].pixel_array.dtype)
    
    print(f'OrigDtype = {OrigDtype}\n')
    
    Dicoms = []
    
    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots(F, 1, figsize=(14, 7*F))
    #n = 1 # sub-plot number
    
    for f in range(F):
        dc = OrigDicoms[f + FirstZind]
        
        # Generate a new SOPInstanceUID:
        dc.SOPInstanceUID = generate_uid()
        dc.file_meta.MediaStorageSOPInstanceUID = dc.SOPInstanceUID
        
        dc.ContentDate = CurrentDate
        dc.ContentTime = CurrentTime
        
        dc.SeriesDescription +=  TxtToAddToSeriesDesc
        
        dc.SeriesInstanceUID = NewSeriesInstanceUID
        
        dc.SeriesNumber = NewSeriesNum
        
        newIPP = CroppedIm.TransformContinuousIndexToPhysicalPoint([0,0,f])
        
        dc.ImagePositionPatient = [str(item) for item in newIPP]
        
        dc.SliceLocation = str(newIPP[2])
        
        dc.Rows = R
        dc.Columns = C
        
        PixArr = sitk.GetArrayFromImage(CroppedIm[:,:,f])
        
        print(f'PixArr.shape = {PixArr.shape}\n')
        print(f'PixArr.dtype = {PixArr.dtype}\n')
        print(f'PixArr.tobytes()[0:50] = {PixArr.tobytes()[0:50]}\n')
        print(f'PixArr.tostring()[0:50] = {PixArr.tostring()[0:50]}\n')
        
        #ax = plt.subplot(F, 1, n)
        #ax.imshow(PixArr, cmap=plt.cm.Greys_r);
        #n += 1 # increment sub-plot number
        
        #""" Convert PixArr to bytes. """
        if True:
            # The following seem to be equivalent but not sure they will 
            # result in required int16 dtype:
            #dc.PixelData = PixArr.tobytes()
            dc.PixelData = PixArr.tostring()
            
            # https://stackoverflow.com/questions/14350675/create-pydicom-file-from-numpy-array:
            validate_file_meta(dc.file_meta, enforce_standard=True)
        
        
        if False:
            print(f'PixArr.dtype = {PixArr.dtype}\n')
            
            # Convert PixArr from PixArr.dtype (e.g. float32) to OrigDtype 
            # (e.g. int16):
            PixArr = ConvertNumpyDtype(PixArr, OrigDtype, LogToConsole)
            
            dc.PixelData = PixArr.tostring()
        
        Dicoms.append(dc)
        
    return Dicoms





def ExportDicoms(Dicoms, RootExportDir):
    """
    Export a list of DICOM Objects to disk.  
    
    Inputs:
    ******
    
    Dicoms : list of Pydicom objects
        A list of DICOM objects to export.
        
    RootExportDir : string
        Root directory path to where the new DICOMs are to be exported. The
        directory name within RootExportDir will have the form 
        "SeriesNumber-SeriesDescription".
    
    
    Outputs:
    *******
    
    None
    
    
    Note:
    ****
    
    The exported DICOMs will have filenames of the form "000001.dcm",
    "000002.dcm", ...
    """
    
    import os
    from pathlib import Path
    
    SeriesNum = Dicoms[0].SeriesNumber
    
    SeriesDesc = Dicoms[0].SeriesDescription
    
    ExportDir = os.path.join(RootExportDir, f'{SeriesNum}-{SeriesDesc}')
    
    if not os.path.isdir(ExportDir):
        #os.mkdir(ExportDir)
        Path(ExportDir).mkdir(parents=True)
    
    
    for i in range(len(Dicoms)):
        
        Fname = f'{i+1:06d}.dcm'
    
        Fpath = os.path.join(ExportDir, Fname)
    
        Dicoms[i].save_as(Fpath)
        
        # https://pydicom.github.io/pydicom/dev/tutorials/dataset_basics.html?highlight=enforce_standard:
        Dicoms[i].save_as(filename=Fpath, write_like_original=False) 
        
        #print(f'\nNew DICOM Object exported to:\n {Fpath}\n')
    
    return