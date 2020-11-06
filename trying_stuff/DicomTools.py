# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:44:06 2020

@author: ctorti
"""

"""
Various DICOM tools
"""


# Import packages:
import os
import pydicom
import SimpleITK as sitk
import numpy as np
#import importlib
import natsort
#import copy




def GetDicomFpaths(DicomDir, SortMethod='slices', LogToConsole=None):
    """
    Get a sorted list of full filepaths for all DICOM files in a directory
    
    Input:
        DicomDir   - (String) Path to directory containing DICOMs
        
        SortMethod - (String) Method to use for sorting of files;
                    Acceptable values are:
                    -- 'none' (or 'None' or any other strings) = no sorting
                    -- 'natural' (or 'Natural') = using natsort
                    -- 'slices' (or 'Slices') = sorted along the slice 
                    direction in ascending order
        
    Returns:
        FilePaths - full filepaths of all Dicom files in DicomDir either:
                    -- without sorting
                    -- sorted in a natural way, 
                    e.g. 3-1, 3-2, ..., 3-9, 3-10, 3-11, ...
                    rather than 3-1, 3-10, 3-11, ... , 3-2, 3-20, ...
                    -- sorted along the scan direction in ascending order
                    
    """
    
    #print('\nSortMethod (input for GetDicomFpaths()) =', SortMethod)
    
    # Create an empty list of all Dicom file paths:
    FilePaths = []

    # Use os to get list of file paths of Dicom-only files:
    for DirName, SubDirList, FileList in os.walk(DicomDir):
        for FileName in FileList:
            if '.dcm' in FileName.lower():  # check for Dicom files only
                FilePaths.append(os.path.join(DirName,FileName))
                
    #print('FilePaths =', FilePaths)
    
    if FilePaths==[]:
        print('No DICOMs found.')
    else:
        if SortMethod in ['natural', 'Natural']:
            # Sort files in natural way:
            FilePaths = natsort.natsorted(FilePaths)
            
        if SortMethod in ['slices', 'Slices']:
            # Create an empty array to store the Positions:
            Positions = []
            
            for FilePath in FilePaths:
                #print('\nFilePath =', FilePath)
                
                # Read in Dicom:
                Dicom = pydicom.read_file(FilePath)
                
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










def ImportDicoms(DicomDir, SortMethod='slices', LogToConsole=False):
    """
    Import DICOM objects from a directory containing DICOM files.
    
    Input:
        DicomDir   - (String) Path to directory containing DICOM files
        
        SortMethod - (String) Method to use for sorting of files using 
                    GetDicomFpaths() with acceptable values:
                    -- 'none' (or 'None' or any other strings) = no sorting
                    -- 'natural' (or 'Natural') = using natsort to sort the 
                    filenames
                    -- 'slices' (or 'Slices') = sorted along the slice direction 
        
    Returns:
        Dicoms    - List of DICOM objects
    """
    
    
    #print('\nSortMethod (input for GetDicoms()) =', SortMethod)
    
    # Get the filepaths of the DICOM files:
    fpaths = GetDicomFpaths(DicomDir, SortMethod, LogToConsole)
    
    # Initialise a dictionary to store everything:
    #dicomDict = {}
    
    # Initialise a list to store the DICOM objects:
    Dicoms = []
    
    # Loop though each DICOM filepath:
    for f in range(len(fpaths)):
        fpath = fpaths[f]
        
        #print('\nfpath =', fpath)
        
        # Read in the DICOM:
        dicom = pydicom.read_file(fpath)
        
        # Append this DICOM object:
        Dicoms.append(dicom)
        
        # Get the SOP Instance UID:
        #uid = dicom.SOPInstanceUID
        
        # Get the Image Position (Patient):
        #position = [float(item) for item in dicom.ImagePositionPatient]
        
        # Get the Image Orientation Patient:
        #orientation = [float(item) for item in dicom.ImageOrientationPatient]
        
        # Get the Slice Location:
        #sliceLocation = float(dicom.SliceLocation)
        
        # Get the Pixel Spacing:
        #pixSpacing = [float(item) for item in dicom.PixelSpacing]
        
        # Get the Pixel Array:
        #pixelArray = dicom.pixel_array
    
        # Store the values in the dictionary using the SOP Instance UID as keys:
        #dicomDict[uid] = {'Dicom slice no':f+1, \
        #                  'Dicom fpath':fpath, \
        #                  'Image pos':position, \
        #                  'Image orient':orientation, \
        #                  'Slice loc':sliceLocation, \
        #                  'Pix spacing':pixSpacing, \
        #                  'Pix array':pixelArray, \
        #                  'Dicom':dicom
        #                 }
    
    #return dicomDict
    #return fpaths, dicoms
    return Dicoms







def ImportImage(DicomDir):
    """
    Import a DICOM series as a SimpleITK multi-dimensional image.
    
    Inputs:
        DicomDir - (String) Directory containing the Source DICOMs
        
    Returns:
        Image    - (SimpleITK object) SimpleITK image
    
    """
    
    
    
    SitkReader = sitk.ImageSeriesReader()
    SitkNames = SitkReader.GetGDCMSeriesFileNames(DicomDir)
    
    SitkReader.SetFileNames(SitkNames)
    Image = SitkReader.Execute()
    Image = sitk.Cast(Image, sitk.sitkFloat32)
    
    
    """
    Above code reads the DICOMs in a random order of filename for some reason.
    https://stackoverflow.com/questions/41037407/itk-simpleitk-dicom-series-loaded-in-wrong-order-slice-spacing-incorrect#41037408
    """
    
    
    
    """
    FileList = sitk.ImageSeriesReader().GetGDCMSeriesFileNames(DicomDir)
    
    Image = sitk.ReadImage(FileList)
    
    Image = sitk.Cast(Image, sitk.sitkFloat32)
    
    The above code didn't work. Try using GetDicomFpaths:
    """
    
    
    
    """
    SortedFilePaths = GetDicomFpaths(DicomDir, 'slices', False)

    SitkReader = sitk.ImageSeriesReader()
    SitkReader.SetFileNames(SortedFilePaths)
    
    Image = SitkReader.Execute()
    Image = sitk.Cast(Image, sitk.sitkFloat32)
    
    
    I think the reason why the images were loaded in the wrong order for a 
    particular series was that there were repeated IPPs and SimpleITK mis-
    interpretted the data (?):
    
    https://github.com/SimpleITK/SimpleITK/issues/1245
    
    So commenting out the section above which uses GetDicomFpaths to provide
    a list of sorted filepaths to SimpleITK.  It's probably fine to keep but
    possibly unnecessary and redundant.
    """
    
    return Image






def ItemsUniqueToWithin(List, epsilon=1e-5):
    
    UniqueList = []
    
    L = len(List)
    
    if L < 2:
        print(f'The list has {L} items (it must have at least 2 items).')
        
        return UniqueList
    
    else:
        UniqueList.append(List[0])
        
        for i in range(L - 1):
            if abs(List[i] - List[0]) > epsilon:
                UniqueList.append(List[i])
                
        return UniqueList
    
    
    




def GetImageAttributes(DicomDir, Package='pydicom', LogToConsole=False):
    """
    Get image size, voxel spacings, image positions (IPPs), and direction 
    cosines either using one of two packages: Pydicom or SimpleITK.
    
    
    Inputs:
        DicomDir   - (String) Directory containing DICOMs  
        
        Package    - (String) Package to use; acceptable inputs are:
                     - 'pydicom' (default)
                     - 'sitk' (SimpleITK)
        
    Returns:
        Size       - (List of integers) The size/dimensions of the 3D image 
                     along x, y and z, e.g. [Columns, Rows, NumOfSlices]
                     
        Spacings   - (List of floats) The pixel spacings along x, y and z 
                     (= SliceThickness appended to PixelSpacing), 
                     e.g. [di, dj, dk]
                     
        Positions  - (List of list of floats) The ImagePositionPatient of all 
                     slices in the DICOM series, 
                     e.g. [[x0_0, y0_0, z0_0], [x1_0, y1_0, z1_0], ...]
        
        Directions - (List of floats) The direction cosine along x (rows), 
                     y (columns) and z (slices) 
                     
                     Notes:
                         1) When using Pydicom the cross product of the x and y
                         direction cosines from ImageOrientationPatient is
                         used to obtain the direction cosine along x.
                         
                         2) There is a sign difference between the cross-terms
                         obtained from Pydicom and SimpleITK. If the Directions
                         obtained from Pydicom are:
                         
                         [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz],
                         
                         the Directions obtained from SimpleITK would be:
                         
                         [Xx, -Xy, -Xz, -Yx, Yy, -Yz, -Zx, -Zy, Zz].
                         
                         3) When using Pydicom the Directions will not be 
                         correct for coronal or sagittal image stacks.
        
    """

    if Package=='sitk':
        # Import the image:
        SitkIm = ImportImage(DicomDir)
        
        Size = SitkIm.GetSize()
        
        Spacings = SitkIm.GetSpacing()
        
        #Origin = SitkIm.GetOrigin() 
        
        Positions = []
        
        for i in range(Size[2]):
            Positions.append(SitkIm.TransformIndexToPhysicalPoint((0, 0, i)))
        
        Directions = SitkIm.GetDirection() 
        
        SliceThick = None
        
        
        
    elif Package=='pydicom':
        """
        NOTE:
            This only gives the correct Directions for an axial stack.
            The components need to be switched around the coronal or sagittal.
        """
        
        # Import the DICOMs:
        Dicoms = ImportDicoms(DicomDir=DicomDir, 
                              SortMethod='slices', 
                              LogToConsole=False)
        
        Size = [int(Dicoms[0].Columns), int(Dicoms[0].Rows), len(Dicoms)]
        
        
        #Origin = [float(item) for item in Dicoms[0].ImagePositionPatient]
        
        Positions = []
    
        for dicom in Dicoms:
            Positions.append([float(item) for item in dicom.ImagePositionPatient])
        
        
        # Check if all IPPs are unique.
        """
        Some scans ("ep2ddiff3scantracep2") have repeated IPPs.
        """
        
        # Get unique positions:
        #UniquePositions = list(set(Positions))
        #UniquePositions = np.unique(np.array(Positions))
        UniquePositions = np.array(list(set(tuple(p) for p in Positions)))
        
        P = len(Positions)
        U = len(UniquePositions)
        
        if U != P:
            print(f'\nWarning:  There are only {U} unique IPPs within the',
                  f'list of {P} IPPs.\n')
        
        
        IOP = [float(item) for item in Dicoms[0].ImageOrientationPatient]
        
        # Get the direction vector along z using the cross product of the x and 
        # y vectors:
        zDir = np.cross(IOP[0:3], IOP[3:])
        
        # Append zDir to IOP:
        Directions = IOP
        Directions.extend(zDir)
        
        
        
        Spacings = [float(item) for item in Dicoms[0].PixelSpacing]
        
        SliceThick = float(Dicoms[0].SliceThickness)
        
        """ 
        Don't use SliceThickness for the z-component of the voxel spacing. 
        Instead use the difference between slices from the IPP.
        """
        # Append SliceThick to Spacings:
        #Spacings.append(SliceThick)
        
        # Compute the slice spacings along the z-direction:
        #Dz = np.diff(np.array(Positions), axis=0)[:,2] # <-- WRONG since it doesn't account for IOP
        
        ##UniqueDz = list(set(Dz))
        #UniqueDz = ItemsUniqueToWithin(Dz)
        
        # Compute the vector lengths of the IPP differences:
        Vlengths = []

        for i in range(len(Positions) - 1):
            vector = np.array(Positions[i + 1]) - np.array(Positions[i])
            
            vectorL = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
            
            Vlengths.append(vectorL)
            
        UniqueVlengths = ItemsUniqueToWithin(Vlengths)
        
        # Append UniqueVlengths to Spacings:
        Spacings.append(UniqueVlengths[0])
        
        if len(UniqueVlengths) > 1:
            print('\nWarning:')
            print('    The voxel spacings along the scan direction are',
                  'non-uniform with the following unique values:\n')
            print('   ', UniqueVlengths, '\n')
        
        
        # Compare UniqueDz to SliceThickness:    
        #epsilon = 1e-5
        #
        #if abs(SliceThick - UniqueDz[0]) > epsilon:
        #    print('\nNote:')
        #    print('    The slice thickness obtained from the IPPs do not',
        #          'agree with the DICOM tag SliceThickness:')
        #    print(f'      d(IPP[2])      = {UniqueDz[0]}')
        #    print(f'      SliceThickness = {SliceThick}')
        
    
        
    if LogToConsole:
        print(f'Size = {Size} \n\nSpacings = {Spacings} \n\n',
              f'SliceThickness = {SliceThick} \n\n',
              f'Positions = {Positions} \n\nDirections = {Directions}')
    
    
    return Size, Spacings, SliceThick, Positions, Directions







def CompareSourceTargetImageAttributes(SrcDcmDir, TrgDcmDir):
#def CompareSrTrgImAttrs(SrcDcmDir, TrgDcmDir):    
    # Get the image attributes using Pydicom:
    SrcPydiSize, SrcPydiSpacing, SrcPydiST,\
    SrcPydiIPP, SrcPydiDir = GetImageAttributes(SrcDcmDir)
    
    TrgPydiSize, TrgPydiSpacing, TrgPydiST,\
    TrgPydiIPP, TrgPydiDir = GetImageAttributes(TrgDcmDir)
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcSitkSize, SrcSitkSpacing, SrcSitkST,\
    SrcSitkIPP, SrcSitkDir = GetImageAttributes(DicomDir=SrcDcmDir,
                                                Package='sitk')
    
    TrgSitkSize, TrgSitkSpacing, TrgSitkST,\
    TrgSitkIPP, TrgSitkDir = GetImageAttributes(DicomDir=TrgDcmDir,
                                                Package='sitk')
    
    
    title = 'Select image attributes for *Source* DICOMs from Pydicom and' \
            + 'SimpleITK:'
    ul = '*' * len(title)
    
    print('\n' + title)
    print(ul + '\n')
    print(f'Pydi Size           = {SrcPydiSize}')
    print(f'Sitk Size           = {list(SrcSitkSize)}\n')
    print(f'Pydi Spacing        = {SrcPydiSpacing}')
    print(f'Sitk Spacing        = {list(SrcSitkSpacing)}\n')
    print(f'Pydi SliceThickness = {SrcPydiST}')
    print(f'Pydi Origin         = {SrcPydiIPP[0]}')
    print(f'Sitk Origin         = {list(SrcSitkIPP[0])}\n\n')
    
    title = 'Select image attributes for *Target* DICOMs from Pydicom and SimpleITK:'
    ul = '*' * len(title)
    
    print('\n' + title)
    print(ul + '\n')
    print(f'Pydi Size           = {TrgPydiSize}')
    print(f'Sitk Size           = {list(TrgSitkSize)}\n')
    print(f'Pydi Spacing        = {TrgPydiSpacing}')
    print(f'Sitk Spacing        = {list(TrgSitkSpacing)}\n')
    print(f'Pydi SliceThickness = {TrgPydiST}')
    print(f'Pydi Origin         = {TrgPydiIPP[0]}')
    print(f'Sitk Origin         = {list(TrgSitkIPP[0])}\n\n')
    
    
    return





def GetDicomSOPuids(DicomDir):
    """
    Get a list of SOP Instance UIDs for a DICOM series.
    
    Inputs:
        DicomDir - (String) Directory containing the DICOMs
        
    Returns:
        SOPuids  - List of strings of the SOP UIDs of the DICOMs 
    """
    
    Dicoms = ImportDicoms(DicomDir=DicomDir, SortMethod='slices', 
                          LogToConsole=False)
    
    SOPuids = [] 
    
    for dicom in Dicoms:
        SOPuids.append(dicom.SOPInstanceUID)
        
    return SOPuids











def InspectDicomStudyDir(StudyDir):
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
    ScanDirNames = natsort.natsorted(ScanDirNames)
    
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
    Dims = []
    FORuids = []
    
    # Loop through each Scan directory:
    for ScanDir in ScanDirs:
    
        # Get the Image Attributes:
        origin, dirs,\
        spacings, dims = GetImageAttributes(DicomDir=ScanDir, 
                                            Package='pydicom')
        
        # Get the DICOMs:
        DicomFpaths, Dicoms = GetDicoms(DirPath=ScanDir, 
                                        SortMethod='slices', 
                                        Debug=False)
        
        FORuid = Dicoms[0].FrameOfReferenceUID
        
        Origins.append(origin)
        Dirs.append(dirs)
        Spacings.append(spacings)
        Dims.append(dims)
        FORuids.append(FORuid)
        
    
    # Create a dictionary:
    StudyDict = {'Scan num':ScanNums,
                 'Scan label':ScanLabels,
                 #'FOR UID':FORuids,
                 'Origin':Origins,
                 'Directions':Dirs,
                 'Dimensions':Dims,
                 'Spacings':Spacings
                 }
    
    #return StudyDict, ScanNums, ScanLabels, FORuids, Origins, Dirs, Dims, Spacings
    return StudyDict
    
    



def InspectDicomSubjectDir(SubjectDir):
    # Get a list of study directories:
    StudyDirs = os.list(SubjectDir)
    
    # Sort using natsort:
    StudyDirs = natsort.natsorted(StudyDirs)
    
    
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