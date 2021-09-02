# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:56:14 2020

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
DICOM IMAGE FUNCTIONS
******************************************************************************
******************************************************************************
"""



def ImportImage(DicomDir, Dtype='float32'):
    """
    06/07/21:
        See image_tools.imports.py
        
    Import a DICOM series as a SimpleITK image.
    
    Inputs:
    ******
        
    DicomDir : string
        Directory containing the DICOMs.
    
    Dtype : string (optional; 'float32' by default)
        The pixel ID type as defined at the link in the Notes section. 
        Acceptable values are:
            - 'SameAsRefImage' (default)
            - 'uint8'
            - 'int8'
            - 'uint16'
            - 'int16'
            - 'uint32'
            - 'int32'
            - 'uint64'
            - 'int64'
            - 'float32'
            - 'float64'
    
    
    Outputs:
    *******
    
    Image : SimpleITK object
        SimpleITK 3D image of the DICOM stack.
    
    
    Notes:
    *****
    
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    import SimpleITK as sitk
    
    if Dtype == 'uint8':
        PixIdType = sitk.sitkUInt8
    elif Dtype == 'int8':
        PixIdType = sitk.sitkInt8
    elif Dtype == 'uint16':
        PixIdType = sitk.sitkUInt16
    elif Dtype == 'int16':
        PixIdType = sitk.sitkInt16
    elif Dtype == 'uint32':
        PixIdType = sitk.sitkUInt32
    elif Dtype == 'int32':
        PixIdType = sitk.sitkInt32
    elif Dtype == 'uint64':
        PixIdType = sitk.sitkUInt64
    elif Dtype == 'int64':
        PixIdType = sitk.sitkInt64
    elif Dtype == 'float32':
        PixIdType = sitk.sitkFloat32
    elif Dtype == 'float64':
        PixIdType = sitk.sitkFloat64
    else:
        msg = f"Dtype '{Dtype}' is not a valid input."
        raise Exception(msg)
    
    Reader = sitk.ImageSeriesReader()
    Fnames = Reader.GetGDCMSeriesFileNames(DicomDir)
    Reader.SetFileNames(Fnames)
    #Reader.ReadImageInformation() # doesn't exist for ImageSeriesReader; works
    # for ImageFileReader
    Image = Reader.Execute()
    #Image = sitk.Cast(Image, sitk.sitkFloat32)
    Image = sitk.Cast(Image, PixIdType)
    
    
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


def ImportNifti(FilePath):
    """
    06/07/21:
        See image_tools.imports.py
        
    Import a NIFTI as a SimpleITK image.
    
    Inputs:
    ******
        
    Filepath : string
        Filepath of NIFTI file.
        
        
    Outputs:
    *******
    
    Image : SimpleITK object
    """
    
    import SimpleITK as sitk
    
    
    Reader = sitk.ImageFileReader()
    Reader.SetImageIO('NiftiImageIO')
    Reader.SetFileName(FilePath)
    Image = Reader.Execute()
    
    return Image


def ExportImageToNifti(Image, FileName, ExportDir='cwd'):
    """ 
    Export a SimpleITK image to a NIFTI file in the current working directory. 
    
    Inputs:
    ******
    
    Image : SimpleITK image 
        The image to be exported.
    
    FileName : string
        The filename to assign to the exported Nifti file. If a .nii file
        extension is not provided it will be added automatically.
    
    ExportDir : string (optional; 'cwd' by default)
        If provided the directory to which the file will be exported to. If not
        provided the file will be exported to the current working directory.
        If the directory doesn't exist it will be created.
    """
    
    import os
    from pathlib import Path
    import SimpleITK as sitk
    
    if not '.nii' in FileName:
        FileName += '.nii'
    
    if ExportDir == 'cwd':
        FilePath = FileName
    else:
        if not os.path.isdir(ExportDir):
            #os.mkdir(ExportDir)
            Path(ExportDir).mkdir(parents=True)
        
        FilePath = os.path.join(ExportDir, FileName)
    
    writer = sitk.ImageFileWriter()
    writer.SetFileName(FilePath)
    writer.Execute(Image)
    
    return


def GetImageAttributes(DicomDir, Package='pydicom', LogToConsole=False):
    """
    06/07/21:
        See image_tools.attributes_info.py
        
    Get image size, voxel spacings, image positions (IPPs), and direction 
    cosines either using one of two packages: Pydicom or SimpleITK.
    
    
    Inputs:
    ******
        
    DicomDir : string
        Directory containing DICOMs.
        
    Package : string (optional; 'pydicom' by default) 
        Package to use; acceptable inputs are:
         - 'pydicom' (default)
         - 'sitk' (SimpleITK)
        
    Outputs:
    *******
        
    Size : list of integers
        The size/dimensions of the 3D image along x, y and z, 
        e.g. [Columns, Rows, NumOfSlices]
                     
    Spacings : list of floats
        The pixel spacings along x, y and z (= SliceThickness appended to 
        PixelSpacing), e.g. [di, dj, dk]
                     
    Positions : list of list of floats
        The ImagePositionPatient of all slices in the DICOM series, 
        e.g. [[x0_0, y0_0, z0_0], [x1_0, y1_0, z1_0], ...]
        
    Directions : list of floats
        The direction cosine along x (rows), y (columns) and z (slices).
        
    ListOfWarnings : list of strings
        List of any warnings.
           
          
    Notes:
    *****
    
    1) When using Pydicom the cross product of the x and y direction cosines 
    from ImageOrientationPatient is used to obtain the direction cosine along x.
    
    2) There is a sign difference between the cross-terms obtained from Pydicom
    and SimpleITK. If the Directions obtained from Pydicom are:
        [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz],
        
    the Directions obtained from SimpleITK would be:
        [Xx, -Xy, -Xz, -Yx, Yy, -Yz, -Zx, -Zy, Zz].
        
    3) When using Pydicom the Directions will not be correct for coronal or 
    sagittal image stacks.    
    """
    
    import numpy as np
    from DicomTools import ImportDicoms
    from GeneralTools import ItemsUniqueToWithin
    
    
    ListOfWarnings = []
    
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
        
        
        """
        Since ImageSeriesReader() (invoked in ImportImage()) doesn't allow for
        calling ReadImageInformation(), it's not straightforward to get 
        metadata using sitk.  So instead get the FOR UID and SliceThickness 
        using Pydicom.
        """
        # Import the DICOMs:
        Dicoms = ImportDicoms(DicomDir=DicomDir, 
                              SortMethod='slices', 
                              LogToConsole=False)
        
        SliceThick = float(Dicoms[0].SliceThickness)
        
        
        
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
            msg = f'\nWarning:  There are only {U} unique IPPs within the '\
                  + f'list of {P} IPPs.\n'
                  
            ListOfWarnings.append(msg)
                  
            print(msg)
        
        
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
            
        
        """ Are the list of unique vector lengths equal to within 1% of the
        maximum value? """
        thresh = max(Vlengths)/100
        
        UniqueVlengths = ItemsUniqueToWithin(Vlengths, thresh)
        
        # Append UniqueVlengths to Spacings:
        Spacings.append(UniqueVlengths[0])
        
        if len(UniqueVlengths) > 1:
            #msg = f'\nWarning:  The voxel spacings along the scan direction '\
            #      + f'are non-uniform with the following unique values:\n'\
            #      + f'{UniqueVlengths}'
            
            msg = '\nWarning:  The voxel spacings along the scan direction '\
                  + f'are not to within 1% of each other:\n{UniqueVlengths}'
            
            ListOfWarnings.append(msg)
                  
            print(msg)
        
        
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
        print(f'\nSize = {Size} \nSpacings = {Spacings}',
              f'\nSliceThickness = {SliceThick} \nPositions = {Positions}',
              f'\nDirections = {Directions}')
    
    
    return Size, Spacings, SliceThick, Positions, Directions, ListOfWarnings


def GetImAttributesFromListOfDcmDirs(ListOfDicomDirs):
    """ 
    Note:  Modelled from GetSegDataFromListOfLabImBySeg.
    
    ListOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.
    """
    
    #import importlib
    #import DicomTools
    #importlib.reload(DicomTools)
    
    #import numpy as np
    from DicomTools import GetDicomFpaths
    #from DicomTools import GetDicomSOPuids
    from ImageTools import GetImageAttributes
    
    ListOfDcmFpaths = []
    ListOfIPPs = []
    ListOfDirections = []
    ListOfSpacings = []
    
    for i in range(len(ListOfDicomDirs)):
        DcmFpaths = GetDicomFpaths(ListOfDicomDirs[i])
        
        Size, Spacings, ST, IPPs, Dirs,\
        ListOfWarnings = GetImageAttributes(ListOfDicomDirs[i])
        
        ListOfDcmFpaths.append(DcmFpaths)
        #ListOfOrigins.append(IPPs[0])
        ListOfIPPs.append(IPPs)
        ListOfDirections.append(Dirs)
        ListOfSpacings.append(Spacings)
    
    return ListOfDcmFpaths, ListOfIPPs, ListOfDirections, ListOfSpacings


def GetImageAttributesFromDict(Dict, ListOfKeywords=[], Package='pydicom', 
                               RemoveFromKeys='', LogToConsole=False):
    """ Get image attributes from a list of directories in a dictionary and
    store in a new dictionary that contain the same keys.  Option to remove 
    character strings from the keys in Dict for the keys for the new 
    dictionary, e.g. if Dict has the key 'MR4_Ser4_DcmDir', the new dictionary
    will contain the key 'MR4_Ser4' if RemoveFromKeys = '_DcmDir'. Also since
    the main use of this function will be to generate a table underscores ('_') 
    will be replaced with spaces (' '), e.g. resulting in 'MR4 Ser4'. """
    
    from GeneralTools import ReduceDictToKeysMatchingKeywords
    from DicomTools import GetDicomFpaths
    
    """ Eliminate any items that aren't directories and whose key doesn't
    contain all keywords in ListOfKeywords. """
    Dict = ReduceDictToKeysMatchingKeywords(Dict, ListOfKeywords, 
                                            KeepDirectories=True,
                                            KeepFilepaths=False)
    
    AttributesDict = {}
    
    for Key, Dir in Dict.items():
        #print('\nDir =', Dir)
        
        """ Use GetDicomFpaths() to determine whether there are any DICOMs in
        Dir. If not, skip to the next item in Dict. """
        Fpaths = GetDicomFpaths(Dir)
        
        if Fpaths:
            Size, Spacings, ST, IPPs, DirCos,\
            ListOfWarnings = GetImageAttributes(Dir, Package, LogToConsole)
            
            """ Skip this key if there are any warnings. """
            if not ListOfWarnings:
                """ Create a dictionary of the attributes that will be a sub-
                dictionary of AttributesDict. """
                SubDict = {}
                
                SubDict['Size'] = Size
                SubDict['Spacings'] = [round(item, 2) for item in Spacings]
                SubDict['Slice thickness'] = ST
                ##SubDict['IPPs'] = IPPs
                SubDict['Origin'] = [round(item, 2) for item in IPPs[0]]
                DirCos = [round(item, 3) for item in DirCos]
                
                #SubDict['Size'] = [str(item) for item in Size]
                #SubDict['Spacings'] = [str(round(item, 3)) for item in Spacings]
                #SubDict['Slice thickness'] = str(ST)
                #SubDict['Origin'] = [str(round(item, 1)) for item in IPPs[0]]
                #DirCos = [str(round(item, 3)) for item in DirCos]
                
                """ Reformat the direction cosines so that each direction 
                cosine is on a new line (makes it easier to read): """
                X = [DirCos[i] for i in range(0,3)]
                Y = [DirCos[i] for i in range(3,6)]
                Z = [DirCos[i] for i in range(6,9)]
                
                #SubDict['Direction cosines'] = f'[{X}, \n{Y}, \n{Z}]'
                SubDict['Direction cosines'] = f'[{X[0]}, {X[1]}, {X[2]},\n'\
                                               + f'{Y[0]}, {Y[1]}, {Y[2]},\n'\
                                               + f'{Z[0]}, {Z[1]}, {Z[2]}]'
                
                if RemoveFromKeys:
                    Key = Key.replace(RemoveFromKeys, '').replace('_', ' ')
                
                AttributesDict[Key] = SubDict
    
    return Dict, AttributesDict


def GetImageInfo(Image, LogToConsole=False):
    """
    06/07/21:
        See image_tools.attributes_info.py
    """
    import numpy as np
    from ConversionTools import Image2PixArr
    #import importlib
    #import GeneralTools
    #importlib.reload(GeneralTools)
    from GeneralTools import UniqueItems
    
    PixArr, F2Sinds = Image2PixArr(Image)
    
    #print(f'\nPixArr.shape = {PixArr.shape}')
    #print(f'F2Sinds = {F2Sinds}')
    #print(f'\nPixArr = {PixArr}')
    #print(f'\ntype(PixArr) = {type(PixArr)}')
    
    UniqueVals = UniqueItems(Items=PixArr, IgnoreZero=False)
    
    #print(f'UniqueVals = {UniqueVals}')
    #print(f'type(UniqueVals) = {type(UniqueVals)}')
    
    PixID = Image.GetPixelID()
    PixIDTypeAsStr = Image.GetPixelIDTypeAsString()
    
    if LogToConsole:
        #print('\nImage info:')
        #print(f'type(Image) = {type(Image)}')
        print(f'    Image PixID = {PixID}')
        print(f'    Image PixIDTypeAsString = {PixIDTypeAsStr}')
        print(f'    Image Size = {Image.GetSize()}')
        print(f'    Image Max = {ImageMax(Image)}')
        print(f'    Image Min = {ImageMin(Image)}')
        print('\n    Conversion of Image to PixArr:')
        print(f'    PixArr shape = {PixArr.shape}')
        if isinstance(UniqueVals, np.ndarray):
            if len(UniqueVals) < 7:
                print('\n    There are {len(UniqueVals)} unique values in PixArr:')
                print(f'    {UniqueVals}')
            else:
                
                print(f'\n    There are {len(UniqueVals)} unique values in PixArr:')
                print(f'    {UniqueVals[:3]}...{UniqueVals[-3:-1]}')
            
            """ Print histogram of UniqueVals
            https://www.machinelearningplus.com/plots/matplotlib-histogram-python-examples/
            """
            
            Vals = PixArr.flatten()
            MinVal = min(Vals)
            MaxVal = max(Vals)
            
            Freq, Bins = np.histogram(Vals, bins=10, range=[MinVal, MaxVal])
            
            print('\n      Distribution of values:')
            #for b, f in zip(Bins[1:], Freq):
            #    #print(round(b, 1), ' '.join(np.repeat('*', f)))
            #    print(f'      {round(b, 2)} - {f}')
            for i in reversed(range(len(Freq))):
                valRange = f'[{round(Bins[i], 2)} - {round(Bins[i+1], 2)}]'
                print(f'      {valRange} : {Freq[i]}')
            
            
            #Vals = PixArr.flatten()
            #MinVal = 0
            ##MaxVal = 0.2
            #MaxVal - 0.1
            #Freq, Bins = np.histogram(Vals, bins=10, range=[MinVal, MaxVal])
            #print('\n      Distribution of values near 0:')
            #for i in reversed(range(len(Freq))):
            #    valRange = f'[{round(Bins[i], 2)} - {round(Bins[i+1], 2)}]'
            #    print(f'      {valRange} : {Freq[i]}')
            
            
        elif UniqueVals == None:
            print('       There are no UniqueVals (= {UniqueVals}) in PixArr:')
        print(f'       \nThere are {len(F2Sinds)} frames with slice indices:')
        print(f'       {F2Sinds}')
    
    return PixID, PixIDTypeAsStr, UniqueVals, F2Sinds


def ChangeImagePixelType(Image, NewPixelType):
    """
    07/07/21:
        See change_im_dtype in image_tools.operations.py
    
    Convert a 3D SimpleITK image to a new pixel type.
    
    Inputs:
    ******                      
        
    Image : SimpleITK image
        The 3D image whose pixel type is to be converted.
        
    NewPixelType : string
        The pixel type to convert Image to.  Acceptable inputs are a subset of
        the SimpleITK pixel types 
        (https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html):
         
        - 'UInt8' 	Unsigned 8-bit integer
        - 'UInt16' 	Unsigned 16-bit integer
        - 'UInt32' 	Unsigned 32-bit integer
        - 'UInt64' 	Unsigned 64-bit integer
        - 'Int8' 	Signed 8-bit integer
        - 'Int16' 	Signed 16-bit integer
        - 'Int32' 	Signed 32-bit integer
        - 'Int64' 	Signed 64-bit integer
        - 'Float32' 	32-bit float
        - 'Float64' 	64-bit float
    
    
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
    
    import SimpleITK as sitk
    
    if NewPixelType == 'UInt8':
        OutputPixelType = sitk.sitkUInt8
    
    elif NewPixelType == 'UInt16':
        OutputPixelType = sitk.sitkUInt16
    
    elif NewPixelType == 'UInt32':
        OutputPixelType = sitk.sitkUInt32
    
    elif NewPixelType == 'UInt64':
        OutputPixelType = sitk.sitkUInt64
    
    elif NewPixelType == 'Int8':
        OutputPixelType = sitk.sitkInt8
    
    elif NewPixelType == 'Int16':
        OutputPixelType = sitk.sitkInt16
    
    elif NewPixelType == 'Int32':
        OutputPixelType = sitk.sitkInt32
    
    elif NewPixelType == 'Int64':
        OutputPixelType = sitk.sitkInt64
        
    elif NewPixelType == 'Float32':
        OutputPixelType = sitk.sitkFloat32
    
    elif NewPixelType == 'Float64':
        OutputPixelType = sitk.sitkFloat64
    
    else:
        msg = "'NewPixelType' must be either:\n 'UInt8'\n 'UInt16'\n 'UInt32'"\
              + "\n 'UInt64'\n 'Int8'\n 'Int16'\n 'Int32'\n 'Int64'"\
              + "\n 'Float32'\n 'Float64'."
        
        raise Exception(msg)
        
    
    Caster = sitk.CastImageFilter()
    
    Caster.SetOutputPixelType(OutputPixelType)
    
    Image = Caster.Execute(Image)
    
    return Image


def GetImageVertices(Image):
    """
    07/07/21:
        See get_im_verts in general_tools.geometry.py
        
    Get the vertices that define the 3D extent of a 3D image.
    
    Inputs:
    ******
    
    Image : SimpleITK image
        A 3D image whose physical extent is to be determined.
    
        
    Outputs:
    *******
    
    pts : list of floats
        A list of length 8 of floats for all vertices that define the 3D extent
        of Image.
    """
    
    W = Image.GetWidth()
    H = Image.GetHeight()
    D = Image.GetDepth()
    
    pts = [Image.TransformIndexToPhysicalPoint((0, 0, 0)), 
           Image.TransformIndexToPhysicalPoint((W, 0, 0)),
           Image.TransformIndexToPhysicalPoint((W, H, 0)),
           Image.TransformIndexToPhysicalPoint((0, H, 0)),
           Image.TransformIndexToPhysicalPoint((0, H, D)),
           Image.TransformIndexToPhysicalPoint((W, H, D)),
           Image.TransformIndexToPhysicalPoint((W, 0, D)),
           Image.TransformIndexToPhysicalPoint((0, 0, D))
           ]
    
    return pts


def GetImageExtent(Image):
    """
    07/07/21:
        See get_im_extent in general_tools.geometry.py
        
    Get the extent of a 3D image.
    
    Inputs:
    ******
    
    Image : SimpleITK image
        A 3D image whose physical extent is to be determined.
    
        
    Outputs:
    *******
    
    Extent : list of floats
        A list (for each dimension) of the max-to-min span along the direction.
    """
    
    Sizes = Image.GetSize()
    Spacings = Image.GetSpacing()
    
    Extent = [Sizes[i]*Spacings[i] for i in range(len(Sizes))]
    
    return Extent


def GetSliceVertices(Image, SliceNum):
    """
    07/07/21:
        See get_im_extent in general_tools.geometry.py
        
    Get the vertices that define the 3D extent of a slice in a 3D image. 
    
    Inputs:
    ******
    
    Image : SimpleITK image
        A 3D image.
        
    SliceNum : integer
        The slice index of interest.
    
        
    Outputs:
    *******
    
    pts : list of floats
        A list of length 8 of floats for all vertices that define the 3D extent
        of the slice.
    """
    
    W = Image.GetWidth()
    H = Image.GetHeight()
    #D = Image.GetDepth()
    
    pts = [Image.TransformContinuousIndexToPhysicalPoint((0, 0, SliceNum - 0.5)), 
           Image.TransformContinuousIndexToPhysicalPoint((W, 0, SliceNum - 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((W, H, SliceNum - 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((0, H, SliceNum - 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((0, H, SliceNum + 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((W, H, SliceNum + 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((W, 0, SliceNum + 0.5)),
           Image.TransformContinuousIndexToPhysicalPoint((0, 0, SliceNum + 0.5))
           ]
    
    return pts


def CheckImageSpacing(DicomDir):
    
    Size, Spacings, ST, IPPs, Dirs,\
    ListOfWarnings = GetImageAttributes(DicomDir, 'sitk')
    
    Image = ImportImage(DicomDir)
    
    vertices = GetImageVertices(Image)
    
    xMin = min([vertices[i][0] for i in range(len(vertices))])
    xMax = max([vertices[i][0] for i in range(len(vertices))])
    yMin = min([vertices[i][1] for i in range(len(vertices))])
    yMax = max([vertices[i][1] for i in range(len(vertices))])
    zMin = min([vertices[i][2] for i in range(len(vertices))])
    zMax = max([vertices[i][2] for i in range(len(vertices))])
    
    W = xMax - xMin
    H = yMax - yMin
    D = zMax - zMin
    
    dx = W/Image.GetWidth()
    dy = H/Image.GetHeight()
    dz = D/Image.GetDepth()
    
    Spacings_calc = [dx, dy, dz]
    
    print(f'Image size      = {Image.GetSize()}')
    print(f'Slice thickness = {ST}\n')
    print(f'Voxel spacings from GetSpacing()        = {Image.GetSpacing()}')
    print(f'Voxel spacings calculated from vertices = {Spacings_calc}')
    diff = [Spacings_calc[i] - Spacings[i] for i in range(3)]
    print(f'Difference in voxel spacings            = {diff}')
    
    return


def CompareImageAttributes(SrcDcmDir, TrgDcmDir):
#def CompareSrTrgImAttrs(SrcDcmDir, TrgDcmDir):    
    # Get the image attributes using Pydicom:
    SrcPydiSize, SrcPydiSpacing, SrcPydiST, SrcPydiIPP, SrcPydiDir,\
    ListOfWarnings = GetImageAttributes(SrcDcmDir)
    
    TrgPydiSize, TrgPydiSpacing, TrgPydiST, TrgPydiIPP, TrgPydiDir,\
    ListOfWarnings = GetImageAttributes(TrgDcmDir)
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcSitkSize, SrcSitkSpacing, SrcSitkST, SrcSitkIPP, SrcSitkDir,\
    ListOfWarnings = GetImageAttributes(DicomDir=SrcDcmDir,
                                                Package='sitk')
    
    TrgSitkSize, TrgSitkSpacing, TrgSitkST, TrgSitkIPP, TrgSitkDir,\
    ListOfWarnings = GetImageAttributes(DicomDir=TrgDcmDir,
                                                Package='sitk')
    
    
    title = 'Image attributes for *Source* DICOMs from Pydicom and SimpleITK:'
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
    
    title = 'Image attributes for *Target* DICOMs from Pydicom and SimpleITK:'
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


def InitialiseImage(RefImage, Dtype='SameAsRefImage'):
    """
    06/07/21:
        See initialise_im in image_tools.operations.py
        
    Initialise an empty SimpleITK image with the same size, spacing, origin and
    direction as another image.  
    
    Inputs:
    ******
    
    RefImage : SimpleITK image
        The reference image whose attributes will be used for the initialised
        image.
    
    Dtype : string (optional; 'SameAsRefImage' by default)
        The pixel ID type as defined at the link in the Notes section. 
        Acceptable values are:
            - 'SameAsRefImage' (default)
            - 'uint8'
            - 'int8'
            - 'uint16'
            - 'int16'
            - 'uint32'
            - 'int32'
            - 'uint64'
            - 'int64'
            - 'float32'
            - 'float64'
        
        
    Outputs:
    *******
    
    NewImage : SimpleITK image
        An empty image with the same PixelID, Size, Spacing, Origin and 
        Direction as RefImage.
        
        
    Notes:
    *****
    
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    import SimpleITK as sitk
    
    if Dtype == 'SameAsRefImage':
        PixIdType = RefImage.GetPixelID()
    elif Dtype == 'uint8':
        PixIdType = sitk.sitkUInt8
    elif Dtype == 'int8':
        PixIdType = sitk.sitkInt8
    elif Dtype == 'uint16':
        PixIdType = sitk.sitkUInt16
    elif Dtype == 'int16':
        PixIdType = sitk.sitkInt16
    elif Dtype == 'uint32':
        PixIdType = sitk.sitkUInt32
    elif Dtype == 'int32':
        PixIdType = sitk.sitkInt32
    elif Dtype == 'uint64':
        PixIdType = sitk.sitkUInt64
    elif Dtype == 'int64':
        PixIdType = sitk.sitkInt64
    elif Dtype == 'float32':
        PixIdType = sitk.sitkFloat32
    elif Dtype == 'float64':
        PixIdType = sitk.sitkFloat64
    else:
        msg = f"Dtype '{Dtype}' is not a valid input."
        raise Exception(msg)
    
    #NewImage = sitk.Image(RefImage.GetSize(), RefImage.GetPixelID())
    NewImage = sitk.Image(RefImage.GetSize(), PixIdType)
    
    NewImage.SetOrigin(RefImage.GetOrigin())
    
    NewImage.SetDirection(RefImage.GetDirection())
    
    NewImage.SetSpacing(RefImage.GetSpacing())
    
    return NewImage


def CopyImage(Im0):
    """
    06/07/21:
        See copy_im in image_tools.operations.py
        
    Copy a SimpleITK image.  
    
    Inputs:
    ******
    
    Im0 : SimpleITK image
        The image to be copied.
        
        
    Outputs:
    *******
    
    Im1 : SimpleITK image
        The copied image.
        
        
    Notes:
    *****
    
    https://simpleitk.org/SimpleITK-Notebooks/01_Image_Basics.html
    """
    
    import SimpleITK as sitk
    
    Im1 = sitk.Image(Im0.GetSize(), Im0.GetPixelID())
    
    #Im1.SetOrigin(Im0.GetOrigin())
    #Im1.SetDirection(Im0.GetDirection())
    #Im1.SetSpacing(Im0.GetSpacing())
    
    Im1 = sitk.Paste(Im0, Im1)
    
    return Im1


def ImageMax(Image):
    """
    06/07/21:
        See im_max in image_tools.operations.py
        
    Find maximum of an SimpleITK image.  
    
    Inputs:
    ******
        
    Image : SimpleITK image
        
        
    Outputs:
    *******
        
    Maximum : float or int
        Maximum value of Image.
    """
    
    import SimpleITK as sitk
    
    MaxImFilt = sitk.MinimumMaximumImageFilter()
    MaxImFilt.Execute(Image)
    
    return MaxImFilt.GetMaximum()


def ImageMin(Image):
    """
    06/07/21:
        See im_min in image_tools.operations.py
        
    Find minimum of an SimpleITK image.  
    
    Inputs:
    ******
        
    Image : SimpleITK image
    
        
    Outputs:
    *******
    
    Minimum : float or int
        Maximum value of Image.
    """
    
    import SimpleITK as sitk
    
    MaxImFilt = sitk.MinimumMaximumImageFilter()
    MaxImFilt.Execute(Image)
    
    return MaxImFilt.GetMinimum()


def ImageMaxByFrame(Image):
    """
    06/07/21:
        See im_max_by_frame in image_tools.operations.py
        
    Return a list of the maximum value of each frame in a SimpleITK image .  
    
    Inputs:
    ******
        
    Image : SimpleITK image
        
        
    Outputs:
    *******
        
    ListOfMax : list of floats or integers
        A list of the maximum value of each frame in Image.
    """
    
    F = Image.GetSize()[2]
    
    ListOfMax = []
    
    for f in range(F):
        ListOfMax.append(ImageMax(Image[:, :, f]))
    
    
    return ListOfMax


def OrImages(Image0, Image1):
    """
    06/07/21:
        See im_or in image_tools.operations.py
        
    Perform pixel-wise OR operation on two SimpleITK images.  
    
    Inputs:
    ******
    
    Image0 : SimpleITK image
        
    Image1 : SimpleITK image
        
    
    Outputs:
    *******
    
    ImageOr : SimpleITK image
        The pixel-wise OR of Image0 and Image1.
    """
    
    import SimpleITK as sitk
            
    OrImageFilt = sitk.OrImageFilter()
    ImageOr = OrImageFilt.Execute(Image0, Image1)
    
    return ImageOr


def AddImages(Image0, Image1):
    """
    06/07/21:
        See im_add in image_tools.operations.py
        
    Perform pixel-wise addition of two SimpleITK images.  
    
    Inputs:
    ******
    
    Image0 : SimpleITK image
        
    Image1 : SimpleITK image
    
    
    Outputs:
    *******
    
    ImageSum : SimpleITK image
        The pixel-wise sum of Image0 and Image1.
    """
    
    import SimpleITK as sitk
    
    AddImFilt = sitk.AddImageFilter()
    ImageSum = AddImFilt.Execute(Image0, Image1)
    
    return ImageSum


def SumAllLabIms(LabIms):
    """
    06/07/21:
        See pixarr_sum_of_ims in image_tools.operations.py
    """
    
    import SimpleITK as sitk
    
    # Initialise labelmap sum:
    LMImsSum = InitialiseImage(LabIms[0])
    
    for i in range(len(LabIms)):
        LMImsSum = AddImages(LMImsSum, LabIms[i])

    # Convert to numpy arrays:
    LMImsSumNpa = sitk.GetArrayFromImage(LMImsSum)
    
    return LMImsSumNpa


def GetPixelAreaRatio(Im0, Im1):
    """
    06/07/21:
        See get_pixel_area_ratio in image_tools.operations.py
        
    Get the pixel area ratio of SimpleITK images Im1 to Im0.
    
    Inputs:
    ******                      
        
    Im0, Im1 : SimpleITK images
    
    
    Outputs:
    *******
    
    Ratio : float
        The ratio of the pixel area of Im1 to Im0.
    """
    
    import numpy as np
    
    Area0 = np.prod(np.array(Im0.GetSpacing())[0:2])
    
    Area1 = np.prod(np.array(Im1.GetSpacing())[0:2])
    
    return Area1/Area0


def GetVoxelVolumeRatio(Im0, Im1):
    """
    06/07/21:
        See get_voxel_volume_ratio in image_tools.operations.py
        
    Get the voxel volume ratio of SimpleITK images Im1 to Im0.
    
    Inputs:
    ******                      
        
    Im0, Im1 : SimpleITK images
    
    
    Outputs:
    *******
    
    Ratio : float
        The ratio of the voxel volume of Im1 to Im0.
    """
    
    import numpy as np
    
    Vol0 = np.prod(np.array(Im0.GetSpacing()))
    
    Vol1 = np.prod(np.array(Im1.GetSpacing()))
    
    return Vol1/Vol0


def FindSuitableThresh(BinaryIm, NonBinaryIm, LogToConsole=False):
    """
    07/07/21:
        See find_thresh in image_tools.operations.py
    
    Find a suitable threshold level that if used to binary threshold a non-
    binary 3D SimpleITK image, NonBinaryIm, would result in a similar 
    distribution of zeros and ones in the binarised image as that of the
    original binary image, BinaryIm.
    
    Inputs:
    ******
    
    BinaryIm : SimpleITK image 
        The binary image.
        
    NonBinaryIm : SimpleITK image
        The non-binary image.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
        
    Outputs:
    *******
    
    Thresh : float
        A threshold level that when applied to NonBinaryIm should result in 
        a similar distribution of zeros and ones as in BinaryIm.
    """
    
    import numpy as np
    from ConversionTools import Image2PixArr
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of FindSuitableThresh():')
        print('\n\n', '-'*120)
    
    
    Ratio = GetVoxelVolumeRatio(BinaryIm, NonBinaryIm)
    
    
    BinPixArr, BinF2Sinds = Image2PixArr(BinaryIm)
    NonBinPixArr, NonBinF2Sinds = Image2PixArr(NonBinaryIm)
    
    NumOfOnes_B = np.count_nonzero(BinPixArr)
            
    Vals_NB = NonBinPixArr.flatten()
    MinVal = min(Vals_NB)
    MaxVal = max(Vals_NB)
    
    Freq, Bins = np.histogram(Vals_NB, bins=1000, range=[MinVal, MaxVal])
    
    if LogToConsole:
        print('\n      Distribution of values:')
        #for b, f in zip(Bins[1:], Freq):
        #    print(f'      {round(b, 8)} - {f}')
        for i in reversed(range(len(Freq))):
                valRange = f'[{round(Bins[i], 3)} - {round(Bins[i+1], 3)}]'
                print(f'      {valRange} : {Freq[i]}')
    
    
    """ The cumulative sum of the frequencies in reverse order: """
    CumSumFreq = np.cumsum(np.flip(Freq))
    
    """ The absolute differences between the cumulative sum of the frequencies
    and the number of ones in the binary image scaled by the ratio of the
    voxel volumes: """
    #AbsDiffs = np.abs(CumSumFreq - NumOfOnes_B)
    AbsDiffs = np.abs(CumSumFreq - NumOfOnes_B/Ratio)
    
    RevInd = np.where(AbsDiffs == np.min(AbsDiffs))[0][0]
    
    Ind = len(CumSumFreq) - RevInd - 1
    
    Thresh = Bins[Ind]
    
    if LogToConsole:
        print(f'\nBinaryIm.GetSpacing() = {BinaryIm.GetSpacing()}')
        print(f'\nNonBinaryIm.GetSpacing() = {NonBinaryIm.GetSpacing()}')
        print(f'\nVoxelVolumeRatio = {Ratio}')
        #print(f'\nTarget number of ones = {NumOfOnes_B}')
        print(f'\nTarget number of ones = {NumOfOnes_B/Ratio}')
        print(f'\nCumSumFreq[0:RevInd+5] = {CumSumFreq[:RevInd+5]}')
        print(f'\nAbsDiffs[0:RevInd+5] = {AbsDiffs[:RevInd+5]}')
        print(f'\nInd = {Ind}')
        print(f'Bins[{Ind}] = {Bins[Ind]}')
    
    BinarisedPixArr = (NonBinPixArr > Thresh).astype(np.int_)
    
    NumOfOnes_BNB = np.count_nonzero(BinarisedPixArr)
    
    if LogToConsole:
        print(f'There are {NumOfOnes_B} ones in the binary pixel array. ',
              f'\nThresholding the non-binary pixel array at {round(Thresh,3)}',
              f'would result in {NumOfOnes_BNB} ones (difference of',
              f'{abs(NumOfOnes_BNB - NumOfOnes_B)}).')
        print('-'*120)
    
    return Thresh


def BinaryThresholdImage(Im, Thresh=0.5):
    """
    07/07/21:
        See binarise_im in image_tools.operations.py
    
    Binary threshold a 3D SimpleITK image.
    
    Inputs:
    ******
    
    Im : SimpleITK image 
        The 3D image to be transformed.
        
    Thresh : float (optional; 0.5 by default)
        Lower limit threshold.
        
        
    Outputs:
    *******
    
    BinaryIm : SimpleITK image
        The 3D binary threshold image.
        
        
    Notes:
    *****
    
    If an empty image (all zeros) is binary thresholded, the resuling image
    has ones everywhere for some unknown reason.
    """
    
    import SimpleITK as sitk
    
    BinThreshImFilt = sitk.BinaryThresholdImageFilter()
    
    #BinThreshImFilt.SetInput(Im)
    
    """
    Setting UpperThreshold (e.g. to 1) removes all pixels that exceed 1, and
    for some reason, transformed binary labelmap images result in non-binary 
    values (including negative values).
    """
    
    BinThreshImFilt.SetLowerThreshold(Thresh)
    #BinThreshImFilt.SetUpperThreshold(1)
    BinThreshImFilt.SetOutsideValue(0)
    BinThreshImFilt.SetInsideValue(1)
    #BinThreshImFilt.SetOutputPixelType(sitk.sitkUInt32) # --> error:
    # 'BinaryThresholdImageFilter' object has no attribute 'SetOutputPixelType'
    #BinThreshImFilt.DebugOn()
    
    #BinThreshImFilt.Execute()
    #BinIm = BinThreshImFilt.GetOutput()
    
    BinIm = BinThreshImFilt.Execute(Im)
    
    #sitk.Show(Im)
    #sitk.Show(BinIm)
    
    return BinIm


def GaussianBlurImage(Im, Variance=(1,1,1), LogToConsole=False):
    """
    07/07/21:
        See gaussian_blur_im in image_tools.operations.py
    
    Gaussian blur a 3D SimpleITK image.
    
    Inputs:
    ******
    
    Im : SimpleITK image 
        The 3D image to be blurred.
        
    Variance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
        
    Outputs:
    *******
    
    BlurredIm : SimpleITK image
        The 3D blurred image.
        
        
    Notes:
    *****
    
    DiscreteGaussianImageFilter:
        Blurs an image by separable convolution with discrete gaussian kernels. 
        This filter performs Gaussian blurring by separable convolution of an 
        image and a discrete Gaussian operator (kernel).
        
    SmoothingRecursiveGaussianImageFilter:
        Computes the smoothing of an image by convolution with the Gaussian 
        kernels implemented as IIR filters
        
    The input image type must be float - otherwise odd results will arise. So
    check will be made prior to Gaussian blurring and the image type will be
    converted to float if applicable.
    """
    
    import SimpleITK as sitk
    
    if LogToConsole:
        print('\n   Image info for Im prior to Gaussian blurring:')
            
    PixID, PixIDTypeAsStr, UniqueValsPostTx,\
    F2Sinds = GetImageInfo(Im, LogToConsole)
    
    
    """ Ensure Im is a 32-bit float (PixID = 8). """
    if PixID != 8: 
        if LogToConsole:
            print(f'\n   Im has PixelID = {PixID} ({PixIDTypeAsStr})). '\
                  'It will be converted to Float32.')
        
        """ Convert Im to Float32: """
        Im = ChangeImagePixelType(Image=Im, NewPixelType='Float32')
        
        if LogToConsole:
            print('\n   Image info for Im after converting to Float32:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(Im, LogToConsole)
    
    
    ImFilt = sitk.DiscreteGaussianImageFilter()
    
    #ImFilt.SetMaximumKernelWidth(Sigma)
    #print(f'   ImFilt.GetMaximumKernelWidth() = {ImFilt.GetMaximumKernelWidth()}')
    
    ImFilt.SetVariance(Variance)
        
    BlurredIm = ImFilt.Execute(Im)
    
    if LogToConsole:
        print('\nApplying Gaussian blur...')
        print(f'   ImFilt.GetMaximumError() = {ImFilt.GetMaximumError()}')
        print(f'   ImFilt.GetMaximumKernelWidth() = {ImFilt.GetMaximumKernelWidth()}')
        print(f'   ImFilt.GetUseImageSpacing() = {ImFilt.GetUseImageSpacing()}')
        print(f'   ImFilt.GetVariance() = {ImFilt.GetVariance()}')
        print(f'\n   Im.GetPixelID() = {Im.GetPixelID()}')
        print(f'   Im.GetPixelIDValue() = {Im.GetPixelIDValue()}')
        print(f'   Im.GetPixelIDTypeAsString() = {Im.GetPixelIDTypeAsString()}')
        print(f'   BlurredIm.GetPixelID() = {BlurredIm.GetPixelID()}')
        print(f'   BlurredIm.GetPixelIDValue() = {BlurredIm.GetPixelIDValue()}')
        print(f'   BlurredIm.GetPixelIDTypeAsString() = {BlurredIm.GetPixelIDTypeAsString()}')
    
    return BlurredIm


def RecursiveGaussianBlurImage(Im, Sigma, Direction, LogToConsole=False):
    """
    07/07/21:
        See recursive_gaussian_blur_im in image_tools.operations.py
    
    Apply a recursive Gaussian blur to a 3D SimpleITK image.
    
    Inputs:
    ******
    
    Im : SimpleITK image 
        The 3D image to be blurred.
        
    Sigma : integer
        The kernel width.
        
    Direction : integer
        The direction to apply blurring to. Direction can be 0, 1, or 2 for the
        x, y or z direction.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
        
    Outputs:
    *******
    
    BlurredIm : SimpleITK image
        The 3D blurred image.
        
        
    Notes:
    *****
    
    http://itk-users.7.n7.nabble.com/ITK-users-sitk-DiscreteGaussianImageFilter-with-different-variances-per-direction-td38185.html
    """
    
    import SimpleITK as sitk
    
    # TypeError: __init__() got an unexpected keyword argument 'sigma':
    #BlurredIm = sitk.RecursiveGaussianImageFilter(Im, sigma=Sigma, 
    #                                              direction=Direction)
    
    ImFilt = sitk.RecursiveGaussianImageFilter()
    ImFilt.SetSigma(Sigma)
    ImFilt.SetDirection(Direction)
    
    if LogToConsole:
        print('\nApplying Gaussian blur...')
        print(f'   ImFilt.GetDirection() = {ImFilt.GetDirection()}')
        print(f'   ImFilt.GetNormalizeAcrossScale() = {ImFilt.GetNormalizeAcrossScale()}')
        print(f'   ImFilt.GetOrder() = {ImFilt.GetOrder()}')
        print(f'   ImFilt.GetSigma() = {ImFilt.GetSigma()}')
    
    BlurredIm = ImFilt.Execute(Im)
    
    if LogToConsole:
        print(f'\n   Im.GetPixelID() = {Im.GetPixelID()}')
        print(f'   Im.GetPixelIDValue() = {Im.GetPixelIDValue()}')
        print(f'   Im.GetPixelIDTypeAsString() = {Im.GetPixelIDTypeAsString()}')
        print(f'   BlurredIm.GetPixelID() = {BlurredIm.GetPixelID()}')
        print(f'   BlurredIm.GetPixelIDValue() = {BlurredIm.GetPixelIDValue()}')
        print(f'   BlurredIm.GetPixelIDTypeAsString() = {BlurredIm.GetPixelIDTypeAsString()}')
    
    return BlurredIm


def ResampleImage(Im, RefIm, Tx='Identity', Interp='Linear', 
                  LogToConsole=False):
    """
    07/07/21:
        See resample_im in image_tools.resampling.py
    
    Resample a 3D SimpleITK image.  The image can be a 3D DICOM image or a
    3D labelmap image.
    
    Inputs:
    ******                      
        
    Im : SimpleITK Image
        The 3D image to be resampled.
                           
    RefIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
        
    Tx : string (optional; 'Identity' by default)
        The transform to be used.  Acceptable inputs are:
        - 'Identity'
        - 'Affine'
    
    Interp : string (optional; 'Linear' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    ResIm : SimpleITK image
        The resampled 3D image.
    
    
    Notes:
    *****
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """

    import SimpleITK as sitk
    
    # Define which transform to use:
    if Tx == 'Identity':    
        sitkTx = sitk.Transform()
        
    elif Interp == 'Affine':
        dim = Im.GetDimension()
        sitkTx = sitk.AffineTransform(dim)
        
    else:
        msg = '"Tx" must be "Identity" or "Affine".'
        
        raise Exception(msg)
    
    # Define which interpolator to use:
    if Interp == 'NearestNeighbor':    
        sitkInterp = sitk.sitkNearestNeighbor
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'LabelGaussian':
        sitkInterp = sitk.sitkLabelGaussian
        #sitkPixType = sitk.sitkUInt64
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Linear':
        sitkInterp = sitk.sitkLinear
        #sitkPixType = sitk.sitkFloat32
        """11/06/21 Should sitkPixType be sitk.sitkUInt32 instead?.. """
        sitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Bspline':
        sitkInterp = sitk.sitkBSpline
        sitkPixType = sitk.sitkFloat32
        
    else:
        msg = '"Interp" must be "Linear", "BSpline" or "LabelGaussian".'
        
        raise Exception(msg)
            
    #print('\nUsing', Interpolation, 'interp\n')
    
    Resampler = sitk.ResampleImageFilter()
    
    Resampler.SetReferenceImage(RefIm)
    
    Resampler.SetInterpolator(sitkInterp)
    
    #Resampler.SetTransform(sitk.Transform())
    Resampler.SetTransform(sitkTx)
    
    #print(f'\nImage.GetPixelIDValue() = {Image.GetPixelIDValue()}')
    
    Resampler.SetOutputSpacing(RefIm.GetSpacing())
    Resampler.SetSize(RefIm.GetSize())
    Resampler.SetOutputDirection(RefIm.GetDirection())
    Resampler.SetOutputOrigin(RefIm.GetOrigin())
    #Resampler.SetDefaultPixelValue(Image.GetPixelIDValue())
    Resampler.SetOutputPixelType(sitkPixType)
    Resampler.SetDefaultPixelValue(0)
    
    #Resampler.LogToConsoleOn() # 'ResampleImageFilter' object has no attribute
    # 'LogToConsoleOn'
    #Resampler.LogToFileOn() # 'ResampleImageFilter' object has no attribute 
    # 'LogToFileOn'
    
    #print('\n')
    #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*Resampler.GetProgress()),end=''))
    #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
    
    ResIm = Resampler.Execute(Im)
    
    #sitk.Show(ResImage)
    
    #ResTx = Resampler.GetTransform()
    
    #if LogToConsole:
    #    print('\nResampling transform:\n', ResTx)

    return ResIm
    #return ResIm, ResTx


def ResampleImage_v2(Im, RefIm, SitkTx, Interp='Linear'):
    """
    07/07/21:
        See resample_im_v2 in image_tools.resampling.py
    
    Resample a 3D SimpleITK image.  The image can be a 3D DICOM image or a
    3D labelmap image.
    
    Inputs:
    ******                      
        
    Im : SimpleITK Image
        The 3D image to be resampled.
                           
    RefIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
        
    SitkTx : SimpleITK Transform
        The SimpleIK Transform to be used.
    
    Interp : string (optional; 'linear' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    ResIm : SimpleITK image
        The resampled 3D image.

    
    Notes:
    *****
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """
    
    import SimpleITK as sitk
    
    if Interp == 'NearestNeighbor':    
        SitkInterp = sitk.sitkNearestNeighbor
        #SitkPixType = sitk.sitkUInt64
        SitkPixType = sitk.sitkUInt32
        
    elif Interp == 'LabelGaussian':
        SitkInterp = sitk.sitkLabelGaussian
        #SitkInterp = sitk.sitkUInt64
        SitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Linear':
        SitkInterp = sitk.sitkLinear
        #SitkPixType = sitk.sitkFloat32
        """11/06/21 Should sitkPixType be sitk.sitkUInt32 instead?.. """
        SitkPixType = sitk.sitkUInt32
        
    elif Interp == 'Bspline':
        SitkInterp = sitk.sitkBSpline
        SitkPixType = sitk.sitkFloat32
        
    else:
        msg = '"Interp" must be "Linear", "BSpline" or "LabelGaussian".'
        
        raise Exception(msg)
    
    ResImFilt = sitk.ResampleImageFilter()
    #ResImFilt.SetTransform(sitk.Transform())
    ResImFilt.SetTransform(SitkTx)
    ResImFilt.SetInterpolator(SitkInterp)
    
    ResImFilt.SetReferenceImage(RefIm)
    #ResImFilt.SetOutputSpacing(RefIm.GetSpacing())
    #ResImFilt.SetSize(RefIm.GetSize())
    #ResImFilt.SetOutputDirection(RefIm.GetDirection())
    #ResImFilt.SetOutputOrigin(RefIm.GetOrigin())
    
    #ResImFilt.SetDefaultPixelValue(RefIm.GetPixelIDValue())
    #ResImFilt.SetOutputPixelType(RefIm.GetPixelID())
    ResImFilt.SetOutputPixelType(SitkPixType)
    ResImFilt.SetDefaultPixelValue(0)
    
    ResIm = ResImFilt.Execute(Im)
    
    return ResIm


def ResampleLabIm(LabIm, F2Sinds, SrcIm, TrgIm, Interp, 
                  PreResVariance=(1,1,1), ApplyPostResBlur=True, 
                  PostResVariance=(1,1,1), LogToConsole=False):
    """
    07/07/21:
        See resample_labim in image_tools.resampling.py
    
    Resample a 3D SimpleITK image representing a binary labelmap using
    either a NearestNeighbor or LabelGaussian interpolator, or by Gaussian
    blurring the labelmap before linearly resampling and binary thresholding.  
    If resampling using a NearestNeighbor or LabelGaussian interpolation,
    aliasing effects can result in a empty resampled labelmap.  If so the input 
    labelmap will be Gaussian blurred, linearly resampled and binary 
    thresholded.  See Notes for more details. 
    
    Inputs:
    ******                      
        
    LabIm : SimpleITK image
        The 3D labelmap image to be resampled.
    
    F2Sinds : list of integers
        List (for each frame) of slice numbers that correspond to each frame in 
        LabIm.
        
    SrcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of LabIm. 
    
    TrgIm : SimpleITK image
        The Target 3D image whose gridspace the LabIm will be resampled to.
    
    Interp : string
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
    
    PreResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior to resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to binarise the labelmap image following 
    #    resampling using a non-label (e.g. linear) interpolator. Example use is
    #    on an image if Gaussian blurring + linearly resampling approach is
    #    taken.  
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        labelmap image is to be Gaussian blurred prior after resampling.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    ResLabIm : SimpleITK images
        The resampled 3D labelmap image.
        
    ResPixArr : Numpy arrays
        The resampled pixel array (converted from ResLabIm).
    
    ResF2Sinds : list of integers
        The list (for each frame) of the frame-to-slice indices in ResPixArr.

    
    Notes:
    *****
    
    A NearestNeighbor or LabelGaussian interpolation are typically applied for
    binary (labelmap) images.  
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. an empty resampled labelmap,
    and hence an empty list of frame-to-slice indices, F2Sinds = []) even if 
    there were contours (i.e. if SrcC2Sinds != []). If this is the case try
    suggestions from Ziv Yaniv (See Link below):
            
    1. Use the sitkLabelGaussian interpolator.
    
    2. Gaussian blur the segmentation, then resample using a linear 
    interpolator, and then threshold so that values in [a<1<b] are mapped to 1.
    You will need to select appropriate values for a, b.
    
    3. Compute the distance map of the segmentation, then resample using a
    linear interpolator, and then threshold the absolute value of the distance
    map (sitk.Abs(lowres_distance_map)<dist_threshold). Again, you need to 
    select an appropriate dist_threshold.


    What has been tried so far:
    
    1. Resulting Gaussian interpolated resampled image had zeros only 
    (the image was converted to float prior to resampling).
    
    2. Works.
    
    
    Comment from Bradley Lowecamp:
        
    The implementation of the sitkLabelGaussian interpolator in the 
    ResampleImageFilter sets the sigma to the pixel spacing of the input image 
    to the resample filter. This is why it does not do a better job in reducing 
    aliasing.
        
    
    Link: 
        
    https://github.com/SimpleITK/SimpleITK/issues/1277
    """
    
    #import numpy as np
    from ConversionTools import Image2PixArr
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of ResampleLabIm():')
        print('\n\n', '-'*120)
        
    if not Interp in ['NearestNeighbor', 'LabelGaussian', 'BlurThenLinear']:
        msg = f'The chosen interpolation, {Interp}, is not one of the '\
              + 'accepted inputs: \'NearestNeighbor\', \'LabelGaussian\', or '\
              + '\'BlurThenLinear\'.'
        
        raise Exception(msg)
    
    
    """ Store the interpolation set as metadata: """
    LabIm.SetMetaData("ResInterpSet", Interp)
    
    
    #SrcSpacings = SrcIm.GetSpacing()
    #TrgSpacings = TrgIm.GetSpacing()
    
    #SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
    #                                                      TrgSpacings))
    
    #SrcSpacings = np.array(SrcIm.GetSpacing())
    #TrgSpacings = np.array(TrgIm.GetSpacing())
        
    #SpacingsRatio = np.divide(TrgSpacings, SrcSpacings)
    
    #VolumeRatio = np.prod(SpacingsRatio)
    
    
    #""" 
    #The FWHM of the Gaussian is:
    #    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
    #    
    #Due to the Nyquist-Shannon sampling theorem, the FWHM should be at
    #least 2x the voxel spacings of TrgIm. For symmetry it makes sense for
    #one Source voxel width to be blurred to three Target voxels widths:
    #    sigma = (3*TrgSpacings)/2.355
    #"""
    
    
    #if LogToConsole:
    #    print(f'\nSrcSpacings = {SrcSpacings}')
    #    print(f'TrgSpacings = {TrgSpacings}')
    #    print(f'\nSpacingsRatio = {SpacingsRatio}')
    
    if Interp in ['NearestNeighbor', 'LabelGaussian']:
        if LogToConsole:
            print('Attempting to resample LabIm using the',
                  f'{Interp} interpolator...\n')
        
        """ Attempt to resample LabIm to the Target image's grid using the
        chosen labelmap interpolator. """
        ResLabIm = ResampleImage(Im=LabIm, RefIm=TrgIm, Tx='Identity',
                                 Interp=Interp)
        
        if LogToConsole:
            print('Image info for ResLabIm after resampling using',
                  f'{Interp} interpolator:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        ResF2Sinds = GetImageInfo(ResLabIm, LogToConsole)
        
        if LogToConsole:
            print('')
        
        
        if ApplyPostResBlur:
            """ Gaussian blur the image. """
            ResLabIm = GaussianBlurImage(Im=ResLabIm, 
                                         Variance=PostResVariance)
            
            if LogToConsole:
                print('Image info for ResLabIm after Gaussian blurring:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            ResF2Sinds = GetImageInfo(ResLabIm, LogToConsole)
            
            if LogToConsole:
                print('')
            
            
            #""" The original binary pixel array: """
            #PixArr_B, F2Sinds = Image2PixArr(LabIm)
            #
            #""" The non-binary pixel array following blur: """
            #PixArr_NB, F2Sinds = Image2PixArr(ResLabIm)
            #
            #Ratio = GetVoxelVolumeRatio(LabIm, ResLabIm)
            #
            #if LogToConsole:
            #    print(f'PixArr_B.shape = {PixArr_B.shape}')
            #    print(f'PixArr_NB.shape = {PixArr_NB.shape}\n')
            
            
            """ Binarise the resampled labelmap if required. """
            if len(UniqueVals) != 2 or sum(UniqueVals) != 1:
                """Find suitable threshold value that approximately preserves
                the number of pre-blurred truth values scaled by VolumeRatio: 
                    """
                Thresh = FindSuitableThresh(BinaryIm=LabIm, 
                                            NonBinaryIm=ResLabIm,
                                            LogToConsole=LogToConsole)
                
                """ Binary threshold the image. """
                #ThreshValue = ThreshLevel*max(UniqueVals)
                #ThreshValue = ThreshLevel
                
                ResLabIm = BinaryThresholdImage(Im=ResLabIm, 
                                                #Thresh=PostResThresh)
                                                Thresh=Thresh) 
                
                if LogToConsole:
                    print('\nImage info for ResLabIm after binary',
                          #f'thresholding at {PostResThresh}:')
                          f'thresholding at {Thresh}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabIm, LogToConsole)
                
                if LogToConsole:
                    print('')
        
        
        #print(f'\n   ResF2Sinds = {ResF2Sinds}')
        
        """ Is ResF2Sinds empty and not expected to be? If so, try the 
        "BlurThenLinear" approach.
        
        Note:
            If F2Sinds isn't empty, ResF2Sinds shouldn't be empty either.
        
            F2Sinds will be empty if there were no segmentations/contours 
            of interest for the r^th ROI. In this case an empty ResF2Sinds 
            is acceptable. """
        
        if ResF2Sinds == []:
            print(f'There are {len(F2Sinds)} non-empty masks in the input',
                  f'labelmap image but {len(ResF2Sinds)} non-empty frames in',
                  f'the resampled labelmap image using {Interp}.',
                  'Try using the \'BlurThenLinear\' approach.\n')
            
            Interp = 'BlurThenLinear'

    
    
    if Interp == 'BlurThenLinear':
        if LogToConsole:
            print('\nResampling LabIm by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        
        #""" The original binary pixel array: """
        #PixArr_B, F2Sinds = Image2PixArr(LabIm)
            
        
        """ Convert LabIm from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        LabIm = ChangeImagePixelType(Image=LabIm, 
                                     NewPixelType='Float32')
        
        if LogToConsole:
            print('\nImage info for LabIm after converting to float:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(LabIm, LogToConsole)
              
        
              
        """ Gaussian blur the image. """
        LabIm = GaussianBlurImage(Im=LabIm, Variance=PreResVariance)
        
        # Use the RecursiveGaussian image filter:
        #ResLabIm = RecursiveGaussianBlurImage(ResLabIm, 
        #                                      Sigma=3,
        #                                      Direction=2)
        
        if LogToConsole:
            print('\nImage info for LabIm after applying Gaussian blur:')
            
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(LabIm, LogToConsole)
        
        if LogToConsole:
            print('\nLabIm prior to resampling:')
            print(f'   LabIm.GetSize() = {LabIm.GetSize()}')
            print(f'   LabIm.GetSpacing() = {LabIm.GetSpacing()}')
            print(f'   TrgIm.GetSize() = {TrgIm.GetSize()}')
            print(f'   TrgIm.GetSpacing() = {TrgIm.GetSpacing()}')
        
        
        """ Linearly resample LabIm to the Target image's grid. """
        ResLabIm = ResampleImage(Im=LabIm, RefIm=TrgIm, Tx='Identity',
                                 Interp='Linear')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after resampling using a',
                  'linear interpolator:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabIm, LogToConsole)
        
        
        if ApplyPostResBlur:
            """ Gaussian blur the image. """
            ResLabIm = GaussianBlurImage(Im=ResLabIm, 
                                         Variance=PostResVariance)
        
        
        #""" The non-binary pixel array following blur + linear resampling: """
        #PixArr_NB, F2Sinds = Image2PixArr(ResLabIm)
        #
        #if LogToConsole:
        #        print(f'PixArr_B.shape = {PixArr_B.shape}')
        #        print(f'PixArr_NB.shape = {PixArr_NB.shape}\n')
            
        """Find suitable threshold value that approximately preserves the
        number of pre-blurred + linear resampled truth values scaled by
        VolumeRatio: """
        Thresh = FindSuitableThresh(BinaryIm=LabIm, 
                                    NonBinaryIm=ResLabIm,
                                    LogToConsole=LogToConsole)
        
        """ Binary threshold the image. """
        #ThreshValue = ThreshLevel*max(UniqueVals)
        #ThreshValue = ThreshLevel
        
        ResLabIm = BinaryThresholdImage(Im=ResLabIm, 
                                        #Thresh=PostResThresh)
                                        Thresh=Thresh) 
        
        if LogToConsole:
            print('\nImage info for ResLabIm after binary thresholding',
                  #f'at {PostResThresh}:')
                  f'at {Thresh}:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabIm, LogToConsole)
    
    
    
    """ Ensure that ResLabIm is a 32-bit unsigned integer (PixID = 5). """
    if PixID != 5: 
        if LogToConsole:
            print(f'\nResLabIm has PixelID = {PixID} ({PixIDTypeAsStr})).')
        
        """ Convert ResLabIm from float to 32-bit unsigned integer. """
        ResLabIm = ChangeImagePixelType(Image=ResLabIm, 
                                        NewPixelType='UInt32')
        
        if LogToConsole:
            print('\nImage info for ResLabIm after converting to 32-bit',
                  'unsigned int:')
            
        #print(f'\nThe metadata keys are:', ResLabIm.GetMetaDataKeys())
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabIm, LogToConsole)
    
    
    """ Convert ResLabIm to a pixel array. """
    ResPixArr, ResF2Sinds = Image2PixArr(ResLabIm)
    
    
    """ Store the interpolation used as metadata (which may be the same or 
    different from the interpolation set): """
    ResLabIm.SetMetaData("ResInterpUsed", Interp)
    
    if LogToConsole:
        print(f'\nThe key "ResInterpUsed" with value "{Interp}" was',
              'added to the metadata of ResLabIm. \nThe metadata keys are:',
              ResLabIm.GetMetaDataKeys())
    
    if Interp == 'BlurThenLinear':
        """ Store the threshold used as metadata: """
        ResLabIm.SetMetaData("PostResThreshUsed", f"{Thresh}")
        
        if LogToConsole:
            print(f'\nThe key "PostResThreshUsed" with value "{Thresh}" was',
                  'added to the metadata of ResLabIm. \nThe metadata keys:',
                  ResLabIm.GetMetaDataKeys())
            
    if LogToConsole:
        print('\nAfter converting ResLabIm to a pixel array:')
        print(f'ResPixArr.shape = {ResPixArr.shape}')
        print(f'ResF2Sinds = {ResF2Sinds}')
        print('-'*120)
        
    return ResLabIm, ResPixArr, ResF2Sinds


def ResampleLabImByRoi(LabImByRoi, F2SindsByRoi, SrcIm, TrgIm, 
                       Interp, PreResVariance, ApplyPostResBlur, 
                       PostResVariance, LogToConsole=False):
    """
    07/07/21:
        See resample_labimByRoi in image_tools.resampling.py
    
    Resample a list 3D SimpleITK images representing binary labelmaps using
    either a NearestNeighbor or LabelGaussian interpolator, or by Gaussian
    blurring the labelmap before linearly resampling and binary thresholding.  
    If resampling using a NearestNeighbor or LabelGaussian interpolation,
    aliasing effects can result in a empty resampled labelmap.  If so the input 
    labelmap will be Gaussian blurred, linearly resampled and binary 
    thresholded.  See Notes in ResampleLabImByRoi() for more details. 
    
    Inputs:
    ******                      
        
    LabImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled.
    
    F2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each frame) of slice numbers that 
        correspond to each frame in the labelmap images.
        
    SrcIm : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of the 
        labelmap images in LabelmapImByRoi.
    
    TrgIm : SimpleITK image
        The Target 3D image whose gridspace the Source labelmap images will be
        resampled to.
    
    Interp : string
        The type of interpolation to be used.  Acceptable inputs are:
        - 'NearestNeighbor'
        - 'LabelGaussian'
        - 'BlurThenLinear' (which represents a Gaussian image blur + resampling
        using a linear interpolator + binary thresholding)
        
    PreResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        Source labelmap image(s) is/are to be Gaussian blurred prior to  
        resampling.
        
    #PostResThresh : float (optional; 0.75 by default)
    #    The threshold level used to binarise the labelmap image following 
    #    resampling using a non-label (e.g. linear) interpolator. Example use is
    #    on an image if Gaussian blurring + linearly resampling approach is
    #    taken.  
    
    ApplyPostResBlur : boolean (optional; True by default)
        If True, the post-resampled labelmap image will be Gaussian blurred.
    
    PostResVariance : tuple of floats (optional; (1,1,1) by default)
        A tuple (for each dimension) of the variance to be applied if the 
        resampled labelmap image(s) is/are to be Gaussian blurred after  
        resampling.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    ResLabImByRoi : list of SimpleITK images
        The list (for each ROI) of resampled 3D labelmap images.
        
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of resampled pixel arrays (converted from each
        labelmap in ResLabImByRoi).
    
    ResF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.
    """
        
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of ResampleLabImByRoi():')
        print('\n\n', '-'*120)
        
    ResLabImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabImByRoi)):
        if LogToConsole:
            print(f'   Resampling of LabImByRoi[{r}]...')
        
        ResLabIm, ResPixArr,\
        ResF2Sinds = ResampleLabIm(LabImByRoi[r], F2SindsByRoi[r], 
                                   SrcIm, TrgIm, Interp, 
                                   PreResVariance, ApplyPostResBlur, 
                                   PostResVariance, LogToConsole)
        
        ResLabImByRoi.append(ResLabIm)
        ResPixArrByRoi.append(ResPixArr)
        ResF2SindsByRoi.append(ResF2Sinds)
        
    if LogToConsole:
            print('-'*120)
        
    return ResLabImByRoi, ResPixArrByRoi, ResF2SindsByRoi


def GetLandmarkPtsFromImageJInds_OLD(Inds, Im):
    """ 
    Get a flat list of points from a list of indices from ImageJ's
    coordinate system.
    
    Inputs:
    ******
    
    Inds : list of a list of integers
        List (for each fiducial) of a list (for each dimension) of indices of
        landmarks (e.g. [[i0, j0, k0], [i1, j1, k1], ...].
    
    Im : SimpleITK image
    
    
    Outputs:
    *******
    
    Pts : list of floats
        Flat list of [x, y, z] coordinates for each index in Indices (e.g.
        [x0, y0, z0, x1, y1, z1, ...]).
    """
    
    #import SimpleITK as sitk
    
    R, C, S = Im.GetSize()
        
    """ Slice indexing in ImageJ is opposite to that of sitk. """
    Inds = [[ind[0], ind[1], S - ind[2]] for ind in Inds]
    
    Pts = [Im.TransformIndexToPhysicalPoint(ind) for ind in Inds]
    
    """ Flatten Pts """
    Pts = [c for p in Pts for c in p]
    
    return Pts


def ImportFiducialsFromTxt(Filepath, FlipK=False, RefIm=None,
                           LogToConsole=False):
    """
    08/07/21:
        See import_fiducials in general_tools.fiducials.py
    
    Import a list of indices or points from a suitably formatted text file.
    
    Inputs:
    ******
        
    Filepath : string
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    
    FlipK : boolean (optional; False by default)
        If the fiducials are indices derived from ImageJ, set to True to flip
        the z-index, i.e. [i, j, S - k], where S = the number of slices in the
        3D image.
    
    RefIm : SimpleITK image (optional unless FlipK is True)
        The 3D SimpleITK image that corresponds to the fiducials.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
        
    Fiducials : list of a list of integers or floats
        A list (for each index/point) of a list (for each dimension) of either
        indices (integers) or points (floats), 
        
        e.g. [[i0, j0, k0], [i1, j1, k1], ...] 
        
        or [[x0, y0, z0], [x1, y1, z1], ...]
    
    FiducialType : string
        Denotes the fiducial type as 'index' or 'point'.
    
    
    Notes:
    *****
    
    The text file must have the format expected of Elastix/SimpleElastix, i.e.:
        
    <index or point>
    <number of indices/points>
    point1 x point1 y [point1 z]
    point2 x point2 y [point2 z]
    
    e.g.
    
    point
    3
    102.8 -33.4 57.0
    178.1 -10.9 14.5
    180.4 -18.1 78.9
    
    https://simpleelastix.readthedocs.io/PointBasedRegistration.html
    """
    
    if FlipK == True and RefIm == None:
        msg = "'FlipK' is True so a reference image must be provided."
        
        raise Exception(msg)
    
    if not '.txt' in Filepath:
        Filepath += '.txt'
    
    Fiducials = []
    
    with open(Filepath, 'r') as file:
        contents = file.readlines()
            
        FiducialType = contents[0].replace('\n', '')
        
        if FiducialType not in ['index', 'point']:
            msg = f"The first line in {Filepath} must be either 'index' or "\
                  + f"'point'.  '{FiducialType}' is not a valid entry."
            raise Exception(msg)
            
        n = int(contents[1].replace('\n', ''))
            
        for i in range(2, len(contents)):
            items = contents[i].replace('\n', '').split(' ')
            
            if FiducialType == 'index':
                fiducial = [int(item) for item in items]
                
                if FlipK:
                    S = RefIm.GetSize()[2]
                    
                    #fiducial[2] = S - fiducial[2] # 26/04/21
                    fiducial[2] = (S - 1) - fiducial[2] # 26/04/21
                    
            else:
                fiducial = [float(item) for item in items]
            
            Fiducials.append(fiducial)
    
    N = len(Fiducials)
    
    if N != n:
        msg = f"Warning: The second line in {Filepath} specifies that there "\
              f"are {n} fiducials but {N} fiducials were parsed."
        
        print(msg)
    
    if LogToConsole:
        print(f'Fiducial type: {FiducialType}')
        print(f'There are {len(Fiducials)} fiducials: \n {Fiducials}\n') 
    
    return Fiducials, FiducialType


def GetFiducialsWithinImExtents(FixFids, MovFids, FixIm, MovIm, Buffer=3,
                                LogToConsole=False):
    """
    08/07/21:
        See get_fiducials_within_im_extents in general_tools.fiducials.py
    
    Return fiducials belonging to two images which lie within its image extent, 
    and do not lie on the boundary or to within Buffer pixels of it.  The 
    fiducials must be indices (not points).
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
    
    FixFids : list of a list of integers
        A list (for each fiducial) of a list of [i, j, k] indices for FixIm.
    
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    MovFids : list of a list of integers
        A list (for each fiducial) of a list of [i, j, k] indices for MovIm.
    
    Buffer : integer (optional; 3 by default)
        Rather than eliminating all fiducials at or beyond the boundary, also
        filter out any fiducials that lie within Buffer pixels of any boundary.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    
    Outputs:
    *******
    
    FixFids : list of a list of integers
        A list (for each fiducial) of a list of [i, j, k] indices for FixIm.
    
    MovFids : list of a list of integers
        A list (for each fiducial) of a list of [i, j, k] indices for MovIm.
    
    CommonInds : list of integers
        A list of indices of the original fiducials that lie within both image
        boundaries.
    
    
    Note:
    *****
    
    For BSpline Transforms it's not enough for the points to lie within the    
    extent of the image - they also must not lie on the boundary:

    https://itk.org/pipermail/insight-users/2012-June/044975.html
    
    It seems that a buffer of 2 (minimum) or 3 pixels is appropriate.
    """
    
    I_f, J_f, K_f = FixIm.GetSize()
    I_m, J_m, K_m = MovIm.GetSize()
    
    # Store the indices of the fiducials that do not lie on any boundary of 
    # either image:
    CommonInds = []
    
    for n in range(len(FixFids)):
        i_f, j_f, k_f = FixFids[n]
        i_m, j_m, k_m = MovFids[n]
        
        # The n^th fiducial lies within the boundaries of both images if all of
        # the following are True:
        b0 = i_f > 0 + Buffer
        b1 = i_f < I_f - 1 - Buffer
        b2 = j_f > 0 + Buffer
        b3 = j_f < J_f - 1 - Buffer
        b2 = k_f > 0 + Buffer
        b3 = k_f < K_f - 1 - Buffer
        b4 = i_m > 0 + Buffer
        b5 = i_m < I_m - 1 - Buffer
        b6 = j_m > 0 + Buffer
        b7 = j_m < J_m - 1 - Buffer
        b8 = k_m > 0 + Buffer
        b9 = k_m < K_m - 1 - Buffer
        
        if b0 & b1 & b2 & b3 & b4 & b5 & b6 & b7 & b8 & b9:
            CommonInds.append(n)
        
    FixFids = [FixFids[i] for i in CommonInds]
    MovFids = [MovFids[i] for i in CommonInds]
    
    if LogToConsole:
        print(f'After rejecting fiducials on or to within {Buffer} pixels of',
              f'the image boundaries there are now {len(FixFids)} fiducials.\n',
              f'FixFids:\n {FixFids}\nMovFids: {MovFids}\n')
    
    return FixFids, MovFids, CommonInds
    

def GetLandmarkPtsFromImageJInds(FiducialFpath, Image):
    """ 
    08/07/21:
        See get_landmark_pts_from_ImageJ_inds in general_tools.fiducials.py
    
    Get a flat list of points from a list of indices from ImageJ's
    coordinate system.
    
    Inputs:
    ******
    
    FiducialFpath : string
        The file path (or file name if the file is present in the current 
        working directory) of the text file. The string need not contain the 
        .txt extension.
    
    Image : SimpleITK image
        The 3D DICOM image for which the fiducials relate to.
    
    
    Outputs:
    *******
    
    Pts : list of floats
        Flat list of [x, y, z] coordinates converted from each index the text
        file, e.g. [x0, y0, z0, x1, y1, z1, ...]).
    """
    
    #import SimpleITK as sitk
    
    R, C, S = Image.GetSize()
    
    Inds, FiducialType = ImportFiducialsFromTxt(FiducialFpath)
    
    if FiducialType != 'index':
        msg = "The text file must contain indices, not points."
        raise Exception(msg)
        
    """ Slice indexing in ImageJ is opposite to that of sitk. """
    Inds = [[ind[0], ind[1], S - ind[2]] for ind in Inds]
    
    Pts = [Image.TransformIndexToPhysicalPoint(ind) for ind in Inds]
    
    """ Flatten Pts """
    Pts = [c for p in Pts for c in p]
    
    return Pts


def GetLandmarkTx_OLD(FixIm, MovIm, FixFidsFpath, MovFidsFpath):
    import SimpleITK as sitk

    FixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath)

    MovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath)

    if FixFidType == 'index':
        FixFids = [FixIm.TransformIndexToPhysicalPoint(ind) for ind in FixFids]

    if MovFidType == 'index':
        MovFids = [MovIm.TransformIndexToPhysicalPoint(ind) for ind in MovFids]
    
    print(f'Fixed fiducials: \n {FixFids}\n')
    print(f'Moving fiducials: \n {MovFids}\n')
    
    """ Flatten the list of points """
    FixPts = [c for p in FixFids for c in p]
    MovPts = [c for p in MovFids for c in p]

    LandmarkTx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
                                                        FixPts, MovPts)
    
    print(f'Landmark-based initial transformation is: \n {LandmarkTx.GetParameters()}')
    
    
    if False:
        # https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/63_Registration_Initialization.ipynb
        LandmarkTx = sitk.VersorRigid3DTransform(sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
                                                                                        FixPts, MovPts))

        # Convert from Versor to Euler, as Versor does not always work well with the optimization. 
        # Internally the optimization sets new parameter values without any constraints, and the versor
        # normalizes its vector component if it is greater than 1-epsilon.
        InitialTx = sitk.Euler3DTransform()
        InitialTx.SetCenter(LandmarkTx.GetCenter())
        InitialTx.SetMatrix(LandmarkTx.GetMatrix())
        InitialTx.SetTranslation(LandmarkTx.GetTranslation())

        print(f'Landmark-based initial transformation is: \n {InitialTx.GetParameters()}')

        return InitialTx
    
    return LandmarkTx


def GetLandmarkTx_v2_OLD(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
                     AlignmentTransform='rigid', NumControlPts=8):
    
    import SimpleITK as sitk
    
    if AlignmentTransform == 'rigid':
        AlignmentTx = sitk.VersorRigid3DTransform()
    elif AlignmentTransform == 'affine':
        AlignmentTx = sitk.AffineTransform(FixIm.GetDimension())
    elif AlignmentTransform == 'bspline':
        AlignmentTx = sitk.BSplineTransform(FixIm.GetDimension())
        #AlignmentTx = sitk.BSplineTransform(FixIm.GetDimension(), 3)
    else:
        msg = "AlignmentTransform must be either 'rigid', 'affine' or "\
              + "'bspline'."
        raise Exception(msg)

    FixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath)

    MovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath)

    if FixFidType == 'index':
        FixFids = [FixIm.TransformIndexToPhysicalPoint(ind) for ind in FixFids]

    if MovFidType == 'index':
        MovFids = [MovIm.TransformIndexToPhysicalPoint(ind) for ind in MovFids]
    
    print(f'Fixed fiducials: \n {FixFids}\n')
    print(f'Moving fiducials: \n {MovFids}\n')
    
    """ Flatten the list of points """
    FixPts = [c for p in FixFids for c in p]
    MovPts = [c for p in MovFids for c in p]
    
    if AlignmentTransform == 'bspline':
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=AlignmentTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm,
                                                            numberOfControlPoints=NumControlPts)
                                                            #BSplineNumberOfControlPoints=NumControlPts)
    else:
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=AlignmentTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm)
    
    
    
    print(f'Landmark-based initial transformation is: \n {LandmarkTx.GetParameters()}')
    
    
    if False:
        # https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/63_Registration_Initialization.ipynb
        LandmarkTx = sitk.VersorRigid3DTransform(sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
                                                                                        FixPts, MovPts))

        # Convert from Versor to Euler, as Versor does not always work well with the optimization. 
        # Internally the optimization sets new parameter values without any constraints, and the versor
        # normalizes its vector component if it is greater than 1-epsilon.
        InitialTx = sitk.Euler3DTransform()
        InitialTx.SetCenter(LandmarkTx.GetCenter())
        InitialTx.SetMatrix(LandmarkTx.GetMatrix())
        InitialTx.SetTranslation(LandmarkTx.GetTranslation())

        print(f'Landmark-based initial transformation is: \n {InitialTx.GetParameters()}')

        return InitialTx
    
    return LandmarkTx


def GetLandmarkTx(Transform, FixFidsFpath, MovFidsFpath, FixIm, MovIm,
                  FlipK=False, Buffer=3, LogToConsole=False):
    """
    08/07/21:
        See get_landmark_tx in general_tools.fiducials.py
    
    Return the landmark-based transform based on fiducials.
    
    Inputs:
    ******
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension.
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be aligned to.
        
    MovIm : SimpleITK image
        The 3D image that will be aligned to FixIm.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ.
    
    Buffer : integer (optional; 3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    LandmarkTx : SimpleITK transform
        The transform used to perform landmark-based alignment.
    
    FixPts : list of floats
        Flat list of landmark coordinates for FixIm.
    
    MovPts : list of floats
        Flat list of landmark coordinates for MovIm.
    
    
    Notes:
    *****
    
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(FixIm.GetDimension())
        sitk.BSplineTransform(FixIm.GetDimension())
        
    NumControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    import SimpleITK as sitk
    from ConversionTools import Indices2Points
    from GeneralTools import FlattenList
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {Transform}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    dim = FixIm.GetDimension()
    
    if Transform == 'rigid':
        SitkTx = sitk.VersorRigid3DTransform()
    elif Transform == 'affine':
        SitkTx = sitk.AffineTransform(dim)
    else:
        SitkTx = sitk.BSplineTransform(dim)
    
    """ Import the fiducials from text files: """
    FixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath, 
                                                 LogToConsole=LogToConsole)

    MovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath,
                                                 LogToConsole=LogToConsole)
    
    
    if Transform == 'bspline':
        """ Reject any fiducials at or near any image boundary: """
        
        print(f'Transform = {Transform} so will reject fiducials at or near',
              'the boundary.\n')
        
        FixFids, MovFids, CommonInds\
            = GetFiducialsWithinImExtents(FixFids=FixFids, MovFids=MovFids, 
                                          FixIm=FixIm, MovIm=MovIm,
                                          Buffer=Buffer,
                                          LogToConsole=LogToConsole)
    
    """ Convert from lists of indices to lists of points: """
    if FixFidType == 'index':
        FixPts = Indices2Points(Indices=FixFids, RefIm=FixIm)
    
    if MovFidType == 'index':
        MovPts = Indices2Points(Indices=MovFids, RefIm=MovIm)
    
    if LogToConsole:
        print(f'There are {len(FixPts)} fixed points: \n {FixPts}\n')
        print(f'There are {len(MovPts)} moving points: \n {MovPts}\n')
    
    """ Flatten the list of points (as required by 
    LandmarkBasedTransformInitializer): """
    FixPts = FlattenList(FixPts)
    MovPts = FlattenList(MovPts)
    
    #""" Additional check (for troubleshooting) - are all points within image
    #extents? """
    #from GeneralTools import IsPointInPolygon
    #
    #FixVerts = GetImageVertices(Image=FixIm)
    #MovVerts = GetImageVertices(Image=MovIm)
    #
    #FixPtsInExtent = []
    #MovPtsInExtent = []
    #
    #for i in range(len(FixPts)):
    #    FixPtsInExtent.append(IsPointIPolygon(FixPts[], FixVerts))
        
    
    if Transform == 'bspline':
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=SitkTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm,
                                                            numberOfControlPoints=4)
                                                            #BSplineNumberOfControlPoints=NumControlPts)
    else:
        LandmarkTx = sitk.LandmarkBasedTransformInitializer(transform=SitkTx, 
                                                            fixedLandmarks=FixPts, 
                                                            movingLandmarks=MovPts,
                                                            referenceImage=FixIm)
    
    # Cast LandmarkTx from sitkTransform to the desired sitk transform:
    # (01/07/21)
    if Transform == 'rigid':
        LandmarkTx = sitk.VersorRigid3DTransform(LandmarkTx)
    elif Transform == 'affine':
        LandmarkTx = sitk.AffineTransform(LandmarkTx)
    else:
        LandmarkTx = sitk.BSplineTransform(LandmarkTx)
        
    if LogToConsole:
        #print(f'Fixed fiducials (points): \n {FixPts}\n')
        #print(f'Moving fiducials (points): \n {MovPts}\n')
        #print(f'Landmark-based initial transformation is: \n {LandmarkTx.GetParameters()}')
        print("-------")
        print('Landmark-based initial transform:')
        print(LandmarkTx)
        print("-------")
    
    #return LandmarkTx # 01/07/21
    return LandmarkTx, FixPts, MovPts


def AlignImagesUsingLandmarks_OLD(FixIm, MovIm, 
                              FixFidsFpath='fixed_fiducials.txt', 
                              MovFidsFpath='moving_fiducials.txt',
                              FlipK=False, Transform='affine',
                              LogToConsole=False):
    """
    Return the landmark-based aligned image and transform based on fiducials.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be aligned to.
        
    MovIm : SimpleITK image
        The 3D image that will be aligned to FixIm.
        
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ.
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    AlignedIm : SimpleITK image
        The 3D aligned image.
    
    LandmarkTx : SimpleITK transform
        The transform used to perform landmark-based alignment.
        
    
    Notes:
    *****
    
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(FixIm.GetDimension())
        sitk.BSplineTransform(FixIm.GetDimension())
    """
    
    import os
    import SimpleITK as sitk
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {Transform}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    if not '.txt' in FixFidsFpath:
        FixFidsFpath += '.txt'
    
    if not '.txt' in MovFidsFpath:
        MovFidsFpath += '.txt'
        
    if not os.path.isfile(FixFidsFpath):
        msg = f"The file {FixFidsFpath} does not exist."
        print(msg)
    
    if not os.path.isfile(MovFidsFpath):
        msg = f"The file {MovFidsFpath} does not exist."
        print(msg)
    
    
    #FixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath) # 26/04/21
    #MovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath) # 26/04/21
    
    FixFids, FixFidType = ImportFiducialsFromTxt(FixFidsFpath,
                                                 FlipK=FlipK,
                                                 #Im=MovIm) # 26/04/21
                                                 Im=FixIm) # 28/04/21
    
    MovFids, MovFidType = ImportFiducialsFromTxt(MovFidsFpath,
                                                 FlipK=FlipK,
                                                 Im=MovIm) # 26/04/21
    
    if FlipK:
        FixR, FixC, FixS = FixIm.GetSize()
        MovR, MovC, MovS = MovIm.GetSize()
        
        """ Slice indexing in ImageJ is opposite to that of sitk. """
        FixFids = [[ind[0], ind[1], FixS - ind[2]] for ind in FixFids]
        MovFids = [[ind[0], ind[1], MovS - ind[2]] for ind in MovFids]
    
    """ LandmarkBasedTransformInitializer requires points. Convert from indices
    if required. """
    if FixFidType == 'index':
        FixFids = [FixIm.TransformIndexToPhysicalPoint(ind) for ind in FixFids]
    
    if MovFidType == 'index':
        MovFids = [MovIm.TransformIndexToPhysicalPoint(ind) for ind in MovFids]
    
        
    """ Flatten the list of points """
    FixPts = [c for p in FixFids for c in p]
    MovPts = [c for p in MovFids for c in p]
    
    dim = FixIm.GetDimension()
    
    ##Tx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
    ##                                            FixPts, MovPts)
    #Tx = sitk.LandmarkBasedTransformInitializer(sitk.AffineTransform(dim), 
    #                                            FixPts, MovPts)
    
    if Transform == 'rigid':
        SitkTx = sitk.VersorRigid3DTransform()
    elif Transform == 'affine':
        SitkTx = sitk.AffineTransform(dim)
    else:
        SitkTx = sitk.BSplineTransform(dim)
    
    LandmarkTx = sitk.LandmarkBasedTransformInitializer(SitkTx, FixPts, MovPts)
    
    if LogToConsole:
        #params = LandmarkTx.GetParameters()
        print(f'Fixed fiducials: \n {FixFids}\n')
        print(f'Moving fiducials: \n {MovFids}\n')
        #print(f'Landmark-based initial transformation is:\n {params}\n')
        print("-------")
        print('LandmarkTx:')
        print(LandmarkTx)
        print("-------")
    
    AlignedIm = sitk.Resample(MovIm, FixIm, LandmarkTx, sitk.sitkLinear, 0.0, 
                              MovIm.GetPixelID())
    
    return AlignedIm, LandmarkTx


def AlignImagesUsingLandmarks(Transform, FixFidsFpath, MovFidsFpath, 
                              FixIm, MovIm, FlipK=False, Buffer=3, 
                              LogToConsole=False):
    """
    08/07/21:
        See align_ims_using_landmarks in image_tools.resampling.py
    
    Return the landmark-based aligned image and transform based on fiducials.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be aligned to.
        
    MovIm : SimpleITK image
        The 3D image that will be aligned to FixIm.
        
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ.
    
    Buffer : integer (optional; 3 by default)
        Fiducials at or within this number of pixels within the image 
        boundaries will be rejected.
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to set up the landmark transform.  
        Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
    
    Outputs:
    *******
    
    AlignedIm : SimpleITK image
        The 3D aligned image.
    
    LandmarkTx : SimpleITK transform
        The transform used to perform landmark-based alignment.
        
    
    Notes:
    *****
    
    Supported 3D transforms for LandmarkBasedTransformInitializer are:
        sitk.VersorRigid3DTransform()
        sitk.AffineTransform(FixIm.GetDimension())
        sitk.BSplineTransform(FixIm.GetDimension())
        
    NumControlPts must be greater than the spline order (=3), i.e. 4.
    But any value greater than 4 results in a highly distorted aligned image.
    """
    
    import SimpleITK as sitk
    
    LandmarkTx, FixPts, MovPts\
        = GetLandmarkTx(Transform=Transform, 
                        FixFidsFpath=FixFidsFpath, 
                        MovFidsFpath=MovFidsFpath, 
                        FixIm=FixIm, MovIm=MovIm, FlipK=FlipK,
                        Buffer=Buffer, LogToConsole=LogToConsole)
    
    
    #print(f'LandmarkTx.GetName() = {LandmarkTx.GetName()}\n')
    
    AlignedIm = sitk.Resample(MovIm, FixIm, LandmarkTx, sitk.sitkLinear, 0.0, 
                              MovIm.GetPixelID())
    
    return AlignedIm, LandmarkTx


#import SimpleITK as sitk
#import matplotlib.pyplot as plt
##from ipywidgets import interact, fixed
#from IPython.display import clear_output


def PlotFixMovImages(FixInd, MovInd, FixPixArr, MovPixArr):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    import matplotlib.pyplot as plt
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr[FixInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Fixed image')
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr[MovInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Moving image')
    plt.axis('off')
    
    plt.show()


def PlotBlendedImage(Ind, alpha, FixIm, ResIm):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. """
    
    import SimpleITK as sitk
    import matplotlib.pyplot as plt

    BlendedIm = (1.0 - alpha)*FixIm[:,:,Ind] + alpha*ResIm[:,:,Ind] 
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.axis('off')
    plt.show()


def StartPlot():
    """ Callback invoked when the StartEvent happens, sets up our new data. """
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []


def EndPlot():
    """ Callback invoked when the EndEvent happens, do cleanup of data and 
    figure. """
    
    import matplotlib.pyplot as plt
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()


def PlotValues(RegMethod):
    """ Callback invoked when the IterationEvent happens, update our data 
    and display new figure. """
    
    import matplotlib.pyplot as plt
    from IPython.display import clear_output
    global metric_values, multires_iterations
    
    metric_values.append(RegMethod.GetMetricValue())                                       
    # Clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    # Plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
 
def UpdateMultiresIterations():
    """ Callback invoked when the sitkMultiResolutionIterationEvent happens, 
    update the index into the metric_values list. """
    
    global metric_values, multires_iterations
    
    multires_iterations.append(len(metric_values))


def RegisterImagesSitk_OLD(FixIm, MovIm, Transform='affine', 
                       InitMethod='moments', NumControlPts=8,
                       FinalInterp='Linear',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       FlipK=False, LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    FinalInterp : string (optional; 'Linear' by default)
        The type of final interpolation to be applied.  Acceptable inputs
        include:
            - 'Linear'
            - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ. This argument is only relevant if 
        InitMethod = 'landmarks'.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
        
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    """
    
    import SimpleITK as sitk
    import time
    import winsound
    
    #import matplotlib.pyplot as plt
    from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    """
    def PlotFixMovImages(FixInd, MovInd, FixPixArr, MovPixArr):
        # Callback invoked by the interact IPython method for scrolling 
        through the image stacks of the two images (moving and fixed).
        
        # Create a figure with two subplots and the specified size.
        plt.subplots(1, 2, figsize=(10,8))
        
        # Draw the fixed image in the first subplot.
        plt.subplot(1, 2, 1)
        plt.imshow(FixPixArr[FixInd,:,:], cmap=plt.cm.Greys_r);
        plt.title('Fixed image')
        plt.axis('off')
        
        # Draw the moving image in the second subplot.
        plt.subplot(1, 2, 2)
        plt.imshow(MovPixArr[MovInd,:,:], cmap=plt.cm.Greys_r);
        plt.title('Moving image')
        plt.axis('off')
        
        plt.show()
    
    
    def PlotBlendedImage(Ind, alpha, FixIm, ResIm):
        # Callback invoked by the IPython interact method for scrolling and 
        modifying the alpha blending of an image stack of two images that 
        occupy the same physical space.
        
        BlendedIm = (1.0 - alpha)*FixIm[:,:,Ind] + alpha*ResIm[:,:,Ind] 
        plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
        plt.axis('off')
        plt.show()
    
    
    def StartPlot():
        # Callback invoked when the StartEvent happens, sets up our new data.
        global metric_values, multires_iterations
        
        metric_values = []
        multires_iterations = []
    
    
    def EndPlot():
        # Callback invoked when the EndEvent happens, do cleanup of data and 
        figure.
        
        global metric_values, multires_iterations
        
        del metric_values
        del multires_iterations
        # Close figure, we don't want to get a duplicate of the plot latter on.
        plt.close()
    
    
    def PlotValues(RegMethod):
        # Callback invoked when the IterationEvent happens, update our data 
        # and display new figure.
        
        global metric_values, multires_iterations
        
        metric_values.append(RegMethod.GetMetricValue())                                       
        # Clear the output area (wait=True, to reduce flickering), and plot current data
        clear_output(wait=True)
        # Plot the similarity metric values
        plt.plot(metric_values, 'r')
        plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
        plt.xlabel('Iteration Number',fontsize=12)
        plt.ylabel('Metric Value',fontsize=12)
        plt.show()
        
     
    def UpdateMultiresIterations():
        # Callback invoked when the sitkMultiResolutionIterationEvent happens, 
        # update the index into the metric_values list.
        
        global metric_values, multires_iterations
        multires_iterations.append(len(metric_values))
    """
    
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {Transform}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    if not InitMethod in ['geometry', 'geometricalcenter', 'moments', 
                          'centerofgravity', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
        
    if not FinalInterp in ['Linear', 'NearestNeighbor']:
        msg = f"The chosen final interpolator (FinalInterp), {FinalInterp}, is"\
              + " not one of the accepted inputs: 'Linear' or 'NearestNeighbor'."
        raise Exception(msg)
    
    
    if InitMethod in ['geometry', 'geometricalcenter']:
        CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
    elif InitMethod in ['moments', 'centerofgravity']:
        CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
    
    
    if False:
        interact(PlotFixMovImages, 
                 FixInd=(0,FixIm.GetSize()[2]-1), 
                 MovInd=(0,MovIm.GetSize()[2]-1), 
                 FixPixArr=fixed(sitk.GetArrayViewFromImage(FixIm)), 
                 MovPixArr=fixed(sitk.GetArrayViewFromImage(MovIm)));
    
    if Transform == 'rigid':
        transf = sitk.Euler3DTransform()
    elif Transform == 'affine':
        transf = sitk.AffineTransform(FixIm.GetDimension())
    elif Transform == 'bspline':
        """ 04/05/21: The variable transf is never used for bspline... """
        transf = sitk.BSplineTransform(FixIm.GetDimension())
    
    if Transform in ['rigid', 'affine']:
        #SamplingPercentage = 0.01
        SamplingPercentage = 0.05
        
        LearningRate = 1.0
        #LearningRate = 5.0
        
    elif Transform == 'bspline':
        #SamplingPercentage = 0.01
        SamplingPercentage = 0.05
        
        LearningRate = 1.0
        #LearningRate = 5.0
    
    if FinalInterp == 'Linear':
        finalinterp = sitk.sitkLinear
    elif FinalInterp == 'NearestNeighbor':
        finalinterp = sitk.sitkNearestNeighbor
        
    
    if Transform == 'bspline':
        #MeshSize = [2] * FixIm.GetDimension() # 16/04/21
        #MeshSize = [8] * MovIm.GetDimension() # 16/04/21
        
        #InitialTx = sitk.BSplineTransformInitializer(FixIm, 
        #                                             MeshSize)
        
        """ Following example here:
        https://simpleitk.org/SPIE2018_COURSE/advanced_registration.pdf 
        """
        # Grid physical spacing:
        GridSpacing = [50.0, 50.0, 50.0] # A control point every 50 mm
        #GridSpacing = [25.0, 25.0, 25.0]
        #GridSpacing = [8.0, 8.0, 8.0] # same values that Elastix uses?
        
        print(f'GridSpacing = {GridSpacing}\n')
        
        # Image physical size:
        ImPhySize = [size*spacing for size, spacing in zip(FixIm.GetSize(), 
                                                           FixIm.GetSpacing())]
    
        print(f'ImPhySize = {ImPhySize}\n')
        
        # Grid mesh size:
        MeshSize = [int(sz/sp + 0.5) for sz, sp in zip(ImPhySize, GridSpacing)]
        
        # The starting mesh size will be 1/4 of the original, it will be  
        # refined by the multi-resolution framework.
        # https://github.com/SimpleITK/TUTORIAL/blob/master/06_advanced_registration.ipynb
        MeshSize = [int(sz/4 + 0.5) for sz in MeshSize] # 04/05/21
    
        print(f'MeshSize = {MeshSize}\n')
        
        if False:
            """ Simple: Use 8 control points. """
            #MeshSize = [8] * MovIm.GetDimension()
            #MeshSize = [3] * MovIm.GetDimension()
            MeshSize = [NumControlPts] * MovIm.GetDimension()
            
            print(f'MeshSize = {MeshSize}\n')
            
            
    """ Initial alignment. """   
    if InitMethod == 'landmarks':
        #InitialTx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
        #                                                   FixPts_flat, 
        #                                                   MovPts_flat)
        
        #print(f'FixFidsFpath = {FixFidsFpath}')
        #print(f'MovFidsFpath = {MovFidsFpath}')
        #print(f'ImageJinds = {ImageJinds}')
        
        AlignedIm, InitialTx = AlignImagesUsingLandmarks(FixIm, MovIm, 
                                                         FixFidsFpath, 
                                                         MovFidsFpath,
                                                         FlipK)
    else:
        if Transform in ['rigid', 'affine']:
            InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                          transf,
                                                          CTIF)
        elif Transform == 'bspline':
            InitialTx\
            = sitk.BSplineTransformInitializer(image1=FixIm, 
                                               transformDomainMeshSize=MeshSize, 
                                               order=3)
    
    
    #ResMovIm = sitk.Resample(MovIm, FixIm, InitialTx, sitk.sitkLinear, 0.0, 
    #                         MovIm.GetPixelID())
    #
    #
    #if False:
    #    interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2]-1), 
    #             alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(ResMovIm));
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    if Transform in ['rigid', 'affine']:
        #RegMethod.SetInitialTransform(InitialTx)
        # Don't optimize in-place, we would possibly like to run this cell multiple 
        # times.
        RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    elif Transform == 'bspline':
        #ScaleFactors = [1,2,5]
        ScaleFactors = [1,2,4] # 04/05/21
        
        RegMethod.SetInitialTransformAsBSpline(InitialTx, inPlace=False,
                                               scaleFactors=ScaleFactors) # 28/04/21

    """ Similarity metric settings. """
    if Transform in ['rigid', 'affine']:
        RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        #RegMethod.SetMetricAsJointHistogramMutualInformation()
    elif Transform == 'bspline':
        RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        #RegMethod.SetMetricAsMeanSquares()
    
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    #RegMethod.SetMetricSamplingPercentage(0.01)
    #RegMethod.SetMetricSamplingPercentage(0.05)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    if Transform in ['rigid', 'affine']:
        ShrinkFactors = [4,2,1]
        SmoothingSigmas = [2,1,1]
        
    elif Transform == 'bspline':
        ##ShrinkFactors = [6,2,1]
        ##SmoothingSigmas = [6,2,1]
        #ShrinkFactors = [6,4,2,1]
        #SmoothingSigmas = [4,2,1,0]
        ShrinkFactors = [4,2,1] # 28/04/21
        #SmoothingSigmas = [4,2,1] # 28/04/21
        SmoothingSigmas = [2,1,1] # 30/04/21
        
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    #if Transform in ['rigid', 'affine']:
    #    RegMethod.SetInterpolator(sitk.sitkLinear)
    #elif Transform == 'bspline':
    #    RegMethod.SetInterpolator(sitk.sitkBSpline)
    
    """ Optimizer settings. """
    if Transform in ['rigid', 'affine']:
        #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
        #                                        numberOfIterations=500, 
        #                                        convergenceMinimumValue=1e-6, 
        #                                        convergenceWindowSize=10)
        RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                                numberOfIterations=500, 
                                                estimateLearningRate=RegMethod.Once)
    elif Transform == 'bspline':
        if False:
            RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
                                                              numberOfIterations=100,
                                                              convergenceMinimumValue=1e-4,
                                                              #convergenceMinimumValue=1e-6,
                                                              convergenceWindowSize=5)
        if False:
            RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                                    numberOfIterations=500, 
                                                    estimateLearningRate=RegMethod.Once)
        # 04/05/21:
        # https://github.com/SimpleITK/TUTORIAL/blob/master/06_advanced_registration.ipynb
        RegMethod.SetOptimizerAsLBFGS2(solutionAccuracy=1e-2, 
                                       numberOfIterations=100, 
                                       deltaConvergenceTolerance=0.01)

        
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    
        
    # Connect all of the observers for plotting during registration:
    RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
    RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
    RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
    RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    """ 28/04/21: Should the MovIm be registered or AlignedIm 
    (for InitMethod = 'landmarks')? """
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print("-------")
    print('InitialTx:')
    print(InitialTx)
    print("-------")
    print("-------")
    print('FinalTx:')
    print(FinalTx)
    print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print(f"Optimizer stopping condition: {RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, finalinterp, 0.0, 
                          MovIm.GetPixelID())
    
    if False:
        interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(RegIm));

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    
    """ Create a dictionary to store the Parameters and FixedParameters, 
    similar to the TxParamMapDict returned with RegisterImagesSelx: """
    TxParamsStr = tuple([str(item) for item in FinalTx.GetParameters()])
    
    TxParamMapDict = {'Transform' : (Transform,),
                      'NumberOfParameters' : (f'{len(FinalTx.GetParameters())}'),
                      'TransformParameters' : TxParamsStr}
    
    if Transform == 'bspline':
        """ 22/04/21:  Need to confirm: """
        TxParamMapDict.update({'GridSize' : tuple([str(item) for item in MeshSize]),
                               'GridSpacing' : tuple([str(item) for item in FinalTx.GetFixedParameters()[0:3]]),
                               'GridOrigin' : tuple([str(item) for item in FinalTx.GetFixedParameters()[3:6]]),
                               """ or 
                               'GridOrigin' : tuple([str(item) for item in FinalTx.GetFixedParameters()[6:9]]),
                               """
                               'GridDirection' : tuple([str(item) for item in FinalTx.GetFixedParameters()[9:]])})
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, FinalTx


def RegisterImagesSitk_OLD2(FixIm, MovIm, Transform='affine', 
                       InitMethod='moments',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       LogToConsole=False):
    """
    Wrapper function for functions SitkNonDefReg(), SitkBsplineReg() and 
    SitkBsplineRegWithFiducials() to register two 3D SimpleITK images using 
    SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    #FinalInterp : string (optional; 'Linear' by default)
    #    The type of final interpolation to be applied.  Acceptable inputs
    #    include:
    #        - 'Linear'
    #        - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    #FlipK : boolean (optional; False by default)
    #    Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
    #    are indices obtained from ImageJ. This argument is only relevant if 
    #    InitMethod = 'landmarks'.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    LandmarkTx : SimpleITK Transform or None 
        The transform used to align MovIm to FixIm prior to registration based
        on fiducials.  If Fiducials were not used to initialise the regitration
        LandmarkTx = None.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    """
    
    if Transform == 'bspline' and InitMethod == 'landmarks':
        RegIm, RegMethod, LandmarkTx, InitialTx,\
        FinalTx = SitkBsplineRegWithFiducials(FixIm, MovIm, 
                                              NumControlPts=8, 
                                              LearningRate=5.0, NumIters=100, 
                                              FixFidsFpath=FixFidsFpath, 
                                              MovFidsFpath=MovFidsFpath, 
                                              LogToConsole=LogToConsole)
    
    elif Transform == 'bspline' and InitMethod != 'landmarks':
        RegIm, RegMethod, InitialTx,\
        FinalTx = SitkBsplineReg(FixIm, MovIm, NumControlPts=8, 
                                 LearningRate=5.0, NumIters=100, 
                                 LogToConsole=LogToConsole)
        
        LandmarkTx = None
    else:
        RegIm, RegMethod, InitialTx,\
        FinalTx = SitkNonDefReg(FixIm, MovIm, Transform=Transform, 
                                InitMethod=InitMethod,
                                FixFidsFpath=FixFidsFpath, 
                                MovFidsFpath=MovFidsFpath, 
                                NumIters=500, LogToConsole=LogToConsole)
        
        """ Even though InitialTx might be a landmark transform (if InitMethod
        = 'landmarks'), it won't be needed. """
        LandmarkTx = None
    
    return RegIm, LandmarkTx, FinalTx





def RegisterImagesSitk(FixIm, MovIm, Transform='affine', 
                       InitMethod='moments',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       SamplingPercentage=5, LogToConsole=False):
    """
    08/07/21:
        See register_ims in image_tools.registering.py
    
    Wrapper function for functions SitkNonDefReg(), SitkBsplineReg() and 
    SitkBsplineRegWithFiducials() to register two 3D SimpleITK images using 
    SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    #FinalInterp : string (optional; 'Linear' by default)
    #    The type of final interpolation to be applied.  Acceptable inputs
    #    include:
    #        - 'Linear'
    #        - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    #FlipK : boolean (optional; False by default)
    #    Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
    #    are indices obtained from ImageJ. This argument is only relevant if 
    #    InitMethod = 'landmarks'.
    
    SamplingPercentage : int or float (optional; 5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
    
    ListOfSitkTxs : list of SimpleITK Transforms
        A list of SimpleITK Transforms used to align MovIm to FixIm. The list
        may be of length 1, which will be the final transform, and will be for 
        most cases. If Transform = 'bspline' and InitMethod = 'landmarks', the 
        first Transform will be the transform used to align the images (based 
        on the fiducials) prior to registration, and the second transform will 
        be the final BSpline Transform.
    
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Possible ideas from:
    https://discourse.itk.org/t/segmentation-fault-error-when-running-rigid-affine-registration/3879
    """
    
    ListOfSitkTxs = []
    
    if Transform == 'bspline' and InitMethod == 'landmarks':
        if LogToConsole:
            print('Running SitkBsplineRegWithFiducials()...\n')
            
        RegIm, RegMethod, LandmarkTx, InitialTx,\
        FinalTx = SitkBsplineRegWithFiducials(FixIm=FixIm, MovIm=MovIm, 
                                              NumControlPts=8, 
                                              LearningRate=5.0, NumIters=100, 
                                              FixFidsFpath=FixFidsFpath, 
                                              MovFidsFpath=MovFidsFpath, 
                                              SamplingPercentage=SamplingPercentage,
                                              LogToConsole=LogToConsole)
        
        """ 01/07/21: Work-around: """
        ListOfSitkTxs = [LandmarkTx, FinalTx]
    
    elif Transform == 'bspline' and InitMethod != 'landmarks':
        if LogToConsole:
            print('Running SitkBsplineReg()...\n')
        
        RegIm, RegMethod, InitialTx,\
        FinalTx = SitkBsplineReg(FixIm=FixIm, MovIm=MovIm, NumControlPts=8, 
                                 LearningRate=5.0, 
                                 SamplingPercentage=SamplingPercentage, 
                                 NumIters=100, LogToConsole=LogToConsole)
        
        #LandmarkTx = None
        ListOfSitkTxs = [FinalTx]
        
    elif Transform in ['rigid', 'affine'] and InitMethod == 'landmarks':
        if LogToConsole:
            print('Running SitkNonDefRegWithFiducials()...\n')
        
        RegIm, RegMethod, LandmarkTx, OptimisedTx,\
        FinalTx = SitkNonDefRegWithFiducials(FixIm=FixIm, MovIm=MovIm, 
                                             FixFidsFpath=FixFidsFpath, 
                                             MovFidsFpath=MovFidsFpath,
                                             RegTxName=Transform, 
                                             SamplingPercentage=SamplingPercentage,
                                             NumIters=500, 
                                             LogToConsole=LogToConsole)
        
        """ Don't need to add LandmarkTx and FinalTx since 
        SitkNonDefRegWithFiducials uses ITKv4 method. """
        ListOfSitkTxs = [FinalTx]
    
    elif Transform in ['rigid', 'affine'] and InitMethod != 'landmarks':
        if LogToConsole:
            print('Running SitkNonDefReg()...\n')
        
        RegIm, RegMethod, InitialTx,\
        FinalTx = SitkNonDefReg(FixIm=FixIm, MovIm=MovIm, RegTxName=Transform, 
                                InitMethod=InitMethod,
                                SamplingPercentage=SamplingPercentage,
                                NumIters=500, LogToConsole=LogToConsole)
        
        #""" Even though InitialTx might be a landmark transform (if InitMethod
        #= 'landmarks'), it won't be needed. """
        #LandmarkTx = None
        ListOfSitkTxs = [FinalTx]
        
        LandmarkTx = None # 01/09/21
    
    return RegIm, LandmarkTx, FinalTx # 01/09/21
    #return RegIm, ListOfSitkTxs # 01/09/21





def command_iteration(method):
    print(f"{method.GetOptimizerIteration():3} = {method.GetMetricValue():10.5f}")
    print("\t#: ", len(method.GetOptimizerPosition()))


def command_multi_iteration(method):
    print("--------- Resolution Changing ---------")




def SitkAffineReg(FixIm, MovIm, InitMethod='moments',
                  FixFidsFpath='fixed_fiducials.txt', 
                  MovFidsFpath='moving_fiducials.txt', 
                  FlipK=False, LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    FinalInterp : string (optional; 'Linear' by default)
        The type of final interpolation to be applied.  Acceptable inputs
        include:
            - 'Linear'
            - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ. This argument is only relevant if 
        InitMethod = 'landmarks'.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
        
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    """
    
    import SimpleITK as sitk
    import time
    import winsound
    
    #import matplotlib.pyplot as plt
    from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    if not InitMethod in ['geometry', 'geometricalcenter', 'moments', 
                          'centerofgravity', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
    
    if InitMethod in ['geometry', 'geometricalcenter']:
        CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
    elif InitMethod in ['moments', 'centerofgravity']:
        CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
    
    #SamplingPercentage = 0.01
    SamplingPercentage = 0.05
    
    LearningRate = 1.0
    #LearningRate = 5.0
    
    SitkTx = sitk.AffineTransform(FixIm.GetDimension())
    
    """ Initial alignment. """   
    if InitMethod == 'landmarks':
        #InitialTx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
        #                                                   FixPts_flat, 
        #                                                   MovPts_flat)
        
        #print(f'FixFidsFpath = {FixFidsFpath}')
        #print(f'MovFidsFpath = {MovFidsFpath}')
        #print(f'ImageJinds = {ImageJinds}')
        
        AlignedIm, InitialTx = AlignImagesUsingLandmarks(FixIm, MovIm, 
                                                         FixFidsFpath, 
                                                         MovFidsFpath,
                                                         FlipK)
    else:
        InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                      SitkTx, CTIF)
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)

    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    ShrinkFactors = [4,2,1]
    SmoothingSigmas = [2,1,1]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
    #                                        numberOfIterations=500, 
    #                                        convergenceMinimumValue=1e-6, 
    #                                        convergenceWindowSize=10)
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=500, 
                                            estimateLearningRate=RegMethod.Once)
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  #convergenceMinimumValue=1e-6,
    #                                                  convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    
    
    # Connect all of the observers for plotting during registration:
    RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
    RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
    RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
    RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    """ 28/04/21: Should the MovIm be registered or AlignedIm 
    (for InitMethod = 'landmarks')? """
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print("-------")
    print('InitialTx:')
    print(InitialTx)
    print("-------")
    print("-------")
    print('FinalTx:')
    print(FinalTx)
    print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print(f"Optimizer stopping condition: {RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                          MovIm.GetPixelID())
    
    if False:
        interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(RegIm));

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    
    """ Create a dictionary to store the Parameters and FixedParameters, 
    similar to the TxParamMapDict returned with RegisterImagesSelx: """
    TxParamsStr = tuple([str(item) for item in FinalTx.GetParameters()])
    
    TxParamMapDict = {'Transform' : ('Affine',),
                      'NumberOfParameters' : (f'{len(FinalTx.GetParameters())}'),
                      'TransformParameters' : TxParamsStr}
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, FinalTx, TxParamMapDict




def SitkNonDefReg_OLD(FixIm, MovIm, Transform='affine', InitMethod='moments',
                  FixFidsFpath='fixed_fiducials.txt', 
                  MovFidsFpath='moving_fiducials.txt', 
                  FlipK=False, SamplingPercentage=5, NumIters=500,
                  LogToConsole=False):
    """
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    Transform : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    FinalInterp : string (optional; 'Linear' by default)
        The type of final interpolation to be applied.  Acceptable inputs
        include:
            - 'Linear'
            - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ. This argument is only relevant if 
        InitMethod = 'landmarks'.
    
    SamplingPercentage : int or float (optional; 5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    NumIters : integer (optional; 500 by default)
        The maximum number of iterations used when optimising.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK Image
        The 3D registered image.
    
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
     
    InitialTx : SimpleITK Transform
        The transform used to align MovIm to FixIm prior to registration.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    """
    
    import SimpleITK as sitk
    import time
    import winsound
    from importlib import reload
    import ConversionTools
    reload(ConversionTools)
    from ConversionTools import Image2BinaryOnesImage
    
    #import matplotlib.pyplot as plt
    from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple
    #PlotToConsole = True
    
    if not InitMethod in ['geometry', 'geometricalcenter', 'moments', 
                          'centerofgravity', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
    
    if InitMethod in ['geometry', 'geometricalcenter']:
        CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
    elif InitMethod in ['moments', 'centerofgravity']:
        CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
    
    #SamplingPercentage = 0.01
    #SamplingPercentage = 0.05 # takes 28.2 s
    #SamplingPercentage = 0.1
    #SamplingPercentage = 0.2
    #SamplingPercentage = 0.2
    
    SamplingPercentage = SamplingPercentage/100
    
    LearningRate = 1.0
    #LearningRate = 5.0
    
    if Transform == 'rigid':
        SitkTx = sitk.Euler3DTransform()
    elif Transform == 'affine':
        SitkTx = sitk.AffineTransform(FixIm.GetDimension())
    else:
        msg = f"Transform must be either 'rigid' or 'affine'. {Transform} is "\
              + "not an accepted value."
        raise Exception(msg)
    
    """ Initial alignment. """   
    if InitMethod == 'landmarks':
        #InitialTx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
        #                                                   FixPts_flat, 
        #                                                   MovPts_flat)
        
        #print(f'FixFidsFpath = {FixFidsFpath}')
        #print(f'MovFidsFpath = {MovFidsFpath}')
        #print(f'ImageJinds = {ImageJinds}')
        
        AlignedIm, InitialTx = AlignImagesUsingLandmarks(FixIm, MovIm, 
                                                         FixFidsFpath, 
                                                         MovFidsFpath,
                                                         FlipK=FlipK,
                                                         LogToConsole=LogToConsole)
    else:
        InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                      SitkTx, CTIF)
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)

    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    #RegMethod.SetMetricSamplingStrategy(RegMethod.REGULAR)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    ShrinkFactors = [4,2,1]
    SmoothingSigmas = [2,1,1]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
    #                                        numberOfIterations=500, 
    #                                        convergenceMinimumValue=1e-6, 
    #                                        convergenceWindowSize=10)
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  #convergenceMinimumValue=1e-6,
    #                                                  convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    
    if LogToConsole:
        if PlotToConsole:
            # Connect all of the observers for plotting during registration:
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    
    #FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ 28/04/21: Should the MovIm be registered or AlignedIm 
    (for InitMethod = 'landmarks')? 
    
    14/06/21: This might be the way to go but not without confining the 
    (randomly selected) sample points to points common to both FixIm and 
    AlignedIm using masks, since without doing so results in error:
    
    itk::ERROR: itk::ERROR: MattesMutualInformationImageToImageMetricv4: 
    All samples map outside moving image buffer. The images do not sufficiently 
    overlap. They need to be initialized to have more overlap before this metric
    will work. For instance, you can align the image centers by translation.
    """
    if InitMethod == 'landmarks':
        # Create a ones mask for FixIm:
        #FixMaskIm = Image2BinaryOnesImage(FixIm)
        #MovMaskIm = Image2BinaryOnesImage(MovIm)
        
        #RegMethod.SetMetricMovingMask(FixMaskIm)
        #RegMethod.SetMetricFixedMask(FixMaskIm)
        
        # https://discourse.itk.org/t/simpleitk-image-registration-4-13-vs-5-1/3649:
        # Recasting avoids problems which can occur for some images
        #RegMethod.SetMetricMovingMask(sitk.Cast(MovMaskIm,
        #                                        FixMaskIm.GetPixelIDValue()))
        #RegMethod.SetMetricMovingMask(sitk.Cast(FixMaskIm,
        #                                        FixMaskIm.GetPixelIDValue()))
        #RegMethod.SetMetricFixedMask(FixMaskIm)
        
        FinalTx = RegMethod.Execute(FixIm, AlignedIm)
    else:
        FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print("Optimizer stopping condition:",
          f"{RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    #RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
    #                      MovIm.GetPixelID())
    
    if InitMethod == 'landmarks':
        RegIm = sitk.Resample(AlignedIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                              AlignedIm.GetPixelID())
    else:
        RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                              MovIm.GetPixelID())
    
    if False:
        interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(RegIm));

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    
    #""" Create a dictionary to store the Parameters and FixedParameters: """
    #TxParamsStr = tuple([str(item) for item in FinalTx.GetParameters()])
    #
    #TxParamMapDict = {'Transform' : ('Affine',),
    #                  'NumberOfParameters' : (f'{len(FinalTx.GetParameters())}'),
    #                  'TransformParameters' : TxParamsStr}
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, RegMethod, InitialTx, FinalTx



def SitkNonDefReg(FixIm, MovIm, RegTxName='affine', InitTx=None, 
                  InitMethod='moments', SamplingPercentage=50, 
                  NumIters=500, LogToConsole=False):
    """
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    RegTxName : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
    
    InitTx : SimpleITK Transform (optional; None by default)
        If provided InitTx will be used to initialise registration. If None,
        initialisation will be using CenteredTransformInitializerFilter with
        the method defined by InitMethod.
    
    InitMethod : string (optional; 'moments' by default)
        The method used for initialisation.  This input will be ignored if
        InitTx is not None.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
    
    FinalInterp : string (optional; 'Linear' by default)
        The type of final interpolation to be applied.  Acceptable inputs
        include:
            - 'Linear'
            - 'NearestNeighbor'
    
    SamplingPercentage : int or float (optional; 50 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    NumIters : integer (optional; 500 by default)
        The maximum number of iterations used when optimising.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK Image
        The 3D registered image.
    
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
     
    InitialTx : SimpleITK Transform
        The transform used to align MovIm to FixIm prior to registration.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    It seems that a higher number of samples (e.g. 50%) is required to get
    results that are not so dependent upon the random sampling.
    
    It's possible to avoid randomness in the random selection of samples for
    repeated registrations by setting the seed in 
    
    SetMetricSamplingPercentage(double percentage, 
                                unsigned int seed = sitkWallClock)
    
    to a constant. i.e.:
    
    seed = 2
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage, seed)
    
    But although the repeated registered images will look the
    same they won't be any better than any other.  It's also not clear what, if
    any, suitable choice in the seed value. 	
    """
    
    import SimpleITK as sitk
    import numpy as np
    import time
    import winsound
    
    #import matplotlib.pyplot as plt
    from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple
    PlotToConsole = True
    
    if not InitMethod in ['geometry', 'geometricalcenter', 'moments', 
                          'centerofgravity', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
    
    if InitMethod in ['geometry', 'geometricalcenter']:
        CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
    elif InitMethod in ['moments', 'centerofgravity']:
        CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
    
    SamplingPercentage = SamplingPercentage/100
    
    LearningRate = 1.0
    #LearningRate = 5.0
    
    if RegTxName == 'rigid':
        SitkTx = sitk.Euler3DTransform()
    elif RegTxName == 'affine':
        SitkTx = sitk.AffineTransform(FixIm.GetDimension())
    else:
        msg = f"RegTxName must be either 'rigid' or 'affine'. {RegTxName} is "\
              + "not an accepted value."
        raise Exception(msg)
    
    """ Initial alignment. """
    if InitTx == None:
        InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                      SitkTx, CTIF)
    else:
        InitialTx = InitTx
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)

    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    #RegMethod.SetMetricSamplingStrategy(RegMethod.REGULAR)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    ShrinkFactors = [4,2,1]
    SmoothingSigmas = [2,1,1]
    #SmoothingSigmas = [2,1,0]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
    #                                        numberOfIterations=500, 
    #                                        convergenceMinimumValue=1e-6, 
    #                                        convergenceWindowSize=10)
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  #convergenceMinimumValue=1e-6,
    #                                                  convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    
    if LogToConsole:
        if PlotToConsole:
            # Connect all of the observers for plotting during registration:
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                                 UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    # Get the number of valid points:
    NumOfValidPts = RegMethod.GetMetricNumberOfValidPoints()
    
    # Get the number of points in FixIm:
    NumOfPtsInFixIm = np.prod(FixIm.GetSize())
    
    # The % of valid points:
    PercentageValidPts = 100*NumOfValidPts/NumOfPtsInFixIm
    
    print(f'NumOfValidPts = {NumOfValidPts} = {PercentageValidPts}%\n')
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print("Optimizer stopping condition:",
          f"{RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                          MovIm.GetPixelID())
    
    
    if False:
        interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(RegIm));

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, RegMethod, InitialTx, FinalTx



def SitkNonDefRegWithFiducials(FixIm, MovIm, 
                               FixFidsFpath='fixed_fiducials.txt', 
                               MovFidsFpath='moving_fiducials.txt', 
                               FlipK=False, 
                               RegTxName='affine',
                               SamplingPercentage=50, NumIters=500,
                               LogToConsole=False):
    """
    16/06/21
    
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    RegTxName : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
    
    #InitMethod : string (optional; 'geometry' by default)
    #    The method used for initialisation.  Acceptable inputs include:
    #        - 'geometry' (or 'geometricalcenter')
    #        - 'moments' (or 'centerofgravity')
    
    #FinalInterp : string (optional; 'Linear' by default)
    #    The type of final interpolation to be applied.  Acceptable inputs
    #    include:
    #        - 'Linear'
    #        - 'NearestNeighbor'
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    FlipK : boolean (optional; False by default)
        Set to True if the fiducials specified in FixFidsFpath and MovFidsFpath
        are indices obtained from ImageJ. This argument is only relevant if 
        InitMethod = 'landmarks'.
    
    SamplingPercentage : int or float (optional; 50 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    NumIters : integer (optional; 500 by default)
        The maximum number of iterations used when optimising.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK Image
        The 3D registered image.
    
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
     
    InitialTx : SimpleITK Transform
        The transform used to align MovIm to FixIm prior to registration.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Since initial alignment is geometrical (use of fiducials), initialising the
    registration with a geometrical centering is suitable (moments is not).
    
    It seems that a higher number of samples (e.g. 50%) is required to get
    results that are not so dependent upon the random sampling.
    """
    
    import SimpleITK as sitk
    import numpy as np
    import time
    import winsound
    
    #import matplotlib.pyplot as plt
    from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple
    PlotToConsole = True
    
    
    SamplingPercentage = SamplingPercentage/100
    
    LearningRate = 1.0
    #LearningRate = 5.0
    
    if RegTxName == 'rigid':
        OptimisedTx = sitk.Euler3DTransform()
    elif RegTxName == 'affine':
        OptimisedTx = sitk.AffineTransform(FixIm.GetDimension())
    else:
        msg = f"RegTxName must be either 'rigid' or 'affine'. {RegTxName} is "\
              + "not an accepted value."
        raise Exception(msg)
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Landmark alignment. """
    AlignedIm, LandmarkTx\
    = AlignImagesUsingLandmarks(FixIm=FixIm, MovIm=MovIm, 
                                FixFidsFpath=FixFidsFpath, 
                                MovFidsFpath=MovFidsFpath,
                                FlipK=FlipK,
                                Transform=RegTxName,
                                LogToConsole=LogToConsole)
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    #RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    
    # Set the initial moving and optimized transforms in ITKv4
    # https://discourse.itk.org/t/multiresolution-registration-with-2d-affine-transformation-on-pairs-of-2d-images/3096
    RegMethod.SetMovingInitialTransform(LandmarkTx)
    RegMethod.SetInitialTransform(OptimisedTx, inPlace=False)
    
    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    #RegMethod.SetMetricSamplingStrategy(RegMethod.REGULAR)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    ShrinkFactors = [4,2,1]
    SmoothingSigmas = [2,1,1]
    #SmoothingSigmas = [2,1,0]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
    #                                        numberOfIterations=500, 
    #                                        convergenceMinimumValue=1e-6, 
    #                                        convergenceWindowSize=10)
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  #convergenceMinimumValue=1e-6,
    #                                                  convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    
    if LogToConsole:
        if PlotToConsole:
            # Connect all of the observers for plotting during registration:
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                                 UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    # Get the number of valid points:
    NumOfValidPts = RegMethod.GetMetricNumberOfValidPoints()
    
    # Get the number of points in FixIm:
    NumOfPtsInFixIm = np.prod(FixIm.GetSize())
    
    # The % of valid points:
    PercentageValidPts = 100*NumOfValidPts/NumOfPtsInFixIm
    
    print(f'NumOfValidPts = {NumOfValidPts} = {PercentageValidPts}%\n')
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    
    #FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    # 16/06/21 Following:
    # https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/61_Registration_Introduction_Continued.ipynb
    #RegMethod.Execute(FixIm, MovIm)
    
    # Try this since Optimised is the default affine transform:
    OptimisedTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Compose the transformations: 
    https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/61_Registration_Introduction_Continued.ipynb
    """
    FinalTx = sitk.CompositeTransform(OptimisedTx)
    FinalTx.AddTransform(LandmarkTx)
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('LandmarkTx:')
        print(LandmarkTx)
        print("-------")
        print("-------")
        print('OptimisedTx:')
        print(OptimisedTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    MetricVal = RegMethod.GetMetricValue()
    StopConditionDesc = RegMethod.GetOptimizerStopConditionDescription()
    #FinalIterNum = RegMethod.GetOptimizerIteration()
    print(f'Final metric value: {MetricVal}')
    print(f"Optimizer stopping condition: {StopConditionDesc}")
    #print(f"Final Iteration: {FinalIterNum}")
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                          MovIm.GetPixelID())
    
    
    if False:
        interact(PlotBlendedImage, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=fixed(FixIm), ResIm=fixed(RegIm));

    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, RegMethod, LandmarkTx, OptimisedTx, FinalTx



def SitkBsplineReg(FixIm, MovIm, NumControlPts=8, LearningRate=5.0, 
                   SamplingPercentage=5, NumIters=100, LogToConsole=False):
    """ 
    07/07/21:
        See bspline_reg in image_tool.registering.py
    
    Register two 3D SimpleITK images using a deformable BSpline transformation
    using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    NumControlPts : integer (optional; 8 by default)
        The number of control points that defines the BSpline grid.
    
    LearningRate : float (optional; 5.0 by default)
        The learning rate used for the optimisation process (using a gradient
        descent approach).
    
    SamplingPercentage : int or float (optional; 5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    NumIters : integers (optional; 100 by default)
        The maximum number of iterations used during optimisation.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK Image
        The 3D registered image.
    
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
     
    InitialTx : SimpleITK Transform
        The transform used to align MovIm to FixIm prior to registration.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
    
    Unlike SimpleElastix, which allows for use of an "energy bending penalty"
    metric, SimpleITK lacks this. So BSplines with "too many" control points
    will result in non-realistic deformations.  8 seems to be a reasonable
    number.  Also, unlike SimpleElastix, which can perform a BSpline 
    registration with 32 control points over a reasonable time frame (several
    minutes), SimpleITK takes several hours.
    
    The parameters follow from the best result from 14/05/21 achieved in 
    20210423_Bspline_registrations.ipynb
    """
    import SimpleITK as sitk
    import time
    import winsound
    
    # Print simple metric values to console or plot of metric value v iteration?
    #PlotToConsole = False # simple
    PlotToConsole = True
    
    #def command_iteration(method):
    #    print(f"{method.GetOptimizerIteration():3} = {method.GetMetricValue():10.5f}")
    #    print("\t#: ", len(method.GetOptimizerPosition()))
    #
    #def command_multi_iteration(method):
    #    print("--------- Resolution Changing ---------")
    
    
    #SamplingPercentage = 0.01
    #SamplingPercentage = 0.05
    SamplingPercentage = SamplingPercentage/100
    
    """ Start timing: """
    times = []
    times.append(time.time())
        
    #MeshSize = [8] * MovIm.GetDimension()
    MeshSize = [NumControlPts] * MovIm.GetDimension()
    
    InitialTx = sitk.BSplineTransformInitializer(image1=FixIm,
                                                 transformDomainMeshSize=MeshSize)
    
    #print("Initial Parameters:")
    #print(InitialTx.GetParameters())
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, True)
    
    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    """ Change to original method: """
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
                                                      numberOfIterations=NumIters,
                                                      convergenceMinimumValue=1e-4,
                                                      convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    RegMethod.SetShrinkFactorsPerLevel([6, 2, 1])
    RegMethod.SetSmoothingSigmasPerLevel([6, 2, 1])
    
    if True:#LogToConsole:
        if PlotToConsole:
            # Connect all of the observers so that we can perform plotting  
            # during registration.
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
        print(f'Final metric value: {RegMethod.GetMetricValue()}')
        print("Optimizer stop condition:", 
              f"{RegMethod.GetOptimizerStopConditionDescription()}")
        #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(FixIm)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0.0)
    resampler.SetTransform(FinalTx)
    
    RegIm = resampler.Execute(MovIm)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    Dtime_min = round(Dtime/60, 1)
    print(f'*Took {Dtime} s ({Dtime_min} min) to perform image registration.\n')
    
    # Convert FinalTx from sitk Transform to sitk BSplineTransform:
    FinalTx = sitk.BSplineTransform(FinalTx)
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    # Took 2283 s = 38 min at iteration 18 (all samples)
    # Took 118.6 s = 2 min (1% samples)
    # Took 214 s = 3.6 min (5% samples)
    
    return RegIm, RegMethod, InitialTx, FinalTx



def SitkBsplineRegWithFiducials(FixIm, MovIm, NumControlPts=8, 
                                LearningRate=5.0, NumIters=100, 
                                FixFidsFpath='fixed_fiducials.txt', 
                                MovFidsFpath='moving_fiducials.txt', 
                                SamplingPercentage=5, LogToConsole=False):
    """ 
    07/07/21:
        See bspline_reg in image_tool.registering.py
    
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
    
    NumControlPts : integer (optional; 8 by default)
        The number of control points that defines the BSpline grid.
    
    LearningRate : float (optional; 5.0 by default)
        The learning rate used for the optimisation process (using a gradient
        descent approach).
    
    NumIters : integers (optional; 100 by default)
        The maximum number of iterations used during optimisation.
    
    FixFidsFpath : string (optional; 'fixed_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    MovFidsFpath : string (optional; 'moving_fiducials.txt' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    
    SamplingPercentage : int or float (optional; 5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK Image
        The 3D registered image.
    
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
    
    LandmarkTx : SimpleITK Transform
        The transform used to align MovIm to FixIm prior to registration based
        on fiducials.
    
    InitialTx : SimpleITK Transform
        The Bspline transform used to initialise a BSpline grid prior to
        registration.
    
    FinalTx : SimpleITK transform
        The final transform used to register/transform MovIm to FixIm.
    
    
    Notes:
    *****
     
    I was unable to successfully initialise the BSpline using a deformable
    use of LandmarkBasedTransformInitializer so instead I am using a 
    global alignment, and performing a BSpline registration on the aligned
    image. This is not ideal since interpolation errors stack up but it's
    a working solution for now. (26/05/21)
    
    The parameters follow from the best result from 14/05/21 achieved in 
    20210423_Bspline_registrations.ipynb
    """
    import SimpleITK as sitk
    #import time
    #import winsound
    
    
    """ Get the rigid landmark-based transform and aligned image: """
    AlignedIm, LandmarkTx = AlignImagesUsingLandmarks(Transform='affine',
                                                      FixFidsFpath=FixFidsFpath, 
                                                      MovFidsFpath=MovFidsFpath,
                                                      FixIm=FixIm, MovIm=MovIm, 
                                                      FlipK=False, Buffer=2,
                                                      LogToConsole=LogToConsole)
    
    print('LandmarkTx:\n')
    print(LandmarkTx, '\n')
    
    RegIm, RegMethod, InitialTx,\
    FinalTx = SitkBsplineReg(FixIm=FixIm, MovIm=AlignedIm, 
                             NumControlPts=NumControlPts, 
                             LearningRate=LearningRate, 
                             SamplingPercentage=SamplingPercentage, 
                             NumIters=NumIters, LogToConsole=LogToConsole)
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('LandmarkTx:')
        print(LandmarkTx)
        print("-------")
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print("Optimizer stop condition:",
          f"{RegMethod.GetOptimizerStopConditionDescription()}")
    #print(f"Final Iteration: {RegMethod.GetOptimizerIteration()}")
    
    
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(FixIm)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0.0)
    resampler.SetTransform(FinalTx)
    
    #RegIm = resampler.Execute(MovIm)
    RegIm = resampler.Execute(AlignedIm)
    
    
    # Convert FinalTx from sitk Transform to sitk BSplineTransform:
    #FinalTx = sitk.BSplineTransform(FinalTx)
    
    
    # Took 2283 s = 38 min at iteration 18 (all samples)
    # Took 118.6 s = 2 min (1% samples)
    # Took 214 s = 3.6 min (5% samples)
    
    #return RegIm, RegMethod, AlignTx, InitialBsplineTx, FinalBSplineTx
    #return RegIm, RegMethod, AlignTx, InitialBsplineTx, FinalBSplineTx, AlignedIm
    return RegIm, RegMethod, LandmarkTx, InitialTx, FinalTx


def Bspline_reg_using_fiducials():
    """ Following is the registration method that finally achieved a bspline reg
    using fiducials, following Zivy's help:
    
    https://github.com/SimpleITK/SimpleITK/issues/1415
    
    following example here: 
    
    https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/65_Registration_FFD.ipynb
    
    See 20210527_Test_CopyXnatRoiCol_and_RunCopyXnatRoiCol.ipynb
    
    roughly 3/5ths the way down as of 08/07/21.
    """
    import ImageTools
    reload(ImageTools)
    from ImageTools import AlignImagesUsingLandmarks
    from ImageTools import GetLandmarkTx, GetFiducialsWithinImExtents
    import registration_utilities as ru
    import registration_callbacks as rc
    import winsound
    
    NumOfControlPts = 8
    
    """ Use scale factors for mesh grid? """
    UseScaleFactors = False
    #UseScaleFactors = True
    
    """ Scale factors: """
    ScaleFactors = [1, 2, 4]
    
    """ Maybe the scale factors should take into account the NumOfControlPts
    that were used in LandmarkTx and the desired NumOfControlPts for the 
    registration?... """
    
    SamplingPercentage = 5/100
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Get landmark transform: """
    #LandmarkTx = GetLandmarkTx(Transform=AlignmentTx, FixFidsFpath=FixFidsFpath, 
    #                           MovFidsFpath=MovFidsFpath, FixIm=FixIm, 
    #                           MovIm=MovIm, FlipK=False, Buffer=Buffer,
    #                           LogToConsole=LogToConsole)
    
    """ Instead of only getting the landmark transform (above), also get the aligned 
    image (useful for debugging): """
    """ Something is wrong - geting 0 fiducials after running this:
    AlignedIm, LandmarkTx = AlignImagesUsingLandmarks(Transform='bspline', 
                                                      FixFidsFpath=TrgFidsFpath, 
                                                      MovFidsFpath=SrcFidsFpath,
                                                      FixIm=TrgIm, MovIm=TrgIm,
                                                      FlipK=False, Buffer=2,
                                                      LogToConsole=True)
    so will run the functions involved separately:
    """
    
    LandmarkTx, FixPts,\
    MovPts = GetLandmarkTx(Transform='bspline', 
                           FixFidsFpath=TrgFidsFpath, 
                           MovFidsFpath=SrcFidsFpath, 
                           FixIm=TrgIm, MovIm=SrcIm, FlipK=False,
                           Buffer=2, LogToConsole=True)
    
    print(f'LandmarkTx is {LandmarkTx.GetName()}\n')
    
    #""" Cast LandmarkTx to BSplineTransform as suggested by zivy:
    #https://github.com/SimpleITK/SimpleITK/issues/1415
    #"""
    #LandmarkTx = sitk.BSplineTransform(LandmarkTx)
    # above was added to GetLandmarkTx
    
    #print(f'LandmarkTx is now {LandmarkTx.GetName()} after casting as BSplineTransform\n')
    
    AlignedIm = sitk.Resample(SrcIm, TrgIm, LandmarkTx, sitk.sitkLinear, 0.0, 
                              SrcIm.GetPixelID())
    
    #print('LandmarkTx:\n')
    #print(LandmarkTx, '\n')
    
    dim = TrgIm.GetDimension()
    
    #MeshSize = [8] * MovIm.GetDimension()
    if UseScaleFactors:
        MeshSize = [NumOfControlPts // ScaleFactors[-1]] * dim
    else:
        MeshSize = [NumOfControlPts] * dim
    
    #OptimisedTx = sitk.BSplineTransformInitializer(image1=TrgIm,
    #                                               transformDomainMeshSize=MeshSize)
    
    #print(f'OptimisedTx is {OptimisedTx.GetName()}\n')
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    #RegMethod.SetInitialTransform(InitialTx, True)
    
    # Set the initial moving and optimized transforms in ITKv4
    # https://discourse.itk.org/t/multiresolution-registration-with-2d-affine-transformation-on-pairs-of-2d-images/3096
    #RegMethod.SetMovingInitialTransform(LandmarkTx)
    ###RegMethod.SetInitialTransform(OptimisedTx, inPlace=False)
    ##RegMethod.SetInitialTransform(OptimisedTx, inPlace=True)
    
    # https://simpleitk.readthedocs.io/en/v1.2.4/Examples/ImageRegistrationMethodBSpline3/Documentation.html
    if UseScaleFactors:
        #RegMethod.SetInitialTransformAsBSpline(OptimisedTx, inPlace=True,
        #                                       scaleFactors=ScaleFactors)
        RegMethod.SetInitialTransformAsBSpline(LandmarkTx, inPlace=True,
                                               scaleFactors=ScaleFactors)
    else:
        #RegMethod.SetInitialTransformAsBSpline(OptimisedTx, inPlace=True)
        #RegMethod.SetInitialTransform(OptimisedTx, inPlace=True)
        RegMethod.SetInitialTransform(LandmarkTx, inPlace=True)
    
    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    """ Change to original method: """
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    # Zivy advised not to use gradient descent line search for bspline:
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=1, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  convergenceWindowSize=5)
    
    # Instead use LBFGSB (without scale factors) or LBFGS2 (with scale factors):
    if UseScaleFactors:
        # Use the LBFGS2 instead of LBFGS. The latter cannot adapt to the changing control grid resolution.
        RegMethod.SetOptimizerAsLBFGS2(solutionAccuracy=1e-2, numberOfIterations=100, deltaConvergenceTolerance=0.01)
    else:
        RegMethod.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-5, numberOfIterations=100)
    
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    #RegMethod.SetShrinkFactorsPerLevel([6, 2, 1])
    #RegMethod.SetSmoothingSigmasPerLevel([6, 2, 1])
    
    RegMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    RegMethod.SetSmoothingSigmasPerLevel([4, 2, 1]) # <-- better result than [2, 1, 0]?..
    #RegMethod.SetSmoothingSigmasPerLevel([2, 1, 0]) # in 65_Registration_FFD.ipynb example
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn() # <-- THIS WAS MISSING until 21:51 on 1/7/21
    
    if True:#LogToConsole:
        if True:
            # Connect all of the observers so that we can perform plotting  
            # during registration.
            RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
            RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
        # Since corresponding points in the fixed and moving image are available  
        # display the similarity metric and the TRE during the registration:
        if False:
            RegMethod.AddCommand(sitk.sitkStartEvent, rc.metric_and_reference_start_plot)
            RegMethod.AddCommand(sitk.sitkEndEvent, rc.metric_and_reference_end_plot)
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: rc.metric_and_reference_plot_values(RegMethod, FixPts, MovPts))
            
    #FinalTx = RegMethod.Execute(FixIm, MovIm)
    FinalTx = RegMethod.Execute(TrgIm, SrcIm)
    
    # 16/06/21 Following:
    # https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/61_Registration_Introduction_Continued.ipynb
    #RegMethod.Execute(FixIm, MovIm)
    
    #""" Try casting OptimisedTx as Transform: """
    #OptimisedTx = sitk.Transform(OptimisedTx)
    #
    #print(f'OptimisedTx is now {OptimisedTx.GetName()} after casting as Transform\n')
    
    #OptimisedTx = RegMethod.Execute(TrgIm, SrcIm)
    
    #""" Compose the transformations: 
    #https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/61_Registration_Introduction_Continued.ipynb
    #"""
    #FinalTx = sitk.CompositeTransform(OptimisedTx)
    #FinalTx.AddTransform(LandmarkTx)
    
    # 17/06/21 Try: Convert FinalTx from sitk Transform to sitk BSplineTransform:
    #FinalTx = sitk.BSplineTransform(FinalTx)
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    Dtime_min = round((times[-1] - times[-2])/60, 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s ({Dtime_min} min) to perform image registration.\n')
    
    """ Post registration analysis. """
    if True:#LogToConsole:
        print("-------")
        print('LandmarkTx:')
        print(LandmarkTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print("Optimizer stopping condition:",
          f"{RegMethod.GetOptimizerStopConditionDescription()}")
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    
def SitkBsplineRegWithFiducials_NOT_WORKING():
    """ Following is method in development, trying to use ITKv4 framework to
    apply bspline registration using fiducials, as was achieved for affine
    registration. 
    
    At line OptimisedTx = RegMethod.Execute(TrgIm, SrcIm)
    
    get the error:
        
    RuntimeError: Exception thrown in SimpleITK ImageRegistrationMethod_Execute: 
    d:\a\1\sitk-build\itk-prefix\include\itk-5.1\itkImageToImageMetricv4GetValueAndDerivativeThreaderBase.hxx:268:
    Exception in GetValueAndDerivativeProcessPoint:
    d:\a\1\sitk-build\itk-prefix\include\itk-5.1\itkBSplineBaseTransform.h:302:
    itk::ERROR: itk::ERROR: BSplineTransform(0000027AA7FCA070): 
    ComputeJacobianWithRespectToPosition not yet implemented for BSplineTransform
    
    See:
    
    https://discourse.itk.org/t/using-landmarkbasedtransforminitializer-for-bspline-registrations/4096/8
    
    https://github.com/SimpleITK/SimpleITK/issues/1415
    """

def RegUsingFiducials_OLD(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
                      SamplingPercentage=0.05, LearningRate=1.0, 
                      NumIters=500, FinalInterpolator='linear'):

    import SimpleITK as sitk
    import time
    #import matplotlib.pyplot as plt
    #from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    if FinalInterpolator == 'linear':
        FinalInterp = sitk.sitkLinear
    elif FinalInterpolator == 'bspline':
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension())
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension(), 3)
        FinalInterp = sitk.sitkBSpline
    elif FinalInterpolator == 'nearestneighbor':
        FinalInterp = sitk.sitkNearestNeighbor
    else:
        msg = "FinalInterpolator must be either 'linear', 'bspline' or "\
              + "'nearestneighbor'."
        raise Exception(msg)
        
    InitialTx, FixPts, MovPts = GetLandmarkTx(FixIm, MovIm, 
                                              FixFidsFpath, MovFidsFpath)
    #InitialTx = GetLandmarkTx_v2(FixIm, MovIm, FixFidsFpath, MovFidsFpath)
    
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    
    RegMethod.SetInterpolator(sitk.sitkLinear)
    #RegMethod.SetInterpolator(sitk.sitkBSpline)

    """ Metric settings """
    #SamplingPercentage = 0.05
        
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    
    """ Optimizer settings """
    #LearningRate = 1.0
    #NumIters = 500
    
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    ShrinkFactors = [4,2,1]
    #SmoothingSigmas = [4,2,1] # worse?
    SmoothingSigmas = [2,1,1] # same as ImageTools code
    #SmoothingSigmas = [2,1,0] # possibly slightly worse?
    #ScaleFactors = [1,2,5]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
    RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
    RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
    RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print('Optimizer\'s stopping condition,',
          f'{RegMethod.GetOptimizerStopConditionDescription()}')
    
    #RegIm = ResampleImage_v2(MovIm, FixIm, ResInterp)
    #RegIm = ResampleImage_v2(MovIm, FixIm, FinalTx, ResInterp)
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, FinalInterp, 0.0, 
                          MovIm.GetPixelID())
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'*Took {Dtime} s to perform image registration.\n')
    
    return RegIm, FinalTx


def RegUsingFiducials_v2_OLD(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
                      AlignmentTransform='rigid', SamplingPercentage=0.05, 
                      LearningRate=1.0, NumIters=500, 
                      RegInterpolator='linear',
                      FinalInterpolator='linear'):

    import SimpleITK as sitk
    import time
    #import matplotlib.pyplot as plt
    #from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    if RegInterpolator == 'linear':
        RegInterp = sitk.sitkLinear
    elif RegInterpolator == 'bspline':
        #RegInterpTx = sitk.BSplineTransform(FixIm.GetDimension())
        #RegInterpTx = sitk.BSplineTransform(FixIm.GetDimension(), 3)
        RegInterp = sitk.sitkBSpline
    elif RegInterpolator == 'nearestneighbor':
        RegInterp = sitk.sitkNearestNeighbor
    else:
        msg = "RegInterpolator must be either 'linear', 'bspline' or "\
              + "'nearestneighbor'."
        raise Exception(msg)
    
    if FinalInterpolator == 'linear':
        FinalInterp = sitk.sitkLinear
    elif FinalInterpolator == 'bspline':
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension())
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension(), 3)
        FinalInterp = sitk.sitkBSpline
    elif FinalInterpolator == 'nearestneighbor':
        FinalInterp = sitk.sitkNearestNeighbor
    else:
        msg = "FinalInterpolator must be either 'linear', 'bspline' or "\
              + "'nearestneighbor'."
        raise Exception(msg)
        
    #InitialTx = GetLandmarkTx(FixIm, MovIm, FixFidsFpath, MovFidsFpath)
    InitialTx = GetLandmarkTx_v2(FixIm, MovIm, FixFidsFpath, MovFidsFpath,
                                 AlignmentTransform)
    
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    
    #RegMethod.SetInterpolator(sitk.sitkLinear)
    #RegMethod.SetInterpolator(sitk.sitkBSpline)
    RegMethod.SetInterpolator(RegInterp)

    """ Metric settings """
    #SamplingPercentage = 0.05
        
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    
    """ Optimizer settings """
    #LearningRate = 1.0
    #NumIters = 500
    
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    ShrinkFactors = [4,2,1]
    #SmoothingSigmas = [4,2,1] # worse?
    SmoothingSigmas = [2,1,1] # same as ImageTools code
    #SmoothingSigmas = [2,1,0] # possibly slightly worse?
    #ScaleFactors = [1,2,5]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
    RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
    RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
    RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print('Optimizer\'s stopping condition,',
          f'{RegMethod.GetOptimizerStopConditionDescription()}')
    
    #RegIm = ResampleImage(MovIm, FixIm, FinalInterpTx)
    #RegIm = ResampleImage(MovIm, FixIm, FinalTx, FinalInterpTx)
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, FinalInterp, 0.0, 
                          MovIm.GetPixelID())
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'*Took {Dtime} s to perform image registration.\n')
    
    return RegIm, FinalTx


def RegUsingFiducials_v3_OLD(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
                      AlignmentTransform='rigid', SamplingPercentage=0.05, 
                      LearningRate=1.0, NumIters=500, NumControlPts=8,
                      RegInterpolator='linear',
                      FinalInterpolator='linear'):

    import SimpleITK as sitk
    import time
    #import matplotlib.pyplot as plt
    #from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    if RegInterpolator == 'linear':
        RegInterp = sitk.sitkLinear
    elif RegInterpolator == 'bspline':
        #RegInterpTx = sitk.BSplineTransform(FixIm.GetDimension())
        #RegInterpTx = sitk.BSplineTransform(FixIm.GetDimension(), 3)
        RegInterp = sitk.sitkBSpline
    elif RegInterpolator == 'nearestneighbor':
        RegInterp = sitk.sitkNearestNeighbor
    else:
        msg = "RegInterpolator must be either 'linear', 'bspline' or "\
              + "'nearestneighbor'."
        raise Exception(msg)
    
    if FinalInterpolator == 'linear':
        FinalInterp = sitk.sitkLinear
    elif FinalInterpolator == 'bspline':
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension())
        #FinalInterp = sitk.BSplineTransform(FixIm.GetDimension(), 3)
        FinalInterp = sitk.sitkBSpline
    elif FinalInterpolator == 'nearestneighbor':
        FinalInterp = sitk.sitkNearestNeighbor
    else:
        msg = "FinalInterpolator must be either 'linear', 'bspline' or "\
              + "'nearestneighbor'."
        raise Exception(msg)
        
    #InitialTx = GetLandmarkTx(FixIm, MovIm, FixFidsFpath, MovFidsFpath)
    InitialTx = GetLandmarkTx_v2(FixIm, MovIm, FixFidsFpath, MovFidsFpath,
                                 AlignmentTransform)
    
    RegMethod = sitk.ImageRegistrationMethod()
    
    if RegInterpolator in ['rigid', 'affine']:
        RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    else:
        """ 
        Attempt 1:
        
        RegMethod.SetInitialTransformAsBSpline(InitialTx, inPlace=False,
                                               scaleFactors=[1, 2, 5])
        
        Error:
        TypeError: in method 'ImageRegistrationMethod_SetInitialTransformAsBSpline', 
        argument 2 of type 'itk::simple::BSplineTransform &'
        """
        
        """ 
        Attempt 2:
        Using the composite transformation we just add them in (stack based, 
        first in - last applied). 
        
        CompositeTx = sitk.CompositeTransform(sitk.BSplineTransform(FixIm.GetDimension()))
        CompositeTx.AddTransform(InitialTx)
        RegMethod.SetInitialTransform(CompositeTx, inPlace=False)
        
        Error:
        module 'SimpleITK' has no attribute 'CompositeTransform'
        """
        
        """ 
        Attempt 3:
        
        MeshSize = [NumControlPts]*FixIm.GetDimension()
        
        CompositeTx = sitk.BSplineTransformInitializer(image1=FixIm, 
                                                       transformDomainMeshSize=MeshSize, 
                                                       order=3)
        CompositeTx.AddTransform(InitialTx)
        
        RegMethod.SetInitialTransform(CompositeTx, inPlace=False)
        
        Error:
        Exception thrown in SimpleITK Transform_AddTransform: 
        C:\...\Common\src\sitkBSplineTransform.cxx:252:
        sitk::ERROR: Transform is not of type BSplineTransform!
        """
        
        # out of ideas... 29/04/21
        
    #RegMethod.SetInterpolator(sitk.sitkLinear)
    #RegMethod.SetInterpolator(sitk.sitkBSpline)
    RegMethod.SetInterpolator(RegInterp)

    """ Metric settings """
    #SamplingPercentage = 0.05
        
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    
    """ Optimizer settings """
    #LearningRate = 1.0
    #NumIters = 500
    
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    ShrinkFactors = [4,2,1]
    #SmoothingSigmas = [4,2,1] # worse?
    SmoothingSigmas = [2,1,1] # same as ImageTools code
    #SmoothingSigmas = [2,1,0] # possibly slightly worse?
    #ScaleFactors = [1,2,5]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    RegMethod.AddCommand(sitk.sitkStartEvent, StartPlot)
    RegMethod.AddCommand(sitk.sitkEndEvent, EndPlot)
    RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, UpdateMultiresIterations) 
    RegMethod.AddCommand(sitk.sitkIterationEvent, lambda: PlotValues(RegMethod))
    
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print('Optimizer\'s stopping condition,',
          f'{RegMethod.GetOptimizerStopConditionDescription()}')
    
    #RegIm = ResampleImage(MovIm, FixIm, FinalInterpTx)
    #RegIm = ResampleImage(MovIm, FixIm, FinalTx, FinalInterpTx)
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, FinalInterp, 0.0, 
                          MovIm.GetPixelID())
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'*Took {Dtime} s to perform image registration.\n')
    
    return RegIm, FinalTx


def TransformImageSitk(MovIm, FixIm, SitkTransform, Interp='Linear', 
                       LogToConsole=False):
    """
    08/07/21:
        This function was redundant. Use resample_im in 
        image_tools.resampling.py
    
    Transform a 3D SimpleITK image using SimpleITK.
    
    Inputs:
    ******
    
    MovIm : SimpleITK image 
        The 3D image to be transformed.
        
    FixIm : SimpleITK image
        The 3D image whose grid MovIm is to be transformed to.
        
    SitkTransform : SimpleITK Transform
        The transform (e.g. image registration transform) that maps Im to TrgIm.
        
    Interp : string (optional; 'Linear' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Linear'
        - 'NearestNeighbor'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    TxIm : SimpleITK image
        The 3D transformed image.
    """
    
    import SimpleITK as sitk
            
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of TransformImageSitk():')
        print('\n\n', '-'*120)
        
    
    if Interp == 'NearestNeighbor':
        #sitkInterp = sitk.sitkLinear # 13/05/21
        sitkInterp = sitk.sitkNearestNeighbor # 13/05/21
    else:
        #sitkInterp = sitk.sitkNearestNeighbor # 13/05/21
        sitkInterp = sitk.sitkLinear # 13/05/21
    
    
    TxIm = sitk.Resample(MovIm, FixIm, SitkTransform, sitkInterp)
    
    if LogToConsole:
        print('-'*120)
    
    return TxIm



def TransformImWithListOfSitkTxs(MovIm, FixIm, ListOfSitkTxs, Interp='Linear', 
                                 LogToConsole=False):
    """
    Note:
        This is a modification of TransformImageSitk() to get around issue
        of trying to initialise a BSpline transform with fiducials.
        
        Performing multiple resampling operations is not ideal due to stacking
        of errors. This is most problematic for nearest-neighbour interpolation
        so to minimise the impact, all but the last resampling step will be
        performed with a nearest-neighbour interpolator - others will be using
        a linear interpolator.
    
    Transform/resample a 3D SimpleITK image by iterating through a list of
    SimpleITK Transforms.
    
    Inputs:
    ******
    
    MovIm : SimpleITK image 
        The 3D image to be transformed.
        
    FixIm : SimpleITK image
        The 3D image whose grid MovIm is to be transformed to.
        
    ListOfSitkTxs : list of SimpleITK Transforms
        A list of SimpleITK Transforms.  The list may be of length 1.
        
    Interp : string (optional; 'Linear' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Linear'
        - 'NearestNeighbor'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    TxIm : SimpleITK image
        The 3D transformed image.
    """
    
    #import SimpleITK as sitk
    
    #OrigMovIm = sitk.CopyInformation(MovIm)
    
    NumOfTx = len(ListOfSitkTxs)
    
    for i in range(NumOfTx):
        SitkTx = ListOfSitkTxs[i]
        
        if Interp != 'Linear':
            if i < NumOfTx - 1:
                ResIm = TransformImageSitk(MovIm, FixIm, SitkTx, 'Linear', 
                                           LogToConsole)
            else:
                ResIm = TransformImageSitk(MovIm, FixIm, SitkTx, Interp, 
                                           LogToConsole)
        else:
            ResIm = TransformImageSitk(MovIm, FixIm, SitkTx, Interp, 
                                       LogToConsole)
        
        #MovIm = sitk.CopyInformation(ResIm)
        MovIm = CopyImage(ResIm)
    
    return ResIm


def TransformLabImByRoi(LabImByRoi, F2SindsByRoi, TrgIm, 
                        SelxImFiltOrSitkTx, Interp='NearestNeighbor', 
                        #ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
                        #ApplyPostTxBin=True, #VoxelVolumeRatio=1,
                        ApplyPostTxBlur=False, PostTxVariance=(1,1,1),
                        ApplyPostTxBin=False, #VoxelVolumeRatio=1,
                        LogToConsole=False):
    """
    08/07/21:
        See resample_labim in image_tools.resampling.py
        
    08/06/21: Default values of ApplyPostTxBlur and ApplyPostTxBin changed
    from True to False due to strang results obtained when testing on CT and MR
    data involving BSpline registration with fiducials alignment.
    
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    
    F2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the frame-to-slice 
        indices of each labelmap image in LabImByRoi.
    
    TrgIm : SimpleITK image
        The 3D image whose grid LabImByRoi is to be transformed to.
        
    SelxImFiltOrSitkTx : Elastix Image Filter or SimpleITK Transform
        The Elastix Image Filter used to perform image registration (if Elastix
        was used); or the final SimpleITK Transform if registration was 
        performed using SimpleITK.
    
    Interp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
        
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image(s) will be Gaussian blurred.
    
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image(s) will be binary thresholded.
    
    #ThreshPostTx : float (optional; 0.75 by default)
    #    The threshold level used to binarise the labelmap image following 
    #    transformation (or following transformation and Gaussian blurring).
    
    #VoxelVolumeRatio : float (optional; 1 by default)
    #    The ratio of the voxel spacings of the Target : Source images.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    TxLabImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabImByRoi
        transformed to the FixImage grid space.
        
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabImByRoi).
    
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    
    
    Notes:
    *****
    
    If the transformed binary labelmap image is blurred, the resulting labelmap
    image will be non-binary, so will need to be binarised. The threshold used
    will be determined such that the resulting (re-) binarised image will have
    a similar number of truth values (ones) as the pre-transformed (binary)
    one.
    """
    
    import time
    #import SimpleITK as sitk
    from ConversionTools import Image2PixArr
    from GeneralTools import UniqueItems
    from SelxTools import TransformImageSelx
    from copy import deepcopy
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running of TransformLabImByRoi():')
        print('\n\n', '-'*120)
    
            
    """ Start timing. """
    times = []
    times.append(time.time())
    
        
    """ Transform LabImByRoi to the Target image's grid to get the resampled
    pixel array. """
    TxLabImByRoi = []
    TxPixArrByRoi = []
    TxF2SindsByRoi = [] 
    
    for r in range(len(LabImByRoi)):
        LabIm = LabImByRoi[r]
        
        if LogToConsole:
            print('\n   Image info for LabIm prior to transformation:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(LabIm, LogToConsole)
        
        OrigF2Sinds = deepcopy(F2Sinds)
        
        """ Transform LabIm. """
        if 'ElastixImageFilter' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSelx(MovIm=LabIm, 
                                         SelxImFilt=SelxImFiltOrSitkTx,
                                         Interp=Interp, 
                                         LogToConsole=LogToConsole)
        elif 'Transform' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSitk(MovIm=LabIm, FixIm=TrgIm, 
                                         SitkTransform=SelxImFiltOrSitkTx, 
                                         Interp=Interp, 
                                         LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'*Took {Dtime} s to transform the {r}^th labelmap image',
              'using the image registration transformation.\n')
        
        if LogToConsole:
            print('\n   Image info for LabIm after transformation:')
            
        PixID, PixIDTypeAsStr, UniqueValsPostTx,\
        F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        if LogToConsole:
            from PlottingTools import PlotTwoImages
            
            PlotTwoImages(Im0=LabIm, Ind0=OrigF2Sinds[0], 
                          PlotLabel0='Original label image', 
                          Im1=TxLabIm, Ind1=F2Sinds[0], 
                          PlotLabel1='Resampled label image')
            
        
        if ApplyPostTxBlur:
            #""" The binary pixel array: """
            #PixArr_B, F2Sinds = Image2PixArr(TxLabIm)
            
            # 11/06/21:
            """ Ensure TxLabIm is a 32-bit float (PixID = 8). """
            if PixID != 8: 
                if LogToConsole:
                    print(f'\n   TxLabIm has PixelID = {PixID} '\
                          f'({PixIDTypeAsStr})). It will be converted to '\
                          'Float32.')
                    
                """ Convert TxLabIm to Float32: """
                TxLabIm = ChangeImagePixelType(Image=TxLabIm, 
                                               NewPixelType='Float32')
                
                if LogToConsole:
                    print('\n   Image info for TxLabIm after converting to Float32:')
                    
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
            
            
            """ Gaussian blur TxLabIm (new 14/02/21) to reduce noise that 
            makes selection of an appropriate threshold for binarisation (below)
            difficult. """
            #TxLabIm = GaussianBlurImage(TxLabIm)
            #TxLabIm = GaussianBlurImage(TxLabIm, Variance=PostTxVariance) # 11/06/21
            BlurTxLabIm = GaussianBlurImage(TxLabIm, Variance=PostTxVariance,
                                            LogToConsole=LogToConsole) # 11/06/21
                
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            print(f'*Took {Dtime} s to Gaussian blur the {r}^th labelmap image.\n')
            
            if LogToConsole:
                print('\n   Image info for BlurTxLabIm after Gaussian blur:')
                
            #PixID, PixIDTypeAsStr, UniqueValsPostBlur,\
            #F2Sinds = GetImageInfo(TxLabIm, LogToConsole) # 11/06/21
            PixID, PixIDTypeAsStr, UniqueValsPostBlur,\
            F2Sinds = GetImageInfo(BlurTxLabIm, LogToConsole) # 11/06/21
            
            #""" The non-binary pixel array following transformation: """
            #PixArr_NB, F2Sinds = Image2PixArr(TxLabIm)
            
        
        
            if ApplyPostTxBin:
                """Find suitable threshold value that approximately preserves
                the number of pre-transformed truth values: """
                #Thresh = FindSuitableThresh(BinaryIm=LabIm, 
                #                            NonBinaryIm=TxLabIm,
                #                            LogToConsole=LogToConsole) # 11/06/21
                Thresh = FindSuitableThresh(BinaryIm=TxLabIm, 
                                            NonBinaryIm=BlurTxLabIm,
                                            LogToConsole=LogToConsole) # 11/06/21
                
                
                """ Binary threshold BlurTxLabIm. """
                #ThreshValue = ThreshLevel*max(UniqueVals)
                #ThreshValue = ThreshLevel                
                
                #TxLabIm = BinaryThresholdImage(Im=TxLabIm, 
                #                               #Thresh=ThreshPostTx)
                #                               Thresh=Thresh)
                BinBlurTxLabIm = BinaryThresholdImage(Im=BlurTxLabIm, 
                                                      #Thresh=ThreshPostTx)
                                                      Thresh=Thresh)
                
                if LogToConsole:
                    print('\n   Image info for BinBlurTxLabIm after binary',
                          #f'thresholding at {ThreshPostTx}:')
                          f'thresholding at {Thresh}:')
                
                #PixID, PixIDTypeAsStr, UniqueVals,\
                #F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(BinBlurTxLabIm, LogToConsole)
                
                #if LogToConsole:
                #    print('\n   Image info for BinBlurTxLabIm after binary',
                #          'thresholding:')
                #    
                #PixID, PixIDTypeAsStr, UniqueVals,\
                #F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
                
                # 11/06/21:
                TxLabIm = BinBlurTxLabIm
            else:
                # 11/06/21:
                TxLabIm = BlurTxLabIm
                
        
        
        """ Ensure TxLabIm is a 32-bit unsigned integer (PixID = 5). """
        if PixID != 5: 
            if LogToConsole:
                print(f'\n   TxLabIm has PixelID = {PixID} '\
                      f'({PixIDTypeAsStr})).')
            
            """ Convert TxLabIm from float to 32-bit unsigned integer. """
            TxLabIm = ChangeImagePixelType(Image=TxLabIm, 
                                           NewPixelType='UInt32')
            
            if LogToConsole:
                print('\n   Image info for TxLabIm after converting to',
                      '32-bit unsigned int:')
                
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        
        
        if ApplyPostTxBlur and ApplyPostTxBin:
            """ Store the threshold used as metadata: """
            TxLabIm.SetMetaData("PostTxThreshUsed", f"{Thresh}")
            
            if LogToConsole:
                print(f'\nThe key "PostTxThreshUsed" with value "{Thresh}"',
                      'was added to the metadata of TxLabIm.')
            
        
        TxLabImByRoi.append(TxLabIm)
        
        
        """ Convert TxLabIm to a pixel array. """
        TxPixArr, TxF2Sinds = Image2PixArr(TxLabIm)
        
        TxPixArrByRoi.append(TxPixArr)
        TxF2SindsByRoi.append(TxF2Sinds)
            
        #PixID, PixIDTypeAsStr, UniqueVals,\
        #F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        if LogToConsole:
            print('\n   Info for conversion of TxLabIm to TxPixArr:')
            PixID, PixIDTypeAsStr, UniqueVals,\
            F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
            unique = UniqueItems(Items=TxPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArr',
                  f'after transforming LabImByRoi[{r}]')
            
            F = TxPixArr.shape[0]
            
            print(f'\nThe segmentation has been registered to {F} frames/',
                  f'contours: {TxF2Sinds}.')
            
        """ Check if there are fewer indices in TxF2Sinds than in 
        F2SindsByRoi[r]. """
        Nin = len(F2SindsByRoi[r])
        Nout = len(TxF2Sinds)
        
        if Nout < Nin and LogToConsole:
            print(f'\nWARNING:  There are {Nin} F2Sinds in the input labelmap',
                  f'but only {Nout} in the transformed labelmap.\n')
            
    if LogToConsole:
        print('-'*120)
        
    return TxLabImByRoi, TxPixArrByRoi, TxF2SindsByRoi



def TransformLabImByRoiWithListOfSitkTxs(LabImByRoi, F2SindsByRoi, TrgIm, 
                                         ListOfSitkTxs, Interp='NearestNeighbor', 
                                         #ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
                                         #ApplyPostTxBin=True, LogToConsole=False):
                                         ApplyPostTxBlur=False, PostTxVariance=(1,1,1),
                                         ApplyPostTxBin=False, LogToConsole=False):
                                       
    """
    08/06/21: Default values of ApplyPostTxBlur and ApplyPostTxBin changed
    from True to False due to strang results obtained when testing on CT and MR
    data involving BSpline registration with fiducials alignment.
    
    Note:
        This is a modification of TransformLabImByRoi() to get around issue
        of trying to initialise a BSpline transform with fiducials.
        
        Performing multiple resampling operations is not ideal due to stacking
        of errors.  This will likely be exacerbated due to multiple operations
        using a NearestNeighbor interpolation.  To minimise the impact, all but 
        the last resampling step will be performed with a nearest-neighbour 
        interpolator - others will be using a linear interpolator.
    
    09/06/21: Results using a linear interpolator for all but the final stage
    looked odd, so instead will by default use NN interpolator for all stages.
    But included the variable ApplySameInterpAllStages to toggle b/t these 2
    approaches.
    
    Resample a list 3D SimpleITK images representing labelmaps by iterating 
    through a list of SimpleITK Transforms. A NearestNeighbor interpolation 
    will be applied initially.  If aliasing effects result in a empty labelmap, 
    the labelmap will be Gaussian blurred, linearly resampled and binary 
    thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    
    F2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the frame-to-slice 
        indices of each labelmap image in LabImByRoi.
    
    TrgIm : SimpleITK image
        The 3D image whose grid LabImByRoi is to be transformed to.
        
    ListOfSitkTxs : list of SimpleITK Transforms
        A list of SimpleITK Transforms.  The list may be of length 1.
    
    Interp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
        
    ApplyPostTxBlur : boolean (optional; True by default)
        If True, the post-transformed labelmap image(s) will be Gaussian blurred.
    
    PostTxVariance : tuple of floats (optional; (1,1,1) by default)
        The variance along all dimensions if Gaussian blurring the post-
        transformed labelmap image(s).
        
    ApplyPostTxBin : boolean (optional; True by default)
        If True, the post-transformed (or post-transformed + Gaussian blurred)
        labelmap image(s) will be binary thresholded.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    TxLabImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabImByRoi
        transformed to the FixImage grid space.
        
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabImByRoi).
    
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    
    
    Notes:
    *****
    
    If the transformed binary labelmap image is blurred, the resulting labelmap
    image will be non-binary, so will need to be binarised. The threshold used
    will be determined such that the resulting (re-) binarised image will have
    a similar number of truth values (ones) as the pre-transformed (binary)
    one.
    
    09/06/21: Option to apply the same interpolator, Interp, for all stages,
    or to apply a linear interpolator for all but the final stage, and apply
    Interp for the final stage only.
    
    11/06/21: Option to apply Gaussian blurring (if ApplyPostTxBlur = True) and
    binarisation (if ApplyPostTxBin = True) at all stages or only the final one.
    
    11/06/21: Transformation of a binary label image with a local transform
    (e.g. BSpline) seems to result in more highly dislocated masks than global
    transforms (e.g. affine). So if ApplyPostTxBlur = True and the transform
    is BSpline set the minimum variance along any dimension to be 3.
    """
    
    from copy import deepcopy
    
    #ApplySameInterpAllStages = False
    ApplySameInterpAllStages = True
    
    #ApplyBlurAllStages = True
    ApplyBlurAllStages = False
    
    #ApplyBinAllStages = True
    ApplyBinAllStages = False
    
    MinVarForBspline = 3
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print('Running TransformLabImByRoiWithListOfSitkTxs():')
        print('\n\n', '-'*120)
        print(f"  ApplySameInterpAllStages = {ApplySameInterpAllStages}")
        print(f"  ApplyBlurAllStages = {ApplyBlurAllStages}")
        print(f"  ApplyBinAllStages = {ApplyBinAllStages}")
        print(f"  MinVarForBspline = {MinVarForBspline}")
    
    NumOfTx = len(ListOfSitkTxs)
    
    for i in range(NumOfTx):
        SitkTx = ListOfSitkTxs[i]
        
        SitkTxName = SitkTx.GetName()
        
        if ApplySameInterpAllStages:
            interp = deepcopy(Interp)
        else:
            if i < NumOfTx - 1:
                interp = 'Linear'
            else:
                interp = deepcopy(Interp)
        
        if ApplyPostTxBlur:
            if ApplyBlurAllStages:
                applyBlur = True
            else:
                if i < NumOfTx - 1:
                    applyBlur = False
                else:
                    applyBlur = True
        else:
            applyBlur = False
        
        if ApplyPostTxBin:
            if ApplyBinAllStages:
                applyBin = True
            else:
                if i < NumOfTx - 1:
                    applyBin = False
                else:
                    applyBin = True
        else:
            applyBin = False
        
        if SitkTxName == 'BSplineTransform':
            var = list(PostTxVariance)
            
            if PostTxVariance[0] < MinVarForBspline:
                var[0] = MinVarForBspline
            if PostTxVariance[1] < MinVarForBspline:
                var[1] = MinVarForBspline
            if PostTxVariance[2] < MinVarForBspline:
                var[2] = MinVarForBspline
            
            var = tuple(var)
        else:
            var = deepcopy(PostTxVariance)
        
        
        if LogToConsole:
            print(f"  Sitk Transform {i+1} of {NumOfTx}:")
            print(f"    SitkTxName = {SitkTxName}")
            print(f"    Interpolation to be applied = {interp}")
            print(f"    Gaussian blurring to be applied = {applyBlur}")
            if applyBlur:
                print(f"    Variance to be applied = {var}")
            print(f"    Binary thresholding to be applied = {applyBin}")
            
        LabImByRoi, PixArrByRoi, F2SindsByRoi\
        = TransformLabImByRoi(LabImByRoi=LabImByRoi, 
                              F2SindsByRoi=F2SindsByRoi, 
                              TrgIm=TrgIm, 
                              SelxImFiltOrSitkTx=SitkTx,
                              Interp=interp, 
                              ApplyPostTxBlur=applyBlur, 
                              PostTxVariance=var,
                              ApplyPostTxBin=applyBin, 
                              LogToConsole=LogToConsole)
    
    if LogToConsole:
        print('-'*120)
        
    return LabImByRoi, PixArrByRoi, F2SindsByRoi