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



def ImportImage(DicomDir):
    """
    Import a DICOM series as a SimpleITK image.
    
    Inputs:
    ******
        
    DicomDir : string
        Directory containing the DICOMs.
        
        
    Outputs:
    *******
    
    Image : SimpleITK object
        SimpleITK 3D image of the DICOM stack.
    """
    
    import SimpleITK as sitk
    
    
    Reader = sitk.ImageSeriesReader()
    Fnames = Reader.GetGDCMSeriesFileNames(DicomDir)
    Reader.SetFileNames(Fnames)
    #Reader.ReadImageInformation() # doesn't exist for ImageSeriesReader; works
    # for ImageFileReader
    Image = Reader.Execute()
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







def ImportNifti(FilePath):
    """
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
            
            msg = f'\nWarning:  The voxel spacings along the scan direction '\
                  + f'are not to within 1% of each other:\n'\
                  + f'{UniqueVlengths}'
            
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
        print(f'      Image PixID = {PixID}')
        print(f'      Image PixIDTypeAsString = {PixIDTypeAsStr}')
        print(f'      Image Size = {Image.GetSize()}')
        print(f'      Image Max = {ImageMax(Image)}')
        print(f'      Image Min = {ImageMin(Image)}')
        print(f'\n      Conversion of Image to PixArr:')
        print(f'      PixArr shape = {PixArr.shape}')
        if isinstance(UniqueVals, np.ndarray):
            if len(UniqueVals) < 7:
                print(f'\n      There are {len(UniqueVals)} unique values',
                      f'in PixArr:')
                print(f'      {UniqueVals}')
            else:
                
                print(f'\n      There are {len(UniqueVals)} unique values',
                      f'in PixArr:')
                print(f'      {UniqueVals[:3]}...{UniqueVals[-3:-1]}')
            
            """ Print histogram of UniqueVals
            https://www.machinelearningplus.com/plots/matplotlib-histogram-python-examples/
            """
            print('\n      Distribution of values:')
            Vals = PixArr.flatten()
            MinVal = min(Vals)
            MaxVal = max(Vals)
            
            Freq, Bins = np.histogram(Vals, bins=10, range=[MinVal, MaxVal])
            for b, f in zip(Bins[1:], Freq):
                #print(round(b, 1), ' '.join(np.repeat('*', f)))
                print(f'      {round(b, 2)} - {f}')
            
            print('\n      Distribution of values near 0:')
            Vals = PixArr.flatten()
            MinVal = 0
            MaxVal = 0.2
            
            Freq, Bins = np.histogram(Vals, bins=10, range=[MinVal, MaxVal])
            for b, f in zip(Bins[1:], Freq):
                #print(round(b, 1), ' '.join(np.repeat('*', f)))
                print(f'      {round(b, 2)} - {f}')
            
            
        elif UniqueVals == None:
            print(f'         There are no UniqueVals (= {UniqueVals}) in',
                  f'PixArr:')
        print(f'         \nThere are {len(F2Sinds)} frames with slice indices:')
        print(f'         {F2Sinds}')
    
    return PixID, PixIDTypeAsStr, UniqueVals, F2Sinds







def GetImageVertices(Image):
    """
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







def InitialiseImage(RefImage, PixType='SameAsRefImage'):
    """
    Initialise an empty SimpleITK image with the same size, spacing, origin and
    direction as another image.  
    
    Inputs:
    ******
    
    RefImage : SimpleITK image
        The reference image whose attributes will be used for the initialised
        image.
    
    PixType : string (optional; 'SameAsRefImage' by default)
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
    
    if PixType == 'SameAsRefImage':
        PixId = RefImage.GetPixelID()
    elif PixType == 'uint8':
        PixId = sitk.sitkUInt8
    elif PixType == 'int8':
        PixId = sitk.sitkInt8
    elif PixType == 'uint16':
        PixId = sitk.sitkUInt16
    elif PixType == 'int16':
        PixId = sitk.sitkInt16
    elif PixType == 'uint32':
        PixId = sitk.sitkUInt32
    elif PixType == 'int32':
        PixId = sitk.sitkInt32
    elif PixType == 'uint64':
        PixId = sitk.sitkUInt64
    elif PixType == 'int64':
        PixId = sitk.sitkInt64
    elif PixType == 'float32':
        PixId = sitk.sitkFloat32
    elif PixType == 'float64':
        PixId = sitk.sitkFloat64
    else:
        msg = f"PixType '{PixType}' is not a valid input."
        raise Exception(msg)
    
    #NewImage = sitk.Image(RefImage.GetSize(), RefImage.GetPixelID())
    NewImage = sitk.Image(RefImage.GetSize(), PixId)
    
    NewImage.SetOrigin(RefImage.GetOrigin())
    
    NewImage.SetDirection(RefImage.GetDirection())
    
    NewImage.SetSpacing(RefImage.GetSpacing())
    
    return NewImage





def ImageMax(Image):
    """
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
        print(f'Results of FindSuitableThresh():')
    
    
    Ratio = GetVoxelVolumeRatio(BinaryIm, NonBinaryIm)
    
    
    BinPixArr, BinF2Sinds = Image2PixArr(BinaryIm)
    NonBinPixArr, NonBinF2Sinds = Image2PixArr(NonBinaryIm)
    
    NumOfOnes_B = np.count_nonzero(BinPixArr)
            
    Vals_NB = NonBinPixArr.flatten()
    MinVal = min(Vals_NB)
    MaxVal = max(Vals_NB)
    
    Freq, Bins = np.histogram(Vals_NB, bins=1000, range=[MinVal, MaxVal])
    
    if LogToConsole:
    #if False:
        for b, f in zip(Bins[1:], Freq):
            print(f'      {round(b, 8)} - {f}')
    
    
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
    """
    
    import SimpleITK as sitk
    
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
    Resample a 3D SimpleITK image.  The image can be a 3D DICOM image or a
    3D labelmap image.
    
    Inputs:
    ******                      
        
    Im : SimpleITK image
        The 3D image to be resampled.
                           
    RefIm : SimpleITK image
        The 3D image reference image whose gridspace Image will be resampled to.
        
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
    
    Resampler : SimpleITK image filter
    

    
    Notes:
    *****
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """

    import SimpleITK as sitk
    
    # Define which tranform to use:
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
        sitkPixType = sitk.sitkFloat32
        
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
    #    print('\nResampling tranform:\n', ResTx)

    return ResIm
    #return ResIm, ResTx








def ResampleLabIm(LabIm, F2Sinds, SrcIm, TrgIm, Interp, 
                  PreResVariance=(1,1,1), ApplyPostResBlur=True, 
                  PostResVariance=(1,1,1), LogToConsole=False):
    """
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
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ResampleLabIm():')
        
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
                  f'Try using the \'BlurThenLinear\' approach.\n')
            
            Interp = 'BlurThenLinear'

    
    
    if Interp == 'BlurThenLinear':
        if LogToConsole:
            print(f'\nResampling LabIm by applying Gaussian blur,',
                  'resampling using a linear interpolator, then binary',
                  'thresholding...')
        
        
        #""" The original binary pixel array: """
        #PixArr_B, F2Sinds = Image2PixArr(LabIm)
            
        
        """ Convert LabIm from 32-bit unsigned integer to float.
        Note:  Result of resampling results in empty labelmap unless
        image is converted to float prior to resampling. """
        LabIm = ConvertImagePixelType(Image=LabIm, 
                                      NewPixelType='sitkFloat32')
        
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
            print(f'\nImage info for LabIm after applying Gaussian blur:')
            
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
        ResLabIm = ConvertImagePixelType(Image=ResLabIm, 
                                         NewPixelType='sitkUInt32')
        
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
        print(f'Results of ResampleLabImByRoi():')
        
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









def RegisterImages_OLD(FixIm, MovIm, Tx='bspline', MaxNumOfIters=512,
                   LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Tx : string (optional; 'bspline' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    
    MaxNumOfIters : integer
        The maximum number of iterations to be used during optimisation.
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    RegImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
    """
    
    import SimpleITK as sitk
    import time
    
    if not Tx in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transformation (Tx), {Tx}, is not one of the '\
              + 'accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        
        raise Exception(msg)
        
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Initiate RegImFilt: """
    RegImFilt = sitk.ElastixImageFilter()
    RegImFilt.LogToConsoleOn() # <-- no output in Jupyter
    #RegImFilt.LogToConsoleOff() 
    RegImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    """ Define the fixed and moving images: """
    RegImFilt.SetFixedImage(FixIm)
    RegImFilt.SetMovingImage(MovIm)
    
    """ Get the default parameter map template for the chosen transformation: """
    #RegImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    #ParamMap = sitk.GetDefaultParameterMap(Tx)
    
    if Tx == 'bspline':
        """ First register using an affine transformation: """
        ParamMap = sitk.GetDefaultParameterMap('affine')
        """ 22/04/2021: I'm not sure this was the correct approach to 
        performing a BSpline registration... """
    else:
        ParamMap = sitk.GetDefaultParameterMap(Tx)
        
    
    """ Re-assign some parameters: """
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
    ParamMap['WriteIterationInfo'] = ['true']
    #ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['MaximumNumberOfIterations'] = [str(MaxNumOfIters)]
    ParamMap['UseDirectionCosines'] = ['true']
    
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    """ 13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70 """
    #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    
    if Tx == 'bspline':
        """ 22/02/21 Try parameters in Appendix A in "Multimodal image 
        registration by edge attraction and regularization using a B-spline 
        grid" by Klein et al 2011; and from the SimpleElastix docs:
        https://simpleelastix.readthedocs.io/ParameterMaps.html """
        
        #ParamMap['MaximumNumberOfIterations'] = ['250']
        #ParamMap['MaximumNumberOfIterations'] = ['512']
        ParamMap['MaximumNumberOfIterations'] = [str(MaxNumOfIters)]
        #ParamMap['FinalGridSpacingInPhysicalUnits'] = ['20']
        #ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
        ParamMap['Registration'] = ['MultiResolutionRegistration']
        ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent']
        ParamMap['Transform'] = ['BSplineTransform']
        ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
        ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
        ParamMap['NumberOfResolutions'] = ['4']
        ParamMap['Metric'] = ['AdvancedMattesMutualInformation']
        ParamMap['NumberOfHistogramBins'] = ['32']
        ParamMap['NumberOfSpatialSamples'] = ['2048']
        ParamMap['NewSamplesEveryIteration'] = ['true']
        ParamMap['ImageSampler'] = ['RandomCoordinate']
        ParamMap['Interpolator'] = ['BSplineInterpolator']
        ParamMap['BSplineInterpolatorOrder'] = ['1']
        #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # typo? (found 22/04/2021)
        ParamMap['FinalBSplineInterpolationOrder'] = ['3'] # 22/04/2021
        ParamMap['Resampler'] = ['DefaultResampler']
        ParamMap['DefaultPixelValue'] = ['0']
        ParamMap['FixedInternalImagePixelType'] = ['float']
        ParamMap['UseDirectionCosines'] = ['true']
        ParamMap['MovingInternalImagePixelType'] = ['float']
        ParamMap['HowToCombineTransforms'] = ['Compose']
    
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    """ Set the parameter map: """
    RegImFilt.SetParameterMap(ParamMap)
    
    if LogToConsole:
        print(f'Performing image registration using {Tx} transform...\n')
        
    #RegImFilt.AddCommand(sitk.sitkProgressEvent, 
    #                     lambda: print("\rProgress: {0:03.1f}%...".format(100*RegImFilt.GetProgress()),end=''))
    # 'ElastixImageFilter' object has no attribute 'AddCommand'
    
    """ Register the 3D images: """
    RegImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = RegImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
        #print('\n', RegImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', RegImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
    
    
    return RegIm, RegImFilt
    #return RegIm, RegImFilt, ParamMap
    
    
    
    
    




def RegisterImagesSelx(FixIm, MovIm, Transform='affine', 
                       MaxNumOfIters='default', InitMethod='centerofgravity',
                       FixFidsFpath='fixed_fiducials.txt', 
                       MovFidsFpath='moving_fiducials.txt',
                       FlipK=False,
                       LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
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
    
    MaxNumOfIters : string (optional; 'default' by default)
        If 'default', the maximum number of iterations used during optimisation
        will be whatever is set by default in the parameter map of choice.  
        If != 'default' it must be a string representation of an integer.
    
    InitMethod : string (optional; 'centerofgravity' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'centerofgravity'
            - 'geometricalcenter'
            - 'origin' (only available for affine transform)
            - 'landmarks'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    SelxImFilt : Elastix image filter
        Elastix image registration transformation filter used to register/
        transform MovIm to FixIm.
    
    TxParamMapDict : dictionary
        A dictionary containing the transform parameter map (including the 
        transform parameters).
    
        
    Notes:
    *****
    
    29/05/20: Trying this instead of trying to change it for Transformix
    ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0']
    
    13/02/21 Try https://github.com/SuperElastix/SimpleElastix/issues/70
    ParamMap['Registration'] = ['MultiMetricResolutionRegistration']
    
    22/02/21 Try parameters in Appendix A in "Multimodal image registration by 
    edge attraction and regularization using a B-spline grid" by Klein et al 
    2011; and from the SimpleElastix docs:
    https://simpleelastix.readthedocs.io/ParameterMaps.html
    """
    
    import SimpleITK as sitk
    import time
    
    if not Transform in ['rigid', 'affine', 'bspline']:
        msg = f'The chosen transform (Transform), {Transform}, is not one of '\
              + 'the accepted inputs: \'rigid\', \'affine\', or \'bspline\'.'
        raise Exception(msg)
    
    if not InitMethod in ['centerofgravity', 'geometricalcenter', 'origins',
                          'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              "is not one of the accepted inputs: 'centerofgravity', "\
              + "'geometricalcenter', 'origins' or 'landmarks'."
        raise Exception(msg)
        
    if Transform != 'affine' and InitMethod == 'origins':
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              "is only allowed for the affine transform."
        raise Exception(msg)
    
    
    """ Initiate RegImFilt: """
    SelxImFilt = sitk.ElastixImageFilter()
    if LogToConsole:
        SelxImFilt.LogToConsoleOn() # no output in Jupyter
    else:
        SelxImFilt.LogToConsoleOff() 
    SelxImFilt.LogToFileOn() # output's elastix.log ONLY ONCE per kernel in Jupyter
    
    SelxImFilt.SetFixedImage(FixIm)
    SelxImFilt.SetMovingImage(MovIm)
    
    """ Get the default parameter map for the chosen transform: """
    #ElxImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    ParamMap = sitk.GetDefaultParameterMap(Transform)
        
    """ 18/03/21:
    
    Kasper said that AutomaticTransformInitialization should be true for rigid,
    and false for affine or bspine: 
    
    https://github.com/SuperElastix/elastix/issues/45
         
    - FixedImagePyramid should be FixedImageSmoothingPyramid for best results 
    (unless memory is really and issue).
    - AutomaticTransformInitialization should be true for rigid and false for 
    affine and bspline.
    - 250 number of iterations is also on the low side. If speed is really an 
    issue, 250 is fine, but if you can afford it use at least 512.
    - Depending on the problem you might want to use at least 3 resolutions. 
    Some places you use 2 which is on the low side.
    - ImageSampler should be RandomSparseMask if users use masks. Otherwise it 
    is fine.
    - If you use a BSplineInterplationOrder of 1 you might as well use 
    LinearInterpolator which produces the same results but is faster.

    """
    
    """ Add/modify parameters: """
    if Transform == 'rigid': # 26/04/21
        ParamMap['AutomaticTransformInitialization'] = ['true']
    else:
        ParamMap['AutomaticTransformInitialization'] = ['false']
    #ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter'] # 19/03/21
    #ParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity'] # 19/03/21
    #ParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
    
    if InitMethod == 'centerofgravity':
        ParamMap['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
    elif InitMethod == 'geometricalcenter':
        ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    elif InitMethod == 'origins':
        ParamMap['AutomaticTransformInitializationMethod'] = ['Origins']
        
    
    ParamMap['AutomaticScalesEstimation'] = ['true'] # 22/02/21
    #ParamMap['DefaultPixelValue'] = ['0']
    """ 05/05/21: Not sure why BSplineInterpolatorOrder was set to 1. """
    if False:
        ParamMap['BSplineInterpolatorOrder'] = ['1']
    #ParamMap['ResampleInterpolator'] = ['FinalBSplineInterpolator'] # recommended in elastix manual (5.3.4 Interpolator)
    #ParamMap['FinalBSplineInterpolatorOrder'] = ['3'] # default for affine and bspline
    #ParamMap['FixedImagePyramid'] = ['FixedSmoothingImagePyramid'] # default for affine and bspline
    #ParamMap['MovingImagePyramid'] = ['MovingSmoothingImagePyramid'] # default for affine and bspline
    #if Transform == 'bspline': # 22/04/2021
    #    ParamMap['FixedImagePyramid'] = ['FixedRecursiveImagePyramid']
    #    ParamMap['MovingImagePyramid'] = ['MovingRecursiveImagePyramid']
    ParamMap['FixedInternalImagePixelType'] = ['float']
    ParamMap['MovingInternalImagePixelType'] = ['float']
    ParamMap['HowToCombineTransforms'] = ['Compose']
    #ParamMap['ImageSampler'] = ['RandomCoordinate'] # default for affine and bspline
    #ParamMap['Interpolator'] = ['LinearInterpolator'] # default for affine
    #ParamMap['Interpolator'] = ['BSplineInterpolator'] # default for bspline
    #ParamMap['MaximumNumberOfIterations'] = ['256'] # default for affine
    #ParamMap['MaximumNumberOfIterations'] = ['512'] # default for bspline
    if MaxNumOfIters != 'default':
        ParamMap['MaximumNumberOfIterations'] = [MaxNumOfIters]
    #ParamMap['Metric'] = ['AdvancedMattesMutualInformation'] # default for affine
    #ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
    #                      'TransformBendingEnergyPenalty'] # default for bspline
    if InitMethod == 'landmarks':
        if Transform in ['rigid', 'affine']:
            ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                  'CorrespondingPointsEuclideanDistanceMetric']
            #ParamMap['Metric0Weight'] = ['1.0']
            ParamMap['Metric0Weight'] = ['0.3']
            ParamMap['Metric1Weight'] = ['1.0']
            
        elif Transform == 'bspline':
            ParamMap['Metric'] = ['AdvancedMattesMutualInformation', 
                                  'TransformBendingEnergyPenalty',
                                  'CorrespondingPointsEuclideanDistanceMetric']
            ParamMap['Metric0Weight'] = ['1.0']
            ParamMap['Metric1Weight'] = ['1.0']
            ParamMap['Metric2Weight'] = ['1.0']
    #ParamMap['NewSamplesEveryIteration'] = ['true'] # default for bspline; recommended in elastix manual (5.3.3 Sampler)
    ParamMap['NumberOfHistogramBins'] = ['32']
    #ParamMap['NumberOfResolutions'] = ['4'] # default for affine and bspline
    #ParamMap['NumberOfSpatialSamples'] = ['2048'] # default for affine and bspline
    #if Transform == 'bspline': # 22/04/2021
    #    ParamMap['NumberOfSpatialSamples'] = ['3000'] # recommended in elastix manual (5.3.3 Sampler)
    #    ParamMap['Metric0Weight'] = ['0.5']
    #    ParamMap['Metric1Weight'] = ['1.0']
    #ParamMap['Optimizer'] = ['AdaptiveStochasticGradientDescent'] # default for affine and bspline
    if InitMethod == 'landmarks':
        """ MultiResolutionRegistration, which is default for affine, only
        allows for 1 metric. """
        ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration']
    #ParamMap['Registration'] = ['MultiResolutionRegistration'] # default for affine
    #ParamMap['Registration'] = ['MultiMetricMultiResolutionRegistration'] # default for bspline
    #ParamMap['Resampler'] = ['DefaultResampler'] # default for affine and bspline
    #ParamMap['Transform'] = ['AffineTransform'] # default for affine
    #ParamMap['Transform'] = ['BSplineTransform'] # default for bspline
    ParamMap['UseDirectionCosines'] = ['true']
    ParamMap['WriteIterationInfo'] = ['true']
    
    
    # Print the parameters:
    #for keys,values in ParamMap.items():
    #    print(keys, '=', values)
        
    """ Set the parameter map: """
    SelxImFilt.SetParameterMap(ParamMap)
    
    if InitMethod == 'landmarks':
        SelxImFilt.SetFixedPointSetFileName(FixFidsFpath)
        SelxImFilt.SetMovingPointSetFileName(MovFidsFpath)
    
    if LogToConsole:
        print(f'Performing image registration using {Transform} transform...\n')
        
    
    """ Start timing: """
    times = []
    times.append(time.time())
    
    """ Register the 3D images: """
    SelxImFilt.Execute()
    
    """ Get the registered image: """
    RegIm = SelxImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
        #print('\n', ElxImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', ElxImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
        
    
    """ Get the items in parameter map and in the transform parameter map
    (both from file and from the Elastix image filter) as a dictionary. """ 

    """ First get the parameter map from SelxImFilt: """
    ParamMapDict = ParamMapFromSelxImFilt(SelxImFilt)
    
    """ Add the transform parameters from SelxImFilt: """
    TxParamMapDict = TxParamMapFromSelxImFilt(SelxImFilt, ParamMapDict)
    
    """ Add the additional transform parameters (and over-write existing ones*)
    from the file TransformParameters.0.txt.
    
    * There are more decimals in the parameters containted from the file than
    those stored in ElxImFilt. """
    TxParamMapDict = TxParamMapFromFile(TxParamMapDict)
    
    
    #return RegIm, SelxImFilt
    return RegIm, SelxImFilt, TxParamMapDict







def ParamMapFromSelxImFilt(SelxImFilt, Dict=None):
    """
    Get the Elastix parameter map from a Elastix image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from 
        .GetTransformParameterMap()) is provided the parameters (from 
        .GetParameterMap()) will be added to it.  Any keys in Dict (e.g. 
        TransformParameterMap) that are also present in ParameterMap will be 
        overwritten. 
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix parameters.
    """
    
    if not Dict:
        Dict = {}
    
    for key, value in SelxImFilt.GetParameterMap()[0].items():
        Dict.update({key : value})
    
    return Dict







def TxParamMapFromSelxImFilt(SelxImFilt, Dict=None):
    """
    Get the Elastix transform parameter map from a Elastix image filter.
    
    Inputs:
    ******
    
    SelxImFilt : Elastix image filter
        Elastix image transformation filter used to perform image registration.
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from .GetParameterMap()) is 
        provided the transform parameters (from .GetTransformParameterMap())
        will be added to it.  Any keys in Dict (e.g. ParameterMap) that are 
        also present in TransformParameterMap will be overwritten. 
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix parameters.
    
    
    Notes:
    *****
    
    1. The transform parameter map can be obtained as follows:
    
    SelxImFilt.GetTransformParameterMap()[0]
    
    and printed as follows:
        
    sitk.PrintParameterMap(SelxImFilt.GetTransformParameterMap()[0])
    
    however printing to the console does not work in Jupyter.  Although it is
    possible to print the value of a specific key to the console as follows:
        
    SelxImFilt.GetTransformParameterMap()[0]['Direction']
    
    or of the entire dictionary:
        
    for key, value in SelxImFilt.GetTransformParameterMap()[0].items():
        print(f'{key} = {value}')
        
    Hence this function is simply a way of getting around the Jupyter notebook
    limitation.
    
    2. Oddly the transform parameters read from file (TransformParameters.0.txt)
    have more significant digits than that read directly from the Elastix image 
    filter.
    
    e.g. 
    
    list(SelxImFilt.GetTransformParameterMap()[0]['TransformParameters'])
    
    returns:
        
    ['1.0046',
     '-0.00831448',
     '-0.0412239',
     '0.0140337',
     '0.997238',
     '0.0593744',
     '0.0420835',
     '-0.0564159',
     '1.00477',
     '-1.70085',
     '-11.4156',
     '-54.8056']
    
    whereas

    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams]
    
    TxParams
    
    returns:
    
    ['1.004604',
     '-0.008314',
     '-0.041224',
     '0.014034',
     '0.997238',
     '0.059374',
     '0.042084',
     '-0.056416',
     '1.004766',
     '-1.700848',
     '-11.415621',
     '-54.805552']

    """
    
    if not Dict:
        Dict = {}
    
    for key, value in SelxImFilt.GetTransformParameterMap()[0].items():
        Dict.update({key : value})
    
    return Dict






def TxParamMapFromFile(Dict=None):
    """
    Get the Elastix transform parameter map from TransformParameters.0.txt in  
    the current working directory.
    
    Inputs:
    ******
    
    Dict : dictionary (optional; None by default)
        If a dictionary (e.g. containing parameters from .GetParameterMap()) is 
        provided the transform parameters (from TransformParameters.0.txt)
        will be added to it.  Any keys in Dict (e.g. ParameterMap) that are 
        also present in TransformParameters.0.txt will be overwritten.
    
    
    Outputs:
    *******
    
    Dict : dictionary
        A dictionary containing Elastix transform parameters (and other 
        parameters if an input dictionary was provided).
        If the file TransformParameters.0.txt is not found in the current 
        working directory an empty (or the input) dictionary will be returned.
    
    
    Notes:
    *****
    
    The transform parameters read from file (TransformParameters.0.txt)
    have more significant digits than that read directly from the Elastix image 
    filter.
    
    e.g. 
    
    list(SelxImFilt.GetTransformParameterMap()[0]['TransformParameters'])
    
    returns:
        
    ['1.0046',
     '-0.00831448',
     '-0.0412239',
     '0.0140337',
     '0.997238',
     '0.0593744',
     '0.0420835',
     '-0.0564159',
     '1.00477',
     '-1.70085',
     '-11.4156',
     '-54.8056']
    
    whereas

    ParamFname = 'TransformParameters.0.txt'
    ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
    TxParams = ParseTransformParameterFile(ParamFpath, 'TransformParameters')

    TxParams = [str(item) for item in TxParams]
    
    TxParams
    
    returns:
    
    ['1.004604',
     '-0.008314',
     '-0.041224',
     '0.014034',
     '0.997238',
     '0.059374',
     '0.042084',
     '-0.056416',
     '1.004766',
     '-1.700848',
     '-11.415621',
     '-54.805552']

    """
    
    import os
    
    if Dict == None:
        Dict = {}
    
    Fname = 'TransformParameters.0.txt'
    
    if os.path.isfile(Fname):
        # Read in the file:
        with open(Fname, 'r') as FileText:
            lines = FileText.readlines()
            
            
            for line in lines:
                if line[0] == '(':
                    # Get the index of the first white space and ')':
                    a = line.index(' ')
                    b = line.index(')')
                    
                    key = line[1:a]
                    value = line[a+1:b].replace('"', '')
                    
                    # Get values as a list of space-delimited strings:
                    values = value.split()
                    
                    Dict.update({key : tuple(values)})

    
    return Dict






def ParseTransformParameterFile(Fpath, ParamName):
    """
    Parse an Elastix transform parameter map from TransformParameters.0.txt for 
    the value(s) of a specific parameter.
            
    Inputs:
    ******
    
    Fpath : string
        The full filepath of TransformParameters text file.
        
    ParamName : string
        Character string of the parameter name to search and parse.
        
        
    Outputs:
    *******
    
    ParamValues : strings, integers or floats (depends on the parameter of 
    interest)
        The value(s) of the parameter of interest.
    
    
    Notes:
    -----
    
    1.  The text file contains data of the following format:
            
        (FixedImageDimension 2)     <-- single value
        (Size 256 256)              <-- two values
        
    2. See further notes below.
    """

    # Determine whether the parameters need to be converted to integers or
    # float.  
    
    # Create a list of parameter names that belong to int and float data types:
    """
    Notes: 
        1. Many of the parameters are to remain as strings.
        2. The lists below only cover the default parameters, so any additional
           parameters will need to be added.
        3. If the parameter is not found to be in ints or floats it will be 
           assumed that the parameter is a string.
    """
    
    IntParams = ['NumberOfParameters', 'FixedImageDimension', 
                 'MovingImageDimension', 'Size', 'Index', 
                 'FinalBSplineInterpolationOrder']
    
    FloatParams = ['TransformParameters', 'Spacing', 'Origin', 'Direction', 
                   'CenterOfRotationPoint', 'DefaultPixelValue']
    
    if ParamName in IntParams:
        dtype = 'int'
    elif ParamName in FloatParams:
        dtype = 'float'
    else:
        dtype = 'string'
        
    
    # Read in the file:
    FileTxt = open(Fpath, 'r').read()
    
    #print('\nFileTxt:\n', FileTxt)
    
    # Get index of string variable parameter:
    a = FileTxt.index(ParamName)
    
    #print(f'\na = {a}')
    
    """ Note:
    The text file has values encoded in the following format:
    (ParameterName ParameterValue1 ParameterValue2 ...)
    
    so the values begin at a+len(ParameterName)+1 and end at the position 
    before ')'.
    """
    StartInd = a + len(ParamName) + 1 # to account for the space after 
    # ParameterName
    
    #print(f'\nStartInd = {StartInd}')
    
    # Find the position of the first ')' after StartInd:
    b = StartInd + FileTxt[StartInd:].index(')')
    
    #print(f'\nb = {b}')
    
    EndInd = b
    
    #print(f'\nEndInd = {EndInd}')
    
    # Crop the text between StartInd and EndInd:
    CroppedTxt = FileTxt[StartInd:EndInd]
    
    #print(f'\nCroppedTxt = {CroppedTxt}')
    
    """ Note:
    CroppedTxt may be single value or may have multiple values. 
    e.g. (FixedImageDimension 2) <-- single value
         (Size 256 256) <-- two values
         
    If CroppedTxt doesn't contain the space character simply return CroppedTxt.
    If CroppedTxt has a space character, slice and append the elements that 
    preceed it to a new array called ParameterValues, then search for another 
    space, repeat until there are no more spaces.
    """
    # Check if there are any space characters in CroppedTxt:
    AreThereSpaces = ' ' in CroppedTxt
    
    #print('\nAreThereSpaces =', AreThereSpaces)
    
    if AreThereSpaces:
        ParamValues = []
        
        # Continue as long as AreThereSpaces is True:
        while AreThereSpaces:
            # Find the first space character:
            ind = CroppedTxt.index(' ')

            # Slice the characters from 0 to ind:
            ParamValues.append(CroppedTxt[0:ind])
            
            # Remove from CroppedTxt the characters that were appended to 
            # ParamValues but not including the space character:
            CroppedTxt = CroppedTxt[ind+1:]
            
            #print(f'\n   CroppedTxt = {CroppedTxt}')
            
            # Check if there are any space characters in CroppedTxt:
            AreThereSpaces = ' ' in CroppedTxt
            
            #print('\n   AreThereSpaces =', AreThereSpaces)
            
        # There are no remaining spaces so simply append what remains of
        # CroppedTxt to ParamValues:
        if CroppedTxt:
            ParamValues.append(CroppedTxt)
            
        #return ParamValues
        if dtype=='int':
            return [int(item) for item in ParamValues]
        elif dtype=='float':
            return [float(item) for item in ParamValues]
        else:
            return ParamValues
        
    else:
        #return CroppedTxt
        if dtype=='int':
            return [int(item) for item in CroppedTxt]
        elif dtype=='float':
            return [float(item) for item in CroppedTxt]
        else:
            return CroppedTxt
        






def GetTxMatrixType(TxParams, LogToConsole=False):
    """
    Determine the Transformation Matrix Type from a list of transformation 
    parameters read from the TransformParameters.0.txt file exported by
    Elastix following an image registration and assumed to be in the current
    working directory.
    
    Inputs:
    ******
    
    TxParams : list of strings
        The registration transform parameters. If TxParams = None, the 
        transform parameters will be parsed from TransformParameters.0.txt
        (see Notes).
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
    
        
    Outputs:
    *******
    
    MatrixType : string
        A string containing the matrix type.  Possible outputs are:
            - 'RIGID'
            - 'RIGID_SCALE'
            - 'AFFINE'
        
        
    Notes:
    -----
    
    1. The transform parameters will be read from the file 
    TransformParameters.0.txt.  The file is generated by Elastix following an
    image registration, found in the current working directory. It will be 
    assumed that the file is located in the current working directory and that
    the filename is as above.
    
    2. The matrix M is the 3x3 matrix within the 4x4 matrix TxParams read from
    TransformParameters.0.txt.  M only includes the elements that describe 
    rotations, scaling and shearing, i.e not translations (hence only the
    elements [0, 1, 2, 4, 5, 6, 8, 9, 10] in TxParams).
    
    3. There are three types of Registration matrices as defined in C.X.1.1.2
    (page 27) in DICOM Standard Supplement 7.3:
        
        RIGID: 
        - a registration involving only translations and rotations
        - the matrix contains 6 degrees of freedom: 3 translations and 3 
        rotations
        - the matrix must be orthnormal
        
        RIGID_SCALE:
        - a registration involving only translations, rotations and scaling
        - the matrix contains 9 degrees of freedom: 3 translations, 3 rotations
        and 3 scales
        - the matrix must be orthogonal
        
        AFFINE:
        - a registration involving translations, rotations, scaling and 
        shearing
        - the matrix contains 12 degrees of freedom: 3 translations, 3 
        rotations, 3 scales and 3 shears
        - there are no constraints on the matrix elements
    """
    
    import os
    import numpy as np
    #from ImageTools import ParseTransformParameterFile
    from GeneralTools import IsMatrixOrthonormal, IsMatrixOrthogonal
    
    #print(f'TxParams = {TxParams}\n')
    
    if TxParams == None:
        ParamFname = 'TransformParameters.0.txt'
        ParamFpath = os.path.join(os.getcwd(), ParamFname)
    
        TxParams = ParseTransformParameterFile(ParamFpath, 
                                               'TransformParameters')
    
    #inds = [0, 1, 2, 4, 5, 6, 8, 9, 10] # 06/05/21
    inds = [0, 1, 2, 3, 4, 5, 6, 7, 8] # 06/05/21

    #M = np.array([TxParams[i] for i in inds]) # 06/05/21
    M = np.array([float(TxParams[i]) for i in inds]) # 06/05/21
    
    M = np.reshape(M, (3, 3))
    
    #print(f'M = {M}\n')
    
    IsOrthonormal = IsMatrixOrthonormal(M, LogToConsole)
    
    if IsOrthonormal:
        MatrixType = 'RIGID'
    else:
        IsOrthogonal = IsMatrixOrthogonal(M, LogToConsole)
        
        if IsOrthogonal:
            MatrixType = 'RIGID_SCALE'
        
        else:
            MatrixType = 'AFFINE'
            
            
    if LogToConsole:
        print('Matrix type is', MatrixType)
    
    return MatrixType






def CompareSelxParameterMaps(SelxImFilt0, SelxImFilt1, DiffsOnly=True):
    ParamMap0 = SelxImFilt0.GetParameterMap()[0]
    ParamMap1 = SelxImFilt1.GetParameterMap()[0]
    
    Keys0 = ParamMap0.keys()
    Keys1 = ParamMap1.keys()
    
    print('SelxImFilt0 parameter v SelxImFilt1 parameter')
    print('*********************************************')
    
    for key0 in Keys0:
        if key0 in Keys1:
            if DiffsOnly:
                if ParamMap0[key0] != ParamMap1[key0]:
                    print(f"\n{key0}:")
                    print(f"{ParamMap0[key0]} v {ParamMap1[key0]}")
                else:
                    print(f"\n{key0}:")
                    print(f"{ParamMap0[key0]} v (key not present)")
            else:
                print(f"\n{key0}:")
                print(f"{ParamMap0[key0]} v {ParamMap1[key0]}")
    
    for key1 in Keys1:
        if not key1 in Keys0:
            if DiffsOnly:
                if ParamMap0[key1] != ParamMap1[key1]:
                    print(f"\n{key1}:")
                    print(f"{ParamMap0[key1]} v {ParamMap1[key1]}")
                else:
                    print(f"\n{key1}:")
                    print(f"(key not present) v {ParamMap1[key1]}")
            else:
                print(f"\n{key1}:")
                print(f"{ParamMap0[key1]} v {ParamMap1[key1]}")

    return




def TransformImageSelx(MovIm, SelxImFilt, Interp='NearestNeighbor', 
                       LogToConsole=False):
    """
    Transform a 3D SimpleITK image using SimpleElastix's Transformix.
    
    Inputs:
    ******
    
    MovIm : SimpleITK image 
        The 3D image to be transformed.
        
    SelxImFilt : Elastix image filter
        SimpleElastix image transformation filter used to perform image 
        registration.
        
    Interp : string (optional; 'NearestNeighbor' by default)
        The interpolator to be used for resampling when transforming the input
        labelmap image(s).  Accepatable values are:
        - 'Default' which leaves unchanged whatever interpolator was used in
        the image registration (i.e. RegImFilt)
        - 'NearestNeighbor'
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
        
    Outputs:
    *******
    
    TxIm : SimpleITK image
        The 3D transformed image.
        
        
    Notes:
    *****
    
    Simple use of Transformix as below:
        
    TxIm = sitk.Transformix(MovIm, ElxImFilt.GetTransformParameterMap())
        
    
    
    https://github.com/SuperElastix/SimpleElastix/issues/409
    """
    
    import SimpleITK as sitk
            
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformImageElx():')
        
    # Initiate TransformixImageFilter:
    TxImFilt = sitk.TransformixImageFilter()
    #TxImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TxImFilt.LogToConsoleOff() # <-- no output in Jupyter
    TxImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    TxImFilt.SetMovingImage(MovIm)
    
    # Get the parameter map used for the registration:
    TxParamMap = SelxImFilt.GetParameterMap()
    
    # Set the parameter map:
    #TxImFilt.SetTransformParameterMap(SelxImFilt.GetTransformParameterMap())
    
    """ Get the transform parameter map that was used for registration: """
    TxMap = SelxImFilt.GetTransformParameterMap()
    
    if LogToConsole:
        print('TxParamMap[0]["Interpolator"] =',
              f'{TxParamMap[0]["Interpolator"]}\n')
        print(f'TxMap[0]["ResampleInterpolator"] =',
              f'{TxMap[0]["ResampleInterpolator"]}\n')
        print(f'TxMap[0]["FinalBSplineInterpolationOrder"] =',
              f'{TxMap[0]["FinalBSplineInterpolationOrder"]}\n')
    
    
    if Interp == 'NearestNeighbor':
        interp = TxMap[0]["ResampleInterpolator"]
        
        TxMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
        
        if LogToConsole:
            print(f"The resample interpolator has been changed from {interp}",
                  "to 'FinalNearestNeighborInterpolator'\n")
    
    
    """ Get the image data type: (new 16/03/21)"""
    dtype = MovIm.GetPixelIDTypeAsString()
    
    if 'integer' in dtype:
        a = TxMap[0]["FinalBSplineInterpolationOrder"]
        
        """ The image is binary so set the FinalBSplineInterpolationOrder to 0
        to avoid "overshoot"-related "garbage" as described in the elastix
        manual section (4.3 The transform parameter file): """
        TxMap[0]["FinalBSplineInterpolationOrder"] = ["0"]
        
        if LogToConsole:
            print("The final BSpline interpolation order has been changed",
                  f"from {a} to '0'\n")
        
    TxImFilt.SetTransformParameterMap(TxMap)
    
    # Set up for computing deformation field:
    #TxImFilt.ComputeDeformationFieldOn()
    
    """ Transform MovIm: """
    TxImFilt.Execute()
    
    TxIm = TxImFilt.GetResultImage()
    
    if LogToConsole:
        print('-'*120)
    
    return TxIm









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








def ImportFiducialsFromTxt(Filepath, FlipK=False, Im=None):
    """
    Import a list of indices or points from a properly formatted text file.
    
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
    
    Im : SimpleITK image (optional unless FlipK is True)
        The 3D SimpleITK image that corresponds to the fiducials.
    
    
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
    
    if FlipK == True and Im == None:
        msg = "'FlipK' is True so a SimpleITK image must be provided."
        
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
                    S = Im.GetSize()[2]
                    
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
    
    return Fiducials, FiducialType






def GetLandmarkPtsFromImageJInds(FiducialFpath, Image):
    """ 
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







def AlignImagesUsingLandmarks(FixIm, MovIm, 
                              FixFidsFpath='fixed_fiducials.txt', 
                              MovFidsFpath='moving_fiducials.txt',
                              FlipK=False):
    """
    Return the landmark-based aligned image and tranform based on fiducials.
    
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
    
    
    Outputs:
    *******
    
    AlignedIm : SimpleITK image
        The 3D aligned image.
    
    Tx : SimpleITK transform
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
        
    print(f'Fixed fiducials: \n {FixFids}\n')
    print(f'Moving fiducials: \n {MovFids}\n')
        
    """ Flatten the list of points """
    FixPts = [c for p in FixFids for c in p]
    MovPts = [c for p in MovFids for c in p]
    
    Tx = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), 
                                                FixPts, MovPts)
    
    print(f'Landmark-based initial transformation is: \n {Tx.GetParameters()}')
    
    AlignedIm = sitk.Resample(MovIm, FixIm, Tx, sitk.sitkLinear, 0.0, 
                              MovIm.GetPixelID())
    
    return AlignedIm, Tx






import SimpleITK as sitk
import matplotlib.pyplot as plt
#from ipywidgets import interact, fixed
from IPython.display import clear_output


def PlotFixMovImages(FixInd, MovInd, FixPixArr, MovPixArr):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
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
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()


def PlotValues(RegMethod):
    """ Callback invoked when the IterationEvent happens, update our data 
    and display new figure. """
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







def RegisterImagesSitk(FixIm, MovIm, Transform='affine', 
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
            - 'geometry'
            - 'moments'
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
    
    if not InitMethod in ['geometry', 'moments', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              "is not one of the accepted inputs: 'geometry', 'moments', or "\
              + "'landmarks'."
        raise Exception(msg)
        
    if not FinalInterp in ['Linear', 'NearestNeighbor']:
        msg = f"The chosen final interpolator (FinalInterp), {FinalInterp}, is"\
              + " not one of the accepted inputs: 'Linear' or 'NearestNeighbor'."
        raise Exception(msg)
    
    
    if InitMethod == 'geometry':
        CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
    elif InitMethod == 'moments':
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
    
    
        
    """ Connect all of the observers so that we can perform plotting during 
    registration. """
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
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print('Optimizer\'s stopping condition,',
          f'{RegMethod.GetOptimizerStopConditionDescription()}')
    
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
        

    
    return RegIm, FinalTx






def GetLandmarkTx(FixIm, MovIm, FixFidsFpath, MovFidsFpath):
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





def GetLandmarkTx_v2(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
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






def ResampleImage_v2(Im, RefIm, ResTx, ResInterp=sitk.sitkLinear):
    import SimpleITK as sitk
    
    ResImFilt = sitk.ResampleImageFilter()
    #ResImFilt.SetTransform(sitk.Transform())
    ResImFilt.SetTransform(ResTx)
    ResImFilt.SetInterpolator(ResInterp)
    
    ResImFilt.SetReferenceImage(RefIm)
    #ResImFilt.SetOutputSpacing(RefIm.GetSpacing())
    #ResImFilt.SetSize(RefIm.GetSize())
    #ResImFilt.SetOutputDirection(RefIm.GetDirection())
    #ResImFilt.SetOutputOrigin(RefIm.GetOrigin())
    
    ResImFilt.SetDefaultPixelValue(RefIm.GetPixelIDValue())
    ResImFilt.SetOutputPixelType(RefIm.GetPixelID())
    
    
    ResIm = ResImFilt.Execute(Im)
    
    return ResIm





def RegUsingFiducials(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
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
        
    InitialTx = GetLandmarkTx(FixIm, MovIm, FixFidsFpath, MovFidsFpath)
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




def RegUsingFiducials_v2(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
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






def RegUsingFiducials_v3(FixIm, MovIm, FixFidsFpath, MovFidsFpath, 
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










def GetTransformFromParams(FixIm, TxMatrix, Transform):
    """
    Return a SimpleITK transform from a list of transform parameters.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image
        The 3D image whose grid is to be transformed to.
    
    TxMatrix : list of float strings
        List of float strings representing the transformation that maps the
        moving (Source) image to the fixed (Target) image (e.g. as parsed from
        the FrameOfReferenceTransformationMatrix of a DRO).
    
    Transform : SimpleITK Transform
        The transform (e.g. image registration transform) that maps Im to TrgIm.
    
        
    Outputs:
    *******
    
    SitkTransform : SimpleITK transform
        The 3D transformed image.
    """
    
    import SimpleITK as sitk
            
    if Transform == 'rigid':
        SitkTransform = sitk.Euler3DTransform()
        
        SitkTransform.SetParameters(TxMatrix[0:6])
        
    elif Transform == 'affine':
        SitkTransform = sitk.AffineTransform(FixIm.GetDimension())
        
        SitkTransform.SetParameters(TxMatrix[0:12])
        
    else:
        msg = "Only 'rigid' and 'affine' transforms are supported at present."\
              + f" Transform = {Transform} is not allowed."
        
        raise Exception(msg)
    
    return SitkTransform






def TransformImageSitk(MovIm, FixIm, Transform, Interp='NearestNeighbor', 
                       LogToConsole=False):
    """
    Transform a 3D SimpleITK image using SimpleITK.
    
    Inputs:
    ******
    
    MovIm : SimpleITK image 
        The 3D image to be transformed.
        
    FixIm : SimpleITK image
        The 3D image whose grid MovIm is to be transformed to.
        
    Transform : SimpleITK Transform
        The transform (e.g. image registration transform) that maps Im to TrgIm.
        
    Interp : string (optional; 'NearestNeighbor' by default)
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
        print(f'Results of TransformImageSitk():')
        
    
    if Interp == 'NearestNeighbor':
        interp = sitk.sitkLinear
    else:
        interp = sitk.sitkNearestNeighbor
    
    
    TxIm = sitk.Resample(MovIm, FixIm, Transform, interp)
    
    if LogToConsole:
        print('-'*120)
    
    return TxIm







def TransformLabImByRoi(LabImByRoi, F2SindsByRoi, TrgIm, SelxImFiltOrSitkTx, 
                       Interp='NearestNeighbor', 
                       ApplyPostTxBlur=True, PostTxVariance=(1,1,1),
                       ApplyPostTxBin=True, #VoxelVolumeRatio=1,
                       LogToConsole=False):
    """
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
        tranformed labelmap image(s).
        
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
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    from GeneralTools import UniqueItems
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformLabImByRoi():')
    
            
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
            print(f'\n   Image info for LabIm prior to transformation:')
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(LabIm, LogToConsole)
        
        """ Transform LabIm. """
        if 'ElastixImageFilter' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSelx(MovIm=LabIm, SelxImFilt=SelxImFiltOrSitkTx,
                                         Interp=Interp, 
                                         LogToConsole=LogToConsole)
        elif 'Transform' in str(type(SelxImFiltOrSitkTx)):
            TxLabIm = TransformImageSitk(MovIm=LabIm, FixIm=TrgIm, 
                                         Tx=SelxImFiltOrSitkTx, Interp=Interp, 
                                         LogToConsole=LogToConsole)
        
        times.append(time.time())
        Dtime = round(times[-1] - times[-2], 1)
        print(f'*Took {Dtime} s to transform the {r}^th labelmap image',
              'using the image registration transformation.\n')
        
        if LogToConsole:
            print(f'\n   Image info for LabIm after transformation:')
            
        PixID, PixIDTypeAsStr, UniqueValsPostTx,\
        F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
        
        
        if ApplyPostTxBlur:
            #""" The binary pixel array: """
            #PixArr_B, F2Sinds = Image2PixArr(TxLabIm)
            
            """ Gaussian blur TxLabIm (new 14/02/21) to reduce noise that 
            makes selection of an appropriate threshold for binarisation (below)
            difficult. """
            #TxLabIm = GaussianBlurImage(TxLabIm)
            TxLabIm = GaussianBlurImage(TxLabIm, Variance=PostTxVariance)
                
            times.append(time.time())
            Dtime = round(times[-1] - times[-2], 1)
            print(f'*Took {Dtime} s to Gaussian blur the {r}^th labelmap image.\n')
            
            if LogToConsole:
                print(f'\n   Image info for LabIm after Gaussian blur:')
                
            PixID, PixIDTypeAsStr, UniqueValsPostBlur,\
            F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
            
            #""" The non-binary pixel array following transformation: """
            #PixArr_NB, F2Sinds = Image2PixArr(TxLabIm)
            
        
        
            if ApplyPostTxBin:
                """Find suitable threshold value that approximately preserves
                the number of pre-transformed truth values: """
                Thresh = FindSuitableThresh(BinaryIm=LabIm, 
                                            NonBinaryIm=TxLabIm,
                                            LogToConsole=LogToConsole)
                
                
                """ Binary threshold TxLabIm. """
                #ThreshValue = ThreshLevel*max(UniqueVals)
                #ThreshValue = ThreshLevel                
                
                TxLabIm = BinaryThresholdImage(Im=TxLabIm, 
                                               #Thresh=ThreshPostTx)
                                               Thresh=Thresh)
                
                if LogToConsole:
                    print('\n   Image info for TxLabIm after binary',
                          #f'thresholding at {ThreshPostTx}:')
                          f'thresholding at {Thresh}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
                
                if LogToConsole:
                    print(f'\n   Image info for LabIm after binary',
                          'thresholding:')
                    
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(TxLabIm, LogToConsole)
                
        
        
        """ Ensure TxLabIm is a 32-bit unsigned integer (PixID = 5). """
        if PixID != 5: 
            if LogToConsole:
                print(f'\n   TxLabIm has PixelID = {PixID}'\
                      f'({PixIDTypeAsStr})).')
            
            """ Convert TxLabIm from float to 32-bit unsigned integer. """
            TxLabIm = ConvertImagePixelType(Image=TxLabIm, 
                                            NewPixelType='sitkUInt32')
            
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







