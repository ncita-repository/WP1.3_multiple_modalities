# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:56:14 2020

@author: ctorti
"""



# Import packages and functions:
#import os
#import SimpleITK as sitk
#import numpy as np
#import time

##import importlib
##import GeneralTools
##importlib.reload(GeneralTools)
##import DicomTools
##importlib.reload(DicomTools)

#from GeneralTools import ItemsUniqueToWithin
#from DicomTools import ImportDicoms





"""
******************************************************************************
******************************************************************************
DICOM IMAGE FUNCTIONS
******************************************************************************
******************************************************************************
"""



def ImportImage(DicomDir):
    """
    Import a DICOM series as a SimpleITK multi-dimensional image.
    
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
        print(f'\nSize = {Size} \nSpacings = {Spacings}',
              f'\nSliceThickness = {SliceThick} \nPositions = {Positions}',
              f'\nDirections = {Directions}')
    
    
    return Size, Spacings, SliceThick, Positions, Directions





def GetImageInfo(Image, LogToConsole=False):
    from ConversionTools import Image2PixArr
    import importlib
    import GeneralTools
    importlib.reload(GeneralTools)
    from GeneralTools import UniqueItems
    
    PixArr, F2Sinds = Image2PixArr(LabmapIm=Image)
    
    #print(f'\nPixArr.shape = {PixArr.shape}')
    #print(f'F2Sinds = {F2Sinds}')
    #print(f'\nPixArr = {PixArr}')
    #print(f'\ntype(PixArr) = {type(PixArr)}')
    
    UniqueVals = UniqueItems(Items=PixArr, IgnoreZero=False)
    
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
        print(f'      Conversion of Image to PixArr:')
        print(f'         PixArr shape = {PixArr.shape}')
        if len(UniqueVals) < 20:
            print(f'         There are {len(UniqueVals)} unique values in',
                  f'PixArr:')
            print(f'         {UniqueVals}')
        else:
            
            print(f'         There are {len(UniqueVals)} unique values in',
                  f'PixArr:')
            print(f'         {UniqueVals[:6]}...{UniqueVals[-6:-1]}')
        print(f'         There are {len(F2Sinds)} frames with slice indices:')
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







def CheckImageSpacing(DicomDir):
    
    Size, Spacings, ST, IPPs, Dirs = GetImageAttributes(DicomDir, 'sitk')
    
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







def InitialiseImage(RefImage):
    """
    Initialise an empty SimpleITK image with the same size, spacing, origin and
    direction as another image.  
    
    Inputs:
    ******
    
    RefImage : SimpleITK image
        The reference image whose attributes will be used for the initialised
        image.
        
        
    Outputs:
    *******
    
    NewImage : SimpleITK image
        An empty image with the same PixelID, Size, Spacing, Origin and 
        Direction as RefImage.
    """
    
    import SimpleITK as sitk
    
    NewImage = sitk.Image(RefImage.GetSize(), RefImage.GetPixelID())
    
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






def SumAllLabmapIms(LabmapIms):
    
    import SimpleITK as sitk
    
    # Initialise labelmap sum:
    LMImsSum = InitialiseImage(LabmapIms[0])
    
    for i in range(len(LabmapIms)):
        LMImsSum = AddImages(LMImsSum, LabmapIms[i])

    # Convert to numpy arrays:
    LMImsSumNpa = sitk.GetArrayFromImage(LMImsSum)
    
    return LMImsSumNpa


    







def ResampleImage(Image, RefImage, Interpolation, LogToConsole=False):
    """
    Resample a 3D SimpleITK image.  The image can be a 3D DICOM image or a
    3D labelmap image.
    
    Inputs:
    ******                      
        
    Image : SimpleITK image
        The 3D image to be resampled.
                           
    RefImage : SimpleITK image
        The 3D image reference image whose gridspace Image will be resampled to.
        
    Interpolation : string
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear' (or 'linear') 
        - 'BSpline' (or 'Bspline' or 'bspline')
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        - 'Gaussian' (or 'gaussian')
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    ResImage : SimpleITK image
        The resampled 3D image.

    
    Notes:
    *****
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """

    import SimpleITK as sitk
    
    # Define which interpolator to use:
    if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:    
        Interpolator = sitk.sitkNearestNeighbor
        
        OutputPixelType = sitk.sitkUInt64
    
    else:
        if 'inear' in Interpolation:
            Interpolator = sitk.sitkLinear
        elif 'pline' in Interpolation:
            Interpolator = sitk.sitkBSpline
        elif 'aussian' in Interpolation:
            Interpolator = sitk.sitkLabelGaussian
        else:
            msg = '"Interpolation" must be "Linear", "BSpline" or "Gaussian".'
            
            raise Exception(msg)
            
        OutputPixelType = sitk.sitkFloat32 
            
    #print('\nUsing', Interpolation, 'interpolation\n')
    
    Resampler = sitk.ResampleImageFilter()
    
    Resampler.SetReferenceImage(RefImage)
    
    Resampler.SetInterpolator(Interpolator)
    
    # Use the Identity transform:
    Resampler.SetTransform(sitk.Transform())
    
    # Use the Affine transform:
    #dim = Image.GetDimension()
    #Resampler.SetTransform(sitk.AffineTransform(dim))
    
    #print(f'\nImage.GetPixelIDValue() = {Image.GetPixelIDValue()}')
    
    Resampler.SetOutputSpacing(RefImage.GetSpacing())
    Resampler.SetSize(RefImage.GetSize())
    Resampler.SetOutputDirection(RefImage.GetDirection())
    Resampler.SetOutputOrigin(RefImage.GetOrigin())
    #Resampler.SetDefaultPixelValue(Image.GetPixelIDValue())
    Resampler.SetOutputPixelType(OutputPixelType)
    Resampler.SetDefaultPixelValue(0)
    
    #Resampler.LogToConsoleOn() # 'ResampleImageFilter' object has no attribute
    # 'LogToConsoleOn'
    #Resampler.LogToFileOn() # 'ResampleImageFilter' object has no attribute 
    # 'LogToFileOn'
    
    #print('\n')
    #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: print("\rProgress: {0:03.1f}%...".format(100*Resampler.GetProgress()),end=''))
    #Resampler.AddCommand(sitk.sitkProgressEvent, lambda: sys.stdout.flush())
    
    ResImage = Resampler.Execute(Image)
    
    #sitk.Show(ResImage)
    
    if LogToConsole:
        print('\nResampler:\n', Resampler.GetTransform())

    return ResImage
    #return ResImage, Resampler







def RegisterImages(FixIm, MovIm, Tx='affine', LogToConsole=False):
    """
    Register two 3D SimpleITK images using SimpleElastix.
    
    Inputs:
    ******
    
    FixIm : SimpleITK image 
        The 3D image that MovIm will be registered to.
        
    MovIm : SimpleITK image
        The 3D image that will be registered to FixIm.
        
    Tx : string (optional; 'affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
        
    LogToConsole : boolean (optional; False by default)
        Denotes whether intermediate results will be logged to the console.
        
        
    Outputs:
    *******
    
    RegIm : SimpleITK image
        The 3D registered image.
        
    RegImFilt : SimpleITK image filter
        Elastix image transformation filter used to transform MovIm to FixIm.
    """
    
    import SimpleITK as sitk
    import time
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Initiate RegImFilt:
    RegImFilt = sitk.ElastixImageFilter()
    RegImFilt.LogToConsoleOn() # <-- no output in Jupyter
    #RegImFilt.LogToConsoleOff() 
    RegImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    # Define the fixed and moving images:
    RegImFilt.SetFixedImage(FixIm)
    RegImFilt.SetMovingImage(MovIm)
    
    # Get the default parameter map template for the chosen transformation:
    #RegImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    #ParamMap = sitk.GetDefaultParameterMap('affine')
    ParamMap = sitk.GetDefaultParameterMap(Tx)
    
    # Re-assign some parameters:
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['WriteIterationInfo'] = ['true']
    ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['UseDirectionCosines'] = ['true']
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0'] 
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    # Set the parameter map:
    RegImFilt.SetParameterMap(ParamMap)
    
    if True:#LogToConsole:
        print('\nPerforming registration...')
        
    #RegImFilt.AddCommand(sitk.sitkProgressEvent, 
    #                     lambda: print("\rProgress: {0:03.1f}%...".format(100*RegImFilt.GetProgress()),end=''))
    # 'ElastixImageFilter' object has no attribute 'AddCommand'
    
    # Register the 3D images:
    RegImFilt.Execute()
    
    # Get the registered image:
    RegIm = RegImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'Took {Dtime} s to register the 3D image stacks.')
        #print('\n', RegImFilt.GetTransform()) # 'ElastixImageFilter' object has 
        ## no attribute 'GetTransform'
        #print('\n', RegImFilt.PrintParameterMap()) # <SimpleITK.SimpleITK.ElastixImageFilter; 
        # proxy of <Swig Object of type 'itk::simple::ElastixImageFilter::Self *' 
        # at 0x00000220072E8B70> >
    
    
    return RegIm, RegImFilt
    #return RegIm, RegImFilt, ParamMap






def TransformImage(Im, RegImFilt, Interpolation='Default'):
    """
    Transform a 3D SimpleITK image using SimpleElastix.
    
    Inputs:
    ******
    
    Im : SimpleITK image 
        The 3D image to be transformed.
        
    RegImFilt : SimpleITK image filter
        Elastix image transformation filter used to perform image registration.
        
    Interpolation : string (optional; 'Default' by default)
        The interpolator to be used for resampling.  At present the only two
        accepatable values are:
        - 'Default' (or 'default') which leaves unchanged whatever interpolator
        has been set in RegImFilt
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour')
        
        
    Outputs:
    *******
    
    TxIm : SimpleITK image
        The 3D transformed image.
        
        
    Notes:
    *****
    
    https://github.com/SuperElastix/SimpleElastix/issues/409
    """
    
    import SimpleITK as sitk
            
            
    # Initiate TransformixImageFilter:
    TxImFilt = sitk.TransformixImageFilter()
    #TxImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TxImFilt.LogToConsoleOff() # <-- no output in Jupyter
    TxImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    TxImFilt.SetMovingImage(Im)
    
    # Get the transform parameter map. Start by getting the transfer parameter 
    # mapused for the registration (= ElastixParamMap):
    #TxParamMap = RegImFilt.GetParameterMap() 
    
    # Set the parameter map:
    #TxImFilt.SetTransformParameterMap(RegImFilt.GetTransformParameterMap())
    TxMap = RegImFilt.GetTransformParameterMap()
    
    #print(f'\nTxMap[0]["ResampleInterpolator"] = {TxMap[0]["ResampleInterpolator"]}')
    
    if ('earest' and 'eighbo' or 'nn' or 'NN') in Interpolation:
        TxMap[0]["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]
        
    TxImFilt.SetTransformParameterMap(TxMap)
    
    # Set up for computing deformation field:
    #TxImFilt.ComputeDeformationFieldOn()
    
    # Transform Im:
    TxImFilt.Execute()
    
    TxIm = TxImFilt.GetResultImage()
    
    return TxIm





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
    BinThreshImFilt.DebugOn()
    
    #BinThreshImFilt.Execute()
    #BinIm = BinThreshImFilt.GetOutput()
    
    BinIm = BinThreshImFilt.Execute(Im)
    
    #sitk.Show(Im)
    #sitk.Show(BinIm)
    
    return BinIm





def GaussianBlurImage(Im, Variance=(1.0, 1.0, 1.0), LogToConsole=False):
    """
    Gaussian blur a 3D SimpleITK image.
    
    Inputs:
    ******
    
    Im : SimpleITK image 
        The 3D image to be blurred.
        
    Variance : tuple of floats (optional; (1.0, 1.0, 1.0) by default)
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








def ResampleLabmapImByRoi_OLD(LabmapImByRoi, C2SindsByRoi, SrcImage, TrgImage, 
                          Variance, ThreshLevel, LogToConsole=False, 
                          Interpolation='NearestNeighbor'):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled.
    
    C2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each contour) of slice numbers that 
        correspond to each contour.
        
    SrcImage : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of the 
        labelmap images in LabelmapImByRoi.
    
    TrgImage : SimpleITK image
        The Target 3D image whose gridspace the Source labelmap images will be
        resampled to.
    
    Variance : tuple of floats (optional; (1.0, 1.0, 1.0) by default)
        The variance along all dimensions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Interpolation : string (optional; 'NearestNeighbor' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear' (or 'linear') 
        - 'BSpline' (or 'Bspline' or 'bspline')
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour'; Default value)
        - 'Gaussian' (or 'gaussian')
    
    
    Outputs:
    *******
    
    ResLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of resampled 3D labelmap images.
        
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of resampled pixel arrays (converted from each
        labelmap in ResLabmapImByRoi).
    
    ResF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.

    
    Notes:
    *****
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. F2Sinds = []) even if there
    were contours (i.e. if SrcC2SindsByRoi[r] != []). If this is the case try 
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
    
    import numpy as np
    from ConversionTools import Image2PixArr, ConvertImagePixelType
    
    """ Which Workaround to use. """
    #Workaround = 1
    Workaround = 2
    #Workaround = 0 # i.e. no workaround
        
        
    if LogToConsole:
        print(f'\nResults of ResampleLabmapImByRoi():')
        
    ResLabmapImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabmapImByRoi)):
        """ Resample LabmapImByRoi to the Target image's grid to get the 
        resampled pixel array. """
        ResLabmapIm = ResampleImage(Image=LabmapImByRoi[r], RefImage=TrgImage, 
                                    Interpolation=Interpolation)
        
        ResLabmapImByRoi.append(ResLabmapIm)
        
        """ Convert ResLabelmapIm to a pixel array. """
        ResPixArr, ResCStoSliceInds = Image2PixArr(LabmapIm=ResLabmapIm)
        
        ResPixArrByRoi.append(ResPixArr)
        
        if LogToConsole:
            print(f'   Image info for ResLabmapImByRoi[{r}] after resampling',
                  f'using {Interpolation} interpolation:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        ResF2SindsByRoi.append(F2Sinds)
        
        
    """ Check for the presence of any empty labelmaps (i.e. empty frame-to-slice
    indices) where not expected (a list of empty contour-to-slice indices, e.g.
    C2SindsByRoi[1], will be empty if there were no contours of interest for 
    that ROI). """
    
    for r in range(len(ResLabmapImByRoi)):
        if C2SindsByRoi[r] != [] and ResF2SindsByRoi[r] == [] and Workaround:
    
            if LogToConsole:
                print(f'   C2SindsByRoi[{r}] = {C2SindsByRoi[r]}')
                print('   but')
                print(f'   ResF2SindsByRoi[{r}] = {ResF2SindsByRoi[r]}')
                print('   i.e. the resampled image has no non-zero frames',
                      f'(e.g. due to aliasing). Try a workaround.')
            
            """ Convert LabmapImByRoi[r] from 32-bit unsigned integer to float.
            Note:  Result of resampling results in empty labelmap unless
            image is converted to float prior to resampling. """
            LabmapIm = ConvertImagePixelType(Image=LabmapImByRoi[r], 
                                             NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print(f'   Image info for ResLabmapImByRoi[{r}] after'
                      'converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
        
            if Workaround == 1:
                if LogToConsole:
                    print(f'   Trying workaround #{Workaround}...')
                
                    
                """ Resample LabmapIm to the Target image's grid using a
                Gaussian interpolation. """
                interp = 'Gaussian'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after'
                          f'resampling using a {interp} interpolation:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                    
                    
                """ Convert back to unsigned int and apply binary threshold. 
                But only proceed if there are non-zero frames in the resampled 
                image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
            
            if Workaround == 2:
                if LogToConsole:
                    print(f'   Trying workaround #{Workaround}...')
                      
                SrcSpacings = SrcImage.GetSpacing()
                TrgSpacings = TrgImage.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
                                                                      TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM 
                should be at least 2x the voxel spacings of TrgIm. For 
                symmetry considerations, it makes sense for one Source 
                voxel width to be blurred to three Target voxels widths.
                    
                    sigma = (3*TrgSpacings)/2.355
                """
                
                if LogToConsole:
                    print(f'   SpacingsRatio = {SpacingsRatio}')
                      
                """ Gaussian blur the image. """
                LabmapIm = GaussianBlurImage(Im=LabmapIm, Variance=Variance)
                
                # Use the RecursiveGaussian image filter:
                #LabmapImToCopy = RecursiveGaussianBlurImage(LabmapIm, 
                #                                            Sigma=3,
                #                                            Direction=2)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'applying Gaussian blur:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                
                if LogToConsole:
                    print(f'   ResLabmapImByRoi[{r}] prior to {interp}',
                          'interpolation:')
                    print(f'      LabmapIm.GetSize() = {LabmapIm.GetSize()}')
                    print(f'      LabmapIm.GetSpacing() = {LabmapIm.GetSpacing()}')
                    print(f'      TrgImage.GetSize() = {TrgImage.GetSize()}')
                    print(f'      TrgImage.GetSpacing() = {TrgImage.GetSpacing()}')
                    
                """ Linearly resample LabmapIm to the Target image's grid. """
                interp = 'Linear'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'linear resampling:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                
                """ Binary threshold the image. """
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                ResLabmapIm = BinaryThresholdImage(Im=ResLabmapIm, 
                                                   Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'   Image info for ResLabmapImByRoi[{r}] after',
                          'binary thresholding at',
                          f'{ThreshLevel}*{max(UniqueVals)} = {ThreshValue}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
    
        
        
        """ Ensure that ResLabmapIm is a 32-bit unsigned integer. """
        if not PixID == 5:
            if LogToConsole:
                print(f'   ResLabmapImByRoi[{r}] has PixelID = {PixID}',
                      f'({PixIDTypeAsStr})). Converting to unsigned 32-bit',
                      'integer (sitkUInt32)..')
            
            """ Convert ResLabmapIm from float to 32-bit unsigned integer. """
            ResLabmapIm = ConvertImagePixelType(Image=ResLabmapIm, 
                                                NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print(f'   Image info for ResLabmapImByRoi[{r}] after',
                      'converting to 32-bit unsigned int:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        """ Convert ResLabmapIm to a pixel array. """
        ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
        
        if LogToConsole:
            print(f'  After converting ResLabmapImByRoi[{r}] to a pixel array:')
            print(f'   ResPixArr.shape = {ResPixArr.shape}')
            print(f'   ResF2Sinds = {ResF2Sinds}')
            
            #PixID, PixIDTypeAsStr, UniqueVals,\
            #F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
            
            #for i in range(ResPixArr.shape[0]):
            #    maxval = np.amax(ResPixArr[i])
            #    minval = np.amin(ResPixArr[i])
            #    
            #    print(f'      Frame {i} has max = {maxval}, min = {minval}')
                
        """ Replace ResLabmapImByRoi[r], ResPixArrByRoi[r] and 
        ResF2SindsByRoi[r] with ResLabelmapIm, ResPixArr and ResC2Sinds. """
        ResLabmapImByRoi[r] = ResLabmapIm
        ResPixArrByRoi[r] = ResPixArr
        #ResF2SindsByRoi[r] = F2Sinds
        ResF2SindsByRoi[r] = ResF2Sinds
        
    return ResLabmapImByRoi, ResPixArrByRoi, ResF2SindsByRoi






def ResampleLabmapImByRoi(LabmapImByRoi, F2SindsByRoi, SrcImage, TrgImage, 
                          Variance, ThreshLevel, LogToConsole=False, 
                          Interpolation='NearestNeighbor'):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled.
    
    F2SindsByRoi : list of a list of integers
        List (for each ROI) of a list (for each frame) of slice numbers that 
        correspond to each frame in the labelmap images.
        
    SrcImage : SimpleITK image
        The Source 3D image whose gridspace matches the gridspace of the 
        labelmap images in LabelmapImByRoi.
    
    TrgImage : SimpleITK image
        The Target 3D image whose gridspace the Source labelmap images will be
        resampled to.
    
    Variance : tuple of floats (optional; (1.0, 1.0, 1.0) by default)
        The variance along all dimensions.
        
    ThreshLevel : float (optional; 0.75 by default)
        The proportion of the maximum value used to perform binary thresholding
        on an image if used.  The threshold value will be ThreshLevel*M, where
        M is the maximum intensity value of the image.
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
    
    Interpolation : string (optional; 'NearestNeighbor' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear' (or 'linear') 
        - 'BSpline' (or 'Bspline' or 'bspline')
        - 'NearestNeighbor' (or 'NearestNeighbour' or 'nearestneighbor' or 
        'nearestneighbour'; Default value)
        - 'Gaussian' (or 'gaussian')
    
    
    Outputs:
    *******
    
    ResLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of resampled 3D labelmap images.
        
    ResPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of resampled pixel arrays (converted from each
        labelmap in ResLabmapImByRoi).
    
    ResF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in ResPixArrByRoi.

    
    Notes:
    *****
    
    If resampling to a smaller grid size aliasing can occur, leading to the
    appearent disappearance of segmentations (i.e. F2Sinds = []) even if there
    were contours (i.e. if SrcC2SindsByRoi[r] != []). If this is the case try 
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
    
    """ Which Workaround to use. """
    #Workaround = 1
    Workaround = 2
    #Workaround = 0 # i.e. no workaround
        
        
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of ResampleLabmapImByRoi():')
        
    ResLabmapImByRoi = []
    ResPixArrByRoi = []
    ResF2SindsByRoi = []
    
    for r in range(len(LabmapImByRoi)):
        LabmapIm = LabmapImByRoi[r]
        F2Sinds = F2SindsByRoi[r]
        
        """ Resample LabmapIm to the Target image's grid to get the resampled
        pixel array. """
        ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                    Interpolation=Interpolation)
        
        """ Convert ResLabmapIm to a pixel array. """
        ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
        
        if LogToConsole:
            print(f'\n   r = {r}')
            print('\n   Image info for ResLabmapIm after resampling using',
                  f'{Interpolation} interpolation:')
        
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        """ Is ResF2Sinds empty and not expected to be? 
        Note:
            If F2Sinds is not empty, ResF2Sinds should not be empty either.
        
            F2Sinds will be empty if there were no segmentations/contours of 
            interest for the r^th ROI. In this case an empty ResF2Sinds is 
            acceptable. """
        
        #if True:
        if F2Sinds != [] and ResF2Sinds == [] and Workaround:
            WorkaroundWasDone = True
            
            if LogToConsole:
                print(f'\n   F2Sinds = {F2Sinds}')
                print('   but')
                print(f'   ResF2Sinds = {ResF2Sinds}')
                print('   i.e. the resampled image has no non-zero frames',
                      f'(e.g. due to aliasing). Try a workaround.')
            
            """ Convert LabmapIm from 32-bit unsigned integer to float.
            Note:  Result of resampling results in empty labelmap unless
            image is converted to float prior to resampling. """
            LabmapIm = ConvertImagePixelType(Image=LabmapIm, 
                                             NewPixelType='sitkFloat32')
            
            if LogToConsole:
                print('\n   Image info for LabmapIm after converting to float:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
        
            if Workaround == 1:
                if LogToConsole:
                    print(f'\n   Trying workaround #{Workaround}...')
                
                    
                """ Resample LabmapIm to the Target image's grid using a
                Gaussian interpolation. """
                interp = 'Gaussian'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print('\n   Image info for ResLabmapIm after resampling'
                          f'using a {interp} interpolation:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                    
                    
                """ Convert back to unsigned int and apply binary threshold. 
                But only proceed if there are non-zero frames in the resampled 
                image. """
                
                if not F2Sinds:
                    Workaround = 2 # try Workaround #2
            
            if Workaround == 2:
                if LogToConsole:
                    print(f'\n   Trying workaround #{Workaround}...')
                      
                SrcSpacings = SrcImage.GetSpacing()
                TrgSpacings = TrgImage.GetSpacing()
                
                SpacingsRatio = tuple(tup1/tup0 for tup0, tup1 in zip(SrcSpacings, 
                                                                      TrgSpacings))
                
                """ 
                The FWHM of the Gaussian is:
                    FWHM = 2*sqrt(2*ln(2))*sigma ~= 2.355*sigma
                    
                Due to the Nyquist-Shannon sampling theorem, the FWHM 
                should be at least 2x the voxel spacings of TrgIm. For 
                symmetry considerations, it makes sense for one Source 
                voxel width to be blurred to three Target voxels widths.
                    
                    sigma = (3*TrgSpacings)/2.355
                """
                
                if LogToConsole:
                    print(f'\n   SpacingsRatio = {SpacingsRatio}')
                      
                """ Gaussian blur the image. """
                LabmapIm = GaussianBlurImage(Im=LabmapIm, Variance=Variance)
                
                # Use the RecursiveGaussian image filter:
                #ResLabmapIm = RecursiveGaussianBlurImage(ResLabmapIm, 
                #                                         Sigma=3,
                #                                         Direction=2)
                
                if LogToConsole:
                    print('\n   Image info for LabmapIm after applying',
                          'Gaussian blur:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(LabmapIm, LogToConsole)
                
                #if not F2Sinds:
                #    Workaround = 3 # try Workaround #3
                
                if LogToConsole:
                    print(f'\n   LabmapIm prior to resampling:')
                    print(f'      LabmapIm.GetSize() = {LabmapIm.GetSize()}')
                    print(f'      LabmapIm.GetSpacing() = {LabmapIm.GetSpacing()}')
                    print(f'      TrgImage.GetSize() = {TrgImage.GetSize()}')
                    print(f'      TrgImage.GetSpacing() = {TrgImage.GetSpacing()}')
                    
                """ Linearly resample LabmapIm to the Target image's grid. """
                interp = 'Linear'
                ResLabmapIm = ResampleImage(Image=LabmapIm, RefImage=TrgImage, 
                                            Interpolation=interp)
                
                if LogToConsole:
                    print(f'\n   Image info for ResLabmapIm after {interp}',
                          'resampling:')
                    
                    PixID, PixIDTypeAsStr, UniqueVals,\
                    F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                
                """ Binary threshold the image. """
                ThreshValue = ThreshLevel*max(UniqueVals)
                
                ResLabmapIm = BinaryThresholdImage(Im=ResLabmapIm, 
                                                   Thresh=ThreshValue)
                
                if LogToConsole:
                    print(f'\n   Image info for ResLabmapIm after binary',
                          f'thresholding at {ThreshLevel}*{max(UniqueVals)} =',
                          f'{ThreshValue}:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        else:
            WorkaroundWasDone = False
    
        
        
        """ Ensure that ResLabmapIm is a 32-bit unsigned integer. """
        if PixID == 5:
            ConversionWasDone = False
        else:
            ConversionWasDone = True
            
            if LogToConsole:
                print(f'\n   ResLabmapIm has PixelID = {PixID}',
                      f'({PixIDTypeAsStr})).')
            
            """ Convert ResLabmapIm from float to 32-bit unsigned integer. """
            ResLabmapIm = ConvertImagePixelType(Image=ResLabmapIm, 
                                                NewPixelType='sitkUInt32')
            
            if LogToConsole:
                print('\n   Image info for ResLabmapIm after converting to',
                      '32-bit unsigned int:')
                
                PixID, PixIDTypeAsStr, UniqueVals,\
                F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
        
        
        
        if WorkaroundWasDone or ConversionWasDone:
            """ Convert ResLabmapIm to a pixel array. """
            ResPixArr, ResF2Sinds = Image2PixArr(LabmapIm=ResLabmapIm)
            
            if LogToConsole:
                print(f'\n  After converting ResLabmapIm to a pixel array:')
                print(f'   ResPixArr.shape = {ResPixArr.shape}')
                print(f'   ResF2Sinds = {ResF2Sinds}')
                
                #PixID, PixIDTypeAsStr, UniqueVals,\
                #F2Sinds = GetImageInfo(ResLabmapIm, LogToConsole)
                
                #for i in range(ResPixArr.shape[0]):
                #    maxval = np.amax(ResPixArr[i])
                #    minval = np.amin(ResPixArr[i])
                #    
                #    print(f'      Frame {i} has max = {maxval}, min = {minval}')
                
        
        ResLabmapImByRoi.append(ResLabmapIm)
        ResPixArrByRoi.append(ResPixArr)
        ResF2SindsByRoi.append(ResF2Sinds)
        
        #if LogToConsole:
        #    print('\n')
    
    if LogToConsole:
            print('-'*120)
        
    return ResLabmapImByRoi, ResPixArrByRoi, ResF2SindsByRoi







def TransformLabmapImByRoi(LabmapImByRoi, FixImage, MovImage, 
                           LogToConsole=False):
    """
    Resample a list 3D SimpleITK images representing labelmaps. A 
    NearestNeighbor interpolation will be applied initially.  If aliasing 
    effects result in a empty labelmap, the labelmap will be Gaussian blurred,
    linearly resampled and binary thresholded.  See Notes for more info.
    
    Inputs:
    ******                      
        
    LabmapImByRoi : list of SimpleITK images
        A list (for each ROI) of 3D labelmap images to be resampled. The 
        labelmaps share the same image attributes as MovImage.
    
    FixImage : SimpleITK image
        The 3D image that MovImage will be registered to (e.g. SrcImage).
    
    MovImage : SimpleITK image
        The 3D image that will be registered to FixImage (e.g. TrgImage).
    
    LogToConsole : boolean (optional; False by default)
        Denotes whether some results will be logged to the console.
        
    
    Outputs:
    *******
    
    TxLabmapImByRoi : list of SimpleITK images
        The list (for each ROI) of the 3D labelmap images in LabmapImByRoi
        transformed to the FixImage grid space.
        
    TxPixArrByRoi : list of Numpy arrays
        The list (for each ROI) of pixel arrays (converted from each labelmap 
        in TxLabmapImByRoi).
    
    TxF2SindsByRoi : list of a list of integers
        The list (for each ROI) of a list (for each frame) of the
        frame-to-slice indices in TxPixArrByRoi.
    """
    
    from ConversionTools import Image2PixArr
    from GeneralTools import UniqueItems
    
    if LogToConsole:
        print('\n\n', '-'*120)
        print(f'Results of TransformLabmapImByRoi():')
    
    # Register the images using an affine transformation:
    RegIm, RegImFilt = RegisterImages(FixIm=FixImage, MovIm=MovImage, 
                                      Tx='affine', LogToConsole=LogToConsole)
        
    # Transform LabmapImByRoi to the Target image's grid to get the resampled
    # pixel array:
    TxLabmapImByRoi = []
    TxPixArrByRoi = []
    TxF2SindsByRoi = [] 
    
    for r in range(len(LabmapImByRoi)):
        # Transform SrcLabmapIm using RegImFilt:
        TxLabmapIm = TransformImage(Im=LabmapImByRoi[r], RegImFilt=RegImFilt,
                                    Interpolation='Nearestneighbor')
        
        TxLabmapImByRoi.append(TxLabmapIm)
        
        # Convert TxLabmapIm to a pixel array
        TxPixArr, TxC2Sinds = Image2PixArr(LabmapIm=TxLabmapIm)
        
        TxPixArrByRoi.append(TxPixArr)
        TxF2SindsByRoi.append(TxC2Sinds)
            
        PixID, PixIDTypeAsStr, UniqueVals,\
        F2Sinds = GetImageInfo(TxLabmapIm, LogToConsole)
        
        #TxF2SindsByRoi.append(F2Sinds)
        
        
        if LogToConsole:
            unique = UniqueItems(Items=TxPixArr, IgnoreZero=False)
            print(f'\nThere are {len(unique)} unique items in TxPixArr',
                  f'after transforming LabmapImByRoi[{r}]')
            
            F = TxPixArr.shape[0]
            
            print(f'\nThe segmentation has been registered to {F} frames/',
                  f'contours: {TxC2Sinds}.')
            
    if LogToConsole:
        print('-'*120)
        
    return TxLabmapImByRoi, TxPixArrByRoi, TxF2SindsByRoi
