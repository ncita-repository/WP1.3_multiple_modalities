# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:03:15 2020

@author: ctorti
"""


"""
Function:
    GetImageAttributes()
    
Purpose:
    Get Origin, Direction cosines and pixel spacings using either pydicom or
    SimpleITK.


Input:
    DicomDir   - Directory containing DICOMs 
    
    Package    - Package to use; acceptable inputs are:
                 - 'pydicom'
                 - 'sitk' (SimpleITK)
    
    
    
Returns:
    Origin     - The ImagePositionPatient (pydicom) or Origin (sitk)
                
                e.g. [x0, y0, z0]
    
    Directions - The direction cosine along x (rows), y (columns) and z (slices)
                 
                e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
                
                SimpleITK returns all returns all direction cosines.
                For pydicom, the cross product of the x and y vectors are used
                to obtain the z vector.
    
    Spacings   - The pixel spacings along x, y and z 
    
                e.g. [di, dj, dk]
                
                SimpleITK returns all pixel spacings.
                For pydicom, the SliceThickness is added to the PixelSpacing.
                
    Dimensions - The dimensions of the image stack
    
"""



def GetImageAttributes(DicomDir, Package):
    
    if Package=='sitk':
        # Import packages:
        import SimpleITK as sitk
    
        """ Get the Image Plane Attributes from SimpleITK image: """
        # Read in the 3D stack of DICOMs:
        SitkReader = sitk.ImageSeriesReader()
        SitkNames = SitkReader.GetGDCMSeriesFileNames(DicomDir)
        SitkReader.SetFileNames(SitkNames)
        SitkIm = SitkReader.Execute()
        SitkIm = sitk.Cast(SitkIm, sitk.sitkFloat32)
        
        Origin = SitkIm.GetOrigin() 
        # e.g. (-104.92378234863281, -152.4906463623047, -9.22148609161377)
        
        Directions = SitkIm.GetDirection() 
        # e.g. (0.9972648738498597, 1.8634260454120088e-10, 0.07391056342109259,
        #       -0.022262997957353512, 0.9535561337886633, 0.3003915089278791,
        #       -0.07047787104598309, -0.301215370979001, 0.9509480374756598)
        
        Spacings = SitkIm.GetSpacing() # e.g. (0.8984375, 0.8984375, 5.0)
        
        Dimensions = SitkIm.GetSize() # e.g. (192, 256, 35)
        
        
    elif Package=='pydicom':
        # Import packages:
        import numpy as np
        from GetDicoms import GetDicoms
        
        """ Get the Image Plane Attributes from Pydicom image: """
        # Get the DICOMs:
        DicomFpaths, Dicoms = GetDicoms(DirPath=DicomDir, 
                                        SortMethod='slices', 
                                        Debug=False)
        
        Origin = [float(item) for item in Dicoms[0].ImagePositionPatient]
        # e.g. [-115.14457762101, -124.83148479744, -10.815127450171]
        
        Directions = [float(item) for item in Dicoms[0].ImageOrientationPatient]
        # e.g. [0.99215282334066, -0.1042593858079, -0.0690127206254, 
        #       0.10603055257536, 0.9941080735439, 0.02250911005349]
        
        # Get the direction vector along z using the cross product of the x and 
        # y vectors:
        zDir = np.cross(Directions[0:3], Directions[3:])
        # e.g. [ 0.06625932 -0.02964993  0.99736181]
        
        # Append zDir to Dirs:
        Directions.append(zDir[0])
        Directions.append(zDir[1])
        Directions.append(zDir[2])
        
        Spacings = [float(item) for item in Dicoms[0].PixelSpacing]
        # e.g. [0.9375, 0.9375]
        
        SliceThick = float(Dicoms[0].SliceThickness)
        # e.g. 5.0
        
        # Append SliceThick to Spacings:
        Spacings.append(SliceThick)
        
        Dimensions = [int(Dicoms[0].Rows), int(Dicoms[0].Columns), len(Dicoms)]
    
    
    return Origin, Directions, Spacings, Dimensions
    
    