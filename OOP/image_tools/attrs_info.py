# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:16:21 2021

@author: ctorti
"""


""" Functions that return attributes or info on a SimpleITK Image. """

from importlib import reload

import image_tools.imports
reload(image_tools.imports)
import dicom_tools.imports
reload(dicom_tools.imports)
import general_tools.general
reload(general_tools.general)


import numpy as np
from image_tools.imports import import_im
#from image_tools.operations import im_min, im_max
from dicom_tools.imports import import_dcms
from general_tools.general import get_items_unique_to_within, get_unique_items
from conversion_tools.pixarrs_ims import im_to_pixarr
from image_tools.operations import im_min, im_max

def get_im_attrs(dicomDir, package='pydicom', p2c=False):
    """
    Get image size, voxel spacings, image positions (IPPs), and direction 
    cosines.
    
    Parameters
    ----------
    dicomDir : str
        Directory containing DICOMs.
    package : str, optional ('pydicom' by default)
        Package to use; acceptable inputs are:
         - 'pydicom' (default)
         - 'sitk' (SimpleITK)
    p2c : bool, optional (False by default)
        If True results will be printed to the console.
    
    Returns
    -------
    size : list of int
        The size/dimensions of the 3D image along x, y and z, 
        e.g. [Columns, Rows, NumOfSlices]
    spacings : list of float
        The pixel spacings along x, y and z (= SliceThickness appended to 
        PixelSpacing), e.g. [di, dj, dk]              
    positions : list of list of float
        The ImagePositionPatient of all slices in the DICOM series, 
        e.g. [[x0_0, y0_0, z0_0], [x1_0, y1_0, z1_0], ...]
    directions : list of float
        The direction cosine along x (rows), y (columns) and z (slices).
    warnings : list of str
        List of any warnings.
          
    Notes
    -----
    
    1) When using Pydicom the cross product of the x and y direction cosines 
    from ImageOrientationPatient is used to obtain the direction cosine along x.
    
    2) There is a sign difference between the cross-terms obtained from Pydicom
    and SimpleITK. If the directions obtained from Pydicom are:
        [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz],
        
    the directions obtained from SimpleITK would be:
        [Xx, -Xy, -Xz, -Yx, Yy, -Yz, -Zx, -Zy, Zz].
        
    3) When using Pydicom the directions will not be correct for coronal or 
    sagittal image stacks.    
    """
    
    warnings = []
    
    if package=='sitk':
        # Import the image:
        im = import_im(dicomDir)
        
        size = im.GetSize()
        
        spacings = im.GetSpacing()
        
        #origin = im.GetOrigin() 
        
        positions = []
        
        for i in range(size[2]):
            positions.append(im.TransformIndexToPhysicalPoint((0, 0, i)))
        
        directions = im.GetDirection() 
        
        """
        Since ImageSeriesReader() (invoked in ImportImage()) doesn't allow for
        calling ReadImageInformation(), it's not straightforward to get 
        metadata using sitk.  So instead get the FOR UID and SliceThickness 
        using Pydicom.
        """
        # Import the DICOMs:
        dicoms = import_dcms(dicomDir=dicomDir, sortMethod='slices', p2c=False)
        
        slcThick = float(dicoms[0].SliceThickness)
        
    elif package=='pydicom':
        """
        NOTE:
            This only gives the correct directions for an axial stack.
            The components need to be switched around the coronal or sagittal.
        """
        
        # Import the DICOMs:
        dicoms = import_dcms(dicomDir=dicomDir, sortMethod='slices', p2c=False)
        
        size = [int(dicoms[0].Columns), int(dicoms[0].Rows), len(dicoms)]
        
        #origin = [float(item) for item in dicoms[0].ImagePositionPatient]
        
        positions = []
    
        for dicom in dicoms:
            positions.append([float(item) for item in dicom.ImagePositionPatient])
        
        # Check if all IPPs are unique.
        """
        Some scans ("ep2ddiff3scantracep2") have repeated IPPs.
        """
        
        # Get unique positions:
        #uniquePositions = list(set(positions))
        #uniquePositions = np.unique(np.array(positions))
        uniquePositions = np.array(list(set(tuple(p) for p in positions)))
        
        P = len(positions)
        U = len(uniquePositions)
        
        if U != P:
            msg = f'\nWarning:  There are only {U} unique IPPs within the '\
                  + f'list of {P} IPPs.\n'
                  
            warnings.append(msg)
                  
            print(msg)
        
        IOP = [float(item) for item in dicoms[0].ImageOrientationPatient]
        
        # Get the direction vector along z using the cross product of the x and 
        # y vectors:
        zDir = np.cross(IOP[0:3], IOP[3:])
        
        # Append zDir to IOP:
        directions = IOP
        directions.extend(zDir)
        
        spacings = [float(item) for item in dicoms[0].PixelSpacing]
        
        slcThick = float(dicoms[0].SliceThickness)
        
        """ 
        Don't use slcThick for the z-component of the voxel spacing. 
        Instead use the difference between slices from the IPP.
        """
        # Append slcThick to spacings:
        #spacings.append(slcThick)
        
        # Compute the slice spacings along the z-direction:
        #Dz = np.diff(np.array(positions), axis=0)[:,2] # <-- WRONG since it doesn't account for IOP
        
        ##UniqueDz = list(set(Dz))
        #UniqueDz = ItemsUniqueToWithin(Dz)
        
        # Compute the vector lengths of the IPP differences:
        vectLengths = []

        for i in range(len(positions) - 1):
            vector = np.array(positions[i + 1]) - np.array(positions[i])
            
            vectorL = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
            
            vectLengths.append(vectorL)
        
        # Are the list of unique vector lengths equal to within 1% of the
        # maximum value?
        thresh = max(vectLengths)/100
        
        uniqueVectLengths = get_items_unique_to_within(vectLengths, thresh)
        
        # Append uniqueVectLengths to spacings:
        spacings.append(uniqueVectLengths[0])
        
        if len(uniqueVectLengths) > 1:
            #msg = f'\nWarning:  The voxel spacings along the scan direction '\
            #      + f'are non-uniform with the following unique values:\n'\
            #      + f'{uniqueVectLengths}'
            
            msg = '\nWarning:  The voxel spacings along the scan direction '\
                  + f'are not to within 1% of each other:\n{uniqueVectLengths}'
            warnings.append(msg)
            print(msg)
        
        # Compare UniqueDz to slcThickness:    
        #epsilon = 1e-5
        #
        #if abs(slcThick - UniqueDz[0]) > epsilon:
        #    print('\nNote:')
        #    print('    The slice thickness obtained from the IPPs do not',
        #          'agree with the DICOM tag slcThickness:')
        #    print(f'      d(IPP[2])      = {UniqueDz[0]}')
        #    print(f'      SliceThickness = {slcThick}')
        
    if p2c:
        print(f'\nsize = {size} \nspacings = {spacings}',
              f'\nslcThick = {slcThick} \npositions = {positions}',
              f'\ndirections = {directions}')
    
    return size, spacings, slcThick, positions, directions, warnings

def get_im_info(image, p2c=False):
    """
    Return some info about an image.

    Parameters
    ----------
    image : SimpleITK Image
    p2c : bool, optional
        If True intermediate results will be printed to the console. The 
        default is False.

    Returns
    -------
    pixID : int
        The pixel type.
    pixIDtypeAsStr : str
        The pixel type as a string.
    uniqueVals : list of ints or floats
        List of values in image that are unique.
    f2sInds : list of ints
        List of integers that indexes each non-zero frame in image, i.e. 
        frame-to-slice indeces.
        
    Notes
    -----
    For pixel type see:
    https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1Image.html#a780923bec6a71b14be35d81069f343f6
    
    For pixel type as a string see:
    https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1Image.html#ad195963d0b257560819b833b9cfe18d6
    """
    
    #from conversion_tools.pixarrs_ims import im_to_pixarr
    #from image_tools.operations import im_min, im_max
    
    pixarr, f2sInds = im_to_pixarr(image)
    
    #print(f'\npixarr.shape = {pixarr.shape}')
    #print(f'f2sInds = {f2sInds}')
    #print(f'\npixarr = {pixarr}')
    #print(f'\ntype(pixarr) = {type(pixarr)}')
    
    uniqueVals = get_unique_items(items=pixarr, ignoreZero=False)
    
    #print(f'uniqueVals = {uniqueVals}')
    #print(f'type(uniqueVals) = {type(uniqueVals)}')
    
    pixID = image.GetPixelID()
    pixIDtypeAsStr = image.GetPixelIDTypeAsString()
    
    if p2c:
        #print('\nImage info:')
        #print(f'type(image) = {type(image)}')
        print(f'    Image PixID = {pixID}')
        print(f'    Image PixIDTypeAsString = {pixIDtypeAsStr}')
        print(f'    Image size = {image.GetSize()}')
        print(f'    Image Max = {im_max(image)}')
        print(f'    Image Min = {im_min(image)}')
        print('\n    Conversion of image to pixarr:')
        print(f'    pixarr shape = {pixarr.shape}')
        if isinstance(uniqueVals, np.ndarray):
            if len(uniqueVals) < 7:
                print('\n    There are {len(uniqueVals)} unique values in pixarr:')
                print(f'    {uniqueVals}')
            else:
                
                print(f'\n    There are {len(uniqueVals)} unique values in pixarr:')
                print(f'    {uniqueVals[:3]}...{uniqueVals[-3:-1]}')
            
            """ Print histogram of uniqueVals
            https://www.machinelearningplus.com/plots/matplotlib-histogram-python-examples/
            """
            
            vals = pixarr.flatten()
            minVal = min(vals)
            maxVal = max(vals)
            
            freq, bins = np.histogram(vals, bins=10, range=[minVal, maxVal])
            
            print('\n      Distribution of values:')
            #for b, f in zip(bins[1:], freq):
            #    #print(round(b, 1), ' '.join(np.repeat('*', f)))
            #    print(f'      {round(b, 2)} - {f}')
            for i in reversed(range(len(freq))):
                valRange = f'[{round(bins[i], 2)} - {round(bins[i+1], 2)}]'
                print(f'      {valRange} : {freq[i]}')
            
            #vals = pixarr.flatten()
            #minVal = 0
            ##maxVal = 0.2
            #maxVal - 0.1
            #freq, bins = np.histogram(vals, bins=10, range=[minVal, maxVal])
            #print('\n      Distribution of values near 0:')
            #for i in reversed(range(len(freq))):
            #    valRange = f'[{round(bins[i], 2)} - {round(bins[i+1], 2)}]'
            #    print(f'      {valRange} : {freq[i]}')
            
        elif uniqueVals == None:
            print('       There are no uniqueVals (= {uniqueVals}) in pixarr:')
        print(f'       \nThere are {len(f2sInds)} frames with slice indices:')
        print(f'       {f2sInds}')
    
    return pixID, pixIDtypeAsStr, uniqueVals, f2sInds