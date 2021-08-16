# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 13:33:52 2021

@author: ctorti
"""


""" old resampling functions.

see image_tools.resampling.py for current functions
"""


def resample_im_v1(Im, RefIm, Tx='Identity', Interp='Linear', LogToConsole=False):
    """
    08/07/21: DELETE THIS
    
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------   
    Im : SimpleITK Image
        The 3D image to be resampled.
    RefIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
    Tx : str, optional ('Identity' by default)
        The transform to be used.  Acceptable inputs are:
        - 'Identity'
        - 'Affine'
    Interp : str, optional ('Linear' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    ResIm : SimpleITK image
        The resampled 3D image.
    
    Note
    ----
    While a linear (or BSpline) interpolator is appropriate for intensity 
    images, only a NearestNeighbor interpolator is appropriate for binary 
    images (e.g. segmentations) so that no new labels are introduced.
    """
    
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

def resample_im_v2(Im, RefIm, SitkTx, Interp='Linear'):
    """
    08/07/21: DELETE THIS
    
    Resample a 3D SimpleITK image.
    
    Parameters
    ----------   
    Im : SimpleITK Image
        The 3D image to be resampled.             
    RefIm : SimpleITK Image
        The 3D image reference image whose gridspace Im will be resampled to.
    SitkTx : SimpleITK Transform
        The SimpleIK Transform to be used.
    Interp : str, optional ('linear' by default)
        The type of interpolation to be used.  Acceptable inputs are:
        - 'Linear'
        - 'Bspline'
        - 'NearestNeighbor'
        - 'LabelGaussian'
    LogToConsole : bool, optional (False by default)
        Denotes whether some results will be logged to the console.
        
    Returns
    -------
    ResIm : SimpleITK image
        The resampled 3D image.
    
    Note
    ----
    While a linear (or BSpline) interpolator is appropriate for intensity 
    images, only a NearestNeighbor interpolator is appropriate for binary 
    images (e.g. segmentations) so that no new labels are introduced.
    """
    
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