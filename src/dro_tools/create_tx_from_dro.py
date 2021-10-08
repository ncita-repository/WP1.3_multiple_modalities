# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:14:44 2021

@author: ctorti
"""


""" Functions that create SimpleITK Transforms from DROs. """

import SimpleITK as sitk
import numpy as np


def create_tx_from_spa_dro(dro, p2c=False):
    """
    Create a SimpleITK Transform based on parameters parsed from a DICOM
    Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    dro : Pydicom Object
        The DICOM Spatial Registration Object.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    sitkTx : SimpleITK Transform
        A SimpleITK Transform based on the parameters stored in dro.
    
    Note
    ----
    It will be assumed that RegistrationSequence has two items: the first for
    the Fixed scan and thev second for the Moving scan.  While the 
    FrameOfReferenceMatrix in the second item is necessary, the first provides
    a convenient means of arciving the FrameOfReferenceUID of the Fixed scan.
    
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    """
    
    matrixType = dro.RegistrationSequence[1]\
                    .MatrixRegistrationSequence[0]\
                    .MatrixSequence[0]\
                    .FrameOfReferenceTransformationMatrixType
    
    matrix = [float(s) for s in dro.RegistrationSequence[1]\
                                   .MatrixRegistrationSequence[0]\
                                   .MatrixSequence[0]\
                                   .FrameOfReferenceTransformationMatrix]
    
    if p2c:
        print(f'matrixType from DRO = {matrixType}')
        print(f'matrix from DRO = {matrix}')
    
    # TODO Finish 'RIGID' and 'RIGID_SCALE' cases 
    if matrixType == 'RIGID':
        sitkTx = sitk.Euler3DTransform()
        
        #sitkTx.SetParameters(txMatrix[0:6]) # 07/05/21
        #txParams = [txMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #sitkTx.SetParameters(txParams) # 07/05/21
    
    elif matrixType == 'RIGID_SCALE':
        sitkTx = sitk.Similarity3DTransform()
        
        #sitkTx.SetParameters(txMatrix[0:12]) # 07/05/21
        #txParams = [txMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #sitkTx.SetParameters(txParams) # 07/05/21
        
    elif matrixType == 'AFFINE':
        sitkTx = sitk.AffineTransform(3)
        
        inds = [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 3, 7, 11] # 01/09/21
        
        #sitkTx.SetParameters(txMatrix[0:12]) # 07/05/21
        #matrix = [matrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        matrix = [matrix[i] for i in inds] # 01/09/21
        
        sitkTx.SetParameters(matrix) # 07/05/21
        
    else:
        msg = f"matrixType = {matrixType} is not a recognised value. Only "\
              + "'RIGID', 'RIGID_SCALE' and 'AFFINE' are acceptable values."
        raise Exception(msg)
    
    if p2c:
        print(f'Parameters in new sitkTx = {sitkTx.GetParameters()}\n')
    
    return sitkTx

def create_pre_tx_from_def_dro(dro, p2c=False):
    """
    Create a SimpleITK Transform based on parameters parsed from
    PreDeformationMatrixRegistrationSequence in a DICOM Deformable Spatial 
    Registration Object (DRO).  
    
    Parameters
    ----------
    dro : Pydicom Object
        The DICOM Deformable Registration Object.
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    preDefTx : SimpleITK Transform
        A SimpleITK Transform based on the parameters stored in dro.
    
    Note
    ----
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    Even thought the FrameOfReferenceTransformationMatrixType may be 'RIGID',
    'RIGID_SCALE' or 'AFFINE', create an affine transform for simplicity.
    
    txMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in txMatrix. Hence
    txParams = [txMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    """
    
    txType = dro.DeformableRegistrationSequence[1]\
                .PreDeformationMatrixRegistrationSequence[0]\
                .FrameOfReferenceTransformationMatrixType
    
    txMatrix = dro.DeformableRegistrationSequence[1]\
                  .PreDeformationMatrixRegistrationSequence[0]\
                  .FrameOfReferenceTransformationMatrix
    
    txParams = [txMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]]
    
    if p2c:
        print(f'txType = {txType}\n')
        print(f'txMatrix = {txMatrix}\n')
        print(f'txParams = {txParams}\n')
    
    #if txType == 'RIGID':
    #    preDefTx = sitk.Euler3DTransform()
    #elif txType == 'RIGID_SCALE':
    #    preDefTx = sitk.Similarity3DTransform()
    #elif txType == 'AFFINE':
    #    preDefTx = sitk.AffineTransform(3)
    #else:
    #    msg = "FrameOfReferenceTransformationMatrixType is not one of the "\
    #          + f"expected values: 'RIGID', 'RIGID_SCALE' or 'AFFINE'."
    #    raise Exception(msg)
    
    preDefTx = sitk.AffineTransform(3)
    preDefTx.SetParameters(txParams)
    
    return preDefTx

def create_tx_from_def_dro(dro, refIm=None, p2c=False):
    """
    Create a SimpleITK BSpline transformation based on parameters parsed from
    a DICOM Deformable Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    dro : Pydicom Object
        The DICOM Deformable Registration Object.
    refIm : SimpleITK Image, optional (None by default)
        A reference image (e.g. fixed image) that shares the same direction 
        cosines as the BSpline grid.  If present, the direction of refIm will
        be used to set that of the grid.  If None, the direction cosines along
        x and y, as stored in 
        
        DeformableRegistrationSequence[1]\
        .DeformableRegistrationGridSequence[0]\
        .ImageOrientationPatient
        
        will be used to determine the direction cosine along the z direction.
    
    p2c : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    sitkTx : SimpleITK Transform
        A transform based on the parameters stored in dro.
    
    Note
    ----
    Sample DRO from:
    https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT#242838679219e971f0494026a216c74aeae636e6
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    In order to set the direction of the coefficient image below, the direction
    cosine along the z-direction is obtained by taking the cross product of the
    direction cosines along the x and y directions.  This is the simplest 
    solution but due to the 16 character limit when storing ImagePositionPatient
    and ImageOrientationPatient there are small differences from the values
    that would be obtained from FixIm.GetDirection().
    
    e.g. the difference between the direction cosines along x, y and z obtained
    from IOP and FixIm.GetDirection() for a sample data was:
        
    [-5.83795068e-09,  2.60326283e-08, -2.79367864e-08, 
    -5.55216258e-09, -1.82105636e-08,  2.78585151e-08,  
    2.62768308e-08, -2.83423783e-08, -2.48949777e-08]
    """
    
    gridOrig = [float(item) for item in dro.DeformableRegistrationSequence[1]\
                                           .DeformableRegistrationGridSequence[0]\
                                           .ImagePositionPatient]
    
    if p2c:
        print(f'type(gridOrig) = {type(gridOrig)}')
        print(f'gridOrig = {gridOrig}\n')
    
    if refIm == None:
        gridDir = [float(item) for item in dro.DeformableRegistrationSequence[1]\
                                              .DeformableRegistrationGridSequence[0]\
                                              .ImageOrientationPatient]
        
        # Add the direction cosine along the z-direction:
        gridDir.extend(list(np.cross(np.array(gridDir[0:3]), np.array(gridDir[3:]))))
    else:
        gridDir = refIm.GetDirection()
    
    if p2c:
        print(f'type(gridDir) = {type(gridDir)}')
        print(f'gridDir = {gridDir}\n')
    
    gridSize = [item for item in dro.DeformableRegistrationSequence[1]\
                                    .DeformableRegistrationGridSequence[0]\
                                    .GridDimensions]
    
    if p2c:
        print(f'type(gridSize) = {type(gridSize)}')
        print(f'gridSize = {gridSize}\n')
    
    gridSpacing = [item for item in dro.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .GridResolution]
    
    if p2c:
        print(f'type(gridSpacing) = {type(gridSpacing)}')
        print(f'gridSpacing = {gridSpacing}\n')
    
    #dim = len(gridSize)
    #
    #meshSize = [item - dim for item in gridSize]
    
    #vectGridData = [item for item in dro.DeformableRegistrationSequence[1]\
    #                                    .DeformableRegistrationGridSequence[0]\
    #                                    .VectorGridData]
    
    #print(f'type(vectGridData) = {type(vectGridData)}')
    #print(f'vectGridData = {vectGridData}\n')
    
    # Convert from bytes to list of floats:
    flatList = list(np.frombuffer(dro.DeformableRegistrationSequence[1]\
                                     .DeformableRegistrationGridSequence[0]\
                                     .VectorGridData, dtype=np.float64))
    
    # Unflatten list to (x, y, z) tuples:
    listOfTuples = zip(*[iter(flatList)]*3)
    
    # Unzip from x, y, z to lists for x, y and z:
    coeffX, coeffY, coeffZ = zip(*listOfTuples)
    
    coeffArrX = np.reshape(np.array(coeffX), gridSize)
    coeffArrY = np.reshape(np.array(coeffY), gridSize)
    coeffArrZ = np.reshape(np.array(coeffZ), gridSize)
    
    if p2c:
        print(f'coeffArrX.shape = {coeffArrX.shape}')
        print(f'coeffArrY.shape = {coeffArrY.shape}')
        print(f'coeffArrZ.shape = {coeffArrZ.shape}\n')
    
    coeffImX = sitk.GetImageFromArray(coeffArrX)
    coeffImX.SetOrigin(gridOrig)
    coeffImX.SetDirection(gridDir)
    coeffImX.SetSpacing(gridSpacing)
    
    coeffImY = sitk.GetImageFromArray(coeffArrY)
    coeffImY.SetOrigin(gridOrig)
    coeffImY.SetDirection(gridDir)
    coeffImY.SetSpacing(gridSpacing)
    
    coeffImZ = sitk.GetImageFromArray(coeffArrZ)
    coeffImZ.SetOrigin(gridOrig)
    coeffImZ.SetDirection(gridDir)
    coeffImZ.SetSpacing(gridSpacing)
    
    sitkTx = sitk.BSplineTransform([coeffImX, coeffImY, coeffImZ])
    
    return sitkTx

