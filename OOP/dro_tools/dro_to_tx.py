# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:14:44 2021

@author: ctorti
"""


""" Functions that create SimpleITK Transforms from DROs. """

import SimpleITK as sitk
import numpy as np


def create_tx_from_spa_dro(Dro, LogToConsole=False):
    """
    09/08/21:
        I probably should also return the pre-reg transform if it's different
        from the identity transform.
    
    Create a SimpleITK Transform based on parameters parsed from a DICOM
    Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    Dro : Pydicom Object
        The DICOM Spatial Registration Object.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    SitkTx : SimpleITK Transform
        A SimpleITK Transform based on the parameters stored in Dro.
    
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
    
    MatrixType = Dro.RegistrationSequence[1]\
                    .MatrixRegistrationSequence[0]\
                    .MatrixSequence[0]\
                    .FrameOfReferenceTransformationMatrixType
    
    Matrix = [float(s) for s in Dro.RegistrationSequence[1]\
                                   .MatrixRegistrationSequence[0]\
                                   .MatrixSequence[0]\
                                   .FrameOfReferenceTransformationMatrix]
    
    if LogToConsole:
        print(f'MatrixType from DRO = {MatrixType}')
        print(f'Matrix from DRO = {Matrix}')
    
    if MatrixType == 'RIGID':
        SitkTx = sitk.Euler3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:6]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
    
    elif MatrixType == 'RIGID_SCALE':
        SitkTx = sitk.Similarity3DTransform()
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        #TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        #SitkTx.SetParameters(TxParams) # 07/05/21
        
    elif MatrixType == 'AFFINE':
        SitkTx = sitk.AffineTransform(3)
        
        #SitkTx.SetParameters(TxMatrix[0:12]) # 07/05/21
        Matrix = [Matrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]] # 07/05/21
        SitkTx.SetParameters(Matrix) # 07/05/21
        
    else:
        msg = f"MatrixType = {MatrixType} is not a recognised value. Only "\
              + "'RIGID', 'RIGID_SCALE' and 'AFFINE' are acceptable values."
        raise Exception(msg)
    
    if LogToConsole:
        print(f'Parameters in new SitkTx = {SitkTx.GetParameters()}\n')
    
    return SitkTx

def create_pre_tx_from_def_dro(Dro, LogToConsole=False):
    """
    Create a SimpleITK Transform based on parameters parsed from
    PreDeformationMatrixRegistrationSequence in a DICOM Deformable Spatial 
    Registration Object (DRO).  
    
    Parameters
    ----------
    Dro : Pydicom Object
        The DICOM Deformable Registration Object.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    PreDefTx : SimpleITK Transform
        A SimpleITK Transform based on the parameters stored in Dro.
    
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
    
    TxMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. Hence
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    """
    
    TxType = Dro.DeformableRegistrationSequence[1]\
                .PreDeformationMatrixRegistrationSequence[0]\
                .FrameOfReferenceTransformationMatrixType
    
    TxMatrix = Dro.DeformableRegistrationSequence[1]\
                  .PreDeformationMatrixRegistrationSequence[0]\
                  .FrameOfReferenceTransformationMatrix
    
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]]
    
    if LogToConsole:
        print(f'TxType = {TxType}\n')
        print(f'TxMatrix = {TxMatrix}\n')
        print(f'TxParams = {TxParams}\n')
    
    #if TxType == 'RIGID':
    #    PreDefTx = sitk.Euler3DTransform()
    #elif TxType == 'RIGID_SCALE':
    #    PreDefTx = sitk.Similarity3DTransform()
    #elif TxType == 'AFFINE':
    #    PreDefTx = sitk.AffineTransform(3)
    #else:
    #    msg = "FrameOfReferenceTransformationMatrixType is not one of the "\
    #          + f"expected values: 'RIGID', 'RIGID_SCALE' or 'AFFINE'."
    #    raise Exception(msg)
    
    PreDefTx = sitk.AffineTransform(3)
    PreDefTx.SetParameters(TxParams)
    
    return PreDefTx

def create_tx_from_def_dro(Dro, RefIm=None, LogToConsole=False):
    """
    Create a SimpleITK BSpline transformation based on parameters parsed from
    a DICOM Deformable Spatial Registration Object (DRO).  
    
    Parameters
    ----------
    Dro : Pydicom Object
        The DICOM Deformable Registration Object.
    RefIm : SimpleITK Image, optional (None by default)
        A reference image (e.g. fixed image) that shares the same direction 
        cosines as the BSpline grid.  If present, the direction of RefIm will
        be used to set that of the grid.  If None, the direction cosines along
        x and y, as stored in 
        
        DeformableRegistrationSequence[1]\
        .DeformableRegistrationGridSequence[0]\
        .ImageOrientationPatient
        
        will be used to determine the direction cosine along the z direction.
    
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    SitkTx : SimpleITK Transform
        A transform based on the parameters stored in Dro.
    
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
    
    GridOrig = [float(item) for item in Dro.DeformableRegistrationSequence[1]\
                                           .DeformableRegistrationGridSequence[0]\
                                           .ImagePositionPatient]
    
    if LogToConsole:
        print(f'type(GridOrig) = {type(GridOrig)}')
        print(f'GridOrig = {GridOrig}\n')
    
    if RefIm == None:
        GridDir = [float(item) for item in Dro.DeformableRegistrationSequence[1]\
                                              .DeformableRegistrationGridSequence[0]\
                                              .ImageOrientationPatient]
        
        # Add the direction cosine along the z-direction:
        GridDir.extend(list(np.cross(np.array(GridDir[0:3]), np.array(GridDir[3:]))))
    else:
        GridDir = RefIm.GetDirection()
    
    if LogToConsole:
        print(f'type(GridDir) = {type(GridDir)}')
        print(f'GridDir = {GridDir}\n')
    
    GridSize = [item for item in Dro.DeformableRegistrationSequence[1]\
                                    .DeformableRegistrationGridSequence[0]\
                                    .GridDimensions]
    
    if LogToConsole:
        print(f'type(GridSize) = {type(GridSize)}')
        print(f'GridSize = {GridSize}\n')
    
    GridSpacing = [item for item in Dro.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .GridResolution]
    
    if LogToConsole:
        print(f'type(GridSpacing) = {type(GridSpacing)}')
        print(f'GridSpacing = {GridSpacing}\n')
    
    #dim = len(GridSize)
    #
    #MeshSize = [item - dim for item in GridSize]
    
    #VectGridData = [item for item in Dro.DeformableRegistrationSequence[1]\
    #                                    .DeformableRegistrationGridSequence[0]\
    #                                    .VectorGridData]
    
    #print(f'type(VectGridData) = {type(VectGridData)}')
    #print(f'VectGridData = {VectGridData}\n')
    
    # Convert from bytes to list of floats:
    FlatList = list(np.frombuffer(Dro.DeformableRegistrationSequence[1]\
                                     .DeformableRegistrationGridSequence[0]\
                                     .VectorGridData, dtype=np.float64))
    
    # Unflatten list to (x, y, z) tuples:
    ListOfTuples = zip(*[iter(FlatList)]*3)
    
    # Unzip from x, y, z to lists for x, y and z:
    CoeffX, CoeffY, CoeffZ = zip(*ListOfTuples)
    
    CoeffArrX = np.reshape(np.array(CoeffX), GridSize)
    CoeffArrY = np.reshape(np.array(CoeffY), GridSize)
    CoeffArrZ = np.reshape(np.array(CoeffZ), GridSize)
    
    if LogToConsole:
        print(f'CoeffArrX.shape = {CoeffArrX.shape}')
        print(f'CoeffArrY.shape = {CoeffArrY.shape}')
        print(f'CoeffArrZ.shape = {CoeffArrZ.shape}\n')
    
    CoeffImX = sitk.GetImageFromArray(CoeffArrX)
    CoeffImX.SetOrigin(GridOrig)
    CoeffImX.SetDirection(GridDir)
    CoeffImX.SetSpacing(GridSpacing)
    
    CoeffImY = sitk.GetImageFromArray(CoeffArrY)
    CoeffImY.SetOrigin(GridOrig)
    CoeffImY.SetDirection(GridDir)
    CoeffImY.SetSpacing(GridSpacing)
    
    CoeffImZ = sitk.GetImageFromArray(CoeffArrZ)
    CoeffImZ.SetOrigin(GridOrig)
    CoeffImZ.SetDirection(GridDir)
    CoeffImZ.SetSpacing(GridSpacing)
    
    SitkTx = sitk.BSplineTransform([CoeffImX, CoeffImY, CoeffImZ])
    
    return SitkTx

