# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 20:19:46 2021

@author: ctorti
"""



def ProportionOfPixArrInExtent_OLD(PixArr, FrameToSliceInds, DicomDir):
    """
    Determine what proportion of voxels that make up a 3D pixel array reside
    within the 3D volume of a 3D image.
    
    Inputs:
    ------
    
    PixArr : Numpy array with dimensions Z x Y x X where Z may be 1
        The pixel array that may represent a mask or labelmap.
    
    FrameToSliceInds : List of list of integers
        A list (for each frame) of the slice numbers that correspond to each 
        frame in PixArr.  This is equivalent to a list of Per-Frame Functional 
        Groups Sequence-to-slice inds (PFFGStoSliceInds). 
        
    DicomDir : string
        Directory containing the DICOMs that relate to PixArr.
        
        
    Outputs:
    -------
    
    p : float
        The fractional proportion (normalised to 1) of voxels in PixArr that 
        intersect the extent of RefImage.
    """
    
    from ConversionTools import PixArr2PtsByContour
    from ImageTools import ImportImage
    
    # Convert the pixel array to a list (for each frame in PixArr) of a list 
    # (for each point) of a list (for each dimension) of physical coordinates:
    PtsByFrame, CntDataByFrame = PixArr2PtsByContour(PixArr, 
                                                     FrameToSliceInds, 
                                                     DicomDir)
    
    Image = ImportImage(DicomDir)
    
    Vertices = GetExtremePoints(Image)
    
    IsInside = []
    
    for PtsInFrame in PtsByFrame:
        for pt in PtsInFrame:
            IsInside.append(IsPointInPolygon(pt, Vertices))
            
    return sum(IsInside)/len(IsInside)





def IsMatrixOrthogonal_OLD(Matrix):
    """
    Determine if a matrix is orthogonal.
    
    Inputs:
    ******
    
    Matrix : Numpy array of square shape
        A Numpy matrix of square shape.
        
        
    Outputs:
    *******
    
    IsOrthogonal : boolean
        True if Matrix is orthogonal, False otherwise.
        
        
    Notes:
    -----
    
    The following tests for orthogonality of a matrix M: 
        
        M is orthogonal if M*M.T = M.T*M, 
        i.e. if np.dot(M, M.T) - np.dot(M.T, M) = 0
    
    
    
    The following tests for orthonomality of a matrix M: 
        
        M is orthonormal if its determinant (np.linalg.det(M)) = +/- 1.
        
        M is orthonormal if M*M.T = I, where I is the identity matrix, 
        i.e. if np.dot(M, M.T) - np.identity(M.shape[0]) = 0
        
        M is orthonormal if M.T = M^-1, 
        i.e. if M.T - np.linalg.inv(M) = 0
        
        
        
    The following test for orthonormality/orthogonality doesn't provide the 
    expected result for either an orthogonal or orthonormal matrix:
        
        M is orthogonal if M.T x M^-1 = I, where I is the identity matrix,
        i.e. if np.cross(M.T, np.linalg.inv(M)) - np.identity(M.shape[0]) = 0
    """
    
    import numpy as np
    
    M = Matrix
    
    #shape = M.shape
    
    if len(M.shape) != 2:
        msg = f'The matrix has {len(M.shape)} dimensions. It must have 2.'
        
        raise Exception(msg)
        
    if M.shape[0] != M.shape[1]:
        msg = f'The matrix has shape {M.shape}. It must be a square matrix.'
        
        raise Exception(msg)
        
    
    Det = np.linalg.det(M)
    
    DotProdResult_1 =  np.dot(M, M.T) - np.identity(M.shape[0])
    
    DotProdResult_2 =  np.dot(M, M.T) - np.dot(M.T, M)
    
    DetInvResult = M.T - np.linalg.inv(M)
    
    CrossProdResult = np.cross(M.T, np.linalg.inv(M)) - np.identity(M.shape[0])
    
    
    IsDet1 = AreItemsEqualToWithinEpsilon(abs(Det), 1, epsilon=1e-06)
    
    SSEofDPR_1 = SumOfSquaresOfMatrix(DotProdResult_1)
    IsDPR0_1 = AreItemsEqualToWithinEpsilon(SSEofDPR_1, 0, epsilon=1e-06)
    
    SSEofDPR_2 = SumOfSquaresOfMatrix(DotProdResult_2)
    IsDPR0_2 = AreItemsEqualToWithinEpsilon(SSEofDPR_2, 0, epsilon=1e-06)
    
    SSEofDIR = SumOfSquaresOfMatrix(DetInvResult)
    IsDIR0 = AreItemsEqualToWithinEpsilon(SSEofDIR, 0, epsilon=1e-06)
    
    SSEofCPR = SumOfSquaresOfMatrix(CrossProdResult)
    IsCPR0 = AreItemsEqualToWithinEpsilon(SSEofCPR, 0, epsilon=1e-06)
    
    
    print(f'\nThe determinant of M = {Det} (must be +/- 1 if M is orthogonal).')
    print(f'\nM*M.T - I = {SSEofDPR_1} (must be 0 if M is orthogonal).')
    print(f'\nM*M.T - M.T*M = {SSEofDPR_2} (must be 0 if M is orthogonal).')
    print(f'\nM.T - M^-1 = {SSEofDIR} (must be 0 if M is orthogonal).')
    print(f'\nM.T x M^-1 - I = {SSEofCPR} (must be 0 if M is orthogonal).')
        
        
    return IsDPR0_1






