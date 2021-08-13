# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:56:06 2021

@author: ctorti
"""


""" Functions that relate to matrices. """

import numpy as np

from general_tools.general import are_items_equal_to_within_eps


def matrix_sum_of_squares(Matrix):
    """
    Compute the root-sum-of-squares of a matrix.
    
    Parameters
    ----------
    Matrix : Numpy array
        A Numpy array (including matrices).
        
    Returns
    -------
    RSS : float
        The root-sum-of-squares of Matrix.
    """
    
    RSS = np.sqrt(np.sum(np.square(Matrix)))
    
    return RSS

def is_matrix_orthogonal(Matrix, LogToConsole=False):
    """
    Determine if a matrix is orthogonal.
    
    Parameters
    ----------
    Matrix : Numpy array of square shape
        A Numpy matrix of square shape.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
        
    Returns
    -------
    IsOrthogonal : bool
        True if Matrix is orthogonal, False otherwise.
    
    Note
    ----
    The following tests for orthogonality of a matrix M: 
        M is orthogonal if M*M.T = M.T*M, 
        i.e. if np.dot(M, M.T) - np.dot(M.T, M) = 0
    """
    
    M = Matrix
    
    if len(M.shape) != 2:
        msg = f'The matrix has {len(M.shape)} dimensions. It must have 2.'
        
        raise Exception(msg)
        
    if M.shape[0] != M.shape[1]:
        msg = f'The matrix has shape {M.shape}. It must be a square matrix.'
        
        raise Exception(msg)
        
    result =  np.dot(M, M.T) - np.dot(M.T, M)
    
    SSE = matrix_sum_of_squares(result)
    
    IsOrthogonal = are_items_equal_to_within_eps(SSE, 0, epsilon=1e-06)
    
    if LogToConsole:
        print(f'\nM*M.T - M.T*M = \n{result}')
        print(f'\nSSE = {SSE}')
        
        if IsOrthogonal:
            print('\nM is orthogonal')
        else:
            print('\nM is not orthogonal')
        
    return IsOrthogonal

def is_matrix_orthonormal(Matrix, LogToConsole=False):
    """
    Determine if a matrix is orthonormal.
    
    Parameters
    ----------
    Matrix : Numpy array of square shape
        A Numpy matrix of square shape.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    IsOrthonormal : bool
        True if Matrix is orthonormal, False otherwise.
        
    Note
    ----
    The following tests for orthonomality of a matrix M: 
        
        M is orthonormal if its determinant (np.linalg.det(M)) = +/- 1.
        
        M is orthonormal if M*M.T = I, where I is the identity matrix, 
        i.e. if np.dot(M, M.T) - np.identity(M.shape[0]) = 0
        
        M is orthonormal if M.T = M^-1, 
        i.e. if M.T - np.linalg.inv(M) = 0
    """
    
    M = Matrix
    
    #print(type(M))
    #print(M)
    
    if len(M.shape) != 2:
        msg = f'The matrix has {len(M.shape)} dimensions. It must have 2.'
        
        raise Exception(msg)
        
    if M.shape[0] != M.shape[1]:
        msg = f'The matrix has shape {M.shape}. It must be a square matrix.'
        
        raise Exception(msg)
        
    """ Either of the three tests below (Det = 1?, M.M_T = 0?, 
    SSE of M.T - M^-1 = 0?) would be suitable. Determinant test is chosen for
    no particular reason. """
    Det = np.linalg.det(M)
    IsOrthonormal = are_items_equal_to_within_eps(abs(Det), 1, epsilon=1e-06)
    
    MdotM_T =  np.dot(M, M.T) - np.identity(M.shape[0])
    SSEofMdotM_T = matrix_sum_of_squares(MdotM_T)
    #IsOrthonormal = AreItemsEqualToWithinEpsilon(SSEofMdotM_T, 0, 
    #                                             epsilon=1e-06)
    
    TequalsInv = M.T - np.linalg.inv(M)
    SSEofTequalsInv = matrix_sum_of_squares(TequalsInv)
    #IsOrthonormal = AreItemsEqualToWithinEpsilon(SSEofTequalsInv, 0, 
    #                                             epsilon=1e-06)
    
    if LogToConsole:
        print(f'\nThe determinant of M = {Det}')
        print(f'\nM*M.T - I = \n{MdotM_T}')
        print(f'\nSSE = {SSEofMdotM_T}')
        print(f'\nM.T - M^-1 = \n{SSEofTequalsInv}')
        print(f'\nSSE = {SSEofTequalsInv}')
        
        if IsOrthonormal:
            print('\nM is orthonormal')
        else:
            print('\nM is not orthonormal')
        
    return IsOrthonormal

def get_tx_matrix_type(TxMatrix, LogToConsole=False):
    """
    Determine the Transformation Matrix Type from a list of transformation
    parameters.
    
    Parameters
    ----------
    TxMatrix : list of strs
        The registration transform parameters.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    MatrixType : str
        A string containing the matrix type.  Possible outputs are:
            - 'RIGID'
            - 'RIGID_SCALE'
            - 'AFFINE'
    
    Note
    ----
    1. The matrix M is the 3x3 matrix within the 4x4 matrix TxParams read from
    TransformParameters.0.txt.  M only includes the elements that describe 
    rotations, scaling and shearing, i.e not translations (hence only the
    elements [0, 1, 2, 4, 5, 6, 8, 9, 10] in TxMatrix).
    
    TxMatrix has additional bottom row [0 0 0 1] (see 
    http://dicom.nema.org/medical/dicom/current/output/chtml/part17/chapter_P.html).
    
    Those additional elements are at indices 3, 7, 11, 15 in TxMatrix. Hence
    TxParams = [TxMatrix[i] for i in [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]].
    
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
    
    #print(f'TxParams = {TxParams}\n')
    
    #inds = [0, 1, 2, 4, 5, 6, 8, 9, 10] # 06/05/21
    inds = [0, 1, 2, 3, 4, 5, 6, 7, 8] # 04/06/21

    #M = np.array([TxParams[i] for i in inds]) # 06/05/21
    #M = np.array([float(TxParams[i]) for i in inds]) # 06/05/21
    #M = np.array([float(TxParams[i]) for i in range(9)]) # 07/05/21
    M = np.array([float(TxMatrix[i]) for i in inds]) # 04/06/21
    
    M = np.reshape(M, (3, 3))
    
    #print(f'M = {M}\n')
    
    IsOrthonormal = is_matrix_orthonormal(M, LogToConsole)
    
    if IsOrthonormal:
        MatrixType = 'RIGID'
    else:
        IsOrthogonal = is_matrix_orthogonal(M, LogToConsole)
        
        if IsOrthogonal:
            MatrixType = 'RIGID_SCALE'
        else:
            MatrixType = 'AFFINE'
    
    if LogToConsole:
        print('Matrix type is', MatrixType)
    
    return MatrixType