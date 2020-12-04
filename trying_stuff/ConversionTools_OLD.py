# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:29:58 2020

@author: ctorti
"""





def Inds2SparseLabmap_OLD(Indices, RefImage): # not useful since Sparse2Labelmap doesn't work
    """
    Note:
        
        The problem with this function is that by using the SimpleITK image 
        format any slice that contains more than one contour (e.g. belonging to 
        more than ROI) will exist on the same slice.  What I want is a Numpy 
        array with a different mask for each contour.  
        
        The other reason for using the SimpleITK
        image format was because I had intended to use sitk filters to convert 
        the sparse mask into a filled one (e.g. using 
        BinaryClosingByReconstructionImageFilter), but since that failed I have
        no reason to use the sitk image format.
        
    
    Convert a list of indices to a sparse 3D Numpy data array.  
    
    Inputs:
        Indices  - (List of a list of integers) List of list of point indices
        
        RefImage - (SimpleITK image) Reference image whose size will be used to
                   define the shape of the sparse data array
        
        
    Returns:
        #Matrix   - (Numpy array, dtype int) The point indices as a Numpy array
        Sparse   - (SimpleITK image) The point indices as a SimpleITK image
    """
    
    import numpy as np
    import SimpleITK as sitk
    
    ImSize = RefImage.GetSize()
    
    # Initialise the sparse matrix:
    """
    Note: .GetSize() provides the size as (Ncols, Nrows, Nslices) but for
    Numpy arrays the order will be (Nslices, Nrows, Ncols).
    """
    Matrix = np.zeros((ImSize[2], ImSize[1], ImSize[0]), dtype=int)
    
    #print(f'Matrix.shape = {Matrix.shape}')
    
    # Loop through all contours:
    for contour in Indices:
        # Ignore any contours that have less than 3 (the minimum number to
        # define a closed contour) indices (as may arise during manual 
        # contouring):
        if len(contour) > 2:
            # Each i,j,k index along x,y,z:
            for i,j,k in contour:
                #print(f'i = {i}, j = {j}, k = {k}')
                
                Matrix[k,j,i] = 1
                
    
    #return Matrix
    return sitk.GetImageFromArray(Matrix)
            
    






def Sparse2Labelmap_OLD(Sparse): # doesn't work
    """
    This doesn't work.
    
    
    Convert a sparse labelmap to a filled labelmap, i.e. close the contour
    indices to result in masks.  
    
    Inputs:
        #Matrix   - (Numpy array, dtype int) The point indices as a Numpy array
        Sparse   - (SimpleITK image) Sparse labelmap (containing only indices
                   that define the contours)
        
        
    Returns:
        Labmap   - (SimpleITK image) Labelmap from closing of contour indices
    """
    
    import SimpleITK as sitk
    
    # First convert from "labelmap" to binary:
    #Labmap2BinaryFilt = sitk.LabelMapToBinaryImageFilter()
    
    #Labmap = Labmap2BinaryFilt.Execute(Sparse)
    
    # Close labelmap:
    ClosingFilt = sitk.BinaryClosingByReconstructionImageFilter()
    ClosingFilt.FullyConnectedOn()
    
    
    #ClosingFilt = sitk.BinaryMorphologicalClosingImageFilter()
    
    
    
    #Labmap = ClosingFilt.Execute(Matrix)
    #Labmap = ClosingFilt.Execute(Labmap)
    Labmap = ClosingFilt.Execute(Sparse)
    
    
    return Labmap









