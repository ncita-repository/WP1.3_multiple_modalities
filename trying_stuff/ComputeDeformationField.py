# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:01:49 2020

@author: ctorti
"""

def ComputeDeformationField(MovingImage, TransformFilter):
    # Import packages:
    import SimpleITK as sitk
    import os

    # Get the current working directory:
    #CurrentWorkingDir = os.getcwd()

    # Initiate TransformixImageFilter:
    TransformixImFilt = sitk.TransformixImageFilter()
    #TransformixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    TransformixImFilt.LogToConsoleOff()
    TransformixImFilt.LogToFileOn() 
    #TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output
    
    # Set up for computing deformation field:
    TransformixImFilt.ComputeDeformationFieldOn() # <-- no output in Jupyter
    
    # Set up other computations:
    TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
    TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'
    
    # Set the transform parameter map:
    TransformixImFilt.SetTransformParameterMap(TransformFilter.GetTransformParameterMap())
    
    # Need to explicitely tell elastix that the image is 3D since by default it
    # will only transform to 2D:
    TransformixImFilt.SetMovingImage(MovingImage)
    
    # Transform the contour points:
    TransformixImFilt.Execute()
    #TransformixIm = TransformixImFilt.Execute()
    
    # Get the deformation field:
    DefField = TransformixImFilt.GetDeformationField()
    
    CWD = os.getcwd()
    
    # Check if deformation file was exported:
    if not os.path.isfile(os.path.join(CWD, 'deformationField.nii')):
        print('\nSeems that the deformation field was not exported.',
              '\nManually exporting file...')
        
        # Convert the deformation field to a numpy array:
        """ Note: 
            Executing this kills the kernel in Anaconda Jupyter. 
        """
        DefFieldNda = sitk.GetArrayFromImage(DefField) 
        
        # Write the deformation field to a .nii file:
        """ Note: 
            Executing this also kills the kernel in Anaconda Jupyter. 
        """
        DefFieldNiiFpath = os.path.join(CWD, 'deformationField.nii')
        
        sitk.WriteImage(DefField, DefFieldNiiFpath)
        
        print('\nDeformation field was exported to', DefFieldNiiFpath)
    
    return DefField