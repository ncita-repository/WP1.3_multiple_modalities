# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:00:33 2020

@author: ctorti
"""











def GetDeformationField():
    # Import packages:
    import SimpleITK as sitk
    #print(sitk.Version())
    import numpy as np
    import time
    import os
    #import copy
    #import importlib
    import json
    #import pydicom
    #import matplotlib.pyplot as plt
    #get_ipython().run_line_magic('matplotlib', 'inline')
    #get_ipython().run_line_magic('matplotlib', 'notebook')
    #from ipywidgets import interact, fixed
    #from IPython.display import clear_output
    
    #import gui
    #import registration_gui as rgui
    
    #from GetDicomFpaths import GetDicomFpaths
    #from GetDicoms import GetDicoms
    #from GetContourPoints3D import GetContourPoints3D
    #from GetAllContourPoints3D import GetAllContourPoints3D
    #from CreateInputFileForElastix import CreateInputFileForElastix
    #from ParseTransformixOutput import ParseTransformixOutput
    #import RegUtilityFuncs as ruf
    
    # Get the current working directory:
    CurrentWorkingDir = os.getcwd()
    
    # Import DataDict:
    ImportDir = r'C:\Temp\2020-05-13 ACRIN MR_4 to MR_12'
    JsonFpath = os.path.join(ImportDir, 'DataDict.json')
    
    file = open(JsonFpath)
    DataDict = json.load(file)
    file.close()
    

    
    # Define FixedDir and MovingDir based on the directories contained in 
    # DataDict:
    FixedDir = DataDict['Source']['DicomDir']
    MovingDir = DataDict['Target']['DicomDir']
    
    SwapFixedMoving = False
    #SwapFixedMoving = True
    
    #if SwapFixedMoving:
    #    print('FixedDir and MovingDir swapped.')
    #else:
    #    print('FixedDir and MovingDir not swapped.')
    
    # Read in the 3D stack of Fixed DICOMs:
    FixedReader = sitk.ImageSeriesReader()
    #FixedReader = sitk.ImageSeriesReader(sitk.sitkFloat32)
    if SwapFixedMoving:
        FixedNames = FixedReader.GetGDCMSeriesFileNames(MovingDir)
    else:
        FixedNames = FixedReader.GetGDCMSeriesFileNames(FixedDir)
    FixedReader.SetFileNames(FixedNames)
    FixedIm = FixedReader.Execute()
    FixedIm = sitk.Cast(FixedIm, sitk.sitkFloat32)
    
    # Read in the 3D stack of Moving DICOMs:
    MovingReader = sitk.ImageSeriesReader()
    #MovingReader = sitk.ImageSeriesReader(sitk.sitkFloat32)
    if SwapFixedMoving:
        MovingNames = MovingReader.GetGDCMSeriesFileNames(FixedDir)
    else:
        MovingNames = MovingReader.GetGDCMSeriesFileNames(MovingDir)
    MovingReader.SetFileNames(MovingNames)
    MovingIm = MovingReader.Execute()
    #MovingIm = MovingReader.Execute(sitk.sitkFloat32)
    MovingIm = sitk.Cast(MovingIm, sitk.sitkFloat32)
    
    
    # Register MovingIm to FixedIm using SimpleElastix:
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Initiate ElastixImageFilter:
    ElastixImFilt = sitk.ElastixImageFilter()
    ElastixImFilt.LogToConsoleOn() # <-- no output
    ElastixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel
    
    # Define the fixed and moving images:
    ElastixImFilt.SetFixedImage(FixedIm)
    ElastixImFilt.SetMovingImage(MovingIm)
    
    # Get the default parameter map template:
    #ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    ParamMap = sitk.GetDefaultParameterMap('affine')
    
    # Set registration parameters:
    ParamMap['AutomaticTransformInitialization'] = ['true']
    ParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ParamMap['WriteIterationInfo'] = ['true']
    ParamMap['MaximumNumberOfIterations'] = ['512']
    ParamMap['UseDirectionCosines'] = ['true']
    
    # Print the parameters:
    #for keys,values in ParamMap.items():
    #    print(keys, '=', values)
    
    # Set the parameter map:
    ElastixImFilt.SetParameterMap(ParamMap)
    
    # Register the 3D images:
    ElastixImFilt.Execute()
    
    # Get the registered image:
    RegIm = ElastixImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\nTook {Dtime} s to register the 3D image stacks.')
    
    
    
    # Compute the deformation field:
    
    """ 
    For backwards compatibility, only one of ComputeDeformationFieldOn() 
    or SetFixedPointSetFileName() can be active at any one time. 
    """
    
    
    
    # Get the transform parameter map:
    #TransformParameterMap = ElastixImFilt.GetParameterMap()
    
    # Initiate TransformixImageFilter:
    TransformixImFilt = sitk.TransformixImageFilter()
    TransformixImFilt.LogToConsoleOn() # <-- no output
    TransformixImFilt.LogToFileOn() 
    TransformixImFilt.SetOutputDirectory(CurrentWorkingDir) # <-- makes no difference, still no output
    # Set up for computing deformation field:
    TransformixImFilt.ComputeDeformationFieldOn() # <-- no output
    # Set up other computations to see if they produce output:
    TransformixImFilt.ComputeSpatialJacobianOn() # <-- produces 'spatialJacobian.nii'
    TransformixImFilt.ComputeDeterminantOfSpatialJacobianOn() # <-- produces 'fullSpatialJacobian.nii'
    # Set the transform parameter map:
    TransformixImFilt.SetTransformParameterMap(ElastixImFilt.GetTransformParameterMap())
    
    #TransformixImFilt.SetFixedPointSetFileName('inputpoints.txt')
    # Need to explicitely tell elastix that the image is 3D since by default it will
    # only transform to 2D:
    TransformixImFilt.SetMovingImage(MovingIm)
    
    # Transform the contour points:
    TransformixImFilt.Execute()
    #TransformixIm = TransformixImFilt.Execute()
    
    # Get the deformation field:
    DeformationField = TransformixImFilt.GetDeformationField()
    
    
    # Convert the deformation field to a numpy array:
    """ Note: 
        Executing this kills the kernel in Anaconda Jupyter. 
    """
    DeformationFieldNda = sitk.GetArrayFromImage(DeformationField) 
    
    #print('\nDeformationFieldNda:\n', DeformationFieldNda)
    
    
    # Export the numpy data array as a .csv:
    DefFieldCsvFpath = os.path.join(CurrentWorkingDir, 'defField.csv')
    print('\nExporting deformation field as CSV to', DefFieldCsvFpath)
    np.savetxt(DefFieldCsvFpath, DeformationFieldNda, delimiter=',')
    
    
    # Also write the image in .nii format:
    """ Note: 
        Executing this also kills the kernel in Anaconda Jupyter. 
    """
    DefFieldNiiFpath = os.path.join(CurrentWorkingDir, 'defField.nii')
    print('\nExporting deformation field as NII to', DefFieldNiiFpath)
    sitk.WriteImage(DeformationField, DefFieldNiiFpath)
    
    
    return DeformationField, DeformationFieldNda