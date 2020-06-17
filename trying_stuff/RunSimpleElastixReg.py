# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:00:14 2020

@author: ctorti
"""

def RunSimpleElastixReg(FixedIm, MovingIm):
    # Import packages:
    import SimpleITK as sitk
    import time
    
    # Start timing:
    times = []
    times.append(time.time())
    
    # Initiate ElastixImageFilter:
    ElastixImFilt = sitk.ElastixImageFilter()
    #ElastixImFilt.LogToConsoleOn() # <-- no output in Jupyter
    ElastixImFilt.LogToConsoleOff() # <-- no output in Jupyter
    ElastixImFilt.LogToFileOn() # <-- output's elastix.log ONLY ONCE per kernel in Jupyter
    
    # Define the fixed and moving images:
    ElastixImFilt.SetFixedImage(FixedIm)
    ElastixImFilt.SetMovingImage(MovingIm)
    
    # Get the default parameter map template for affine transformation:
    #ElastixImFilt.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    ElastixParamMap = sitk.GetDefaultParameterMap('affine')
    
    # Re-assign some parameters:
    ElastixParamMap['AutomaticTransformInitialization'] = ['true']
    ElastixParamMap['AutomaticTransformInitializationMethod'] = ['GeometricalCenter']
    ElastixParamMap['WriteIterationInfo'] = ['true']
    ElastixParamMap['MaximumNumberOfIterations'] = ['512']
    ElastixParamMap['UseDirectionCosines'] = ['true']
    """ 29/05: Trying this instead of trying to change it for Transformix """
    #ElastixParamMap['FinalBSplineInterpolationOrder'] = ['0'] 
    
    # Print the parameters:
    #for keys,values in ElastixParamMap.items():
    #    print(keys, '=', values)
        
    # Set the parameter map:
    ElastixImFilt.SetParameterMap(ElastixParamMap)
    
    print('Performing registration...')
    
    # Register the 3D images:
    ElastixImFilt.Execute()
    # Get the registered image:
    RegIm = ElastixImFilt.GetResultImage()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    print(f'\nTook {Dtime} s to register the 3D image stacks.')
    
    return RegIm, ElastixImFilt