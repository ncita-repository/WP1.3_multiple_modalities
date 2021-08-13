# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 09:18:36 2021

@author: ctorti
"""


""" Registering functions for SimpleITK images. """

import SimpleITK as sitk
import numpy as np
import time
import winsound
import matplotlib.pyplot as plt
from IPython.display import clear_output
from ipywidgets import interact, fixed

from general_tools.general import check_file_ext
from general_tools.fiducials import get_landmark_tx
import image_tools.registration_utilities as ru
import image_tools.registration_callbacks as rc


def command_iteration(method):
    print(f"{method.GetOptimizerIteration():3} = {method.GetMetricValue():10.5f}")
    print("\t#: ", len(method.GetOptimizerPosition()))

def command_multi_iteration(method):
    print("--------- Resolution Changing ---------")
    
def start_plot():
    """ Callback invoked when the StartEvent happens, sets up our new data. """
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []

def end_plot():
    """ Callback invoked when the EndEvent happens, do cleanup of data and 
    figure. """
    
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()

def plot_values(RegMethod):
    """ Callback invoked when the IterationEvent happens, update our data 
    and display new figure. """
    
    global metric_values, multires_iterations
    
    metric_values.append(RegMethod.GetMetricValue())                                       
    # Clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    # Plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
def update_multires_iterations():
    """ Callback invoked when the sitkMultiResolutionIterationEvent happens, 
    update the index into the metric_values list. """
    
    global metric_values, multires_iterations
    
    multires_iterations.append(len(metric_values))

def plot_fix_mov_ims(FixInd, MovInd, FixPixArr, MovPixArr):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(FixPixArr[FixInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Fixed image')
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(MovPixArr[MovInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Moving image')
    plt.axis('off')
    
    plt.show()

def plot_blended_im(Ind, alpha, FixIm, ResIm):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. """

    BlendedIm = (1.0 - alpha)*FixIm[:,:,Ind] + alpha*ResIm[:,:,Ind] 
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.axis('off')
    plt.show()

def non_def_reg_ims(FixIm, MovIm, RegTxName='affine', InitMethod='landmarks', 
                    FixFidsFpath='', MovFidsFpath='', SamplingPercentage=50, 
                    NumIters=500, LearningRate=1.0, LogToConsole=False):
    """  
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Parameters
    ----------
    FixIm : SimpleITK Image 
        The 3D image that MovIm will be registered to
    MovIm : SimpleITK Image
        The 3D image that will be registered to FixIm.
    RegTxName : str, optional ('affine' by default)
        Transformation to use for registration.  Acceptable values include:
        - 'rigid'
        - 'affine'
    InitMethod : str, optional ('landmarks' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
        If InitMethod = 'landmarks' but file paths to txt files containing
        fiducials are not provided the next default value will be 'moments'.
    FixFidsFpath : str, optional ('' by default)
        The file path of the text file containing fiducials for FixIm. 
        This parameter is only relevant if InitMethod = 'landmarks'.
    MovFidsFpath : str, optional ('' by default)
        The file path of the text file containing fiducials for MovIm. 
        This argument is only relevant if InitMethod = 'landmarks'.
    SamplingPercentage : int or float, optional (50 by default)
        The percentage of the images that are sampled for evaluating the metric.
    NumIters : int, optional (500 by default)
        The maximum number of iterations used when optimising.
    LearningRate : float, optional (1.0 by default)
        The learning rate used for the gradient descent optimiser.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
        
    Returns
    -------
    RegIm : SimpleITK Image
        The 3D registered image.
    InitialTx : SimpleITK Transform
        The SimpleITK Transform used to initialise the registration.
    FinalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform MovIm to FixIm.
    
    Note
    ----
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Since initial alignment is geometrical (use of fiducials), initialising the
    registration with a geometrical centering is suitable (moments is not).
    
    It seems that a higher number of samples (e.g. 50%) is required to get
    results that are not so dependent upon the random sampling.
    """
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple text only
    PlotToConsole = True # plots
    
    SamplingPercentage = SamplingPercentage/100
    
    if RegTxName == 'rigid':
        SitkTx = sitk.Euler3DTransform()
    elif RegTxName == 'affine':
        SitkTx = sitk.AffineTransform(FixIm.GetDimension())
    else:
        msg = f"RegTxName must be either 'rigid' or 'affine'. {RegTxName} is "\
              + "not an accepted value."
        raise Exception(msg)
    
    # Start timing:
    times = []
    times.append(time.time())
    
    if InitMethod=='landmarks':
        if (check_file_ext(FixFidsFpath, '.txt') and
            check_file_ext(MovFidsFpath, '.txt')):
            print('File paths to fiducials were provided. Attempting to',
                  'create a landmark transform..\n')
            
            # Get landmark transform:
            InitialTx, FixPts, MovPts\
                = get_landmark_tx(Transform=RegTxName, 
                                  FixFidsFpath=FixFidsFpath, 
                                  MovFidsFpath=MovFidsFpath, 
                                  FixIm=FixIm, MovIm=MovIm, FlipK=False,
                                  LogToConsole=True)
        else:
            print('File paths to fiducials were not provided/found. Initialising',
                  'using CenteredTransformInitializerFilter.MOMENTS..\n')
        
            CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
            
            InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                          SitkTx, CTIF)
    else:
        if InitMethod in ['geometry', 'geometricalcenter']:
            CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
        elif InitMethod in ['moments', 'centerofgravity']:
            CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
        
        InitialTx = sitk.CenteredTransformInitializer(FixIm, MovIm,
                                                      SitkTx, CTIF)
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    RegMethod.SetInitialTransform(InitialTx, inPlace=False)
    
    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #RegMethod.SetMetricAsJointHistogramMutualInformation()
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    #RegMethod.SetMetricSamplingStrategy(RegMethod.REGULAR)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Setup for the multi-resolution framework. """
    ShrinkFactors = [4,2,1]
    SmoothingSigmas = [2,1,1]
    #SmoothingSigmas = [2,1,0]
    
    RegMethod.SetShrinkFactorsPerLevel(shrinkFactors=ShrinkFactors)
    RegMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=SmoothingSigmas)
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    
    """ Interpolator settings. """
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    #RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
    #                                        numberOfIterations=500, 
    #                                        convergenceMinimumValue=1e-6, 
    #                                        convergenceWindowSize=10)
    RegMethod.SetOptimizerAsGradientDescent(learningRate=LearningRate, 
                                            numberOfIterations=NumIters, 
                                            estimateLearningRate=RegMethod.Once)
    #RegMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=LearningRate, 
    #                                                  numberOfIterations=100,
    #                                                  convergenceMinimumValue=1e-4,
    #                                                  #convergenceMinimumValue=1e-6,
    #                                                  convergenceWindowSize=5)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    if LogToConsole:
        if PlotToConsole:
            # Connect all of the observers for plotting during registration:
            RegMethod.AddCommand(sitk.sitkStartEvent, start_plot)
            RegMethod.AddCommand(sitk.sitkEndEvent, end_plot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                                 update_multires_iterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: plot_values(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
    # Get the number of valid points:
    NumOfValidPts = RegMethod.GetMetricNumberOfValidPoints()
    
    # Get the number of points in FixIm:
    NumOfPtsInFixIm = np.prod(FixIm.GetSize())
    
    # The % of valid points:
    PercentageValidPts = 100*NumOfValidPts/NumOfPtsInFixIm
    
    print(f'NumOfValidPts = {NumOfValidPts} = {PercentageValidPts}%\n')
    
    #FinalTx = RegMethod.Execute(sitk.Cast(FixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(MovIm, sitk.sitkFloat32))
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    MetricVal = RegMethod.GetMetricValue()
    StopConditionDesc = RegMethod.GetOptimizerStopConditionDescription()
    #FinalIterNum = RegMethod.GetOptimizerIteration()
    print(f'Final metric value: {MetricVal}')
    print(f"Optimizer stopping condition: {StopConditionDesc}")
    #print(f"Final Iteration: {FinalIterNum}")
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                          MovIm.GetPixelID())
    
    if False:
        interact(plot_blended_im, Ind=(0,FixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), FixIm=(FixIm), ResIm=fixed(RegIm));
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s to perform image registration.\n')
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    return RegIm, InitialTx, FinalTx

def bspline_reg_ims(FixIm, MovIm, FixFidsFpath='', MovFidsFpath='',
                    NumControlPts=8, SamplingPercentage=5, NumIters=100, 
                    LearningRate=5.0, LogToConsole=False):
    """ 
    Register two 3D SimpleITK images using a B-spline transformation in
    SimpleITK.
    
    Parameters
    ----------
    FixIm : SimpleITK Image 
        The 3D image that MovIm will be registered to.
    MovIm : SimpleITK Image
        The 3D image that will be registered to FixIm.
    FixFidsFpath : str, optional ('' by default)
        The file path of the txt file containing fiducials for FixIm. If '' the
        BSpline will be initialised using BSplineTransformInitializer.
    MovFidsFpath : str, optional ('' by default)
        The file path of the txt file containing fiducials for MovIm. If '' the
        BSpline will be initialised using BSplineTransformInitializer.
    NumControlPts : int, optional (8 by default)
        The number of control points that defines the BSpline grid.
    SamplingPercentage : int or float, optional (5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    NumIters : int, optional (100 by default)
        The maximum number of iterations used during optimisation.
    LearningRate : float, optional (5.0 by default)
        The learning rate used for the gradient descent optimiser.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
        
    Returns
    -------
    RegIm : SimpleITK Image
        The 3D registered image.
    RegMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
    InitialTx : SimpleITK Transform
        The Bspline SimpleITK Transform used to initialise the BSpline grid
        prior to registration.
    FinalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform MovIm to FixIm.
    
    Note
    ----
    Zivy advised not to use gradient descent line search for bspline:
    https://github.com/SimpleITK/SimpleITK/issues/1415
    """
    
    # Print simple metric values to console or plot of metric value v iteration?
    PlotToConsole = False # simple text only
    PlotToConsole = True # plots
    
    # Use scale factors for mesh grid?
    # Note 08/07/21: Not sure it works properly using scale factors.
    UseScaleFactors = False
    #UseScaleFactors = True
    
    """ Scale factors: 
    
    Maybe the scale factors should take into account the NumOfControlPts
    that were used in LandmarkTx and the desired NumOfControlPts for the 
    registration?
    """
    ScaleFactors = [1, 2, 4]
    
    dim = FixIm.GetDimension()
    
    if UseScaleFactors:
        MeshSize = [NumControlPts // ScaleFactors[-1]] * dim
    else:
        MeshSize = [NumControlPts] * dim
    
    SamplingPercentage = SamplingPercentage/100
    
    # Start timing:
    times = []
    times.append(time.time())
    
    if (check_file_ext(FixFidsFpath, '.txt') and
            check_file_ext(MovFidsFpath, '.txt')):
        print('File paths to fiducials were provided. Attempting to create a',
              'landmark transform..\n')
        
        # Get landmark transform:
        InitialTx, FixPts,\
            MovPts = get_landmark_tx(Transform='bspline', 
                                     FixFidsFpath=FixFidsFpath, 
                                     MovFidsFpath=MovFidsFpath, 
                                     FixIm=FixIm, MovIm=MovIm, FlipK=False,
                                     Buffer=2, LogToConsole=True)
    else:
        print('File paths to fiducials were not provided/found. Initialising',
              'the Bspline using BSplineTransformInitializer..\n')
        
        InitialTx\
            = sitk.BSplineTransformInitializer(image1=FixIm,
                                               transformDomainMeshSize=MeshSize)
    
    print(f'InitialTx is {InitialTx.GetName()}\n')
    
    """ Registration. """
    RegMethod = sitk.ImageRegistrationMethod()
    
    if UseScaleFactors:
        RegMethod.SetInitialTransformAsBSpline(InitialTx, inPlace=True,
                                               scaleFactors=ScaleFactors)
    else:
        RegMethod.SetInitialTransform(InitialTx, inPlace=True)
    
    """ Similarity metric settings. """
    RegMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    RegMethod.SetMetricSamplingStrategy(RegMethod.RANDOM)
    RegMethod.SetMetricSamplingPercentage(SamplingPercentage)
    
    """ Optimiser settings. Use LBFGSB (without scale factors) or LBFGS2 (with 
    scale factors).
    """
    if UseScaleFactors:
        """ Use the LBFGS2 instead of LBFGS. The latter cannot adapt to the 
        changing control grid resolution. """
        RegMethod.SetOptimizerAsLBFGS2(solutionAccuracy=1e-2, 
                                       numberOfIterations=NumIters, 
                                       deltaConvergenceTolerance=0.01)
    else:
        RegMethod.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-5, 
                                       numberOfIterations=NumIters)
    RegMethod.SetOptimizerScalesFromPhysicalShift()
    
    RegMethod.SetInterpolator(sitk.sitkLinear)
    
    #RegMethod.SetShrinkFactorsPerLevel([6, 2, 1])
    #RegMethod.SetSmoothingSigmasPerLevel([6, 2, 1])
    RegMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    RegMethod.SetSmoothingSigmasPerLevel([4, 2, 1]) # <-- better result than [2, 1, 0]?..
    #RegMethod.SetSmoothingSigmasPerLevel([2, 1, 0]) # in 65_Registration_FFD.ipynb example
    RegMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn() # <-- THIS WAS MISSING until 21:51 on 1/7/21
    
    if LogToConsole:
        if PlotToConsole:
            # Connect all of the observers so that we can perform plotting  
            # during registration.
            RegMethod.AddCommand(sitk.sitkStartEvent, start_plot)
            RegMethod.AddCommand(sitk.sitkEndEvent, end_plot)
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                                 update_multires_iterations) 
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: plot_values(RegMethod))
        else:
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: command_iteration(RegMethod))
            RegMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                                 lambda: command_multi_iteration(RegMethod))
    
        # Since corresponding points in the fixed and moving image are available  
        # display the similarity metric and the TRE during the registration:
        """ TRY THIS - COULD BE USEFUL: """
        if False:
            RegMethod.AddCommand(sitk.sitkStartEvent, rc.metric_and_reference_start_plot)
            RegMethod.AddCommand(sitk.sitkEndEvent, rc.metric_and_reference_end_plot)
            RegMethod.AddCommand(sitk.sitkIterationEvent, 
                                 lambda: rc.metric_and_reference_plot_values(RegMethod, FixPts, MovPts))
            
    FinalTx = RegMethod.Execute(FixIm, MovIm)
    
    RegIm = sitk.Resample(MovIm, FixIm, FinalTx, sitk.sitkLinear, 0.0, 
                          MovIm.GetPixelID())
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 1)
    Dtime_min = round((times[-1] - times[-2])/60, 1)
    if True:#LogToConsole:
        print(f'*Took {Dtime} s ({Dtime_min} min) to perform image registration.\n')
    
    """ Post registration analysis. """
    if LogToConsole:
        print("-------")
        print('InitialTx:')
        print(InitialTx)
        print("-------")
        print("-------")
        print('FinalTx:')
        print(FinalTx)
        print("-------")
    print(f'Final metric value: {RegMethod.GetMetricValue()}')
    print("Optimizer stopping condition:",
          f"{RegMethod.GetOptimizerStopConditionDescription()}")
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)

    return RegIm, InitialTx, FinalTx

def register_ims(FixIm, MovIm, RegTxName='affine', InitMethod='landmarks',
                 FixFidsFpath='', MovFidsFpath='', SamplingPercentage=5,
                 LogToConsole=False):
    """
    Wrapper function for functions SitkNonDefReg(), SitkBsplineReg() and 
    SitkBsplineRegWithFiducials() to register two 3D SimpleITK images using 
    SimpleITK.
    
    Parameters
    ----------
    FixIm : SimpleITK Image 
        The 3D image that MovIm will be registered to.
    MovIm : SimpleITK Image
        The 3D image that will be registered to FixIm.
    RegTxName : str, optional ('affine' by default)
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
    InitMethod : str, optional ('landmarks' by default)
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
        If InitMethod = 'landmarks' but file paths to txt files containing
        fiducials are not provided the next default value will be 'moments'.
    FixFidsFpath : str, optional ('' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for FixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    MovFidsFpath : str, optional ('' by default)
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for MovIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if InitMethod = 'landmarks'.
    SamplingPercentage : int or float, optional (5 by default)
        The percentage of the images that are sampled for evaluating the metric.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
        
    Returns
    -------
    RegIm : SimpleITK Image
        The 3D registered image.
    InitialTx : SimpleITK Transform
        The SimpleITK Transform used to initialise the registration.
    FinalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform MovIm to FixIm.
    
    Note
    ----
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Possible ideas from:
    https://discourse.itk.org/t/segmentation-fault-error-when-running-rigid-affine-registration/3879
    """
    
    if not RegTxName in ['rigid', 'affine', 'bspline']:
        msg = f"The chosen transform, {RegTxName}, is not one of the "\
              + "accepted inputs: 'rigid', 'affine', or 'bspline'."
        raise Exception(msg)
    
    if not InitMethod in ['geometry', 'geometricalcenter', 'moments', 
                          'centerofgravity', 'landmarks']:
        msg = f"The chosen initialisation method (InitMethod), {InitMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
    
    if RegTxName == 'bspline':
        if LogToConsole:
            print('Running bspline_reg()...\n')
            
        RegIm, InitialTx, FinalTx\
            = bspline_reg_ims(FixIm=FixIm, MovIm=MovIm,
                              FixFidsFpath=FixFidsFpath, 
                              MovFidsFpath=MovFidsFpath,
                              NumControlPts=8,
                              SamplingPercentage=SamplingPercentage,
                              NumIters=100, LearningRate=5.0, 
                              LogToConsole=LogToConsole)
    else:
        if LogToConsole:
            print('Running non_def_reg_ims()...\n')
        
        RegIm, InitialTx, FinalTx,\
            = non_def_reg_ims(FixIm=FixIm, MovIm=MovIm, RegTxName=RegTxName,
                              FixFidsFpath=FixFidsFpath, 
                              MovFidsFpath=MovFidsFpath,
                              SamplingPercentage=SamplingPercentage,
                              NumIters=500, LearningRate=1.0,
                              LogToConsole=LogToConsole)
    
    return RegIm, InitialTx, FinalTx