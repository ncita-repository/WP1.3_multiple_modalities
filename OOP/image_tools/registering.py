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

#from general_tools.general import check_file_ext
from general_tools.fiducials import get_landmark_tx
#import image_tools.registration_utilities as ru
import image_tools.registration_callbacks as rc


def command_iteration(method):
    print(f"{method.GetOptimizerIteration():3} = {method.GetMetricValue():10.5f}")
    #print("\t#: ", len(method.GetOptimizerPosition()))

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
    
    #""" 31/08/21
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()
    #"""

def plot_values(regMethod):
    """ Callback invoked when the IterationEvent happens, update the data 
    and display new figure. """
    
    global metric_values, multires_iterations
    
    metric_values.append(regMethod.GetMetricValue())                                       
    # Clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True) # <-- this clears all output
    #plt.clf() # doesn't have intended result
    # Plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, 
             [metric_values[index] for index in multires_iterations], 
             'b*')
    plt.xlabel('Iteration Number', fontsize=12)
    plt.ylabel('Metric Value', fontsize=12)
    plt.show()

def update_multires_iterations():
    """ Callback invoked when the sitkMultiResolutionIterationEvent happens, 
    update the index into the metric_values list. """
    
    global metric_values, multires_iterations
    
    multires_iterations.append(len(metric_values))

def update_metric_values(regMethod):
    """ Callback invoked when the IterationEvent happens, update the data. 
    """
    
    global metric_values, multires_iterations
    
    metric_values.append(regMethod.GetMetricValue())

def plot_fix_mov_ims(fixInd, movInd, fixPixarr, movPixarr):
    """ Callback invoked by the interact IPython method for scrolling 
    through the image stacks of the two images (moving and fixed). """
    
    # Create a figure with two subplots and the specified size.
    plt.subplots(1, 2, figsize=(10,8))
    
    # Draw the fixed image in the first subplot.
    plt.subplot(1, 2, 1)
    plt.imshow(fixPixarr[fixInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Fixed image')
    plt.axis('off')
    
    # Draw the moving image in the second subplot.
    plt.subplot(1, 2, 2)
    plt.imshow(movPixarr[movInd,:,:], cmap=plt.cm.Greys_r);
    plt.title('Moving image')
    plt.axis('off')
    
    plt.show()

def plot_blended_im(Ind, alpha, fixIm, resIm):
    """ Callback invoked by the IPython interact method for scrolling and 
    modifying the alpha blending of an image stack of two images that 
    occupy the same physical space. """

    BlendedIm = (1.0 - alpha)*fixIm[:,:,Ind] + alpha*resIm[:,:,Ind] 
    plt.imshow(sitk.GetArrayViewFromImage(BlendedIm), cmap=plt.cm.Greys_r);
    plt.axis('off')
    plt.show()

def plot_values_and_export(metric_values, multires_iterations, regPlotFpath):
    
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, 
             [metric_values[ind] for ind in multires_iterations], 'b*')
    plt.xlabel('Iteration number', fontsize=12)
    plt.ylabel('Metric value', fontsize=12)
    
    plt.savefig(regPlotFpath, bbox_inches='tight')
        
    print(f'Plot exported to:\n {regPlotFpath}\n')

def rigid_reg_im(
        fixIm, movIm, regTxName='affine', initMethod='landmarks', 
        fixFidsFpath='', movFidsFpath='', samplingPercentage=50, 
        numIters=500, learningRate=1.0, p2c=False, regPlotFpath='',
        ):
    """  
    Register two 3D SimpleITK images using a non-deformable transformation 
    using SimpleITK.
    
    Parameters
    ----------
    fixIm : SimpleITK Image 
        The 3D image that movIm will be registered to
    movIm : SimpleITK Image
        The 3D image that will be registered to fixIm.
    regTxName : str, optional
        Transformation to use for registration.  Acceptable values include:
        - 'rigid'
        - 'affine'
        The default value is 'affine'.
    initMethod : str, optional
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
        If initMethod = 'landmarks' but file paths to txt files containing
        fiducials are not provided the next default value will be 'moments'.
        The default value is 'landmarks'.
    fixFidsFpath : str, optional
        The file path of the text file containing fiducials for fixIm. 
        This parameter is only relevant if initMethod = 'landmarks'. 
        The default value is ''.
    movFidsFpath : str, optional
        The file path of the text file containing fiducials for movIm. 
        This argument is only relevant if initMethod = 'landmarks'.
        The default value is ''.
    samplingPercentage : int or float, optional
        The percentage of the images that are sampled for evaluating the metric.
        The default value is 50.
    numIters : int, optional
        The maximum number of iterations used when optimising. The default 
        value is 500.
    learningRate : float, optional
        The learning rate used for the gradient descent optimiser. The default 
        value is 1.0.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False
    regPlotFpath : str, optional
        The file path to be assigned to the Metric v Iteration number plot
        following optimisation if regPlotFpath is not ''. 
        The default value is ''.
        
    Returns
    -------
    regIm : SimpleITK Image
        The 3D registered image.
    initialTx : SimpleITK Transform
        The SimpleITK Transform used to initialise the registration.
    finalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform movIm to fixIm.
    
    Note
    ----
    It appears that a high sampling percentage (e.g. 50%) is required to get 
    good (and repeatable) results. Low values (e.g. 5%) result in poorer 
    results and are highly variable from run to run (due to the use of clock
    to generate samples).
    
    An indication that the sampling percentage is too low is looking at the
    Metric value v Iteration Number during optimisation. A highly erratic plot
    may suggest that the value is too low. Typically the part of the plot 
    from the first resolution (of the multi-resolution strategy) looks very 
    "noisy" (the part of the plot between the first and second asterixes). 
    
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Since initial alignment is geometrical (use of fiducials), initialising the
    registration with a geometrical centering is suitable (moments is not).
    
    It seems that a higher number of samples (e.g. 50%) is required to get
    results that are not so dependent upon the random sampling.
    """
    
    # Print simple metric values to console or plot of metric value v iteration?
    plotToConsole = False # simple text only
    plotToConsole = True # plots
    
    """ metric values and multiresIters don't return expected values when
    plotToConsole = False!
    """
    # TODO correct code so that metric values and iters returned when plotToConsole=False
    
    samplingPercentage = float(samplingPercentage/100)
    
    if regTxName == 'rigid':
        sitkTx = sitk.Euler3DTransform()
    elif regTxName == 'rigid_scale':
        sitkTx = sitk.Similarity3DTransform()
    elif regTxName == 'affine':
        sitkTx = sitk.AffineTransform(fixIm.GetDimension())
    else:
        msg = f"regTxName must be either 'rigid' or 'affine'. {regTxName} is "\
              + "not an accepted value."
        raise Exception(msg)
    
    # Start timing:
    times = []
    times.append(time.time())
    
    if initMethod == 'landmarks' and fixFidsFpath and movFidsFpath:
        """
        if (check_file_ext(fixFidsFpath, '.txt') and
                check_file_ext(movFidsFpath, '.txt')):
        """
        # TODO tidy this up
        if True: # 07/09 above returning false when true expected
            print('File paths to fiducials were provided. Attempting to',
                  'create a landmark transform..\n')
            
            # Get landmark transform:
            initialTx, fixPts, movPts = get_landmark_tx(
                txName=regTxName, 
                fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath, 
                fixIm=fixIm, movIm=movIm, flipK=False, p2c=True
                )
        else:
            print('File paths to fiducials were not provided/found.',
                  'Initialising using geometry..\n')
        
            #CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
            CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
            
            initialTx = sitk.CenteredTransformInitializer(
                fixIm, movIm, sitkTx, CTIF)
    else:
        if initMethod in ['geometry', 'geometricalcenter']:
            CTIF = sitk.CenteredTransformInitializerFilter.GEOMETRY
        elif initMethod in ['moments', 'centerofgravity']:
            CTIF = sitk.CenteredTransformInitializerFilter.MOMENTS
        
        initialTx = sitk.CenteredTransformInitializer(
            fixIm, movIm, sitkTx, CTIF
            )
    
    
        
    """ Registration. """
    regMethod = sitk.ImageRegistrationMethod()
    
    regMethod.SetInitialTransform(initialTx, inPlace=False)
    
    """ Similarity metric settings. """
    regMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    #regMethod.SetMetricAsJointHistogramMutualInformation()
    #print('\n\n\n*** Changing MetricSamplingStrategy from RANDOM to NONE (08/09/21)\n\n\n')
    regMethod.SetMetricSamplingStrategy(regMethod.RANDOM) # 08/09
    #regMethod.SetMetricSamplingStrategy(regMethod.REGULAR)
    #regMethod.SetMetricSamplingStrategy(regMethod.NONE) # 08/09
    regMethod.SetMetricSamplingPercentage(samplingPercentage) # 08/09
    #regMethod.SetMetricSamplingPercentage(
    #    samplingPercentage, sitk.sitkWallClock
    #    )
    
    """ Setup for the multi-resolution framework. """
    shrinkFactors = [4,2,1]
    smoothingSigmas = [2,1,1]
    #smoothingSigmas = [2,1,0]
    #print('\n\n\n*** Changed smoothingSigmas from [2,1,1] to [4,2,1] on',
    #      '08/09/21\n\n\n')
    #smoothingSigmas = [4,2,1]
    
    # 06/09/21
    #shrinkFactors = [8,4,2,1]
    #smoothingSigmas = [4,2,1,1]
    
    #""" 06/09/21
    regMethod.SetShrinkFactorsPerLevel(shrinkFactors=shrinkFactors)
    regMethod.SetSmoothingSigmasPerLevel(smoothingSigmas=smoothingSigmas)
    regMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    #"""
    
    """ Interpolator settings. """
    regMethod.SetInterpolator(sitk.sitkLinear)
    
    """ Optimizer settings. """
    
    """
    regMethod.SetOptimizerAsGradientDescent(
        learningRate=learningRate, numberOfIterations=500, 
        convergenceMinimumValue=1e-6, convergenceWindowSize=10
        )
    """
    
    #""" 06/09/21
    regMethod.SetOptimizerAsGradientDescent(
        learningRate=learningRate, numberOfIterations=numIters, 
        estimateLearningRate=regMethod.Once
        )
    #"""
    # 06/09/21:
    #regMethod.SetOptimizerAsRegularStepGradientDescent(1.0, .001, 200)
    #regMethod.SetOptimizerAsRegularStepGradientDescent(1.0, 1e-6, 200)
    #regMethod.SetOptimizerAsRegularStepGradientDescent(1.0, 1e-8, 200)
    
    """
    regMethod.SetOptimizerAsGradientDescentLineSearch(
        learningRate=learningRate, numberOfIterations=100,
        convergenceMinimumValue=1e-4, #convergenceMinimumValue=1e-6,
        convergenceWindowSize=5
        )
    """
    
    regMethod.SetOptimizerScalesFromPhysicalShift()
    
    if plotToConsole:
        # Connect all of the observers for plotting during registration:
        regMethod.AddCommand(sitk.sitkStartEvent, start_plot)
        # Don't want to delete metric_values or multires_iterations since
        # I want to return them:
        #regMethod.AddCommand(sitk.sitkEndEvent, end_plot)
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                             update_multires_iterations) 
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             lambda: plot_values(regMethod))
    else:
        # Need to run start_plot to initialise metric_values and
        # multires_iterations, even if not plotting:
        regMethod.AddCommand(sitk.sitkStartEvent, start_plot)
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                             update_multires_iterations)
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             update_metric_values(regMethod))
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             lambda: command_iteration(regMethod))
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                             lambda: command_multi_iteration(regMethod))
    
    # Get the number of valid points:
    numOfValidPts = regMethod.GetMetricNumberOfValidPoints()
    
    # Get the number of points in fixIm:
    numOfPtsInFixIm = np.prod(fixIm.GetSize())
    
    # The % of valid points:
    percentageValidPts = 100*numOfValidPts/numOfPtsInFixIm
    
    print(f'numOfValidPts = {numOfValidPts} = {percentageValidPts}%\n')
    
    #finalTx = regMethod.Execute(sitk.Cast(fixIm, sitk.sitkFloat32), 
    #                            sitk.Cast(movIm, sitk.sitkFloat32))
    finalTx = regMethod.Execute(fixIm, movIm)
    
    print('\nParameters used during registration:')
    print(f'regTxName = {regTxName}')
    print(f'initMethod = {initMethod}')
    print(f'samplingPercentage = {samplingPercentage}')
    print(f'learningRate = {learningRate}')
    print(f'numIters = {numIters}\n')
    
    """ Post registration analysis. """
    if p2c:
        if plotToConsole:
            plot_values_and_export(
                metric_values, multires_iterations, regPlotFpath
                )
        print("-------")
        print('initialTx:')
        print(initialTx)
        print("-------")
        print("-------")
        print('finalTx:')
        print(finalTx)
        print("-------")
    metricVal = regMethod.GetMetricValue()
    stopConditionDesc = regMethod.GetOptimizerStopConditionDescription()
    #finalIterNum = regMethod.GetOptimizerIteration()
    print(f'Final metric value: {metricVal}')
    print(f"Optimizer stopping condition: {stopConditionDesc}")
    #print(f"Final Iteration: {finalIterNum}")
    
    # Resample movIm using finalTx to get the registered image:
    regIm = sitk.Resample(
        movIm, fixIm, finalTx, sitk.sitkLinear, 0.0, movIm.GetPixelID()
        )
    
    if False:
        interact(plot_blended_im, Ind=(0,fixIm.GetSize()[2] - 1), 
                 alpha=(0.0,1.0,0.05), fixIm=(fixIm), resIm=fixed(regIm));
    
    times.append(time.time())
    #dTime = round(times[-1] - times[-2], 1)
    dTime = times[-1] - times[-2]
    if True:#p2c:
        print(f'*Took {dTime:.1f} s ({dTime/60:.1f} min) to perform image',
              'registration.\n')
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    # Resample movIm using intialTx to get the pre-registration aligned image 
    # (for info only - not required):
    alignedIm = sitk.Resample(
        movIm, fixIm, initialTx, sitk.sitkLinear, 0.0, movIm.GetPixelID()
        )
    
    return initialTx, alignedIm, finalTx, regIm, metric_values,\
        multires_iterations

def bspline_reg_im(
        fixIm, movIm, fixFidsFpath='', movFidsFpath='', numControlPts=8,
        samplingPercentage=5, numIters=100, learningRate=5.0, p2c=False,
        regPlotFpath=''
        ):
    """ 
    Register two 3D SimpleITK images using a B-spline transformation in
    SimpleITK.
    
    Parameters
    ----------
    fixIm : SimpleITK Image 
        The 3D image that movIm will be registered to.
    movIm : SimpleITK Image
        The 3D image that will be registered to fixIm.
    fixFidsFpath : str, optional
        The file path of the txt file containing fiducials for fixIm. If '' the
        BSpline will be initialised using BSplineTransformInitializer.
        The default value is ''.
    movFidsFpath : str, optional
        The file path of the txt file containing fiducials for movIm. If '' the
        BSpline will be initialised using BSplineTransformInitializer.
        The default value is ''.
    numControlPts : int, optional
        The number of control points that defines the BSpline grid. The default
        value is 8.
    samplingPercentage : int or float, optional
        The percentage of the images that are sampled for evaluating the metric.
        The default value is 5.
    numIters : int, optional
        The maximum number of iterations used during optimisation. The default
        value is 100.
    learningRate : float, optional
        The learning rate used for the gradient descent optimiser. The default
        value is 5.0.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False.
    regPlotFpath : str, optional
        The file path to be assigned to the Metric v Iteration number plot
        following optimisation if regPlotFpath is not ''. 
        The default value is ''.
        
    Returns
    -------
    regIm : SimpleITK Image
        The 3D registered image.
    regMethod : SimpleITK ImageRegistrationMethod
        The SimpleITK ImageRegistrationMethod used to perform registration.
    initialTx : SimpleITK Transform
        The Bspline SimpleITK Transform used to initialise the BSpline grid
        prior to registration.
    finalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform movIm to fixIm.
    
    Note
    ----
    Unlike rigid registrations it seems that a low sampling percentage (e.g.
    5%) is sufficient to get good results.
    
    Zivy advised not to use gradient descent line search for bspline:
    https://github.com/SimpleITK/SimpleITK/issues/1415
    """
    
    # Print simple metric values to console or plot of metric value v iteration?
    plotToConsole = False # simple text only
    plotToConsole = True # plots
    
    # Use scale factors for mesh grid?
    # Note 08/07/21: Not sure it works properly using scale factors.
    useScaleFactors = False
    #useScaleFactors = True
    
    """ Scale factors: 
    
    Maybe the scale factors should take into account the NumOfControlPts
    that were used in LandmarkTx and the desired NumOfControlPts for the 
    registration?
    """
    scaleFactors = [1, 2, 4]
    
    dim = fixIm.GetDimension()
    
    if useScaleFactors:
        meshSize = [numControlPts // scaleFactors[-1]] * dim
    else:
        meshSize = [numControlPts] * dim
    
    samplingPercentage = samplingPercentage/100
    
    # Start timing:
    times = []
    times.append(time.time())
    
    if fixFidsFpath and movFidsFpath:
        """
        if (check_file_ext(fixFidsFpath, '.txt') and
                check_file_ext(movFidsFpath, '.txt')):
            print('File paths to fiducials were provided. Attempting to create a',
                  'landmark transform..\n')
        """
        # TODO tidy this up
        if True: # 07/09 above returning false when true expected
            
            # Get landmark transform:
            initialTx, fixPts, movPts = get_landmark_tx(
                txName='bspline', fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath, fixIm=fixIm, movIm=movIm, 
                flipK=False, buffer=2, p2c=True
                )
    else:
        print('File paths to fiducials were not provided/found. Initialising',
              'the Bspline using BSplineTransformInitializer..\n')
        
        initialTx = sitk.BSplineTransformInitializer(
            image1=fixIm, transformDomainMeshSize=meshSize
            )
    
    print(f'initialTx is {initialTx.GetName()}\n')
    
    """ Registration. """
    regMethod = sitk.ImageRegistrationMethod()
    
    if useScaleFactors:
        regMethod.SetInitialTransformAsBSpline(
            initialTx, inPlace=True, scaleFactors=scaleFactors
            )
    else:
        regMethod.SetInitialTransform(initialTx, inPlace=True)
    
    """ Similarity metric settings. """
    regMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    regMethod.SetMetricSamplingStrategy(regMethod.RANDOM)
    regMethod.SetMetricSamplingPercentage(samplingPercentage)
    
    """ Optimiser settings. Use LBFGSB (without scale factors) or LBFGS2 (with 
    scale factors).
    """
    if useScaleFactors:
        """ Use the LBFGS2 instead of LBFGS. The latter cannot adapt to the 
        changing control grid resolution. """
        regMethod.SetOptimizerAsLBFGS2(
            solutionAccuracy=1e-2, numberOfIterations=numIters, 
            deltaConvergenceTolerance=0.01
            )
    else:
        regMethod.SetOptimizerAsLBFGSB(
            gradientConvergenceTolerance=1e-5, numberOfIterations=numIters
            )
    regMethod.SetOptimizerScalesFromPhysicalShift()
    
    regMethod.SetInterpolator(sitk.sitkLinear)
    
    #regMethod.SetShrinkFactorsPerLevel([6, 2, 1])
    #regMethod.SetSmoothingSigmasPerLevel([6, 2, 1])
    regMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    regMethod.SetSmoothingSigmasPerLevel([4, 2, 1]) # <-- better result than [2, 1, 0]?..
    #regMethod.SetSmoothingSigmasPerLevel([2, 1, 0]) # in 65_Registration_FFD.ipynb example
    regMethod.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn() # <-- THIS WAS MISSING until 21:51 on 1/7/21
    
    if plotToConsole:
        # Connect all of the observers so that we can perform plotting  
        # during registration.
        regMethod.AddCommand(sitk.sitkStartEvent, start_plot)
        # Don't want to delete metric_values or multires_iterations since
        # I want to return them:
        #regMethod.AddCommand(sitk.sitkEndEvent, end_plot)
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                             update_multires_iterations) 
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             lambda: plot_values(regMethod))
    else:
        # Need to run start_plot to initialise metric_values and
        # multires_iterations, even if not plotting:
        regMethod.AddCommand(sitk.sitkStartEvent, start_plot)
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent, 
                             update_multires_iterations)
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             update_metric_values(regMethod))
        regMethod.AddCommand(sitk.sitkIterationEvent, 
                             lambda: command_iteration(regMethod))
        regMethod.AddCommand(sitk.sitkMultiResolutionIterationEvent,
                             lambda: command_multi_iteration(regMethod))
    
        # Since corresponding points in the fixed and moving image are available  
        # display the similarity metric and the TRE during the registration:
        """ TRY THIS - COULD BE USEFUL: """
        if False:
            regMethod.AddCommand(
                sitk.sitkStartEvent, rc.metric_and_reference_start_plot
                )
            regMethod.AddCommand(
                sitk.sitkEndEvent, rc.metric_and_reference_end_plot
                )
            regMethod.AddCommand(
                sitk.sitkIterationEvent, 
                lambda: rc.metric_and_reference_plot_values(
                    regMethod, fixPts, movPts
                    )
                )
            
    finalTx = regMethod.Execute(fixIm, movIm)
    
    regIm = sitk.Resample(
        movIm, fixIm, finalTx, sitk.sitkLinear, 0.0, movIm.GetPixelID()
        )
    
    times.append(time.time())
    dTime = round(times[-1] - times[-2], 1)
    dTime_min = round((times[-1] - times[-2])/60, 1)
    if True:#p2c:
        print(f'*Took {dTime} s ({dTime_min} min) to perform image registration.\n')
    
    """ Post registration analysis. """
    if p2c:
        if plotToConsole:
            plot_values_and_export(
                metric_values, multires_iterations, regPlotFpath
                )
        print("-------")
        print('initialTx:')
        print(initialTx)
        print("-------")
        print("-------")
        print('finalTx:')
        print(finalTx)
        print("-------")
    print(f'Final metric value: {regMethod.GetMetricValue()}')
    print("Optimizer stopping condition:",
          f"{regMethod.GetOptimizerStopConditionDescription()}")
    
    frequency = 2500  # Hz
    duration = 1500  # ms
    winsound.Beep(frequency, duration)
    
    # Resample movIm using intialTx to get the pre-registration aligned image 
    # (for info only - not required):
    alignedIm = sitk.Resample(
        movIm, fixIm, initialTx, sitk.sitkLinear, 0.0, movIm.GetPixelID()
        )
    
    return initialTx, alignedIm, finalTx, regIm, metric_values,\
        multires_iterations

def register_im(
        fixIm, movIm, regTxName='affine', initMethod='landmarks',
        fixFidsFpath='', movFidsFpath='', p2c=False, regPlotFpath=''
        ):
    """
    Wrapper function for functions rigid_reg_im() and bspline_reg_im() 
    to register two 3D SimpleITK images using SimpleITK.
    
    Parameters
    ----------
    fixIm : SimpleITK Image 
        The 3D image that movIm will be registered to.
    movIm : SimpleITK Image
        The 3D image that will be registered to fixIm.
    regTxName : str, optional
        Denotes type of transformation to use for registration.  Acceptable 
        values include:
        - 'rigid'
        - 'affine'
        - 'bspline'
        The default value is 'affine'.
    initMethod : str, optional
        The method used for initialisation.  Acceptable inputs include:
            - 'geometry' (or 'geometricalcenter')
            - 'moments' (or 'centerofgravity')
            - 'landmarks'
        The default value is 'landmarks'. If initMethod = 'landmarks' but file 
        paths to txt files containing fiducials are not provided the default 
        value will be 'moments'.
    fixFidsFpath : str, optional unless initMethod=='landmarks'
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for fixIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if initMethod = 'landmarks'. The default value is ''.
    movFidsFpath : str, optional unless initMethod=='landmarks'
        The file path (or file name if the file is present in the current 
        working directory) of the text file containing fiducials for movIm. 
        The string need not contain the .txt extension. This argument is only
        relevant if initMethod = 'landmarks'. The default value is ''.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False.
    regPlotFpath : str, optional
        The file path to be assigned to the Metric v Iteration number plot
        following optimisation if regPlotFpath is not ''. 
        The default value is ''.
        
    Returns
    -------
    regIm : SimpleITK Image
        The 3D registered image.
    initialTx : SimpleITK Transform
        The SimpleITK Transform used to initialise the registration.
    finalTx : SimpleITK Transform
        The final SimpleITK Transform used to register/transform movIm to fixIm.
    
    Note
    ----
    https://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/60_Registration_Introduction.html
    
    https://simpleitk.readthedocs.io/en/master/link_ImageRegistrationMethodBSpline3_docs.html
    
    Possible ideas from:
    https://discourse.itk.org/t/segmentation-fault-error-when-running-rigid-affine-registration/3879
    """
    
    # Deal with invalid value for regTxName:
    if not regTxName in ['rigid', 'rigid_scale', 'affine', 'bspline']:
        msg = f"The chosen transform, {regTxName}, is not one of the "\
              + "accepted inputs: 'rigid', 'rigid_scale', 'affine', or "\
              + "'bspline'."
        raise Exception(msg)
    
    # Deal with invalid value for initMethod:
    if not initMethod in [
            'geometry', 'geometricalcenter', 'moments', 'centerofgravity', 
            'landmarks'
            ]:
        msg = f"The chosen initialisation method (initMethod), {initMethod}, "\
              + "is not one of the accepted inputs: 'geometry'/'geometricalcenter',"\
              + "  'moments'/'centerofgravity', or 'landmarks'."
        raise Exception(msg)
    
    # Deal with missing fiducials file (i.e. if one or both fiducials file 
    # paths are '' or None) for initMethod = 'landmarks':
    if initMethod == 'landmarks' and not (fixFidsFpath or movFidsFpath):
        msg = f"The chosen initialisation method (initMethod), {initMethod}, "\
              + "requires both fixed and moving fiducial file paths."\
                  + f"\nfixFidsFpath = {fixFidsFpath}"\
                      + f"\nmovFidsFpath = {movFidsFpath}"
        raise Exception(msg)
    
    if regTxName == 'bspline':
        if p2c:
            print('Running bspline_reg()...\n')
            
        initialTx, alignedIm, regIm, finalTx, metric_values,\
            multires_iterations = bspline_reg_im(
                fixIm=fixIm, movIm=movIm,
                fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath,
                numControlPts=8, samplingPercentage=5,
                numIters=100, learningRate=5.0, p2c=p2c, regPlotFpath=regPlotFpath
                )
    else:
        if p2c:
            print('Running rigid_reg_im()...\n')
        
        """
        learningRate = 1
        learningRate = 5
        #learningRate = 10
        #learningRate = 0.5
        if learningRate != 1:
            print('\n* On 06/09 changed learning rate from 1 to',
                  f'{learningRate}\n\n')
        """
        
        """
        samplingPercentage = 50
        #samplingPercentage = 100
        if samplingPercentage != 50:
            print('\n* On 06/09 changed samplingPercentage from 50% to',
                  f'{samplingPercentage}%\n\n')
        """
        
        """
        initMethod = 'moments'
        #initMethod = 'geometry'
        if initMethod != 'moments':
            print('\n* On 06/09 changed initMethod from "moments" to',
                  '"geometry"\n\n')
        """
        
        learningRate = 1.0
        samplingPercentage = 50
        
        initialTx, alignedIm, finalTx, regIm, metric_values,\
            multires_iterations = rigid_reg_im(
                fixIm=fixIm, movIm=movIm, 
                regTxName=regTxName, initMethod=initMethod,
                fixFidsFpath=fixFidsFpath, 
                movFidsFpath=movFidsFpath,
                samplingPercentage=samplingPercentage,
                numIters=500, learningRate=learningRate, 
                p2c=p2c, regPlotFpath=regPlotFpath
                )
    
    return initialTx, alignedIm, finalTx, regIm, metric_values,\
        multires_iterations