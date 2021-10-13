# WP1.3_multiple_modalities

### NCITA Repository Unit Work Package 1.3: XNAT-OHIF viewer - multiple modalities
Specific repositories and projects of the form WP1.x-y are for specific second-level, viewer-related work packages in the NCITA Repository Unit strategy.

The aim of this work package of the NCITA Repository Unit will be to develop approaches within the OHIF viewer and XNAT to work with images from connected studies. The initial focus will be on image registration and the copying of regions-of-interest between:
* different DICOM series within the same imaging session (e.g., draw ROI on DWI and see where it is on the corresponding T1w image);
* corresponding images from different time points for the same patient (registration required and we will investigate rigid and non-rigid);
* corresponding images from different modalities (start with the simple case of inherently co-registered PET-CT and work towards more difficult scenarios).

The NCITA Repository Unit site as a whole is at an early stage of the life cycle and does not currently have a formal structure of documentation releases. Please bear with us: things will become more organised as the NCITA project gets into gear! Note that the ICR also releases code related to the XNAT-OHIF image viewer at https://bitbucket.org/icrimaginginformatics/. 

### Basic instructions on running the code

A *source* ROI Collection (RTSTRUCT or SEG) is copied (or propagated) to a *target* DICOM series.  One such single operation will be referred to as a *run*.  When executing the code, a `runID` must be provided, and the directory path of the configuration file that contains the metadata attributed to that run, `cfgDir`.  The steps required to perform a run are as follows:

1. Create a configuration file(s) for the run(s) you wish to carry out.

2. Run the code on each run, providing a `runID` and `cfgDir`.

#### 1. Creating a configuration file

The package [*config.py*](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/config.py) is used to generate configuration files.  Within the *config.py* module is the function *create_config_files()*.  This contains pre-populated data specific to the XNAT used to develop the code. You can use the pre-populated code as a template for your own XNAT.

There are several variables to be defined, including:

- Where the fetched data is to be downloaded to (`rootDownloadDir`)
- Where any outputs should be exported (e.g. `rtsExportDir`)
- Where to find sample DROs and fiducials (e.g. `sampleDroDir` and `fidsDir`)
- Whether or not to export various data (e.g. `exportRoicol`)
- Whether or not to print verbose results to the console (`p2c`)

A dictionary `cfg` is initialised containing the above variables across the board, but there are additional variables that are set within the dictionary elements themselves.  The key-dictionary pairs that are defined within the dictionary include:

- `url`           : XNAT address
- `username`      : XNAT username
- `password`      : XNAT password
- `proID`         : XNAT Project ID
- `subjLab`       : XNAT Subject Label
- `srcExpLab`     : Source dataset's XNAT Experiment Label
- `srcScanID`     : Source dataset's XNAT Scan ID
- `srcSlcNum`     : Slice number of the Source ROI Collection to be copied. Leave as `None` if all contours/segmentations are to be copied.
- `srcRoicolName` : Source ROI Collection's Name
- `srcRoiName`    : Name of the ROI/segment within `srcRoicolName` to be copied. Leave as `None` if all ROIs/segments are to be copied.
- `roicolMod`     : Source ROI Collection's Modality
- `trgExpLab`     : Target dataset's XNAT Experiment Label
- `trgScanID`     : Target dataset's XNAT Scan ID
- `trgSlcNum`     : Slice number destination where the contour/segmentation is to be copied to. Leave as `None` if all contours/segmentations are to be copied.
- `trgRoicolName` : Target ROI Collection's Name. Leave as `None` if one does not exist (i.e. one does not wish to copy a specific contour/segmentation within the Source ROI Collection to an existing Target ROI Collection). Such an operation is only valid for "direct" copies (see below for further details).
- `trgRoiName`    : Name of the ROI/segment within `trgRoicolName` where a contour/segmentation is to be copied to. Leave as `None` if all ROIs/segments are to be copied.
- `srcFidsFpath`  : Full filepath to the .txt containing fiducials for the Source DICOM scan.
- `trgFidsFpath`  : Full filepath to the .txt containing fiducials for the Target DICOM scan.
- `forceReg`      : Set to `True` if you wish to apply image registration on a run for which the Source and Target scans have the same *FrameOfReferenceUID*.
- `regTxName`     : Defines the transformation to apply during registration. Acceptable values include "rigid", "affine" and "bspline".
- `initMethod`    : Defines how a registration is initialised.  Acceptable values include "geometry" (aligning to image centres), "moments" (align to centre-of-mass of image intensity), and "landmarks" (`srcFidsFpath` and `trgFidsFpath` must be defined, otherwise will default to "geometry").
- `maxIters`      : The maximum number of iterations to allow during optimisation for image registation.
- `useDroForTx`   : If `True` and a suitable DRO is found on XNAT, the transformation parameters stored in the DRO will be used, by-passing image registration.
- `applyPreResBlur` : If `True` the Source label image will be Gaussian blurred prior to resampling or registration. Doing so helps to avoid aliasing effects.
- `preResVar`       : The variance of the Gaussian blur that will be applied if `applyPreResBlur` is `True`.
- `resInterp`       : Defines the interpolation to use when resampling label images.  Acceptable values include "NearestNeighbor", "LabelGaussian" and `BlurThenLinear`.  *NearestNeighbor* and *LabelGaussian* are binary interpolators, whereas *BlurThenLinear* applies a Gaussian blur (defined by `preResVar`), followed by linear resampling, then binary thresholding to retain a binary label image. This option helps to avoid aliasing effects.
- `applyPostResBlur` : If `True` the resampled label image will be Gaussian blurred.
- `postResVar`       : The variance of the Gaussian blur that will be applied if `applyPostResBlur` is `True`. 

The remaining keys are simply the variables that were assigned at the top of the module (e.g. `sampleDroDir`).

Once the dictionaries have been defined, the configuration files can be created in one of two ways:

1. In a command shell, enter (change path to *src* directory accordingly):

		cd C:\Code\WP1.3_multiple_modalities\src
		python config.py configs

2. Import the module and execute the function:
	
		from config import create_config_files
		config_files("C:\Code\WP1.3_multiple_modalities\src\config")

A list of .json files should be generated in the config directory - one file for each key-dictionary pair within `cfg`.

#### 2. Performing a run

The main script for performing a run is contained in the module *app.py* in the *src* directory.  The run can be executed either as a function within a python environment, i.e.

	from app import main
	main("C:\Code\WP1.3_multiple_modalities\src\configs", "NCITA_TEST_RR2")

or within a command shell as a script:

	python app.py C:\Code\WP1.3_multiple_modalities\src\configs NCITA_TEST_RR2

### Relationship-preserving v Non-relationship-preserving, and Copies v Propagations

A *non-relationship preserving* copy (also referred to as a "direct" copy) is analogous to a copy-and-paste function - e.g. "Copy the *source* contour/segmentation corresponding to slice `srcSlcNum` in ROI/segment `srcRoiName` in ROI Collection `srcRoicol` to a contour/segmentation that will overlay onto *target* slice `trgSlcNum`.  "Direct" copies do not maintain spatial relationships between the *source* contour/segmentation and the *target* contour/segmentation's location in space.  One might want to perform a "direct" copy to use an existing contour/segmentation as a "starting point" for a new one that can sculpted to match the anatomical features on the *target* scan.  By contrast, a *relationship-preserving* does respect spatial relationships and hence the contour(s)/segmentation(s) end up where they "should".

*Copies* are made between two 3D images that have the same voxel resolution (*Pixel Spacing* and difference between adjacent *ImagePositionPatient* values along the scan direction), same patient position (*ImagePositionPatient*), orientation (*ImageOrientationPatient*) and frame-of-reference (*FrameOfReferenceUID*).  Propagations require a bit more work - either a resampling is required (if the two 3D images do not have the same voxel spacings, or same extent, but have the same FOR), or image registration (if they do not have the same FOR).  

*Non-relationship preserving* copies/propagations always involve the copy/propagation of a single contour/segmentation to a single contour/segmentation.  If it is the case that resampling/registration of the source contour/segmentation led to multiple contours/segmentations, the ROI/segment is collapsed to a single contour/segmentation to maintain the expected behaviour.  On the contrary, *relationship-preserving* copies/propagations are broader, in that a single contour/segmentation may map to multiple contours/segmentations (depending on the relative voxel spacings and FOR of the two 3D images). 

When making relationship-preserving copies/propagations, in addition to copies/propagations of a single contour/segmentation, an entire ROI/segment, consisting of any number of contour(s)/segment(s) may be copied, or an entire ROI Collection, consisting of any number of ROI(s)/segment(s) containing any number of contour(s)/segmentation(s).  The behaviour entirely depends on the user-defined configuration parameters and relationships between the two image domains.

The code used to copy or propagate an ROI Collection for one image session to another relies heavily on the use of [*SimpleITK*](https://simpleitk.org/) [1-3].

[1]: R. Beare, B. C. Lowekamp, Z. Yaniv, “Image Segmentation, Registration and Characterization in R with SimpleITK”, J Stat Softw, 86(8), https://doi.org/10.18637/jss.v086.i08, 2018.

[2]: Z. Yaniv, B. C. Lowekamp, H. J. Johnson, R. Beare, “SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research”, J Digit Imaging., https://doi.org/10.1007/s10278-017-0037-8, 31(3): 290-303, 2018.

[3]: B. C. Lowekamp, D. T. Chen, L. Ibáñez, D. Blezek, “The Design of SimpleITK”, Front. Neuroinform., 7:45. https://doi.org/10.3389/fninf.2013.00045, 2013.
