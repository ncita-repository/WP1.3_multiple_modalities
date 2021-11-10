# NCITA Repository Unit Work Package 1.3: XNAT-OHIF viewer - multiple modalities
Specific repositories and projects of the form WP1.x-y are for specific second-level, viewer-related work packages in the NCITA Repository Unit strategy.

The aim of this work package of the NCITA Repository Unit will be to develop approaches within the OHIF viewer and XNAT to work with images from connected studies. The initial focus will be on image registration and the copying of regions-of-interest between:
* different DICOM series within the same imaging session (e.g., draw ROI on DWI and see where it is on the corresponding T1w image);
* corresponding images from different time points for the same patient (registration required and we will investigate rigid and non-rigid);
* corresponding images from different modalities (start with the simple case of inherently co-registered PET-CT and work towards more difficult scenarios).

The NCITA Repository Unit site as a whole is at an early stage of the life cycle and does not currently have a formal structure of documentation releases. Please bear with us: things will become more organised as the NCITA project gets into gear! Note that the ICR also releases code related to the XNAT-OHIF image viewer at https://bitbucket.org/icrimaginginformatics/. 

# Tools covered in this repository

This repository contains three main tools:

1. A "stand-alone" *ROI copy/propagation* tool

2. The "back-end" of the *ROI propagation* portion of a wider *ROI copy/propagation* tool to be run in the OHIF-Viewer

3. An XNAT "snapshot" tool

All code is written in Python 3.  The *ROI* in *ROI copy/propagation* refers to a region-of-interest in the broadest sense, and is not to be confused with a collection of contours within an RTSTRUCT-based *ROI Collection*.  The definitions at the end of this page will be useful in understanding the terminology used to describe the tools.  In brief, a "copy" will refer to an operation that does not preserve the spatial coordinates of the contour or segmentation, whereas a "propagation" implies that the resulting contour or segmentation upon pasting maintains spatial coordinates.

## "Stand-alone" *ROI copy/propagation* tool

The "stand-alone" ROI copy/propagation tool can be used to make "non-relationship-preserving" copies of a single *contour* or *segmentation* (think copy and paste function), or "relationship-preserving" copies of a single *contour* or *segmentation*, a collection of *contour*s (= "ROI") or *segmentation*s (= "segment"), or an entire ROI Collection (consisting of any number of *contour*s or *segmentation*s within any number of *ROI*s or *segment*s).

This tool relies on the user providing a configuration JSON file that contains the metadata that identifies the "source" and "target" DICOM series ("scans" in XNAT parlance), *source* ROI Collection, and if applicable, *target* ROI Collection.  Based on the user-inputed data, the tool will use XNAT REST API calls to fetch the required data, download it to the local filesystem, import the data and depending on the relationship between the *source* and *target* DICOM series and other user-provided metadata, will either perform a non-relationship-preserving copy or relationship-preserving propagation of the entity of interest - entity here means *contour*, *segmentation*, *ROI*, *segment* or *ROI Collection*.

## Back-end of the *ROI propagation* portion of a wider *ROI copy/propagation* tool to be run in the OHIF-Viewer

Current development of the OHIF-Viewer consists of a *contour* copy-and-paste feature (which will be extended to *segmentation* copying).  This is a non-relationship preserving copy in that the contour that is pasted will occupy the same location within the image extent (interpolation will be used in the *source* and *target* DICOM slices have different resolutions), but their coordinates within space will not be preserved.  Hence a *contour* that segmented the ventricle in the brain in the *source* DICOM slice will not necessarily coincide with the ventricle in the *target* slice upon the pasting operation.  Such non-relationship-preserving copies will be handled entirely within the OHIF-Viewer front-end.  Relationship-preserving propagations will be handled by the "back-end" in Python, for example, using image resampling or registration to preserve spatial relationships.





## Install required packages

After cloning the repository, *pip* install the required packages listed in the *requirements.txt* file:

    pip3 install -r requirements.txt


## Basic instructions on running the code

A *source* ROI Collection (DICOM-RTSTRUCT or DICOM-SEG) is copied (or propagated) to a *target* DICOM series.  One such single operation will be referred to as a *run*.  When executing the code, a `runID` must be provided, and the directory path of the configuration file that contains the metadata attributed to that run, `cfgDir`.  The steps required to perform a run are as follows:

1. Create a configuration file(s) for the run(s) you wish to carry out.

2. Run the code on each run, providing a `runID` and `cfgDir`.

### 1. Creating a configuration file

The package [*config.py*](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/config.py) is used to generate configuration files.  Within the *config.py* module is the function *create_config_files()*.  This contains pre-populated data specific to the XNAT that was used to develop the code. You can use the pre-populated code as a template for your own XNAT.

There are several variables to be defined, including:

- Where the fetched data is to be downloaded to (`rootDownloadDir`)
- Where any outputs should be exported (e.g. `rtsExportDir`)
- Where to find sample DROs and fiducials (e.g. `sampleDroDir` and `fidsDir`)
- Whether or not to export various data (e.g. `exportRoicol`)
- Whether or not to print (highly) verbose results to the console (`p2c`) 

*Note: Some comments/results will be printed to the console even if `p2c = False`.*

A dictionary of dictionaries, `cfg`, is populated with dictionaries that containing the above variables across the board, since they are variables that are less likely to be modified from run to run.  Each sub-dictionary in `cfg` has a `runID` that matches the key of the sub-dictionary within the `cfg` dictionary.  In addition to the above variables there are additional variables that are set within each sub-dictionary.  The key-value pairs that define each sub-dictionary includes:

- `url`           : XNAT address
- `username`      : XNAT username
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

There are a few ways to go about creating a configuration file(s):

#### Option 1 - Copy and modify an existing configuration file

This is by far the quickest way to go about implementing your first run of the code.  

Simply copy an existing configuration file from within the *.../src/configs* directory (e.g. *NCITA_test_RD1.json*), rename the file, open it in your preferred text editor and modify the dictionary's values accordingly.  Note that the the name you assign to the file must agree with the value given for the `runID` key.  

#### Option 2 - Define a new dictionary within the function *create_config_files()* in *config.py*

This option is more suitable if you wish to run the code on multiple datasets and you don't want to manually create confguration files one by one.  

Within the function *create_config_files()* there are cases where a dictionary has been fully defined (by scratch) - see for example the definition of the dictionary with key *NCITA_test_RR1* (lines 317-374).  In this case a `runID` is assigned (line 316), and a dictionary is appended to the dictionary `cfg` with that `runID`.  

Further down on line 379, a new dictionary is added to `cfg` using the key *NCITA_test_RR2* (line 378), and the previously defined dictionary is copied (line 379).  Then the `runID` (line 380) and other selected keys are modified (lines 381-385).

Use a previously fully-defined dictionary as a template for new dictionaries as you go along.  Remember to re-run the command *create_config_files()* to update your list of config files.

Once the dictionaries have been defined, the configuration files can be created in one of two ways:

##### Option 2a - Run *config.py* as a script in a command shell

In a command shell, cd to the *src* directory and run the Python module and provide the (relative or absolute) directory path containing the configuration file (change path to *src* directory accordingly):

	cd C:\Code\WP1.3_multiple_modalities\src
	python config.py configs

##### Option 2b - Run *create_config_files()* as a function within a Python shell

In a Python shell, import the function and provide the (absolute) directory path:
	
	import sys
	sys.path.append("C:\Code\WP1.3_multiple_modalities\src")
	from config import create_config_files
	create_config_files("C:\Code\WP1.3_multiple_modalities\src\configs")

A list of .json files should be generated and/or overwritten in the configs directory - one file for each key-dictionary pair within the dictionary `cfg`.  

### 2. Performing a run

Now that the configuration file has been created it's time to run the code.  As before you have two options - you can execute it in a command shell, or in a Python shell.

#### Option 1 - Run *app.py* as a script in a command shell

In a command shell, run the Python module and provide the (relative or absolute) directory path containing the configuration file and the `runID`:

	python app.py configs NCITA_TEST_RR2
	
#### Option 2 - Run the function *main()* in a Python shell

In a Python shell, run *main()* providing as input arguments the (absolute) directory path containing the configuration file and the `runID`:

	from app import main
	main("C:\Code\WP1.3_multiple_modalities\src\configs", "NCITA_TEST_RR2")

Regardless of which option you take, you should be prompted to enter an XNAT password.  If the XNAT `url`, `username` and password you just entered were correct, an XNAT session will be established and an XNAT alias token generated.  The alias token will be saved to the same *configs* directory.  If the code is executed again using the same `url` and authentication credentials, and within the lifetime of the alias token, the token will be used, avoiding the need to re-enter your password.  




# XNAT Snapshots

A basic tool that fetches high-level metadata from an XNAT and produces a "snapshot", exported as an XLSX file is a stand-alone feature separate from the business of copying/propagating ROIs.

At present two XLSX files are exported - one that provides a snapshot broken down by projects, and one that is XNAT-wide.  The former includes the principle investigator, project description, a list of investigators, number of subjects, experiments, MR/CT/PET sessions, etc.  The XNAT-wide table includes totals and averages per project of subjects, experiments, image sessions, DICOM image files, etc.

![An example of an XNAT-by-project snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_by_project_snapshot.png "An example of an XNAT-by-project snapshot")

![An example of an XNAT-wide snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_wide_snapshot.png "An example of an XNAT-wide snapshot")

The tool can either be run as a module in a command shell or as a Python function.

#### Option 1 - Run *snapshots.py* as a script in a command shell

In a command shell, run the Python module providing the XNAT url:

	python snapshots.py configs http://10.1.1.20
	
#### Option 2 - Run the function *get_xnat_snapshot()* in a Python shell

In a Python shell, run *get_xnat_snapshot()* providing the XNAT url as an input argument:

	from snapshots import get_xnat_snapshot
	get_xnat_snapshot("http://10.1.1.20")



# Definitions

## Relationship-preserving v Non-relationship-preserving, and Copies v Propagations

A *non-relationship preserving* copy (also referred to as a "direct" copy) is analogous to a copy-and-paste function - e.g. "Copy the *source* contour/segmentation corresponding to slice `srcSlcNum` in ROI/segment `srcRoiName` in ROI Collection `srcRoicol` to a contour/segmentation that will overlay onto *target* slice `trgSlcNum`.  "Direct" copies do not maintain spatial relationships between the *source* contour/segmentation and the *target* contour/segmentation's location in space.  One might want to perform a "direct" copy to use an existing contour/segmentation as a "starting point" for a new one that can sculpted to match the anatomical features on the *target* scan.  By contrast, a *relationship-preserving* does respect spatial relationships and hence the contour(s)/segmentation(s) end up where they "should".

*Copies* are made between two 3D images that have the same voxel resolution (*Pixel Spacing* and difference between adjacent *ImagePositionPatient* values along the scan direction), same patient position (*ImagePositionPatient*), orientation (*ImageOrientationPatient*) and frame-of-reference (*FrameOfReferenceUID*).  Propagations require a bit more work - either a resampling is required (if the two 3D images do not have the same voxel spacings, or same extent, but have the same FOR), or image registration (if they do not have the same FOR).  

*Non-relationship preserving* copies/propagations always involve the copy/propagation of a single contour/segmentation to a single contour/segmentation.  If it is the case that resampling/registration of the source contour/segmentation led to multiple contours/segmentations, the ROI/segment is collapsed to a single contour/segmentation to maintain the expected behaviour.  On the contrary, *relationship-preserving* copies/propagations are broader, in that a single contour/segmentation may map to multiple contours/segmentations (depending on the relative voxel spacings and FOR of the two 3D images). 

When making relationship-preserving copies/propagations, in addition to copies/propagations of a single contour/segmentation, an entire ROI/segment, consisting of any number of contour(s)/segment(s) may be copied, or an entire ROI Collection, consisting of any number of ROI(s)/segment(s) containing any number of contour(s)/segmentation(s).  The behaviour entirely depends on the user-defined configuration parameters and relationships between the two image domains.

The code used to copy or propagate an ROI Collection for one image session to another relies heavily on the use of [SimpleITK](https://simpleitk.org/) [1-3].


[1]: R. Beare, B. C. Lowekamp, Z. Yaniv, “Image Segmentation, Registration and Characterization in R with SimpleITK”, J Stat Softw, 86(8), https://doi.org/10.18637/jss.v086.i08, 2018.

[2]: Z. Yaniv, B. C. Lowekamp, H. J. Johnson, R. Beare, “SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research”, J Digit Imaging., https://doi.org/10.1007/s10278-017-0037-8, 31(3): 290-303, 2018.

[3]: B. C. Lowekamp, D. T. Chen, L. Ibáñez, D. Blezek, “The Design of SimpleITK”, Front. Neuroinform., 7:45. https://doi.org/10.3389/fninf.2013.00045, 2013.
