# NCITA Repository Unit Work Package 1.3: XNAT-OHIF viewer - multiple modalities
Specific repositories and projects of the form WP1.x-y are for specific second-level, viewer-related work packages in the NCITA Repository Unit strategy.

The aim of this work package of the NCITA Repository Unit will be to develop approaches within the OHIF viewer and XNAT to work with images from connected studies. The initial focus will be on image registration and the copying of regions-of-interest between:
* different DICOM series within the same imaging session (e.g., draw ROI on DWI and see where it is on the corresponding T1w image);
* corresponding images from different time points for the same patient (registration required and we will investigate rigid and non-rigid);
* corresponding images from different modalities (start with the simple case of inherently co-registered PET-CT and work towards more difficult scenarios).

The NCITA Repository Unit site as a whole is at an early stage of the life cycle and does not currently have a formal structure of documentation releases. Please bear with us: things will become more organised as the NCITA project gets into gear! Note that the ICR also releases code related to the XNAT-OHIF image viewer at https://bitbucket.org/icrimaginginformatics/. 

# Tools covered in this repository

This repository contains two tools:

1. A *ROI copy/propagation* tool

2. An XNAT "snapshot" tool

All code is written in Python 3.  

The *ROI* in *ROI copy/propagation* refers to a region-of-interest in the broadest sense, which may refer to any "entity-of-interest", or simply "entity", within a *DICOM-RTSTRUCT* or *DICOM-SEG* *ROI Collection*, and is not to be confused with a collection of contours within an *RTSTRUCT*-based ROI Collection.  The definitions at the end of this page will be useful in understanding the terminology used within this *README*.  In brief, an *entity* may refer to single *contour*, a single *segmentation*, a collection of *contour*s (= *ROI*) or *segmentation*s (= *segment*), or an entire *ROI Collection* (consisting of any number of *contour*s or *segmentation*s within any number of *ROI*s or *segment*s).

A *copy* will refer to an operation that does not preserve the spatial coordinates of the entity, whereas a *propagation* implies that the propagated entity shares the same spatial coordinates as the original.  Hence a contour that segments the ventricle in the brain in a *source* DICOM slice will not necessarily coincide with the ventricle in a *target* slice upon making the *copy*, whereas it will be expected to overlay when performing a *propagation*.

## *ROI copy/propagation* tool

The purpose of this tool is to make a "copy" of, or to "propagate" an entity from a *source* ROI Collection so that the entity overlays onto a *target* DICOM series.  The copy may be a "non-relationship-preserving copy" (think copy and paste function) or a "relationship-preserving propagation" of the entity.

This tool requires two JSON files to run:  a JSON containing global variables and a JSON containing XNAT configuration parameters.  The global variables are stored *src/global_variables.json* as a dictionary.  The parameters relate to specifics on how image resampling and registration is carried out, as well as default directories amongst other things.  The module used to generate *global_variables.json* is *src/create_global_file.py*.  See the *Definitions* section for a listing of all parameters.


The XNAT config file contains the metadata that identifies the *source* and *target* DICOM series ("scans" in XNAT parlance), *source* ROI Collection, and if applicable, *target* ROI Collection.  See the *Definitions* section for a listing of all parameters.  Based on the user-inputed data, the tool will use XNAT REST API calls to fetch the required data, download it to *src/xnat_downloads* (by default), import the data and depending on the relationship between the *source* and *target* DICOM series and other metadata, will either perform a *non-relationship-preserving copy* or *relationship-preserving propagation* of the entity of interest.

By default the XNAT config JSON file has the filename *xnatCfg*, but the user can define an alternative file name. The configuration file that is imported upon executing the code is *src/cfgDict.json*.  Ultimately this file will be provided by the XNAT Container Service.  Since a container does not yet exist for this tool, the following work-around exists:  The user selects an XNAT config file (e.g. *src/xnat_configs/runID.json*), and the file is copied to *src/xnatCfg.json* (file name by default). That way when the tool is run *xnatCfg.json* is imported along with *global_variables.json*.

The reason for splitting the variables in this way was to differentiate between variables that the user is not expected to need to modify readily and those that will. The global variables will likely be stored in a (not-yet-existing) XNAT container, whilst the other variables will be provided by the XNAT Container Service.

The option to assign a different file name (other than "xnatCfg.json"), e.g. "xnatCfg_241js23.json", is so that when *app.py* is run with the optional input argument *--xnatCfgFname=cfgDict_241js23*, the desired XNAT config file will be imported, allowing for the possibility of concurrent calls to *app.py* with unique XNAT config files. This feature was added in anticipation of the need for a scalable solution within XNAT.

Upon the successful copying or propagation of an entity, a new ROI Collection will be created and exported (by default) to the *src/outputs/new_roicols* directory.  The new ROI Collection will also be automatically uploaded to XNAT using the same credentials used to fetch the data from XNAT.

If an image registration is required to propagate the entity, a search will be done on XNAT for a suitable DRO to use instead (more on this later).  If one is not found, image registration will be performed and a new *DICOM Registration Object (DRO)* will be created and exported (by default) to *src/outputs/new_DRO*.  A rigid registration will yield a *Spatial Registration Object (SRO)*, while a deformable one will result in the creation of a *Deformable Spatial Registration Object (DSRO)*.  As well as the exported file, the DRO will automatically be uploaded to XNAT as a DICOM resource at the subject level.  A sample SRO and DSRO, stored in *src/inputs/sample_DROs*, are used as templates to generate the new DRO. 

Once a suitable DRO exists on XNAT, future propagation calls that involve the same *source* and *target* DICOM series and the same registration type (i.e. "rigid", "affine" or "bspline"), will result in the use of the transformation matrix stored within the DRO, thus by-passing the computationally expensive registration process.  While a run using image registration might take over 200 s, use of a DRO will reduce the execution time to 30 s, for example.  

At present all subject DICOM resources are parsed and scanned for a match for the run in progress.  Since this is not an efficient way of doing things in future work the DRO will be stored as a subject assessor with a corresponding XML file so that matches can be made without needing to parse all resource files.  This should shave on the order of 10 s for parsing of ~30 DICOM resources, for example.

The algorithm makes use of XNAT Alias Tokens to avoid the creation of multiple user sessions.  The alias token is stored in the *src/tokens* directory.

## XNAT "snapshot" tool

This is a basic tool that fetches high-level metadata from an XNAT and produces a "snapshot", exported as an XLSX file to *src/xnat_snapshots*, and is a stand-alone feature completely separate from the business of copying/propagating ROIs.  However this tool does share some common modules used for *ROI* copying/propagating, hence its inclusion here.

At present two XLSX files are exported - one that provides a snapshot broken down by projects, and one that is XNAT-wide.  The former includes the principle investigator, project description, a list of investigators, number of subjects, experiments, MR/CT/PET sessions, etc.  The XNAT-wide table includes totals and averages per project of subjects, experiments, image sessions, DICOM image files, etc.  Examples can be found in *src/xnat_snapshots_examples*.

![An example of an XNAT-by-project snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_by_project_snapshot.png "An example of an XNAT-by-project snapshot")

![An example of an XNAT-wide snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_wide_snapshot.png "An example of an XNAT-wide snapshot")

# Using the tools

## Install required packages

After cloning the repository, *pip* install the required packages listed in the *requirements.txt* file:

    pip3 install -r requirements.txt


## Running the *ROI copy/propagation* tool

A *source ROI Collection* (DICOM-RTSTRUCT or DICOM-SEG) is copied (or propagated) to a *target* DICOM series.  One such single operation will be referred to as a *run*.  As described above, running of the tool requires there to be a config JSON file (with default file name *cfgDict.json*).  Hence a bit of prep work is required before trying to run the code:

0. (Optional) Modify *global_variables.json*.

1. Create a configuration file(s) for the run(s) you wish to carry out. This only has to be done once, unless you wish to add more XNAT config files or modify existing ones.

2. Copy the XNAT config file to be run from *src/xnat_configs/{runID}.json* to *src/xnatCfg.json*. 

3. Run the code on *src/xnatCfg.json*.


### 0. Modify *global_variables.json*

### 1. Creating a configuration file(s)

The package *config.py* is used to generate configuration files.  Within the *config.py* module is the function *create_config_files()*.  This contains pre-populated data specific to the XNAT that was used to develop the code. You can use the pre-populated code as a template for your own XNAT.

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


## Running the "OHIF-integration" version

This tool will not be used directly by an end-user but details are provided for posterity.

The choice of which variables to treat as global has been determined in *global_var_file_ohif.py*.  If changes are to be made, the easiest way to do this is:

1. Modify *global_var_file_ohif.py* and save the changes

2. Run in a command shell:

	`python global_var_file_ohif.py`

or in a Python terminal: 

	`from global_var_file_ohif import create_global_var_file`
	`create_global_var_file()`

This will update *global_variables.json* in *src/*.

Next create a dictionary of dictionaries in *paths_from_ohif.py* that contains the required metadata that fully defines the datasets and their file or directory paths, the resampling/registration parameters, output directory names, etc. As before:

1. Modify *paths_from_ohif.py* and save the changes

2. Run in a command shell:

	`python paths_from_ohif.py`

or in a Python terminal:

	`from paths_from_ohif import create_cfgDict_json`
	`create_cfgDict_json()`

This will update *cfgDict_ohif.json*.  Ensure that the paths defined in the JSON link to data in *src/inputs/scans*, *src/inputs/roicols*, and *src/inputs/dro* and *src/inputs/fiducials* if applicable.


Next run *app_ohif.py* providing as input arguments the location of *cfgDict_ohif.json* (which should be in the *src/* directory) and the key (= `runID`) in *cfgDict_ohif.json* to be run, i.e. run in a command shell:

	`python app_ohif.py . RR3_contour`

or in a Python terminal:

	`import os`
	`from app_ohif import main`
	
	`main(os.getcwd(), 'RR3_contour')`

It should be re-iterated that this is not a tool to be used for end-users, but merely demonstrates that if the metadata, paths to data and other variables are passed in place of *global_variables.json* and *cfgDict_ohif.json*, the code should run as expected.  


# XNAT Snapshots

The tool can either be run as a module in a command shell or as a Python function.

#### Option 1 - Run *snapshots.py* as a script in a command shell

In a command shell, run the Python module providing the XNAT url:

	python snapshots.py configs http://10.1.1.20
	
#### Option 2 - Run the function *get_xnat_snapshot()* in a Python shell

In a Python shell, run *get_xnat_snapshot()* providing the XNAT url as an input argument:

	from snapshots import get_xnat_snapshot
	get_xnat_snapshot("http://10.1.1.20")



# Definitions

## 3D images

Every DICOM series (or "scan") consists of any number of 2D pixel arrays - i.e. 2D images.  The entire series or scan make up a 3D image.  When referring to 3D images, this will usually imply the DICOM series/scan.  3D images also crop up when referring to 3D "label images".  These are a collection of binary masks (i.e. segmentations) that occupy the same extent, voxel resolution, origin and direction as the 3D image representation of the DICOM series.

## Entities

The term *entity* is used to refer to a region-of-interest in the most general sense, and may refer to single contour, a collection of contours (= "ROI"), or an entire ROI Collection consisting of any number of contours within any number of ROIs, all originating from a DICOM-RTSTRUCT ROI Collection.  Likewise, *entity* may also refer to a single segmentation, a collection of segmentations (= "segment"), or an entire ROI Collection consisting of any number of segmentations within any number of segments, all originating from a DICOM-SEG ROI Collection.
 
## Relationship-preserving propagations v Non-relationship-preserving copies

A *non-relationship preserving copy* (sometimes referred to as a "direct" copy) is analogous to a copy-and-paste function - e.g. "Copy the *source* contour/segmentation corresponding to slice `srcSlcNum` in ROI/segment `srcRoiName` in ROI Collection `srcRoicol` to a contour/segmentation that will overlay onto *target* slice `trgSlcNum`.  "Direct" copies do not maintain spatial relationships between the *source* contour/segmentation and the *target* contour/segmentation's location in space.  One might want to perform a "direct" copy to use an existing contour/segmentation as a template for a new one that can sculpted to match the anatomical features on the *target* scan.  By contrast, a *relationship-preserving propagation* does respect spatial relationships and hence the contour(s)/segmentation(s) end up where they "should anatomically".

*Copies* are made between two 3D images that have the same voxel resolution (*Pixel Spacing* and difference between adjacent *ImagePositionPatient* values along the scan direction), same patient position (*ImagePositionPatient*), orientation (*ImageOrientationPatient*) and frame-of-reference (*FrameOfReferenceUID*).  Propagations require a bit more work - either a resampling is required (if the two 3D images do not have the same voxel spacings, or same extent, but have the same FOR), or image registration (if they do not have the same FOR).  

A *non-relationship preserving copy* always involves the copying of a single contour/segmentation to a single contour/segmentation.  If the *source* and *target* 3D images have different voxel resolutions the copied entity will be interpolated (i.e. image resampling).  If resampling of the *source* contour/segmentation led to multiple contours/segmentations, the *ROI/segment* is collapsed to a single *contour/segmentation* to maintain the expected behaviour.  On the contrary, a *relationship-preserving propagation* are broader, in that a single contour/segmentation may map to multiple contours/segmentations (depending on the relative voxel spacings and FOR of the two 3D images). 

When making a *relationship-preserving propagation*, in addition to propagations of a single contour/segmentation, an entire ROI/segment, consisting of any number of contour(s)/segment(s) may be propagated, or an entire ROI Collection, consisting of any number of ROI(s)/segment(s) containing any number of contour(s)/segmentation(s).  The behaviour entirely depends on the user-defined configuration parameters and relationships between the two 3D image domains.

## The *key:default value* pairs that define *global_variables.json*

- `forceReg` : Set to `True` if you wish to apply image registration on a run for which the source and target scans have the same *FrameOfReferenceUID* (image registration would not normally be used in this case).  The default value is `False`.
- `useDroForTx` : If `True` and a suitable DRO is found on XNAT, the transformation parameters stored in the DRO will be used, by-passing image registration. The default value is `True`.
- `regTxName` : Defines the transformation to apply during registration. Acceptable values include `"rigid"`, `"affine"` and `"bspline"`. The default value is `"affine"`.
- `initMethod` : Defines how a registration is initialised.  Acceptable values include `"geometry"` (aligning to image centres), `"moments"` (align to centre-of-mass of image intensity), and `"landmarks"` (`srcFidsFpath` and `trgFidsFpath` must be defined, otherwise `"geometry"` will be used). The default value is `"geometry"`.
- `maxIters` : The maximum number of iterations to allow during optimisation for image registation. The default value is `512`.
- `applyPreResBlur` : If `True` the source label image will be Gaussian blurred prior to resampling or registration. Doing so helps to avoid aliasing effects. The default value is `False` (since pre-registration Gaussian blurring is performed when `resInterp = "BlurThenLinear"` (default value)).
- `preResVar` : The variance of the Gaussian blur that will be applied if `applyPreResBlur` is `True`. The default value is `(1, 1, 1)`.
- `resInterp` : Defines the interpolation to use when resampling label images.  Acceptable values include `"NearestNeighbor"`, `"LabelGaussian"` and `"BlurThenLinear"`.  `"NearestNeighbor"` and `"LabelGaussian"` are binary interpolators, whereas `"BlurThenLinear"` applies a Gaussian blur (defined by `preResVar`), followed by linear resampling, then binary thresholding to retain a binary label image. This option helps to avoid aliasing effects.
- `applyPostResBlur` : If `True` the resampled label image will be Gaussian blurred. The default value is `True`.
- `postResVar` : The variance of the Gaussian blur that will be applied if `applyPostResBlur` is `True`. The default value is `(1, 1, 1)`.
- `whichSrcRoicol` : There may be multiple ROI Collection that match the source metadata defined in the XNAT config file. If multiple hits occur, this parameter determines how to proceed.  Acceptable inputs include `"oldest"` (select the oldest file), `"newest"` (select the newest file) and `"user"` (present a list of all matching files and prompt the user for a selection). The default value is `"oldest"`.
- `addToRoicolLab`: This optional parameter allows for text to be added to the *StructureSetLabel* (for DICOM-RTSTRUCT) or *SeriesDescription* (DICOM-SEG). The default value is `""`.
- `p2c` : If `True` more verbose output will be printed to the console. The default value is `False`.
- `cwd` : The current working directory during the time *global_variables.json* was created.
- `xnatCfgDir` : The name of the directory containing XNAT configuration files (relative to *src/*). The default value is `"xnat_configs"`.
- `inputsDir` : The name of the directory where data files are to be imported from (relative to *src/*). The default value is `"inputs"`.
- `outputsDir` : The name of the directory where data files are to be exported to (relative to *src/*). The default value is `"outputs"`.
- `sampleDroDir` : The name of the directory where sample DRO files are to be imported from (relative to *src/*). The default value is `{inputsDir}"/sample_dros"`.
- `fidsDir` : The name of the directory where fiducials files (*TXT*) are to be imported from (relative to *src/*). The default value is `{inputsDir}"/fiducials"`.
- `rtsExportDir` : The name of the directory where DICOM-RTSTRUCT ROI Collection files (*DCM*) are to be exported to (relative to *src/*). The default value is `{outputsDir}"/roicols"`.
- `segExportDir` : The name of the directory where DICOM-SEG ROI Collection files (*DCM*)are to be exported to (relative to *src/*). The default value is `{outputsDir}"/roicols"`.
- `droExportDir` : The name of the directory where DICOM-DRO files (*DCM*)are to be exported to (relative to *src/*). The default value is `{outputsDir}"/dros"`.
- `txExportDir` : The name of the directory where (*SimpleITK Transform*) registration transforms (*TFM* and *HDF*) files are to be exported to (relative to *src/*). The default value is `{outputsDir}"/transforms"`.
- `imExportDir` : The name of the directory where (*SimpleITK Image*) source and target images (*HDF*) files are to be exported to (relative to *src/*). The default value is `{outputsDir}"/images"`.
- `labimExportDir` : The name of the directory where (*SimpleITK Image*) source, target and resampled/registered binary label (if applicable) images (*HDF*) files are to be exported to (relative to *src/*). The default value is `{outputsDir}"/label_images"`.
- `logsExportDir` : The name of the directory where log files (*TXT*) files are to be exported to (relative to *src/*). The default value is `{outputsDir}"/logs"`.
- `rtsPlotsExportDir` : The name of the directory where figures of DICOM-RTSTRUCT copy/propgation results overlaid onto DICOM images (*JPG*) are to be exported to (relative to *src/*). The default value is `{outputsDir}"/plots_rts"`.
- `segPlotsExportDir` : The name of the directory where figures of DICOM-SEG copy/propgation results overlaid onto DICOM images (*JPG*) are to be exported to (relative to *src/*). The default value is `{outputsDir}"/plots_seg"`.
- `resPlotsExportDir` : The name of the directory where figures of resampling/registering (if applicable) results (*JPG*) are to be exported to (relative to *src/*). The default value is `{outputsDir}"/plots_res"`.
- `exportRoicol` : If `True` a new ROI Collection will be exported to `rtsExportDir` (for DICOM-RTSTRUCTs) or `segExportDir` (for DICOM-SEGs). The default value is `True`.
- `exportDro` : If `True` a new DICOM-DRO will be exported to `droExportDir`. The default value is `True`.
- `exportTx` : If `True` the *SimpleITK Transform* will be exported to `txExportDir`. The default value is `False`.
- `exportIm` : If `True` the *SimpleITK Image*s for source and target will be exported to `imExportDir`. The default value is `False`.
- `exportLabim` : If `True` the *SimpleITK Image*s for source, target and resampled/registered (if applicable) binary label images will be exported to `labimExportDir`. The default value is `False`.
- `exportPlots` : If `True` plots will be exported to `rtsPlotsExportDir` or `segPlotsExportDir`. The default value is `False`.
- `exportLogs` : If `True` logs will be exported to `logsExportDir`. The default value is `False`.
- `uploadDro` : If `True` the new DICOM-DRO will be uploaded to XNAT. The default value is `True`.
- `overwriteDro` : If `True` an existing DICOM-DRO on XNAT will be overwritten. The default value is `False`.


# Credits

The code used to copy or propagate an entity for one image session to another relies heavily on the use of [SimpleITK](https://simpleitk.org/) [1-3].


[1]: R. Beare, B. C. Lowekamp, Z. Yaniv, “Image Segmentation, Registration and Characterization in R with SimpleITK”, J Stat Softw, 86(8), https://doi.org/10.18637/jss.v086.i08, 2018.

[2]: Z. Yaniv, B. C. Lowekamp, H. J. Johnson, R. Beare, “SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research”, J Digit Imaging., https://doi.org/10.1007/s10278-017-0037-8, 31(3): 290-303, 2018.

[3]: B. C. Lowekamp, D. T. Chen, L. Ibáñez, D. Blezek, “The Design of SimpleITK”, Front. Neuroinform., 7:45. https://doi.org/10.3389/fninf.2013.00045, 2013.
