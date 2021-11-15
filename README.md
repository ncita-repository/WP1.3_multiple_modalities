# NCITA Repository Unit Work Package 1.3: XNAT-OHIF viewer - multiple modalities
Specific repositories and projects of the form WP1.x-y are for specific second-level, viewer-related work packages in the NCITA Repository Unit strategy.

The aim of this work package of the NCITA Repository Unit will be to develop approaches within the OHIF viewer and XNAT to work with images from connected studies. The initial focus will be on image registration and the copying of regions-of-interest between:
* different DICOM series within the same imaging session (e.g., draw ROI on DWI and see where it is on the corresponding T1w image);
* corresponding images from different time points for the same patient (registration required and we will investigate rigid and non-rigid);
* corresponding images from different modalities (start with the simple case of inherently co-registered PET-CT and work towards more difficult scenarios).

The NCITA Repository Unit site as a whole is at an early stage of the life cycle and does not currently have a formal structure of documentation releases. Please bear with us: things will become more organised as the NCITA project gets into gear! Note that the ICR also releases code related to the XNAT-OHIF image viewer at https://bitbucket.org/icrimaginginformatics/. 

## Tools covered in this repository

This repository contains two tools:

1. *ROI copy/propagation* tool (*app.py*)

2. *XNAT snapshot* tool (*snapshots.py*)

All code is written in Python 3.  

The definitions at the end of this page will be useful in understanding the terminology used within this README.  In brief, an *entity-of-interest* (or simply *entity*) may refer to a single *contour*, a single *segmentation*, a collection of *contour*s (= *ROI*) or *segmentation*s (= *segment*), or an entire *ROI Collection* (consisting of any number of *contour*s or *segmentation*s within any number of *ROI*s or *segment*s).

The "ROI" in *ROI copy/propagation* refers to a region-of-interest in the broadest sense, which may refer to any *entity* within a *DICOM-RTSTRUCT* or *DICOM-SEG* *ROI Collection*, and is not to be confused with a collection of contours within an *RTSTRUCT*-based ROI Collection.   

A *copy* will refer to an operation that does not preserve the spatial coordinates of the entity, whereas a *propagation* implies that the propagated entity shares the same spatial coordinates as the original.  Hence a contour that segments the ventricle in the brain in a *source* DICOM slice will not necessarily coincide with the ventricle in a *target* slice upon making the *copy*, whereas it will be expected to overlay when performing a *propagation*.

# *ROI copy/propagation* tool (*app.py*)

The purpose of this tool is to make a "copy" of, or to "propagate" an entity from a *source* ROI Collection so that the entity overlays onto a *target* DICOM series.  The copy may be a "non-relationship-preserving copy" (think copy and paste function) or a "relationship-preserving propagation" of the entity.

Upon the successful copying or propagation of an entity, a new ROI Collection will be created and exported to *src/outputs/roicols/*.  The new ROI Collection will also be automatically uploaded to XNAT as an ROI Collection for the *target* DICOM series.

If image registration is required to propagate the entity, a search will be done on XNAT for a suitable DRO to use instead (more on this later).  If one is not found, image registration will be performed and a new *DICOM Registration Object (DRO)* will be created and exported to *src/outputs/dros/*.  

A rigid registration will yield a *Spatial Registration Object (SRO)*, while a deformable one will result in the creation of a *Deformable Spatial Registration Object (DSRO)*.  As well as the exported file, the DRO will automatically be uploaded to XNAT as a DICOM resource at the subject level.  A sample SRO and DSRO (stored in *src/inputs/sample_dros/*) are used as templates to generate the new DRO. 

Once a suitable DRO exists on XNAT, future propagation calls that involve the same *source* and *target* DICOM series and the same registration type (i.e. "rigid", "affine" or "bspline"), will result in the use of the transformation matrix stored within the DRO, thus by-passing the computationally expensive registration process.  While a run using image registration might take over 200 s, use of a DRO will reduce the execution time to ~30 s for the same datasets (of which ~10 s is spent performing the search).  

At present all subject DICOM resources are parsed for a match.  Since this is inefficient, in future work the DRO will be stored as a subject assessor with an accompanying XML file so that matches can be made without needing to parse all resource files.  This should shave on the order of 10 s for parsing of ~30 DICOM resources, for example.

The tool utilises XNAT Alias Tokens to avoid the creation of multiple user sessions.  When running *app.py* a search for an alias token will be made in *src/xnat_tokens/*. If a token is found and if the token is valid, it will be used, avoiding the need to enter a password. If a token does not exist, or if the token is invalid, the user will be prompted for a password, and the token file will be overritten for future use.

Running of the tool requires two JSON files:  a JSON containing global variables and a JSON containing XNAT configuration parameters.  The reason for splitting the variables in this way was to differentiate between variables that the user is not expected to need to modify readily and those that will be specific to each use of the tool.  In future work, the global variables will likely be stored in a (not-yet-existing) XNAT container, whilst the other variables will be provided by the XNAT Container Service.  **Note: Any parameters defined in *xnatCfg.json* will override those defined in *global_variables.json*.**

## Global variables (*global_variables.json*)

The global variables are stored in *src/global_variables.json* as a dictionary.  The parameters relate to specifics on how image resampling and registration is to be carried out, and other default parameters.  The module used to generate *global_variables.json* is *src/create_global_file.py*.  See the **Definitions** section for a listing of all parameters.

## XNAT configuration file (*xnatCfg.json*)

The XNAT config file contains the metadata that identifies the *source* and *target* DICOM series ("scans" in XNAT parlance), *source* ROI Collection, and if applicable, *target* ROI Collection.  See the **Definitions** section for a listing of all parameters.  

Based on the metadata contained in *xnatCfg.json*, the tool will use XNAT REST API calls to fetch the required data, download it to *src/xnat_downloads/*, import the data and depending on the relationship between the *source* and *target* DICOM series and other metadata, will either perform a *non-relationship-preserving copy* or *relationship-preserving propagation* of the entity of interest.

By default the XNAT config JSON file has the default filename *xnatCfg*, but the user can define an alternative file name, e.g. "xnatCfg_241js23", is so that when *app.py* is run with the optional input argument *--xnatCfgFname=cfgDict_241js23*, that is the file that will be imported.  

This will allow for the possibility of concurrent calls to *app.py* with distinct XNAT config files. This feature was added in anticipation of the need for a scalable solution within XNAT. In future work, the XNAT config file will be provided by the XNAT Container Service.  Since a container does not yet exist for this tool, a work-around does, and is covered next.


## XNAT configuration files in *src/xnat_configs/* and *select_xnat_config.py* (i.e. "work-around")

The user creates any number of run-specific XNAT configuration files stored in *src/xnat_configs/*. The contents of the run-specific XNAT config file, e.g. *src/xnat_configs/runID.json*, is copied to *xnatCfg.json* (or a custom file name). That way when the tool is run *xnatCfg.json* (or the custom file name) is imported along with *global_variables.json*.

The module *config.py* is used to generate configuration files.  Within it is the function *create_config_files()*, which contains pre-populated data specific to the XNAT that was used to develop the code. Hence the contents will only be relevant to you if your XNAT contains the same data and same XNAT metadata.  You can use the pre-populated code as a template for your own XNAT.  Note that the `runID` key in the dictionary is expected to match the file name of the JSON.

Once there exists at least one XNAT configuration file to be run, the next step is to select the desired file and copy the dictionary to *xnatCfg.json* (or other custom file name). This is done using *select_xnat_config.py*.

The commands used to execute the above steps are covered in the section **Using the tools**.


# XNAT "snapshot" tool (*snapshots.py*)

This is a basic tool that fetches high-level metadata from an XNAT and produces a "snapshot", exported as an XLSX file to *src/xnat_snapshots/*, and is a stand-alone feature completely separate from the business of copying/propagating ROIs.  However this tool does share some common modules used for ROI copying/propagating, hence its inclusion here.

At present two XLSX files are exported - one that provides a snapshot broken down by projects, and one that is XNAT-wide.  The former includes the principle investigator, project description, a list of investigators, number of subjects, experiments, MR/CT/PET sessions, etc.  The XNAT-wide table includes totals and averages per project of subjects, experiments, image sessions, DICOM image files, etc.  Examples can be found in *src/xnat_snapshots_examples/*.

![An example of an XNAT-by-project snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_by_project_snapshot.png "An example of an XNAT-by-project snapshot")

![An example of an XNAT-wide snapshot](https://github.com/ncita-repository/WP1.3_multiple_modalities/blob/master/src/xnat_snapshots_examples/example_XNAT_wide_snapshot.png "An example of an XNAT-wide snapshot")

# Using the tools

## Install required packages

After cloning the repository, *pip* install the required packages listed in the *requirements.txt* file:

    pip3 install -r requirements.txt


## Running the *ROI copy/propagation* tool (*app.py*)

A *source ROI Collection* (DICOM-RTSTRUCT or DICOM-SEG) is copied (or propagated) to a *target* DICOM series.  One such single operation will be referred to as a *run*.  As described above, running of the tool requires there to be an XNAT config JSON file (with default file name *xnatCfg.json*) and a global variables JSON file (*global_variables.json*).  Hence a bit of prep work is required before trying to run the code:

0. (Optional) Modify *global_variables.json*.  Alternatively set non-default parameters in the XNAT config file, since the parameters stored there take priority over those in *global_variables.json*.

1. Create a configuration file(s) for the run(s) you wish to carry out. This only has to be done once, unless you wish to add more XNAT config files or modify existing ones.

2. Copy the XNAT config file to be run from *src/xnat_configs/* to *src/xnatCfg.json*. 

3. Run the tool.

All commands below are to be made from the *src* directory.


### 0. (Optional) Modify *global_variables.json*

Modify *global_var_file.py* as desired and save changes.  Then in a command shell run:

	python global_var_file.py

This will refresh *global_variables.json*

(Optional) If running in a Python terminal:

	from global_var_file import create_global_vars
	create_global_vars()


### 1. Creating a XNAT configuration file(s)

There are a few ways to go about creating an XNAT configuration file(s):

#### Option 1 - Copy and modify an existing XNAT configuration file

This is by far the quickest way to go about implementing your first run of the code.  

Simply copy an existing XNAT config file in */src/xnat_configs/*, rename the file, open it in your preferred text editor and modify the dictionary's values accordingly.  Note that the the name you assign to the file must agree with the value given for the `runID` key.  

#### Option 2 - Define a new dictionary within the function *create_xnat_config_files()* in *xnat_config_files.py*

This option is more suitable if you wish to run the code on multiple datasets and you don't want to manually create files one at a time.  

Within the function *create_xnat_config_files()* there are cases where a dictionary has been fully defined (by scratch) - see for example the definition of the dictionary with key *NCITA_test_RR1* (lines 97-116).  In this case a `runID` is assigned (line 97), and a dictionary is appended to the dictionary `cfg` with `runID`.  

Further down a new dictionary is added to `cfg` using the key *NCITA_test_RR2* (line 120), and the previously defined dictionary is copied (line 121).  Then the `runID` (line 122) and other selected keys are modified (lines 123-127).

Use a previously fully-defined dictionary as a template for new dictionaries as you go along.  Remember to re-run the command *create_xnat_config_files()* to update your list of XNAT config files.

Once the dictionaries have been defined in the script, the XNAT config files can be created by running the following in a command shell:

	python xnat_config_files.py

This will generate and/or overwrite a list of JSONs in *src/xnat_configs/* - one file for each key-dictionary pair within the dictionary `cfg`.

(Optional) If running in a Python terminal:

	from xnat_config_files import create_xnat_config_files
	create_xnat_config_files()

### 2. Copy the XNAT config file to be run from *src/xnat_configs/* to *src/xnatCfg.json*

In a command shell run:

	python select_xnat_config.py runID

where `runID` is the file name of the XNAT config file in *src/xnat_configs/* you wish to run.

This will copy the dictionary in *src/xnat_configs/runID.json* to *src/xnatCfg.json*.

(Optional) If you wish to assign the XNAT config file a custom name (e.g. to allow for concurrent runs of *app.py* with distinct XNAT config files), use the optional input parameter `xnatCfgFname`:

	python select_xnat_config.py runID --xnatCfgFname=custom_fname

(Optional) If running in a Python terminal:

	from select_xnat_config import create_xnat_config_files
	create_xnat_config_files("runID")` or `create_xnat_config_files("runID", "custom_fname")
	

### 3. Run the tool

In a command shell run:

	python app.py

(Optional) or if running the tool on an XNAT config file with custom file name:

	python app.py --xnatCfgFname=custom_fname

(Optional) If running in a Python terminal:

	from app import main
	main()` or `main("custom_fname")

You should be prompted to enter an XNAT password.  If the XNAT `url`, `username` and password you just entered were correct, an XNAT session will be established and an XNAT alias token generated.  The alias token will be saved to *src/xnat_tokens/*.  If the code is executed again using the same `url` and authentication credentials, and within the lifetime of the existing alias token, the token will be used, avoiding the need to re-enter your password, and the creation of another XNAT user session.


# XNAT Snapshots (*snapshots.py*)

The tool can either be run by executing the following in a command shell:

	python snapshots.py configs http://10.1.1.20 --username --password
	
(Optional) In a Python terminal run:

	from snapshots import get_xnat_snapshot
	get_xnat_snapshot("http://10.1.1.20")` or `get_xnat_snapshot("http://10.1.1.20", "username", "password")

The optional argument `username` can be provided without also providing `password` (the user will be prompted to enter a password).

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

## The minimal *key:value* pairs that define *xnatCfg.json*

- `runID`         : An identifier for the run
- `url`           : XNAT address
- `username`      : XNAT username
- `proID`         : XNAT Project ID
- `subjLab`       : XNAT Subject Label
- `srcExpLab`     : XNAT Experiment Label for the *source* dataset
- `srcScanID`     : XNAT Scan ID for the *source* dataset
- `srcSlcNum`     : Slice number of the *source* ROI Collection to be copied. Leave as `None` if all contours/segmentations are to be copied.
- `srcRoicolName` : Name of the *source* ROI Collection
- `srcRoiName`    : Name of the ROI/segment within `srcRoicolName` to be copied. Leave as `None` if all ROIs/segments are to be copied.
- `roicolMod`     : Modality of the *source* ROI Collection
- `trgExpLab`     : Target dataset's XNAT Experiment Label for the *target* dataset
- `trgScanID`     : XNAT Scan ID for the *target* dataset
- `trgSlcNum`     : Slice number where the contour/segmentation is to be copied to for the *target* dataset. Leave as `None` if all contours/segmentations are to be copied.
- `trgRoicolName` : Name of the *target* ROI Collection. Leave as `None` if one does not exist (i.e. one does not wish to copy a specific contour/segmentation within the Source ROI Collection to an existing Target ROI Collection). Such an operation is only valid for "direct" copies (see below for further details).
- `trgRoiName`    : Name of the ROI/segment within `trgRoicolName` where a contour/segmentation is to be copied to. Leave as `None` if all ROIs/segments are to be copied.
- `srcFidsFname`  : File name to the .txt containing fiducials for the Source DICOM scan.
- `trgFidsFname`  : File name to the .txt containing fiducials for the Target DICOM scan.

The above list is of the *minimal* key:value pairs expected in *xnatCfg.json*.  The dictionary may also contain any non-default variable defined in *global_variables.json*, since the value in *xnatCfg.json* will take priority over that in *global_variables.json*.  Hence the global variables will not necessarily need to be changed, as run-specific changes can be made in *xnatCfg.json*.

For example, many runs will not require the use of fiducials, hence `srcFidsFname` and `trgFidsFname` will both be `""`, and by default `initMethod` will be `"geometry"`. If fiducials are to be used, then `srcFidsFname` and `trgFidsFname` will not be empty strings, and `initMethod` should be `"landmarks"`.

## Fiducials


# Credits

The code used to copy or propagate an entity for one image session to another relies heavily on the use of [SimpleITK](https://simpleitk.org/) [1-3].


[1]: R. Beare, B. C. Lowekamp, Z. Yaniv, “Image Segmentation, Registration and Characterization in R with SimpleITK”, J Stat Softw, 86(8), https://doi.org/10.18637/jss.v086.i08, 2018.

[2]: Z. Yaniv, B. C. Lowekamp, H. J. Johnson, R. Beare, “SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research”, J Digit Imaging., https://doi.org/10.1007/s10278-017-0037-8, 31(3): 290-303, 2018.

[3]: B. C. Lowekamp, D. T. Chen, L. Ibáñez, D. Blezek, “The Design of SimpleITK”, Front. Neuroinform., 7:45. https://doi.org/10.3389/fninf.2013.00045, 2013.
