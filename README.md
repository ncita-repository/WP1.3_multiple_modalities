# WP1.3_multiple_modalities

### NCITA Repository Unit Work Package 1.3: XNAT-OHIF viewer - multiple modalities
Specific repositories and projects of the form WP1.x-y are for specific second-level, viewer-related work packages in the NCITA Repository Unit strategy.

The aim of this work package of the NCITA Repository Unit will be to develop approaches within the OHIF viewer and XNAT to work with images from connected studies. The initial focus will be on image registration and the copying of regions-of-interest between:
* different DICOM series within the same imaging session (e.g., draw ROI on DWI and see where it is on the corresponding T1w image);
* corresponding images from different time points for the same patient (registration required and we will investigate rigid and non-rigid);
* corresponding images from different modalities (start with the simple case of inherently co-registered PET-CT and work towards more difficult scenarios).

The NCITA Repository Unit site as a whole is at an early stage of the life cycle and does not currently have a formal structure of documentation releases. Please bear with us: things will become more organised as the NCITA project gets into gear! Note that the ICR also releases code related to the XNAT-OHIF image viewer at https://bitbucket.org/icrimaginginformatics/. 

The code used to copy or propagate an ROI Collection for one image session to another relies heavily on the use of SimpleITK [1-3].


References:

[1] R. Beare, B. C. Lowekamp, Z. Yaniv, “Image Segmentation, Registration and Characterization in R with SimpleITK”, J Stat Softw, 86(8), https://doi.org/10.18637/jss.v086.i08, 2018.

[2] Z. Yaniv, B. C. Lowekamp, H. J. Johnson, R. Beare, “SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research”, J Digit Imaging., https://doi.org/10.1007/s10278-017-0037-8, 31(3): 290-303, 2018.

[3] B. C. Lowekamp, D. T. Chen, L. Ibáñez, D. Blezek, “The Design of SimpleITK”, Front. Neuroinform., 7:45. https://doi.org/10.3389/fninf.2013.00045, 2013.
