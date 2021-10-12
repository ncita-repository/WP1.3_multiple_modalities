# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:51:00 2020

@author: ctorti
"""

"""
Function:
    CopyRoi()
    
Purpose:
    Copy a ROI from one DICOM scan (collection) to another.

Input:
    fromRoiProjLab     - Project Label of the ROI Collection
    
    fromRoiSubjLab     - Subject Label of the ROI Collection
    
    fromRoiExpLab      - Experiment Label of the ROI Collection
    
    fromRoiLabDateTime - DateTime within the ROI Label of the ROI Collection
    
    fromDicomScanLab   - Scan Label of the DICOM scan (collection) with
                        ROI Collection
    
    toDicomProjLab     - Project Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomSubjLab     - Subject Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomExpLab      - Experiment Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
    
    toDicomScanLab     - Scan Label of the DICOM scan (collection) the
                        ROI Collection is to be copied to
                        
    searchBy           - How to define the closest slice with possible values:
                        -- 'x' 
                        -- 'y'
                        -- 'z'
                        -- 'r'
                        to find closest slice along x/y/z direction or r, where
                        r = sqrt(dx^2 + dy^2 + dz^2)
                        
    scanPosCorr        - Scan position correction that needs to be applied to
                        the DICOM scan (collection) to account for different
                        positional reference points (as is the case for 
                        different imaging modalities). The correction will be
                        applied to the same dimension as defined by searchBy.
                        Accepatable values are any float.  Note that if
                        searchBy='r' the correction will not be applied.
                        
    regMethod          - Type of registration to perform on the DICOM scan
                        (collection) that the ROI Collection is to be copied
                        to.  Acceptable values include:
                            -- 'rigid'
                            -- 'affine'
                            -- 'non-rigid'
                            -- 'none' (or 'None' or '')
    
    delDicoms         - Delete downloaded DICOMs if True
    
    debug             - Print results if True

    
Returns:
    roi2               - ROI object that contains the contours in roi1 but 
                        modified for overlaying over the 2nd DICOM scan;
                       - A new DICOM-RTSTRUCT file will be saved to the location
                       specified
                   
    xnatDict           - Dictionary containing XNAT login details, REST path 
                        variables and directory paths
                        (see DownloadFilesForRoiCopy() for details)
                       - The dictionary is exported as a JSON file 
"""



def CopyRoi(fromRoiProjLab, fromRoiSubjLab, fromRoiExpLab, fromRoiLabDateTime,\
            fromDicomScanLab,\
            toDicomProjLab, toDicomSubjLab, toDicomExpLab, toDicomScanLab,\
            searchBy, scanPosCorr, regMethod, delDicoms, debug):
    
    # IMPORT PACKAGES AND FUNCTIONS:
    
    #import importlib
    #import DicomHelperFuncs
    #importlib.reload(DicomHelperFuncs)
    #from DicomHelperFuncs import DownloadFilesForRoiCopy
    from DownloadFilesForRoiCopy import DownloadFilesForRoiCopy
    #from DicomHelperFuncs import GetDicomData
    from GetDicomData import GetDicomData
    #from DicomHelperFuncs import GetRoiData
    #from DicomHelperFuncs import GetClosestSlices
    from GetClosestSlices import GetClosestSlices
    #from DicomHelperFuncs import ModifyRoi
    from ModifyRoi import ModifyRoi
    #from DicomHelperFuncs import RegisterSlices
    from RegisterSlices import RegisterSlices
    import time, os, shutil
    #import SimpleITK as sitk
    #import copy
    
    #import GetClosestSlices
    #import importlib
    #from GetClosestSlices import GetClosestSlices
    
    
    # Keep track of the time it takes to execute certain functions:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    # Get the current date:
    todaysDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    # Assign the root directory:
    rootDir = os.path.join(r'C:\Temp', todaysDate)
    
    """
    If registration is to be applied to the second DICOM scan, the ROIs from
    the first DICOM scan will be copied and modified so that they will overlay
    onto the registered series. To help make things clear, the exported 
    RTSTRUCT file will be stored in a folder with the e.g. directory structure:
        10 zPosCorr affineReg ROIs
        
    where 10 is the Series Number, 'zPosCorr' indicates that the Image 
    Positions of the second scan were corrected, and affineReg indicates that 
    an affine registration method was applied.
    
    If no position correction needs to be made, or no registration applied, the
    structure will be simpler, e.g.:
        10 ROIs
    """
    
    
    # Create a dictionary to store various variables:
    xnatDict = {# Define the XNAT login details:
                'xnatAddress':'http://10.1.1.17',\
                'username':'admin',\
                'password':'admin',\
        
                # Define the REST variables that point to the ROIs to be copied:
                'fromRoiProjLab':fromRoiProjLab,\
                'fromRoiSubjLab':fromRoiSubjLab,\
                'fromRoiExpLab':fromRoiExpLab,\
                #'fromRoiLab':'AIM_20200304_110009',\
                #'fromRoiLab':'AIM_20200306_114104',\
                'fromRoiLabDateTime':fromRoiLabDateTime,\
        
                # Define the REST variables that point to the DICOM files that the ROI applies to:
                'fromDicomScanLab':fromDicomScanLab,\
        
                # Define the REST variables that point to the DICOM files that the ROI is to be applied to:
                'toDicomProjLab':toDicomProjLab,\
                'toDicomSubjLab':toDicomSubjLab,\
                'toDicomExpLab':toDicomExpLab,\
                'toDicomScanLab':toDicomScanLab
                }
    
    # First define the folder structures of the DICOMs and DICOM-RTSTRUCT of
    # the scan whose ROIs are to be copied:
    
    # Where the DICOM-RTSTRUCT to be copied is to be saved to:
    fromRoiDir = os.path.join(rootDir, fromDicomScanLab + ' ROIs')
    
    # Where the DICOM scan that relate to the ROIs to be copies are to be saved 
    # to:
    fromDicomDir = os.path.join(rootDir, fromDicomScanLab + ' DICOMs')
    
    
    # Now deal with the folder structures for the DICOMs and DICOM-RTSTRUCT of
    # the scan that the ROIs are to be copied to:
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    #toDicomDir = os.path.join(rootDir, toDicomScanLab + ' DICOMs')
    
    # For the directory names to export the position-corrected DICOMs, the 
    # position-corrected-and-registered DICOMs, and the DICOM-RTSTRUCT for the
    # latter DICOMs, start by creating a empty character string that will be 
    # appended to:
    dirName = ''
    
    # Add the Scan Label:
    dirName = dirName + toDicomScanLab
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    toDicomDir = os.path.join(rootDir, toDicomScanLab + ' DICOMs')
    
    # If scanPosCorr is non-zero:
    if scanPosCorr:
        dirName = dirName + ' ' + searchBy + 'PosCorr'
      
    # The DICOMs of the series that the ROIs are to be copied to AFTER a change
    # to the Image Position but PRIOR to image registration:
    # If scanPosCorr is non-zero:
    if scanPosCorr:
        toPosCorrDicomDir = os.path.join(rootDir, dirName + ' DICOMs')
    else:
        toPosCorrDicomDir = ''
        
    # The DICOMs that relate to the new ROIs (which may be resulting from 
    # a registration) are to be stored at:
    # If registration is to be applied:
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        dirName = dirName + ' ' + regMethod + 'Reg'
        
        toRegDicomDir = os.path.join(rootDir, dirName + ' DICOMs')
    else:
        toRegDicomDir = ''
    
    # The new (modified copy) ROIs are to be stored at:
    toRoiDir = os.path.join(rootDir, dirName + ' ROIs')
    
    
    #*************************************************************************
    # March 21:  Expand on above but take slightly different approach:
    
    
    """
    Firstly use the DICOM scan label (i.e. Series Number) and create new
    scan labels for DICOMs whose scan positions will be modified, or for
    DICOMs that will be the result of image registration.
    Note that the original scan label is toDicomScanLab. 
    """
    # Convert toScanLab to a number:
    toDicomScanNo = int(toDicomScanLab)
    
    """
    To help makes things clear, the scan numbers for scan position corrected,
    rigid-registered, affine-registered and non-rigid-registered datasets 
    will be the original scan number plus:
     100 if the positions were changed only
     1000, 2000 or 3000 if rigid-, affine- or non-rigid-registered only
     1100, 2100 or 3100 if the positions were changed and then rigid-, affine- 
     or non-rigid-registered
     
    The directory names for the associated different datasets will consist of
    the above scan numbers (converted to strings) plus strings that relate to
    the scan position correction or type of registration if applicable.
    """
    # The scan no and scan label for position-corrected data only:
    toPosCorrScanNo = toDicomScanNo + 100
    toPosCorrScanLab = str(toPosCorrScanNo) + '_' + searchBy + 'PosCorr'
    
    # The scan no for registered-only and for position-corrected and registered
    # data:
    if regMethod=='rigid':
        toRegScanNo = toDicomScanNo + 1000
        
        toPosCorrRegScanNo = toDicomScanNo + 1100
        
    elif regMethod=='affine':
        toRegScanNo = toDicomScanNo + 2000
        
        toPosCorrRegScanNo = toDicomScanNo + 2100
        
    elif regMethod=='non-rigid':
        toRegScanNo = toDicomScanNo + 3000
        
        toPosCorrRegScanNo = toDicomScanNo + 3100
    else:
        #toRegScanNo = copy.deepcopy(toDicomScanNo)
        toRegScanNo = []
        toPosCorrRegScanNo = []
        
    # The scan label for registered-only data: 
    toRegScanLab = str(toRegScanNo) + '_' + regMethod + 'Reg'
    
    # The scan label for position-corrected and registered data:  
    toPosCorrRegScanLab = str(toPosCorrRegScanNo) + '_' + searchBy + 'PosCorr'\
                          + '_' + regMethod + 'Reg'
    
    
    # The directory names for the DICOMs of the series that the ROIs are to be 
    # copied to (prior to any change to Image Position or image registration)
    # and the corresponding new ROIs (the ROIs will use the same scan no. as 
    # the DICOMs:
    toDicomDir = os.path.join(rootDir, toDicomScanLab + ' DICOMs')
    toRoiDir = os.path.join(rootDir, toDicomScanLab + ' ROIs')
    
    # The directory names for position corrected data only:
    toPosCorrDicomDir = os.path.join(rootDir, toPosCorrScanLab + ' DICOMs')
    toPosCorrRoiDir = os.path.join(rootDir, toPosCorrScanLab + ' ROIs')
    
    # The directory names for registered data only:
    toRegDicomDir = os.path.join(rootDir, toRegScanLab + ' DICOMs')
    toRegRoiDir = os.path.join(rootDir, toRegScanLab + ' ROIs')
    
    # The directory names for position corrected and registered data:
    toPosCorrRegDicomDir = os.path.join(rootDir, toPosCorrRegScanLab + ' DICOMs')
    toPosCorrRegRoiDir = os.path.join(rootDir, toPosCorrRegScanLab + ' ROIs')


    
    # Add the above variables + searchBy and scanPosCorr to xnatDict:
    xnatDict.update({
            # Define the directories to download files to:
            'rootDir':rootDir, \
            'fromDicomDir':fromDicomDir,\
            'fromRoiDir':fromRoiDir,\
            'toDicomDir':toDicomDir,\
            'toRoiDir':toRoiDir,\
            'toPosCorrDicomDir':toPosCorrDicomDir,\
            'toPosCorrRoiDir':toPosCorrRoiDir,\
            'toRegDicomDir':toRegDicomDir,\
            'toRegRoiDir':toRegRoiDir,\
            'toPosCorrRegDicomDir':toPosCorrRegDicomDir,\
            'toPosCorrRegRoiDir':toPosCorrRegRoiDir,\
                
            # Store some of the input variables:
            'searchBy':searchBy,
            'scanPosCorr':scanPosCorr
            })
        
  
    
    # DOWNLOAD THE DICOM-RTSTRUCT AND DICOM FILES FROM XNAT:
    xnatDict = DownloadFilesForRoiCopy(xnatDict, debug=debug)
    
    # times[1] = time after downloading of files from XNAT:
    times.append(time.time())
    
    print(f'\nTook {round(times[-1] - times[-2], 1)} s to download data from ',\
                    'XNAT')
    
    # Parse the ROI Collection to be copied:
    """ This is not needed following changes to the helper functions """
    #roiDict1 = GetRoiData(xnatDict['fromRoiFpath'])


    # PARSE THE DICOM DATA THAT CORRESPONDS TO ROI COLLECTION TO BE COPIED:
    dicomDict1 = GetDicomData(xnatDict['fromDicomDir'], sortMethod='slices')
    
    
    # PARSE THE DICOM DATA THAT THE ROI COLLECTION IS TO BE COPIED TO:
    dicomDict2 = GetDicomData(xnatDict['toDicomDir'], sortMethod='slices')

    
    # GET THE MAPPING FROM THE 1st TO 2nd DICOM SCANS 
    # (using the z-coordinate to find the nearest slice): 
    """
    dicomDict1_to_2, dicomDict2posCorr = GetClosestSlices(dicomDict1,\
                                                          dicomDict2,\
                                                          searchBy=searchBy,\
                                                          scanPosCorr=scanPosCorr,\
                                                          debug=debug)
    """
    
    dicomDict1_to_2, dicomDict2posCorr = GetClosestSlices(dicomDict1, \
                                                          dicomDict2, \
                                                          xnatDict, \
                                                          debug=debug)
    
    
    """
    GetClosestSlices() works with dictionaries - so the positional changes
    are applied to the dictionaries and not the DICOMs themselves. So to 
    """
    
    # Print some results for diagnostics purposes:
    """
    There are problems with this:
        KeyError: 'Dicom frame no'

    So for now commenting it out.
    """
    if False: #debug:
        from DicomHelperFuncs import IndexDictInTier2
        
        key1 = IndexDictInTier2(dicomDict1, tier2key='Dicom frame no',\
                                tier2val=5)


        print('Contour copied from Image Position =', \
              dicomDict1[key1]['Image pos'], '\n')
        
        # loop through each key in dicomDict2:
        #keys2 = list(dicomDict2.keys())
        keys2 = list(dicomDict2posCorr.keys())
        
        for i in range(len(keys2)):
            print(f'Slice {i+1} Image Position =', \
                  #dicomDict2[keys2[i]]['Image pos'])
                  dicomDict2posCorr[keys2[i]]['Image pos'])
            
    
    
    # APPLY REGISTRATION IF regMethod INDICATES TO:
    if regMethod in ['rigid', 'affine', 'non-rigid']:
        dicomDict1_to_2, xnatDict = RegisterSlices(dicomDict1_to_2,\
                                                   regMethod=regMethod,\
                                                   xnatDict=xnatDict,\
                                                   debug=debug)
        
        # times[2] = time after performing registration:
        times.append(time.time())
        
        print(f'\nTook {round(times[-1] - times[-2], 1)} s to register images')
        
        
    # COPY AND MODIFY THE ROI COLLECTION 
    # ((from the original DICOM scan) so they may be overlaid on the new DICOM 
    # scan):
    
    # Added 21/03:  It would also be useful to have the ROI copied and modified
    # for the non-registered images too:
    roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2, regMethod='') # 21/03/2020
    
    #roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2)
    #roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2) # 10/03/2020
    roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2, regMethod) # 13/03/2020
    
    # times[3] = time after parsing files, finding closest slice locations 
    # and modifying the ROIs:
    times.append(time.time())
    
    print(f'\nTook {round(times[-1] - times[-2], 1)} s to copy and modify the ',\
                    'ROIs')
    
    # Do some housekeeping - delete the DICOM files in fromDicomDir and 
    # toDicomDir:
    if delDicoms:
        print('\nDeleting temporary DICOM files..')
        for directory in [xnatDict['fromDicomDir'], xnatDict['toDicomDir']]:
            for root, dirs, files in os.walk(directory):
                for f in files:
                    os.unlink(os.path.join(root, f))
                for d in dirs:
                    shutil.rmtree(os.path.join(root, d))
                    
    
    print(f'\nTook {round(times[-1] - times[0], 1)} s to do everything')
    
    #return roi2, xnatDict
    #return dicomDict1, roi1, dicomDict2, roi2, dicomDict1_to_2, xnatDict # 10/03/2020
    return dicomDict1, roi1, dicomDict2posCorr, roi2, dicomDict1_to_2, xnatDict # 15/03/2020
