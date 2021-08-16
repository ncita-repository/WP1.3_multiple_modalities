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
    SourceRoiProjLabel     - Project Label of the source ROI Collection
    
    SourceRoiSubjLabel     - Subject Label of the source ROI Collection
    
    SourceRoiExpLabel      - Experiment Label of the source ROI Collection
    
    SourceRoiLabelDateTime - DateTime within the ROI Label of the source ROI 
                           Collection
    
    SourceSeriesNo         - Series number of the DICOM scan (collection) with
                           source ROI Collection
    
    SourceSeriesDesc       - Series description of the DICOM scan (collection) 
                           with source ROI Collection (Note: This is only used
                           for directory names of downloaded and new created
                           DICOM and ROI Collections)
    
    TargetDicomProjLabel   - Project Label of the target DICOM scan 
                             (collection) the ROI Collection is to be copied to
    
    TargetDicomSubjLabel   - Subject Label of the target DICOM scan 
                            (collection) the ROI Collection is to be copied to
    
    TargetDicomExpLabel    - Experiment Label of the target DICOM scan 
                            (collection) the ROI Collection is to be copied to
    
    TargetSeriesNo         - Series number of the target DICOM scan (collection) 
                            the ROI Collection is to be copied to
                            
    TargetSeriesDesc       - Series description of the target DICOM scan 
                           (collection) the ROI Collection is to be copied to
                           (Note: This is only used for directory names of 
                           downloaded and new created DICOM and ROI Collections)
                        
    SearchBy               - How to define the closest slice with possible 
                            values:
                             -- 'x' 
                             -- 'y'
                             -- 'z'
                             -- 'r'
                            to find closest slice along x/y/z direction or r, 
                            where r = sqrt(dx^2 + dy^2 + dz^2)
                        
    Shift                  - Scan position shift that needs to be applied to
                           the DICOM scan (collection) to account for different
                           positional reference points (as is the case for 
                           different imaging modalities). The correction will 
                           be applied to the same dimension as defined by 
                           SearchBy.
                           Acceptable values are any float.  Note that if
                           SearchBy='r' the correction will not be applied.
                        
    RegMethod              - Type of registration to perform on the DICOM scan
                           (collection) that the ROI Collection is to be copied
                           to.  Acceptable values include:
                            -- 'rigid'
                            -- 'affine'
                            -- 'non-rigid'
                            -- '' (actually any string other than those above)
    
    DelDicoms              - Delete downloaded DICOMs if True
    
    Debug                  - Print results if True

    
Returns:
    SourceDicoms           - List of source DICOM objects
    
    SourceRois             - ROI object belonging to SourceDicoms
    
    TargetDicoms           - List of target DICOM objects
    
    TargetRois             - ROI object belonging to TargetDicoms
    
    ShiftedDicoms          - List of DICOM objects from TargetDicoms whose scan 
                            positions have been shifted by Shift if Shift is 
                            non-zero. If Shift = 0, ShiftedDicoms = [].
                            
    ShiftedRois            - ROI object belonging to ShiftedDicoms if Shift is
                            non-zero. If Shift = 0, ShiftedRois = [].
    
    RemappedDicoms         - List of DICOM objects from SourceDicoms (if S > T),
                            or from TargetDicoms (if S < T and Shift = 0), or
                            from ShiftedDicoms (if S < T and Shift != 0).
                            If S = T, ShiftedDicoms = [].
                            
                            where S = len(SourceDicoms)
                                  T = len(TargetDicoms)
                                  P = len(ShiftedDicoms)
                                  
    RemappedRois           - ROI object belonging to RemappedDicoms
    
    RegisteredDicoms       - List of DICOM objects resulting from image 
                            registration of a list of "fixed" DICOMs and a list
                            of "moving" DICOMs if RegMethod is 'rigid', 'affine' 
                            or 'non-rigid'.  Otherwise, RegisteredRois = []. 
                            The two lists will be of equal length. 
                            Several possible combinations of "fixed" and
                            "moving" DICOMs are possible, e.g.:
                            --- TargetDicoms registered to SourceDicoms if 
                            T = S and Shift = 0;
                            --- TargetDicoms registered to RemappedDicoms if
                            S > T and Shift = 0;
                            --- RemappedDicoms registered to SourceDicoms if
                            S < T and Shift = 0;
                            --- ShiftedDicoms registered to SourceDicoms if 
                            T = S and Shift != 0;
                            --- TargetDicoms registered to RemappedDicoms if 
                            S > T and Shift != 0 (Note that RemappedDicoms 
                            would have followed from ShiftedDicoms);
                            --- RemappedDicoms registered to SourceDicoms if
                            S < T and Shift = 0 (Note that RemappedDicoms 
                            would have followed from ShiftedDicoms);
                            
    RegisteredRois         - ROI object belonging to RegisteredDicoms if 
                            RegMethod is 'rigid', 'affine' or 'non-rigid'.
                            Otherwise, RegisteredRois = [].
                   
    DataDict               - Dictionary containing XNAT login details, REST path 
                            variables, directory paths, file paths, and other
                            details organised in keys (e.g. 'Source', 'Target',
                            'Shifted', 'Remapped', 'Registered')
                           - The dictionary is exported as a JSON file 
"""



#def CopyRoi(SourceRoiProjLabel, SourceRoiSubjLabel, SourceRoiExpLabel, \
#            SourceRoiLabelDateTime, SourceScanLabel,\
#            TargetDicomProjLabel, TargetDicomSubjLabel, TargetDicomExpLabel, TargetScanLabel,\
#            SearchBy, Shift, RegMethod, DelDicoms, Debug):
    
def CopyRoi(SourceRoiProjLabel, SourceRoiSubjLabel, SourceRoiExpLabel, \
            SourceRoiLabelDateTime, SourceSeriesNo, SourceSeriesDesc, \
            TargetDicomProjLabel, TargetDicomSubjLabel, TargetDicomExpLabel, \
            TargetSeriesNo, TargetSeriesDesc, \
            SearchBy, Shift, RegMethod, DelDicoms, Debug, RerunReg):
    
    # IMPORT PACKAGES AND FUNCTIONS:
    import time, os, shutil
    import importlib
    import copy
    
    #import CreateDirs
    #importlib.reload(CreateDirs)
    #from CreateDirs import CreateDirs
    
    #import AddToDataDict
    #importlib.reload(AddToDataDict)
    #from AddToDataDict import AddToDataDict
    
    import CreateDir
    importlib.reload(CreateDir)
    from CreateDir import CreateDir
    
    import CreateDirsForDownloads
    importlib.reload(CreateDirsForDownloads)
    from CreateDirsForDownloads import CreateDirsForDownloads
    
    import DownloadFilesForRoiCopy
    importlib.reload(DownloadFilesForRoiCopy)
    from DownloadFilesForRoiCopy import DownloadFilesForRoiCopy
    
    import GetDicoms
    importlib.reload(GetDicoms)
    from GetDicoms import GetDicoms
    
    import CreateDirsForShifts
    importlib.reload(CreateDirsForShifts)
    from CreateDirsForShifts import CreateDirsForShifts
    
    import ShiftDicoms
    importlib.reload(ShiftDicoms)
    from ShiftDicoms import ShiftDicoms
    
    import GetClosestSlices
    importlib.reload(GetClosestSlices)
    from GetClosestSlices import GetClosestSlices
    
    import CreateDirsForRemap
    importlib.reload(CreateDirsForRemap)
    from CreateDirsForRemap import CreateDirsForRemap
    
    import RemapDicoms
    importlib.reload(RemapDicoms)
    from RemapDicoms import RemapDicoms
    
    import CreateDirsForReg
    importlib.reload(CreateDirsForReg)
    from CreateDirsForReg import CreateDirsForReg
    
    import RegisterSlices
    importlib.reload(RegisterSlices)
    from RegisterSlices import RegisterSlices
    
    import ModifyRoi
    importlib.reload(ModifyRoi)
    from ModifyRoi import ModifyRoi
    
    import PlotRegResults
    importlib.reload(PlotRegResults)
    from PlotRegResults import PlotRegResults
    
    
    
    # Keep track of the time it takes to execute certain functions:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    # Get the current date:
    todaysDate = time.strftime("%Y-%m-%d", time.gmtime())
    
    """
    Define the root directory (RootDir) to be a sub-directory of C:\Temp with
    todaysDate + SourceRoiProjLabel + SourceRoiSubjLabel + SourceRoiExpLabel 
    + ' copied to ' + TargetDicomExpLabel.
    """
    # Assign the root directory:
    #RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel \
    #                       + ' - ' + SourceRoiExpLabel + ' - ' + TargetDicomExpLabel)
    RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel)
    
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
    #"""
    #Only run the following if DataDict doesn't already exist:
    #"""
    #if not 'DataDict' in locals(): 
    ##if not 'DataDict' in globals(): 
    
    DataDict = {# Define the XNAT login details:
                'Login':{'xnatAddress':'http://10.1.1.17',\
                         'username':'admin',\
                         'password':'admin'
                         }, \
                
                # Add the root directory:
                'RootDir':RootDir, \
                
                # Add the other variables:
                'Shift':Shift, \
                'SearchBy':SearchBy, \
                'RegMethod':RegMethod, \
                'Debug':Debug, \
                'RerunReg':RerunReg, \
                
                # Define the REST variables that point to the ROIs to be 
                # copied and the DICOM files that the ROI applies to:
                'Source':{'RoiProjLabel':SourceRoiProjLabel,\
                        'RoiSubjLabel':SourceRoiSubjLabel,\
                        'RoiExpLabel':SourceRoiExpLabel,\
                        'RoiLabelDateTime':SourceRoiLabelDateTime,\
                        #'ScanLabel':SourceScanLabel
                        'SeriesNo':SourceSeriesNo,\
                        'SeriesDesc':SourceSeriesDesc
                        }, \
        
                # Define the REST variables that point to the DICOM files  
                # that the ROI is to be applied to:
                'Target':{'DicomProjLabel':TargetDicomProjLabel,\
                      'DicomSubjLabel':TargetDicomSubjLabel,\
                      'DicomExpLabel':TargetDicomExpLabel,\
                      #'ScanLabel':TargetScanLabel
                      'SeriesNo':TargetSeriesNo,\
                      'SeriesDesc':TargetSeriesDesc
                      }
                }

    
    # Initialise variables that will be updated below:
    SourceRois = []
    TargetRois = []
    ShiftedRois = []
    RemappedRois = []
    FixedDicoms = []
    FixedRois = []
    MovingDicoms = []
    MovingRois = []
    RegisteredDicoms = []
    RegisteredRois = []
    
  
    """ CREATE DIRECTORIES WHERE FILES WILL BE DOWNLOADED FROM XNAT: """   
    #DataDict = AddToDataDict(DataDict)
    DataDict = CreateDirsForDownloads(DataDict)
    
    
    
    """ DOWNLOAD THE DICOM-RTSTRUCT AND DICOM FILES FROM XNAT: """
    DataDict = DownloadFilesForRoiCopy(DataDict)

    #return 0, 1, 2, 3, 4, dataDict
    
    # times[1] = time after downloading of files from XNAT:
    times.append(time.time())
    
    # The time taken to download the data from XNAT:
    Dtime = round(times[-1] - times[-2], 1)
    
    print(f'\nTook {Dtime} s to download data from XNAT.')
    
    
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    # Parse the ROI Collection to be copied:
    """ This is not needed following changes to the helper functions """
    #roiDict1 = GetRoiData(xnatDict['fromRoiFpath'])

    
    
    
    """ GET THE DICOM OBJECTS AND FILEPATHS THAT CORRESPOND TO THE ROI
    COLLECTION TO BE COPIED: """
    #print('\n\nGetting the DICOM objects from', \
    #      DataDict['From']['RoiExpLabel'], '...')
    print('\n\nGetting the DICOM objects from', \
          DataDict['Source']['SeriesDesc'], '...')
    
    SourceDicomFpaths, SourceDicoms = GetDicoms(DataDict['Source']['DicomDir'], \
                                            SortMethod='slices')

    # ADD SourceDicomFpaths TO dataDict:
    DataDict['Source'].update({'DicomFpaths':SourceDicomFpaths})
    
    
    # times[2] = time after reading in all DICOMs that corresond to the ROI
    # Collection:
    times.append(time.time())
    
    # The time taken to load all DICOMs in the "Source" scan:
    Dtime = round(times[-1] - times[-2], 1)
    
    #print(f'\nTook {Dtime} s to load all DICOMs in scan ' + SourceScanLabel + '.')
    print(f'\n   Took {Dtime} s to load all DICOMs in scan ' \
          + SourceSeriesNo + '.')
    

    
    
    """ GET THE DICOM OBJECTS AND FILEPATHS THAT CORRESPOND TO THE SCAN THAT
    THE ROI COLLECTION IS TO BE COPIED TO: """
    #"""
    #Only run GetDicoms if Rerun=True and if 'ToDicoms' doesn't already exist:
    #"""
    #if Rerun and not 'ToDicoms' in locals():
    ##if Rerun and not 'ToDicoms' in globals():
    
    #print('\n\nGetting the DICOM objects from', \
    #      DataDict['To']['DicomExpLabel'], '...')
    print('\n\nGetting the DICOM objects from', \
          DataDict['Target']['SeriesDesc'], '...')
    
    TargetDicomFpaths, TargetDicoms = GetDicoms(DataDict['Target']['DicomDir'], SortMethod='slices')
        
    # ADD TargetDicomFpaths TO dataDict:
    DataDict['Target'].update({'DicomFpaths':TargetDicomFpaths})
                
    
    # times[3] = time after reading in all DICOMs that the ROI Collection is to
    # be copied to:
    times.append(time.time())
    
    # The time taken to load all DICOMs in the "to" scan:
    Dtime = round(times[-1] - times[-2], 1)
    
    #print(f'\nTook {Dtime} s to load all DICOMs in scan ' + TargetScanLabel + '.')
    print(f'\n   Took {Dtime} s to load all DICOMs in scan ' \
          + TargetSeriesNo + '.')
    
    
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
   
    """ Commenting this section b/c the mapping between the two DICOM series
    needs to be done first, by running GetClosestSlices()! """
    if False:
        """ COPY THE ROI COLLECTION FOR SourceDicos TO TargetDicoms: """
        """ (This is only useful as a test of the ROI copying algorithm) """
        
        # Initialise variables that will be updated below:
        SourceRois = []
        TargetRois = []
        ShiftedRois = []
        RemappedRois = []
        
        # The key of the series whose ROIs are to be copied:
        CopyFromKey = 'Source'
        
        # The key of the series to copy the ROIs to:
        CopyToKey = 'Target'
        
        print('\n\nCopying the ROI Collection from', \
              DataDict[CopyFromKey]['SeriesDesc'], 'to', \
              DataDict[CopyToKey]['SeriesDesc'], '...')
        
        """ COPY THE ROI COLLECTION FROM SourceDicoms TO RemappedDicoms: """
        SourceRois, TargetRois, DataDict = ModifyRoi(SourceDicoms, TargetDicoms, DataDict, CopyFromKey, CopyToKey)       
        
        
        return SourceDicoms, SourceRois, TargetDicoms, TargetRois, 0, 0, \
               0, 0, 0, 0, DataDict        
            
            
    
    
    

    
    
    
    
    
    
    
    """ PROCEED WITH SHIFTS ONLY IF Shift IS NON-ZERO: """
    if Shift:
        if Debug:
            print('\nKeys in DataDict prior to running CreateDirsForShifts:\n', \
                  DataDict.keys())
        
        """ CREATE DIRECTORIES FOR SHIFTED DICOMS AND ROIS: """  
        
        # The key of the series to be shifted:
        ToShiftKey = 'Target'
        
        DataDict = CreateDirsForShifts(DataDict, ToShiftKey)
        
        if Debug:
            print('\nKeys in DataDict after running CreateDirsForShifts:\n', \
                  DataDict.keys())
        
            print('\nSeries Description of ShiftedTarget after CreateDirsForShifts() =', \
                  DataDict['ShiftedTarget']['SeriesDesc'])
        
        """ APPLY SHIFT TO POSITIONS in TargetDicoms: """
        
        ShiftedDicoms, DataDict = ShiftDicoms(TargetDicoms, DataDict)
        
        if Debug:
            print('\nSeries Description of ShiftedTarget after ShiftDicoms() =', \
                  DataDict['ShiftedTarget']['SeriesDesc'])
    else:
        ShiftedDicoms = []

    
    
    #return FromDicoms, 0, ToDicoms, 0, ShiftedDicoms, 0, \
    #       0, 0, 0, 0, DataDict
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')



    """ GET THE MAPPING FROM SourceDicoms TO TargetDicoms (if Shift=0) OR TO
    ShiftedDicoms (IF Shift!=0)
    (using the z-coordinate to find the nearest slice): """
    #print('\n\nGetting the indeces that map each slice in', \
    #      DataDict['From']['RoiExpLabel'], 'to each slice in ', \
    #      DataDict['To']['DicomExpLabel'], '...')
    ##print('\n\nGetting the indeces that map each slice in', \
    ##      DataDict['From']['SeriesDesc'], 'to each slice in ', \
    ##      DataDict['To']['SeriesDesc'], '...')
    #PosCorrDicoms, \
    #MappedDicoms, \
    #DataDict = GetClosestSlices(FromDicoms, ToDicoms, DataDict)
    
    print('\n\nGetting the indeces that map each slice in', \
              DataDict['Target']['DicomExpLabel'], 'to each slice in ', \
              DataDict['Source']['RoiExpLabel'], '...')
    
    if Shift:
        DataDict = GetClosestSlices(SourceDicoms, ShiftedDicoms, DataDict)
        
    else:
        DataDict = GetClosestSlices(SourceDicoms, TargetDicoms, DataDict)
    
    #return SourceDicoms, 0, TargetDicoms, 0, PosCorrDicoms, MappedDicoms, 0, DataDict
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
    
    """ It's more sensible to copy from SourceDicoms to MappedDicoms rather than
    having to use MappingInds within ModifyRoi()... """
    if False:
        """ COPY THE ROI COLLECTION OF SourceDicoms TO ShiftedDicoms: """
        
        # Initialise variables that will be updated below:
        SourceRois = []
        TargetRois = []
        ShiftedRois = []
        RemappedRois = []
        
        # The key of the series whose ROIs are to be copied:
        CopyFromKey = 'Source'
        
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['ShiftedKey']
        
        print('\n\nCopying the ROI Collection from', \
                  DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                  DataDict[CopyToKey]['SeriesDesc'], '...')
            
        """ COPY THE ROI COLLECTION FROM SourceDicoms TO ShiftedDicoms: """
        SourceRois, ShiftedRois, DataDict = ModifyRoi(SourceDicoms, ShiftedDicoms, DataDict, CopyFromKey, CopyToKey)  
        
        
        return SourceDicoms, SourceRois, TargetDicoms, 0, ShiftedDicoms, ShiftedRois, \
               0, 0, 0, 0, DataDict        
               
               
               
    
        
        
    
    """ USE THE MAPPING INDECES TO REMAP THE DICOMS OF GREATER LENGTH: """
    """ (This will be useful for performing registrations) """
    if len(SourceDicoms) != len(TargetDicoms):
        if len(SourceDicoms) > len(TargetDicoms):
            # The key of the series to be remapped:
            ToRemapKey = 'Source'
            
            print('\nThere are more slices SourceDicoms than in TargetDicoms', \
                  'so', ToRemapKey, 'will be remapped.')
            
            """ CREATE DIRECTORIES FOR REMAPPED DICOMS AND ROIS: """   
            DataDict = CreateDirsForRemap(DataDict, ToRemapKey)
            
            """ RE-MAP SourceDicoms: """
            RemappedDicoms, DataDict = RemapDicoms(SourceDicoms, DataDict)
        
        elif len(SourceDicoms) < len(TargetDicoms):
            if Shift:
                # The key of the series to be remapped:
                ToRemapKey = DataDict['ShiftedKey']
                
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so', ToRemapKey, 'will be remapped.')
                
                """ CREATE DIRECTORIES FOR REMAPPED DICOMS AND ROIS: """   
                DataDict = CreateDirsForRemap(DataDict, ToRemapKey)
            
                """ RE-MAP ShiftedTargetDicoms: """
                RemappedDicoms, DataDict = RemapDicoms(ShiftedDicoms, DataDict)
            else:
                # The key of the series to be remapped:
                ToRemapKey = 'Target'
                
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so', ToRemapKey, 'will be remapped.')
                
                """ CREATE DIRECTORIES FOR REMAPPED DICOMS AND ROIS: """   
                DataDict = CreateDirsForRemap(DataDict, ToRemapKey)
                
                """ RE-MAP TargetDicoms: """
                RemappedDicoms, DataDict = RemapDicoms(TargetDicoms, DataDict)
    else:
        print('\nThere are equal numbers of slices in SourceDicoms and', \
              'TargetDicoms so no DICOMs will be remapped.')
                
    
    #print(f'\nlen(FromDicoms) = {len(FromDicoms)}')        
        
    
    #return FromDicoms, 0, ToDicoms, 0, ShiftedDicoms, 0, \
    #       RemappedDicoms, 0, 0, 0, DataDict
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
    
    """ April 6: Section below replaced for different cases below... """
    if False:
        """ COPY THE ROI COLLECTION OF SourceDicoms TO ShiftedDicoms: """
        
        """ 
        THIS HAS TO BE DONE IN TWO STEPS:
            1. Copy ROIs from SourceDicoms to RemappedDicoms (= SourceDicoms remapped
            so that they have the same number of DICOMs as ShiftedDicoms)
            
            2. Copy the ROIs from RemappedDicoms to ShiftedDicoms
        """
        
        # Initialise variables that will be updated below:
        SourceRois = []
        TargetRois = []
        ShiftedRois = []
        RemappedRois = []
        
        """ STEP 1:  Copy ROIs from SourceDicoms to RemappedDicoms: """
    
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['RemappedKey']
        
        print('\n\nCopying the ROI Collection from', \
                  DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                  DataDict[CopyToKey]['SeriesDesc'], '...')
        
        """ COPY THE ROI COLLECTION FROM SourceDicoms TO RemappedDicoms: """
        SourceRois, RemappedRois, DataDict = ModifyRoi(SourceDicoms, RemappedDicoms, DataDict, CopyFromKey, CopyToKey)  
    
    
        """ STEP 2:  Copy ROIs from RemappedDicoms to ShiftedDicoms: """
    
        # The key of the series whose ROIs are to be copied:
        CopyFromKey = DataDict['RemappedKey']
        
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['ShiftedKey']
        
        print('\n\nCopying the ROI Collection from', \
                  DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                  DataDict[CopyToKey]['SeriesDesc'], '...')
            
        """ COPY THE ROI COLLECTION FROM RemappedDicoms TO ShiftedDicoms: """
        RemappedRois, ShiftedRois, DataDict = ModifyRoi(RemappedDicoms, ShiftedDicoms, DataDict, CopyFromKey, CopyToKey)  
    
    
    
    
    
    
    
    """ Modified from above on April 6: """
    
    """
    Note:
        If the length of SourceDicoms is less than or equal to the length of
        TargetDicoms, this can be done in a single step (as in older version
        of code above).
        
        If the length of SourceDicoms is greater than that of TargetDicoms, it
        might be necessary to do this in two steps: 
            1. Copy ROIs from SourceDicoms to RemappedDicoms (= Remapped 
            SourceDicoms)
            2. Copy the ROIs from RemappedDicoms to ShiftedDicoms (= Shifted
            TargetDicoms)
    """
    
    
    """ COPY THE ROI COLLECTION TO THE PRE-REGISTERED SERIES: """
    
    """ APRIL 7:  I think the following is incorrect, since regardless of the
    relative lengths of SourceDicoms and TargetDicoms, the ROIs will always be
    copied from SourceRois! So changes made below this section: """
      
    if False:
        # The key of the series whose ROIs are to be copied:
        if len(SourceDicoms) != len(TargetDicoms):
            if len(SourceDicoms) > len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = DataDict['RemappedKey']
                CopyFromDicoms = copy.deepcopy(RemappedDicoms)
                
                # The key of the series to copy the ROIs to:
                if Shift:
                    CopyToKey = DataDict['ShiftedKey']
                    CopyToDicoms = copy.deepcopy(ShiftedDicoms)
                else:
                    CopyToKey = 'Target'
                    CopyToDicoms = copy.deepcopy(TargetDicoms)
                    
                print('\nThere are more slices SourceDicoms than in TargetDicoms', \
                      'so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                if Shift:
                    RemappedRois, ShiftedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
                else:
                    RemappedRois, TargetRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
            
            elif len(SourceDicoms) < len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = 'Source'
                CopyFromDicoms = copy.deepcopy(SourceDicoms)
                    
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['RemappedKey']
                CopyToDicoms = copy.deepcopy(RemappedDicoms)
                    
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                SourceRois, RemappedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
                    
                    
        else:
            # The key of the series whose ROIs are to be copied:
            CopyFromKey = 'Source'
            CopyFromDicoms = copy.deepcopy(SourceDicoms)
            
            # The key of the series to copy the ROIs to:
            if Shift:
                CopyToKey = DataDict['ShiftedKey']
                CopyToDicoms = copy.deepcopy(ShiftedDicoms)
            else:
                CopyToKey = 'Target'
                CopyToDicoms = copy.deepcopy(TargetDicoms)
                
            print('\nThere are equal numbers of slices in SourceDicoms and', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
            
            """ COPY THE ROI COLLECTION: """
            if Shift:
                SourceRois, ShiftedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
            else:
                SourceRois, TargetRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
            
           
            
    # The key of the series whose ROIs are to be copied:
    CopyFromKey = 'Source'
    CopyFromDicoms = copy.deepcopy(SourceDicoms)
    
    # Get list of keys in DataDict:
    keys = list(DataDict.keys())
    
    # If 'RemappedShiftedTarget' exists, copy for those DICOMs:
    if 'RemappedShiftedTarget' in keys:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['RemappedKey']
        CopyToDicoms = copy.deepcopy(RemappedDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, RemappedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)
        
    # If 'ShiftedTarget' exists, copy for those DICOMs:
    if 'ShiftedTarget' in keys:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['ShiftedKey']
        CopyToDicoms = copy.deepcopy(ShiftedDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, ShiftedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
        
    # Copy ROIs for Target DICOMs:
    # The key of the series to copy the ROIs to:
    CopyToKey = 'Target'
    CopyToDicoms = copy.deepcopy(TargetDicoms)
        
    print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
    
    """ COPY THE ROI COLLECTION: """
    SourceRois, TargetRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
    
    
    """ April 7:  This is fine up to here: """
    return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
           ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
           FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
           RegisteredDicoms, RegisteredRois, DataDict   
           
   
        
        

           
        
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
    
    """ Print some results for diagnostics purposes: """
    """
    There are problems with this:
        KeyError: 'Dicom frame no'

    Plus I'm now working on arrays of DICOMs rather than dictionaries.
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
            
    
    
    """ April 7: Section below replaced for different cases below... """
    if False:
        """ PROCEED WITH REGISTRATION ONLY IF RegMethod IS ONE OF THE EXPECTED
        REGISTRATION METHODS: """
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
            # Define the keys within DataDict that determine which filepaths will
            # be used in RegisterSlices(), and where the registered images should
            # be exported to:
            FixedKey = 'Source'
            #MovingKey = 'MapPosCorr'
            #MovingKey = 'Mapped'
            # If Shift is non-zero, ShiftedDicoms will be registered to FromDicoms.
            # If Shift = 0, ToDicoms will be registered to FromDicoms.
            if Shift:
                MovingKey = 'Shifted'
            else:
                MovingKey = 'Target'
            #RegKey = 'RegMapPosCorr'
            RegKey = 'Registered'
            
            
            print('\nKeys in DataDict prior to running CreateDirsForReg:\n', \
                      DataDict.keys())
            
            
            """ CREATE DIRECTORIES FOR SHIFTED DICOMS AND ROIS: """  
            DataDict = CreateDirsForReg(DataDict, FixedKey, MovingKey)
            
            
            print('\nKeys in DataDict after running CreateDirsForReg:\n', \
                      DataDict.keys())
            
            
            """ APPLY REGISTRATION: """
            
            #print('\n\nRegistering images in', DataDict[MovingKey]['DicomExpLabel'], \
            #       'to', DataDict[FixedKey]['DicomExpLabel'], '...')
            print('\n\nRegistering images in', DataDict[MovingKey]['SeriesDesc'], \
                   'to', DataDict[FixedKey]['SeriesDesc'], '...')
            
            RegDicoms, DataDict = RegisterSlices(DataDict, FixedKey, MovingKey, RegKey)
            
            # times[2] = time after performing registration:
            times.append(time.time())
            
            # The time taken to register images:
            Dtime = round(times[-1] - times[-2], 1)
            
            print(f'\n   Took {Dtime} s to register images.')
            
            #return FromDicoms, 0, ToDicoms, 0, PosCorrDicoms, MappedDicoms, RegDicoms, DataDict
        
        
        
        
    """ Modified from above on April 7: """
    
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        """ PERFORM REGISTRATION: """
        
        # The key name of the registered series:
        RegKey = 'Registered'
    
        # Define FixedKey and MovingKey:
        if len(SourceDicoms) != len(TargetDicoms):
            if len(SourceDicoms) > len(TargetDicoms):
                # The key of the series to be fixed:
                FixedKey = DataDict['RemappedKey']
                
                # The key of the series to be registered to the fixed series:
                if Shift:
                    MovingKey = DataDict['ShiftedKey']
                else:
                    MovingKey = 'Target'
                    
                print('\nThere are more slices SourceDicoms than in', \
                      'TargetDicoms so FixedKey =', FixedKey, \
                      ', and MovingKey =', MovingKey, '.')
            
            elif len(SourceDicoms) < len(TargetDicoms):
                # The key of the series to be fixed:
                FixedKey = 'Source'
                    
                # The key of the series to be registered to the fixed series:
                MovingKey = DataDict['RemappedKey']
                    
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so FixedKey =', FixedKey, \
                  ', and MovingKey =', MovingKey, '.')
                
        else:
            # The key of the series to be fixed:
            FixedKey = 'Source'
            
            # The key of the series to be registered to the fixed series:
            if Shift:
                MovingKey = DataDict['ShiftedKey']
            else:
                MovingKey = 'Target'
                
            print('\nThere are equal numbers of slices in SourceDicoms and', \
                  'TargetDicoms so FixedKey =', FixedKey, \
                  ', and MovingKey =', MovingKey, '.')
            
            
        """ CREATE DIRECTORIES FOR REGISTERED DICOMS AND ROIS: """  
        DataDict = CreateDirsForReg(DataDict, FixedKey, MovingKey, RegKey)
        
        
        print('\nKeys in DataDict after running CreateDirsForReg:\n', \
                  DataDict.keys())
        
        
        """ APPLY REGISTRATION: """
        
        #print('\n\nRegistering images in', DataDict[MovingKey]['DicomExpLabel'], \
        #       'to', DataDict[FixedKey]['DicomExpLabel'], '...')
        print('\n\nRegistering images in', DataDict[MovingKey]['SeriesDesc'], \
               'to', DataDict[FixedKey]['SeriesDesc'], '...')
        
        FixedDicoms, MovingDicoms, RegisteredDicoms, DataDict = RegisterSlices(DataDict, FixedKey, MovingKey, RegKey)
        
        # times[2] = time after performing registration:
        times.append(time.time())
        
        # The time taken to register images:
        Dtime = round(times[-1] - times[-2], 1)
        
        print(f'\n   Took {Dtime} s to register images.')
        
        #return FromDicoms, 0, ToDicoms, 0, PosCorrDicoms, MappedDicoms, RegDicoms, DataDict 
        #return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
        #       ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
        #       FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
        #       RegisteredDicoms, RegisteredRois, DataDict  
        

        
            
        print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
        
        
        
        """ COPY THE ROI COLLECTION FOR FixedDicoms: """
        
        # The key of the series whose ROIs are to be copied:
        if len(SourceDicoms) != len(TargetDicoms):
            if len(SourceDicoms) > len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = DataDict['RemappedKey']
                CopyFromDicoms = copy.deepcopy(RemappedDicoms)
                
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['FixedKey']
                CopyToDicoms = copy.deepcopy(FixedDicoms)
                    
                print('\nThere are more slices SourceDicoms than in TargetDicoms', \
                      'so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                RemappedRois, FixedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)   
            
            elif len(SourceDicoms) < len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = 'Source'
                CopyFromDicoms = copy.deepcopy(SourceDicoms)
                    
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['FixedKey']
                CopyToDicoms = copy.deepcopy(FixedDicoms)
                    
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                SourceRois, FixedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
                    
                    
        else:
            # The key of the series whose ROIs are to be copied:
            CopyFromKey = 'Source'
            CopyFromDicoms = copy.deepcopy(SourceDicoms)
            
            # The key of the series to copy the ROIs to:
            CopyToKey = DataDict['FixedKey']
            CopyToDicoms = copy.deepcopy(FixedDicoms)
                
            print('\nThere are equal numbers of slices in SourceDicoms and', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
            
            """ COPY THE ROI COLLECTION: """
            SourceRois, FixedRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
        
        
        """ April 7: Seems fine up to here:  """
        #return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
        #       ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
        #       FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
        #       RegisteredDicoms, RegisteredRois, DataDict 
    
    
    
        print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
        
        
        
        
        """ COPY THE ROI COLLECTION FOR MovingDicoms: """
        
        # The key of the series whose ROIs are to be copied:
        if len(SourceDicoms) != len(TargetDicoms):
            if len(SourceDicoms) > len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = DataDict['RemappedKey']
                CopyFromDicoms = copy.deepcopy(RemappedDicoms)
                
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['MovingKey']
                CopyToDicoms = copy.deepcopy(MovingDicoms)
                    
                print('\nThere are more slices SourceDicoms than in TargetDicoms', \
                      'so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                RemappedRois, MovingRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)   
            
            elif len(SourceDicoms) < len(TargetDicoms):
                # The key of the series to be copied:
                CopyFromKey = 'Source'
                CopyFromDicoms = copy.deepcopy(SourceDicoms)
                    
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['MovingKey']
                CopyToDicoms = copy.deepcopy(FixedDicoms)
                    
                print('\nThere are fewer slices in SourceDicoms than in', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
                
                """ COPY THE ROI COLLECTION: """
                SourceRois, MovingRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
                    
                    
        else:
            # The key of the series whose ROIs are to be copied:
            CopyFromKey = 'Source'
            CopyFromDicoms = copy.deepcopy(SourceDicoms)
            
            # The key of the series to copy the ROIs to:
            CopyToKey = DataDict['MovingKey']
            CopyToDicoms = copy.deepcopy(FixedDicoms)
                
            print('\nThere are equal numbers of slices in SourceDicoms and', \
                  'TargetDicoms so', CopyFromKey, 'will be copied to', CopyToKey, '.')
            
            """ COPY THE ROI COLLECTION: """
            SourceRois, MovingRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)  
        
        
        """ April 7: Seems fine up to here:  """
        return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
               ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
               FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
               RegisteredDicoms, RegisteredRois, DataDict 
               
               
        
        
        print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
        
        
        
    
        
        """ Plot registration results: """
        if True:#Debug:
            # Import function:
            from PlotRegResults import PlotRegResults
            
            # Get the slice numbers:
            SourceSliceNos = DataDict[FixedKey]['SliceNos']
            #MappedSliceNos = DataDict[MovingKey]['SliceNos']
            ShiftedSliceNos = DataDict[MovingKey]['SliceNos']
            RegSliceNos = DataDict[RegKey]['SliceNos']
            
            # Get the current date-time:
            CurrentDateTime = time.strftime("%Y%m%d_%H%M%H", time.gmtime())
            
            # Create the filename for the exported plot:
            PlotFname = 'Registration results (' + CurrentDateTime + ').jpg'
            
            # Define the directory to export the plot:
            PlotDir = DataDict[RegKey]['RoiDir']
            
            # Create the filepath for the exported plot:
            PlotFpath = os.path.join(PlotDir, PlotFname)
            
            #ExportFpath = '' # if empty, not plot will be exported
            
            """
            March 31:  Above code has a problem - not sure what, so for now:
            """
            #ExportFpath = r'C:\Temp\2020-03-31 ACRIN - ACRIN-FMISO-Brain-011_MR_2 - ACRIN-FMISO-Brain-011_CT_PET_1\30112 ROIs (3D_Brain miso- zPosCorr - Mapped - affine-Reg) (3D_Brain miso- zPosCorr - Mapped - affine-Reg)\Registration results.jpg'
            
            # Get the mappping indeces:
            #Inds = DataDict[MovingKey]['Inds']
            
            # Use the indeces from MovingKey to index the mapped slices in
            # MappedDicoms:
            #print(f'\nThere are {len(MappedDicoms)} DICOMs in MappedDicoms, ', \
            #      f'and {len(Inds)} indeces in Inds.')
            #IndMappedDicoms = [MappedDicoms[i] for i in Inds]
            
            #PlotRegResults(FromDicoms, IndMappedDicoms, RegDicoms, \
            #               FromSliceNos, MappedSliceNos, RegSliceNos)
            
            #PlotRegResults(FromDicoms, MappedDicoms, RegDicoms, \
            #               FromSliceNos, MappedSliceNos, RegSliceNos, \
            #               PlotFpath)
            
            PlotRegResults(SourceDicoms, ShiftedDicoms, RegDicoms, \
                           SourceSliceNos, ShiftedSliceNos, RegSliceNos, \
                           PlotFpath)
    else:
        RegDicoms = []

            
    
    """
    March 31:  I don't think it makes sense to run this section below because
    the "Target" series is not mapped.  So commenting this section out..
    """
    #""" COPY THE ROI COLLECTION FROM FromDicoms TO ToDicoms: """
    
    
    # Define the keys within DataDict that determine which filepaths will
    # be used in ModifyRoi(), and where the new DICOM-RTSTRUCT file should be
    # exported to::
    #FromKey = 'From'
    #ToKey = 'To'
    
    ##print('\n\nCopying the ROI Collection from', \
    ##      DataDict[FromKey]['DicomExpLabel'], 'to', \
    ##      DataDict[ToKey]['DicomExpLabel'], '...')
    #print('\n\nCopying the ROI Collection from', \
    #      DataDict[FromKey]['SeriesDesc'], 'to', \
    #      DataDict[ToKey]['SeriesDesc'], '...')
 
    # Use the above keys to get the ROI filepath, ROI directory and mapped
    # indeces:
    ###CopyFromRoiFpath = DataDict[CopyFromKey]['RoiFpath']
    ###CopyToRoiDir = DataDict[CopyToKey]['RoiDir']
    ##Inds = DataDict[ToKey]['Inds'] 
    
    # Use the indeces from ToKey to index the mapped slices in MappedDicoms:
    ##CopyToDicoms = [MappedDicoms[i] for i in Inds]
    
    # Copy and modify the DICOM-RTSTRUCT file:
    ##roi1, roi2, xnatDict = ModifyRoi(CopyFromRoiFpath, CopyToRoiDir, FromDicoms, IndMappedDicoms, Inds, DataDict)
    #FromRois, ToRois, DataDict = ModifyRoi(FromDicoms, ToDicoms, DataDict, FromKey, ToKey)
    
    #return FromDicoms, FromRoi, ToDicoms, ToRoi, PosCorrDicoms, MappedDicoms, RegDicoms, DataDict
    
    
    """ April 2: Commmenting this out for the version below that handles more cases """
    if False:
        """ COPY THE ROI COLLECTION FROM FromDicoms TO THE PRE-REGISTERED SERIES: """
        
        # Define the keys within DataDict that determine which filepaths will
        # be used in ModifyRoi(), and where the new DICOM-RTSTRUCT file should be
        # exported to::
        FromKey = 'From'
        #ToKey = 'MapPosCorr'
        #ToKey = 'Mapped'
        if Shift:
            ToKey = 'Shifted'
        else:
            ToKey = 'To'
        
        #print('\n\nCopying the ROI Collection from', \
        #      DataDict[FromKey]['DicomExpLabel'], 'to', \
        #      DataDict[ToKey]['DicomExpLabel'], '...')
        print('\n\nCopying the ROI Collection from', \
              DataDict[FromKey]['SeriesDesc'], 'to', \
              DataDict[ToKey]['SeriesDesc'], '...')
     
        # Copy and modify the DICOM-RTSTRUCT file:
        #FromRois, MapPosCorrRois, DataDict = ModifyRoi(FromDicoms, MappedDicoms, DataDict, FromKey, ToKey)
        if Shift:
            FromRois, ShiftedRois, DataDict = ModifyRoi(FromDicoms, ShiftedDicoms, DataDict, FromKey, ToKey)
            
            ToRois = []
        else:
            FromRois, ToRois, DataDict = ModifyRoi(FromDicoms, ToDicoms, DataDict, FromKey, ToKey)
            
            ShiftedRois = []
        
        #return FromDicoms, FromRoi, ToDicoms, MapPosCorrRoi, PosCorrDicoms, MappedDicoms, RegDicoms, DataDict
    
    
    
    
    
    if False:
        """ COPY THE ROI COLLECTION TO THE PRE-REGISTERED SERIES: """
        """ commented for simpler version below """
        # Initialise variables that will be updated below:
        FromRois = []
        ToRois = []
        ShiftedRois = []
        RemappedRois = []
        
        
        if len(FromDicoms) != len(ToDicoms):
            if len(FromDicoms) > len(ToDicoms):
                # The key of the series whose ROIs are to be copied:
                CopyFromKey = DataDict['RemappedKey']
                
                if Shift:
                    # The key of the series to copy the ROIs to:
                    CopyToKey = DataDict['ShiftedKey']
                    
                    print('\n\nCopying the ROI Collection from', \
                          DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                          DataDict[CopyToKey]['SeriesDesc'], '...')
                    
                    """ COPY THE ROI COLLECTION FROM RemappedDicoms TO ShiftedDicoms: """
                    RemappedRois, ShiftedRois, DataDict = ModifyRoi(RemappedDicoms, ShiftedDicoms, DataDict, CopyFromKey, CopyToKey)
                    
                else:
                    # The key of the series to copy the ROIs to:
                    CopyToKey = 'To'
                    
                    print('\n\nCopying the ROI Collection from', \
                          DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                          DataDict[CopyToKey]['SeriesDesc'], '...')
                    
                    """ COPY THE ROI COLLECTION FROM RemappedDicoms TO ToDicoms: """
                    RemappedRois, ToRois, DataDict = ModifyRoi(RemappedDicoms, ToDicoms, DataDict, CopyFromKey, CopyToKey)
            
            elif len(FromDicoms) < len(ToDicoms):
                # The key of the series whose ROIs are to be copied:
                CopyFromKey = 'From'
    
                if Shift:
                    # The key of the series to copy the ROIs to:
                    CopyToKey = DataDict['ShiftedKey']
                    
                    print('\n\nCopying the ROI Collection from', \
                          DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                          DataDict[CopyToKey]['SeriesDesc'], '...')
                    
                    """ COPY THE ROI COLLECTION FROM FromDicoms TO ShiftedDicoms: """
                    FromRois, ShiftedRois, DataDict = ModifyRoi(FromDicoms, ShiftedDicoms, DataDict, CopyFromKey, CopyToKey)
                else:
                    # The key of the series to copy the ROIs to:
                    CopyToKey = 'To'
                    
                    print('\n\nCopying the ROI Collection from', \
                          DataDict[CopyFromKey]['SeriesDesc'], 'to', \
                          DataDict[CopyToKey]['SeriesDesc'], '...')
                    
                    """ COPY THE ROI COLLECTION FROM FromDicoms TO ShiftedDicoms: """
                    RemappedRois, ToRois, DataDict = ModifyRoi(FromDicoms, ToDicoms, DataDict, CopyFromKey, CopyToKey)
                    
        else:
            # The key of the series whose ROIs are to be copied:
            CopyFromKey = 'From'
        
            if Shift:
                # The key of the series to copy the ROIs to:
                CopyToKey = 'Shifted'
                
                """ COPY THE ROI COLLECTION FROM FromDicoms TO ShiftedDicoms: """
                FromRois, ShiftedRois, DataDict = ModifyRoi(FromDicoms, ShiftedDicoms, DataDict, CopyFromKey, CopyToKey)
            else:
                # The key of the series to copy the ROIs to:
                CopyToKey = 'To'
                
                """ COPY THE ROI COLLECTION FROM FromDicoms TO ToDicoms: """
                FromRois, ToRois, DataDict = ModifyRoi(FromDicoms, ToDicoms, DataDict, CopyFromKey, CopyToKey)
    
    
    
    
    """ COPY THE ROI COLLECTION TO THE PRE-REGISTERED SERIES: """
    """ simpler version than that above... """
    
    # Initialise variables that will be updated below:
    SourceRois = []
    TargetRois = []
    ShiftedRois = []
    RemappedRois = []
    
    
    """ 
    The ROIs will be copied from SourceDicoms unless len(SourceDicoms) > len(TargetDicoms),
    in which case, the ROIs will be copied from RemappedDicoms (i.e. remapped 
    SourceDicoms).
    """
    
    # The key of the series whose ROIs are to be copied:
    CopyFromKey = 'Source'
    
    # And the DICOMs that belong to that series:
    CopyFromDicoms = copy.deepcopy(SourceDicoms)
    
    """
    The ROIs will be copied to TargetDicoms or ShiftedDicoms:
    """
    if Shift:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['ShiftedKey']
        
        # And the DICOMs that belong to that series:
        CopyToDicoms = copy.deepcopy(ShiftedDicoms)
        
    else:
        # The key of the series to copy the ROIs to:
        CopyToKey = 'Target'
        
        # And the DICOMs that belong to that series:
        CopyToDicoms = copy.deepcopy(TargetDicoms)
                
                
    # Modify CopyFromKey and CopyFromDicoms only if 
    # len(FromDicoms) > len(ToDicoms):
    if len(SourceDicoms) != len(TargetDicoms):
        if len(SourceDicoms) > len(TargetDicoms):
            """ Need to make a two-stage copy: First need to copy ROIs from 
            SourceDicoms to RemappedSourceDicoms; Then need to copy ROIs from
            RemappedSourceDicoms to ShiftedDicoms (if Shift !=0) or to TargetDicoms
            (if Shift=0). """
            # First copy ROIs from SourceDicoms to RemappedSourceDicoms:
            # Note:  CopyFromKey = 'Source'
            # and CopyFromDicoms = SourceDicoms
            
            # The key of the series to copy the ROIs to:
            CopyToKey = DataDict['RemappedKey']
            
            # And the DICOMs that belong to that series:
            CopyToDicoms = copy.deepcopy(RemappedDicoms)
            
            """ COPY THE ROI COLLECTION FROM SourceDicoms TO RemappedDicoms: """
            CopyFromRois, CopyToRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)       
            
            # Copy the ROIs before ModifyRoi is run again:
            SourceRois = copy.deepcopy(CopyFromRois)
            RemappedRois = copy.deepcopy(CopyToRois)
            
            # The next key of the series whose ROIs are to be copied:
            CopyFromKey = DataDict['RemappedKey']
            
            # And the DICOMs that belong to that series:
            CopyFromDicoms = copy.deepcopy(RemappedDicoms)
            
            if Shift:
                # The key of the series to copy the ROIs to:
                CopyToKey = DataDict['ShiftedKey']
                
                # And the DICOMs that belong to that series:
                CopyToDicoms = copy.deepcopy(ShiftedDicoms)
                
            else:
                # The key of the series to copy the ROIs to:
                CopyToKey = 'Target'
                
                # And the DICOMs that belong to that series:
                CopyToDicoms = copy.deepcopy(TargetDicoms)
            
    
    print('\n\nCopying the ROI Collection from', \
          DataDict[CopyFromKey]['SeriesDesc'], 'to', \
          DataDict[CopyToKey]['SeriesDesc'], '...')
    
    """ COPY THE ROI COLLECTION FROM RemappedDicoms TO ShiftedDicoms: """
    CopyFromRois, CopyToRois, DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, CopyFromKey, CopyToKey)        
           
    # Copy the ROIs before ModifyRoi is run again:
    if 'Shifted' in CopyToKey:
        ShiftedRois = copy.deepcopy(CopyToRois)
    else:
        TargetRois = copy.deepcopy(CopyToRois)
    
    
    # Create the key for the registered series:
    #RegisteredKey = 'Reg' + CopyToKey
    
    
    return SourceDicoms, SourceRois, TargetDicoms, TargetRois, ShiftedDicoms, ShiftedRois, \
           RemappedDicoms, RemappedRois, RegDicoms, RegRois, DataDict
           
    
    
    """ April 2: Commented out for version below that handles more cases """
    if False:
        """ COPY THE ROI COLLECTION FROM SourceDicoms TO THE REGISTERED SERIES: """
        
        # Initialise variables that will be updated below:
        SourceRois = []
        TargetRois = []
        ShiftedRois = []
        RemappedRois = []
        
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
            # Define the keys within DataDict that determine which filepaths will
            # be used in ModifyRoi():
            FromKey = 'Source'
            #ToKey = 'RegMapPosCorr'
            ToKey = 'Registered'
            
            #print('\n\nCopying the ROI Collection from', \
            #      DataDict[FromKey]['DicomExpLabel'], 'to', \
            #      DataDict[ToKey]['DicomExpLabel'], '...')
            print('\n\nCopying the ROI Collection from', \
                  DataDict[FromKey]['SeriesDesc'], 'to', \
                  DataDict[ToKey]['SeriesDesc'], '...')
            
            #roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2)
            #roi1, roi2, xnatDict = ModifyRoi(xnatDict, dicomDict1_to_2) # 10/03/2020
            SourceRois, RegRois, DataDict = ModifyRoi(SourceDicoms, RegDicoms, DataDict, FromKey, ToKey)
            
            #return SourceDicoms, SourceRois, ToDicoms, ToRois, PosCorrDicoms, [], \
            #       MappedDicoms, [], RegDicoms, RegRois, DataDict
            
            # times[3] = time after parsing files, finding closest slice locations 
            # and modifying the ROIs:
            times.append(time.time())
            
            # The time taken to modify the ROIs:
            Dtime = round(times[-1] - times[-2], 1)
            
            print(f'\n   Took {Dtime} s to copy and modify the ROIs.')
        else:
            RegRois = []
        
        
        
        
    """ COPY THE ROI COLLECTION TO THE REGISTERED SERIES: """    
        
        
        
    
    
    # Do some housekeeping - delete the DICOM files in SourceDicomDir and 
    # TargetDicomDir:
    if DelDicoms:
        print('\nDeleting temporary DICOM files..')
        for directory in [DataDict['Source']['DicomDir'], DataDict['Target']['DicomDir']]:
            for root, dirs, files in os.walk(directory):
                for f in files:
                    os.unlink(os.path.join(root, f))
                for d in dirs:
                    shutil.rmtree(os.path.join(root, d))
                    
    
    # The time taken to do everything:
    Dtime = round(times[-1] - times[-2], 1)
    
    print(f'\n\nTook {Dtime} s to do everything.')
    
    #return roi2, xnatDict
    #return dicomDict1, roi1, dicomDict2, roi2, dicomDict1_to_2, xnatDict # 10/03/2020
    #return dicomDict1, roi1, dicomDict2posCorr, roi2, dicomDict1_to_2, xnatDict # 15/03/2020
    #return FromDicoms, FromRois, ToDicoms, ToRois, PosCorrDicoms, [], \
    #       MappedDicoms, [], RegDicoms, RegRois, DataDict
    #return FromDicoms, FromRois, ToDicoms, [], PosCorrDicoms, [], \
    #       MappedDicoms, MapPosCorrRois, RegDicoms, RegRois, DataDict
    #return FromDicoms, FromRois, ToDicoms, [], ShiftedDicoms, [], \
    #       MappedDicoms, MappedRois, RegDicoms, RegRois, DataDict
    return SourceDicoms, SourceRois, TargetDicoms, TargetRois, ShiftedDicoms, ShiftedRois, \
           RemappedDicoms, RemappedRois, RegDicoms, RegRois, DataDict
