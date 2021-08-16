# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:36:47 2020

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
                           
    ShiftFeature           - The feature that was used to manually determine
                           the required Shift (e.g. 'lens of eye', 
                           'peak of nose', 'butterfly-like structure in brain')
                        
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


    
def CopyRoi(XnatAddress, XnatUsername, XnatPassword, RootDir, CopyFrom, CopyTo, \
            SourceRoiProjLabel, SourceRoiSubjLabel, SourceRoiExpLabel, \
            SourceRoiLabelDateTime, SourceSeriesNo, SourceSeriesDesc, \
            TargetDicomProjLabel, TargetDicomSubjLabel, TargetDicomExpLabel, \
            TargetSeriesNo, TargetSeriesDesc, \
            SearchBy, Shift, ShiftFeature, RegMethod, DelDicoms, RerunReg, Debug):
    
    # IMPORT PACKAGES AND FUNCTIONS:
    import time, os, shutil
    import importlib
    import copy
    
    import CreateDir
    importlib.reload(CreateDir)
    from CreateDir import CreateDir
    
    import CreateDirsForDownloads
    importlib.reload(CreateDirsForDownloads)
    from CreateDirsForDownloads import CreateDirsForDownloads
    
    import DownloadFilesForRoiCopy
    importlib.reload(DownloadFilesForRoiCopy)
    from DownloadFilesForRoiCopy import DownloadFilesForRoiCopy
    
    import GetDicomFpaths
    importlib.reload(GetDicomFpaths)
    from GetDicomFpaths import GetDicomFpaths
    
    import GetDicoms
    importlib.reload(GetDicoms)
    from GetDicoms import GetDicoms
    
    import GetSopUids
    importlib.reload(GetSopUids)
    from GetSopUids import GetSopUids
    
    import GetScanPositions
    importlib.reload(GetScanPositions)
    from GetScanPositions import GetScanPositions
    
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
    
    import CheckDicomPixValAndDtypeCompatibility
    importlib.reload(CheckDicomPixValAndDtypeCompatibility)
    from CheckDicomPixValAndDtypeCompatibility import CheckDicomPixValAndDtypeCompatibility
    
    import ModifyRoi
    importlib.reload(ModifyRoi)
    from ModifyRoi import ModifyRoi
    
    import CheckRoiToDicomsCompatibility
    importlib.reload(CheckRoiToDicomsCompatibility)
    from CheckRoiToDicomsCompatibility import CheckRoiToDicomsCompatibility
    
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
    #RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel)
    RootDir = os.path.join(r'C:\Temp', todaysDate + ' ' + SourceRoiProjLabel \
                           + ' ' + CopyFrom + ' to ' + CopyTo)
    
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
                'XnatAddress':XnatAddress,\
                'XnatUsername':XnatUsername,\
                'XnatPassword':XnatPassword,\
                
                # Add the other variables:
                'RootDir':RootDir, \
                'Shift':Shift, \
                'ShiftFeature':ShiftFeature, \
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
    RemappedSourceRois = []
    RemappedTargetRois = []
    RemappedShiftedTargetRois = []
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
    

    
    
    """ GET THE DICOM OBJECTS AND FILEPATHS THAT CORRESPOND TO THE ROI
    COLLECTION TO BE COPIED: """
    #print('\n\nGetting the DICOM objects from', \
    #      DataDict['From']['RoiExpLabel'], '...')
    print('\n\nGetting the DICOM objects from', \
          DataDict['Source']['SeriesDesc'], '...')
    
    SourceDicomFpaths, SourceDicoms = GetDicoms(DataDict['Source']['DicomDir'], \
                                            SortMethod='slices')

    # ADD SourceDicomFpaths and SliceNos TO dataDict:
    DataDict['Source'].update({'DicomFpaths':SourceDicomFpaths})
    DataDict['Source'].update({'SliceNos':list(range(len(SourceDicomFpaths)))})
    
    
    
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
    
    #print('\n\nGetting the DICOM objects from', \
    #      DataDict['To']['DicomExpLabel'], '...')
    print('\n\nGetting the DICOM objects from', \
          DataDict['Target']['SeriesDesc'], '...')
    
    TargetDicomFpaths, TargetDicoms = GetDicoms(DataDict['Target']['DicomDir'], SortMethod='slices')
        
    # ADD TargetDicomFpaths and SliceNos TO dataDict:
    DataDict['Target'].update({'DicomFpaths':TargetDicomFpaths})
    DataDict['Target'].update({'SliceNos':list(range(len(TargetDicomFpaths)))})
                
    
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
    
   


    
    
    """ GET THE SOP INSTANCE UIDs AND SCAN POSITIONS OF THE SOURCE DICOMs: """
    
    print('\n\nGetting the SOP Instance UIDs and scan positions for', \
          DataDict['Source']['SeriesDesc'], '...')
    
    SourceSopUids = GetSopUids(SourceDicoms)
    
    SourceScanPos = GetScanPositions(SourceDicoms, SearchBy)
        
    # ADD SourceSopUids and SourceScanPositions TO dataDict:
    DataDict['Source'].update({'SopUids':SourceSopUids})
    DataDict['Source'].update({'ScanPos':SourceScanPos})
                

    
    
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
   

    
    
    """ GET THE SOP INSTANCE UIDs AND SCAN POSITIONS OF THE TARGET DICOMs: """
    
    print('\n\nGetting the SOP Instance UIDs and scan positions for', \
          DataDict['Target']['SeriesDesc'], '...')
    
    TargetSopUids = GetSopUids(TargetDicoms)
    
    TargetScanPos = GetScanPositions(TargetDicoms, SearchBy)
        
    # ADD TargetSopUids TO dataDict:
    DataDict['Target'].update({'SopUids':TargetSopUids})
    DataDict['Target'].update({'ScanPos':TargetScanPos})
                

    
    
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
    
    
    
    """ APPLY SCAN POSITIONAL SHIFT TO TargetDicoms: """
    if Shift:
        if Debug:
            print('\nKeys in DataDict prior to running CreateDirsForShifts:\n', \
                  DataDict.keys())
        
        """ CREATE DIRECTORIES FOR SHIFTED DICOMS AND ROIS: """  
        
        # The key of the series to be shifted:
        ToShiftKey = 'Target'
        
        # Store ToShiftKey in DataDict:
        DataDict.update({'ToShiftKey':ToShiftKey})
        
        DataDict = CreateDirsForShifts(DataDict)
        
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
              
         
        S = len(DataDict['Source']['ScanPos'])
        Tpre = len(DataDict['Target']['ScanPos'])
        Tpost = len(DataDict['ShiftedTarget']['ScanPos'])
        
        print(f'\nThe Source scan positions are (N = {S}):\n\n', \
              [round(n,2) for n in DataDict['Source']['ScanPos']])
            
        print(f'\nThe pre-shifted Target scan positions were (N = {Tpre}):\n\n', \
              [round(n,2) for n in DataDict['Target']['ScanPos']])
        
        print(f'\nThe post-shifted Target scan positions are (N = {Tpost}):\n\n', \
              [round(n,2) for n in DataDict['ShiftedTarget']['ScanPos']])
        
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
    
    print('\n\nGetting the indeces that map each slice in', \
          DataDict['Target']['DicomExpLabel'], '\nto each slice in', \
          DataDict['Source']['RoiExpLabel'], '...')
    
    if Shift:
        #print('\n\nGetting the indeces that map each slice in', \
        #      DataDict['ShiftedTarget']['DicomExpLabel'], 'to each slice in ', \
        #      DataDict['Source']['RoiExpLabel'], '...')
        
        DataDict = GetClosestSlices(SourceDicoms, ShiftedDicoms, DataDict)
        
    else:
        #print('\n\nGetting the indeces that map each slice in', \
        #      DataDict['Target']['DicomExpLabel'], 'to each slice in ', \
        #      DataDict['Source']['RoiExpLabel'], '...')
        
        DataDict = GetClosestSlices(SourceDicoms, TargetDicoms, DataDict)
    
    #return SourceDicoms, 0, TargetDicoms, 0, PosCorrDicoms, MappedDicoms, 0, DataDict
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
        
  
    
    """ USE THE MAPPING INDECES TO REMAP THE DICOMS OF GREATER LENGTH 
    AND DEFINE FixedKey, MovingKey, FixedDicoms and MovingDicom REQUIRED
    FOR IMAGE REGISTRATION: """

    
    """ Note:
        The indeces that map each DICOM in SourceDicoms to TargetDicoms, and
        from TargetDicoms to SourceDicoms are 
        DataDict['SourceToTargetNearest']['MappingInds'] and
        DataDict['TargetToSourceNearest']['MappingInds'].
        
        But the nearest are not necessarily close enough to be a relevant
        mapping.  So a subset of those indeces above was created by omitting
        any index for which the scan position difference was greater than some
        threshold distance (e.g. 1 mm).  Those corresponding indeces are
        DataDict['SourceToTargetNearOnly']['MappingInds'] and
        DataDict['TargetToSourceNearOnly']['MappingInds'].
        
        One of those two lists of indeces will be shorter and the longer list
        will need to be resliced (reduced in length) to match the length of the
        shorter one.
    """
    
    #S = len(SourceDicoms)
    #T = len(TargetDicoms)
    S2T = len(DataDict['SourceToTargetNearOnly']['MappingInds'])
    T2S = len(DataDict['TargetToSourceNearOnly']['MappingInds'])
    
    """ April 27th revision to code below: """
    
    if S2T != T2S:
        print('\nThere are an unequal number of indeces in', \
              f'SourceToTargetNearestInds (N = {S2T}) than in', \
              f'TargetToSourceNearestInds (N = {T2S}).  Both Source and', \
              'Target DICOMs will be remapped.')
        
        """ FIRST REMAP SourceDicoms: """
        KeyToRemap = 'Source'
        DicomsToRemap = copy.deepcopy(SourceDicoms)
        
        # Create directories for remapped KeyToRemap:   
        DataDict = CreateDirsForRemap(DataDict, KeyToRemap)
            
        # Remap DicomsToRemap:    
        RemappedSourceDicoms, DataDict = RemapDicoms(DicomsToRemap, DataDict, \
                                                     KeyToRemap)
    
        """ NEXT REMAP TargetDicoms or ShiftedDicoms: """
        if Shift:
            KeyToRemap = DataDict['ShiftedKey']
            DicomsToRemap = copy.deepcopy(ShiftedDicoms)
        else:
            KeyToRemap = 'Target'
            DicomsToRemap = copy.deepcopy(TargetDicoms)
    
        # Create directories for remapped KeyToRemap:   
        DataDict = CreateDirsForRemap(DataDict, KeyToRemap)
            
        # Remap DicomsToRemap: 
        if Shift:
            RemappedShiftedTargetDicoms, DataDict = RemapDicoms(DicomsToRemap, \
                                                                DataDict, \
                                                                KeyToRemap)
        else:
            RemappedTargetDicoms, DataDict = RemapDicoms(DicomsToRemap, \
                                                         DataDict, \
                                                         KeyToRemap)
        
        
        """ DEFINE FixedKey, MovingKey, FixedDicoms AND MovingDicoms: """
        # The key of the series to be fixed is the remapped Source DICOMs:
        FixedKey = 'RemappedSource'
            
        # And the corresponding DICOMs:
        FixedDicoms = copy.deepcopy(RemappedSourceDicoms)
        
        # The key of the series to be registered to the fixed series and
        # the corresponding DICOMs:
        if Shift:
            MovingKey = 'ShiftedKey'
            MovingDicoms = copy.deepcopy(RemappedShiftedTargetDicoms)
        else:
            MovingKey = 'Target'
            MovingDicoms = copy.deepcopy(RemappedTargetDicoms)
            
        
            
        print('\nFixedKey = ' + FixedKey + ' and MovingKey = ' + MovingKey + '.')     
        
        
        
        """ PRINT FOR DEBUGGING: """
        
        # The length of pre- and post-remapped SourceDicoms are:
        SpreR = len(SourceDicoms)
        SpostR = len(RemappedSourceDicoms)
        
        # The length of pre- and post-remapped Target/ShiftedDicoms are:
        if Shift:
            TpreR = len(ShiftedDicoms)
            TpostR = len(RemappedShiftedTargetDicoms)
        else:
            TpreR = len(TargetDicoms)
            TpostR = len(RemappedTargetDicoms)
        
        
        
        """ PRINT FOR DEBUGGING: """
        print(f'\nThe pre-remapped SourceDicoms had {SpreR} slices with', \
              'scan positional differences of:\n\n', \
              [round(n,2) for n in DataDict['SourceToTargetNearest']['ScanPosDiffs']])
        
        print(f'\nThe post-remapped SourceDicoms has {SpostR} slices with', \
              'scan positional differences of:\n\n', \
              [round(n,2) for n in DataDict['RemappedSource']['ScanPosDiffs']])
        
        print(f'\nThe pre-remapped TargetDicoms had {TpreR} slices with', \
              'scan positional differences of:\n\n', \
              [round(n,2) for n in DataDict['TargetToSourceNearest']['ScanPosDiffs']])
        
        if Shift:
            print(f'\nThe post-remapped TargetDicoms has {TpostR} slices with', \
                  'scan positional differences of:\n\n', \
                  [round(n,2) for n in DataDict['RemappedShiftedTarget']['ScanPosDiffs']])
        else:
            print(f'\nThe post-remapped TargetDicoms has {TpostR} slices with', \
                  'scan positional differences of:\n\n', \
                  [round(n,2) for n in DataDict['RemappedTarget']['ScanPosDiffs']])
        

                
    else: # i.e. S2T == T2S
        #print('\nThere are equal numbers of slices in SourceDicoms and', \
        #      'TargetDicoms so no DICOMs will be remapped.')
        print('\nThere are an equal number of indeces in', \
              f'SourceToTargetNearestInds (N = {S2T}) and', \
              f'TargetToSourceNearestInds (N = {T2S}) so no remapping required.')
        
        
        """ DEFINE FixedKey, MovingKey, FixedDicoms AND MovingDicoms: """
        # The key of the series to be fixed:
        FixedKey = 'Source'
        
        # And the corresponding DICOMs:
        FixedDicoms = copy.deepcopy(SourceDicoms)
        
        # The key of the series to be registered to the fixed series and 
        # the corresponding DICOMs:
        if Shift:
            MovingKey = DataDict['ShiftedKey']
            
            MovingDicoms = copy.deepcopy(ShiftedDicoms)
        else:
            MovingKey = 'Target'
            
            MovingDicoms = copy.deepcopy(TargetDicoms)
        
        
    
    """ If statement below was pre-April 27: """
    if False:
        if S2T != T2S:
            if S2T > T2S:
                # The key of the series to be remapped:
                KeyToRemap = 'Source'
                    
                # The key of the series that KeyToRemap maps to:
                if Shift:
                    KeyForRemap = DataDict['ShiftedKey']
                else:
                    KeyForRemap = 'Target'
                    
                # Store ToRemapKey and ForRemapKey in DataDict:
                DataDict.update({'KeyToRemap':KeyToRemap})
                DataDict.update({'KeyForRemap':KeyForRemap})
                
                #print(f'\nThere are more slices in SourceDicoms (N = {S})', \
                #      f'than in TargetDicoms (N = {T}) so', ToRemapKey, \
                #      'will be remapped.')
                
                print('\nThere are more indeces in SourceToTargetNearestInds',\
                      f'(N = {S2T}) than in TargetToSourceNearestInds', \
                      f'(N = {T2S})\nso', KeyToRemap, 'will be remapped.')
                
                """ CREATE DIRECTORIES FOR REMAPPED DICOMS AND ROIS: """   
                DataDict = CreateDirsForRemap(DataDict, KeyToRemap)
                    
                    
                """ RE-MAP SourceDicoms: """
                RemappedDicoms, DataDict = RemapDicoms(SourceDicoms, DataDict, \
                                                       KeyToRemap)
                
                
                """ DEFINE FixedKey, MovingKey, FixedDicoms AND MovingDicoms: """
                # The key of the series to be fixed is the remapped Source DICOMs:
                FixedKey = DataDict['RemappedKey']
                    
                # And the corresponding DICOMs:
                FixedDicoms = copy.deepcopy(RemappedDicoms)
                
                # The key of the series to be registered to the fixed series and
                # the corresponding DICOMs:
                if Shift:
                    MovingKey = DataDict['ShiftedKey']
                    
                    MovingDicoms = copy.deepcopy(ShiftedDicoms)
                else:
                    MovingKey = 'Target'
                    
                    #MovingDicoms = copy.deepcopy(MovingDicoms) # <- was before April 24
                    MovingDicoms = copy.deepcopy(TargetDicoms) # correction
                    
                print('\nFixedKey = ' + FixedKey + ' and MovingKey = ' + MovingKey + '.')     
                
            
            elif S2T < T2S:
                # The key of the series that KeyToRemap maps to:
                KeyForRemap = 'Source'
                
                # The key of the series to be remapped:
                if Shift:  
                    KeyToRemap = DataDict['ShiftedKey']
                else:
                    KeyToRemap = 'Target'    
                    
                # Store KeyToRemap and KeyForRemap in DataDict:
                DataDict.update({'KeyToRemap':KeyToRemap})
                DataDict.update({'KeyForRemap':KeyForRemap})
                
                #print(f'\nThere are fewer slices in SourceDicoms (N = {S})', \
                #      f'than in TargetDicoms (N = {T}) so', ToRemapKey, \
                #      'will be remapped.')
                
                print('\nThere are fewer indeces in SourceToTargetNearestInds',\
                      f'(N = {S2T}) than in TargetToSourceNearestInds', \
                      f'(N = {T2S}) \nso', KeyToRemap, 'will be remapped.')
                
                """ CREATE DIRECTORIES FOR REMAPPED DICOMS AND ROIS: """   
                DataDict = CreateDirsForRemap(DataDict, KeyToRemap)
                
                
                if Shift:
                    """ RE-MAP ShiftedTargetDicoms: """
                    RemappedDicoms, DataDict = RemapDicoms(ShiftedDicoms, \
                                                           DataDict, \
                                                           KeyToRemap)
                else:
                    """ RE-MAP TargetDicoms: """
                    RemappedDicoms, DataDict = RemapDicoms(TargetDicoms, \
                                                           DataDict, \
                                                           KeyToRemap)
                    
                
                """ DEFINE FixedKey, MovingKey, FixedDicoms AND MovingDicoms: """
                # The key of the series to be fixed:
                FixedKey = 'Source'
                
                # And the corresponding DICOMs:
                FixedDicoms = copy.deepcopy(SourceDicoms)
                    
                # The key of the series to be registered to the fixed series:
                MovingKey = DataDict['RemappedKey']
                
                # And the corresponding DICOMs:
                MovingDicoms = copy.deepcopy(RemappedDicoms)
    
              
            # Print for debugging:
            N = len(DataDict[KeyToRemap]['ScanPos'])
            
            print('\nThe pre-remapped', KeyToRemap, \
                  f'scan positions were (N = {N}):\n\n', \
                  DataDict[KeyToRemap]['ScanPos'])
            
            N = len(DataDict[DataDict['RemappedKey']]['ScanPos'])
            
            print('\nThe post-remapped', KeyToRemap, \
                  f'scan positions are (N = {N}):\n\n', \
                  DataDict[DataDict['RemappedKey']]['ScanPos'])
            
            N = len(DataDict[KeyForRemap]['ScanPos'])
            
            print('\nThe scan positions of the DICOMs that the remapping was', \
                  'for (' + KeyForRemap + f') are (N = {N}):\n\n', \
                  DataDict[KeyForRemap]['ScanPos'])
                    
        else: # i.e. S2T == T2S
            #print('\nThere are equal numbers of slices in SourceDicoms and', \
            #      'TargetDicoms so no DICOMs will be remapped.')
            print('\nThere are an equal number of indeces in', \
                  f'SourceToTargetNearestInds (N = {S2T}) and', \
                  f'TargetToSourceNearestInds (N = {T2S}) so no remapping required.')
            
            
            """ DEFINE FixedKey, MovingKey, FixedDicoms AND MovingDicoms: """
            # The key of the series to be fixed:
            FixedKey = 'Source'
            
            # And the corresponding DICOMs:
            FixedDicoms = copy.deepcopy(SourceDicoms)
            
            # The key of the series to be registered to the fixed series and 
            # the corresponding DICOMs:
            if Shift:
                MovingKey = DataDict['ShiftedKey']
                
                MovingDicoms = copy.deepcopy(ShiftedDicoms)
            else:
                MovingKey = 'Target'
                
                MovingDicoms = copy.deepcopy(TargetDicoms)
    
    
    
    
    print('\nFixedKey = ' + FixedKey + f'. \nFixedDicoms has {len(FixedDicoms)}' \
          ' slices. \nMovingKey = ' + MovingKey + '. \nMovingDicoms has ' \
          f'{len(MovingDicoms)} slices.') 
    
    # Define the key to be assigned to the registered series:
    RegisteredKey = 'Registered'
    
    # Store FixedKey, MovingKey and RegisteredKey in DataDict:
    DataDict.update({'FixedKey':FixedKey})
    DataDict.update({'MovingKey':MovingKey})
    DataDict.update({'RegisteredKey':RegisteredKey})
    
    #print(f'\nlen(FromDicoms) = {len(FromDicoms)}')        
        
    
    #return FromDicoms, 0, ToDicoms, 0, ShiftedDicoms, 0, \
    #       RemappedDicoms, 0, 0, 0, DataDict
    

    """ APRIL 23:  Return before Registration for debugging: """
    if Shift:
        return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
               ShiftedDicoms, ShiftedRois, \
               RemappedSourceDicoms, RemappedSourceRois, \
               RemappedShiftedTargetDicoms, RemappedShiftedTargetRois, \
               FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
               RegisteredDicoms, RegisteredRois, DataDict
    else:
        return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
               ShiftedDicoms, ShiftedRois, \
               RemappedSourceDicoms, RemappedSourceRois, \
               RemappedTargetDicoms, RemappedTargetRois, \
               FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
               RegisteredDicoms, RegisteredRois, DataDict
           
    

       
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    

    
    """ APRIL 24:  CAN REMOVE ALL OF THE FOLLOWING SECTION IN THE if False
    STATEMENT BELOW (SINCE THIS WAS ALL ADDED TO THE SECTION ABOVE): """
    if False:
        """ DEFINE FixedDicoms and MovingDicoms FOR REGISTRATION: """
        
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
        
            # The key name of the registered series:
            RegisteredKey = 'Registered'
        
            # Define FixedKey, MovingKey, FixedDicoms and MovingDicoms:
            if len(SourceDicoms) != len(TargetDicoms):
                if len(SourceDicoms) > len(TargetDicoms):
                    # The key of the series to be fixed:
                    FixedKey = DataDict['RemappedKey']
                    
                    # And the corresponding DICOMs:
                    FixedDicoms = copy.deepcopy(RemappedDicoms)
                    
                    if Shift:
                        # The key of the series to be registered to the fixed series:
                        MovingKey = DataDict['ShiftedKey']
                        
                        # And the corresponding DICOMs:
                        MovingDicoms = copy.deepcopy(ShiftedDicoms)
                    else:
                        # The key of the series to be registered to the fixed series:
                        MovingKey = 'Target'
                        
                        # And the corresponding DICOMs:
                        MovingDicoms = copy.deepcopy(MovingDicoms)
                        
                    print('\nThere are more slices in SourceDicoms than in ' \
                          + 'TargetDicoms, \nso FixedKey = ' + FixedKey \
                          + ' and MovingKey = ' + MovingKey + '.')
                
                elif len(SourceDicoms) < len(TargetDicoms):
                    # The key of the series to be fixed:
                    FixedKey = 'Source'
                    
                    # And the corresponding DICOMs:
                    FixedDicoms = copy.deepcopy(SourceDicoms)
                        
                    # The key of the series to be registered to the fixed series:
                    MovingKey = DataDict['RemappedKey']
                    
                    # And the corresponding DICOMs:
                    MovingDicoms = copy.deepcopy(RemappedDicoms)
                        
                    print('\nThere are fewer slices in SourceDicoms than in ' \
                      + 'TargetDicoms, \nso FixedKey = ' + FixedKey \
                      + ', and MovingKey = ' + MovingKey + '.')
                    
            else:
                # The key of the series to be fixed:
                FixedKey = 'Source'
                
                # And the corresponding DICOMs:
                FixedDicoms = copy.deepcopy(SourceDicoms)
                
                if Shift:
                    # The key of the series to be registered to the fixed series:
                    MovingKey = DataDict['ShiftedKey']
                    
                    # And the corresponding DICOMs:
                    MovingDicoms = copy.deepcopy(ShiftedDicoms)
                else:
                    # The key of the series to be registered to the fixed series:
                    MovingKey = 'Target'
                    
                    # And the corresponding DICOMs:
                    MovingDicoms = copy.deepcopy(TargetDicoms)
                    
                print('\nThere are equal numbers of slices in SourceDicoms and ' \
                      + 'TargetDicoms, so FixedKey = ' + FixedKey \
                      + ', and MovingKey = ' + MovingKey + '.')
                
                
            # Store FixedKey, MovingKey and RegisteredKey in DataDict:
            DataDict.update({'FixedKey':FixedKey})
            DataDict.update({'MovingKey':MovingKey})
            DataDict.update({'RegisteredKey':RegisteredKey})
                
            """ CREATE DIRECTORIES FOR REGISTERED DICOMS AND ROIS: """  
            DataDict = CreateDirsForReg(DataDict)
            
            
            print('\nKeys in DataDict after running CreateDirsForReg:\n\n', \
                      DataDict.keys())
            
            
            if False:
                """ GET THE POSITIONAL DIFFERENCES BETWEEN THE SLICES IN FixedDicoms
                AND MovingDicoms: """
                
                N = len(FixedDicoms)
                
                dScanPos = [DataDict[MovingKey]['ScanPos'][i] 
                            - DataDict[FixedKey]['ScanPos'][i] for i in range(N)]
                            
                
                # Store dScanPos in DataDict:
                #DataDict.update({'FixedToMovingScanPosDiffs':dScanPos}) 
                
                
                """ RESLICE FixedDicoms and MovingDicoms, OMITTING SLICES WHERE THE 
                SCAN POSITIONAL DIFFERENCES ARE GREATER THAN SOME THRESHOLD DISTANCE: 
                """
                
                ThreshDist = 1 # mm
                
                FixedDicoms = [FixedDicoms[i] for i in range(N) \
                               if abs(dScanPos[i]) < ThreshDist]
                
                MovingDicoms = [MovingDicoms[i] for i in range(N) \
                               if abs(dScanPos[i]) < ThreshDist]
                
                # Finally reslice dScanPos:
                dScanPos = [dScanPos[i] for i in range(N) if abs(dScanPos[i]) < ThreshDist]
                
                # Store dScanPos in DataDict:
                DataDict.update({'FixedToMovingScanPosDiffs':dScanPos}) 
            
            
            """ APRIL 23:  Return before Registration for debugging: """
            return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
               ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
               FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
               RegisteredDicoms, RegisteredRois, DataDict 
     
        """ END OF if False (section above can be deleted) """
        
        
        
    """ PERFORM IMAGE REGISTRATION: """
    
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        
        """ CREATE DIRECTORIES FOR REGISTERED DICOMS AND ROIS: """  
        DataDict = CreateDirsForReg(DataDict)
            
            
        #print('\n\nRegistering images in', DataDict[MovingKey]['DicomExpLabel'], \
        #       'to', DataDict[FixedKey]['DicomExpLabel'], '...')
        print('\n\nRegistering images in', DataDict[MovingKey]['SeriesDesc'], \
               'to', DataDict[FixedKey]['SeriesDesc'], '...')
        
        FixedDicoms, \
        MovingDicoms, \
        RegisteredDicoms, \
        DataDict = RegisterSlices(DataDict, FixedKey, MovingKey, RegisteredKey)
        
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
        



        """ 
        CHECK THE COMPATIBILITY BETWEEN THE REGISTERED DICOM'S PIXEL VALUES AND  
        THE CORRESPONDING DATATYPE:
        """
        
        print('\nChecking the compatibility between the registered DICOM\'s', \
              'pixel values and its data type:\n')
        
        # Define the directories:
        FixedDir = DataDict[FixedKey]['DicomDir']
        MovingDir = DataDict[MovingKey]['DicomDir']
        RegDir = DataDict[RegisteredKey]['DicomDir']
        
        # Check the ROI and DICOMs for each key in CheckKeys:
        CheckDicomPixValAndDtypeCompatibility(FixedDir, MovingDir, RegDir)
        
        
        
        
        print('\n***********************************************************' \
              + '************************************************************' \
              '\n***********************************************************' \
              + '************************************************************' \
              '\n***********************************************************' \
              + '************************************************************')
    

    
    
    """ COPY THE ROI COLLECTION FROM SOURCE TO VARIOUS TARGET SERIES: """
            
    # The key of the series whose ROIs are to be copied:
    CopyFromKey = 'Source'
    CopyFromDicoms = copy.deepcopy(SourceDicoms)
    
    """ Copy ROIs for Target DICOMs: """
    # The key of the series to copy the ROIs to:
    CopyToKey = 'Target'
    CopyToDicoms = copy.deepcopy(TargetDicoms)
        
    print('\nROIs will be copied from ' + CopyFromKey + ' to ' + CopyToKey \
          + '...')
    
    """ COPY THE ROI COLLECTION: """
    SourceRois, \
    TargetRois, \
    DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, \
                         CopyFromKey, CopyToKey)  
    
    # Get list of keys in DataDict:
    keys = list(DataDict.keys())
    
    """ If 'ShiftedTarget' exists, copy for those DICOMs: """
    if 'ShiftedTarget' in keys:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['ShiftedKey']
        CopyToDicoms = copy.deepcopy(ShiftedDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, \
        ShiftedRois, \
        DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, \
                             CopyFromKey, CopyToKey)  
        
    """ If 'RemappedSource' exists, copy for those DICOMs: """
    if 'RemappedSource' in keys:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['RemappedKey']
        CopyToDicoms = copy.deepcopy(RemappedDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, \
        RemappedRois, \
        DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, \
                             CopyFromKey, CopyToKey)
        
    """ If 'RemappedShiftedTarget' exists, copy for those DICOMs: """
    if 'RemappedShiftedTarget' in keys:
        # The key of the series to copy the ROIs to:
        CopyToKey = DataDict['RemappedKey']
        CopyToDicoms = copy.deepcopy(RemappedDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, \
        RemappedRois, \
        DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, \
                             CopyFromKey, CopyToKey)
        
    """ If 'Registered' exists, copy for those DICOMs: """
    if 'Registered' in keys:
        # Copy the ROIs for the RegisteredDicoms:
        CopyToKey = DataDict['RegisteredKey']
        CopyToDicoms = copy.deepcopy(RegisteredDicoms)
        
        print('\nROIs will be copied from', CopyFromKey, 'to', CopyToKey, '.')
        
        """ COPY THE ROI COLLECTION: """
        SourceRois, \
        RegisteredRois, \
        DataDict = ModifyRoi(CopyFromDicoms, CopyToDicoms, DataDict, \
                             CopyFromKey, CopyToKey)
        
    
    
    
    """ April 7:  This is fine up to here: """
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
    
    
    
    """ 
    CHECK THE COMPATIBILITY BETWEEN THE ROI COLLECTIONS AND THE 
    CORRESPONDING DICOMS:
    """
    
    if True: 
        print('\nChecking the compatibility between the DICOMs and the ROI', \
              'Collections that were created for them:\n')
        
        # Define the keys in DataDict to check:
        CheckKeys = ['Source', 'Target', DataDict['ShiftedKey'], \
                     DataDict['RemappedKey'], 'Registered']
        
        # Check the ROI and DICOMs for each key in CheckKeys:
        for CheckKey in CheckKeys:
            print('\n\nChecking', CheckKey, 'data...')
        
            CheckRoiToDicomsCompatibility(DataDict, CheckKey)
    
    
    
    
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
      
    
    """ DEFINE FixedRois AND MovingRois FOR PLOTTING PURPOSES: """
    
    
    # Define FixedRois and MovingRois:
    if len(SourceDicoms) != len(TargetDicoms):
        if len(SourceDicoms) > len(TargetDicoms):
            FixedRois = copy.deepcopy(RemappedRois)
            
            if Shift:
                MovingRois = copy.deepcopy(ShiftedRois)
            else:
                MovingRois = copy.deepcopy(MovingRois)
        
        elif len(SourceDicoms) < len(TargetDicoms):
            FixedRois = copy.deepcopy(SourceRois)
                
            MovingRois = copy.deepcopy(RemappedRois)
            
    else:
        FixedRois = copy.deepcopy(SourceRois)
        
        if Shift:
            MovingRois = copy.deepcopy(ShiftedRois)
        else:
            MovingRois = copy.deepcopy(TargetRois)
        

           
        
    
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
            
    

               
               
        
        
        print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
        
        
     
        
    """ Plot registration results: """    
      
    """ April 8:  I'll do this outside of this function. """
    
    if False:
        if RegMethod in ['rigid', 'affine', 'non-rigid']:
                
            print('\nPlotting registration results:')
            
            # Plot the registration results:            
            DataDict = PlotRegResults(FixedDicoms, MovingDicoms, RegisteredDicoms,\
                                      FixedRois, MovingRois, RegisteredRois,\
                                      DataDict)
        



            
    print('\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************' \
          '\n***********************************************************' \
          + '************************************************************')
    
    
        
        
    
    
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
    
 
    return SourceDicoms, SourceRois, TargetDicoms, TargetRois, \
           ShiftedDicoms, ShiftedRois, RemappedDicoms, RemappedRois, \
           FixedDicoms, FixedRois, MovingDicoms, MovingRois, \
           RegisteredDicoms, RegisteredRois, DataDict 
