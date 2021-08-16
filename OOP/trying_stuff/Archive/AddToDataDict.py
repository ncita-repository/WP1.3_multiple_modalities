# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:08:26 2020

@author: ctorti
"""





def AddToDataDict(DataDict):
    # Import packages:
    import os
    
    
    # Get some necessary info from dataDict:
    RootDir = DataDict['RootDir']
    #FromScanLabel = DataDict['From']['ScanLabel']
    #ToScanLabel = DataDict['To']['ScanLabel']
    FromSeriesNo = DataDict['From']['SeriesNo']
    ToSeriesNo = DataDict['To']['SeriesNo']
    Shift = DataDict['Shift']
    #searchBy = dataDict['SearchBy']
    RegMethod = DataDict['RegMethod']
    Debug = DataDict['Debug']
    
    
    """
    DEFINE THE SCAN LABELS AND DIRECTORY NAMES:
    """
    
    """
    For the DICOMs and DICOM-RTSTRUCT of the scan whose ROIs are to be copied:
    """
    
    # Where the DICOM-RTSTRUCT to be copied is to be saved to:
    #FromRoiDir = os.path.join(RootDir, FromScanLabel + ' ROIs')
    FromRoiDir = os.path.join(RootDir, FromSeriesNo + ' ROIs')
    
    # Where the DICOM scan that relate to the ROIs to be copied are to be saved 
    # to:
    #FromDicomDir = os.path.join(RootDir, FromScanLabel + ' DICOMs')
    FromDicomDir = os.path.join(RootDir, FromSeriesNo + ' DICOMs')
    
    # Add the above directories to dataDict:
    DataDict['From'].update({'RoiDir':FromRoiDir,\
                             'DicomDir':FromDicomDir,\
                             #'ScanLabel':FromScanLabel,\
                             'SeriesNo':FromSeriesNo,\
                             'SeriesDesc':'',\
                             'SliceNos':[]
                             })
    
    """
    Note:  The entry 'SeriesDesc' is a place-holder that will be added to all
    of the main keys in DataDict.  
    It will later be populated with the DICOMs' Series Description. 
    """
    
    
    """ 
    For the DICOMs and DICOM-RTSTRUCT of the scan that the ROIs are to be 
    copied to:
    """
    
    # The DICOMs of the series that the ROIs are to be copied to (prior to any
    # change to Image Position or image registration):
    #ToDicomDir = os.path.join(RootDir, ToScanLabel + ' DICOMs')
    ToDicomDir = os.path.join(RootDir, ToSeriesNo + ' DICOMs')
    
    # The new (modified copy) ROIs are to be stored at:
    #ToRoiDir = os.path.join(RootDir, ToScanLabel + ' ROIs')
    ToRoiDir = os.path.join(RootDir, ToSeriesNo + ' ROIs')
    
    
    # Add the above directories to dataDict:
    DataDict['To'].update({#'ScanLabel':ToScanLabel,\
                           'RoiDir':ToRoiDir,\
                           'DicomDir':ToDicomDir,\
                           'SeriesNo':ToSeriesNo,\
                           'SeriesDesc':'',\
                           'SliceNos':[]
                           })
    
        
    
    
    # This was moved down:
    #dataDict.update({'Mapped':{'DicomScanLab':toMapDicomScanLab,\
    #                           'DicomDir':toMapDicomDir,\
    #                           'RoiDir':toMapRoiDir
    #                           }
    #                 })
    
    """
    For efficiency, create a list of all directories that need to be created,
    AllDirs, and add the directories defined above as standard.  More
    directories will be appended within the if statements below.
    """
    #AllDirs = [RootDir, FromRoiDir, FromDicomDir, ToRoiDir, ToDicomDir, \
    #           MapRoiDir, MapDicomDir]
    AllDirs = [RootDir, FromRoiDir, FromDicomDir, ToRoiDir, ToDicomDir]
    
    
    """
    For the DICOMs that result from applying a "scan-position-shift" to 
    the Image Positions of the DICOMs in ToDicoms and ToMapDicoms.
    """
    
    if Shift:
        ##toPosCorrScanNo = toDicomScanNo + 100
        ##toPosCorrDicomScanNo = toDicomScanNo*10
        ##toPosCorrDicomScanLab = str(toPosCorrDicomScanNo)
        #toPosCorrDicomScanLab = toDicomScanLab + a + '1' + '0' + '0'
        ##toPosCorrDicomScanNo = int(toPosCorrDicomScanLab)
        #toMapPosCorrDicomScanLab = toDicomScanLab + a + '1' + '0' + '0'
        
        #toPosCorrDirName = toDirName + ' ' + searchBy + '-corr'
        
        
        # Position-corrected:
        #PosCorrScanLabel = ToScanLabel + '_posCorr'
        ShiftedSeriesNo = ToSeriesNo + '0' + '1'
        
        if Debug:
            print('\nShiftedSeriesNo =', ShiftedSeriesNo)
        
        """
        Note: 
            Format of Series Numbers will be the ToSeriesNo + '0' + a number of
            binaries which indicate whether position-correction was done,
            mapping, registration, etc.
        """
      
        #ShiftedDicomDir = os.path.join(RootDir, PosCorrScanLabel + ' DICOMs')
        #ShiftedRoiDir = os.path.join(RootDir, PosCorrScanLabel + ' ROIs')
        ShiftedDicomDir = os.path.join(RootDir, ShiftedSeriesNo + ' DICOMs')
        ShiftedRoiDir = os.path.join(RootDir, ShiftedSeriesNo + ' ROIs')
        
        # Append directories:
        #AllDirs.append(ShiftedRoiDir)
        #AllDirs.append(ShiftedDicomDir)
        
        
        # Mapped position-corrected:
        #MapPosCorrScanLabel = ToScanLabel + '_posCorr_mapped-to-' \
        #                      + FromScanLabel
        MappedSeriesNo = ShiftedSeriesNo + '1'
        
        if Debug:
            print('\nMappedSeriesNo =', MappedSeriesNo)
         
        #MapShiftedDicomDir = os.path.join(RootDir, MapPosCorrScanLabel + ' DICOMs')
        #MapShiftedRoiDir = os.path.join(RootDir, MapPosCorrScanLabel + ' ROIs')
        MappedDicomDir = os.path.join(RootDir, MappedSeriesNo + ' DICOMs')
        MappedRoiDir = os.path.join(RootDir, MappedSeriesNo + ' ROIs')
        
        # Append directories:
        #AllDirs.append(MapShiftedRoiDir)
        #AllDirs.append(MapShiftedDicomDir)
        
        
        # Add the above directories to dataDict:                     
        DataDict.update({'Shifted':{#'ScanLabel':PosCorrScanLabel,\
                                    'RoiDir':ShiftedRoiDir,\
                                    'DicomDir':ShiftedDicomDir,\
                                    'SeriesNo':ShiftedSeriesNo,\
                                    'SeriesDesc':'',\
                                    'SliceNos':[]
                                     }
                         })
    

        DataDict.update({'Mapped':{#'ScanLabel':MapPosCorrScanLabel,\
                                   'RoiDir':MappedRoiDir,\
                                   'DicomDir':MappedDicomDir,\
                                   'SeriesNo':MappedSeriesNo,\
                                   'SeriesDesc':'',\
                                   'SliceNos':[]
                                   }
                         })
    
    else:
        """
        For the DICOMs that result from mapping the DICOMs in ToDicoms to the
        DICOMs in FromDicoms.
        """
    
        """
        Notes:
        The slices in the "from" and "to" scans were not necessarily taken in 
        the same location. If so, the slices in the "to" series will need to be 
        matched to the slices in the "from" scan. Furthermore, there may be a
        different number of slices in the two scans.
        
        At this stage it's not known whether a mapping will need to be done, so
        each definition of ScanLab, DirName, DicomDir and RoiDir for registered
        data for the case of Shift or not will be duplicated with the
        addition of "Map" to allow the the possibility that DICOMs will need to
        be mapped.
        """
        
        #toMapDicomScanLab = toDicomScanLabel + a + '0' + '1' + '0'
        #MapScanLabel = ToScanLabel + '_mapped-to-' + FromScanLabel
        MappedSeriesNo = ToSeriesNo + '0' + '0' + '1'
        
        if Debug:
            print('\nMappedSeriesNo =', MappedSeriesNo)
        
        #MapDicomDir = os.path.join(RootDir, MapScanLabel + ' DICOMs')
        #MapRoiDir = os.path.join(RootDir, MapScanLabel + ' ROIs')
        MappedDicomDir = os.path.join(RootDir, MappedSeriesNo + ' DICOMs')
        MappedRoiDir = os.path.join(RootDir, MappedSeriesNo + ' ROIs')
        
        DataDict.update({'Mapped':{#'ScanLabel':MapScanLabel,\
                                   'RoiDir':MappedRoiDir,\
                                   'DicomDir':MappedDicomDir,\
                                   'SeriesNo':MappedSeriesNo,\
                                   'SeriesDesc':'',\
                                   'SliceNos':[]
                                   }
                         })
    
        
        
        
    """
    For the DICOMs that result from applying image registration to DICOMs that
    have either had scan-position-correction or not, and mapping from the 
    DICOMs in ToDicoms to those in FromDicoms.
    """     
    
    # The scan no and scan label for registered-only and for position-corrected 
    # and registered data:
    if RegMethod in ['rigid', 'affine', 'non-rigid']:
        if Shift:
            # Registered position-corrected:
            #RegPosCorrScanLabel = ToScanLabel + '_posCorr_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegShiftedSeriesNo = ShiftedSeriesNo + '1'
            elif RegMethod=='affine':
                RegShiftedSeriesNo = ShiftedSeriesNo + '2'
            elif RegMethod=='non-rigid':
                RegShiftedSeriesNo = ShiftedSeriesNo + '3'
                
                                                              
            #RegShiftedDicomDir = os.path.join(RootDir, RegPosCorrScanLabel + ' DICOMs')
            #RegShiftedRoiDir = os.path.join(RootDir, RegPosCorrScanLabel + ' ROIs')
            RegShiftedDicomDir = os.path.join(RootDir, RegShiftedSeriesNo + ' DICOMs')
            RegShiftedRoiDir = os.path.join(RootDir, RegShiftedSeriesNo + ' ROIs')
            
            # Append directories:
            #AllDirs.append(RegShiftedRoiDir)
            #AllDirs.append(RegShiftedDicomDir)
            
                            
            # Registered mapped position-corrected:
            #RegMapPosCorrScanLabel = ToScanLabel + '_posCorr_mapped-to-' \
            #                         + FromScanLabel + '_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegMappedSeriesNo = MappedSeriesNo + '1'
            elif RegMethod=='affine':
                RegMappedSeriesNo = MappedSeriesNo + '2'
            elif RegMethod=='non-rigid':
                RegMappedSeriesNo = MappedSeriesNo + '3'
            
            if Debug:
                print('\nRegMappedSeriesNo =', RegMappedSeriesNo)
            
            #RegMapShiftedDicomDir = os.path.join(RootDir, RegMapPosCorrScanLabel + ' DICOMs')
            #RegMapShiftedRoiDir = os.path.join(RootDir, RegMapPosCorrScanLabel + ' ROIs')
            RegMappedDicomDir = os.path.join(RootDir, RegMappedSeriesNo + ' DICOMs')
            RegMappedRoiDir = os.path.join(RootDir, RegMappedSeriesNo + ' ROIs')
            
            # Append directories:
            #AllDirs.append(RegMapShiftedRoiDir)
            #AllDirs.append(RegMapShiftedDicomDir)
            
            
            # Add the above directories to dataDict:                     
            DataDict.update({'RegShifted':{#'ScanLabel':RegPosCorrScanLabel,\
                                           'RoiDir':RegShiftedRoiDir,\
                                           'DicomDir':RegShiftedDicomDir,\
                                           'SeriesNo':RegShiftedSeriesNo,\
                                           'SeriesDesc':'',\
                                           'SliceNos':[]
                                           }
                             })
    
            DataDict.update({'RegMapped':{#'ScanLabel':RegMapPosCorrScanLabel,\
                                              'RoiDir':RegMappedRoiDir,\
                                              'DicomDir':RegMappedDicomDir,\
                                              'SeriesNo':RegMappedSeriesNo,\
                                              'SeriesDesc':'',\
                                              'SliceNos':[]
                                              }
                             })
            
                                     
        else:
            # Registered:
            #RegScanLabel = ToScanLabel + '_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegSeriesNo = MapSeriesNo + '0' + '0' + '0' + '1'
            elif RegMethod=='affine':
                RegSeriesNo = MapSeriesNo + '0' + '0' + '0' + '2'
            elif RegMethod=='non-rigid':
                RegSeriesNo = MapSeriesNo + '0' + '0' + '0' + '3'
                
            if Debug:
                print('\nRegSeriesNo =', RegSeriesNo)
            
            #RegDicomDir = os.path.join(RootDir, RegScanLabel + ' DICOMs')
            #RegRoiDir = os.path.join(RootDir, RegScanLabel + ' ROIs')
            RegDicomDir = os.path.join(RootDir, RegSeriesNo + ' DICOMs')
            RegRoiDir = os.path.join(RootDir, RegSeriesNo + ' ROIs')
            
            # Append directories:
            #AllDirs.append(RegRoiDir)
            #AllDirs.append(RegDicomDir)
        
            
            # Registered mapped:
            #RegMapScanLabel = ToScanLabel + '_mapped-to-' + FromScanLabel \
            #                  + '_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegMapSeriesNo = MapSeriesNo + '0' + '0' + '1' + '1'
            elif RegMethod=='affine':
                RegMapSeriesNo = MapSeriesNo + '0' + '0' + '1' + '2'
            elif RegMethod=='non-rigid':
                RegMapSeriesNo = MapSeriesNo + '0' + '0' + '1' + '3'
                                        
            #RegMapDicomDir = os.path.join(RootDir, RegMapScanLabel + ' DICOMs')
            #RegMapRoiDir = os.path.join(RootDir, RegMapScanLabel + ' ROIs')
            RegMapDicomDir = os.path.join(RootDir, RegMapSeriesNo + ' DICOMs')
            RegMapRoiDir = os.path.join(RootDir, RegMapSeriesNo + ' ROIs')
        
            # Append directories:
            #AllDirs.append(RegMapRoiDir)
            #AllDirs.append(RegMapDicomDir)
            
            
            # Add the above directories to dataDict:                     
            DataDict.update({'Reg':{#'ScanLabel':RegScanLabel,\
                                    'RoiDir':RegRoiDir,\
                                    'DicomDir':RegDicomDir,\
                                    'SeriesNo':RegSeriesNo,\
                                    'SeriesDesc':'',\
                                    'SliceNos':[]
                                    }
                             })
    
            DataDict.update({'RegMap':{#'ScanLab':RegMapScanLabel,\
                                       'RoiDir':RegMapRoiDir,\
                                       'DicomDir':RegMapDicomDir,\
                                       'SeriesNo':RegMapSeriesNo,\
                                       'SeriesDesc':'',\
                                       'SliceNos':[]
                                       }
                             })
            
        
    # CREATE THE DIRECTORIES IF THEY DON'T EXIST:
    for dirPath in AllDirs:
        if Debug:
            print('\nChecking if directory', dirPath, 'exists')
        try:
            os.mkdir(dirPath)
            if Debug:
                print('\nDirectory', dirPath, 'created')
        except FileExistsError:
            if Debug:
                print('\nDirectory', dirPath, 'already exists')
    
    
    
    return DataDict