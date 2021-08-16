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
    ScanPosCorr = DataDict['ScanPosCorr']
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
    For the DICOMs that result from applying a "scan-position-correction" to 
    the Image Positions of the DICOMs in ToDicoms and ToMapDicoms.
    """
    
    if ScanPosCorr:
        ##toPosCorrScanNo = toDicomScanNo + 100
        ##toPosCorrDicomScanNo = toDicomScanNo*10
        ##toPosCorrDicomScanLab = str(toPosCorrDicomScanNo)
        #toPosCorrDicomScanLab = toDicomScanLab + a + '1' + '0' + '0'
        ##toPosCorrDicomScanNo = int(toPosCorrDicomScanLab)
        #toMapPosCorrDicomScanLab = toDicomScanLab + a + '1' + '0' + '0'
        
        #toPosCorrDirName = toDirName + ' ' + searchBy + '-corr'
        
        
        # Position-corrected:
        #PosCorrScanLabel = ToScanLabel + '_posCorr'
        PosCorrSeriesNo = ToSeriesNo + '0' + '1'
        
        if Debug:
            print('\nPosCorrSeriesNo =', PosCorrSeriesNo)
        
        """
        Note: 
            Format of Series Numbers will be the ToSeriesNo + '0' + a number of
            binaries which indicate whether position-correction was done,
            mapping, registration, etc.
        """
      
        #PosCorrDicomDir = os.path.join(RootDir, PosCorrScanLabel + ' DICOMs')
        #PosCorrRoiDir = os.path.join(RootDir, PosCorrScanLabel + ' ROIs')
        PosCorrDicomDir = os.path.join(RootDir, PosCorrSeriesNo + ' DICOMs')
        PosCorrRoiDir = os.path.join(RootDir, PosCorrSeriesNo + ' ROIs')
        
        # Append directories:
        #AllDirs.append(PosCorrRoiDir)
        #AllDirs.append(PosCorrDicomDir)
        
        
        # Mapped position-corrected:
        #MapPosCorrScanLabel = ToScanLabel + '_posCorr_mapped-to-' \
        #                      + FromScanLabel
        MapPosCorrSeriesNo = PosCorrSeriesNo + '1'
        
        if Debug:
            print('\nMapPosCorrSeriesNo =', MapPosCorrSeriesNo)
         
        #MapPosCorrDicomDir = os.path.join(RootDir, MapPosCorrScanLabel + ' DICOMs')
        #MapPosCorrRoiDir = os.path.join(RootDir, MapPosCorrScanLabel + ' ROIs')
        MapPosCorrDicomDir = os.path.join(RootDir, MapPosCorrSeriesNo + ' DICOMs')
        MapPosCorrRoiDir = os.path.join(RootDir, MapPosCorrSeriesNo + ' ROIs')
        
        # Append directories:
        #AllDirs.append(MapPosCorrRoiDir)
        #AllDirs.append(MapPosCorrDicomDir)
        
        
        # Add the above directories to dataDict:                     
        DataDict.update({'PosCorr':{#'ScanLabel':PosCorrScanLabel,\
                                    'RoiDir':PosCorrRoiDir,\
                                    'DicomDir':PosCorrDicomDir,\
                                    'SeriesNo':PosCorrSeriesNo,\
                                    'SeriesDesc':'',\
                                    'SliceNos':[]
                                     }
                         })
    

        DataDict.update({'MapPosCorr':{#'ScanLabel':MapPosCorrScanLabel,\
                                       'RoiDir':MapPosCorrRoiDir,\
                                       'DicomDir':MapPosCorrDicomDir,\
                                       'SeriesNo':MapPosCorrSeriesNo,\
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
        data for the case of scanPosCorr or not will be duplicated with the
        addition of "Map" to allow the the possibility that DICOMs will need to
        be mapped.
        """
        
        #toMapDicomScanLab = toDicomScanLabel + a + '0' + '1' + '0'
        #MapScanLabel = ToScanLabel + '_mapped-to-' + FromScanLabel
        MapSeriesNo = ToSeriesNo + '0' + '0' + '1'
        
        if Debug:
            print('\nMapSeriesNo =', MapSeriesNo)
        
        #MapDicomDir = os.path.join(RootDir, MapScanLabel + ' DICOMs')
        #MapRoiDir = os.path.join(RootDir, MapScanLabel + ' ROIs')
        MapDicomDir = os.path.join(RootDir, MapSeriesNo + ' DICOMs')
        MapRoiDir = os.path.join(RootDir, MapSeriesNo + ' ROIs')
        
        DataDict.update({'Mapped':{#'ScanLabel':MapScanLabel,\
                                   'RoiDir':MapRoiDir,\
                                   'DicomDir':MapDicomDir,\
                                   'SeriesNo':MapSeriesNo,\
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
        if ScanPosCorr:
            # Registered position-corrected:
            #RegPosCorrScanLabel = ToScanLabel + '_posCorr_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegPosCorrSeriesNo = PosCorrSeriesNo + '1'
            elif RegMethod=='affine':
                RegPosCorrSeriesNo = PosCorrSeriesNo + '2'
            elif RegMethod=='non-rigid':
                RegPosCorrSeriesNo = PosCorrSeriesNo + '3'
                
                                                              
            #RegPosCorrDicomDir = os.path.join(RootDir, RegPosCorrScanLabel + ' DICOMs')
            #RegPosCorrRoiDir = os.path.join(RootDir, RegPosCorrScanLabel + ' ROIs')
            RegPosCorrDicomDir = os.path.join(RootDir, RegPosCorrSeriesNo + ' DICOMs')
            RegPosCorrRoiDir = os.path.join(RootDir, RegPosCorrSeriesNo + ' ROIs')
            
            # Append directories:
            #AllDirs.append(RegPosCorrRoiDir)
            #AllDirs.append(RegPosCorrDicomDir)
            
                            
            # Registered mapped position-corrected:
            #RegMapPosCorrScanLabel = ToScanLabel + '_posCorr_mapped-to-' \
            #                         + FromScanLabel + '_' + RegMethod + 'Reg'
            if RegMethod=='rigid':
                RegMapPosCorrSeriesNo = MapPosCorrSeriesNo + '1'
            elif RegMethod=='affine':
                RegMapPosCorrSeriesNo = MapPosCorrSeriesNo + '2'
            elif RegMethod=='non-rigid':
                RegMapPosCorrSeriesNo = MapPosCorrSeriesNo + '3'
            
            if Debug:
                print('\nRegMapPosCorrSeriesNo =', RegMapPosCorrSeriesNo)
            
            #RegMapPosCorrDicomDir = os.path.join(RootDir, RegMapPosCorrScanLabel + ' DICOMs')
            #RegMapPosCorrRoiDir = os.path.join(RootDir, RegMapPosCorrScanLabel + ' ROIs')
            RegMapPosCorrDicomDir = os.path.join(RootDir, RegMapPosCorrSeriesNo + ' DICOMs')
            RegMapPosCorrRoiDir = os.path.join(RootDir, RegMapPosCorrSeriesNo + ' ROIs')
            
            # Append directories:
            #AllDirs.append(RegMapPosCorrRoiDir)
            #AllDirs.append(RegMapPosCorrDicomDir)
            
            
            # Add the above directories to dataDict:                     
            DataDict.update({'RegPosCorr':{#'ScanLabel':RegPosCorrScanLabel,\
                                           'RoiDir':RegPosCorrRoiDir,\
                                           'DicomDir':RegPosCorrDicomDir,\
                                           'SeriesNo':RegPosCorrSeriesNo,\
                                           'SeriesDesc':'',\
                                           'SliceNos':[]
                                           }
                             })
    
            DataDict.update({'RegMapPosCorr':{#'ScanLabel':RegMapPosCorrScanLabel,\
                                              'RoiDir':RegMapPosCorrRoiDir,\
                                              'DicomDir':RegMapPosCorrDicomDir,\
                                              'SeriesNo':RegMapPosCorrSeriesNo,\
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