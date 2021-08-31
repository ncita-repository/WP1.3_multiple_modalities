# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:50:51 2021

@author: ctorti
"""

from importlib import reload
import image_tools.attrs_info
reload(image_tools.attrs_info)

def get_seg_data_from_list_of_segs(
        listOfSegs, listOfDicomDirs, p2c=False):
    """ 
    Note:
    
    listOfSegs is a list (which could be of length 1) of SEG (Pydicom) objects.
    
    listOfDicomDirs is a list (which could be of length 1) of strings 
    containing the directory containing DICOMs that correspond to each SEG.

    """
    
    from dicom_tools.metadata import get_dcm_uids
    from dicom_tools.imports import get_dcm_fpaths
    from seg_tools.metadata import get_p2sIndsBySeg
    from seg_tools.pixarrs import get_pixarrBySeg
    from image_tools.attrs_info import get_im_attrs
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of get_seg_data_from_list_of_segs():')
        print('\n\n', '-'*120)
    
    listOfPixarrBySeg = []
    listOfF2SindsBySeg = []
    #listOfSegNums = []
    listOfImSizes = []
    listOfImSpacings = []
    listOfImSlcThicks = []
    listOfImIPPs = []
    listOfImDirections = []
    listOfDicomFpaths = []
    
    for i in range(len(listOfSegs)):
        if p2c:
            print(f'\n\nlistOfSegs[{i}]:')
            
        studyUID, seriesUID, FORUID, SOPUIDs = get_dcm_uids(listOfDicomDirs[i])
        
        if listOfSegs[i]:
            f2sIndsBySeg = get_p2sIndsBySeg(listOfSegs[i], SOPUIDs)
            
            pixarrBySeg = get_pixarrBySeg(listOfSegs[i], f2sIndsBySeg, p2c)
            
            if p2c:      
                print(f'f2sIndsBySeg = {f2sIndsBySeg}')
            
        else:
            pixarrBySeg = None
            f2sIndsBySeg = None
            
            if p2c:      
                print('No segmentations for this dataset')
        
        listOfPixarrBySeg.append(pixarrBySeg)
        listOfF2SindsBySeg.append(f2sIndsBySeg)
            
        if p2c:      
            print(f'\nlen(listOfF2SindsBySeg) = {len(listOfF2SindsBySeg)}')
            #print(f'\nPixArr.shape = {PixArr.shape}')
            #[print(f'np.amax(PixArr[{i}]) = {np.amax(PixArr[i])}') for i in range(PixArr.shape[0])]
        
        imSize, imSpacings, imSlcThick, imIPPs, imDirections,\
            warnings = get_im_attrs(listOfDicomDirs[i], p2c=p2c)
        
        listOfImSizes.append(imSize)
        listOfImSpacings.append(imSpacings)
        listOfImSlcThicks.append(imSlcThick)
        listOfImIPPs.append(imIPPs)
        listOfImDirections.append(imDirections)
        listOfDicomFpaths.append(get_dcm_fpaths(listOfDicomDirs[i]))
        
    #ListOfDcmFpaths, ListOfIPPs, ListOfDirections, ListOfSpacings\
    #    = GetImAttributesFromListOfDcmDirs(listOfDicomDirs)
    
    if p2c:
        print('-'*120)
        
    return listOfPixarrBySeg, listOfF2SindsBySeg, listOfImSizes,\
        listOfImSpacings, listOfImIPPs, listOfImDirections, listOfDicomFpaths

