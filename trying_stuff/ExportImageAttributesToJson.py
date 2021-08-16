# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:18:32 2020

@author: ctorti
"""




def ExportImageAttributesToJson(DicomDir, FilePrefix):
    """
    Export the Origin, Direction cosines and pixel spacings of a collection of
    DICOM images to JSON files in the current working directory.
    
    
    Input:
        DicomDir   - String of the path to the directory containing DICOMs  
        
        FilePrefix - String of the prefix to be added to the JSON file names
                     (e.g. The prefix 'Fix' would result in the file name
                     'FixOrigin.json')
        
        
    Returns:
        Nothing returned;
        
        The following JSON files will be exported to the current working 
        directory:
            
        '_Origin.json'
            List of the ImagePositionPatient, e.g. [x0, y0, z0]
        
        '_Directions.json'
            List of the direction cosine along x (rows), y (columns) and z 
            (slices), e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
        
        '_Spacings.json'
            List of the pixel spacings along x, y and z, e.g. [di, dj, dk]
                    
        '_Dimensions.json'
            List of the dimensions of the image stack
            
        where the _ represents the FilePrefix
        
    """

    import json
    from GetImageAttributes import GetImageAttributes
    
    Origin, Directions,\
    Spacings, Dimensions = GetImageAttributes(DicomDir=DicomDir, 
                                              Package='pydicom')
    
    
    
    # Group all variables into a list:
    Vars = [Origin, Directions, Spacings, Dimensions]
    
    Names = ['Origin', 'Dirs', 'Spacings', 'Dims']
        
    
    for i in range(len(Vars)):
        # Complete the file name:
        FileName = FilePrefix + Names[i] + '.json'
    
        with open(FileName, 'w') as file:
            json.dump(Vars[i], file)
    
    return


    