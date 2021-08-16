# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:52:56 2020

@author: ctorti
"""




def ImportImageAttributesFromJson(FilePrefix):
    """
    Import the Origin, Direction cosines and pixel spacings of a collection of
    DICOM images from JSON files in the current working directory.
    
    
    Input:    
        FilePrefix - String of the prefix that was added to the JSON file names
                     (e.g. 'Fix' or 'Mov' before '_Origin.json')
        
        
    Returns:
        Attributes - List of the following attributes:
            
            Origin     - List of the ImagePositionPatient, e.g. [x0, y0, z0]
            
            Directions - List of the direction cosine along x (rows), y   
                         (columns) and z (slices), 
                         e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
            
            Spacings   - List of the pixel spacings along x, y and z, 
                         e.g. [di, dj, dk]
                        
            Dimensions - List of the dimensions of the image stack
            
        
    """

    import json
    
    
    Names = ['Origin', 'Dirs', 'Spacings', 'Dims']
        
    Attributes = []
    
    for i in range(len(Names)):
        # Complete the file name:
        FileName = FilePrefix + Names[i] + '.json'
    
        with open(FileName, 'r') as file:
            Attributes.append(json.load(file))
    
    return Attributes
