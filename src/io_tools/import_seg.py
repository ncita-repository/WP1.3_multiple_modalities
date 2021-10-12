# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:10:28 2021

@author: ctorti
"""

from pydicom import dcmread


class SegImporter:
    """
    This class creates a segData Object from a SEG downloaded from XNAT.
    
    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
    
    Returns
    -------
    
    """
    
    def __init__(self, params):
        
        # Import SEG Pydicom Objects for Source and Target:
        self.import_seg(params, 'src')
        self.import_seg(params, 'trg')
        
    def import_seg(self, params, srcORtrg):
        """
        Import SEG as a Pydicom Object.
        
        Parameters
        ----------
        params : Params Object
            Object containing various parameters.
        srcORtrg : str
            'src' or 'trg' for Source or Target.
        
        Returns
        -------
        None.
        
        """
        
        if srcORtrg == 'src':
            filepath = params.srcXnatParams['roicolFpath']
            
            if filepath != None:
                self.srcSeg = dcmread(filepath)
            else:
                self.srcSeg = None
        else:
            filepath = params.trgXnatParams['roicolFpath']
            
            if filepath != None:
                self.trgSeg = dcmread(filepath)
            else:
                self.trgSeg = None
        