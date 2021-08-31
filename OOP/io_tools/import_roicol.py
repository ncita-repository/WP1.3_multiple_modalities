# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:10:28 2021

@author: ctorti
"""

from pydicom import dcmread


class RoicollectionImporter:
    """
    05/08/21:
        Rather than use a separate class for ROI collections I've decided to 
        import them within the Dataset class.
    
    THIS IS OBSOLETE because the ROI Collection is imported along with other 
    data using the DataImporter class (io_tools.import_data.py).
        
    This class creates an roicolData Object from a SEG downloaded from XNAT.
    
    ROI Collection data refers to the ROI Collection itself (RTS/SEG) and
    associated data.
    
    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
    
    Returns
    -------
    
    """
    
    def __init__(self, params):
        
        # Import RTS/SEG Pydicom Objects for Source and Target:
        self.import_roicol(params, 'src')
        self.import_roicol(params, 'trg')
        
    def import_roicol(self, params, srcORtrg):
        """
        Import an ROI Collection as a Pydicom Object.
        
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
                self.srcRoicol = dcmread(filepath)
            else:
                self.srcRoicol = None
        else:
            # i.e.: srcORtrg == 'trg'
            filepath = params.trgXnatParams['roicolFpath']
            
            if filepath != None:
                self.trgRoicol = dcmread(filepath)
            else:
                self.trgRoicol = None
        