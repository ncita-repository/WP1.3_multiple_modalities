# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:19:21 2021

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
TESTING FUNCTIONS
******************************************************************************
******************************************************************************
"""  


def CreateCopyRoiConfigFile():
    """ 
    Create a configuration file ('CopyRoiDefaults.py') assigning default values 
    for some inputs required for CopyRoiCol() (in RoiCopyTools.py).
    
    The file will be exported to the current working directory.
    """
    
    import os
    from pypref import Preferences
    
    ConfigFname = 'CopyRoiDefaults.py'
    
    cwd = os.getcwd()
    
    Dict = {'TrgRoiColMod' : None, 
            'TrgRoiColName' : None, 
            'TrgSliceNum' : None, 
            'TxtToAddToTrgRoiColName' : '', 
            'XnatSession' : None, 
            'PathsDict' : None, 
            'ExportRootDir' : 'cwd', 
            'LogToConsole' : False, 
            'ExportLogFiles' : False,
            'ResInterp' : 'BlurThenLinear', 
            'PreResVariance' : (1,1,1),
            'ApplyPostResBlur' : False, 
            'PostResVariance' : (1,1,1),
            'ForceReg' : False,
            'TxMatrix' : None,
            'SelxOrSitk' : 'Selx',
            'Transform' : 'affine',
            'MaxIters' : '512', 
            'TxInterp' : 'NearestNeighbor', 
            'ApplyPostTxBin' : True,
            'ApplyPostTxBlur' : True, 
            'PostTxVariance' : (1,1,1),
            'ExportNewTrgRoiCol' : True,
            'RtsExportDir' : os.path.join(cwd, 'new_RTS'),
            'SegExportDir' : os.path.join(cwd, 'new_SEG'),
            'ExportNewDro' : True,
            'DroExportDir' : os.path.join(cwd, 'new_DRO'),
            'UploadNewDro' : True,
            'PlotAllSlices' : False,
            'RtsPlotExportDir' : os.path.join(cwd, 'plots_RTS'),
            'SegPlotExportDir' : os.path.join(cwd, 'plots_SEG'),
            'LogExportDir' : os.path.join(cwd, 'logs')
            }
    
    pref = Preferences(directory=cwd, filename=ConfigFname)

    pref.set_preferences(Dict)
    
    return




def CreateTestConfigFile():
    """ 
    Create a configuration file ('CopyRoiTestConfig.py') assigning values for  
    inputs required for CopyRoiCol() (in RoiCopyTools.py) based on the tests in 
    https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g.
    
    The file will be exported to the current working directory.
    """
    
    import os
    from pypref import Preferences
    
    ConfigFname = 'CopyRoiTestConfig.py'
    
    Dict = {'XnatUrl' : 'http://10.1.1.20', 
            'ProjId' : 'NCITA_TEST', 
            'SubjLabel' : 'TCGA-BB-A5HY',
            
            'RR1': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '4', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 4'
                    },
            
            'RR2': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '10', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 10'
                    },
            
            'RR3': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 5'
                    },
            
            'RR4': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '11', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 11'
                    },
            
            'RR5': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '14', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 14'
                    },
                    
            'RR6': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '10', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 10'
                    },
            
            'RR7': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '12', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 12'
                    },
            
            'RR8': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '13', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 13'
                    },
            
            'RR9': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '14', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 14'
                    },
            
            'RR10': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '11', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'RR11': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '13', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'RR12': {'SrcExpLabel' : 'Session5', 
                     'SrcScanId' : '12', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'RR13': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '11', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            'RR14': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '13', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            'RR15': {'SrcExpLabel' : 'Session5', 
                     'SrcScanId' : '12', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            #'RR10r': {'SrcExpLabel' : 'Session3',
            #          'SrcScanId' : '6', 
            #          'SrcRoiColMod' : 'RTSTRUCT', 
            #          'SrcRoiColName' : 'Ventricles', 
            #          'SrcRoiName' : None,
            #          'SrcSliceNum' : None, 
            #          'TrgExpLabel' : 'Session4', 
            #          'TrgScanId' : '11', 
            #          'TrgRoiColMod' : 'RTSTRUCT', 
            #          'TrgRoiColName' : None, 
            #          'TrgRoiName' : None,
            #          'TrgSliceNum' : None,
            #          'TxtToAddToTrgRoiColName' : ' 4_11'
            #         },
            
            'RD1': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '2', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgRoiColName' : ' 2'
                    },
            
            'RD2': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '4', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgRoiColName' : ' 4'
                    },
            
            'RD3': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 65, # Source slice 5 maps to Target slices 55-64
                    'TxtToAddToTrgRoiColName' : ' 5'
                    },
            
            'RD4': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session6',
                    'TrgScanId' : '3', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 8, # Source slice 5 maps to Target slice 7
                    'TxtToAddToTrgRoiColName' : ' 6_3'
                    },
            
            'RD5': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session6',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 78, # Source slice 5 maps to Target slices 67-77
                    'TxtToAddToTrgRoiColName' : ' 6_5'
                    },
            
            'RD6': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'RTSTRUCT', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 18, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '10', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 21, 
                    'TxtToAddToTrgRoiColName' : ' 10'
                    },
                    
                    
            'SR1': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '4', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 4'
                    },
            
            'SR2': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '10', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 10'
                    },
            
            'SR3': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 5'
                    },
            
            'SR4': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '11', 
                    'TrgRoiColMod' : 'RTSTRUCT', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 11'
                    },
            
            'SR5': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '13', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '14', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 14'
                    },
                    
            'SR6': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '10', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 10'
                    },
            
            'SR7': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '12', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 12'
                    },
            
            'SR8': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '13', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 13'
                    },
            
            'SR9': {'SrcExpLabel' : 'Session4', 
                    'SrcScanId' : '11', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgExpLabel' : 'Session4',
                    'TrgScanId' : '14', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgRoiColName' : ' 14'
                    },
            
            'SR10': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '11', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'SR11': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '13', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'SR12': {'SrcExpLabel' : 'Session5', 
                     'SrcScanId' : '12', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session3',
                     'TrgScanId' : '6', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 3_6'
                     },
            
            'SR13': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '11', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            'SR14': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '13', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            'SR15': {'SrcExpLabel' : 'Session5', 
                     'SrcScanId' : '12', 
                     'SrcRoiColMod' : 'SEG', 
                     'SrcRoiColName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session1',
                     'TrgScanId' : '5', 
                     'TrgRoiColMod' : 'SEG', 
                     'TrgRoiColName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' 1_5'
                     },
            
            'SD1': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '2', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgRoiColName' : ' 2'
                    },
            
            'SD2': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '4', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgRoiColName' : ' 4'
                    },
            
            'SD3': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session2',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 65, # Source slice 5 maps to Target slices 55-64
                    'TxtToAddToTrgRoiColName' : ' 5'
                    },
            
            'SD4': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session6',
                    'TrgScanId' : '3', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 8, # Source slice 5 maps to Target slice 7
                    'TxtToAddToTrgRoiColName' : ' 6_3'
                    },
            
            'SD5': {'SrcExpLabel' : 'Session2', 
                    'SrcScanId' : '2', 
                    'SrcRoiColMod' : 'SEG', 
                    'SrcRoiColName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgExpLabel' : 'Session6',
                    'TrgScanId' : '5', 
                    'TrgRoiColMod' : 'SEG', 
                    'TrgRoiColName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 78, # Source slice 5 maps to Target slices 67-77
                    'TxtToAddToTrgRoiColName' : ' 6_5'
                    },
            
            'RR16': {'SrcExpLabel' : 'Session4', 
                     'SrcScanId' : '13', 
                     'SrcRoiColMod' : 'RTSTRUCT', 
                     'SrcRoiColName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgExpLabel' : 'Session4',
                     'TrgScanId' : '11', 
                     'TrgRoiColMod' : 'RTSTRUCT', 
                     'TrgRoiColName' : 'Ventricles', 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgRoiColName' : ' Ser11 2021-03-26'
                     }
            }
    
    pref = Preferences(directory=os.getcwd(), filename=ConfigFname)

    pref.set_preferences(Dict)
    
    return











""" OLD FUNCTIONS BELOW """





def CreateDictOfPaths():
    """
    DEFINE FILE PATHS TO SOURCE AND TARGET (IF APPLICABLE) RTS/SEG FILES, 
    TO SOURCE AND TARGET DICOM DIRECTORIES, AND OTHER INPUTS REQUIRED.
    
    WARNING:
        Don't forget to add new directories/filepaths to the dictionary Dict
        below!
    """

    import os
    
    Dict = {}
    
    
    """ Subject TCGA-BB-A5HY in XNAT_TEST """
    
    Sub1_Dir = r"C:\Data\XNAT_TEST\TCGA-BB-A5HY"
    
    Dict['Sub1_Dir'] = Sub1_Dir
    
    """ Session 1 """
    
    Dict['Ses1_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\5_AX T2 BRAIN  PROPELLER\DICOM")
    
    Dict['Ses1_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\6_AX GRE BRAIN\DICOM")
    
    Dict['Ses1_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\7_AX T1 BRAIN\DICOM")
    
    Dict['Ses1_Ser8_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\8_Ax T2 FLAIR BRAIN POST\DICOM")
    
    Dict['Ses1_Ser9_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session1\scans\9_AX T1 POST BRAIN\DICOM")
    
    
    """ Session 2 """
    
    Dict['Ses2_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\2_Head Routine  5.0  H30s\DICOM")
    
    Dict['Ses2_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\3_Head Routine  5.0  J30s  3\DICOM")
    
    Dict['Ses2_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\4_Head Routine  5.0  H60s\DICOM")
    
    Dict['Ses2_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session2\scans\5_Head Routine  0.75  H20s\DICOM")

    Dict['Ses2_Ser2_LeftEyeCT_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hkhGnxQ_2kX4w5b39s0\RTSTRUCT",
                   r"AIM_20210201_000731.dcm")
    
    Dict['Ses2_Ser2_LeftEyeCT_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hlDxpyQ_7cJDZotLuP_icr_roiCollectionData\SEG",
                   r"SEG_20210213_140215.dcm") 
    
    Dict['Ses2_Ser2_RightEyeCT_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hkhGnxQ_2NZvk8B9pwA\RTSTRUCT",
                   r"AIM_20210201_000842.dcm") 
    
    Dict['Ses2_Ser2_RightEyeCT_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session2\assessors",
                   r"RoiCollection_hlDxpyQ_AqEfVVmmysA_icr_roiCollectionData\SEG",
                   r"SEG_20210213_135857.dcm") 
    
    
    """ Session 3 """
    
    Dict['Ses3_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"6_T2 FLAIR Ax\DICOM")
    
    Dict['Ses3_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"7_Diff Ax\DICOM")
    
    Dict['Ses3_Ser13_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"13_T1 Ax Post fs IAC\DICOM")
    
    Dict['Ses3_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session3\scans", r"14_T1 Ax Post Brain\DICOM")
    
    
    """ Session 4 """
    
    Dict['Ses4_Ser10_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans", r"10_T1  AXIAL 2MM POST repeat_S7_DIS3D\DICOM")
    
    Dict['Ses4_Ser11_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\11_T1 3D AX POST_S5_DIS3D\DICOM")
    
    Dict['Ses4_Ser12_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\12_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    
    Dict['Ses4_Ser13_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\13_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    
    Dict['Ses4_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session4\scans\14_T2_FLAIR 3MM_S2_DIS3D\DICOM")
    
    Dict['Ses4_Ser11_Ventricles_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                    r"RoiCollection_hkhGnxQ_BNths5fNmA\RTSTRUCT",
                     r"AIM_20210131_230845.dcm") 
    
    Dict['Ses4_Ser11_Ventricles_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                    r"RoiCollection_hlDxpyQ_5Dz9OulRH6n_icr_roiCollectionData\SEG",
                     r"SEG_20210213_144855.dcm") 
    
    Dict['Ses4_Ser13_NasalCavity_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                   r"RoiCollection_hkhGnxQ_6daB7S8Eckp\RTSTRUCT",
                   r"AIM_20210131_214821.dcm") 
    
    Dict['Ses4_Ser13_NasalCavity_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session4\assessors",
                   r"RoiCollection_hlDxpyQ_AJ47Vd6ikeE_icr_roiCollectionData\SEG",
                   r"SEG_20210213_152014.dcm") 
    
    
    """ Session 5"""
    
    Dict['Ses5_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\6_T2 FLAIR Ax\DICOM")
    
    Dict['Ses5_Ser7_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\7_CISS Ax\DICOM")
    
    Dict['Ses5_Ser12_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\12_T2 Ax fs Post Brain\DICOM")
    
    Dict['Ses5_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\14_T1 Ax Post fs IAC\DICOM")
    
    Dict['Ses5_Ser16_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session5\scans\16_T1 Ax Post fs IAC repeat\DICOM")
    
    Dict['Ses5_Ser12_LeftEye_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hkhGnxQ_4rFYFKE1r6r\RTSTRUCT",
                   r"AIM_20210131_212509.dcm")
    
    Dict['Ses5_Ser12_LeftEye_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hlDxpyQ_8CVLpZmvRD8_icr_roiCollectionData\SEG",
                   r"SEG_20210213_153745.dcm")
    
    Dict['Ses5_Ser16_Tumour_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hkhGnxQ_79YdSO5eYPo\RTSTRUCT",
                   r"AIM_20210131_211622.dcm")
    
    Dict['Ses5_Ser16_Tumour_SegFpath']\
    = os.path.join(Sub1_Dir, r"Session5\assessors",
                   r"RoiCollection_hlDxpyQ_18XInHEwVjJ_icr_roiCollectionData\SEG",
                   r"SEG_20210213_155524.dcm")
    
    
    """ Session 6"""
    
    Dict['Ses6_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\3_Head Routine  5.0  H30f\DICOM")
    
    Dict['Ses6_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Dict['Ses6_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session6\scans\5_Head Routine  0.75  H20f\DICOM")
    
    
    """ Session 7"""
    
    Dict['Ses7_Ser21_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session7\scans\21_T1  AXIAL 2MM POST_S6_DIS3D\DICOM")
    
    Dict['Ses7_Ser22_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\22_T1 3D AX POST_S5_DIS3D\DICOM")
    
    Dict['Ses7_Ser23_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\23_T2 AX 3MM STRAIGHT (post)_S4_DIS3D\DICOM")
    
    Dict['Ses7_Ser24_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\24_T1  AXIAL 2MM_S3_DIS3D\DICOM")
    
    Dict['Ses7_Ser25_DcmDir']\
     = os.path.join(Sub1_Dir, r"Session7\scans\25_T2_FLAIR 3MM_S2_DIS3D\DICOM")

    
    """ Session 8"""
    
    Dict['Ses8_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\2_Head Routine  5.0  H30f\DICOM")
    
    Dict['Ses8_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\4_Head Routine  5.0  J30f  3\DICOM")
    
    Dict['Ses8_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session8\scans\5_Head Routine  0.75  H20f\DICOM")
    
    
    """ Session 9"""
    
    Dict['Ses9_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\3_T2 FLAIR Ax\DICOM")
    
    Dict['Ses9_Ser4_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\4_Diff Ax\DICOM")
    
    Dict['Ses9_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\5_Diff Ax_ADC\DICOM")
    
    Dict['Ses9_Ser14_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\14_T2spc Ax\DICOM")
    
    Dict['Ses9_Ser15_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session9\scans\15_T1 3D Ax Post\DICOM")
    
    
    """ Session 10"""
    
    Dict['Ses10_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\3\DICOM")
    
    Dict['Ses10_Ser6_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\6\DICOM")
    
    Dict['Ses10_Ser19_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\19\DICOM")
    
    Dict['Ses10_Ser22_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session10\22\DICOM")
    
    Dict['Ses10_Ser3_RightEye_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session10", r"AIM_20210131_183753\RTSTRUCT", 
                   r"AIM_20210131_183753.dcm")
    
    Dict['Ses10_Ser3_LowerBrain_RtsFpath']\
    = os.path.join(Sub1_Dir, r"Session10", r"AIM_20210131_232127\RTSTRUCT",
                   r"AIM_20210131_232127.dcm")
    
    
    """ Session 11"""
    
    Dict['Ses11_Ser2_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\2_Head Routine  0.75  H20s\DICOM")
    
    Dict['Ses11_Ser3_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\3_Head Routine  5.0  H30s\DICOM")
    
    Dict['Ses11_Ser5_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\5_Head Routine  5.0  J30s  3\DICOM")
    
    Dict['Ses11_Ser9_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\9_TEMPORAL BONES RIGHT  5.0  H30s\DICOM")
    
    Dict['Ses11_Ser10_DcmDir']\
    = os.path.join(Sub1_Dir, r"Session11\scans\10_TEMPORAL BONES RIGHT  0.6  H30s\DICOM") 
    
    
    
    """ Subject ACRIN-FMISO-Brain-011 in ACRIN-FMISO-Brain """
    
    Sub2_Dir = r"C:\Data\Cancer imaging archive\ACRIN-FMISO-Brain\ACRIN-FMISO-Brain-011"
    Sub2_RoiDir = r"C:\Data\Data copied to Google Drive\ACRIN-FMISO-Brain-011"
    Sub2_RtsDir = os.path.join(Sub2_RoiDir, "Fixed_RTS_Collections")
    Sub2_SegDir = os.path.join(Sub2_RoiDir, "Fixed_SEG_Collections")
    
    Dict['Sub2_Dir'] = Sub2_Dir
    Dict['Sub2_RoiDir'] = Sub2_RoiDir
    Dict['Sub2_RtsDir'] = Sub2_RtsDir
    Dict['Sub2_SegDir'] = Sub2_SegDir

    MR4_Dir = r"04-10-1960-MRI Brain wwo Contrast-69626"
    # IPP =
    # [1.00, 0.00, 0.07,
    # -0.02, 0.95, 0.30,
    # -0.07, -0.30, 0.95]
    
    MR12_Dir = r"06-11-1961-MRI Brain wwo Contrast-79433"
    # IPP =
    # [0.99, 0.11, 0.07,
    # -0.10, 0.99, -0.03,
    # -0.07, 0.02, 1.00]
    
    Dict['MR4_Dir'] = MR4_Dir
    
    Dict['MR12_Dir'] = MR12_Dir
    
    
    """ DICOM directories for MR4 """
    
    Dict['MR4_Ser3_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"3-T2 TSE AXIAL-07507") 
    # 378, 448, 21
    # 0.51, 0.51, 7.5 mm
    # 5.0 mm
    
    Dict['MR4_Ser4_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"4-T2 AXIAL FLAIR DARK FL-48830") 
    # 416, 512, 21 pix
    # 0.45, 0.45, 5.0 mm
    # 5.0 mm
    
    Dict['MR4_Ser5_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"5-T1 SE AXIAL 3MM-81246")  
    # 208, 256, 50 
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    Dict['MR4_Ser7_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"7-ep2ddiff3scantracep2ADC-22197")
    # 192, 192, 21 pix
    # 1.3, 1.3, 7.5 mm
    # 5.0 mm
    
    Dict['MR4_Ser8_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"8-T1 SE AXIAL POST FS FC-59362") 
    # 212, 256, 30 pix
    # 0.9, 0.9, 5.0 mm
    # 5.0 mm
    
    Dict['MR4_Ser9_DcmDir']\
    = os.path.join(Sub2_Dir, MR4_Dir, r"9-T1 SE AXIAL POST 3MM 5 MIN DELAY-07268") 
    # 208, 256, 50 pix
    # 0.9, 0.9, 3.0 mm
    # 3.0 mm
    
    
    """ RTS/SEG filepaths for MR4 """
    
    Dict['MR4_Ser3_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm")
    #= os.path.join(Sub2_RtsDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201112_162427.dcm") # tumour &
    # brain are not separate ROIs!
    
    Dict['MR4_Ser3_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") 
    #= os.path.join(Sub2_SegDir, "MR4_S3_brain_SEG_20201030_133110.dcm")
    #= os.path.join(Sub2_SegDir,"MR4_S3_tumour_and_brain_SEG_20201105_090359.dcm") # one ROI!
    #= os.path.join(Sub2_SegDir, "MR4_S3_tumour_and_brain_RTS_AIM_20201125_133250.dcm") # two ROIs
    #= os.path.join(Sub2_SegDir, "MR4_S3_tumour_SEG_20201216_085421.dcm") # single frame on slice 10
    
    Dict['MR4_Ser5_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    
    Dict['MR4_Ser5_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S5_tumour_SEG_20201216_084441.dcm") 
    #= os.path.join(Sub2_SegDir, "Seg_from_MR4_S5_tumour_RTS_AIM_20201021_151904.dcm") 
    
    Dict['MR4_Ser8_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    
    Dict['MR4_Ser8_SegFpath']\
    = os.path.join(Sub2_SegDir, "Seg_from_MR4_S8_tumour_RTS_AIM_20200511_073405.dcm") 
    
    Dict['MR4_Ser9_RtsFpath']\
    = os.path.join(Sub2_RtsDir, "MR4_S9_tumour_and_brain_RTS_AIM_20201116_122858.dcm")
    #= os.path.join(Sub2_RtsDir, "MR4_S9_brain_RTS_AIM_20201112_142559.dcm")
    
    Dict['MR4_Ser9_SegFpath']\
    = os.path.join(Sub2_SegDir, "MR4_S9_tumour_and_brain_SEG_20201105_115939.dcm")
    #= os.path.join(Sub2_SegDir, "Seg_from_MR4_S9_tumour_RTS_AIM_20201014_163721.dcm") 
    #= os.path.join(Sub2_SegDir, "MR4_S9_tumour_SEG_20201105_115408.dcm")
    
    
    
    """ DICOM directories for MR12 """
    
    Dict['MR12_Ser3_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"3-T2 TSE AX IPAT2-08659") 
    
    Dict['MR12_Ser4_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"4-T2 AXIAL FLAIR DARK FL-94212")
    
    Dict['MR12_Ser5_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"5-T1 SE AXIAL-43742")
    
    Dict['MR12_Ser8_DcmDir']\
    = os.path.join(Sub2_Dir, MR12_Dir, r"8-T1 SE AXIAL POST FS FC-81428") 
    
    return Dict





def GetPathInputs(TestNum, UseOrigTrgRoi=False):
    """
    Get paths from a dictionary.
    
    Inputs:
    ******                      
        
    TestNum : string
        The character string that denotes the test number to be run.
        
    UseOrigTrgRoi : boolean (optional; False by default)
        If True, if the value for the RTS/SEG filepath key is not None, and if
        a Direct copy is to be made, any existing contours/segmentations in the
        original Target RTS/SEG that match the desired RoiLabelToCopy input to
        CopyRoi() will be preserved.  If False, the value of the filepath key 
        to the original Target RTS/SEG will be over-written with None, so that 
        a new Target RTS/SEG file is created that does not contain any of the 
        original contours/segmentations even if they exist for the ROI of 
        interest.
        

    Outputs:
    *******
    
    SrcRoiFpath : string
        The full filepath to the Source ROI Collection.
    
    SrcSliceNum : integer (0-indexed) or None
        The slice number containing the contour/segmentation of interest, or
        None if all contours/segmentations and/or all ROIs/segments are of
        interest.
    
    SrcRoiLabel : string
        The ROI Collection name.
    
    TrgSliceNum
    
    SrcDcmDir : string
        The Source DICOM directory.
    
    TrgDcmDir : string
        The Target DICOM directory.
    
    TrgRoiFpath : string or None
        The full filepath to the Source ROI Collection, or None.
    
    TxtToAddToRoiLabel
    
    SrcExpLabel : string
        The Source DICOM series / XNAT experiment label.
    
    TrgExpLabel : string
        The Target DICOM series / XNAT experiment label.

    
    Note:
    ----
        While a linear (or BSpline) interpolator is appropriate for intensity 
        images, only a NearestNeighbor interpolator is appropriate for binary 
        images (e.g. segmentations) so that no new labels are introduced.
    """
    
    #import os
    #from copy import deepcopy
    
    Dict = CreateDictOfPaths()
    
    """ ******************************************************************* """
    """ TestNums for Testing                                                """
    """ ******************************************************************* """
    
    """ Contours """
    
    if TestNum == 'RR1':
        """ Rts Relationship-preserving copy test 1 """
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        # Alternative inputs to test other possible combinations:
        #SrcRoiLabel = None # copy all ROIs
        #SrcSliceNum = 10 # there are no contours on slice 10
        #SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1 # direct copy to next slice
        #TrgSliceNum = 0
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'RR2':
        """ Rts Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR3':
        """ Rts Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RR4':
        """ Rts Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'RR5':
        """ Rts Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'RR6':
        """ Rts Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR7':
        """ Rts Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'RR8':
        """ Rts Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'RR9':
        """ Rts Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    #elif TestNum == 'RR10':
    #    """ Rts Relationship-preserving copy test 10 """
    #    
    #    SrcDcmKey = 'Ses2_Ser2_DcmDir'
    #    SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
    #    
    #    TrgDcmKey = 'Ses3_Ser6_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    SrcRoiLabel = 'Left eye'
    #    SrcSliceNum = None # copy contours on all slices
    #    TrgSliceNum = None # relationship-preserving copy
    #    
    #    SrcExpLabel = 'Ses2_Ser2'
    #    TrgExpLabel = 'Ses3_Ser6'
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' 3_6'
    #
    #elif TestNum == 'RR11':
    #    """ Rts Relationship-preserving copy test 11 """
    #    
    #    SrcDcmKey = 'Ses2_Ser2_DcmDir'
    #    SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
    #    
    #    TrgDcmKey = 'Ses3_Ser6_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    SrcRoiLabel = 'Right eye'
    #    SrcSliceNum = None # copy contours on all slices
    #    TrgSliceNum = None # relationship-preserving copy
    #    
    #    SrcExpLabel = 'Ses2_Ser2'
    #    TrgExpLabel = 'Ses3_Ser6'
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'RR10':
        """ Rts Relationship-preserving copy test 10 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR11':
        """ Rts Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'RR12':
        """ Rts Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses5_Ser12'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR13':
        """ Rts Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR14':
        """ Rts Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR15':
        """ Rts Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses5_Ser12'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
        
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmKey = 'Ses4_Ser13_DcmDir'
    #    SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
    #    
    #    TrgDcmKey = 'MR12_Ser3_DcmDir'
    #    TrgDcmKey = 'MR12_Ser4_DcmDir'
    #    TrgDcmKey = 'MR12_Ser8_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    SrcRoiLabel = 'Nasal cavity'
    #    SrcSliceNum = None # copy contours on all slices
    #    TrgSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
    #    TxtToAddToRoiLabel = ' MR12 S8'
    
    #elif TestNum == 'RR16':
    # """ Rts Relationship-preserving copy test 16 """
    #
    #    SrcDcmKey = 'Ses4_Ser11_DcmDir'
    #    SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
    #    
    #    TrgDcmKey = 'MR12_Ser3_DcmDir'
    #    TrgDcmKey = 'MR12_Ser4_DcmDir'
    #    TrgRoiKey = None
    #    
    #    # Inputs required for this test:
    #    SrcRoiLabel = 'Ventricles'
    #    SrcSliceNum = None # copy contours on all slices
    #    TrgSliceNum = None # relationship-preserving copy
    #    
    #    # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
    #    TxtToAddToRoiLabel = ' MR12 S3'
    #    TxtToAddToRoiLabel = ' MR12 S4'
     
    if TestNum == 'RD1':
        """ Rts Direct copy test 1 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser2_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        #SrcSliceNum = 6
        SrcSliceNum = 5
        TrgSliceNum = SrcSliceNum + 1
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'RD2':
        """ Rts Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        TrgSliceNum = SrcSliceNum + 1
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'RD3':
        """ Rts Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RD4':
        """ Rts Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'RD5':
        """ Rts Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses6_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
    elif TestNum == 'RD6':
        """ Rts Direct copy test 6 - Direct copy analogue of RR2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = 18 
        TrgSliceNum = 21
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    
    """ Segmentations """
    
    
    if TestNum == 'SR1':
        """ Seg Relationship-preserving copy test 1 """
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy segmentations on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'SR2':
        """ Seg Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR3':
        """ Seg Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SR4':
        """ Seg Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'SR5':
        """ Seg Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'SR6':
        """ Seg Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR7':
        """ Seg Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'SR8':
        """ Seg Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'SR9':
        """ Seg Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    
    elif TestNum == 'SR10':
        """ Seg Relationship-preserving copy test 10 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser11'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR11':
        """ Seg Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'SR12':
        """ Seg Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses5_Ser12'
        TrgExpLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR13':
        """ Seg Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Ventricles'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR14':
        """ Seg Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Nasal cavity'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses4_Ser13'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR15':
        """ Seg Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Left eye'
        SrcSliceNum = None # copy contours on all slices
        TrgSliceNum = None # relationship-preserving copy
        
        SrcExpLabel = 'Ses5_Ser12'
        TrgExpLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
        
     
    if TestNum == 'SD1':
        """ Seg Direct copy test 1 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser2_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        #SrcSliceNum = 6 # <-- no segmentation on slice 6?
        SrcSliceNum = 5
        TrgSliceNum = SrcSliceNum + 1
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'SD2':
        """ Seg Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        TrgSliceNum = SrcSliceNum + 1
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'SD3':
        """ Seg Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SD4':
        """ Seg Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'SD5':
        """ Seg Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'Right eye'
        SrcSliceNum = 5
        #TrgSliceNum = SrcSliceNum + 1
        TrgSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcExpLabel = 'Ses2_Ser2'
        TrgExpLabel = 'Ses6_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
    
    
    
    
    
    
    
        """ ************************************************************** """
        """ ************************************************************** """
        """ TestNums for Development                                       """
        """ ************************************************************** """
        """ ************************************************************** """
    
    elif TestNum in ['1R', '1S']:
        """ UseCase 1 (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 23
        TrgSliceNum = 29
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}_s{TrgSliceNum}'
        
    elif TestNum in ['2aR', '2aS']:
        """ UseCase 2a (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser5_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser5_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser5_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 23
        TrgSliceNum = 22
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}_s{TrgSliceNum}'
        
    elif TestNum in ['2bR', '2bS']:
        """ UseCase 2b (relationship-preserving copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser5_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser5_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser5_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 27
        TrgSliceNum = None
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}'
    
    elif TestNum in ['3aiR', '3aiS']:
        """ UseCase 3ai (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser3_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser3_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser3_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 25
        TrgSliceNum = 20
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}_s{TrgSliceNum}'
    
    elif TestNum in ['3aiiR', '3aiiS']:
        """ UseCase 3aii (direct copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser3_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser3_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser3_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 11
        TrgSliceNum = 26
        
        SrcExpLabel = 'MR4_Ser3'
        TrgExpLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}_s{TrgSliceNum}'
    
    elif TestNum in ['3biR', '3biS']:
        """ UseCase 3bi (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR4_Ser3_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser3_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser3_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 26
        TrgSliceNum = None
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}'
    
    elif TestNum in ['3biiR', '3biiS']:
        """ UseCase 3bii (relationship-preserving copy) Rts or Seg """
        
        SrcDcmKey = 'MR4_Ser3_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser3_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser3_SegFpath'
        
        TrgDcmKey = 'MR4_Ser9_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                TrgRoiKey = 'MR4_Ser9_RtsFpath'
            else:
                TrgRoiKey = 'MR4_Ser9_SegFpath'
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 11
        TrgSliceNum = None
        
        SrcExpLabel = 'MR4_Ser3'
        TrgExpLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}'
    
    elif TestNum in ['5bR', '5bS']:
        """ UseCase 5b (relationship-preserving copy) Rts or Seg"""
        
        SrcDcmKey = 'MR4_Ser9_DcmDir'
        if 'R' in TestNum:
            SrcRoiKey = 'MR4_Ser9_RtsFpath'
        else:
            SrcRoiKey = 'MR4_Ser9_SegFpath'
        
        TrgDcmKey = 'MR12_Ser8_DcmDir'
        if UseOrigTrgRoi:
            if 'R' in TestNum:
                raise Exception('A RTS does not yet exist for MR12 Ser8.')
            else:
                raise Exception('A SEG does not yet exist for MR12 Ser8.')
        else:
            TrgRoiKey = None
        
        # Inputs required for this test:
        SrcRoiLabel = 'tumour'
        SrcSliceNum = 23
        TrgSliceNum = None
        
        SrcExpLabel = 'MR4_Ser9'
        TrgExpLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{SrcSliceNum}_to_{TrgExpLabel}'
    
    
    
    """ Use the keys to select the directories and file paths from Dict. """
    SrcDcmDir = Dict[SrcDcmKey]
    SrcRoiFpath = Dict[SrcRoiKey]
    
    TrgDcmDir = Dict[TrgDcmKey]
    if TrgRoiKey:
        TrgRoiFpath = Dict[TrgRoiKey]
    else:
        TrgRoiFpath = None
    
    
    
    return SrcRoiFpath, SrcSliceNum, SrcRoiLabel, TrgSliceNum, SrcDcmDir,\
           TrgDcmDir, TrgRoiFpath, TxtToAddToRoiLabel, SrcExpLabel, TrgExpLabel