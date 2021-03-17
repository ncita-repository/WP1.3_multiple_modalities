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
    Create a configuration file assigning default values for some inputs 
    required for CopyRoi().
    
    The file will be exported to the current working directory.
    """
    
    import os
    from pypref import Preferences
    
    cwd = os.getcwd()
    
    Dict = {'TrgAsrMod' : None, 
            'TrgAsrName' : None, 
            'TrgSliceNum' : None, 
            'TxtToAddToTrgAsrName' : '', 
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
            'Tx' : 'affine', 
            'TxMaxIters' : '512', 
            'TxInterp' : 'NearestNeighbor', 
            'ApplyPostTxBin' : True,
            'ApplyPostTxBlur' : True, 
            'PostTxVariance' : (1,1,1),
            'ExportNewTrgAsr' : True,
            'RtsExportDir' : os.path.join(cwd, 'new_RTS'),
            'SegExportDir' : os.path.join(cwd, 'new_SEG'),
            'PlotAllSlices' : False,
            'RtsPlotExportDir' : os.path.join(cwd, 'plots_RTS'),
            'SegPlotExportDir' : os.path.join(cwd, 'plots_SEG'),
            'LogExportDir' : os.path.join(cwd, 'logs')
            }
    
    pref = Preferences(directory=cwd, filename='CopyRoiDefaults.py')

    pref.set_preferences(Dict)
    
    return




def CreateTestConfigFile():
    """ 
    Create a configuration file assigning values for inputs required for 
    CopyRoi() based on the tests in 
    https://docs.google.com/document/d/1_LxC1FUdQME4xuXZFGjjGBGYRb1cJX8YiEADjtZTn_g/edit#.
    
    The file will be exported to the current working directory.
    """
    
    import os
    from pypref import Preferences
    
    Dict = {'XnatUrl' : 'http://10.1.1.20', 
            'ProjId' : 'NCITA_TEST', 
            'SubjLabel' : 'TCGA-BB-A5HY',
            
            'RR1': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '4', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 4'
                    },
            
            'RR2': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '10', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 10'
                    },
            
            'RR3': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 5'
                    },
            
            'RR4': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '11', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 11'
                    },
            
            'RR5': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '14', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 14'
                    },
                    
            'RR6': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '10', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 10'
                    },
            
            'RR7': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '12', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 12'
                    },
            
            'RR8': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '13', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 13'
                    },
            
            'RR9': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '14', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 14'
                    },
            
            'RR10': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '11', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'RR11': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '13', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'RR12': {'SrcStudyLabel' : 'Session5', 
                     'SrcSeriesLabel' : '12', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'RR13': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '11', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            'RR14': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '13', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            'RR15': {'SrcStudyLabel' : 'Session5', 
                     'SrcSeriesLabel' : '12', 
                     'SrcAsrMod' : 'RTSTRUCT', 
                     'SrcAsrName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'RTSTRUCT', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            #'RR10r': {'SrcStudyLabel' : 'Session3',
            #          'SrcSeriesLabel' : '6', 
            #          'SrcAsrMod' : 'RTSTRUCT', 
            #          'SrcAsrName' : 'Ventricles', 
            #          'SrcRoiName' : None,
            #          'SrcSliceNum' : None, 
            #          'TrgStudyLabel' : 'Session4', 
            #          'TrgSeriesLabel' : '11', 
            #          'TrgAsrMod' : 'RTSTRUCT', 
            #          'TrgAsrName' : None, 
            #          'TrgRoiName' : None,
            #          'TrgSliceNum' : None,
            #          'TxtToAddToTrgAsrName' : ' 4_11'
            #         },
            
            'RD1': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '2', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgAsrName' : ' 2'
                    },
            
            'RD2': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '4', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgAsrName' : ' 4'
                    },
            
            'RD3': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 65, # Source slice 5 maps to Target slices 55-64
                    'TxtToAddToTrgAsrName' : ' 5'
                    },
            
            'RD4': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session6',
                    'TrgSeriesLabel' : '3', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 8, # Source slice 5 maps to Target slice 7
                    'TxtToAddToTrgAsrName' : ' 6_3'
                    },
            
            'RD5': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session6',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 78, # Source slice 5 maps to Target slices 67-77
                    'TxtToAddToTrgAsrName' : ' 6_5'
                    },
            
            'RD6': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'RTSTRUCT', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 18, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '10', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 21, 
                    'TxtToAddToTrgAsrName' : ' 10'
                    },
                    
                    
            'SR1': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '4', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 4'
                    },
            
            'SR2': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '10', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 10'
                    },
            
            'SR3': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Left eye CT', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 5'
                    },
            
            'SR4': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '11', 
                    'TrgAsrMod' : 'RTSTRUCT', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 11'
                    },
            
            'SR5': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '13', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Nasal cavity', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '14', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 14'
                    },
                    
            'SR6': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '10', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 10'
                    },
            
            'SR7': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '12', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 12'
                    },
            
            'SR8': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '13', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 13'
                    },
            
            'SR9': {'SrcStudyLabel' : 'Session4', 
                    'SrcSeriesLabel' : '11', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Ventricles', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : None, 
                    'TrgStudyLabel' : 'Session4',
                    'TrgSeriesLabel' : '14', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : None,
                    'TxtToAddToTrgAsrName' : ' 14'
                    },
            
            'SR10': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '11', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'SR11': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '13', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'SR12': {'SrcStudyLabel' : 'Session5', 
                     'SrcSeriesLabel' : '12', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session3',
                     'TrgSeriesLabel' : '6', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 3_6'
                     },
            
            'SR13': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '11', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Ventricles', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            'SR14': {'SrcStudyLabel' : 'Session4', 
                     'SrcSeriesLabel' : '13', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Nasal cavity', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            'SR15': {'SrcStudyLabel' : 'Session5', 
                     'SrcSeriesLabel' : '12', 
                     'SrcAsrMod' : 'SEG', 
                     'SrcAsrName' : 'Left eye', 
                     'SrcRoiName' : None,
                     'SrcSliceNum' : None, 
                     'TrgStudyLabel' : 'Session1',
                     'TrgSeriesLabel' : '5', 
                     'TrgAsrMod' : 'SEG', 
                     'TrgAsrName' : None, 
                     'TrgRoiName' : None,
                     'TrgSliceNum' : None,
                     'TxtToAddToTrgAsrName' : ' 1_5'
                     },
            
            'SD1': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '2', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgAsrName' : ' 2'
                    },
            
            'SD2': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '4', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 6,
                    'TxtToAddToTrgAsrName' : ' 4'
                    },
            
            'SD3': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session2',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 65, # Source slice 5 maps to Target slices 55-64
                    'TxtToAddToTrgAsrName' : ' 5'
                    },
            
            'SD4': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session6',
                    'TrgSeriesLabel' : '3', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 8, # Source slice 5 maps to Target slice 7
                    'TxtToAddToTrgAsrName' : ' 6_3'
                    },
            
            'SD5': {'SrcStudyLabel' : 'Session2', 
                    'SrcSeriesLabel' : '2', 
                    'SrcAsrMod' : 'SEG', 
                    'SrcAsrName' : 'Right eye CT (single slice)', 
                    'SrcRoiName' : None,
                    'SrcSliceNum' : 5, 
                    'TrgStudyLabel' : 'Session6',
                    'TrgSeriesLabel' : '5', 
                    'TrgAsrMod' : 'SEG', 
                    'TrgAsrName' : None, 
                    'TrgRoiName' : None,
                    'TrgSliceNum' : 78, # Source slice 5 maps to Target slices 67-77
                    'TxtToAddToTrgAsrName' : ' 6_5'
                    }
            }
    
    pref = Preferences(directory=os.getcwd(), filename='CopyRoiTestConfig.py')

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
    
    Image : SimpleITK image
        The 3D image with modified pixel type.

    
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
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Alternative inputs to test other possible combinations:
        #FromRoiLabel = None # copy all ROIs
        #FromSliceNum = 10 # there are no contours on slice 10
        #FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1 # direct copy to next slice
        #ToSliceNum = 0
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'RR2':
        """ Rts Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR3':
        """ Rts Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RR4':
        """ Rts Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'RR5':
        """ Rts Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'RR6':
        """ Rts Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'RR7':
        """ Rts Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'RR8':
        """ Rts Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'RR9':
        """ Rts Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
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
    #    FromRoiLabel = 'Left eye'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    SrcLabel = 'Ses2_Ser2'
    #    TrgLabel = 'Ses3_Ser6'
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
    #    FromRoiLabel = 'Right eye'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
    #    
    #    SrcLabel = 'Ses2_Ser2'
    #    TrgLabel = 'Ses3_Ser6'
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
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR11':
        """ Rts Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'RR12':
        """ Rts Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'RR13':
        """ Rts Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR14':
        """ Rts Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'RR15':
        """ Rts Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_RtsFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses1_Ser5'
        
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
    #    FromRoiLabel = 'Nasal cavity'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
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
    #    FromRoiLabel = 'Ventricles'
    #    FromSliceNum = None # copy contours on all slices
    #    ToSliceNum = None # relationship-preserving copy
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
        FromRoiLabel = 'Right eye'
        #FromSliceNum = 6
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'RD2':
        """ Rts Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'RD3':
        """ Rts Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'RD4':
        """ Rts Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'RD5':
        """ Rts Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_RtsFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_5'
        
    elif TestNum == 'RD6':
        """ Rts Direct copy test 6 - Direct copy analogue of RR2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_RtsFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = 18 
        ToSliceNum = 21
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
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
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy segmentations on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
    
    elif TestNum == 'SR2':
        """ Seg Relationship-preserving copy test 2 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR3':
        """ Seg Relationship-preserving copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_LeftEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SR4':
        """ Seg Relationship-preserving copy test 4 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser11_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser11'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 11'
    
    elif TestNum == 'SR5':
        """ Seg Relationship-preserving copy test 5 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses4_Ser14'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    elif TestNum == 'SR6':
        """ Seg Relationship-preserving copy test 6 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser10_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser10'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 10'
        
    elif TestNum == 'SR7':
        """ Seg Relationship-preserving copy test 7 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser12_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser12'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 12'
    
    elif TestNum == 'SR8':
        """ Seg Relationship-preserving copy test 8 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser13_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser13'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 13'
    
    elif TestNum == 'SR9':
        """ Seg Relationship-preserving copy test 9 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses4_Ser14_DcmDir'
        TrgRoiKey = None
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses4_Ser14'
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 14'
    
    
    elif TestNum == 'SR10':
        """ Seg Relationship-preserving copy test 10 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser11'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR11':
        """ Seg Relationship-preserving copy test 11 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
    
    elif TestNum == 'SR12':
        """ Seg Relationship-preserving copy test 12 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses3_Ser6_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses3_Ser6'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 3_6'
        
    elif TestNum == 'SR13':
        """ Seg Relationship-preserving copy test 13 """
        
        SrcDcmKey = 'Ses4_Ser11_DcmDir'
        SrcRoiKey = 'Ses4_Ser11_Ventricles_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Ventricles'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR14':
        """ Seg Relationship-preserving copy test 14 """
        
        SrcDcmKey = 'Ses4_Ser13_DcmDir'
        SrcRoiKey = 'Ses4_Ser13_NasalCavity_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Nasal cavity'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses4_Ser13'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
    
    elif TestNum == 'SR15':
        """ Seg Relationship-preserving copy test 15 """
        
        SrcDcmKey = 'Ses5_Ser12_DcmDir'
        SrcRoiKey = 'Ses5_Ser12_LeftEye_SegFpath'
        
        TrgDcmKey = 'Ses1_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Left eye'
        FromSliceNum = None # copy contours on all slices
        ToSliceNum = None # relationship-preserving copy
        
        SrcLabel = 'Ses5_Ser12'
        TrgLabel = 'Ses1_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 1_5'
        
     
    if TestNum == 'SD1':
        """ Seg Direct copy test 1 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser2_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        #FromSliceNum = 6 # <-- no segmentation on slice 6?
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser2'
    
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 2'
    
    elif TestNum == 'SD2':
        """ Seg Direct copy test 2 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser4_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        ToSliceNum = FromSliceNum + 1
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser4'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 4'
        
    elif TestNum == 'SD3':
        """ Seg Direct copy test 3 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses2_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 65 # Source slice 5 maps to Target slices 55-64
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses2_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 5'
    
    elif TestNum == 'SD4':
        """ Seg Direct copy test 4 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser3_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 8 # Source slice 5 maps to Target slice 7
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = ' 6_3'
        
    elif TestNum == 'SD5':
        """ Seg Direct copy test 5 """
        
        SrcDcmKey = 'Ses2_Ser2_DcmDir'
        SrcRoiKey = 'Ses2_Ser2_RightEyeCT_SegFpath'
        
        TrgDcmKey = 'Ses6_Ser5_DcmDir'
        TrgRoiKey = None
        
        # Inputs required for this test:
        FromRoiLabel = 'Right eye'
        FromSliceNum = 5
        #ToSliceNum = FromSliceNum + 1
        ToSliceNum = 78 # Source slice 5 maps to Target slices 67-77
        
        SrcLabel = 'Ses2_Ser2'
        TrgLabel = 'Ses6_Ser5'
        
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 29
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = 22
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
        
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 27
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser5'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 25
        ToSliceNum = 20
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = 26
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}_s{ToSliceNum}'
    
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 26
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR4_Ser3'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 11
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser3'
        TrgLabel = 'MR4_Ser9'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
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
        FromRoiLabel = 'tumour'
        FromSliceNum = 23
        ToSliceNum = None
        
        SrcLabel = 'MR4_Ser9'
        TrgLabel = 'MR12_Ser8'
        
        # Text to add to RTS ROIStructureSetLabel or SEG SeriesDescription:
        TxtToAddToRoiLabel = f'_s{FromSliceNum}_to_{TrgLabel}'
    
    
    
    """ Use the keys to select the directories and file paths from Dict. """
    SrcDcmDir = Dict[SrcDcmKey]
    SrcRoiFpath = Dict[SrcRoiKey]
    
    TrgDcmDir = Dict[TrgDcmKey]
    if TrgRoiKey:
        TrgRoiFpath = Dict[TrgRoiKey]
    else:
        TrgRoiFpath = None
    
    
    
    return SrcRoiFpath, FromSliceNum, FromRoiLabel, ToSliceNum, SrcDcmDir,\
           TrgDcmDir, TrgRoiFpath, TxtToAddToRoiLabel, SrcLabel, TrgLabel