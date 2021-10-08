# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 19:07:06 2021

@author: ctorti
"""

from importlib import reload
import dro_tools.create_spa_dro
import dro_tools.create_def_dro
reload(dro_tools.create_spa_dro)
reload(dro_tools.create_def_dro)

from dro_tools.create_spa_dro import create_spa_dro_from_tx
from dro_tools.create_def_dro import create_def_dro_from_tx

def create_dro(srcDataset, trgDataset, newDataset, params):
    # TODO update docstrings
    """
    Create a spatial or deformable DICOM Registration Object if applicable.
    
    If a transform exists as a result of image registration a DRO will be
    created. Otherwise None will be returned. The DRO will be added to 
    propObj.
    
    Parameters
    ----------
    srcDataset : Importer Object, optional
        Object containing various data for the source DICOM series.
    trgDataset : Importer Object, optional
        Object containing various data for the target DICOM series.
    newDataset : Propagator Object
        Object containing various data for the source-to-target propagated (or 
        copied) dataset.
    params : DataDownloader Object
        Contains parameters (cfgDict), file paths (pathsDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    dro : Pydicom Object or None
        DICOM (spatial or deformable) Registration Object if applicable, None
        otherwise.
    """
    
    """ 
    Consider adding this as an input variable for the user to specify.
    Empty for now.
    """
    # ContentDescription:
    contentDesc = ''
        
    useCaseToApply = params.cfgDict['useCaseToApply']
    # Name of the registration transform:
    regTxName = params.cfgDict['regTxName']
    # The registration and pre-registration transforms:
    regTx = newDataset.finalTx
    preRegTx = newDataset.preRegTx
    sampleDroDir = params.cfgDict['sampleDroDir']
    p2c = params.cfgDict['p2c']
    
    if useCaseToApply in ['5a', '5b']:
        timingMsg = "Creating the DICOM Registration Object...\n"
        params.add_timestamp(timingMsg)
        
        if p2c:
            print('regTx =', print(regTx), '\n')
        
        if regTxName in ['rigid', 'affine']:
            # Create a Spatial Registration DRO.
            dro = create_spa_dro_from_tx(
                srcDataset, 
                trgDataset,
                newDataset, 
                params, 
                description=contentDesc
                )
            # TODO Consider adding preRegTx (but need to return it from create_spa_dro_from_tx..) (09/07/21)
        else:
            # Create a Deformable Spatial Registration DRO.
            dro = create_def_dro_from_tx(
                srcDataset, 
                trgDataset,
                newDataset, 
                params, 
                description=contentDesc
                )
        
        timingMsg = "Took [*] s to create the DICOM Registration Object.\n"
        params.add_timestamp(timingMsg)
    else:
        dro = None
    
    return dro