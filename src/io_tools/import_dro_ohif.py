# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 21:11:08 2021

@author: ctorti
"""


""" 
This replaces import_dro.py for the OHIF-integration, which expects a DRO in
src/inputs/dro rather than searching for one on XNAT using the REST API.

THIS HAS NOT BEEN TESTED YET!!! (08/11/2021)
"""

import os
from pydicom import dcmread
from datetime import datetime
#from xnat_tools.FOR_uid import get_FOR_uid
#from xnat_tools.dicom_metadata import get_dicom_metadata
from dicom_tools.dcm_metadata import get_dcm_uids
from image_tools.attrs_info import get_im_attrs
    

class DroImporterOhif:
    """
    This class creates a DroImporter Object by importing a DICOM Registration 
    Object (DRO) from src/inputs/dro.
    
    Parameters
    ----------
    cfgObj : ConfigFromOhif Object
        Contains parameters and directory/file paths (cfgDict), timestamps
        (timings) and timing messages (timingMsgs).
    
    Returns
    -------
    self.dro : Pydicom Object or None
        Pydicom representation of the DRO if a suitable DRO was found, and None
        otherwise.
    self.txMatrix : list of strs or None
        List of float strings representing the non-deformable transformation 
        that transforms the moving (Source) image to the fixed (Target) image 
        (if there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    self.gridDims : list of ints or None
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image (if 
        there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    self.gridRes : list of floats or None
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image (if there 
        exists a DRO (as a subject assessor) that matches the requirements); 
        None if not.
    self.vectGridData : list of floats or None
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image (if there exists a
        DRO (as a subject assessor) that matches the requirements); None if not.
    """
    
    """
    def __init__(self, params):
        useCaseToApply = params.useCaseToApply
        
        # The following were moved to fetch_dro:
        ## Initialise values:
        #self.dro = None
        #self.txMatrix = None
        #self.gridDims = None
        #self.gridRes = None
        #self.vectGridData = None
        
        #useCaseToApply = params.useCaseToApply
        
        # If image registration is required (useCaseToApply = '5a' or '5b')
        # search XNAT for a suitable DRO for this set of Source and Target:
        if '5' in useCaseToApply:
            self.fetch_dro(params)
            
        print(f'type(self.dro) = {type(self.dro)}\n')
    """
    
    def import_dro(self, cfgObj):
        """
        Import DRO from src/inputs/dro if a suitable one is found.

        Parameters
        ----------
        cfgObj : ConfigFromOhif Object
            Contains parameters and directory/file paths (cfgDict), timestamps
            (timings) and timing messages (timingMsgs).

        Returns
        -------
        self.dro : Pydicom Object or None
            Pydicom representation of the DRO if a suitable DRO was found, and
            None otherwise.
        self.txMatrix : list of strs or None
            List of float strings representing the non-deformable  
            transformation that transforms the moving (Source) image to the 
            fixed (Target) image if dro != None; None otherwise.
        self.gridDims : list of ints or None
            List of integers representing the dimensions of the BSpline grid  
            used to deform the moving (Source) image to the fixed (Target)  
            image if dro != None; None otherwise.
        self.gridRes : list of floats or None
            List of floats representing the resolution of the BSpline grid used 
            to deform the moving (Source) image to the fixed (Target) image if
            dro != None; None otherwise.
        self.vectGridData : list of floats or None
            List of floats representing the vector deformations that deform the
            moving (Source) image to the fixed (Target) image if dro != None; 
            None otherwise.
        
        Note
        ----
        Modified from fetch_dro() in io_tools.import_dro.py
        """
        
        cfgDict = cfgObj.cfgDict
        droDir = cfgDict['droDir']
        srcDicomDir = cfgDict['srcDicomDir']
        trgDicomDir = cfgDict['trgDicomDir']
        #projID = cfgDict['projID']
        #subjLab = cfgDict['subjLab']
        srcExpLab = cfgDict['srcExpLab']
        srcScanID = cfgDict['srcScanID']
        trgExpLab = cfgDict['trgExpLab']
        trgScanID = cfgDict['trgScanID']
        regTxName = cfgDict['regTxName']
        p2c = cfgDict['p2c']
        
        if p2c:
            print('\n\n', '-'*120)
            print('Running of parse_dro():')
            print('\n\n', '-'*120)
        
        """
        srcFORuid = get_FOR_uid(
            url, projID, subjLab, srcExpLab, srcScanID, xnatSession
            )
        
        trgFORuid = get_FOR_uid(
            url, projID, subjLab, trgExpLab, trgScanID, xnatSession
            )
        """
        
        """
        srcStudyUID_req, srcSeriesUID_req, srcFORuid_req, srcIPPs_req,\
            srcIOP_req = get_dicom_metadata(
                url, projID, subjLab, srcExpLab, srcScanID, xnatSession
                )
        
        trgStudyUID_req, trgSeriesUID_req, trgFORuid_req, trgIPPs_req,\
            trgIOP_req = get_dicom_metadata(
                url, projID, subjLab, trgExpLab, trgScanID, xnatSession
                )
        """
        
        srcStudyUID_req, srcSeriesUID_req, srcFORuid_req,\
            srcSOPuids = get_dcm_uids(srcDicomDir)
        
        trgStudyUID_req, trgSeriesUID_req, trgFORuid_req,\
            trgSOPuids = get_dcm_uids(trgDicomDir)
        
        srcImSize, srcImSpacings, srcSlcThick, srcIPPs_req, srcDirection,\
            srcWarnings = get_im_attrs(srcDicomDir)
        
        trgImSize, trgImSpacings, trgSlcThick, trgIPPs_req, trgDirection,\
            trgWarnings = get_im_attrs(trgDicomDir)
        
        # srcDirection and trgDirection are the direction cosines along x, y
        # and z, so take the first 6 elements for the IOPs:
        srcIOP_req = srcDirection[0:6]
        trgIOP_req = trgDirection[0:6]
        
        if regTxName == 'bspline':
            droType_req = 'Deformable Spatial Registration Storage'
        else:
            droType_req = 'Spatial Registration Storage'
        
        if p2c:
            print('Items to match:')
            print(f'   droType_req = {droType_req}')
            print(f'   srcFORuid_req = {srcFORuid_req}')
            print(f'   trgFORuid_req = {trgFORuid_req}\n')
            print('Items expected to match:')
            print(f'   srcStudyUID_req = {srcStudyUID_req}')
            print(f'   trgStudyUID_req = {trgStudyUID_req}\n')
            #print(f'   srcIOP_req = {srcIOP_req}')
            #print(f'   trgIOP_req = {trgIOP_req}\n')
        
        
        # Get list of files in droDir containing .dcm extension:
        fnames = [
            file for file in os.listdir(droDir) if file.endswith('.dcm')
            ]
        
        if p2c:
            print(f'There were {len(fnames)} DICOMs found in {droDir}:\n')
        
        # Initialise lists:
        txMatrices = []
        txMatrixTypes = []
        gridDimss = [] # note extra s since a list of gridDims
        gridRess = [] # note extra s since a list of gridRes
        vectGridDatas = []
        
        # The indices of any DROs whose FORuids match the Source and Target FORuids
        # and whose SOPClassUID is compatible with regTxName: 
        inds = []
        
        # files will be a list of any DROs that match; contentDateTimes are their
        # corresponding ContentDate and ContentTime DICOM tags combined:
        contentDateTimes = []
        files = []
        
        for i in range(len(fnames)):
            # Read in the file (which may or may not be a DRO):
            file = dcmread(os.join.path(droDir, fnames[i]))
            
            files.append(file)
            
            mod = f'{file.Modality}'
            
            # Only proceed if the DICOM file is a registration object:
            if mod == 'REG': # 07/05/21
                """ Get meta data """
                
                contentDateTime = f'{file.ContentDate}{file.ContentTime}'
                seriesDesc = f'{file.SeriesDescription}'
                contentDesc = f'{file.ContentDescription}'
                
                #droType = f'{file.SOPClassUID}' # this is a number 
                # (e.g. 1.2.840.10008.5.1.4.1.1.66.3)
                droType = f'{file[0x0008, 0x0016].repval}' # this is the 
                # representation of the element's value
                
                # Use f-strings instead of deepcopy:
                if 'Deformable' in droType:
                    trgFORuid = f'{file.DeformableRegistrationSequence[0].SourceFrameOfReferenceUID}'
                    srcFORuid = f'{file.DeformableRegistrationSequence[1].SourceFrameOfReferenceUID}'
                else:
                    trgFORuid = f'{file.RegistrationSequence[0].FrameOfReferenceUID}'
                    srcFORuid = f'{file.RegistrationSequence[1].FrameOfReferenceUID}'
                
                trgStudyUID = f'{file.StudyInstanceUID}'
                
                if 'Deformable' in droType:
                    gridDims = file.DeformableRegistrationSequence[1]\
                                   .DeformableRegistrationGridSequence[0]\
                                   .GridDimensions
                    
                    gridRes = file.DeformableRegistrationSequence[1]\
                                  .DeformableRegistrationGridSequence[0]\
                                  .GridResolution
                    
                    vectGridData = file.DeformableRegistrationSequence[1]\
                                       .DeformableRegistrationGridSequence[0]\
                                       .VectorGridData
                else:
                    txMatrix = file.RegistrationSequence[1]\
                                   .MatrixRegistrationSequence[0]\
                                   .MatrixSequence[0]\
                                   .FrameOfReferenceTransformationMatrix
                    
                    txMatrixType = file.RegistrationSequence[1]\
                                       .MatrixRegistrationSequence[0]\
                                       .MatrixSequence[0]\
                                       .FrameOfReferenceTransformationMatrixType
                
                
                """ Does meta data match requirements? """
                
                """
                if 'Deformable' in droType:
                    if regTxName == 'bspline':
                        matchOnDroType = True
                    else:
                        matchOnDroType = False
                else:
                    if regTxName in ['rigid', 'affine']:
                        matchOnDroType = True
                    else:
                        matchOnDroType = False
                """
                
                if droType == droType_req:
                    matchOnDroType = True
                else:
                    matchOnDroType = False
                
                if srcFORuid == srcFORuid_req:
                    matchOnSrcFOR = True
                else:
                    matchOnSrcFOR = False
                
                if trgFORuid == trgFORuid_req:
                    matchOnTrgFOR = True
                else:
                    matchOnTrgFOR = False
                
                if trgStudyUID == trgStudyUID_req:
                    matchOnTrgStudyUID = True
                else:
                    matchOnTrgStudyUID = False
                
                if p2c:
                    print(f'  File {i+1} of {len(fnames)}')
                    print(f'    fname = {fnames[i]}')
                    print(f'    contentDateTime = {contentDateTime}')
                    print(f'    seriesDesc = {seriesDesc}')
                    print(f'    contentDesc = {contentDesc}')
                    print(f'    droType = {droType}')
                    if 'Deformable' in droType:
                        print(f'    gridDims = {gridDims}')
                        print(f'    gridRes = {gridRes}')
                        print(f'    len(vectGridData) = {len(vectGridData)}')
                    else:
                        print(f'    txMatrixType = {txMatrixType}')
                        """
                        print(f'    txMatrix = [{float(txMatrix[0])},',
                              f'{float(txMatrix[1])}, {float(txMatrix[2])},',
                              f'{float(txMatrix[3])},')
                        print(f'                {float(txMatrix[4])},',
                              f'{float(txMatrix[5])}, {float(txMatrix[6])},',
                              f'{float(txMatrix[7])},')
                        print(f'                {float(txMatrix[8])},',
                              f'{float(txMatrix[9])}, {float(txMatrix[10])},',
                              f'{float(txMatrix[11])},')
                        print(f'                {float(txMatrix[12])},',
                              f'{float(txMatrix[13])}, {float(txMatrix[14])},',
                              f'{float(txMatrix[15])}]')
                        """
                        txt = '    txMatrix = ['
                        for j in range(4):
                            txt += f'{txMatrix[4*j]}, {txMatrix[4*j + 1]}, '\
                                + f'{txMatrix[4*j + 2]}, {txMatrix[4*j + 3]}'
                            if j < 3:
                                txt += ',\n                '
                            else:
                                txt += ']'
                        print(txt)
                    print(f'    matchOnDroType = {matchOnDroType}')
                    print(f'    matchOnSrcFOR = {matchOnSrcFOR}')
                    print(f'    matchOnTrgFOR = {matchOnTrgFOR}')
                    print(f'    matchOnTrgStudyUID = {matchOnTrgStudyUID}\n')
                
                
                """ Store meta data that match requirements in lists """
                
                if (matchOnDroType and matchOnSrcFOR and matchOnTrgFOR):
                    inds.append(i)
                    
                    if 'Deformable' in droType:
                        gridDimss.append(gridDims)
                        
                        gridRess.append(gridRes)
                        
                        vectGridDatas.append(vectGridData)
                    else:
                        txMatrices.append(txMatrix)
                        
                        txMatrixTypes.append(txMatrixType)
                    
                    # Convert contentDateTime to datetime object:
                    try:
                        contentDateTime = datetime.strptime(
                            contentDateTime, '%Y%m%d%H%M%S.%f'
                            )
                    except ValueError:
                        contentDateTime = datetime.strptime(
                            contentDateTime, '%Y%m%d%H%M%S'
                            )
                    contentDateTimes.append(contentDateTime)
        
        if p2c:
            print(f'  The indices of the files that match = {inds}')
            #print(f'  contentDateTimes = {contentDateTimes}')
            txt = '  contentDateTimes = ['
            for j in range(len(contentDateTimes)):
                txt += f'{contentDateTimes[j]}'
                if j < len(contentDateTimes) - 1:
                    txt += ',\n                      '
                else:
                    txt += ']'
            print(txt)
        
        if contentDateTimes:#time_diffs:
            # The index of the most recent contentDateTime:
            newInd = contentDateTimes.index(max(contentDateTimes))
            
            # The index of the most recent file that matches:
            droInd = inds[newInd]
            
            if p2c:
                print(f'\n  The index of the most recent file that matches = {droInd}')
            
            dro = files[droInd]
            
            if regTxName == 'bspline':
                txMatrix = None
                gridDims = gridDimss[newInd]
                gridRes = gridRess[newInd]
                vectGridData = vectGridDatas[newInd]
                
                gridResRounded = [round(item, 2) for item in gridRes]
                
                msg = f"* There were {len(fnames)} subject assessors found,"+\
                    f" of which {len(contentDateTimes)} were a {droType}, " +\
                    "with\n  FrameOfReferenceUIDs matching those of the " +\
                    "Source and Target image series, matching the transform "+\
                    "'{regTxName}', \n  the most recent of which contains " +\
                    f"the grid dimensions {gridDims}, grid resolution " +\
                    f"{gridResRounded}, and vector grid data containing " +\
                    f"{len(vectGridData)} elements.\n"
            
            if regTxName in ['rigid', 'affine']:
                txMatrix = txMatrices[newInd]
                gridDims = None
                gridRes = None
                vectGridData = None
                
                txMatrixRounded = [round(item, 3) for item in txMatrix]
                
                msg = f"* There were {len(fnames)} subject assessors found," +\
                    f" of which {len(contentDateTimes)} were a {droType}, " +\
                    "with\n  FrameOfReferenceUIDs matching those of the " +\
                    "Source and Target image series, matching the transform "+\
                    f"'{regTxName}', \n  the most recent of which contains " +\
                    f"the transformation matrix: \n  {txMatrixRounded}\n"
        else:
            dro = None
            txMatrix = None
            gridDims = None
            gridRes = None
            vectGridData = None
            
            msg = "* There were no DROs with FrameOfReferenceUIDs matching " +\
                f"those of the Source ({srcFORuid_req}) \nand Target " +\
                f"({trgFORuid_req}) image series, and whose SOP Class " +\
                f"matched the desired transform type ({regTxName}).\n"
        print(msg)
        
        self.dro = dro
        self.txMatrix = txMatrix
        self.gridDims = gridDims
        self.gridRes = gridRes
        self.vectGridData = vectGridData
        
        if p2c:
            print('-'*120)