# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 12:33:15 2021

@author: ctorti
"""


from pydicom.filebase import DicomBytesIO
from pydicom import dcmread
from datetime import datetime
from xnat_tools.FOR_uid import get_FOR_uid
    

class DroImporter:
    """
    This class creates a droData Object by fetching a DRO from XNAT.
    
    Parameters
    ----------
    params : Params Object
        Object containing various parameters.
    
    Returns
    -------
    
    """
    
    def __init__(self, params):
        useCaseToApply = params.useCaseToApply
        
        # Initialise values:
        self.dro = None
        self.txMatrix = None
        self.gridDims = None
        self.gridRes = None
        self.vectGridData = None
        
        useCaseToApply = params.useCaseToApply
        
        # If image registration is required (useCaseToApply = '5a' or '5b')
        # search XNAT for a suitable DRO for this set of Source and Target:
        if '5' in useCaseToApply:
            self.fetch_dro(params)
            
        print(f'type(self.dro) = {type(self.dro)}\n')
    
    def fetch_dro(self, params):
        """
        Search XNAT for a suitable DRO.
        
        Modified from search_dro().

        Parameters
        ----------
        params : Params Object
            Object containing various parameters.

        Returns
        -------
        None.

        """
        
        xnatSession = params.xnatSession
        url = xnatSession.url
        projID = params.srcXnatParams['projID']
        subjLab = params.srcXnatParams['subjLab']
        srcExpLab = params.srcXnatParams['expLab']
        srcScanID = params.srcXnatParams['scanID']
        trgExpLab = params.trgXnatParams['expLab']
        trgScanID = params.trgXnatParams['scanID']
        regTxName = params.resParams['regTxName']
        p2c = params.genParams['p2c']
        
        srcFORUID = get_FOR_uid(url, projID, subjLab, srcExpLab, srcScanID,
                                xnatSession)
        
        trgFORUID = get_FOR_uid(url, projID, subjLab, trgExpLab, trgScanID,
                                xnatSession)
        
        # Get a listing of all resource files for the subject:
        uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/files'
        
        request = xnatSession.get(uri)
        
        # Raise status error if not None:
        if request.raise_for_status() != None:
            print(request.raise_for_status())
        
        fnames = []
    
        for result in request.json()['ResultSet']['Result']:
            #fnames.append(result['Name']) # 07/05/21
            fname = result['Name'] # 07/05/21
            if '.dcm' in fname: # 07/05/21
                fnames.append(fname) # 07/05/21
        
        if p2c:
            print(f'  There were {len(fnames)} DICOM subject assessors found.')
        
        # Initialise lists:
        txMatrices = []
        #txMatrixTypes = []
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
            uri = f'{url}/data/projects/{projID}/subjects/{subjLab}/'\
                  + f'files/{fnames[i]}'
            
            request = xnatSession.get(uri)
            
            raw = DicomBytesIO(request.content)
            
            # Read in the file (which may or may not be a DRO):
            file = dcmread(raw)
            
            files.append(file)
            
            mod = f'{file.Modality}'
            
            # Only proceed if the DICOM file is a registration object:
            if mod == 'REG': # 07/05/21
                #droType = f'{file.SOPClassUID}' # this is a number 
                # (e.g. 1.2.840.10008.5.1.4.1.1.66.3)
                droType = f'{file[0x0008, 0x0016].repval}' # this is the 
                # representation of the element's value
                
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
                
                if p2c:
                    print(f'  File {i+1} of {len(fnames)}')
                    print(f'    droType = {droType}\n')
                    print(f'    matchOnDroType = {matchOnDroType}\n')
                
                # Use f-strings instead of deepcopy:
                if 'Deformable' in droType:
                    FORUID0 = f'{file.DeformableRegistrationSequence[0].SourceFrameOfReferenceUID}'
                    FORUID1 = f'{file.DeformableRegistrationSequence[1].SourceFrameOfReferenceUID}'
                else:
                    FORUID0 = f'{file.RegistrationSequence[0].FrameOfReferenceUID}'
                    FORUID1 = f'{file.RegistrationSequence[1].FrameOfReferenceUID}'
                
                #if FORUID0 == trgFORUID and FORUID1 == srcFORUID and Type == regTxName\
                #and mod == 'REG': # 07/05/21
                if FORUID0 == trgFORUID and FORUID1 == srcFORUID and matchOnDroType:
                    inds.append(i)
                    
                    if 'Deformable' in droType:
                        gridDimss.append(
                            file.DeformableRegistrationSequence[1]\
                                .DeformableRegistrationGridSequence[0]\
                                .GridDimensions)
                        
                        gridRess.append(
                            file.DeformableRegistrationSequence[1]\
                                .DeformableRegistrationGridSequence[0]\
                                .GridResolution)
                        
                        vectGridDatas.append(
                            file.DeformableRegistrationSequence[1]\
                                .DeformableRegistrationGridSequence[0]\
                                .VectorGridData)
                    else:
                        txMatrices.append(
                            file.RegistrationSequence[1]\
                                .MatrixRegistrationSequence[0]\
                                .MatrixSequence[0]\
                                .FrameOfReferenceTransformationMatrix)
                
                    """
                    txMatrixTypes.append(
                        file.RegistrationSequence[1]\
                            .MatrixRegistrationSequence[0]\
                            .MatrixSequence[0]\
                            .FrameOfReferenceTransformationMatrixType)
                    """
        
                    contentDateTime = f'{file.ContentDate}{file.ContentTime}'
                    
                    if p2c:
                        print(f'  contentDateTime = {contentDateTime}')
                    
                    # Convert to datetime object:
                    contentDateTime = datetime.strptime(contentDateTime, 
                                                        '%Y%m%d%H%M%S.%f')
                    contentDateTimes.append(contentDateTime)
        
        if p2c:
            print(f'  The indices of the files that match = {inds}')
            print(f'  contentDateTimes = {contentDateTimes}')
        
        if contentDateTimes:#time_diffs:
            # The index of the most recent contentDateTime:
            newInd = contentDateTimes.index(max(contentDateTimes))
            
            # The index of the most recent file that matches:
            droInd = inds[newInd]
            
            if p2c:
                print(f'  The index of the most recent file that matches = {droInd}')
            
            dro = files[droInd]
            
            if regTxName == 'bspline':
                txMatrix = None
                gridDims = gridDimss[newInd]
                gridRes = gridRess[newInd]
                vectGridData = vectGridDatas[newInd]
                
                gridResRounded = [round(item, 2) for item in gridRes]
                
                msg = f"  There were {len(fnames)} subject assessors found, of "\
                      + f"which {len(contentDateTimes)} was a {droType}, with "\
                      + "\n  FrameOfReferenceUIDs matching those of the Source "\
                      + "and Target image series, matching the transform type, "\
                      + "\n  the most recent of which contains the grid "\
                      + f"dimensions {gridDims}, grid resolution {gridResRounded},"\
                      + f" and vector grid data containing {len(vectGridData)} "\
                      + "elements.\n"
            
            if regTxName in ['rigid', 'affine']:
                txMatrix = txMatrices[newInd]
                gridDims = None
                gridRes = None
                vectGridData = None
                
                txMatrixRounded = [round(item, 3) for item in txMatrix]
                
                msg = f"  There were {len(fnames)} subject assessors found, of "\
                      + f"which {len(contentDateTimes)} was a {droType}, with "\
                      + "\n  FrameOfReferenceUIDs matching those of the Source and"\
                      + " Target image series, matching the transform type, "\
                      + "\n  the most recent of which contains the transformation"\
                      + f" matrix \n  {txMatrixRounded}\n"
        else:
            dro = None
            txMatrix = None
            gridDims = None
            gridRes = None
            vectGridData = None
            
            msg = "  There were no DROs with FrameOfReferenceUIDs matching those "\
                  + "of the Source ({srcFORUID}) \nand Target ({trgFORUID}) "\
                  + "image series, and whose SOP Class matched the desired "\
                  + f"transform type ({regTxName}).\n"
        print(msg)
        
        self.dro = dro
        self.txMatrix = txMatrix
        self.gridDims = gridDims
        self.gridRes = gridRes
        self.vectGridData = vectGridData