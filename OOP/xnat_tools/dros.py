# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:58:34 2021

@author: ctorti
"""


def search_dro(xnatSession, projID, subjLab, srcExpLab, srcScanID, trgExpLab, 
               trgScanID, regTxName, p2c=False):
    """
    Search XNAT for a DRO (as a subject assessor) whose FrameOfReferenceUIDs   
    matches the Source (moving) and Target (fixed) FrameOfReferenceUIDs, and 
    whose SOP Class UID matches the transformation type specified by regTxName. 
    
    If there is one (or more) match(es), return the (newest) DRO, the
    FrameOfReferenceTransformation (for non-deformable transforms), or the 
    GridDimensions, GridResolution and VectorGridData (for deformable 
    transforms) of the moving (Source) sequence.
    
    Parameters
    ----------
    xnatSession : Requests session
        A Requests session for XNAT.
    projID : str
        The project ID of interest.
    subjLab : str
        The subject label of interest. 
    srcExpLab : str
        The Source DICOM study / XNAT experiment label of interest. 
    srcScanID : str
        The Source DICOM series label / XNAT scan ID of interest.
    trgExpLab : str
        The Target DICOM study / XNAT experiment label of interest. 
    trgScanID : str
        The Target DICOM series label / XNAT scan ID of interest.
    regTxName : str
        Denotes type of registration transformation. Acceptable values include:
        - 'rigid'
        - 'affine'
        - 'bspline' (i.e. deformable)
    p2c : bool, optional
        Denotes whether some results will be logged to the console.
    
    Returns
    -------
    dro : Pydicom Object or None
        A DICOM Registration Object (if there exists a DRO, as a subject 
        assessor) that matches the requirements); None if not.
    txMatrix : list of float str or None
        List of float strings representing the non-deformable transformation 
        that transforms the moving (Source) image to the fixed (Target) image 
        (if there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    gridDims : list of int or None
        List of integers representing the dimensions of the BSpline grid used 
        to deform the moving (Source) image to the fixed (Target) image (if 
        there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    gridRes : list of float or None
        List of floats representing the resolution of the BSpline grid used to 
        deform the moving (Source) image to the fixed (Target) image (if there 
        exists a DRO (as a subject assessor) that matches the requirements); 
        None if not.
    vectGridData : list of float str or None
        List of floats representing the vector deformations that deform the
        moving (Source) image to the fixed (Target) image (if there exists a
        DRO (as a subject assessor) that matches the requirements); None if not.
    """
    
    from pydicom.filebase import DicomBytesIO
    from pydicom import dcmread
    #import time
    from datetime import datetime
    from xnat_tools.FOR_uid import get_FOR_uid
    
    if p2c:
        print('\n\n', '-'*120)
        print('Running of search_dro():')
        print('\n\n', '-'*120)
        
    #print(xnatSession)
    
    url = xnatSession.url
    
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
                
    # Get a list of the FrameOfReferenceUID as a tuple of pairs, 
    # e.g. (SrcForUid, TrgForUid), for all items:
    #ForUidPairs = []
    
    # Initialise lists:
    txMatrices = []
    #txMatrix_types = []
    gridDimss = [] # note extra s since a list of gridDims
    gridRess = [] # note extra s since a list of gridRes
    vectGridDatas = []
    
    #contentDateTimes = []
    
    # The indices of any DROs whose FORuids match the Source and Target FORuids
    # and whose SOPClassUID is compatible with regTxName: 
    inds = []
    
    # Find the time difference between now and each contentDateTime:
    #now = time.time()
    #now = datetime.now()
    #now = datetime.strftime(Now, '%Y%m%d%H%M%S.%f')
    
    #if p2c:
    #    print(f'  now = {now}')

    #time_diffs = []
    
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
        
        #print('file =', Dro, '\n')
        #print(f'file.SOPClassUID = {file.SOPClassUID}\n')
        #print(f'file[0x0008, 0x0016].name = {file[0x0008, 0x0016].name}\n')
        #print(f'file[0x0008, 0x0016].value = {file[0x0008, 0x0016].value}\n')
        
        mod = f'{file.Modality}'
        
        # Only proceed if the DICOM file is a registration object:
        if mod == 'REG': # 07/05/21
            #droType = f'{file.SOPClassUID}' # this is a number (e.g. 1.2.840.10008.5.1.4.1.1.66.3)
            droType = f'{file[0x0008, 0x0016].repval}' # the representation of the element's value
            
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
            
            #ForUidPairs.append((Dro.RegistrationSequence[0].FrameOfReferenceUID,
            #                    Dro.RegistrationSequence[1].FrameOfReferenceUID))
                       
            # Use f-strings instead of deepcopy:
            if 'Deformable' in droType:
                FORUID0 = f'{file.DeformableRegistrationSequence[0].SourceFrameOfReferenceUID}'
                FORUID1 = f'{file.DeformableRegistrationSequence[1].SourceFrameOfReferenceUID}'
            else:
                FORUID0 = f'{file.RegistrationSequence[0].FrameOfReferenceUID}'
                FORUID1 = f'{file.RegistrationSequence[1].FrameOfReferenceUID}'
            
            #Type = f'{file.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrixType}'
            #Type = Type.lower()
            
            #mod = f'{file.Modality}' # 07/05/21
                        
            #print(Type)
            
            #if FORUID0 == trgFORUID and FORUID1 == srcFORUID and Type == regTxName\
            #and mod == 'REG': # 07/05/21
            if FORUID0 == trgFORUID and FORUID1 == srcFORUID and matchOnDroType:
                inds.append(i)
                
                if 'Deformable' in droType:
                    gridDimss.append(file.DeformableRegistrationSequence[1]\
                                          .DeformableRegistrationGridSequence[0]\
                                          .GridDimensions)
                    
                    gridRess.append(file.DeformableRegistrationSequence[1]\
                                         .DeformableRegistrationGridSequence[0]\
                                         .GridResolution)
                    
                    vectGridDatas.append(file.DeformableRegistrationSequence[1]\
                                               .DeformableRegistrationGridSequence[0]\
                                               .VectorGridData)
                else:
                    txMatrices.append(file.RegistrationSequence[1]\
                                           .MatrixRegistrationSequence[0]\
                                           .MatrixSequence[0]\
                                           .FrameOfReferenceTransformationMatrix)
            
                #txMatrix_types.append(file.RegistrationSequence[1]\
                #                           .MatrixRegistrationSequence[0]\
                #                           .MatrixSequence[0]\
                #                           .FrameOfReferenceTransformationMatrixType)
    
                contentDateTime = f'{file.ContentDate}{file.ContentTime}'
                
                if p2c:
                    print(f'  contentDateTime = {contentDateTime}')
                
                # Convert to datetime object:
                contentDateTime = datetime.strptime(contentDateTime, 
                                                    '%Y%m%d%H%M%S.%f')
                contentDateTimes.append(contentDateTime)
                
                #print(f'contentDateTime = {contentDateTime}')
                #print(f'now = {now}')
                
                #t0 = time.mktime(time.strptime(contentDateTime, 
                #                               '%Y%m%d%H%M%S.%f'))
                #t0 = datetime.strptime(contentDateTime, '%Y%m%d%H%M%S.%f')
                
                #if p2c:
                #    print(f'  t0 = {t0}')
                
                #diff = now - t0
                
                #print(f't0 = {t0}')
                #print(f'diff = {diff}\n')
    
                #time_diffs.append(diff)
    
    if p2c:
        print(f'  The indices of the files that match = {inds}')
        #print(f'  time_diffs = {time_diffs}')
        print(f'  contentDateTimes = {contentDateTimes}')
    
    if contentDateTimes:#time_diffs:
        # The index of the most recent contentDateTime:
        #newInd = time_diffs.index(min(time_diffs))
        newInd = contentDateTimes.index(max(contentDateTimes))
        
        # The index of the most recent file that matches:
        droInd = inds[newInd]
        
        if p2c:
            print(f'  The index of the most recent file that matches = {droInd}')
        
        dro = files[droInd]
        
        #if 'Deformable' in droType:
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
        #else:
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
        msg = "  There were no DROs with FrameOfReferenceUIDs matching those "\
              + "of the Source ({srcFORUID}) \nand Target ({trgFORUID}) "\
              + "image series, and whose SOP Class matched the desired "\
              + f"transform type ({regTxName}).\n"
        
        dro = None
        txMatrix = None
        gridDims = None
        gridRes = None
        vectGridData = None
    
    print(msg)
    
    if p2c:
        print('-'*120)
    
    return dro, txMatrix, gridDims, gridRes, vectGridData