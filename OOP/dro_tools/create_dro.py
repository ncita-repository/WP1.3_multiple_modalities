# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 19:07:06 2021

@author: ctorti
"""

from importlib import reload
import dro_tools.create_spa_dro
import dro_tools.create_def_dro
import dro_tools.matrices
import general_tools.general
reload(general_tools.general)
reload(dro_tools.create_spa_dro)
reload(dro_tools.create_def_dro)
reload(dro_tools.matrices)
import xnat_tools.subject_assessors
reload(xnat_tools.subject_assessors)

import os
from pathlib import Path
#import time
from datetime import datetime
from pydicom import dcmread
from pydicom.uid import generate_uid
import SimpleITK as sitk
import numpy as np

from dro_tools.general import (
    modify_ref_ser_seq, modify_stu_con_oth_ref_ins_seq, modify_ref_im_seq
    )
#from dro_tools.matrices import is_matrix_orthonormal, is_matrix_orthogonal
from dro_tools.matrices import (get_txMatrix_from_tx, get_tx_matrix_type)
from general_tools.general import (
    reduce_list_of_str_floats_to_16, generate_reg_fname
    )
from xnat_tools.subject_assessors import upload_subj_asr

#from dro_tools.create_spa_dro import create_spa_dro_from_tx
#from dro_tools.create_def_dro import create_def_dro_from_tx




class DroCreator:
    # TODO update docstrings
    """
    This class creates a DroCreator Object to create a DICOM Registration 
    Object (DRO).
    
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
    self.sampleDro : Pydicom Object
        The DICOM Registration Object that will be used as a template for
        the new DRO.
    self.txParams : tuple of floats
        The registration transformation parameters.
    self.txMatrix : list of strs or None
        List of float strings representing the non-deformable transformation 
        that transforms the moving (Source) image to the fixed (Target) image 
        (if there exists a DRO (as a subject assessor) that matches the 
        requirements); None if not.
    self.gridOrig : list of strs, conditional
        The origin of the bspline grid along x, y and z. Only if the 
        registration transform is a bspline.
    self.gridDir : list of strs, conditional
        The direction of the bspline grid. Only if the registration 
        transform is a bspline.
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
    self.dro : Pydicom Object
        The new DICOM Spatial Registration Object.
        
    Note
    ----
    Sample DROs from:
    https://www.insight-journal.org/browse/publication/923
    
    Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html
    
    Example Spatial DRO:
    https://docs.google.com/document/d/1tMUWGe4kw6yLC2j7-y9WC-gY6LoiRL8XbG-_LeRU2_U
    
    Deformable Spatial Registration Object module:
    https://dicom.innolitics.com/ciods/deformable-spatial-registration
    http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.3.html
    
    Example Deformable DROs:
    https://docs.google.com/document/d/1i1C_bJNfWkqfxyCctzSxCd_y4lCAwiiPxSWa3QKkrFQ
    https://docs.google.com/document/d/1MDTQfAr3Ws-TDPZi_Bj58X2YwBpOvEv4mb_12U_NXTw
    
    VectorGridData has VR = OF:
    
    A stream of 32-bit IEEE 754:1985 floating point words. OF is a VR that 
    requires byte swapping within each 32-bit word when changing byte ordering 
    (see Section 7.3).
    
    http://dicom.nema.org/medical/dicom/current/output/chtml/part05/sect_6.2.html
    """
    
    
    def __init__(self, newDataset, params):
        # TODO update docstrings
        """
        Initialise the object.
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.sampleDro : Pydicom Object
            The DICOM Registration Object that will be used as a template for
            the new DRO.
        #self.resTxParams : tuple of floats, conditional
        #    The registration transformation parameters if the transform is a
        #    rigid transform; None otherwise.
        self.resTxMatrix : list of floats, conditional
            List form of the transformation matrix if the transform is a
            rigid transform; None otherwise.
        self.preRegTxMatrix : list of floats, conditional
            List form of the transformation matrix if the transform is a
            bspline transform; None otherwise.
        self.gridOrig : list of strs, conditional
            The origin of the bspline grid along x, y and z if the transform is
            a bspline transform; None otherwise.
        self.gridDir : list of strs, conditional
            The direction of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.gridDims : list of floats, conditional
            The dimensions of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.gridRes : list of floats, conditional
            The spacings of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.vectGridData : list of floats, conditional
            The deformations of the bspline if the transform is a bspline
            transform, e.g.: 
            [def0_x, def0_y, def0_z, def1_x, def1_y, def1_z, ...];
            None otherwise.
        """
        
        """
        Decide whether to reference all SOPInstanceUIDs in srcDicoms and 
        trgDicoms within ReferencedInstanceSequence (not required), or just one 
        (for compactness) to record the SeriesInstanceUIDs of the two series.
        """
        #self.refAllSOPs = True
        self.refAllSOPs = False
        
        """ 
        Consider adding this as an input variable for the user to specify.
        Empty for now.
        """
        # ContentDescription:
        #self.contentDesc = ''
        
        # Initialise the transform data that will be updated for rigid 
        # transforms:
        self.resTxMatrix = None
        # and bspline transforms:
        self.preRegTxMatrix = None
        self.gridOrig = None
        self.gridDir = None
        self.gridDims = None
        self.gridRes = None
        self.vectGridData = None
        
        
        # Get the transformation matrix and grid data (if the registration
        # transformation is a bspline):
        #self.get_txMatrix_from_tx(newDataset, params)
        self.get_data_from_tx(newDataset, params)
        
        # Import the sample DRO:
        self.import_sample_dro(params)
    
    def import_sample_dro(self, params):
        # TODO update docstrings
        """
        Import the sample DRO appropriate for the registration transform.  
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        description : str, optional ('' by default)
            Description to be added to ContentDescription, and used in the file
            name of the exported DRO.
        
        Returns
        -------
        self.sampleDro : Pydicom Object
            The DICOM Registration Object that will be used as a template for
            the new DRO.
        """
        
        regTxName = params.cfgDict['regTxName']
        sampleDroDir = params.cfgDict['sampleDroDir']
        
        if regTxName == 'bspline':
            sampleDroFname = 'sample_deformable_dro.dcm'
        else:
            sampleDroFname = 'sample_spatial_dro.dcm'
            
        sampleDroFpath = os.path.join(sampleDroDir, sampleDroFname)
        # Read in the DRO template:
        sampleDro = dcmread(sampleDroFpath)
        
        self.sampleDro = sampleDro
    
    def get_data_from_rigid_tx(self, newDataset, params):
        # TODO update docstrings
        """
        Get the transform parameters and transform matrix.  
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.resTxMatrix : list of floats
            List form of the transformation matrix.
        """
        
        #regTxName = params.cfgDict['regTxName']
        p2c = params.cfgDict['p2c']
        
        # The registration transform:
        #regTx = newDataset.finalTx # 04/09/21
        regTx = newDataset.resTx # 04/09/21
        
        self.resTxMatrix = get_txMatrix_from_tx(regTx, p2c)
    
    def get_data_from_bspline_tx(self, newDataset, params):
        # TODO update docstrings
        """
        Get the registration transform matrix as a list, and grid parameters if
        the registration transform is non-rigid.  
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        self.preRegTxMatrix : list of floats
            List form of the pre-registration transformation matrix.
        self.gridOrig : list of strs
            The origin of the bspline grid along x, y and z.
        self.gridDir : list of strs
            The direction of the bspline grid.
        self.gridDims : list of floats
            The dimensions of the bspline grid.
        self.gridRes : list of floats
            The spacings of the bspline grid.
        self.vectGridData : list of floats
            The deformations of the bspline, e.g.: 
            [def0_x, def0_y, def0_z, def1_x, def1_y, def1_z, ...].
        """
        
        #regTxName = params.cfgDict['regTxName']
        p2c = params.cfgDict['p2c']
        
        # The registration and pre-registration transforms:
        #regTx = newDataset.finalTx # 04/09/21
        regTx = newDataset.resTx # 04/09/21
        preRegTx = newDataset.preRegTx
        
        # The transform parameters:
        #txParams = regTx.GetParameters()
        
        """ 
        Could also get the tuple of the registration transform parameters
        directly from?:
        txParams = newDataset.txParams
        """
        # TODO tidy above
        
        """ 
        The final transform after executing the registration is of class
        SimpleITK.SimpleITK.Transform but needs to be of class
        SimpleITK.SimpleITK.BSplineTransform. Check and change if needed:
        """
        if type(regTx) == sitk.SimpleITK.Transform:
            regTx = sitk.BSplineTransform(regTx)
        
        # Get the coefficient images (deformations along x, y and z):
        coeffImX, coeffImY, coeffImZ = regTx.GetCoefficientImages()
        
        # Zip the vectors to a list of x, y, z tuples:
        listOfTuples = list(zip(
            list(sitk.GetArrayViewFromImage(coeffImX).flatten()),
            list(sitk.GetArrayViewFromImage(coeffImY).flatten()),
            list(sitk.GetArrayViewFromImage(coeffImZ).flatten())
            ))
        
        flatListOfVects = [item for tup in listOfTuples for item in tup]
        
        vectGridData = np.array(flatListOfVects).tobytes()
        
        gridOrig = [str(item) for item in coeffImX.GetOrigin()]
        gridDir = [str(item) for item in coeffImX.GetDirection()]
        gridDims = list(coeffImX.GetSize())
        gridRes = list(coeffImX.GetSpacing())
        
        if len(gridDir) > 6:
            # Ignore the direction cosines along the z-direction:
            gridDir = gridDir[0:6]
        
        # Get the transformation matrix (txMatrix) from the pre-registration
        # transform parameters: """
        if preRegTx == None:
            #preRegTxParams = None
            
            txMatrix = ['1.0', '0.0', '0.0', '0.0',
                        '0.0', '1.0', '0.0', '0.0',
                        '0.0', '0.0', '1.0', '0.0',
                        '0.0', '0.0', '0.0', '1.0']
        else:
            """
            txParams = preRegTx.GetParameters()
            
            # Add the row vector [0, 0, 0, 1] to txParams to complete the Frame of
            # Reference Transformation Matrix 
            # (http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.20.2.html#sect_C.20.2.1.1). 
            #txParams = [str(item) for item in txParams] + ['0.0', '0.0', '0.0', '1.0'] # 06/05/21
            txMatrix = [
                str(txParams[0]), str(txParams[1]), str(txParams[2]), '0.0',
                str(txParams[3]), str(txParams[4]), str(txParams[5]), '0.0',
                str(txParams[6]), str(txParams[7]), str(txParams[8]), '0.0',
                str(txParams[9]), str(txParams[10]), str(txParams[11]), '1.0'
                ]
            """
            
            txMatrix = get_txMatrix_from_tx(preRegTx, p2c)
        
        if p2c:
            print(f'type(gridOrig) = {type(gridOrig)}')
            print(f'gridOrig = {gridOrig}\n')
            
            print(f'type(gridDir) = {type(gridDir)}')
            print(f'gridDir = {gridDir}\n')
            
            print(f'type(gridDims) = {type(gridDims)}')
            print(f'gridDims = {gridDims}\n')
            
            print(f'type(gridRes) = {type(gridRes)}')
            print(f'gridRes = {gridRes}\n')
            
            print(f'type(vectGridData) = {type(vectGridData)}\n')
            #print(f'vectGridData = {vectGridData}\n')
            
            #print(f'txParams = {txParams}\n')
            print(f'preRegTxMatrix = {txMatrix}\n')
        
        self.gridOrig = gridOrig
        self.gridDir = gridDir
        self.gridDims = gridDims
        self.gridRes = gridRes
        self.vectGridData = vectGridData
        self.preRegTxMatrix = txMatrix  # 07/09/21
    
    def get_data_from_tx(self, newDataset, params):
        # TODO update docstrings
        """
        Get the data from the transform that will be required for creating a
        DRO.  
        
        Parameters
        ----------
        newDataset : Propagator Object
            Object containing various data for the source-to-target propagated (or 
            copied) dataset.
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        #self.resTxParams : tuple of floats, conditional
        #    The registration transformation parameters if the transform is a
        #    rigid transform; None otherwise.
        self.resTxMatrix : list of floats, conditional
            List form of the transformation matrix if the transform is a
            rigid transform; None otherwise.
        self.preRegTxMatrix : list of floats, conditional
            List form of the transformation matrix if the transform is a
            bspline transform; None otherwise.
        self.gridOrig : list of strs, conditional
            The origin of the bspline grid along x, y and z if the transform is
            a bspline transform; None otherwise.
        self.gridDir : list of strs, conditional
            The direction of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.gridDims : list of floats, conditional
            The dimensions of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.gridRes : list of floats, conditional
            The spacings of the bspline grid if the transform is a bspline
            transform; None otherwise.
        self.vectGridData : list of floats, conditional
            The deformations of the bspline if the transform is a bspline
            transform, e.g.: 
            [def0_x, def0_y, def0_z, def1_x, def1_y, def1_z, ...];
            None otherwise.
        """
        
        regTxName = params.cfgDict['regTxName']
        
        if regTxName == 'bspline':
            self.get_data_from_bspline_tx(newDataset, params)
        else:
            self.get_data_from_rigid_tx(newDataset, params)
    
    def create_spa_dro(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        """
        Create a DICOM Spatial Registration Object (DRO).  
        
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
        description : str, optional ('' by default)
            Description to be added to ContentDescription, and used in the file
            name of the exported DRO.
        
        Returns
        -------
        self.dro : Pydicom Object
            The new DICOM Spatial Registration Object.
        """
        
        # The tuple of the registration transform parameters, and the list of
        # parameters that define the transform matrix:
        #txParams = self.txParams # 04/09/21
        #txMatrix = self.txMatrix # 04/09/21
        #resTxParams = self.resTxParams # 04/09/21
        resTxParams = newDataset.resTxParams # 07/09/21
        resTxMatrix = self.resTxMatrix # 04/09/21
        
        refAllSOPs = self.refAllSOPs
        #contentDesc = self.contentDesc
        dro = self.sampleDro
        p2c = params.cfgDict['p2c']
        
        srcDicoms = srcDataset.dicoms
        trgDicoms = trgDataset.dicoms
        
        srcExpLab = params.cfgDict['srcExpLab']
        srcScanID = params.cfgDict['srcScanID']
        trgExpLab = params.cfgDict['trgExpLab']
        trgScanID = params.cfgDict['trgScanID']
        #regTxName = params.cfgDict['regTxName']
        
        if p2c:
            if refAllSOPs:
                print('Note: All SOPInstanceUIDs of both series will be referenced',
                      'in ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
            else:
                print('Note: Only the first SOPInstanceUID of both series will be',
                      'referenced in',
                      'ReferencedSeriesSequence[i].ReferencedInstanceSequence.')
        
        # Start timing:
        #times = []
        #times.append(time.time())
        
        #timingMsg = "Creating the DICOM Spatial Registration Object...\n"
        #params.add_timestamp(timingMsg)
        
        #currentDate = time.strftime("%Y%m%d", time.gmtime())
        #currentTime = time.strftime("%H%M%S", time.gmtime())
        ##currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        timeNow = datetime.now()
        currentDate = timeNow.strftime('%Y%m%d')
        currentTime = timeNow.strftime('%H%M%S.%f')
        
        dro.InstanceCreationDate = currentDate
        dro.InstanceCreationTime = currentTime
        
        # Generate a new SOPInstanceUID:
        dro.SOPInstanceUID = generate_uid()
        dro.file_meta.MediaStorageSOPInstanceUID = dro.SOPInstanceUID
        
        dro.StudyDate = trgDicoms[0].StudyDate
        
        try:
            dro.SeriesDate = trgDicoms[0].SeriesDate
        except AttributeError:
            pass
        
        dro.StudyTime = trgDicoms[0].StudyTime
        
        try:
            dro.SeriesTime = trgDicoms[0].SeriesTime
        except AttributeError:
            pass
        
        #dro.ContentDate = trgDicoms[0].ContentDate
        #dro.ContentTime = trgDicoms[0].ContentTime
        dro.ContentDate = currentDate
        dro.ContentTime = currentTime
        dro.ContentDescription = f'ExpLab {srcExpLab} ScanID {srcScanID} '\
            + f'registered to ExpLab {trgExpLab} ScanID {trgScanID}'
        
        #dro.Manufacturer = trgDicoms[0].Manufacturer
        dro.Manufacturer = 'NCITA'
        """ Consider modifying StudyDescription (= '' in the template DRO. """
        try:
            dro.SeriesDescription = trgDicoms[0].SeriesDescription
        except AttributeError:
            #pass
            dro.SeriesDescription = ''
        #dro.ManufacturerModelName = trgDicoms[0].ManufacturerModelName
        dro.ManufacturerModelName = ''
        dro.PatientName = trgDicoms[0].PatientName
        dro.PatientID = trgDicoms[0].PatientID
        dro.PatientBirthDate = trgDicoms[0].PatientBirthDate
        dro.PatientSex = trgDicoms[0].PatientSex
        try:
            dro.PatientAge = trgDicoms[0].PatientAge
        except AttributeError:
            pass
        dro.SoftwareVersions = ''
        dro.StudyInstanceUID = trgDicoms[0].StudyInstanceUID
        # Generate a new UID for SeriesInstanceUID:
        dro.SeriesInstanceUID = generate_uid()
        dro.StudyID = trgDicoms[0].StudyID
        dro.SeriesNumber = trgDicoms[0].SeriesNumber
        dro.InstanceNumber = trgDicoms[0].InstanceNumber
        dro.FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
        dro.PositionReferenceIndicator = trgDicoms[0].PositionReferenceIndicator
        """ Keep ContentLabel as 'REGISTRATION'. """
        """ Keep ContentCreatorName as ''. """
        
        # Modify the RegistrationSequence for the fixed domain:
        dro.RegistrationSequence[0]\
           .FrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
        
        """
        Since the first sequence in RegistrationSequence will be the fixed image
        domain, FrameOfReferenceTransformationMatrix will be the identity matrix
        ['1.0', '0.0', '0.0', '0.0', 
         '0.0', '1.0', '0.0', '0.0', 
         '0.0', '0.0', '1.0', '0.0', 
         '0.0', '0.0', '0.0', '1.0'].
        Hence no change is required FrameOfReferenceTransformationMatrix for this
        sequence.
        """
        
        """
        Since the FrameOfReferenceTransformationMatrix is the identity matrix, the 
        FrameOfReferenceTransformationMatrixType will be 'RIGID'.  Hence no change
        is required to FrameOfReferenceTransformationMatrixType for this sequence.
        """
        
        """ 
        Consider adding FrameOfReferenceTransformationComment:
        dro.RegistrationSequence[0]\
           .MatrixRegistrationSequence[0]\
           .FrameOfReferenceTransformationComment = 
        """
        
        # Modify the RegistrationSequence for the moving domain:
        dro.RegistrationSequence[1]\
           .FrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
        
        # Modify the FrameOfReferenceTransformationMatrix:
        dro.RegistrationSequence[1]\
           .MatrixRegistrationSequence[0]\
           .MatrixSequence[0]\
           .FrameOfReferenceTransformationMatrix = resTxMatrix # was txParams (06/05/21)
           
        # Modify the FrameOfReferenceTransformationMatrixType.  
        """
        Acceptable values:
         'RIGID' (= value in template DRO)
         'RIGID_SCALE' (Similarity transform)
         'AFFINE' 
         Note:  Need to verify that get_tx_matrix_type covers all 3 options.
        """
        dro.RegistrationSequence[1]\
           .MatrixRegistrationSequence[0]\
           .MatrixSequence[0]\
           .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(
               resTxParams, p2c
               )
        
        #times.append(time.time())
        #dtime = round(times[-1] - times[-2], 3)
        #if p2c:
        #    print(f'\n   *Took {dtime} s to get the transform matrix type.')
            
        """ 
        Consider adding FrameOfReferenceTransformationComment:
        dro.RegistrationSequence[1]\
           .MatrixRegistrationSequence[0]\
           .FrameOfReferenceTransformationComment = ''
        """
         
        # Modify the ReferencedSeriesSequence:
        """
        Note:
            ReferencedSeriesSequence will have two items: The first for the fixed 
            series, and second for the moving series.
        """
        # Modify ReferencedSeriesSequence for the fixed series:
        dro = modify_ref_ser_seq(
            dro, seqNum=0, dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
        
        # Modify ReferencedSeriesSequence for the moving series:
        dro = modify_ref_ser_seq(
            dro, seqNum=1, dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
         
        # Modify the StudiesContainingOtherReferencedInstancesSequence:
        """
        Note:
            StudiesContainingOtherReferencedInstancesSequence will have two items: 
            The first for the fixed series, and second for the moving series.
        """
        # Modify StudiesContainingOtherReferencedInstancesSequence for the fixed 
        # series:
        dro = modify_stu_con_oth_ref_ins_seq(
            dro, seqNum=0, dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
        
        # Modify StudiesContainingOtherReferencedInstancesSequence for the moving 
        # series:
        dro = modify_stu_con_oth_ref_ins_seq(
            dro, seqNum=1, dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
        
        """ > Consider adding UsedFiducialsSequence (optional). 
        http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
        
        Spatial Fiducials Series Module:
        http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
        """
        
        
        """
        https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
        """
        
        if p2c:
            print('Deleting GroupLength tags.')
        
        """
        The SRO template contains 4 retired tags with name "GroupLength". 
        Delete them.
        """
        del dro[0x08, 0x00]
        del dro[0x10, 0x00]
        del dro[0x20, 0x00]
        del dro[0x70, 0x00]
        
        #times.append(time.time())
        #dtime = round(times[-1] - times[-2], 1)
        #if p2c:
        #    print(f'\n   *Took {dtime} s to delete all GroupLength tags.')
            
        #dtime = round(times[-1] - times[0], 1)
        #if p2c:
        #    print(f'\n   *Took a total of {dtime} s to create the DRO.')
        
        #timingMsg = "Took [*] s to create the DICOM Spatial Registration "\
        #    + "Object.\n"
        #params.add_timestamp(timingMsg)
        
        self.dro = dro

    def create_def_dro(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        """
        Create a DICOM Deformable Spatial Registration Object (DRO).  
        
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
        description : str, optional ('' by default)
            Description to be added to ContentDescription, and used in the file
            name of the exported DRO.
        
        Returns
        -------
        self.dro : Pydicom object
            The DICOM Deformable Spatial Registration Object.
        """
        
        # The tuple of the registration transform parameters, the list of
        # parameters that define the transform matrix, and the bspline grid
        # parameters:
        #txParams = self.txParams # 04/09/21
        #txMatrix = self.txMatrix # 04/09/21
        #resTxParams = self.resTxParams # 04/09/21
        #resTxParams = newDataset.resTxParams # 07/09/21
        #resTxMatrix = self.resTxMatrix # 04/09/21
        preRegTxMatrix = self.preRegTxMatrix # 07/09/21
        gridOrig = self.gridOrig
        gridDir = self.gridDir
        gridDims = self.gridDims
        gridRes = self.gridRes
        vectGridData = self.vectGridData
        
        refAllSOPs = self.refAllSOPs
        #contentDesc = self.contentDesc
        dro = self.sampleDro
        p2c = params.cfgDict['p2c']
        
        srcDicoms = srcDataset.dicoms
        trgDicoms = trgDataset.dicoms
        
        srcExpLab = params.cfgDict['srcExpLab']
        srcScanID = params.cfgDict['srcScanID']
        trgExpLab = params.cfgDict['trgExpLab']
        trgScanID = params.cfgDict['trgScanID']
        
        if p2c:
            if refAllSOPs:
                print('Note: All SOPInstanceUIDs of both series will be',
                      'referenced in',
                      'DeformableRegistrationSequence[0].ReferencedImageSequence.')
            else:
                print('Note: Only the first SOPInstanceUID of both series will',
                      'be referenced in',
                      'DeformableRegistrationSequence[0].ReferencedImageSequence.')
        
        # Start timing:
        #times = []
        #times.append(time.time())
        
        #timingMsg = "Creating the DICOM Deformable Spatial Registration "\
        #    + "Object...\n"
        #params.add_timestamp(timingMsg)
        
        #currentDate = time.strftime("%Y%m%d", time.gmtime())
        #currentTime = time.strftime("%H%M%S", time.gmtime())
        ##currentDateTime = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        
        timeNow = datetime.now()
        currentDate = timeNow.strftime('%Y%m%d')
        currentTime = timeNow.strftime('%H%M%S.%f')
        
        dro.InstanceCreationDate = currentDate
        dro.InstanceCreationTime = currentTime
        
        # Generate a new SOPInstanceUID:
        dro.SOPInstanceUID = generate_uid()
        dro.file_meta.MediaStorageSOPInstanceUID = dro.SOPInstanceUID
        
        dro.StudyDate = srcDicoms[0].StudyDate
        try:
            seriesDate = srcDicoms[0].SeriesDate
            try:
                dro.SeriesDate = seriesDate
            except AttributeError:
                pass
        except AttributeError:
            pass
        
        dro.StudyTime = srcDicoms[0].StudyTime
        
        try:
            seriesTime = srcDicoms[0].SeriesTime
            try:
                dro.SeriesTime = seriesTime
            except AttributeError:
                pass
        except AttributeError:
            pass
        
        #dro.ContentDate = trgDicoms[0].ContentDate
        #dro.ContentTime = trgDicoms[0].ContentTime
        dro.ContentDate = currentDate
        dro.ContentTime = currentTime
        dro.ContentDescription = f'ExpLab {srcExpLab} ScanID {srcScanID} '\
            + f'registered to ExpLab {trgExpLab} ScanID {trgScanID}'
        
        #dro.Manufacturer = trgDicoms[0].Manufacturer
        dro.Manufacturer = 'NCITA'
        
        dro.AccessionNumber = ''
        #dro.Manufacturer = srcDicoms[0].Manufacturer
        dro.Manufacturer = 'NCITA'
        """ Consider modifying StudyDescription (= '' in the template DRO. """
        try:
            studyDesc = srcDicoms[0].StudyDescription
            try:
                dro.StudyDescription = studyDesc
            except AttributeError:
                pass
        except AttributeError:
            pass
        try:
            seriesDesc = srcDicoms[0].SeriesDescription
            try:
                dro.SeriesDescription = seriesDesc
            except AttributeError:
                #pass
                dro.SeriesDescription = ''
        except AttributeError:
            pass
        #dro.ManufacturerModelName = srcDicoms[0].ManufacturerModelName
        dro.ManufacturerModelName = ''
        dro.PatientName = srcDicoms[0].PatientName
        dro.PatientID = srcDicoms[0].PatientID
        dro.PatientBirthDate = srcDicoms[0].PatientBirthDate
        dro.PatientSex = srcDicoms[0].PatientSex
        try:
            patientAge = srcDicoms[0].PatientAge
            try:
                dro.PatientAge = patientAge
            except AttributeError:
                pass
        except AttributeError:
            pass
        """ 
        Keep De-identificationMethod as it is, but consider modifying to 
        represent the actual methods used. 
    
        Delete PrivateCreator and PrivateTagData tags (not in the DICOM 
        standard): 
        """
        del dro[0x0013, 0x0010]
        del dro[0x0013, 0x1010]
        del dro[0x0013, 0x1013]
        dro.SoftwareVersions = ''
        dro.StudyInstanceUID = srcDicoms[0].StudyInstanceUID
        # Generate a new UID for SeriesInstanceUID:
        dro.SeriesInstanceUID = generate_uid()
        dro.StudyID = srcDicoms[0].StudyID
        dro.SeriesNumber = srcDicoms[0].SeriesNumber
        dro.InstanceNumber = srcDicoms[0].InstanceNumber
        dro.FrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
        dro.PositionReferenceIndicator = srcDicoms[0].PositionReferenceIndicator
        """ Keep LongitudinalTemporalInformationModified as 'MODIFIED'. """
        
        """ 
        Modify the ReferencedSeriesSequence. 
        
        Note:
            ReferencedSeriesSequence will reference the moving (source) series.
        """
        #dro = ModifyRefSerSeq(dro, SeqNum=0, Dicoms=srcDicoms, 
        #                      refAllSOPs=refAllSOPs, p2c=p2c)
        
        seq = modify_ref_im_seq(
            seq=dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence,
            dicoms=srcDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
        
        dro.ReferencedSeriesSequence[0].ReferencedInstanceSequence = seq
        
        dro.ReferencedSeriesSequence[0]\
           .SeriesInstanceUID = srcDicoms[0].SeriesInstanceUID
        
        
        """ 
        Modify the StudiesContainingOtherReferencedInstancesSequence. 
        
        Note:
            StudiesContainingOtherReferencedInstancesSequence will reference the
            fixed (target) series.
        """
        #dro = modify_stu_con_oth_ref_ins_seq(dro, seqNum=0, dicoms=trgDicoms, 
        #                               refAllSOPs=refAllSOPs, 
        #                               p2c=p2c)
        
        seq = modify_ref_im_seq(
            seq=dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value,
            dicoms=trgDicoms, refAllSOPs=refAllSOPs, p2c=p2c
            )
        
        dro[0x08, 0x1200][0][0x08, 0x1115][0][0x08, 0x114a].value = seq
        
        dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .ReferencedSeriesSequence[0]\
           .SeriesInstanceUID = trgDicoms[0].SeriesInstanceUID
        
        dro.StudiesContainingOtherReferencedInstancesSequence[0]\
           .StudyInstanceUID = trgDicoms[0].StudyInstanceUID
        
        
        """ Modify the DeformableRegistrationSequence. 
        
        Note:
        The sample DRO contains two sequences in DeformableRegistrationSequence
        - the first relates to the fixed image series, the second to the moving
        image series. 
        
        The DICOM standard seems to suggest that only the series to be
        registered is to be referenced in DeformableRegistrationSequence:
        
        https://dicom.innolitics.com/ciods/deformable-spatial-registration/deformable-spatial-registration/00640002
        
        The first sequence even contains MatrixRegistrationSequence, which is
        not a child of DeformableRegistrationSequence.  It seems that they 
        wanted to reference both moving and fixed series, but it's not clear 
        why since the moving series was already fully referenced in 
        ReferencedSeriesSequence, and the fixed in 
        StudiesContainingOtherReferencedInstanceSequence.
        
        But since I've only found examples of DROs (spatial and deformable)
        that contain two sequences in DeformableRegistrationSequence (or
        RegistrationSequence for spatial DROs) I will maintain this structure
        in case there are best-practices that are not necessarily to the 
        standard.
        """
        if p2c:
            print('\nModify the DeformableRegistrationSequence.')
        
        #""" Remove the first sequence in DeformableRegistrationSequence (which
        #references the fixed image series). """
        #dro.DeformableRegistrationSequence.pop(0)
        
        
        """ 
        Modify the first ReferencedImageSequence. 
        
        Note:
            ReferencedImageSequence[0] will reference the fixed (target) series and
            ReferencedImageSequence[1] will reference the moving (source) series.
        """
        #dro = modify_ref_ser_seq(dro, seqNum=0, Dicoms=trgDicoms, 
        #                      refAllSOPs=refAllSOPs, p2c=p2c)
        seq = modify_ref_im_seq(
            seq=dro.DeformableRegistrationSequence[0].ReferencedImageSequence,
            dicoms=trgDicoms, refAllSOPs=True, p2c=p2c
            )
        
        dro.DeformableRegistrationSequence[0].ReferencedImageSequence = seq
        
        
        """ > Modify the SourceFrameOfReferenceUID.
        Note:
        
        The DICOM standard defines registration from Source to Registered RCS 
        (Referenced Coordinate System).
        
        http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_A.39.html
        """
        dro.DeformableRegistrationSequence[0]\
           .SourceFrameOfReferenceUID = trgDicoms[0].FrameOfReferenceUID
        
        # Delete MatrixRegistrationSequence.
        del dro[0x64, 0x02][0][0x70, 0x0309]
        
        
        """ 
        Modify the second ReferencedImageSequence. 
        
        Note:
            ReferencedImageSequence[0] will reference the fixed (target) series and
            ReferencedImageSequence[1] will reference the moving (source) series.
        """
        #dro = modify_ref_ser_seq(dro, seqNum=1, Dicoms=srcDicoms, 
        #                      refAllSOPs=refAllSOPs, p2c=p2c)
        seq = modify_ref_im_seq(
            seq=dro.DeformableRegistrationSequence[1].ReferencedImageSequence,
            dicoms=srcDicoms, refAllSOPs=True, p2c=p2c
            )
        
        dro.DeformableRegistrationSequence[1].ReferencedImageSequence = seq
        
        
        """ > Modify the SourceFrameOfReferenceUID. """
        dro.DeformableRegistrationSequence[1]\
           .SourceFrameOfReferenceUID = srcDicoms[0].FrameOfReferenceUID
        
        
        """ > Modify the DeformableRegistrationGridSequence. """
        dro.DeformableRegistrationSequence[1]\
           .DeformableRegistrationGridSequence[0]\
           .ImagePositionPatient = reduce_list_of_str_floats_to_16(gridOrig)
        
        if len(gridDir) > 6:
            # Ignore the direction cosines along the z-direction:
            gridDir = gridDir[0:6]
            
        dro.DeformableRegistrationSequence[1]\
           .DeformableRegistrationGridSequence[0]\
           .ImageOrientationPatient = reduce_list_of_str_floats_to_16(gridDir)
        
        dro.DeformableRegistrationSequence[1]\
           .DeformableRegistrationGridSequence[0]\
           .GridDimensions = gridDims
        
        dro.DeformableRegistrationSequence[1]\
           .DeformableRegistrationGridSequence[0]\
           .GridResolution = gridRes
        
        dro.DeformableRegistrationSequence[1]\
           .DeformableRegistrationGridSequence[0]\
           .VectorGridData = vectGridData
        
        
        """ > Modify PreDeformationMatrixRegistrationSequence. """
        dro.DeformableRegistrationSequence[1]\
           .PreDeformationMatrixRegistrationSequence[0]\
           .FrameOfReferenceTransformationMatrix = preRegTxMatrix
        
        """ > Modify the FrameOfReferenceTransformationMatrixType.  Acceptable 
        values:
        'RIGID' (= value in template DRO)
        'RIGID_SCALE' (Similarity transform)
        'AFFINE' 
        
        Note:  Need to verify that get_tx_matrix_type covers all 3 options.
        """
        
        #if txParams != None: # 04/06/21
        #    print(f'txParams = {txParams}')
        #    
        #    dro.DeformableRegistrationSequence[1]\
        #       .PreDeformationMatrixRegistrationSequence[0]\
        #       .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(txParams, 
        #                                                                   p2c)
        #
        #    times.append(time.time())
        #    dtime = round(times[-1] - times[-2], 3)
        #    if p2c:
        #        print(f'\n   *Took {dtime} s to get the transform matrix type.')
        
        # 04/06/21:
        dro.DeformableRegistrationSequence[1]\
           .PreDeformationMatrixRegistrationSequence[0]\
           .FrameOfReferenceTransformationMatrixType = get_tx_matrix_type(
               preRegTxMatrix, p2c
               )
               
        """ > Consider adding details for the FrameOfReferenceTransformationComment
        and RegistrationTypeCodeSequence (optional). 
        dro.DeformableRegistrationSequence[1]\
           .FrameOfReferenceTransformationComment = ''
        
        dro.DeformableRegistrationSequence[1]\
           .RegistrationTypeCodeSequence = ''
        """
        
        
        """ > Leave PostDeformationMatrixRegistrationSequence unchanged. """
        
        
        """ Keep ContentLabel as 'REGISTRATION'. """
        dro.ContentCreatorName = 'NCITA'
        # Delete PrivateCreator tag (not in the DICOM standard).
        del dro[0x3773, 0x01]
        
        
        """ > Consider adding UsedFiducialsSequence (optional). 
        
        http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.20.html#sect_C.20.3
        
        Spatial Fiducials Series Module:
        http://dicom.nema.org/dicom/2013/output/chtml/part03/sect_C.21.html
        """
        
        """
        https://dicom.innolitics.com/ciods/spatial-registration/common-instance-reference/00081115/0008114a/00081155
        """
        
        #dtime = round(times[-1] - times[0], 1)
        #if p2c:
        #    print(f'\n   *Took a total of {dtime} s to create the DRO.')
        
        #timingMsg = "Took [*] s to create the DICOM Deformable Spatial "\
        #    + "Registration Object.\n"
        #params.add_timestamp(timingMsg)
        
        self.dro = dro

    def create_dro(self, srcDataset, trgDataset, newDataset, params):
        # TODO update docstrings
        """
        Create a Spatial or Deformable Spatial DICOM Registration Object.
        
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
        self.dro : Pydicom Object
            DICOM Spatial or Deformable Spatial Registration Object.
        """
        
        regTxName = params.cfgDict['regTxName']
        
        timingMsg = "Creating the DICOM Registration Object...\n"
        params.add_timestamp(timingMsg)
            
        if regTxName == 'bspline':
            self.create_def_dro(srcDataset, trgDataset, newDataset, params)
        else:
            self.create_spa_dro(srcDataset, trgDataset, newDataset, params)
        
        timingMsg = "Took [*] s to create the DICOM Registration Object.\n"
        params.add_timestamp(timingMsg)
    
    def export_dro(self, params):
        
        dro = self.dro
        cfgDict = params.cfgDict
        
        exportDro = cfgDict['exportDro']
        exportDir = cfgDict['droExportDir']
        srcExpLab = cfgDict['srcExpLab']
        srcScanID = cfgDict['srcScanID']
        trgExpLab = cfgDict['trgExpLab']
        trgScanID = cfgDict['trgScanID']
        regTxName = cfgDict['regTxName']
        
        if exportDro:
            # Create directory if it doesn't already exist:
            if not os.path.isdir(exportDir):
                #os.mkdir(exportDir)
                Path(exportDir).mkdir(parents=True)
            
            #cdt = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
            #
            #fname = f'{cdt}_ExpLab_{srcExpLab}_ScanID_{srcScanID}_{regTxName}'\
            #    + f'_reg_to_ExpLab_{trgExpLab}_ScanID_{trgScanID}.dcm'
            
            fname = generate_reg_fname(
                srcExpLab, srcScanID, trgExpLab, trgScanID, regTxName
                ) + '.dcm'
            
            fpath = os.path.join(exportDir, fname)
            
            dro.save_as(fpath)
            
            self.droFpath = fpath
        
            print(f'\nNew DICOM Registration Object exported to:\n {fpath}\n')
    
    def upload_dro(self, params):
        """
        Upload the DRO to XNAT if the DRO is to be uploaded, the use case to
        apply was one that potentially required image registration, and if an
        existing DRO was not used to by-pass registration.
        
        Parameters
        ----------
        params : DataDownloader Object
            Contains parameters (cfgDict), file paths (pathsDict), timestamps
            (timings) and timing messages (timingMsgs).
        
        Returns
        -------
        None.
        """
        xnatSession = params.xnatSession
        cfgDict = params.cfgDict
        
        useCaseToApply = cfgDict['useCaseToApply']
        useDroForTx = cfgDict['useDroForTx']
        regTxName = cfgDict['regTxName']
        uploadDro = cfgDict['uploadDro']
        url = cfgDict['url']
        projID = cfgDict['projID']
        subjLab = cfgDict['subjLab']
        
        if uploadDro and '5' in useCaseToApply and not useDroForTx:
            if regTxName in ['rigid', 'affine']:
                droType = 'SRO_DRO'
            else:
                droType = 'DSRO_DRO'
            
            xnatSession = upload_subj_asr(
                subj_asr_fpath=self.droFpath, 
                url=url, proj_id=projID, subj_label=subjLab, 
                content_label=droType, session=xnatSession
                )
        
    def create_dro_OLD_DELETE(srcDataset, trgDataset, newDataset, params):
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