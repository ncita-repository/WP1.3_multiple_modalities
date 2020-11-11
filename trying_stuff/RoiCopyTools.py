# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:12:23 2020

@author: ctorti
"""


"""
ROI copying tools
"""



# Import packages and functions:
import pydicom
import copy
import numpy as np
import SimpleITK as sitk

#import importlib
#import DicomTools
#importlib.reload(DicomTools)

from HelperTools import AddItemToListAndSort

from DicomTools import ImportDicoms
from DicomTools import ImportImage
from DicomTools import GetImageAttributes
from DicomTools import GetDicomSOPuids
from DicomTools import CompareSourceTargetImageAttributes

from ImageTools import ImageMin
from ImageTools import ImageMax
from ImageTools import OrImages
from ImageTools import ResampleImage

from RtsTools import GetCIStoDcmInds
from RtsTools import GetCStoDcmInds
from RtsTools import VerifyCISandCS
from RtsTools import GetPFFGStoDcmInds
from RtsTools import AddToRtsSequences
from RtsTools import ModifyRtsTagVals
#from RtsTools import ExportRtsRoi

from SegTools import GetSegLabels
from SegTools import GetDIVs
from SegTools import GroupListBySegment
from SegTools import GetRIStoDcmInds
from SegTools import AddToSegSequences
from SegTools import ModifySegTagVals
from SegTools import PixArr2Labmap
from SegTools import PixArr2Image
from SegTools import PixArrForOneFrame
from SegTools import Image2PixArr
#from SegTools import ExportSegRoi




"""
******************************************************************************
******************************************************************************
CONTOUR-RELATED FUNCTIONS
******************************************************************************
******************************************************************************
"""


def CopyContourWithinSeries(SrcDcmDir, SrcRtsFpath, FromSliceNum, ToSliceNum, 
                            LogToConsole=False):
    """
    Copy a single contour on a given slice to another slice within the same
    DICOM series.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source (= Target) 
                       DICOMs
                            
        SrcRtsFpath  - (String) Full path of the Source (= Target) 
                       DICOM-RTSTRUCT file 
                             
        FromSliceNum - (Integer) Slice index of the Source (= Target) DICOM 
                       stack corresponding to the contour to be copied 
                       (counting from 0)
                             
        ToSliceNum   - (Integer) Slice index of the Target (= Source) DICOM 
                       stack where the contour will be copied to (counting from 
                       0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        TrgRtsRoi    - Modified Target (= Source) ROI object 
    
    """
    
    # Import the Source RTSTRUCT ROI:
    SrcRtsRoi = pydicom.dcmread(SrcRtsFpath)
    
    # Import the Source (= Target) DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir)
    
    # Get the Source (= Target) DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices, 
    # and ContourSequence-to-DICOM slice indices for the Source DICOMs:
    SrcCIStoDcmInds = GetCIStoDcmInds(SrcRtsRoi, SrcSOPuids)
    
    SrcCStoDcmInds = GetCStoDcmInds(SrcRtsRoi, SrcSOPuids)
    
    # Both indices should be the same.  Return if not:
    if not VerifyCISandCS(SrcCIStoDcmInds, SrcCStoDcmInds, LogToConsole):
        return
        
    
    # Verify that a contour exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcCIStoDcmInds:
        print(f'\nA contour does not exist on slice {FromSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in SrcCIStoDcmInds:
        # Use SrcRtsRoi as a template for TrgRtsRoi: 
        TrgRtsRoi = copy.deepcopy(SrcRtsRoi)
        
        # Create a new list of CIStoDcmInds:
        TrgCIStoDcmInds = copy.deepcopy(SrcCIStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        TrgRtsRoi = AddToRtsSequences(SrcRtsRoi)
        
        # Create a new list of CIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to SrcCIStoDcmInds:
        TrgCIStoDcmInds = AddItemToListAndSort(SrcCIStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoDcmInds =', SrcCIStoDcmInds)
        print('\nTrgCIStoDcmInds =', TrgCIStoDcmInds)
    
    # Modify select tags in TrgRtsRoi:
    TrgRtsRoi = ModifyRtsTagVals(SrcRtsRoi, TrgRtsRoi, FromSliceNum, ToSliceNum, 
                                 SrcDicoms, SrcCIStoDcmInds, TrgCIStoDcmInds,
                                 LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    TrgRtsRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    TrgRtsRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return TrgRtsRoi





def DirectCopyContourAcrossSeries(SrcDcmDir, SrcRtsFpath, FromSliceNum, 
                                  TrgDcmDir, TrgRtsFpath, ToSliceNum, 
                                  LogToConsole=False):
    """
    Make a direct copy of a single contour on a given slice to another slice in 
    a different DICOM series.  A direct copy implies that the z-coordinates of
    the Target points (copy of Source points) do not relate to the 
    z-coordinates of the DICOM slice, since the voxels corresponding to the 
    Target ROI volume do not relate to the voxels of the Source ROI volume.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcRtsFpath  - (String) Full path of the Source DICOM-RTSTRUCT file 
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the contour to be copied (counting from 
                       0)
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgRtsFpath  - (String) Full path of the Target DICOM-RTSTRUCT file 
                             
        ToSliceNum   - (Integer) Slice index of the Target DICOM stack where
                       the contour will be copied to (counting from 0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgRtsRoi  - Modified Target ROI object 
    
    """
    
    # Import the Source and Target RTSTRUCT ROI:
    SrcRtsRoi = pydicom.dcmread(SrcRtsFpath)
    TrgRtsRoi = pydicom.dcmread(TrgRtsFpath)
    
    # Import the Target DICOMs:
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices, 
    # and ContourSequence-to-DICOM slice indices for Source and Target:
    SrcCIStoDcmInds = GetCIStoDcmInds(SrcRtsRoi, SrcSOPuids)
    
    SrcCStoDcmInds = GetCStoDcmInds(SrcRtsRoi, SrcSOPuids)
    
    TrgCIStoDcmInds = GetCIStoDcmInds(TrgRtsRoi, TrgSOPuids)
    
    TrgCStoDcmInds = GetCStoDcmInds(TrgRtsRoi, TrgSOPuids)
    
    # Both indices should be the same.  Return if not:
    if not VerifyCISandCS(SrcCIStoDcmInds, SrcCStoDcmInds, LogToConsole):
        return
    
    if not VerifyCISandCS(TrgCIStoDcmInds, TrgCStoDcmInds, LogToConsole):
        return
        
    
    # Verify that a contour exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcCIStoDcmInds:
        print(f'\nA contour does not exist on slice {FromSliceNum} in the',
              'Source DICOM series.')
        
        return
    
    # Verify that the slice that a contour is to be copied to exists (this 
    # check won't be necessary in the OHIF-Viewer but a bad input could occur
    # here):
    if not ToSliceNum in list(range(len(TrgDicoms))):
        print(f'\nThe Target DICOM series does not have slice {ToSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in TrgCIStoDcmInds:
        # Use TrgRtsRoi as a template for NewTrgRtsRoi: 
        NewTrgRtsRoi = copy.deepcopy(TrgRtsRoi)
        
        # Create a new list of NewTrgCIStoDcmInds:
        NewTrgCIStoDcmInds = copy.deepcopy(TrgCIStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgRtsRoi = AddToRtsSequences(TrgRtsRoi)
        
        # Create a new list of TrgCIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to TrgCIStoDcmInds:
        NewTrgCIStoDcmInds = AddItemToListAndSort(TrgCIStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoDcmInds    =', SrcCIStoDcmInds)
        print('\nTrgCIStoDcmInds    =', TrgCIStoDcmInds)
        print('\nNewTrgCIStoDcmInds =', NewTrgCIStoDcmInds)
    
    # Modify select tags in TrgRtsRoi:
    NewTrgRtsRoi = ModifyRtsTagVals(TrgRtsRoi, NewTrgRtsRoi, FromSliceNum, 
                                    ToSliceNum, TrgDicoms, 
                                    TrgCIStoDcmInds, NewTrgCIStoDcmInds,
                                    LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    NewTrgRtsRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRtsRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return NewTrgRtsRoi





def MappedCopyContourAcrossSeries_OLD(SrcDcmDir, SrcRoiFpath, FromSliceNum, 
                                 TrgDcmDir, TrgRoiFpath, LogToConsole):
    """
    I need to work with labelmaps!!!  See MappedCopySegmentAcrossSeries...
    
    Make a relationship-preserving copy of a single contour on a given slice to  
    a different DICOM series.  Since the relationship is preserved, the voxels 
    corresponding to the points copied to the Target ROI volume spatially
    relate to the voxels corresponding to the points copied from the Source ROI
    volume.
    
    There are four possible cases that will need to be accommodated*:
        
        2. The Source and Target images have the same FOR, IOP, IPP, PS and ST 
        
        3. The Source and Target images have the same FOR and IOP but have
           different PS and/or ST (if different ST they will likely have 
           different IPP)
        
        4. The Source and Target images have the same FOR but have different
           IOP and IPP, and possibly different PS and/or ST
           
        5. The Source and Target images have different FOR (hence different IPP
           and IOP), and possibly different ST and/or PS
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, ST = SliceThickness, PS = PixelSpacing
    
    * The numbering deliberately starts at 2 for consistency with the 5 overall
    cases.  The first case is covered by another function.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcRoiFpath  - (String) Full path of the Source DICOM RTSTRUCT/SEG file
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the contour/segment to be copied 
                       (counting from 0)
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgRoiFpath  - (String) Full path of the Target DICOM RTSTRUCT/SEG file 
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgRoi  - Modified Target ROI object 
    
    """
    
    # Import the Source and Target RTSTRUCT/SEG ROI:
    SrcRoi = pydicom.dcmread(SrcRoiFpath)
    TrgRoi = pydicom.dcmread(TrgRoiFpath)
    
    # Get the modality of the ROI:
    SrcModality = SrcRoi.Modality
    TrgModality = TrgRoi.Modality
    
    # Both modalities should be the same.  Return if not:
    if SrcModality == TrgModality:
        Modality = SrcModality
        
    else:
        print(f'\nThe Source ({SrcModality}) and Target ({TrgModality})',
              'modalities are different.')
        
        return
    
    # Import the Source and Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir, SortMethod='slices', Debug=False)
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir, SortMethod='slices', Debug=False)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ContourImageSequence-to-DICOM slice indices and
    # ContourSequence-to-DICOM slice indices, or the 
    # ReferencedImageSequence-to-DICOM slice indices and 
    # PerFrameFunctionalGroupsSequence-to-DICOM slice indices for Source and 
    # Target:
    if Modality == 'RTSTRUCT':
        SrcCIStoDcmInds = GetCIStoDcmInds(SrcRoi, SrcSOPuids)
    
        SrcCStoDcmInds = GetCStoDcmInds(SrcRoi, SrcSOPuids)
        
        TrgCIStoDcmInds = GetCIStoDcmInds(TrgRoi, TrgSOPuids)
        
        TrgCStoDcmInds = GetCStoDcmInds(TrgRoi, TrgSOPuids)
        
        # Both indices should be the same.  Return if not:
        if not VerifyCISandCS(SrcCIStoDcmInds, SrcCStoDcmInds, LogToConsole):
            return
        
        if not VerifyCISandCS(TrgCIStoDcmInds, TrgCStoDcmInds, LogToConsole):
            return
        
        # Verify that a contour exists for FromSliceNum (this check won't be 
        # necessary in the OHIF-Viewer but a bad input could occur here):
        if not FromSliceNum in SrcCIStoDcmInds:
            print(f'\nA contour does not exist on slice {FromSliceNum} in the',
                  'Source DICOM series.')
            
            return
        
    if Modality == 'SEG':
        SrcRIStoDcmInds = GetRIStoDcmInds(SrcRoi, SrcSOPuids)
    
        SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcRoi, SrcSOPuids)
        
        TrgRIStoDcmInds = GetRIStoDcmInds(TrgRoi, TrgSOPuids)
        
        TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgRoi, TrgSOPuids)
        
        # Verify that a segment exists for FromSliceNum (this check won't be 
        # necessary in the OHIF-Viewer but a bad input could occur here):
        if not FromSliceNum in SrcPFFGStoDcmInds:
            print(f'\nA segment does not exist on slice {FromSliceNum} in the',
                  'Source DICOM series.')
            
            return
    
    
    
    
    # Get the Source and Target FrameOfReferenceUID:
    SrcFOR = SrcDicoms[0].FrameOfReferenceUID
    TrgFOR = TrgDicoms[0].FrameOfReferenceUID
    
    # Get the image attributes:
    SrcIPPs, SrcDirs, SrcSpacings, SrcDims = GetImageAttributes(SrcDcmDir)
    TrgIPPs, TrgDirs, TrgSpacings, TrgDims = GetImageAttributes(TrgDcmDir)
    
    
    # Check which case follows:
    if SrcFOR == TrgFOR:
        if SrcDirs == TrgDirs:
            if SrcSpacings == TrgSpacings:
                Case = 2
                
            else:
                Case = 3
            
        else:
            Case = 4
    
    else:
        Case = 5
    
    
    
    
    
    
    
    # Check if ToSliceNum already has a contour:
    """
    Note:
        If slice number ToSliceNum already has a contour, it will be 
        over-written with the contour data of slice number FromSliceNum.
        If it doesn't have a contour, ContourImageSequence and 
        ContourSequence need to be extended by one.
    """
    if ToSliceNum in TrgCIStoDcmInds:
        # Use TrgRtsRoi as a template for NewTrgRtsRoi: 
        NewTrgRtsRoi = copy.deepcopy(TrgRtsRoi)
        
        # Create a new list of NewTrgCIStoDcmInds:
        NewTrgCIStoDcmInds = copy.deepcopy(TrgCIStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgRtsRoi = AddToRtsSequences(TrgRtsRoi)
        
        # Create a new list of TrgCIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ContourImageSequence and to
        # ContourSequence) to TrgCIStoDcmInds:
        NewTrgCIStoDcmInds = AddItemToListAndSort(TrgCIStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nSrcCIStoDcmInds    =', SrcCIStoDcmInds)
        print('\nTrgCIStoDcmInds    =', TrgCIStoDcmInds)
        print('\nNewTrgCIStoDcmInds =', NewTrgCIStoDcmInds)
    
    # Modify select tags in TrgRtsRoi:
    NewTrgRtsRoi = ModifyRtsTagVals(TrgRtsRoi, NewTrgRtsRoi, FromSliceNum, 
                                    ToSliceNum, TrgDicoms, 
                                    TrgCIStoDcmInds, NewTrgCIStoDcmInds,
                                    LogToConsole)
    
    
    
    # Generate a new SOP Instance UID:
    NewTrgRtsRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgRtsRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return NewTrgRtsRoi







"""
******************************************************************************
******************************************************************************
SEGMENTATION-RELATED FUNCTIONS
******************************************************************************
******************************************************************************
"""


def CopySegmentWithinSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, ToSliceNum, 
                            LogToConsole=False):
    """
    Copy a single segment on a given slice to another slice within the same
    DICOM series.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcSegFpath  - (String) Full path of the Source DICOM SEG file 
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the segment to be copied (counting from 
                       0)
                             
        ToSliceNum   - (Integer) Slice index in the Target DICOM stack where
                       the segment will be copied to (counting from 0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewSrcSegRoi - Modified Source SEG ROI object
    
    """
    
    # Import the Source SEG ROI:
    SegRoi = pydicom.dcmread(SrcSegFpath)
    
    # Import the Source DICOMs:
    Dicoms = ImportDicoms(DicomDir=SrcDcmDir)
    
    # Get the DICOM SOP UIDs:
    SOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, and the
    # Per-frameFunctionalGroupsSequence-to-DICOM slice indices:
    RIStoDcmInds = GetRIStoDcmInds(SegRoi, SOPuids)
    
    PFFGStoDcmInds = GetPFFGStoDcmInds(SegRoi, SOPuids)
    
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in PFFGStoDcmInds:
        print(f'\nA segment does not exist on slice {FromSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a segment:
    """
    Note:
        If slice number ToSliceNum already has a segment, it will be 
        over-written with the segment data of slice number FromSliceNum.
        If it doesn't have a segment, ReferencedImageSequence and 
        PerFrameFunctionalGroupsSequence need to be extended by one.
    """
    if ToSliceNum in RIStoDcmInds:
        # Use SegRoi as a template for NewSegRoi: 
        NewSegRoi = copy.deepcopy(SegRoi)
        
        # Create a new list of RIStoDcmInds: <-- not needed?
        #NewRIStoDcmInds = copy.deepcopy(RIStoDcmInds)
        
        # Create a new list of PFFGStoDcmInds:
        NewPFFGStoDcmInds = copy.deepcopy(PFFGStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewSegRoi = AddToSegSequences(SegRoi)
        
        # Create a new list of RIStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to ReferencedImageSequence)
        # to RIStoDcmInds: <-- not needed?
        #NewRIStoDcmInds = AddItemToListAndSort(RIStoDcmInds, ToSliceNum)
        
        # Create a new list of PFFGStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to 
        # PerFrameFunctionalGroupsSequence) to PFFGStoDcmInds:
        NewPFFGStoDcmInds = AddItemToListAndSort(PFFGStoDcmInds, ToSliceNum)
    
    
    if LogToConsole:
        print('\nPFFGStoDcmInds    =', PFFGStoDcmInds)
        print('\nNewPFFGStoDcmInds =', NewPFFGStoDcmInds)
        
    # Modify select tags in NewSegRoi:
    NewSegRoi = ModifySegTagVals(SrcSeg=SegRoi, SrcDicoms=Dicoms, 
                                 SrcPFFGStoDcmInds=PFFGStoDcmInds, 
                                 FromSliceNum=FromSliceNum, 
                                 TrgSeg=SegRoi, TrgDicoms=Dicoms, 
                                 TrgPFFGStoDcmInds=PFFGStoDcmInds, 
                                 ToSliceNum=ToSliceNum, 
                                 NewTrgSeg=NewSegRoi, 
                                 NewTrgPFFGStoDcmInds=NewPFFGStoDcmInds, 
                                 LogToConsole=LogToConsole)
                                 
    
    
    
    # Generate a new SOP Instance UID:
    NewSegRoi.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewSegRoi.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return NewSegRoi





def DirectCopySegmentAcrossSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, 
                                  TrgDcmDir, TrgSegFpath, ToSliceNum, 
                                  LogToConsole=False):
    """
    Copy a single segment on a given slice to another slice in a different 
    DICOM series.
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcSegFpath  - (String) Full path of the Source DICOM SEG file 
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the segment to be copied (counting from 
                       0)
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgSegFpath  - (String) Full path of the Target DICOM SEG file 
                             
        ToSliceNum   - (Integer) Slice index of the Target DICOM stack where
                       the segment will be copied to (counting from 0)
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgSeg    - Modified Target SEG object 
    
    """
    
    # Import the Source and Target SEG ROI:
    SrcSeg = pydicom.dcmread(SrcSegFpath)
    TrgSeg = pydicom.dcmread(TrgSegFpath)
    
    # Import the Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir) # <-- not needed but a 
    # required input to ModifySegTagVals 
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, and the
    # Per-frameFunctionalGroupsSequence-to-DICOM slice indices for Source and 
    # Target:
    SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
    
    TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        print(f'\nA segment does not exist on slice {FromSliceNum} in the',
              'Source DICOM series.')
        
        return
    
    # Verify that the slice that a contour is to be copied to exists (this 
    # check won't be necessary in the OHIF-Viewer but a bad input could occur
    # here):
    if not ToSliceNum in list(range(len(TrgDicoms))):
        print(f'\nThe Target DICOM series does not have slice {ToSliceNum}.')
        
        return
    
    
    # Check if ToSliceNum already has a segment:
    """
    Note:
        If slice number ToSliceNum already has a segment, it will be 
        over-written with the segment data of slice number FromSliceNum.
        If it doesn't have a segment, ReferencedImageSequence and 
        PerFrameFunctionalGroupsSequence need to be extended by one.
    """
    if ToSliceNum in TrgRIStoDcmInds:
        # Use TrgSegRoi as a template for NewTrgSeg:
        NewTrgSeg = copy.deepcopy(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = copy.deepcopy(TrgRIStoDcmInds)
        
        # Create a new list of PFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = copy.deepcopy(TrgPFFGStoDcmInds)
        
    else:
        # Increase the length of the sequences:
        NewTrgSeg = AddToSegSequences(TrgSeg)
        
        # Create a new list of TrgRIStoDcmInds by adding the index ToSliceNum 
        # (representing the contour to be added to ReferencedImageSequence) to 
        # TrgRIStoDcmInds: <-- not needed?
        #NewTrgRIStoDcmInds = AddItemToListAndSort(TrgRIStoDcmInds, ToSliceNum)
        
        # Create a new list of TrgPFFGStoDcmInds by adding the index ToSliceNum 
        # (representing the segment to be added to 
        # PerFrameFunctionalGroupsSequence) to TrgPFFGStoDcmInds:
        NewTrgPFFGStoDcmInds = AddItemToListAndSort(TrgPFFGStoDcmInds, ToSliceNum)
        
    
    if LogToConsole:
        print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
        print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
    
    # Modify select tags in NewTrgSeg:
    NewTrgSeg = ModifySegTagVals(SrcSeg, SrcDicoms, SrcPFFGStoDcmInds, 
                                 FromSliceNum, TrgSeg, TrgDicoms, 
                                 TrgPFFGStoDcmInds, ToSliceNum, NewTrgSeg, 
                                 NewTrgPFFGStoDcmInds, LogToConsole)
    
    
    # Generate a new SOP Instance UID:
    NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
    
    # Generate a new Series Instance UID:
    NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
    
    
    return NewTrgSeg







def MappedCopySegmentAcrossSeries(SrcDcmDir, SrcSegFpath, FromSliceNum, 
                                  FromSegLabel,
                                  TrgDcmDir, TrgSegFpath, ToSegLabel, 
                                  LogToConsole=False):
                                
    """
    Make a relationship-preserving copy of a single segment on a given slice to  
    a different DICOM series.  Since the relationship is preserved, the voxels 
    corresponding to the voxels copied to the Target ROI volume spatially
    relate to the voxels corresponding to the voxels copied from the Source ROI
    volume.
    
    There are four possible cases that will need to be accommodated*:
        
        2b. The Source and Target images have the same FOR, IOP, IPP, PS and ST 
        
        3b. The Source and Target images have the same FOR and IOP but have
           different PS and/or ST (if different ST they will likely have 
           different IPP)
        
        4. The Source and Target images have the same FOR but have different
           IOP and IPP, and possibly different PS and/or ST
           
        5. The Source and Target images have different FOR (hence different IPP
           and IOP), and possibly different ST and/or PS
        
    FOR = FrameOfReferenceUID, IPP = ImagePositionPatient, 
    IOP = ImageOrientationPatient, ST = SliceThickness, PS = PixelSpacing
    
    * The numbering deliberately starts at 2b for consistency with the 5 
    overall cases.  Cases 1, 2a and 3a (=Direct copies) are covered by other 
    functions (CopySegmentWithinSeries for Case 1, and 
    DirectCopySegmentAcrossSeries for Cases 2a and 3a).
    
    
    Inputs:
        SrcDcmDir    - (String) Directory containing the Source DICOMs
                            
        SrcSegFpath  - (String) Full path of the Source DICOM SEG file
                             
        FromSliceNum - (Integer) Slice index of the Source DICOM stack 
                       corresponding to the contour/segment to be copied 
                       (counting from 0)
                       
        FromSegLabel - (String) All or part of the segment label of the ROI
                       in the Source SEG containing the segment to be copied
                       
        TrgDcmDir    - (String) Directory containing the Target DICOMs
                            
        TrgSegFpath  - (String) Full path of the Target DICOM SEG file;
                       May also be empty string ('') if there is no DICOM SEG
                       for the Target Series
                       
        ToSegLabel   - (String) All or part of the segment label of the ROI
                       in the Target SEG that the segment is to be copied to
                            
        LogToConsole - (Boolean) Denotes whether or not some intermediate
                       results will be logged to the console
                            
    Returns:
        NewTrgSeg    - Modified Target ROI object 
    
    """
    
    # Import the Source and Target SEG ROI:
    SrcSeg = pydicom.dcmread(SrcSegFpath)
    if TrgSegFpath: # if TrgSegFpath is not empty
        TrgSeg = pydicom.dcmread(TrgSegFpath)
    
    # Get the modality of the ROIs:
    SrcModality = SrcSeg.Modality
    if TrgSegFpath:
        TrgModality = TrgSeg.Modality
    
    # Both modalities should be SEG.  Return if not:
    if TrgSegFpath:
        if SrcModality != 'SEG' or TrgModality != 'SEG':
            print(f'\nThe Source ({SrcModality}) and Target ({TrgModality})',
                  'modalities are different.')
            
            return
    
    # Import the Source and Target DICOMs:
    SrcDicoms = ImportDicoms(DicomDir=SrcDcmDir)
    TrgDicoms = ImportDicoms(DicomDir=TrgDcmDir)
    
    # Get the Source and Target FrameOfReferenceUID:
    SrcFOR = SrcDicoms[0].FrameOfReferenceUID
    TrgFOR = TrgDicoms[0].FrameOfReferenceUID
    
    # Get the Image Attributes for Source and Target using SimpleITK:
    SrcSitkSize, SrcSitkSpacing, SrcSitkST,\
    SrcSitkIPP, SrcSitkDir = GetImageAttributes(DicomDir=SrcDcmDir,
                                                Package='sitk')
    
    TrgSitkSize, TrgSitkSpacing, TrgSitkST,\
    TrgSitkIPP, TrgSitkDir = GetImageAttributes(DicomDir=TrgDcmDir,
                                                Package='sitk')
    
    if LogToConsole:
        CompareSourceTargetImageAttributes(SrcDcmDir, TrgDcmDir)
    
    
    """
    I'm not sure if the items in DimensionIndexSequence,
    e.g. DimensionDescriptionLabel = 'ReferencedSegmentNumber'
                                   = 'ImagePositionPatient'
                                   = 'StackID'
                                   = 'InStackPositionNumber'
                                   
    are important.  While the label 'ReferenedIndexSequence' are common to SEGs
    both generated by the OHIF-Viewer and RTS-to-SEG conversions, it seems that
    the label 'ImagePositionPatient' is only found in OHIF-generated SEGs, and
    the labels 'StackID' and 'InStackPositionNumber' from RTS-to-SEG 
    conversions.
    """
    
    # Get the number of ROIs in Source and Target:
    SrcNumOfROIs = len(SrcSeg.SegmentSequence)
    if TrgSegFpath:
        TrgNumOfROIs = len(TrgSeg.SegmentSequence)
    
    # Get the labels of all ROIs in the Source and Target SEGs:
    SrcSegLabels = GetSegLabels(SrcSeg)
    #SrcNumOfROIs = len(SrcSegLabels)
    
    if TrgSegFpath:
        TrgSegLabels = GetSegLabels(TrgSeg)
        #TrgNumOfROIs = len(TrgSegLabels)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print(f'\nThe Source SEG has {SrcNumOfROIs} ROIs with label(s)',
                  f'{SrcSegLabels} \nand the Target SEG has {TrgNumOfROIs} ROIs',
                  f'with label(s) {TrgSegLabels}')
        else:
            print(f'\nThe Source SEG has {SrcNumOfROIs} ROIs with label(s)',
                  f'{SrcSegLabels}')
    
    
    
    # Get the DimensionIndexValues (which relate the segments to the item in 
    # the Referenced Instance Sequence):
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVs = GetDIVs(SrcSeg)
    if TrgSegFpath:
        TrgDIVs = GetDIVs(TrgSeg)
    
    # Group the DIVs by ROI:
    """
    Note:
        The DimensionalIndexValues are integers that start from 1 (unlike the
        RIS or PFFGS which begin at 0)!
    """
    SrcDIVsByRoi = GroupListBySegment(ListToGroup=SrcDIVs, DIVs=SrcDIVs)
    if TrgSegFpath:
        TrgDIVsByRoi = GroupListBySegment(ListToGroup=TrgDIVs, DIVs=TrgDIVs)
    
    if True:#LogToConsole:
        if TrgSegFpath:
            print('\nThe Source SEG DimensionIndexValues grouped by ROI are',
                  f'\n{SrcDIVsByRoi} \n\n and the Target SEG',
                  f'DimensionIndexValues grouped by ROI are \n{TrgDIVsByRoi}')
        else:
            print('\nThe Source SEG DimensionIndexValues grouped by ROI are',
                  f'\n{SrcDIVsByRoi}')
        
    # Get the Source and Target DICOM SOP UIDs:
    SrcSOPuids = GetDicomSOPuids(DicomDir=SrcDcmDir)
    TrgSOPuids = GetDicomSOPuids(DicomDir=TrgDcmDir)
    
    # Get the ReferencedImageSequence-to-DICOM slice indices, 
    # the same grouped by ROI,
    # the PerFrameFunctionalGroupsSequence-to-DICOM slice indices, 
    # the same grouped by ROI, and 
    # the frame numbers grouped by ROI for Source and Target:
    """
    Note:
        The RIS and PFFGS are integers that start from 0 (unlike the DIVs which
        begin at 1)!
    """
    SrcRIStoDcmInds = GetRIStoDcmInds(SrcSeg, SrcSOPuids)

    SrcPFFGStoDcmInds = GetPFFGStoDcmInds(SrcSeg, SrcSOPuids)
    
    SrcFrameNums = list(range(len(SrcPFFGStoDcmInds)))
    
    # Group the above by ROIs:
    #SrcRIStoDcmIndsByRoi = GroupListBySegment(ListToGroup=SrcRIStoDcmInds, 
    #                                      DIVs=SrcDIVs) # <-- RIS has more elements than DIV
    SrcPFFGStoDcmIndsByRoi = GroupListBySegment(ListToGroup=SrcPFFGStoDcmInds, 
                                                DIVs=SrcDIVs)
    
    SrcFrameNumsByRoi = GroupListBySegment(ListToGroup=SrcFrameNums, 
                                           DIVs=SrcDIVs)
    
    if TrgSegFpath:
        TrgRIStoDcmInds = GetRIStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgPFFGStoDcmInds = GetPFFGStoDcmInds(TrgSeg, TrgSOPuids)
        
        TrgFrameNums = list(range(len(TrgPFFGStoDcmInds)))
        
        # Group the above by ROIs:
        #TrgRIStoDcmIndsByRoi = GroupListBySegment(ListToGroup=TrgRIStoDcmInds, 
        #                                      DIVs=TrgDIVs) # <-- RIS has more elements than DIV
        TrgPFFGStoDcmIndsByRoi = GroupListBySegment(ListToGroup=TrgPFFGStoDcmInds, 
                                                    DIVs=TrgDIVs)
        
        TrgFrameNumsByRoi = GroupListBySegment(ListToGroup=TrgFrameNums, 
                                               DIVs=TrgDIVs)
        
        if True:#LogToConsole:
            print('\nThe Source series has segments on slices',
                  f'{SrcPFFGStoDcmIndsByRoi}, \nand the Target series on',
                  f'slices {TrgPFFGStoDcmIndsByRoi} \n(grouped by ROI)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}, \nand for the Target series is',
            #      f'\n{TrgRIStoDcmInds}')
    
    else:
        if True:#LogToConsole:
            print('The Source series has segments on slices',
                  f'{SrcPFFGStoDcmIndsByRoi} \n(grouped by ROI)')
            
            #print('\nThe RIS-to-DICOM slice num for the Source series is',
            #      f'\n{SrcRIStoDcmInds}')
            
    
    # Verify that a segment exists for FromSliceNum (this check won't be 
    # necessary in the OHIF-Viewer but a bad input could occur here):
    if not FromSliceNum in SrcPFFGStoDcmInds:
        print('\nThe Source series does not have any segments on slice',
              f'{FromSliceNum}')
        
        return
    
    
    # The 3D pixel array:
    #SrcSegPixArr = SrcSeg.pixel_array
    #TrgSegPixArr = TrgSeg.pixel_array
    
    if False:
        # Create 3D labelmaps (= pixel arrays in same frame position as 
        # corresponding DICOM slices, with 2D zero masks in between):
        SrcLabmap = PixArr2Labmap(PixelArray=SrcSeg.pixel_array,
                                  NumOfSlices=SrcSitkSize[2],
                                  FrameToSliceInds=SrcPFFGStoDcmInds)
        
        if TrgSegFpath:
            TrgLabmap = PixArr2Labmap(PixelArray=TrgSeg.pixel_array,
                                      NumOfSlices=TrgSitkSize[2],
                                      FrameToSliceInds=TrgPFFGStoDcmInds)
        else:
            # Initialise the Target LabelMap:
            TrgLabmap = np.zeros((TrgSitkSize[2], TrgSitkSize[1], TrgSitkSize[0]), 
                                 dtype='uint')
        
        """ More useful to convert Pixel Array to sitk image.. """
    
    # Read in the Source and Target SimpleITK 3D images:
    SrcImage = ImportImage(SrcDcmDir)
    TrgImage = ImportImage(TrgDcmDir)
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Source*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    # Create SimpleITK images of the 3D labelmaps (= pixel arrays in same frame 
    # position as corresponding DICOM slices, with 2D zero masks in between)
    # for each ROI:
    SrcLabmapIms = []
    
    for RoiNum in range(SrcNumOfROIs):
        # The starting (= a) and ending (= b) frame numbers:
        a = SrcFrameNumsByRoi[RoiNum][0]
        b = SrcFrameNumsByRoi[RoiNum][-1]
        
        # The (partial) pixel array for this ROI:
        PixArr = SrcSeg.pixel_array[a:b+1]
        
        Inds = SrcPFFGStoDcmIndsByRoi[RoiNum]
        
        if LogToConsole:
            print(f'len(PixArr) = {len(PixArr)}')
            print(f'len(Inds) = {len(Inds)}')
        
        SrcLabmapIm = PixArr2Image(PixelArray=PixArr,
                                   Image=SrcImage,
                                   FrameToSliceInds=Inds)
        
        SrcLabmapIms.append(SrcLabmapIm)
        
    
    if LogToConsole:
        title = 'Creating 3D labelmaps for *Target*...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    if TrgSegFpath:
        TrgLabmapIms = []
        
        for RoiNum in range(TrgNumOfROIs):
            # The starting (= a) and ending (= b) frame numbers:
            a = TrgFrameNumsByRoi[RoiNum][0]
            b = TrgFrameNumsByRoi[RoiNum][-1]
            
            # The (partial) pixel array for this ROI:
            PixArr = TrgSeg.pixel_array[a:b+1]
            
            Inds = TrgPFFGStoDcmIndsByRoi[RoiNum]
            
            TrgLabmapIm = PixArr2Image(PixelArray=PixArr,
                                       Image=TrgImage,
                                       FrameToSliceInds=Inds)
            
            TrgLabmapIms.append(TrgLabmapIm)
            
        
    else:
        # Initialise the Target LabelMap:
        TrgLabmap = np.zeros((TrgSitkSize[2], TrgSitkSize[1], TrgSitkSize[0]), 
                             dtype='uint')
        
        TrgLabmapIm = sitk.GetImageFromArray(TrgLabmap)
        
        # Set the Origin, Spacing and Direction of TrgLabmapIm to that of 
        # TrgImage:
        TrgLabmapIm.SetOrigin(TrgImage.GetOrigin())
        TrgLabmapIm.SetSpacing(TrgImage.GetSpacing())
        TrgLabmapIm.SetDirection(TrgImage.GetDirection())
        
        # For consistency with other IF case:
        TrgLabmapIms = []
        TrgLabmapIms.append(TrgLabmapIm)
        
        
    """ 
    SrcLabmapIm is the 3D labelmap of all segments within SrcSeg.
    
    But when looking to copy a single segment from Source to Target, I'll need
    to know how that single segment is resampled.  So create a new labelmap
    containing the segment to be copied only.
    """
    
    if LogToConsole:
        title = 'Creating 3D labelmap for *Source* containing segment to copy only...'
        ul = '*' * len(title)
        print('\n' + title)
        print(ul + '\n')
    
    # Create SimpleITK image of the 3D labelmap for the segment to be copied 
    # only:
    """
    There may be more than one frame number in PixelArray that corresponds to 
    FromSliceNum...
    # The frame number in PixelArray that corresponds to FromSliceNum:
    FromFrameNum = SrcPFFGStoDcmInds.index(FromSliceNum)
    """
    
    # Determine which Source ROI contains the segment to be copied:
    """
    What I call RoiNum is equivalent to Segment Number in the Segment Sequence.
    """
    FromRoiNum = [i for i, s in enumerate(SrcSegLabels) if FromSegLabel in s]
    FromRoiNum = FromRoiNum[0]
    
    if TrgSegFpath:
        # Determine which Target ROI the Source segment is to be copied to:
        ToRoiNum = [i for i, s in enumerate(TrgSegLabels) if ToSegLabel in s]
        ToRoiNum = ToRoiNum[0]
    
        print(f'\nFromRoiNum = {FromRoiNum} \nToRoiNum = {ToRoiNum}')
        
    else:
        print(f'\nFromRoiNum = {FromRoiNum}')
    
    # Determine which frame number in PixelArray corresponds to FromSliceNum:
    """
    FromFrameNum = SrcPFFGStoDcmIndsByRoi[FromRoiNum].index(FromSliceNum)  
    
    The problem with the above is that it's the frame number within the subset
    of frames pertaining to the ROI given by FromRoiNum - it doesn't relate to 
    the frame number in PixelArray (which contains all frames for all ROIs).
    So need something a bit more clever...
    """
    n = 0 # initialise frame number counter
    
    for RoiNum in range(SrcNumOfROIs):
        if RoiNum == FromRoiNum:
            # This ROI contains the frame to be copied, whose position is given
            # by the index of FromSliceNum is in SrcPFFGStoDcmIndsByRoi[RoiNum]
            # plus any frames that preceeded it (i.e. the frame counter n):
            FromFrameNum = n + SrcPFFGStoDcmIndsByRoi[RoiNum].index(FromSliceNum)
            
        else:
            # Add to the frame counter the number of frames in this ROI:
            n += len(SrcPFFGStoDcmIndsByRoi[RoiNum])
    
    
    print(f'\nFromSliceNum = {FromSliceNum} relates to frame number',
          f'{FromFrameNum} in PixelArray')
    
    # Determine which frame number(s) in PixelArray correspond to FromSliceNum: 
    #FromFrameNums = []
    #
    #for RoiNum in range(SrcNumOfROIs):
    #    if FromSliceNum in SrcPFFGStoDcmIndsByRoi[RoiNum]:
    #        ind = SrcPFFGStoDcmIndsByRoi[RoiNum].index(FromSliceNum)
    #        
    #        FromFrameNums.append(ind)
    #     
    #print(f'\nFromSliceNum = {FromSliceNum} relates to frame number(s)',
    #      f'{FromFrameNums} in PixelArray')
    
    # Knowing which ROI contains the frame to be copied and the frame number
    # defines the pixel array to be copied:
    FromPixelArr = SrcSeg.pixel_array[FromFrameNum]
    
    # Modify the pixel array so that all but the frame to be copied has zeros:
    """
    SrcPixArrToCopy = PixArrForOneFrame(PixelArray=FromPixelArr, 
                                        FrameNum=FromFrameNum)
    
    I either need to pass the entire PixelArray and specify the FrameNum to be
    copied (as was my original intention with the function PixArrForOneFrame),
    or pass the frame to be copied (FromPixelArr) and modify PixArrForOneFrame
    accordingly!
    """
    SrcPixArrToCopy = PixArrForOneFrame(PixelArray=SrcSeg.pixel_array, 
                                        FrameNum=FromFrameNum)
    
    if LogToConsole:
        print(f'\nSrcPixArrToCopy.max() = {SrcPixArrToCopy.max()}')
    
    SrcLabmapToCopyIm = PixArr2Image(PixelArray=SrcPixArrToCopy,
                                     Image=SrcImage,
                                     FrameToSliceInds=[FromSliceNum])
    
    if LogToConsole:
        print(f'\nMax value of SrcLabmapToCopyIm = {ImageMax(SrcLabmapToCopyIm)}')
    
    
    # Check which case follows:
    if SrcFOR == TrgFOR:
        if LogToConsole:
            print('\nThe Source and Target FrameOfReferenceUID are the same')
        
        # Compare the Source and Target directions. Get the maximum value of 
        # their absolute differences:
        AbsDiffs = abs(np.array(SrcSitkDir) - np.array(TrgSitkDir))
        MaxAbsDiff = max(AbsDiffs)
        
        # Consider the directions different if any of the vector differences 
        # are greater than epsilon:
        epsilon = 1e-06
        
        #if SrcDirs == TrgDirs: 
        """ The line above will result in insignificant differences (e.g. due
        to quantisation/pixelation) resulting in False, so instead consider
        MaxAbsDiff and epsilon. 
        """
        if MaxAbsDiff < epsilon:
            if LogToConsole:
                print('The Source and Target directions are the same')
            
            if SrcSitkSpacing == TrgSitkSpacing:
                if LogToConsole:
                    print('The Source and Target voxel sizes are the same')
                
                CaseNum = '2b'
                
            else:
                if LogToConsole:
                    print('The Source and Target voxel sizes are different:\n',
                          f'   Source: {SrcSitkSpacing} \n',
                          f'   Target: {TrgSitkSpacing}')
                    
                CaseNum = '3b'
            
        else:
            if LogToConsole:
                print('The Source and Target directions are different:\n',
                      f'   Source: [{SrcSitkDir[0]}, {SrcSitkDir[1]}, {SrcSitkDir[2]},',
                      f'            {SrcSitkDir[3]}, {SrcSitkDir[4]}, {SrcSitkDir[5]},',
                      f'            {SrcSitkDir[6]}, {SrcSitkDir[7]}, {SrcSitkDir[8]}]\n',
                      f'   Target: [{TrgSitkDir[0]}, {TrgSitkDir[1]}, {TrgSitkDir[2]},',
                      f'            {TrgSitkDir[3]}, {TrgSitkDir[4]}, {TrgSitkDir[5]},',
                      f'            {TrgSitkDir[6]}, {TrgSitkDir[7]}, {TrgSitkDir[8]}]')
                
                print(f'AbsDiffs = {AbsDiffs}')
            
            
            
            CaseNum = '4'
    
    else:
        if LogToConsole:
            print('The Source and Target FrameOfReferenceUID are different')
        
        CaseNum = '5'
        
    
    if LogToConsole:    
        print('The Source and Target origins:\n',
                          f'   Source: {SrcSitkIPP[0]} \n',
                          f'   Target: {TrgSitkIPP[0]}')
        
        print('The Source and Target image dims:\n',
                          f'   Source: {SrcSitkSize} \n',
                          f'   Target: {TrgSitkSize}')
    
    
    
    if LogToConsole:
        print(f'\nCase {CaseNum} applies.')
    
    
    if CaseNum == '2b':
        """
        The origins for Source and Target are the same since they have the same
        FrameOfReferenceUID.  And since the PixelSpacing and SliceThickness are
        also the same, ToSliceNum = FromSliceNum, and the relationships will be
        preserved.
        
        26/10: It's seems that having the same FOR does not guarantee the same
        origin.  Series 9 and 6 for Subject 011, MR4 have different origins!
        """
        
        print('\nApplying Case 2b..')
    
        ToSliceNum = copy.deepcopy(FromSliceNum)
        
        if LogToConsole:
            print(f'\nCopying slice {FromSliceNum} from Source SEG to slice',
                  f'{ToSliceNum} in Target SEG..')
    
        # Check if ToSliceNum already has a segment:
        """
        Note:
            If slice number ToSliceNum already has a segment, it will be 
            over-written with the segment data of slice number FromSliceNum.
            If it doesn't have a segment, ReferencedImageSequence and 
            PerFrameFunctionalGroupsSequence need to be extended by one.
        """
        if ToSliceNum in TrgRIStoDcmInds:
            # Use TrgSeg as a template for NewTrgSegRoi:
            NewTrgSeg = copy.deepcopy(TrgSeg)
            
            # Create a new list of TrgRIStoDcmInds: <-- not needed?
            #NewTrgRIStoDcmInds = copy.deepcopy(TrgRIStoDcmInds)
            
            # Create a new list of PFFGStoDcmInds:
            NewTrgPFFGStoDcmInds = copy.deepcopy(TrgPFFGStoDcmInds)
            
        else:
            # Increase the length of the sequences:
            NewTrgSeg = AddToSegSequences(TrgSeg)
            
            # Create a new list of TrgRIStoDcmInds by adding the index  
            # ToSliceNum (representing the contour to be added to 
            # ReferencedImageSequence) to TrgRIStoDcmInds: <-- not needed?
            #NewTrgRIStoDcmInds = AddItemToListAndSort(TrgRIStoDcmInds, 
            #                                          ToSliceNum)
            
            # Create a new list of TrgPFFGStoDcmInds by adding the index 
            # ToSliceNum (representing the segment to be added to 
            # PerFrameFunctionalGroupsSequence) to TrgPFFGStoDcmInds:
            NewTrgPFFGStoDcmInds = AddItemToListAndSort(TrgPFFGStoDcmInds, 
                                                        ToSliceNum)
            
        
        if LogToConsole:
            print('\nTrgPFFGStoDcmInds    =', TrgPFFGStoDcmInds)
            print('\nNewTrgPFFGStoDcmInds =', NewTrgPFFGStoDcmInds)
        
        # Modify select tags in NewTrgSeg:
        #NewTrgSeg = ModifySegTagVals(TrgSeg, NewTrgSeg, FromSliceNum, 
        #                             ToSliceNum, TrgDicoms, 
        #                             #TrgPFFGStoDcmInds, NewTrgPFFGStoDcmInds,
        #                             SrcPFFGStoDcmInds, NewTrgPFFGStoDcmInds,
        #                             LogToConsole)
        
        NewTrgSeg = ModifySegTagVals(SrcSeg, SrcDicoms, FromSliceNum, 
                                     TrgSeg, NewTrgSeg, TrgDicoms, ToSliceNum, 
                                     SrcPFFGStoDcmInds, 
                                     TrgPFFGStoDcmInds, NewTrgPFFGStoDcmInds, 
                                     LogToConsole)
        
        
        
        """
        23/10:  Still need to verify that Case 2 provides the expected results.
        
        26/10:  Resampling will need to be performed if the origins are not the
        same...
        """
        
        
    if CaseNum == '3b':
        print(f'\nApplying Case 3b..')
        
        # Read in the Source and Target SimpleITK 3D images:
        #SrcImage = ImportImage(SrcDcmDir)
        #TrgImage = ImportImage(TrgDcmDir)
        
        #SrcSitkOrigin = SrcSitkIm.GetOrigin() 
        #TrgSitkOrigin = SrcSitkIm.GetOrigin() 
        #
        #SrcSitkDirs = SrcSitkIm.GetDirection() 
        #
        #SrcSitkSpacings = SrcSitkIm.GetSpacing() 
        #
        #SrcSitkDims = SrcSitkIm.GetSize() 
        
        
        # Get the Image Attributes for Source and Target using SimpleITK:
        #SrcSitkSize, SrcSitkSpacing,\
        #SrcSitkIPPs, SrcSitkDirs = GetImageAttributes(DicomDir=SrcDcmDir,
        #                                              Package='sitk')
        #
        #TrgSitkSize, TrgSitkSpacing,\
        #TrgSitkIPPs, TrgSitkDirs = GetImageAttributes(DicomDir=TrgDcmDir,
        #                                              Package='sitk')
        
        
        """
        Although it's the Source Labelmap that needs to be resampled, first
        prove that I can resample the image.
        """
        
        
        # Resample the Source image:
        ResSrcImage = ResampleImage(Image=SrcImage, RefImage=TrgImage, 
                                    Interpolation='Linear')
        
        if LogToConsole:                            
            print('Results after resampling Source image:\n')
            
            print('The Source and Target image size:\n',
                          f'   Original Source:  {SrcImage.GetSize()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSize()} \n',
                          f'   Target:           {TrgImage.GetSize()}')
                
            print('The Source and Target voxel spacings:\n',
                          f'   Original Source:  {SrcImage.GetSpacing()} \n',
                          f'   Resampled Source: {ResSrcImage.GetSpacing()} \n',
                          f'   Target:           {TrgImage.GetSpacing()}')
                
            print('The Source and Target Origin:\n',
                          f'   Original Source:  {SrcImage.GetOrigin()} \n',
                          f'   Resampled Source: {ResSrcImage.GetOrigin()} \n',
                          f'   Target:           {TrgImage.GetOrigin()}')
                
            print('The Source and Target Direction:\n',
                          f'   Original Source:  {SrcImage.GetDirection()} \n',
                          f'   Resampled Source: {ResSrcImage.GetDirection()} \n',
                          f'   Target:           {TrgImage.GetDirection()}')
                                   
                                    
        # Resample the Source labelmap images:
        ResSrcLabmapIms = []
        
        for SrcLabmapIm in SrcLabmapIms:
            ResSrcLabmapIm = ResampleImage(Image=SrcLabmapIm, RefImage=TrgImage, 
                                           Interpolation='NearestNeighbor')
            
            ResSrcLabmapIms.append(ResSrcLabmapIm)
            
        
        """
        Now need to work with the labelmap only containing the segment to be 
        copied.
        """
        # Resample the Source labelmap:
        ResSrcLabmapToCopyIm = ResampleImage(Image=SrcLabmapToCopyIm, 
                                             RefImage=TrgImage, 
                                             Interpolation='NearestNeighbor')
        
        # Perform a pixel-wise OR of ResSrcLabmapToCopyIm and the labelmap of
        # Target that the segment is to be copied to:
        ORLabmapIm = OrImages(ResSrcLabmapToCopyIm, TrgLabmapIms[ToRoiNum])
        
        """
        I can't tell if OR has worked since my ROIs (tumour and brain) overlap!
        So rather than performing OR operation, add the labelmaps.
        """
        # Perform a pixel-wise addition of ResSrcLabmapToCopyIm and the 
        # labelmap of Target that the segment is to be copied to: 
        #ADDLabmapIm = AddImages(ResSrcLabmapToCopyIm, TrgLabmapIms[ToRoiNum])
        
        # Create a new list of Target Labelmaps to include the new/modified 
        # ROI/labelmap as well as any others (unmodified):
        NewTrgLabmapIms = []
        
        for RoiNum in range(len(TrgLabmapIms)):
            if RoiNum == ToRoiNum:
                NewTrgLabmapIms.append(ORLabmapIm)
                #NewTrgLabmapIms.append(ADDLabmapIm) # <-- useful for visualisation but not correct for PixelData
            
            else:
                NewTrgLabmapIms.append(TrgLabmapIms[RoiNum])
                
                
        """
        Now need to convert NewTrgLabmapIms into a pixel array.
        """
        
        # Convert each labelmap image in NewTrgLabmapIms to a pixel array:
        NewTrgPixArrByRoi = []
        NewTrgPFFGStoDcmIndsByRoi = []
        NewTrgPFFGStoDcmInds = []
        NewTrgRoiNumsByRoi = []
        NewTrgRoiNums = []
        
        for RoiNum in range(len(NewTrgLabmapIms)):
            PixArr, FTSInds = Image2PixArr(LabmapImage=NewTrgLabmapIms[RoiNum])
            
            NewTrgPixArrByRoi.append(PixArr)
            
            NewTrgPFFGStoDcmIndsByRoi.append(FTSInds)
            
            NewTrgPFFGStoDcmInds.extend(FTSInds)
            
            NewTrgRoiNumsByRoi.append([RoiNum]*len(FTSInds))
            
            NewTrgRoiNums.extend([RoiNum]*len(FTSInds))
            
        
        # Combine all frames PixArrByRoi:
        NewTrgPixArr = np.zeros((len(NewTrgPFFGStoDcmInds), 
                                 PixArr.shape[1], 
                                 PixArr.shape[2]))
        
        n = 0
        
        for RoiNum in range(len(NewTrgLabmapIms)):
            for i in range(len(NewTrgPFFGStoDcmIndsByRoi[RoiNum])):
                NewTrgPixArr[n] = NewTrgPixArrByRoi[RoiNum][i]
                
                n += 1
        
        
        
        print(f'\nThere are {NewTrgPixArr.shape[0]} frames in NewTrgPixArr')
        print(f'NewTrgPFFGStoDcmIndsByRoi: {NewTrgPFFGStoDcmIndsByRoi}')
        print(f'NewTrgRoiNumsByRoi: {NewTrgRoiNumsByRoi}')
        
        
        
        
        # Create a list of frame numbers for NewTrg that indicate the frame
        # numbers that the masks originated from, and a similar list indicating
        # where they came from:
        #NewTrgFrameNums = []
        #NewTrgFrameSrcs = []
        
        # Create new Target SEG object and modify the tags accordingly:
        NewTrgSeg = ModifySegTagVals(TrgSeg, TrgDicoms, NewTrgPixArr, 
                                     NewTrgRoiNums, NewTrgPFFGStoDcmInds, 
                                     FromRoiNum, FromSliceNum, SrcSeg, 
                                     ToRoiNum, LogToConsole=False)
        
        
        
        
        if LogToConsole:
            print('\n\n')
            title = 'Results after resampling Source LabelMap image and ' \
                    + 'Source-to-copy LabelMap image:'
            ul = '*' * len(title)
            
            print('\n' + title)
            print(ul)
            
            print('\n*The Source and Target labelmap image size:*\n',
                          f'   Source:                   {SrcLabmapIm.GetSize()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetSize()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetSize()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSize()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSize()}')
                
            print('\n*The Source and Target labelmap image voxel spacings:*\n',
                          f'   Source:                   {SrcLabmapIm.GetSpacing()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetSpacing()} \n',
                          f'   Source-to-copy            {SrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetSpacing()} \n',
                          f'   Target:                   {TrgLabmapIm.GetSpacing()}')
                
            print('\n*The Source and Target labelmap image Origin:*\n',
                          f'   Source:                   {SrcLabmapIm.GetOrigin()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetOrigin()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetOrigin()} \n',
                          f'   Target:                   {TrgLabmapIm.GetOrigin()}')
                
            print('\n*The Source and Target labelmap image Direction:*\n',
                          f'   Source:                   {SrcLabmapIm.GetDirection()} \n',
                          f'   Resampled Source:         {ResSrcLabmapIm.GetDirection()} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyIm.GetDirection()} \n',
                          f'   Target:                   {TrgLabmapIm.GetDirection()}')
        
        
        
        if LogToConsole:
            SrcLabmapImMax = ImageMax(SrcLabmapIm)
            SrcLabmapImMin = ImageMin(SrcLabmapIm)
            
            ResSrcLabmapImMax = ImageMax(ResSrcLabmapIm)
            ResSrcLabmapImMin = ImageMin(ResSrcLabmapIm)
            
            SrcLabmapToCopyImMax = ImageMax(SrcLabmapToCopyIm)
            SrcLabmapToCopyImMin = ImageMin(SrcLabmapToCopyIm)
            
            ResSrcLabmapToCopyImMax = ImageMax(ResSrcLabmapToCopyIm)
            ResSrcLabmapToCopyImMin = ImageMin(ResSrcLabmapToCopyIm)
            
            TrgLabmapImMax = ImageMax(TrgLabmapIm)
            TrgLabmapImMin = ImageMin(TrgLabmapIm)
            
            print('\n*The Source and Target labelmap image minimum:*\n',
                          f'   Source:                   {SrcLabmapImMin} \n',
                          f'   Resampled Source:         {ResSrcLabmapImMin} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMin} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMin} \n',
                          f'   Target:                   {TrgLabmapImMin}')
            
            print('\n*The Source and Target labelmap image maximum:*\n',
                          f'   Source:                   {SrcLabmapImMax} \n',
                          f'   Resampled Source:         {ResSrcLabmapImMax} \n',
                          f'   Source-to-copy:           {SrcLabmapToCopyImMax} \n',
                          f'   Resampled Source-to-copy: {ResSrcLabmapToCopyImMax} \n',
                          f'   Target:                   {TrgLabmapImMax}')
        
        
        
        
        
        # Plot results:
        #PlotSitkImages_Src_ResSrc_Trg(SrcIm=SrcImage, 
        #                              ResSrcIm=ResSrcImage, 
        #                              TrgIm=TrgImage, 
        #                              LogToConsole=False)
        
        
        #return SrcImage, ResSrcImage, TrgImage
        #return SrcLabmapIm, ResSrcLabmapIm, TrgLabmapIm
        #return SrcLabmapIm, ResSrcLabmapIm, TrgLabmapIm,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm
        #return SrcImage, ResSrcImage, TrgImage,\
        #       SrcLabmapIm, ResSrcLabmapIm,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm, TrgLabmapIm
        #return SrcImage, ResSrcImage, TrgImage,\
        #       SrcLabmapIms, ResSrcLabmapIms,\
        #       SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
        #       TrgLabmapIms, NewTrgLabmapIms,\
        #       NewTrgSeg
         
         
        
        
    if CaseNum == '4':
        print(f'\nApplying Case 4..')
        
        
        
    if CaseNum == '5':
        print(f'\nApplying Case 5..')
        
        
    
        
        
    """
    The following was moved to ModifySegTagVals:
    # Check if NewTrgSeg exists (temporary check):
    if 'NewTrgSeg' in locals():
    
        # Generate a new SOP Instance UID:
        NewTrgSeg.SOPInstanceUID = pydicom.uid.generate_uid()
        
        # Generate a new Series Instance UID:
        NewTrgSeg.SeriesInstanceUID = pydicom.uid.generate_uid()
        
        
        return NewTrgSeg
    
    #else:
        #return []
        #return ResSrcImage
    """
        
        
        
    return SrcImage, ResSrcImage, TrgImage,\
           SrcLabmapIms, ResSrcLabmapIms,\
           SrcLabmapToCopyIm, ResSrcLabmapToCopyIm,\
           TrgLabmapIms, NewTrgLabmapIms,\
           NewTrgSeg