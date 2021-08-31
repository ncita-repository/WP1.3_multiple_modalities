# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:38:51 2021

@author: ctorti
"""


""" Helper functions for creating DROs. """

import time
from copy import deepcopy

def modify_ref_im_seq(seq, dicoms, refAllSOPs=False, p2c=False):
    """
    Modify a ReferencedImageSequence (or ReferencedInstanceSequence).  
    
    Parameters
    ----------
    seq : Sequence of a Pydicom Object
        A ReferencedImageSequence to be modified.
    SeqNum : int
        The (zero-indexed) integer of the sequence to modify.
    dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    refAllSOPs : bool, optional
        If True, all SOPInstanceUIDs in dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced. The default value is False.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False.
    
    Returns
    -------
    seq : Sequence of a Pydicom Object
        The modified ReferencedImageSequence.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
       
    # The number of sequences:
    N0 = len(seq)
    
    # The required number of sequences:
    if refAllSOPs:
        N1 = len(dicoms)
        
        if p2c:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence,',
                  f'and {N1} (= no. of DICOMs) required sequences.')
    else:
        N1 = 1
        
        if p2c:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence',
                  f'and {N1} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if N1 > N0:
        for i in range(N1 - N0):           
            seq.append(seq[-1])          
    else:
        for i in range(N0 - N1):
            seq.pop()
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n   *Took {dtime} s to add sequences to ReferencedImageSequence.')
        
    # Update the sequences:
    for i in range(N1):
        seq[i].ReferencedSOPClassUID = dicoms[i].SOPClassUID
           
        seq[i].ReferencedSOPInstanceUID = dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n   *Took {dtime} s to update ReferencedImageSequence.')
    
    return seq

def modify_ref_ser_seq(dro, seqNum, dicoms, refAllSOPs=False, p2c=False):
    """
    Modify the ReferencedSeriesSequence in a DICOM Spatial or Deformable
    Spatial Registration Object (DRO) for a chosen sequence number.  
    
    Parameters
    ----------
    dro : Pydicom Object
        The DICOM Registration Object to be modified.
    seqNum : int
        The (zero-indexed) integer of the sequence to modify.
    dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    refAllSOPs : bool, optional
        If True, all SOPInstanceUIDs in dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced. The default value is False.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False.
    
    Returns
    -------
    dro : Pydicom Object
        The modified DICOM Registration Object.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
    
    dro.ReferencedSeriesSequence[seqNum]\
       .SeriesInstanceUID = dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    origRIS = len(dro.ReferencedSeriesSequence[seqNum]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if refAllSOPs:
        reqRIS = len(dicoms)
        
        if p2c:
            print(f'\nThere are {origRIS} sequences in',
                  f'ReferencedSeriesSequence[{seqNum}].ReferencedInstanceSequence,',
                  f'and {reqRIS} (= no. of DICOMs) required sequences.')
    else:
        reqRIS = 1
        
        if p2c:
            print(f'\nThere are {origRIS} sequences in',
                  f'ReferencedSeriesSequence[{seqNum}].ReferencedInstanceSequence,',
                  f'and {reqRIS} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if reqRIS > origRIS:
        for i in range(reqRIS - origRIS):
            lastItem = deepcopy(dro.ReferencedSeriesSequence[seqNum]\
                                   .ReferencedInstanceSequence[-1])
            
            dro.ReferencedSeriesSequence[seqNum]\
               .ReferencedInstanceSequence.append(lastItem)          
    else:
        for i in range(origRIS - reqRIS):
            dro.ReferencedSeriesSequence[seqNum]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n   *Took {dtime} s to add sequences to',
              f'ReferencedSeriesSequence[{seqNum}].ReferencedInstanceSequence.')
        
    # Update the sequences:
    for i in range(reqRIS):
        dro.ReferencedSeriesSequence[seqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = dicoms[i].SOPClassUID
           
        dro.ReferencedSeriesSequence[seqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        print(f'\n   *Took {dtime} s to update',
              f'ReferencedSeriesSequence[{seqNum}].ReferencedInstanceSequence.')
    
    return dro

def modify_stu_con_oth_ref_ins_seq(
        dro, seqNum, dicoms, refAllSOPs=False, p2c=False
        ):
    """
    Modify the StudiesContainingOtherReferencedInstancesSequence in a DICOM 
    Spatial or Deformable Spatial Registration Object (DRO) for a chosen 
    sequence number.  
    
    Parameters
    ----------
    dro : Pydicom Object
        The DICOM Registration Object to be modified.
    seqNum : int
        The (zero-indexed) integer of the sequence to modify.
    dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    refAllSOPs : bool, optional
        If True, all SOPInstanceUIDs in dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced. The default value is False.
    p2c : bool, optional
        Denotes whether intermediate results will be logged to the console.
        The default value is False.
    
    Returns
    -------
    dro : Pydicom Object
        The modified DICOM Registration Object.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
    
    dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
       .StudyInstanceUID = dicoms[0].StudyInstanceUID
    
    dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    origRIS = len(
        dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
            .ReferencedSeriesSequence[0]\
            .ReferencedInstanceSequence
            )
    
    # The required number of sequences:
    if refAllSOPs:
        reqRIS = len(dicoms)
        
        if p2c:
            msg = f'\nThere are {origRIS} sequences in '\
                + f'StudiesContainingOtherReferencedInstancesSequence[{seqNum}]'\
                + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                + f' and {reqRIS} (= no. of DICOMs) required sequences.'
            print(msg)
    else:
        reqRIS = 1
        
        if p2c:
            msg = f'\nThere are {origRIS} sequences in '\
                + f'StudiesContainingOtherReferencedInstancesSequence[{seqNum}]'\
                + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                + f' and {reqRIS} (the first DICOM) required sequence.'
            print(msg)
    
    # Add/remove sequences:
    if reqRIS > origRIS:
        for i in range(reqRIS - origRIS):
            lastItem = deepcopy(
                dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
                    .ReferencedSeriesSequence[0]\
                    .ReferencedInstanceSequence[-1]
                    )
            
            dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(lastItem)          
    else:
        for i in range(origRIS - reqRIS):
            dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        msg = f'\n   *Took {dtime} s to add sequences to '\
            + f'StudiesContainingOtherReferencedInstancesSequence[{seqNum}]'\
            + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
        
    # Update the sequences:
    for i in range(reqRIS):
        dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = dicoms[i].SOPClassUID
           
        dro.StudiesContainingOtherReferencedInstancesSequence[seqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    dtime = round(times[-1] - times[-2], 3)
    if p2c:
        msg = f'\n   *Took {dtime} s to update '\
            + f'StudiesContainingOtherReferencedInstancesSequence[{seqNum}]'\
            + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
    
    return dro