# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 20:38:51 2021

@author: ctorti
"""


""" Helper functions for creating DROs. """

import time
from copy import deepcopy

def modify_ref_im_seq(Sequence, Dicoms, RefAllSOPs=False, LogToConsole=False):
    """
    Modify a ReferencedImageSequence (or ReferencedInstanceSequence).  
    
    Parameters
    ----------
    Sequence : Sequence of a Pydicom Object
        A ReferencedImageSequence to be modified.
    SeqNum : int
        The (zero-indexed) integer of the sequence to modify.
    Dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    RefAllSOPs : bool, optional (False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    Sequence : Sequence of a Pydicom Object
        The modified ReferencedImageSequence.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
       
    # The number of sequences:
    N0 = len(Sequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        N1 = len(Dicoms)
        
        if LogToConsole:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence,',
                  f'and {N1} (= no. of DICOMs) required sequences.')
    else:
        N1 = 1
        
        if LogToConsole:
            print(f'\nThere are {N0} sequences in ReferencedImageSequence',
                  f'and {N1} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if N1 > N0:
        for i in range(N1 - N0):           
            Sequence.append(Sequence[-1])          
    else:
        for i in range(N0 - N1):
            Sequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to ReferencedImageSequence.')
        
    # Update the sequences:
    for i in range(N1):
        Sequence[i].ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Sequence[i].ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update ReferencedImageSequence.')
    
    return Sequence

def modify_ref_ser_seq(Dro, SeqNum, Dicoms, RefAllSOPs=False, 
                       LogToConsole=False):
    """
    Modify the ReferencedSeriesSequence in a DICOM Spatial or Deformable
    Spatial Registration Object (DRO) for a chosen sequence number.  
    
    Parameters
    ----------
    Dro : Pydicom Object
        The DICOM Registration Object to be modified.
    SeqNum : int
        The (zero-indexed) integer of the sequence to modify.
    Dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    RefAllSOPs : bool, optional (False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    Dro : Pydicom Object
        The modified DICOM Registration Object.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
    
    Dro.ReferencedSeriesSequence[SeqNum]\
       .SeriesInstanceUID = Dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.ReferencedSeriesSequence[SeqNum]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(Dicoms)
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (= no. of DICOMs) required sequences.')
    else:
        ReqRIS = 1
        
        if LogToConsole:
            print(f'\nThere are {OrigRIS} sequences in',
                  f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence,',
                  f'and {ReqRIS} (the first DICOM) required sequence.')
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.ReferencedSeriesSequence[SeqNum]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.ReferencedSeriesSequence[SeqNum]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.ReferencedSeriesSequence[SeqNum]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to add sequences to',
              f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence.')
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.ReferencedSeriesSequence[SeqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Dro.ReferencedSeriesSequence[SeqNum]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        print(f'\n   *Took {Dtime} s to update',
              f'ReferencedSeriesSequence[{SeqNum}].ReferencedInstanceSequence.')
    
    return Dro

def modify_stu_con_oth_ref_ins_seq(Dro, SeqNum, Dicoms, RefAllSOPs=False,
                                   LogToConsole=False):
    """
    Modify the StudiesContainingOtherReferencedInstancesSequence in a DICOM 
    Spatial or Deformable Spatial Registration Object (DRO) for a chosen 
    sequence number.  
    
    Parameters
    ----------
    Dro : Pydicom Object
        The DICOM Registration Object to be modified.
    SeqNum : int
        The (zero-indexed) integer of the sequence to modify.
    Dicoms : list of Pydicom Objects
        The DICOMs that related to the sequence.
    RefAllSOPs : bool, optional (False by default)
        If True, all SOPInstanceUIDs in Dicoms will be referenced in the
        sequence.  If False, only the SOPInstanceUID of the first DICOM will be
        referenced.
    LogToConsole : bool, optional (False by default)
        Denotes whether intermediate results will be logged to the console.
    
    Returns
    -------
    Dro : Pydicom Object
        The modified DICOM Registration Object.
    """
    
    # Start timing:
    times = []
    times.append(time.time())
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
       .StudyInstanceUID = Dicoms[0].StudyInstanceUID
    
    Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
       .ReferencedSeriesSequence[0]\
       .SeriesInstanceUID = Dicoms[0].SeriesInstanceUID
       
    # The number of sequences:
    OrigRIS = len(Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
                     .ReferencedSeriesSequence[0]\
                     .ReferencedInstanceSequence)
    
    # The required number of sequences:
    if RefAllSOPs:
        ReqRIS = len(Dicoms)
        
        if LogToConsole:
            msg = f'\nThere are {OrigRIS} sequences in '\
                  + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
                  + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                  + f' and {ReqRIS} (= no. of DICOMs) required sequences.'
            print(msg)
    else:
        ReqRIS = 1
        
        if LogToConsole:
            msg = f'\nThere are {OrigRIS} sequences in '\
                  + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
                  + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence,'\
                  + f' and {ReqRIS} (the first DICOM) required sequence.'
            print(msg)
    
    # Add/remove sequences:
    if ReqRIS > OrigRIS:
        for i in range(ReqRIS - OrigRIS):
            LastItem = deepcopy(Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
                                   .ReferencedSeriesSequence[0]\
                                   .ReferencedInstanceSequence[-1])
            
            Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.append(LastItem)          
    else:
        for i in range(OrigRIS - ReqRIS):
            Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
               .ReferencedSeriesSequence[0]\
               .ReferencedInstanceSequence.pop()
    
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        msg = f'\n   *Took {Dtime} s to add sequences to '\
              + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
              + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
        
    # Update the sequences:
    for i in range(ReqRIS):
        Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPClassUID = Dicoms[i].SOPClassUID
           
        Dro.StudiesContainingOtherReferencedInstancesSequence[SeqNum]\
           .ReferencedSeriesSequence[0]\
           .ReferencedInstanceSequence[i]\
           .ReferencedSOPInstanceUID = Dicoms[i].SOPInstanceUID
                      
    times.append(time.time())
    Dtime = round(times[-1] - times[-2], 3)
    if LogToConsole:
        msg = f'\n   *Took {Dtime} s to update '\
              + f'StudiesContainingOtherReferencedInstancesSequence[{SeqNum}]'\
              + '.ReferencedSeriesSequence[0].ReferencedInstanceSequence.'
        print(msg)
    
    return Dro