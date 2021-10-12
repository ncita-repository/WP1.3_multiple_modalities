# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 10:26:04 2020

@author: ctorti
"""


"""
Note that RoiTags may be:

- a single attribute name at the base level, 

    e.g. 'StudyInstanceSOP' 

         = RoiObject['StudyInstanceSOP']

- a single attribute name at a nested level, 

    e.g. ['StructureSetROISequence', 0, 'ReferencedFrameOfReferenceUID']

         = RoiObject['StructureSetROISequence'][0]['ReferencedFrameOfReferenceUID']

         Note that in this case the 0 indexes a single-item attribute.

- an attribute that needs to be changed N number of times for an attribute with N items,

    e.g. ['ReferencedFrameOfReferenceSequence', 0, 'RTReferencedStudySequence', 0, \
          'RTReferencedSeriesSequence', 0, 'ContourImageSequence', 1, 'ReferencedSOPInstanceUID']   

          = RoiObject['ReferencedFrameOfReferenceSequence'][0]['RTReferencedStudySequence'][0]\
            ['RTReferencedSeriesSequence'][0]['ContourImageSequence'][N]['ReferencedSOPInstanceUID'] 

            Note that the index [1] will denote an attribute that has more than 1 item.

To accommodate the fact that there may be one or more attribute names in RoiTags, each RoiTags 
must be itself a list (even if the list includes only one item, e.g. ['StudyInstanceUID']).
"""

def WalkThroughAttributes(DicomObject, RoiObject, TagPair, Debug):
    # Import packages:
    import copy
    
    # Get the DicomTag and RoiTags:
    DicomTag = copy.deepcopy(TagPair[0])
    RoiTags = copy.deepcopy(TagPair[1])
    
    if Debug:
        print('\nTagPair = ', TagPair)
        print('\nDicomTag = ', DicomTag)
        print('\nRoiTags = ', RoiTags)
      
    
    """ 
    Note:  Since some attributes are optional, it's possible that DicomTag
    doesn't exist in DicomObject, so use try/except.
    """
    try:
        # The DICOM value is:
        DicomVal = DicomObject[DicomTag].value

        # Use pop to get the first item in RoiTags:
        RoiTag = RoiTags.pop(0)
    
        if Debug:
            print('\nPopped item (=RoiTag) = ', RoiTag)
            
        # Update TagPair with the remaining items in RoiTags:
        TagPair[1] = copy.deepcopy(RoiTags)
    
        """
        If what remains in RoiTags after pop'ing the first item has a length greater 
        than 1, the desired attribute is nested, in which case, and the next item in 
        RoiTags is the index that the first item is to be indexed with.  If the 
        remaining RoiTags has length equal to 1, the attribute at this level has only
        1 item, so no indexing is required.  If the remaining RoiTags has length equal
        to 0, the attribute is at the final nested level where the desired attribute is.
        """
        if len(RoiTags) > 1:
            if Debug:
                print('\nlen(RoiTags) > 1, so next level has more than 1 item..')
            
            # Pop again to get the tag index:
            TagInd = RoiTags.pop(0)
            
            # Update TagPair with the remaining items in RoiTags:
            TagPair[1] = copy.deepcopy(RoiTags)
            
            if Debug:
                print('\nPopped item (=TagInd) = ', TagInd)
                print('\nAfter pop (=RoiTags) =', RoiTags)
                print('\nRoiObject[RoiTag=' + RoiTag + '][TagInd =' \
                      + str(TagInd) + '] = ', RoiObject[RoiTag][TagInd])
            
            return WalkThroughAttributes(DicomObject, \
                                         RoiObject[RoiTag][TagInd], \
                                         TagPair, \
                                         Debug)
        
        elif len(RoiTags) == 1: # so at next level has only 1 item
            if Debug:
                print('\nlen(RoiTags) = 1, so next level has only 1 item..')
                print('\nPopped item (=TagInd) = ', TagInd)
                print('\nAfter pop (=RoiTags) =', RoiTags)
                print('\nRoiObject[RoiTag=' + RoiTag + '] = ', \
                      RoiObject[RoiTag][TagInd])
            
            return WalkThroughAttributes(DicomObject, \
                                         RoiObject[RoiTag], \
                                         TagPair, \
                                         Debug)
            
            
        else: # len(RoiTags) = 0, so at highest nested level
            if Debug:
                print('\nlen(RoiTags) = 0, so at final nested level..')     
                print('\nRoiObject[RoiTag=' + RoiTag + '] = ', \
                      RoiObject[RoiTag])
                
            # The original ROI value is:
            OldRoiVal = copy.deepcopy(RoiObject[RoiTag].value)
            
            # The DICOM value is:
            #DicomVal = DicomObject[DicomTag].value
            
            # No changes to the ROI object need to be made if the ROI and DICOM values 
            # are the same:
            if OldRoiVal==DicomVal:
                if Debug:
                    print('\nThe value of the DICOM and ROI attribute are the same.')
                    
                return OldRoiVal, OldRoiVal, RoiObject
                   
            else:
                if Debug:
                    print('\nThe value of the DICOM attribute that will ' \
                          + 'replace the ROI attribute value = ', DicomVal)
    
                # Modify the attribute value of RoiObject to the attribute 
                # value of DicomObject:
                RoiObject[RoiTag].value = copy.deepcopy(DicomObject[DicomTag].value)
    
                # The new ROI value is:
                NewRoiVal = copy.deepcopy(RoiObject[RoiTag].value)
    
                #return oldValue, newValue, RoiObject[RoiTag]
                return OldRoiVal, NewRoiVal, RoiObject
            
    except ValueError:
        if Debug:
            print('\n' + DicomTag + ' does not exist in the DICOM object.')