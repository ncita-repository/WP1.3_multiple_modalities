# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 19:04:35 2020

@author: ctorti
"""



"""
March 19:
I was trying to figure out how to make changes to the RTSTRUCT's metadata in
a more compact way but I didn't manage to figure out how to make the change...
The following is the relevant snipped from ModifyRoi():
"""



    # Create a dictionary to store the changes made to the ROI that will 
    # later be added to xnatDict:
    changesDict = {}
    
    # Make changes to the metadata in roi2 to reflect the corresponding 
    # metadata contained in the DICOM. 
    
    # The list of metadata to be altered with are the following in pairs with
    # the DICOM tag followed by the DICOM-RTSTRUCT tag in each pair:
    tagPairs = [
                ['StudyInstanceUID', 'StudyInstanceUID'],\
                
                ['StudyInstanceUID', ['ReferencedFrameOfReferenceSequence',\
                                      'RTReferencedStudySequence',\
                                      'ReferencedSOPInstanceUID'
                                      ]],\
    
                ['SeriesInstanceUID', ['ReferencedFrameOfReferenceSequence',\
                                      'RTReferencedStudySequence',\
                                       'RTReferencedSeriesSequence',\
                                       'SeriesInstanceUID'
                                       ]],\
    
                ['FrameOfReferenceUID', 'FrameOfReferenceUID'],\
    
                ['FrameOfReferenceUID', ['StructureSetROISequence',\
                                         'ReferencedFrameOfReferenceUID'
                                        ]],\
                                         
                ['StudyDate', 'StudyDate'],\
                
                ['StudyTime', 'StudyTime'],\
                
                ['PatientName', 'PatientName'],\
                
                ['PatientID', 'PatientID'],\
                
                ['PatientSex', 'PatientSex']
               ]
    
    # Make copy of roi2, roi2copy:
    roi2copy = copy.deepcopy(roi2)
    
    # Loop for each tagPair in tagPairs:
    for tagPair in tagPairs:
        dicomTag = tagPair[0]
        roiTags = tagPair[1] 
        
        """
        Note that roiTags may be a list with a single tag or a list of tags
        that point to a nested metadata.
        """
        
        # Getting the DICOM tag value is easy (if it exists). It may not exist
        # so use try/exception:
        try:
            dicomVal = dicom2[dicomTag].value
            
            # Getting the DICOM-RTSTRUCT tag value requires dealing with tags
            # that are nested:
            """
            e.g. if roiTags = ['ReferencedFrameOfReferenceSequence',\
                               'RTReferencedStudySequence',\
                               'ReferencedSOPInstanceUID']
            then the ultimate tag value would be:
                roi['ReferencedFrameOfReferenceSequence'][0]
                .['RTReferencedStudySequence'][0]
                .['ReferencedSOPInstanceUID'].value
            """
            while len(roiTags) > 1:
                roi2copy = roi2copy[roiTags.pop(0)][0]
      
            # The ROI tag value is:
            """
            At this point roiTags will be a single-item list, so need to index
            the 0th item.
            """
            roiVal = roi2copy[roiTags[0]].value
            
            # Finally, compare dicomVal with roiVal:
            if dicomVal==roiVal:
                print('\n' + dicomTag + ' are the same:\n',\
                      roi2copy[roiTags[0]].value)
            else:
                print('\nChanging '+ dicomTag + ' from:\n',\
                      roi2copy[roiTags[0]].value, \
                      '\nto:\n', dicom2[dicomTag].value)
        
                # Update the dictionary changesDict:
                changesDict.update({'from' + dicomTag : roi2copy[roiTags[0]].value})
                
                # Modify the RTSTUCT's metadata:
                """
                This works easily for tags at the first tier but need to figure
                out how to do this for nested tags. 
                THE NEXT TWO LINES ARE INCORRECT!
                """
                roi2cpy[???].value = copy.deepcopy(dicom2[dicomTag].value)
                
                # Update the dictionary changesDict:
                changesDict.update({ 'to' + dicomTag : roi2copy[???].value})
            
            
        except AttributeError:
            print('\n' + dicomTag + ' does not exist in the DICOMs')
    
        
        