# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:48:35 2020

@author: ctorti
"""


"""
Function:
    RegisterSlices()
    
Purpose:
    Apply registration between two lists of DICOM objects.

Input:
    FixedFpaths   - list of full filepaths of fixed DICOM objects
    
    MovingFpaths  - list of full filepaths of DICOM objects to be registered to 
                    fixedDicoms
                        
    RegDicomDir   - directory to export registered DICOMs     
                        
    DataDict      - dictionary containing several variables (including
                    RegMethod)
    
    RegKey     - key in DataDict in which the registration filepaths are to
                    be added to (e.g. 'Reg' or 'RegPosCorr' or 'RegMapPosCorr')
                      

    
Returns:
    RegDicoms     - list of registered DICOM objects
    
    DataDict      - dictionary containing several variables
    
    

Note 1: 
    RegMethod defines the type of registration to perform on the DICOM scan
    (collection) that the ROI Collection is to be copied to.
    Acceptable values include:
    -- 'rigid' = perform rigid registration
    -- 'affine' = perform affine registration
    -- 'non-rigid' = perform non-rigid registration
    -- any other character string (including '') = do not perform any 
    registration
    
    
Note 2:
    Although the DICOMs are available it seems that SimpleITK needs to read 
    in the image directly from a DICOM file.
    
                  
"""



#def RegisterSlices(FixedFpaths, MovingFpaths, RegDicomDir, \
#                   DataDict, MovingKey, RegKey):
    
def RegisterSlices(DataDict, FixedKey, MovingKey, RegisteredKey):
    
    # Import packages and functions:
    
    import os
    import SimpleITK as sitk
    import copy
    import pydicom
    import time
    import numpy as np
    import importlib
    import ConvertNumpyDtype
    importlib.reload(ConvertNumpyDtype)
    from ConvertNumpyDtype import ConvertNumpyDtype
    from CreateDir import CreateDir
    
    
    # Get some necessary info from DataDict:
    RerunReg = DataDict['RerunReg']
    RegMethod = DataDict['RegMethod']
    SearchBy = DataDict['SearchBy']
    Debug = DataDict['Debug']
    #MappingInds = DataDict['Mapping']['MappingInds']
    
    # The key-dependent DICOM filepaths and SOP Instance UIDs:
    FixedFpaths = DataDict[FixedKey]['DicomFpaths']
    MovingFpaths = DataDict[MovingKey]['DicomFpaths']
    FixedSopUids = DataDict[FixedKey]['SopUids']
    MovingSopUids = DataDict[MovingKey]['SopUids']
    
    # Get the series numbers:
    FixedSeriesNo = DataDict[FixedKey]['SeriesNo']
    MovingSeriesNo = DataDict[MovingKey]['SeriesNo']
    
    # Get the series descriptions:
    FixedSeriesDesc = DataDict[FixedKey]['SeriesDesc']
    MovingSeriesDesc = DataDict[MovingKey]['SeriesDesc']
    
    # Get the slice numbers: 
    FixedSliceNos = DataDict[FixedKey]['SliceNos']
    MovingSliceNos = DataDict[MovingKey]['SliceNos']
    
    # If either the SourceDicoms or TargetDicoms were remapped, modify the 
    # appropriate list of filepaths:
    if 'Remapped' in [FixedKey, MovingKey]:
        # The key of the series that was remapped:
        RemappedKey = DataDict['RemappedKey']
    
        # Get the mapping indeces:
        if 'Source' in RemappedKey:
            MappingInds = DataDict['TargetToSourceMapping']['MappingInds']
            
            # Modify the list FixedFpaths:
            FixedFpaths = [FixedFpaths[ind] for ind in MappingInds]
        else:
            MappingInds = DataDict['SourceToTargetMapping']['MappingInds']
            
            # Modify the list MovingFpaths:
            MovingFpaths = [MovingFpaths[ind] for ind in MappingInds]
    
    
    # The Series Number allocated for the registered series:
    RegSeriesNo = DataDict[RegisteredKey]['SeriesNo']
    
    # And the directory where the registered DICOMs are to be saved:
    RegDicomDir = DataDict[RegisteredKey]['DicomDir']
    

    
    """
    The Study Instance UID and Frame of Reference UID of RegDicoms will be kept 
    unchanged from those of MovingDicoms, but a new Series Instance UID will 
    need to be generated.  
    
    Unique SOP Instance UIDs will be generated for each RegDicom.  
    
    I'm not sure if the Frame of Reference UID of RegDicoms should be changed
    from that of MovingDicoms.  For now I'll assume not.
    """
    
    
    # Create a dictionary to store the changes made to the DICOM of the  
    # registered slices:
    #changesDict = {}
    
    # Read in the first DICOM in MovingDicoms:
    MovingDicom = pydicom.read_file(MovingFpaths[0])    
    
    # Get the details of the pixel array of the first DICOM in MovingDicoms
    # (not necessary but useful for Debugging):
    MovingPixSpacing = MovingDicom.PixelSpacing
    #preRegRows = dicom2.Rows
    #preRegColumns = dicom2.Columns
    #MovingPixArray = MovingDicom.pixel_array
    MovingPixArrayShape = np.shape(MovingDicom.pixel_array)
    
    #print(f'\nPre-registered pixel array dims = {np.shape(MovingPixArray)}')
    print(f'\n\n--- Pre-registered pixel array dims = {MovingPixArrayShape}')
    #print(f'Pre-registered no of rows = {preRegRows}')
    #print(f'Pre-registered no of columns = {preRegColumns}')
    if Debug:
        print(f'\n---Pre-registered pixel spacing = {MovingPixSpacing}')
    
    # Get the Series Description of MovingDicom:
    #MovingSeriesDesc = MovingDicom.SeriesDescription
    
    # Create a new Series Description by appending RegMethod to the 
    # Series Description of MovingSeriesDesc:
    #RegSeriesDesc = MovingSeriesDesc + ' - ' + RegMethod + '-Reg'
    RegSeriesDesc = MovingSeriesDesc + ' ' + RegMethod + '-registered to ' \
                    + FixedSeriesDesc
    
    # Modify RegDicomDir to include the Series Description in parentheses:
    #RegDicomDir = RegDicomDir + ' (' + RegSeriesDesc + ')'
    
    # Modify RegDicomDir in DataDict:
    #DataDict[RegKey]['DicomDir'] = RegDicomDir
    
    # Create the directory:
    CreateDir(RegDicomDir, Debug)
    
    # Create a new Series Instance UID:
    RegSeriesUid = pydicom.uid.generate_uid()
    
    # Store the original and modified SOP Instance UIDs of the registered
    # sliced:
    #preRegSOPUIDs = []
    #postRegSOPUIDs = []
    
    # Initialise a list to store the filepaths and SOP Instance UIDs of the 
    # exported registered DICOMs:
    RegFpaths = []
    RegSopUids = []
    RegScanPos = []
    
    # Initialise a list to store the fixed, moving and registered DICOMs:
    FixedDicoms = []
    MovingDicoms = []
    RegDicoms = []
    
    # Keep track of the time it takes to execute each loop:
    times = []
    
    # times[0] = start time:
    times.append(time.time())
    
    F = len(MovingFpaths)
    
    # Iterate for each DICOM in MovingDicoms:
    for f in range(F):
        print(f'\n--- Performing registration on slice {f+1}/{F}..')
        
        # Get the filepaths:
        FixedFpath = FixedFpaths[f]
        MovingFpath = MovingFpaths[f]
        
        #print('\nMovingFpath =', MovingFpath)
        
        # Read in the Fixed and Moving DICOMs:
        FixedDicom = pydicom.read_file(FixedFpath) 
        MovingDicom = pydicom.read_file(MovingFpath) 
        
        # Store FixedDicom and MovingDicom in FixedDicoms and MovingDicoms:
        FixedDicoms.append(FixedDicom)
        MovingDicoms.append(MovingDicom)
        
        """
        Before doing anything serious (like registrations) check to see
        if the files already exist (there may be reasons to run this 
        function but without over-writting the already registered images).
        """
        
        # Get current filename:
        MovingFname = os.path.basename(MovingFpath)
        
        # Get filename (with extension):
        MovingFname = os.path.splitext(MovingFname)
        
        # Create a new filename, adding some info about the registration
        # method used:
        RegFname = MovingFname[0] + '-' + RegMethod + '-reg.dcm'
        
        if Debug:
            print('\n--- Changing the file name from:\n', \
              MovingFname[0], '\nto:\n', RegFname)
        
        # Create new filepath:
        RegFpath = os.path.join(RegDicomDir, RegFname)
        
        # Append RegFpath to RegFpaths:
        RegFpaths.append(RegFpath)
        
        
        
        # Only run the registration if the registered DICOM file doesn't 
        # already exist or if it exists and RerunReg=True:
        if os.path.isfile(RegFpath) and not RerunReg:
            if True:#Debug:
                print('\n--- The registered DICOM file already exists', \
                      'and will not be overwritten.')
                
                # Read in the DICOM:
                RegDicom = pydicom.read_file(RegFpath)
                
                # Append this DICOM object:
                RegDicoms.append(RegDicom)
            
        else:
            if not os.path.isfile(RegFpath):
                if Debug:
                    print('\n--- A registered DICOM file does not', \
                          'already exist.')
                    
            elif os.path.isfile(RegFpath) and RerunReg:
                if True:#Debug:
                    print('\n--- The registered DICOM file will be', \
                          'overwritten.')
            
        
            # Read in the 3D images using SimpleITK:
            FixedIm3D = sitk.ReadImage(FixedFpath)
            MovingIm3D = sitk.ReadImage(MovingFpath)
            
            """
            SimpleITK treats the images as 3D so need to extract the 2D
            images.
            """
            # Extract the 2D images from the 3D images:
            """
            Note: 
                
            When loading DICOMs with 512x512 pixel arrays the resulting
            2D images end up with 384x384 pixels for some unknown reason.
            
            If I replace (fixed3D.GetWidth(), fixed3D.GetHeight(), 0) 
            with (512, 512, 0) I get the following error:
                
            RuntimeError: Exception thrown in SimpleITK Extract: 
                    
            C:\SimpleElastix\build_20200108\ITK\Modules\Core\Common\src\itkDataObject.cxx:393:
                    
            Requested region is (at least partially) outside the largest
            possible region.
            
            So for now will have to accept that the image resolutions will
            be different.
            """
            FixedIm = sitk.Extract(FixedIm3D, \
                                   (FixedIm3D.GetWidth(), \
                                    FixedIm3D.GetHeight(), 0), \
                                    (0, 0, 0)
                                    )
                                   
            MovingIm = sitk.Extract(MovingIm3D, \
                                    (MovingIm3D.GetWidth(), \
                                     MovingIm3D.GetHeight(), 0), \
                                     (0, 0, 0)
                                     )
            
            
            # Perform registration depending on RegMethod:
            if RegMethod in ['rigid', 'affine']:
                # Perform rigid or affine registration:
                imFilter = sitk.ElastixImageFilter()
                imFilter.SetFixedImage(FixedIm)
                imFilter.SetMovingImage(MovingIm)
                imFilter.SetParameterMap(sitk.GetDefaultParameterMap(RegMethod))
                imFilter.LogToConsoleOn()
                imFilter.LogToFileOn()
                RegIm = imFilter.Execute()
            
            else:
                # Perform non-rigid registration:
                RegIm = sitk.Elastix(FixedIm, MovingIm)
                
                
            # Convert the SimpleITK image RegIm to a numpy data array:
            RegImNda = sitk.GetArrayFromImage(RegIm)
            
            # Get the data type:
            RegImDtype = str(RegImNda.dtype)
            
            """
            Since registered SimpleITK objects seem to return floats after
            converting to numpy, convert to the same data type of the pre-
            registered DICOM.
            
            April 9: 
                The code below forced the registered DICOM to have the same
                datatype as the fixed DICOM. This is probably the cause for the
                issues I sometimes get with pixel values that don't match the
                datatype of the DICOM, and the visually saturated pixel arrays.
                
                Since the registered DICOM starts off as a copy of the moving
                DICOM, it makes sense to force the registered DICOM to have the
                same datatype as the moving DICOM, rather than the fixed one!
            """
            
            # Convert the pre-registered image to numpty array:
            #FixedImNda = sitk.GetArrayFromImage(FixedIm) # April 9
            MovingImNda = sitk.GetArrayFromImage(MovingIm)
            
            # Get the data type:
            #FixedImDtype = str(FixedImNda.dtype) # April 9
            MovingImDtype = str(MovingImNda.dtype)
            
            """ 
            April 15: Converting the registered image to the same data type of
            MovingIm is ok if it is int16 or uint16 (since DICOMs need to be
            in either format). But I've already seen a case where the reported
            data type of MovingIm is float64 (might be a bug), so force the 
            conversion to be uint16:
            """
            ConvertTo = 'uint16'
            
            if True: #Debug:
                print('\n--- Data type of pre-registered image is ', \
                      MovingImDtype) # FixedImDtype) April 9
                
                print('\n--- Data type of registered image is ', \
                      RegImDtype) # April 15
                
                print('\n--- Data type of registered image will be', \
                      'converted to', ConvertTo) # April 15
                      
            
            # Convert the numpy data array from float to FixedImDtype:
            #RegImNdaConverted = ConvertNumpyDtype(RegImNda,\
            #                                      ConvertTo=FixedImDtype,\
            #                                      Debug=Debug)  April 9
            
            # Convert the numpy data array from float to MovingImDtype:
            #RegImNdaConverted = ConvertNumpyDtype(RegImNda,\
            #                                      ConvertTo=MovingImDtype,\
            #                                      Debug=Debug) # April 15
            
            RegImNdaConverted = ConvertNumpyDtype(RegImNda,\
                                                  ConvertTo=ConvertTo,\
                                                  Debug=Debug)
            
            
                
            
            
            # Read in the DICOM using pydicom:
            MovingDicom = pydicom.read_file(MovingFpath)
        
            # Create RegDicom, starting as a copy of MovingDicom (whose pixel 
            # array will be modified).
            RegDicom = copy.deepcopy(MovingDicom)
            
            # Change the pixel array to imRegNdaConverted:
            #dicomReg.pixel_array = imRegNda # <- AttributeError: can't set attribute
            #dicomReg['PixelArray'] = imRegNda # <-TypeError: Dataset contents must be DataElement instances
            #dicomReg.pixel_array = preRegPixelArray # <- AttributeError: can't set attribute
            #dicomReg.pixel_array = dicomReg.pixel_array # <- same error
            #dicomReg.PixelData = imRegNda.tostring()
            RegDicom.PixelData = RegImNdaConverted.tostring()
            
            """
            Update the other metadata: 
            """
            
            # Update the Image Position Patient and Image Orientation Patient:
            """ Note that since the SimpleITK object RegIm is a 2D image, the
            Origin (equivalent to Image Position Patient) only has two items 
            (whereas IPP has three).  Likewise, the Direction (equivalent to 
            Image Orientation Patient) only has four items (whereas IOP has
            six).  
            So modify the first two items in IPP to be the two items in the 
            Origin tuple, and modify the first, second, fourth and fifth items 
            in IOP to be the four items in the Orientation tuple. """
            RegDicom.ImagePositionPatient[0:2] = [str(item) for item in RegIm.GetOrigin()]
            
            print(f'\n--- Moving Image Position Patient =', MovingDicom.ImagePositionPatient)
            print(f'--- Registered Origin             = {RegIm.GetOrigin()}')
            
            # Store the scan position:
            if SearchBy == 'x':
                RegScanPos.append(float(RegDicom.ImagePositionPatient[0]))
            if SearchBy == 'y':
                RegScanPos.append(float(RegDicom.ImagePositionPatient[1]))
            if SearchBy == 'z':
                RegScanPos.append(float(RegDicom.ImagePositionPatient[2]))
            
            
            """ The line below doesn't work (TypeError: float() argument must 
            be a string or a number, not 'list'): """
            #RegDicom.ImageOrientationPatient[0:2,3:5] = [str(item) for item in RegIm.GetDirection()]
            
            # Create a list of indeces from 0 to 4:
            inds = list(range(5))
            
            # Remove the 3rd index:
            inds.pop(2)
            
            # Use inds to slice RegDicom.ImageOrientationPatient:
            #RegDicom.ImageOrientationPatient[inds] = [str(item) for item in RegIm.GetDirection()]
            
            """ The above line doesn't work either! So do it in a less elegant
            way: """
            RegDicom.ImageOrientationPatient[0] = str(RegIm.GetDirection()[0])
            RegDicom.ImageOrientationPatient[1] = str(RegIm.GetDirection()[1])
            RegDicom.ImageOrientationPatient[3] = str(RegIm.GetDirection()[2])
            RegDicom.ImageOrientationPatient[4] = str(RegIm.GetDirection()[3])
            
            print(f'\n--- Moving Image Orientation Patient =', MovingDicom.ImageOrientationPatient)
            print(f'--- Registered Direction             = {RegIm.GetDirection()}')
            
            
            # Update the pixel spacing:
            RegDicom.PixelSpacing = [str(item) for item in RegIm.GetSpacing()]
            
            print(f'\n--- Moving pixel spacing     = {MovingDicom.PixelSpacing}')
            print(f'--- Registered pixel spacing = {RegIm.GetSpacing()}')
            
            
            # Since the array dimensions have changed, update the Rows and 
            # Columns metadata:
            ##dicomReg.Rows, dicomReg.Columns = imRegNda.shape
            #RegPixArrayShape = RegImNdaConverted.shape
            #RegDicom.Rows, RegDicom.Columns = RegPixArrayShape
            #print(f'\n--- Registered pixel array dims = {RegPixArrayShape}')
            
            RegDicom.Rows = RegIm.GetWidth()
            RegDicom.Columns = RegIm.GetHeight()
            
            print(f'\n--- Moving pixel array dims     = ({MovingDicom.Rows}, {MovingDicom.Columns})')
            print(f'--- Registered pixel array dims = {RegIm.GetSize()}')
            
            
            # Get the original (pre-registered) SOP Instance UID:
            ##preRegSOPUID = dicomDict1_to_2[key]['Closest Dicom SOP UID']
            OrigSopUid = RegDicom.SOPInstanceUID
            
            ## Append to postRegSOPUIDs:
            #preRegSOPUIDs.append(preRegSOPUID)
            
            
            # Generate a new SOP Instance UID:
            """
            Note that there may be slices that are repeated in MovingDicoms 
            due to differences in slice thickness, e.g. slice no 1 in 
            FixedDicoms may be closest to slices 1 and 2 in MovingDicoms, so  
            the repeated slice will need to have a unique SOP Instance UID. 
            Hence each slice will need a unique SOP Instance UID to avoid 
            throwing up errors when it comes to uploading DICOMs to XNAT.
            """
            NewSopUid = pydicom.uid.generate_uid() 
            
            # Append to RegSopUids:
            RegSopUids.append(NewSopUid)
            
            if Debug:
                print('\n--- Changing the SOP Instance UID from:\n', \
                  OrigSopUid, '\nto:\n', NewSopUid)
            

        
            # Make changes to the tags:
            RegDicom.SOPInstanceUID = NewSopUid
            RegDicom.SeriesInstanceUID = RegSeriesUid
            RegDicom.SeriesDescription = RegSeriesDesc
            RegDicom.SeriesNumber = RegSeriesNo
            #RegDicom.StudyInstanceUID = postRegStudyUID
            #RegDicom.SeriesInstanceUID = postRegSeriesUID
            #RegDicom.FrameOfReferenceUID = postRegFrameOfRefUID
            
            
            if Debug:
                print('\n--- Exporting the registered DICOM to:\n', \
                      RegFpath)
            
            # Export the DICOM file:
            RegDicom.save_as(RegFpath)
            
            # Append RegDicom to RegDicoms:
            RegDicoms.append(RegDicom)

            
            # times[k] = this loop time:
            times.append(time.time())
            
            # Time taken to complete this loop:
            dt = times[-1] - times[-2]
            
            # Total time taken so far:
            Dt = times[-1] - times[0]
            
            # Average time taken per loop so far:
            period = Dt/(f+1)
            
            # Number of iterations remaining:
            rem = F - (f+1)
            
            # Anticipated time remaining to completion:
            t_rem = period*rem
            
            print(f'\n--- Took {round(dt, 1)} s for this registration.')
            if Dt < 60:
                print(f'\n--- Total time so far {round(Dt, 1)} s.')
            else:
                print(f'\n--- Total time so far {round(Dt/60, 1)} min.')
            if t_rem < 60:
                print(f'\n--- Expected time to completion {round(t_rem, 1)} s.')
            else:
                print(f'\n--- Expected time to completion {round(t_rem/60, 1)} min.')
            print('')
    
    
        
        """
        # Display changes:    
        print('\nChanging the Series Description from:\n', \
          preRegSeriesDesc, '\nto:\n', postRegSeriesDesc)
        print('\nChanging the Series Number from:\n', \
          preRegSeriesNo, '\nto:\n', postRegSeriesNo)
        print('\nChanging the Study Instance UID from:\n', \
          preRegStudyUID, '\nto:\n', postRegStudyUID)
        print('\nChanging the Series Instance UID from:\n', \
          preRegSeriesUID, '\nto:\n', postRegSeriesUID)
        print('\nChanging the Frame of Reference UID from:\n', \
          preRegFrameOfRefUID, '\nto:\n', postRegFrameOfRefUID)
        
        # Update the dictionary changesDict:
        changesDict.update({'preRegSeriesDescription':preRegSeriesDesc})
        changesDict.update({'postRegSeriesDescription':postRegSeriesDesc})
        changesDict.update({'preRegSeriesNumber':preRegSeriesNo})
        changesDict.update({'postRegSeriesNumber':postRegSeriesNo})
        changesDict.update({'preRegStudyInstanceUID':preRegStudyUID})
        changesDict.update({'postRegStudyInstanceUID':postRegStudyUID})
        changesDict.update({'preRegSeriesInstanceUID':preRegSeriesUID})
        changesDict.update({'postRegSeriesInstanceUID':postRegSeriesUID})
        changesDict.update({'preRegFrameOfReferenceUID':preRegFrameOfRefUID})
        changesDict.update({'postRegFrameOfReferenceUID':postRegFrameOfRefUID})
        changesDict.update({'preRegSOPInstanceUIDs':preRegSOPUIDs})
        changesDict.update({'postRegSOPInstanceUIDs':postRegSOPUIDs})
        """
        
    # Update the dictionary DataDict:
    DataDict['Fixed'] = {'DicomFpaths':FixedFpaths, \
                         'SopUids':FixedSopUids, \
                         'SeriesNo':FixedSeriesNo, \
                         'SeriesDesc':FixedSeriesDesc, \
                         'SliceNos':FixedSliceNos
                         }
    
    DataDict['Moving'] = {'DicomFpaths':MovingFpaths, \
                          'SopUids':MovingSopUids, \
                          'SeriesNo':FixedSeriesNo, \
                          'SeriesDesc':FixedSeriesDesc, \
                          'SliceNos':MovingSliceNos
                          }

    DataDict[RegisteredKey].update({'DicomFpaths':RegFpaths, \
                                    'SopUids':RegSopUids, \
                                    # Set SliceNos = MovingSliceNos:
                                    'SliceNos':MovingSliceNos
                                    })

    
    
    # April 1:  Stuff below not needed so commented out
    if False:
        # Get keys in DataDict[MovingKey]:
        keys = list(DataDict[MovingKey].keys())
        
        # If 'MappingInds' key exists in MovingKey, copy them over to RegKey,
        # since the registered slices will have the same indeces as the 
        # pre-registered series:
        if 'MappingInds' in keys:
            # Get the indeces:
            MappingInds = DataDict[MovingKey]['MappingInds']
            
            # Copy MappingInds to DataDict[RegisteredKey]:
            DataDict[RegisteredKey].update({'MappingInds':MappingInds})
        
        
    
    
    return FixedDicoms, MovingDicoms, RegDicoms, DataDict 

