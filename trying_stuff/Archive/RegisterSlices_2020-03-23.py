# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:50:25 2020

@author: ctorti
"""


"""
Function:
    RegisterSlices()
    
Purpose:
    Apply registration to slices from a dictionary resulting from calling the
    function GetClosestSlices().

Input:
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- SOP Instance UID
                      -- filepath
                      -- Image Position
                      -- Image Orientation
                      -- Slice Location
                      -- etc.
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series
                        
    regMethod       - Type of registration to perform on the DICOM scan
                      (collection) that the ROI Collection is to be copied
                      to.  Acceptable values include:
                        -- 'rigid' = perform rigid registration
                        -- 'affine' = perform affine registration
                        -- 'non-rigid' = perform non-rigid registration
                        -- any other character string (including '') = do not
                            perform any registration
                        
    xnatDict        - Dictionary containing XNAT login details, REST path 
                      variables and directory paths
                      (see DownloadFilesForRoiCopy() for details)
                      
    debug           - Prints some useful results if True

    
Returns:
    dicomDict1_to_2 - Dictionary dicomDict1 with added entries:   
                      -- Registration method used
                      -- Registered closest Dicom SOP Instance UID
                      -- Registered closest Dicom filepath
                      -- Registered closest pixel array 
                      -- Registered closest Dicom 
                      of the slice in the 2nd series that is closest to each
                      slice in the 1st series 
    
    xnatDict        - Dictionary containing XNAT login details, REST path 
                      variables and directory paths
                      (see DownloadFilesForRoiCopy() for details)
                  
"""



#def RegisterSlices(dicomDict1_to_2, regMethod, xnatDict, debug):
def RegisterSlices(dicomDict1_to_2, dataDict):
    
    # Import packages and functions:
    
    #import importlib
    #import DicomHelperFuncs
    #importlib.reload(DicomHelperFuncs)
    import os
    import SimpleITK as sitk
    import copy
    import pydicom
    import time
    import numpy as np
    #from DicomHelperFuncs import ConvertNumpyDtype
    from ConvertNumpyDtype import ConvertNumpyDtype
    
    
    # Get some necessary info from dataDict:
    regMethod = dataDict['RegMethod']
    debug = dataDict['Debug']
    
    
    #print(regMethod)
    
    # If regMethod is not 'rigid', 'affine' or 'non-rigid', simply return an
    # empty array for imReg:
    """ This probably isn't necessary as why would the function be called if
    regMethod is ''?... """
    if not regMethod in ['rigid', 'affine', 'non-rigid']:
        print('Not performing registration')
        return [], dataDict
    
    else: 
        # Get the directory to export the registered DICOMs:
        toRegDicomDir = dataDict['To']['RegDicomDir']
        
        ## Create the directory to store the registered DICOMs if they don't 
        ## exist:
        #try:
        #    os.mkdir(toRegDicomDir)
        #    if debug:
        #        print('Directory', toRegDicomDir, 'created')
        #except FileExistsError:
        #    if debug:
        #        print('Directory', toRegDicomDir, 'already exists')
                
        
        # Create a dictionary to store the changes made to the DICOM of the  
        # registered slices:
        changesDict = {}
        
        #print('\ndicomDict1_to_2:\n', dicomDict1_to_2, '\n\n')
        
        # Get the keys (SOP UIDs) in dicomDict1_to_2:
        keys = list(dicomDict1_to_2.keys())
        K = len(keys)
        
        
        """
        The new DICOMS will need some new (unique) tags. New SOP Instance UIDs
        will be generated for each key in dicomDict1_to_2, in the for loop 
        below.  Other tags only need generating once, so will use the values 
        from the first key (keys[0]).
        """
        # Load the first DICOM in the set to be registered:
        dicom2 = dicomDict1_to_2[keys[0]]['Closest Dicom']
        
        # Add some details of the type of registration used in the Series
        # Description:
        preRegSeriesDesc = dicom2.SeriesDescription
        postRegSeriesDesc = preRegSeriesDesc + ' ' + regMethod + ' registered'
        
        # Add some details of the type of registration used in the Series
        # Number:
        preRegSeriesNo = dicom2.SeriesNumber
        postRegSeriesNo = str(preRegSeriesNo) + '-' + regMethod + 'Reg'
        """
        Maybe Series Number cannot have strings, so just change the number.
        """
        postRegSeriesNo = str(preRegSeriesNo) + '-2'
        """
        That was a problem too.
        """
        postRegSeriesNo = preRegSeriesNo + 10
        
        # Create a different series number for different registration methods
        # (so that multiple different registrations can be uploaded to XNAT):
        if 'rigid' in regMethod:
            postRegSeriesNo = preRegSeriesNo + 100
        elif 'affine' in regMethod:
            postRegSeriesNo = preRegSeriesNo + 200
        else:
            postRegSeriesNo = preRegSeriesNo + 300
        
        # Create a new Study Instance UID:
        preRegStudyUID = dicom2.StudyInstanceUID
        postRegStudyUID = pydicom.uid.generate_uid()   
        
        # Create a new Series Instance UID:
        preRegSeriesUID = dicom2.SeriesInstanceUID
        postRegSeriesUID = pydicom.uid.generate_uid()  
        
        # Create a new Frame of Reference UID:
        preRegFrameOfRefUID = dicom2.FrameOfReferenceUID
        postRegFrameOfRefUID = pydicom.uid.generate_uid() 
        
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
            
        
        # Get the details of the pixel array of the first DICOM from the set 
        # to be registered (not necessary but useful for debugging):
        preRegPixelSpacing = dicom2.PixelSpacing
        #preRegRows = dicom2.Rows
        #preRegColumns = dicom2.Columns
        preRegPixelArray = dicom2.pixel_array
        
        print(f'\nPre-registered pixel array dims = {np.shape(preRegPixelArray)}')
        #print(f'Pre-registered no of rows = {preRegRows}')
        #print(f'Pre-registered no of columns = {preRegColumns}')
        print(f'\nPre-registered pixel spacing = {preRegPixelSpacing}')
        
        
        # Store the original and modified SOP Instance UIDs of the registered
        # sliced:
        preRegSOPUIDs = []
        postRegSOPUIDs = []
        
        # Keep track of the time it takes to execute each loop:
        times = []
        
        # times[0] = start time:
        times.append(time.time())
        
        # Iterate for each key:
        for k in range(K):
            key = keys[k]
            
            print('\nkey =', key)
            
            print(f'\nPerforming registration on slice {k+1}/{K}..')
            
            """
            Although the pixel arrays are already contained within 
            dicomDict1_to_2, it seems that SimpleITK needs to read and image 
            directly from a DICOM file.
            """
            # So get the filepaths:
            fpath1 = dicomDict1_to_2[key]['Dicom fpath']
            fpath2 = dicomDict1_to_2[key]['Closest Dicom fpath']
            
            print('\nfpath2 =', fpath2)
            
            """
            Before doing anything serious (like registrations) check to see
            if the files already exist (there may be reasons to run this 
            function but without over-writting the already registered images).
            """
            
            # Get current filename:
            preRegFname = os.path.basename(fpath2)
            
            # Get filename (with extension):
            preRegFname = os.path.splitext(preRegFname)
            
            # Create a new filename, adding some info about the registration
            # method used:
            postRegFname = preRegFname[0] + '-' + regMethod + '-registered.dcm'
            
            print('\n   Changing the file name from:\n', \
              preRegFname[0], '\nto:\n', postRegFname)
            
            # Create new filepath:
            postRegFpath = os.path.join(toRegDicomDir, postRegFname)
            
            
            """
            # Decide whether to over-write the registered files (if they 
            # already exist) or not:
            overWrite = False
            
            # Check if the registered DICOM file already exists:
            if os.path.isfile(postRegFpath):
                print('\nThe registered DICOM file already exists.')
                
                if not overWrite:
                    print('\nThe registered DICOM file will not be overwritten.')
                else:
                    print('\nThe registered DICOM file will be overwritten.')
            """
        
            
            # Read in the 3D images using SimpleITK:
            im3D1 = sitk.ReadImage(fpath1)
            im3D2 = sitk.ReadImage(fpath2)
            
            """
            SimpleITK treats the images as 3D so need to extract the 2D images.
            """
            # Extract the 2D images from the 3D images:
            """
            Note: When loading DICOMs with 512x512 pixel arrays the resulting
            2D images end up with 384x384 pixels for some unknown reason.
            If I replace (im3D1.GetWidth(), im3D1.GetHeight(), 0) with
            (512, 512, 0) I get the following error:
                RuntimeError: Exception thrown in SimpleITK Extract: 
                    C:\SimpleElastix\build_20200108\ITK\Modules\Core\Common\src\itkDataObject.cxx:393:
                    Requested region is (at least partially) outside the 
                    largest possible region.
            So for now will have to accept that the image resolutions will be
            different.
            """
            im1 = sitk.Extract(im3D1, (im3D1.GetWidth(), im3D1.GetHeight(), 0), (0, 0, 0))
            im2 = sitk.Extract(im3D2, (im3D2.GetWidth(), im3D2.GetHeight(), 0), (0, 0, 0))
            
            
            # Perform registration depending on regMethod:
            if regMethod in ['rigid', 'affine']:
                # Perform rigid or affine registration:
                imFilter = sitk.ElastixImageFilter()
                imFilter.SetFixedImage(im1)
                imFilter.SetMovingImage(im2)
                imFilter.SetParameterMap(sitk.GetDefaultParameterMap(regMethod))
                imReg = imFilter.Execute()
            
            else:
                # Perform non-rigid registration:
                imReg = sitk.Elastix(im1, im2)
                
                
            # Convert the SimpleITK image imReg to a numpy data array:
            imRegNda = sitk.GetArrayFromImage(imReg)
            
            # Since registered SimpleITK objects seem to return floats after
            # converting to numpy, convert to the same data type of the pre-
            # registered DICOM.
            
            # Convert the pre-registered image to numpty array:
            im2Nda = sitk.GetArrayFromImage(im2)
            
            # Get the data type:
            im2dtype = str(im2Nda.dtype)
            
            if debug:
                print('\nData type of pre-registered image was ', im2dtype)
            
            # Convert the numpy data array from float to im2dtype:
            imRegNdaConverted = ConvertNumpyDtype(imRegNda,\
                                                  convertTo=im2dtype,\
                                                  debug=debug)
                
            
            # Create a new DICOM to modify with the registered pixel array:
            dicomReg = copy.deepcopy(dicomDict1_to_2[key]['Closest Dicom'])
            
            # Change the pixel array to imRegNdaConverted:
            #dicomReg.pixel_array = imRegNda # <- AttributeError: can't set attribute
            #dicomReg['PixelArray'] = imRegNda # <-TypeError: Dataset contents must be DataElement instances
            #dicomReg.pixel_array = preRegPixelArray # <- AttributeError: can't set attribute
            #dicomReg.pixel_array = dicomReg.pixel_array # <- same error
            #dicomReg.PixelData = imRegNda.tostring()
            dicomReg.PixelData = imRegNdaConverted.tostring()
            
            """
            Since the array dimensions have changed, update the relevant 
            metadata
            """
            #dicomReg.Rows, dicomReg.Columns = imRegNda.shape
            dicomReg.Rows, dicomReg.Columns = imRegNdaConverted.shape
            
            # Get the original (pre-registered) SOP Instance UID:
            #preRegSOPUID = dicomDict1_to_2[key]['Closest Dicom SOP UID']
            preRegSOPUID = dicomReg.SOPInstanceUID
            
            # Append to postRegSOPUIDs:
            preRegSOPUIDs.append(preRegSOPUID)
            
            
            # Generate a new SOP Instance UID:
            """
            Note that there may be slices that are repeated in the second set 
            due to differences in slice thickness, e.g. slice no 1 in the first 
            set may be closest to slices 1 and 2 in the second set.  Hence the
            same slice in the second set will have a unique SOP Instance UID
            if it is repeated.  In any case, each slice will need a unique SOP
            Instance UID to avoid throwing up errors when it comes to uploading
            DICOMs to XNAT.
            """
            postRegSOPUID = pydicom.uid.generate_uid() 
            
            # Append to postRegSOPUIDs:
            postRegSOPUIDs.append(postRegSOPUID)
            
            print('\n   Changing the SOP Instance UID from:\n', \
              preRegSOPUID, '\nto:\n', postRegSOPUID)
            

        
            # Make changes to the tags:
            dicomReg.SOPInstanceUID = postRegSOPUID
            dicomReg.SeriesDescription = postRegSeriesDesc
            dicomReg.SeriesNumber = postRegSeriesNo
            dicomReg.StudyInstanceUID = postRegStudyUID
            dicomReg.SeriesInstanceUID = postRegSeriesUID
            dicomReg.FrameOfReferenceUID = postRegFrameOfRefUID
            
            
            """
            # Get current filename:
            preRegFname = os.path.basename(fpath2)
            
            # Get filename (with extension):
            preRegFname = os.path.splitext(preRegFname)
            
            # Create a new filename, adding some info about the registration
            # method used:
            postRegFname = preRegFname[0] + '-' + regMethod + '-registered.dcm'
            
            print('\n   Changing the file name from:\n', \
              preRegFname[0], '\nto:\n', postRegFname)
            
            # Create new filepath:
            postRegFpath = os.path.join(toDicomRegDir, postRegFname)
            """
            
            
            # Export the DICOM file:
            print('\nExporting the registered DICOM to:\n', postRegFpath)
            dicomReg.save_as(postRegFpath)

        
            # Update the dictionary dicomDict1_to_2:
            dicomDict1_to_2[key].update(
                    {
                    'Registration method':regMethod,\
                    'Registered closest Dicom SOP UID':postRegSOPUID,\
                    'Registered closest Dicom fpath':postRegFpath,\
                    #'Registered closest Pix array':imRegNda,\
                    'Registered closest Pix array':imRegNdaConverted,\
                    'Registered closest Dicom':dicomReg
                    }
                    )
            
            # times[k] = this loop time:
            times.append(time.time())
            
            # Time taken to complete this loop:
            dt = times[-1] - times[-2]
            
            # Total time taken so far:
            Dt = times[-1] - times[0]
            
            # Average time taken per loop so far:
            period = Dt/(k+1)
            
            # Number of iterations remaining:
            rem = K - (k+1)
            
            # Anticipated time remaining to completion:
            t_rem = period*rem
            
            print(f'\n   Took {round(dt, 1)} s for this registration.')
            if Dt < 60:
                print(f'\n   Total time so far {round(Dt, 1)} s.')
            else:
                print(f'\n   Total time so far {round(Dt/60, 1)} min.')
            if t_rem < 60:
                print(f'\n   Expected time to completion {round(t_rem, 1)} s.')
            else:
                print(f'\n   Expected time to completion {round(t_rem/60, 1)} min.')
            print('')
        
        
            
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
        
        return dicomDict1_to_2, xnatDict 

