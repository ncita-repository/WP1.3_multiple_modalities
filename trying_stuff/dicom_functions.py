# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 07:23:50 2019

@author: ctorti


****************************************************************************
****************************************************************************
****************************************************************************


Function:
    dc_filepaths()
    
Purpose:
    To get list of full file paths for all DICOM files in a given directory.

Input:
    dirPath - path to directory containing DICOM files
    
Returns:
    filePaths - full file paths of all DICOM files in dirPath sorted in a 
                natural way (e.g. 3-1, 3-2, ..., 3-10, 3-11, ... rather than
                3-1, 3-10, 3-11, ... , 3-2, 3-20, ...)
"""

def dc_filepaths(dirPath):
    # Import packages:
    import os
    import natsort as ns
    
    filePaths = []  # create an empty list of all DICOM file paths

    # Use os to get list of file paths of DICOM-only files:
    for dirName, subdirList, fileList in os.walk(dirPath):
        for fileName in fileList:
            if '.dcm' in fileName.lower():  # check for DICOM files only
                filePaths.append(os.path.join(dirName,fileName))
                 
    # Sort files in natural way:
    filePaths = ns.natsorted(filePaths)
    
    return filePaths


"""

****************************************************************************
****************************************************************************
****************************************************************************



Function:
    dc_to_np()
    
Purpose:
    Load DICOM files from a list of file paths and store in numpy array.

Input:
    filePaths - full file paths of DICOM files to be loaded
    
Returns:
    x - array of x-coordinates for frame
    y - array of y-coordinates for frame
    z - array of z-coordinates for frame
    frames - DICOM files stored in 3D numpy array (frames along third dim),
             converted to uint8 dtype
"""


def dc_to_np(filePaths):
    # Load packages:
    import pydicom as dc
    import numpy as np
    
    # Load in first frame to obtain parameters required to create numpy array:
    frame = dc.read_file(filePaths[0])
    
    # Need to create an empty array so need dims of frames:

    # The 3D array will have dims Columns x Rows x Slices (X x Y x Z):
    dimX = frame.Columns
    dimY = frame.Rows
    dimZ = len(filePaths)
    
    # To display the images in the correct physiological scalings, 
    # need to know the pixel spacing (dpix) along x, y and z:
    dpixY = float(frame.PixelSpacing[0])
    dpixX = float(frame.PixelSpacing[1])
    dpixZ = float(frame.SliceThickness)
    
    # Create x, y and z axes (useful for plotting):
    x = np.arange(0.0, dimX*dpixX, dpixX)
    y = np.arange(0.0, dimY*dpixY, dpixY)
    z = np.arange(0.0, dimZ*dpixZ, dpixZ)
    
    # Create empty 3D array using same datatype as pixel_array (conversion
    # to uint8 will be done later):
    frames = np.zeros((dimY, dimX, dimZ), dtype=frame.pixel_array.dtype)
        
    # Load image data into data array by looping through each slice:
    for f in range(dimZ):
        frame = dc.read_file(filePaths[f])
        
        frames[:, :, f] = frame.pixel_array # store in array
    
    # Flip frames left-right and up-down to match orientation shown in
    # assignment instructions:    
    frames = np.fliplr(np.flipud(frames))    
    
    # Now can determine max pix value of array:
    maxPixVal = np.max(frames[:,:,:])
    
    # Use maxPixVal to convert from uint16 to uint8:
    frames = np.uint8(255.0*frames/maxPixVal)
    
    return x, y, z, frames



"""

****************************************************************************
****************************************************************************
****************************************************************************


Function:
    left_fill()
    
Purpose:
    Due to low pix values in body regions it's difficult to isolate background
    from body+breast without ignoring darker body regions. Starting with a 
    black & white mask, "left-filling" will involve looping along rows of 
    the B&W mask to determine the "centre-of-mass" of all white pixels along
    that row, then looping along columns, turn a black pixel white if it is 
    located left of the centre-of-mass.
    
Note:
    I accept that highly tailoring an algorithm to a particular dataset will
    mean that it does not generalise well to other datasets. I reasoned that 
    this would be an acceptable way of dealing with the challenge given 1) the
    time constraints; and 2) Since the task was to segment the data provided, 
    with no explicit instruction to ensure that the algorithm would generalise 
    to other data sets.
    

Input:
    mask - black and white mask with 0s and 1s
    
Returns:
    mask - the input mask left-filled as described above

"""




def left_fill(mask):
    # Load packages:
    import numpy as np
    
    # Get dims of array:
    dimR, dimC = mask.shape
    
    # Get maxPixVal of mask (this will allow what follows to work equally
    # with 0/1 or 0/255 values:
    maxPixVal = np.max(mask[:,:])
    
    for r in range(dimR):
        # Find the indexes of the white pixels for this row:
        inds = np.where(mask[r,:]==maxPixVal)[0] # this will return [] if 
        # no white pixels!

        if inds.size > 0: # if inds not []
            # Find the centre-of-mass:
            com = int(np.average([inds[0], inds[-1]]))

            for c in range(dimC):
                if c < com: # the pixel is left of com
                    mask[r,c] = 1

    return mask




"""

****************************************************************************
****************************************************************************
****************************************************************************


Function:
    segment_masks()
    
Purpose:
    To segment images in data set into three regions: body, breast and 
    background. Returns binary masks for each region.
    

Input:
    images - Images stored in a 3D numpy array
    
Returns:
    bodyMasks, breastMasks, backgroundMasks:
    
    bodyMasks - A 3D numpy array of binary masks for body regions
    breastMasks - A 3D numpy array of binary masks for breast regions
    backgroundMasks - A 3D numpy array of binary masks for background regions
    
    


Details:
    1. Body+breast is isolated from background by thresholding the image with
    "thresh_background".  A B&W mask for body+breast results.
    
    2. The mask is eroded to help eliminate speckles
    
    3. The mask is dilated to help restore mask to physiological size.
    
    4. The mask is "left-filled" as described in comments for function
    left_fill().
    
    5. A closing operation is performed to remove holes from mask.
    
    6. Resulting mask is "bodyAndBreastMask".
    
    7. From "bodyAndBreastMask" the mask for background ("backgroundMask") is
    easily determined.
    
    8. In order to isolate body from breast the images are blurred.  This
    reduces the challenges that arise from vasculature in the breast.
    
    9. The image is thresholded with "thresh_body". If too high, low intensity
    parts of the body get excluded. If too low, significant regions of breast
    with low pixel values remain. So despite a compromising level, there will
    remain a boundary of pixels surrounding the breast region that needs to
    be eroded...but first...
    
    10. Since background is included with the above thresholding, the 
    background can be eliminated using bodyAndBreastMask.
    
    11. The mask is eroded to try to eliminate low level pixels that surround
    the breast region.  Quite aggressive erosion is required to achieve this.
    
    12. Dilation is performed on the mask to try to restore the regions of the
    body that were eroded from the previous step.
    
    13. A closing operation is performed on the mask to close holes.
    
    14. What remains is a good guess at the body regions. However due to the 
    aggressive erosion and resulting high dilation in the latest steps, 
    what remain are regions that extend beyond the body+breast mask
    bodyAndBreastMask.  A simple correction is to set to 0 any pixel value
    of 1 in breastMask that does not have a pixel value of 1 in 
    bodyAndBreastMask.
    
    15. Since bodyAndBreastMask and bodyMask exist, breastMask can be obtained
    using the XOR of bodyAndBreastMask and bodyMask.
    
    16. The result for breastMask can be improved by dilating it. This helps
    to recover breast regions that were lost in the process of trying to 
    remove the body beneath breast.
    
    17. The latest dilation not only returned some pixels to breastMask that 
    belonged to the breast where it meets the body, it also increased the
    mask beyond where there is any breast (i.e. into the background regions).
    So once again use bodyAndBreastMask to correct breastMask by setting to 0
    any pixel with value 1 in breastMask that isn't also 1 in 
    bodyAndBreastMask.
    
    
"""


def segment_masks(images):
    # Load packages:
    import numpy as np
    import cv2
    import dicom_functions as dcf
    
    """
    *****************************************************
    Parameters for isolating body+breast from background:
    *****************************************************
    """
    
    
    """
    When thresholding the images to isolate body+breast from background 
    it seems reasonable to assign a threshold value as a % of the average
    pixel values averaged across all frames.
    """
    
    # The average pixel values averaged across all frames:
    avePixVal = np.average(images[:,:,:]) 
    
    # The threshold value for isolating body+breast from background,
    # expressed as a % of avePixVal:
    thresh_background = 0.4*avePixVal # i.e. 40% of avePixVal 
    
    
    
    """
    *****************************************************
    Parameters for isolating breast from body+breast:
    *****************************************************
    """
    
    """ 
    When thresholding the images to isolate breast from body+breast
    it seems reasonable to assign a threshold value as a % of the maximum
    pixel values averaged across all frames.
    """
    
    # Calculate maxPixVal across all frames:
    maxPixVal = np.max(images[:,:,:])
    
    # The threshold value for isolating the body from body+breast,
    # expressed as a % of maxPixVal:
    thresh_body = 0.2*maxPixVal # i.e. 20% of maxPixVal
    
    
    # Blur kernel for blurring of frames:
    blurKernelSize = 15
    
    
    """
    ***********************************************************
    Create kernels used for erosion, dilating and closing operations:
    ***********************************************************
    """
    Kernel3 = np.ones((3,3),np.uint8) # used for eroding and dilating
    Kernel5 = np.ones((5,5),np.uint8) # used for dilating
    Kernel25 = np.ones((25,25),np.uint8) # used for closing holes
    Kernel50 = np.ones((50,50),np.uint8) # used for closing holes
    
    
    # Get dims of frames for creation of empty arrays:
    dimY, dimX, dimZ = images.shape
    
    """
    ***********************************************************
    Create empty 3D arrays for bodyAndBreastMasks, backgroundMasks,
    bodyMasks and breastMasks:
    ***********************************************************
    """
    bodyAndBreastMasks = np.zeros((dimY, dimX, dimZ), dtype='uint8')
    backgroundMasks = np.zeros((dimY, dimX, dimZ), dtype='uint8')
    bodyMasks = np.zeros((dimY, dimX, dimZ), dtype='uint8')
    breastMasks = np.zeros((dimY, dimX, dimZ), dtype='uint8')

    
    """
    ***********************************************************
    Loop for each image in images:
    ***********************************************************
    """
    
    for i in range(dimZ):
        #print('Working on frame', str(f+1), 'of', str(dimZ), '...')
        
        image = images[:,:,i]
        
        """ 1) Isolate the body and breast using thresh_background: """
    
        # Create empty array for mask:
        mask = np.zeros((dimY, dimX), dtype='uint8')
        
        # Convert to B&W using thresh:
        for r in range(dimY):
            for c in range(dimX):
                if image[r,c] > thresh_background:
                    mask[r,c] = 1
        
            
        """ 2) Erode the mask: """

        mask = cv2.erode(mask, Kernel3, iterations = 3)
    
        
        """ 3) Dilate the mask: """    
        
        mask = cv2.dilate(mask, Kernel5, iterations = 3)
    
         
        """ 4) Left-fill the mask: """
                    
        # Left-fill holes by turning black pixels white if they are
        # left of the right-most white pixel at that same row (y):
        mask = dcf.left_fill(mask)
    
    
        """ 5) Close holes: """
                    
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, Kernel25)
    
        
        """ 6) Define bodyAndBreastMask: """
                    
        # Convert mask to binary:
        bodyAndBreastMask = cv2.threshold(mask, 0.5, 1, cv2.THRESH_BINARY)[1]
        
        
        """ 7) Define backgroundMask: """
        
        # Background mask is the logical compliment of bodyAndBreastMask:
        backgroundMask = np.logical_not(bodyAndBreastMask)
    
    
        # Store masks in 3D arrays:
        bodyAndBreastMasks[:,:,i] = bodyAndBreastMask
        backgroundMasks[:,:,i] = backgroundMask
    
        
        """ 8) Isolate body from breast in bodyAndBreastMask. 
               First blur image twice: """
      
        image = cv2.GaussianBlur(image, (blurKernelSize, blurKernelSize), 0)
        image = cv2.GaussianBlur(image, (blurKernelSize, blurKernelSize), 0) 
        image = cv2.GaussianBlur(image, (blurKernelSize, blurKernelSize), 0)
        
    
        """ 9) Use "thresh_body" to threshold blurred image.
               Note that background will remain for the time being. """
        
        # Create empty array for mask:
        mask = np.zeros((dimY, dimX), dtype='uint8')
        
        # Convert to B&W using thresh:
        for r in range(dimY):
            for c in range(dimX):
                if image[r,c] < thresh_body:
                    mask[r,c] = 1
        
        
        """ 10) Use bodyAndBreastMask to remove the background: """
        
        mask = mask*bodyAndBreastMask
        
        
        """ 11) Erode the mask: """    
        
        mask = cv2.erode(mask, Kernel5, iterations = 10)

        
        """ 12) Dilate the mask: """
        
        mask = cv2.dilate(mask, Kernel5, iterations = 35) # gets rid of body beneath breast
        
        
        """ 13) Close holes: """
            
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, Kernel50)
    
        
        """ 14) Correct mask by setting to black any white pixel that is not
        also white in bodyAndBreastMask: """
        
        for r in range(dimY):
            for c in range(dimX):
                if bodyAndBreastMask[r,c]==0:
                    mask[r,c] = 0
        
        
        # Convert mask to binary:
        bodyMask = cv2.threshold(mask, 0.5, 1, cv2.THRESH_BINARY)[1]
        
        # Store masks in 3D arrays:
        bodyMasks[:,:,i] = bodyMask
    
        
        """ 15) Define breastMask using bodyAndBreastMask and bodyMask: """    
            
        breastMask = np.bitwise_xor(bodyAndBreastMask, bodyMask)
        
        
        """ 16) Dilate breastMask to recover breast regions lost in attempts
        to remove body beneath breast: """    
                    
        breastMask = cv2.dilate(breastMask, Kernel5, iterations = 10)
        
        
        """ 17) Use bodyAndBreastMask to correct breastMask: """ 

        # Set to 0 any 1 in breastMask that isn't also 1 in bodyAndBreastMask:
        for r in range(dimY):
            for c in range(dimX):
                if bodyAndBreastMask[r,c]==0:
                    breastMask[r,c] = 0
        
        
        # Store final result for breastMask to 3D array:
        breastMasks[:,:,i] = breastMask
    
    return bodyMasks, breastMasks, backgroundMasks
