# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:32:41 2020

@author: ctorti
"""

"""
Function:
    CreateInputFileForElastix()
    
Purpose:
    Create a text file containing points for Elastix's Transformix to deform.
    
    A text file will be generated in the current working directory with the 
    format required of Elastix, i.e.:
        
    point
    number_of_points
    point1_x point1_y point1_z
    point2_x point2_y point2_z
    ...

    The contour points (ContourPts) will be obtained by running a function 
    called GetContourPoints3D(), that uses tags in the DICOMs (Dicoms) and the 
    DICOM-RTSTRUCT (ROI Collection generated by OHIF-Viewer) file to obtain the 
    contour points in image space.  Another function called GetDicoms() reads 
    in the DICOMs.
    
    FixedContourPts are in the form:
        
    [
     [],
     [],
     ...
     [
      [point1_x_contourN, point1_y_contourN, point1_z_contourN],
      [point2_x_contourN, point2_y_contourN, point2_z_contourN],
      ...
      [pointM_x_contourN, pointM_y_contourN, pointM_z_contourN]
      ], 
     [],
     ]
     
    where each array corresponds to a slice in Dicoms. 
    
    If any given array is empty, i.e. [], there is no contour for that slice.  
    If there is a contour, the points 
    
    [point1_x_contourN, point1_y_contourN, point1_z_contourN],
    [point2_x_contourN, point2_y_contourN, point2_z_contourN],
    ...
    [pointM_x_contourN, pointM_y_contourN, pointM_z_contourN]
      
    are the M points that make up a contour for the nth slice in Dicoms.


Input:
    DicomDir        - Directory containing DICOMs 
    
    RoiFpath        - Filepath of the DICOM-RTSTRUCT ROI Collection
    
    Origin          - The ImagePositionPatient (pydicom) or Origin (sitk)
                
                     e.g. [x0, y0, z0]
    
    Directions     - The direction cosine along x (rows), y (columns) and z 
                     (slices)
                 
                     e.g. [Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz]
                
                     SimpleITK returns all returns all direction cosines.
                     For pydicom, the cross product of the x and y vectors are
                     used to obtain the z vector.
    
    Spacings        - The pixel spacings along x, y and z 
    
                     e.g. [di, dj, dk]
                
                     SimpleITK returns all pixel spacings.
                     For pydicom, the SliceThickness is added to the PixelSpacing.
    
    CoordSys        - The coordinate system of the contour points returned;
                      acceptable values are:
                       - 'PCS' for Patient Coordinate System (i.e. ContourData
                         in Roi is in PCS)
                       - 'ICS' for Image Coordinate System (accounts for the
                          Image Plane attributes)
    
Returns:
    inputpoints.pts - Text file containing points in ROI Collection in format
                      required of Elastix
    
"""


def CreateInputFileForElastix(DicomDir, RoiFpath, 
                              Origin, Directions, Spacings, CoordSys):
    # Import packages:
    import os
    import importlib
    import GetInputPoints
    importlib.reload(GetInputPoints)
    from GetInputPoints import GetInputPoints
    
    
    # Get the contour points:    
    Pts_PCS, PtsArr_PCS, PtsDict_PCS,\
    Pts_ICS, PtsArr_ICS, PtsDict_ICS = GetInputPoints(DicomDir=DicomDir, 
                                                      RoiFpath=RoiFpath,
                                                      Origin=Origin,
                                                      Directions=Directions,
                                                      Spacings=Spacings)
    
    """ 
    Should use the contour points in the Patient Coordinate System since it
    seems that Elastix already performs a conversion to the Image Coordinate 
    System.
    """
    
    if CoordSys=='PCS':
        Points = Pts_PCS
    
    if CoordSys=='ICS':
        Points = Pts_ICS
        

    # Define the filename for the exported file:
    TextFname = 'inputpoints.txt'
    #TextFname = 'inputpoints.pcs'
    
    # Open a text file:
    TextFile = open(TextFname, 'w')

    # Write 'point' to the first line:
    TextFile.write('point')
    
    # Write the number of points to the second line:
    TextFile.write(f'\n{len(Points)}')
    
    # Write all points in ContourPts to each subsequent line in the required 
    # format:
    for Point in Points:
        TextFile.write(f'\n{Point[0]} {Point[1]} {Point[2]}')
                
    TextFile.close()
    
    # Get working directory:
    CWD = os.getcwd()
    
    print('\nContour points exported to', os.path.join(CWD, TextFname), '.')
    
    #return Pts_PCS, PtsArr_PCS, PtsDict_PCS,\
    #       Pts_ICS, PtsArr_ICS, PtsDict_ICS
    return