# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:09:22 2020

@author: ctorti
"""

def GetInputPoints(DicomDir, RoiFpath):
    # Import packages:
    #import os
    import SimpleITK as sitk
    #import pydicom
    #from GetDicoms import GetDicoms
    import importlib
    from PCStoICS import PCStoICS
    #import RegUtilityFuncs as ruf
    import GetContourPtsForDicoms
    importlib.reload(GetContourPtsForDicoms)
    from GetContourPtsForDicoms import GetContourPtsForDicoms
    
    
    # Get the contour points:    
    Pts_PCS, PtsArr_PCS, PtsDict_PCS = GetContourPtsForDicoms(DicomDir, 
                                                              RoiFpath)
    
    # Read in the 3D stack of DICOMs:
    SitkReader = sitk.ImageSeriesReader()
    SitkNames = SitkReader.GetGDCMSeriesFileNames(DicomDir)
    SitkReader.SetFileNames(SitkNames)
    SitkIm = SitkReader.Execute()
    SitkIm = sitk.Cast(SitkIm, sitk.sitkFloat32)

    # Convert points in Pts_PCS from PCS to ICS:
    Pts_ICS = PCStoICS(Pts_PCS=Pts_PCS, SitkIm=SitkIm)
    
    # Convert points in PtsArr_PCS from PCS to ICS:
    PtsArr_ICS = []
    
    for array in PtsArr_PCS:
        # If array is not empty:
        if array:
            PtsArr_ICS.append(PCStoICS(Pts_PCS=array, SitkIm=SitkIm))
        else:
            PtsArr_ICS.append([])
    
    
    # Convert points in Pts_ICS to a dictionary:
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in Pts_ICS and append the x, y 
    # and z coordinates to X, Y and Z:
    for point in Pts_ICS:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    
    PtsDict_ICS = {'x':x, 'y':y, 'z':z} 
    
    return Pts_PCS, PtsArr_PCS, PtsDict_PCS,\
           Pts_ICS, PtsArr_ICS, PtsDict_ICS