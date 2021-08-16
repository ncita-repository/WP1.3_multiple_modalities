# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:00:31 2020

@author: ctorti
"""



#def GetContourPtsForDicoms(DicomDir, RoiFpath, CoordSys):
def GetContourPtsForDicoms(DicomDir, RoiFpath):
    # Import packages:
    from GetDicomFpaths import GetDicomFpaths
    import importlib
    import GetContourPtsForDicom
    importlib.reload(GetContourPtsForDicom)
    from GetContourPtsForDicom import GetContourPtsForDicom
    #import GetContourPtsForDicom_v2
    #importlib.reload(GetContourPtsForDicom_v2)
    #from GetContourPtsForDicom_v2 import GetContourPtsForDicom_v2
    
    # Get the DICOM filepaths:
    DicomFpaths = GetDicomFpaths(DirPath=DicomDir, 
                                 SortMethod='slices', 
                                 Debug=False)

    
    # Initialise the array of contour points (AllContourPts) and array of array 
    # of contour points (AllContourPtsArr):
    AllPts_PCS = []
    AllPtsArr_PCS = []
    AllPts_ICS = []
    AllPtsArr_ICS = []

    
    # Loop through each DICOM filepath and get the array of contour points 
    # (ContourPts), array of array of contour points (ContourPtsArr) and 
    # dictionary of contour points (ContourPtsDict) that relate to it:
    for DicomFpath in DicomFpaths:
        #ContourPts, ContourPtsArr, ContourPtsDict = GetContourPtsForDicom(DicomFpath=DicomFpath,
        #                                                                  RoiFpath=RoiFpath,
        #                                                                  CoordSys=CoordSys)
        
        Pts_PCS, PtsArr_PCS, PtsDict_PCS, Pts_ICS, PtsArr_ICS, PtsDict_ICS = GetContourPtsForDicom(DicomFpath=DicomFpath,
                                                                                                   RoiFpath=RoiFpath)
        
        #Pts_PCS, PtsArr_PCS, PtsDict_PCS, Pts_ICS, PtsArr_ICS, PtsDict_ICS = GetContourPtsForDicom_v2(DicomFpath=DicomFpath,
        #                                                                                              RoiFpath=RoiFpath)
        
        
        AllPts_PCS.extend(Pts_PCS)
        AllPtsArr_PCS.append(PtsArr_PCS)
        
        AllPts_ICS.extend(Pts_ICS)
        AllPtsArr_ICS.append(PtsArr_ICS)
        
    
        
    # Create a dictionary:
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in AllContourPts and append the x, y 
    # and z coordinates to X, Y and Z:
    for Point in AllPts_PCS:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    AllPtsDict_PCS = {'x':x, 'y':y, 'z':z}      
    
    # Initialise an array to store the x, y and z coordinates:
    x = []
    y = []
    z = []  
    
    # Loop through each array of points in AllContourPts and append the x, y 
    # and z coordinates to X, Y and Z:
    for Point in AllPts_ICS:
        x.append(Point[0])
        y.append(Point[1])
        z.append(Point[2])
    
    AllPtsDict_ICS = {'x':x, 'y':y, 'z':z} 
        
        
    #return AllContourPts, AllContourPtsArr, AllContourPtsDict
    return AllPts_PCS, AllPtsArr_PCS, AllPtsDict_PCS, AllPts_ICS, AllPtsArr_ICS, AllPtsDict_ICS
        