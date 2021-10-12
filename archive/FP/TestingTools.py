# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:19:21 2021

@author: ctorti
"""



"""
******************************************************************************
******************************************************************************
EXAMPLE CODE FOR TESTING/RUNNING OF ROI COPYING FUNCTIONS
******************************************************************************
******************************************************************************
"""  





"""
******************************************************************************
Re-generate config files
******************************************************************************
"""
import importlib
import ConfigTools
importlib.reload(ConfigTools)
from ConfigTools import CreateCopyRoiConfigFile, CreateTestConfigFile

CreateCopyRoiConfigFile()

CreateTestConfigFile()






"""
******************************************************************************
Establish connection to XNAT
******************************************************************************
""" 
from XnatTools import CreateXnatSession

XnatSession = CreateXnatSession('http://10.1.1.20')

print(XnatSession)




"""
******************************************************************************
Test RunCopyXnatRoi() with pre-configured Test inputs
******************************************************************************
""" 
import importlib
import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import RunCopyXnatRoiCol
#from pypref import Preferences
#import os

TestNum = 'RR1'
#TestNum = 'RR2'
#TestNum = 'RR3'
TestNum = 'RR6'
#TestNum = 'RD1'
#TestNum = 'RD3'
#TestNum = 'RD5'
#TestNum = 'SR1'
#TestNum = 'SR2'
#TestNum = 'SR3'
#TestNum = 'SR4'
#TestNum = 'SR5'
#TestNum = 'RD5'
TestNum = 'RR15'

LogToConsole = False
PlotResults = True
ExportPlot = True
DevOutputs = True

""" Minimum inputs (XnatSession, TestNum): """
#RunCopyXnatRoi(XnatSession, TestNum)

""" Additional inputs (LogToConsole, PlotResults, ExportPlot): """
#RunCopyXnatRoiCol(XnatSession, TestNum, LogToConsole, PlotResults, ExportPlot)

""" Additional inputs and returning of outputs (i.e. DevOutputs = True): """
XnatSession, SrcRoiCol, TrgRoiCol, NewTrgRoiCol, Dro,\
PathsDict, DictOfInputs, ListOfInputs, TimingMsgs, Times\
= RunCopyXnatRoiCol(XnatSession, TestNum, LogToConsole, PlotResults, ExportPlot, DevOutputs)






"""
******************************************************************************
Test CopyXnatRoi() with pre-configured Test inputs
******************************************************************************
"""  
import importlib
import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopyXnatRoiCol
from pypref import Preferences
import os

TestNum = 'RR1'
TestNum = 'RR2'
TestNum = 'RR3'
TestNum = 'RD1'
TestNum = 'RD3'
TestNum = 'RD5'
TestNum = 'SR1'
TestNum = 'SR2'
TestNum = 'SR3'
TestNum = 'SR4'
TestNum = 'SR5'
#TestNum = 'RD5'


#XnatSession = None
PathsDict = None
#XnatDownloadDir = 'default'
XnatDownloadDir = r'C:\Users\ctorti\Downloads\XnatDownloads\test3' 

# Get the inputs for this TestNum:
TestCfg = Preferences(directory=os.getcwd(), filename='CopyRoiTestConfig.py')

XnatUrl = TestCfg.get('XnatUrl')
ProjId = TestCfg.get('ProjId')
SubjLabel = TestCfg.get('SubjLabel')
SrcExpLabel = TestCfg.get(TestNum)['SrcExpLabel']
SrcScanId = TestCfg.get(TestNum)['SrcScanId']
SrcSliceNum = TestCfg.get(TestNum)['SrcSliceNum']
SrcRoiColMod = TestCfg.get(TestNum)['SrcRoiColMod']
SrcRoiColName = TestCfg.get(TestNum)['SrcRoiColName']
SrcRoiName = TestCfg.get(TestNum)['SrcRoiName']
TrgExpLabel = TestCfg.get(TestNum)['TrgExpLabel']
TrgScanId = TestCfg.get(TestNum)['TrgScanId']
TrgSliceNum = TestCfg.get(TestNum)['TrgSliceNum']
TrgRoiColMod = TestCfg.get(TestNum)['TrgRoiColMod']
TrgRoiColName = TestCfg.get(TestNum)['TrgRoiColName']
TrgRoiName = TestCfg.get(TestNum)['TrgRoiName']
TxtToAddToTrgRoiColName= TestCfg.get(TestNum)['TxtToAddToTrgRoiColName']

# Get other default inputs:
Defaults = Preferences(directory=os.getcwd(), filename='CopyRoiDefaults.py')

ResInterp = Defaults.get('ResInterp')
PreResVar = Defaults.get('PreResVariance')
ApplyPostResBlur = Defaults.get('ApplyPostResBlur')
PostResVar = Defaults.get('PostResVariance')
ForceReg = Defaults.get('ForceReg')
TxMatrix = Defaults.get('TxMatrix')
SelxOrSitk = Defaults.get('SelxOrSitk')
Transform = Defaults.get('Transform')
MaxIters = Defaults.get('MaxIters')
TxInterp = Defaults.get('TxInterp')
ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
PostTxVar = Defaults.get('PostTxVariance')
ExportNewTrgRoiCol = Defaults.get('ExportNewTrgRoiCol')
RtsExportDir = Defaults.get('RtsExportDir')
SegExportDir = Defaults.get('SegExportDir')
ExportNewDro = Defaults.get('ExportNewDro')
DroExportDir = Defaults.get('DroExportDir')
UploadNewDro = Defaults.get('UploadNewDro')
PlotAllSlices = Defaults.get('PlotAllSlices')
RtsPlotExportDir = Defaults.get('RtsPlotExportDir')
SegPlotExportDir = Defaults.get('SegPlotExportDir')
LogToConsole = Defaults.get('LogToConsole')
ExportLogFiles = Defaults.get('ExportLogFiles')
LogExportDir = Defaults.get('LogExportDir')

print(f'\nPre-configured settings for TestNum {TestNum}:')
print('*****************************************')
print(f'ProjId = {ProjId}')
print(f'SubjLabel = {SubjLabel}')
print(f'SrcExpLabel = {SrcExpLabel}')
print(f'SrcScanId = {SrcScanId}')
print(f'SrcSliceNum = {SrcSliceNum}')
print(f'SrcRoiColMod = {SrcRoiColMod}')
print(f'SrcRoiColName = {SrcRoiColName}')
print(f'SrcRoiName = {SrcRoiName}')
print(f'TrgExpLabel = {TrgExpLabel}')
print(f'TrgScanId = {TrgScanId}')
print(f'TrgSliceNum = {TrgSliceNum}')
print(f'TrgRoiColMod = {TrgRoiColMod}')
print(f'TrgRoiColName = {TrgRoiColName}')
print(f'TrgRoiName = {TrgRoiName}')

print('\nDefault settings for ROI copying algorithm:')
print('*******************************************')
print(f'ResInterp = {ResInterp}')
print(f'PreResVar = {PreResVar}')
print(f'ApplyPostResBlur = {ApplyPostResBlur}')
print(f'PostResVar = {PostResVar}')
print(f'ForceReg = {ForceReg}')
print(f'TxMatrix = {TxMatrix}')
print(f'SelxOrSitk = {SelxOrSitk}')
print(f'Transform = {Transform}')
print(f'MaxIters = {MaxIters}')
print(f'TxInterp = {TxInterp}')
print(f'ApplyPostTxBin = {ApplyPostTxBin}')
print(f'ApplyPostTxBlur = {ApplyPostTxBlur}')
print(f'PostTxVar = {PostTxVar}')

print(f'ExportNewTrgRoiCol = {ExportNewTrgRoiCol}')
print(f'RtsExportDir = {RtsExportDir}')
print(f'SegExportDir = {SegExportDir}')
print(f'ExportNewDro = {ExportNewDro}')
print(f'DroExportDir = {DroExportDir}')
print(f'UploadNewDro = {UploadNewDro}')
print(f'PlotAllSlices = {PlotAllSlices}')
print(f'RtsPlotExportDir = {RtsPlotExportDir}')
print(f'SegPlotExportDir = {SegPlotExportDir}')
print(f'LogToConsole = {LogToConsole}')
print(f'ExportLogFiles = {ExportLogFiles}')
print(f'LogExportDir = {LogExportDir}')

TrgRoi, Dro, PathsDict, XnatSession, DictOfInputs, ListOfInputs,\
TimingMsgs = CopyXnatRoiCol(XnatUrl, XnatSession, ProjId, SubjLabel, 
                            SrcExpLabel, SrcScanId, SrcSliceNum, 
                            SrcRoiColMod, SrcRoiColName, SrcRoiName, 
                            TrgExpLabel, TrgScanId, TrgSliceNum, 
                            TrgRoiColMod, TrgRoiColName, TrgRoiName, 
                            TxtToAddToTrgRoiColName, 
                            PathsDict, XnatDownloadDir, LogToConsole, ExportLogFiles,
                            ResInterp, PreResVar, ApplyPostResBlur, ForceReg,  
                            TxMatrix, Transform, MaxIters, TxInterp, 
                            ApplyPostTxBin, ApplyPostTxBlur, PostTxVar)







"""
******************************************************************************
Test CopyXnatRoi() with explicit inputs
******************************************************************************
"""
import importlib
import RoiCopyTools
importlib.reload(RoiCopyTools)
from RoiCopyTools import CopyXnatRoi

XnatUrl = 'http://10.1.1.20'
ProjId = 'NCITA_TEST'
SubjLabel = 'TCGA-BB-A5HY'

SrcStudyLabel = 'Session2'
SrcSeriesLabel = '2'
SrcAsrMod = 'RTSTRUCT'
SrcAsrName = 'Left eye CT'
SrcRoiName = None
SrcSliceNum = None
TrgStudyLabel = 'Session2'
TrgSeriesLabel = '4'
TrgAsrMod = 'RTSTRUCT'
TrgAsrName = None
TrgRoiName = None
TrgSliceNum = None
TxtToAddToTrgAsrName = ' 4'

#if not 'XnatSession' in locals():
#    XnatSession = None
PathsDict = None
#XnatDownloadDir = 'default'
XnatDownloadDir = r'C:\Users\ctorti\Downloads\XnatDownloads\test3' 
LogToConsole = False
ExportLogFiles = False
ResInterp = 'BlurThenLinear'
PreResVar = (1,1,1)
ForceReg = False,   
Tx = 'affine'
TxMaxIters = '512'
TxInterp = 'NearestNeighbor'
ApplyPostTxBlur = True
PostTxVar = (2,2,2)
ApplyPostTxBin = True

#print(f'\nTrgSliceNum = {TrgSliceNum}')

TrgRoi, PathsDict, XnatSession, DictOfInputs, ListOfInputs,\
ListOfTimings = CopyXnatRoi(XnatUrl, XnatSession, ProjId, SubjLabel, SrcStudyLabel,  
                            SrcSeriesLabel, SrcAsrMod, SrcAsrName, SrcRoiName, SrcSliceNum,   
                            TrgStudyLabel, TrgSeriesLabel, TrgAsrMod, TrgAsrName, TrgRoiName,
                            TrgSliceNum, TxtToAddToTrgAsrName, 
                            PathsDict, XnatDownloadDir, LogToConsole, ExportLogFiles,
                            ResInterp, PreResVar, ForceReg, Tx, TxMaxIters, TxInterp, 
                            ApplyPostTxBlur, PostTxVar, ApplyPostTxBin)








"""
******************************************************************************
Run subset of code in CopyXnatRoi() with explicit inputs
******************************************************************************
"""
import time
import os
from pypref import Preferences
from XnatTools import CreateXnatSession, DownloadScan, DownloadImAsr
from XnatTools import GetFnameAndIdFromName
from ImageTools import ImportImage, RegisterImagesSelx, RegisterImagesSitk

TestNum = 'RR15'

TestCfg = Preferences(directory=os.getcwd(), filename='CopyRoiTestConfig.py')

XnatUrl = TestCfg.get('XnatUrl')
ProjId = TestCfg.get('ProjId')
SubjLabel = TestCfg.get('SubjLabel')
SrcExpLabel = TestCfg.get(TestNum)['SrcExpLabel']
SrcScanId = TestCfg.get(TestNum)['SrcScanId']
SrcSliceNum = TestCfg.get(TestNum)['SrcSliceNum']
SrcRoiColMod = TestCfg.get(TestNum)['SrcRoiColMod']
SrcRoiColName = TestCfg.get(TestNum)['SrcRoiColName']
SrcRoiName = TestCfg.get(TestNum)['SrcRoiName']
TrgExpLabel = TestCfg.get(TestNum)['TrgExpLabel']
TrgScanId = TestCfg.get(TestNum)['TrgScanId']
TrgSliceNum = TestCfg.get(TestNum)['TrgSliceNum']
TrgRoiColMod = TestCfg.get(TestNum)['TrgRoiColMod']
TrgRoiColName = TestCfg.get(TestNum)['TrgRoiColName']
TrgRoiName = TestCfg.get(TestNum)['TrgRoiName']
TxtToAddToTrgRoiColName= TestCfg.get(TestNum)['TxtToAddToTrgRoiColName']

# Get other default inputs:
Defaults = Preferences(directory=os.getcwd(), filename='CopyRoiDefaults.py')

ResInterp = Defaults.get('ResInterp')
PreResVar = Defaults.get('PreResVariance')
ApplyPostResBlur = Defaults.get('ApplyPostResBlur')
PostResVar = Defaults.get('PostResVariance')
ForceReg = Defaults.get('ForceReg')
TxMatrix = Defaults.get('TxMatrix')
SelxOrSitk = Defaults.get('SelxOrSitk')
Transform = Defaults.get('Transform')
MaxIters = Defaults.get('MaxIters')
TxInterp = Defaults.get('TxInterp')
ApplyPostTxBin = Defaults.get('ApplyPostTxBin')
ApplyPostTxBlur = Defaults.get('ApplyPostTxBlur')
PostTxVar = Defaults.get('PostTxVariance')
ExportNewTrgRoiCol = Defaults.get('ExportNewTrgRoiCol')
RtsExportDir = Defaults.get('RtsExportDir')
SegExportDir = Defaults.get('SegExportDir')
ExportNewDro = Defaults.get('ExportNewDro')
DroExportDir = Defaults.get('DroExportDir')
UploadNewDro = Defaults.get('UploadNewDro')
PlotAllSlices = Defaults.get('PlotAllSlices')
RtsPlotExportDir = Defaults.get('RtsPlotExportDir')
SegPlotExportDir = Defaults.get('SegPlotExportDir')
LogToConsole = Defaults.get('LogToConsole')
ExportLogFiles = Defaults.get('ExportLogFiles')
LogExportDir = Defaults.get('LogExportDir')

print(f'\nPre-configured settings for TestNum {TestNum}:')
print('*****************************************')
print(f'ProjId = {ProjId}')
print(f'SubjLabel = {SubjLabel}')
print(f'SrcExpLabel = {SrcExpLabel}')
print(f'SrcScanId = {SrcScanId}')
print(f'SrcSliceNum = {SrcSliceNum}')
print(f'SrcRoiColMod = {SrcRoiColMod}')
print(f'SrcRoiColName = {SrcRoiColName}')
print(f'SrcRoiName = {SrcRoiName}')
print(f'TrgExpLabel = {TrgExpLabel}')
print(f'TrgScanId = {TrgScanId}')
print(f'TrgSliceNum = {TrgSliceNum}')
print(f'TrgRoiColMod = {TrgRoiColMod}')
print(f'TrgRoiColName = {TrgRoiColName}')
print(f'TrgRoiName = {TrgRoiName}')

print('\nDefault settings for ROI copying algorithm:')
print('*******************************************')
print(f'ResInterp = {ResInterp}')
print(f'PreResVar = {PreResVar}')
print(f'ApplyPostResBlur = {ApplyPostResBlur}')
print(f'PostResVar = {PostResVar}')
print(f'ForceReg = {ForceReg}')
print(f'TxMatrix = {TxMatrix}')
print(f'SelxOrSitk = {SelxOrSitk}')
print(f'Transform = {Transform}')
print(f'MaxIters = {MaxIters}')
print(f'TxInterp = {TxInterp}')
print(f'ApplyPostTxBin = {ApplyPostTxBin}')
print(f'ApplyPostTxBlur = {ApplyPostTxBlur}')
print(f'PostTxVar = {PostTxVar}')

print(f'ExportNewTrgRoiCol = {ExportNewTrgRoiCol}')
print(f'RtsExportDir = {RtsExportDir}')
print(f'SegExportDir = {SegExportDir}')
print(f'ExportNewDro = {ExportNewDro}')
print(f'DroExportDir = {DroExportDir}')
print(f'UploadNewDro = {UploadNewDro}')
print(f'PlotAllSlices = {PlotAllSlices}')
print(f'RtsPlotExportDir = {RtsPlotExportDir}')
print(f'SegPlotExportDir = {SegPlotExportDir}')
print(f'LogToConsole = {LogToConsole}')
print(f'ExportLogFiles = {ExportLogFiles}')
print(f'LogExportDir = {LogExportDir}')

# Changes to default settings:
Transform = 'bspline'
ExportNewTrgRoiCol = False
ExportNewDro = False
UploadNewDro = False

print('\nChanges to default settings:')
print('****************************')
print(f'Transform = {Transform}')
print(f'ExportNewTrgRoiCol = {ExportNewTrgRoiCol}')
print(f'ExportNewDro = {ExportNewDro}')
print(f'UploadNewDro = {UploadNewDro}')

times = []

""" Create XNAT session. """
if XnatSession == None:
    XnatSession = CreateXnatSession(XnatUrl)


""" Start timing. """
times.append(time.time())

""" Download data from XNAT. """
PathsDict, XnatSession\
= DownloadScan(XnatUrl, ProjId, SubjLabel, SrcExpLabel, SrcScanId, 
               XnatSession, XnatDownloadDir, PathsDict)

PathsDict, XnatSession\
= DownloadScan(XnatUrl, ProjId, SubjLabel, TrgExpLabel, TrgScanId, 
               XnatSession, XnatDownloadDir, PathsDict)

PathsDict, XnatSession\
= DownloadImAsr(XnatUrl, ProjId, SubjLabel, SrcExpLabel, 
                SrcScanId, SrcRoiColMod, SrcRoiColName, 
                XnatSession, XnatDownloadDir, PathsDict)

if TrgRoiColName != None:
    PathsDict, XnatSession\
    = DownloadImAsr(XnatUrl, ProjId, SubjLabel, TrgExpLabel, 
                    TrgScanId, TrgRoiColMod, TrgRoiColName, 
                    XnatSession, XnatDownloadDir, PathsDict)

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
if True:#LogToConsole:
    msg = f'Took {Dtime} s to download data from XNAT.\n'
    #TimingMsgs.append(msg)
    print(f'*{msg}')


#print('\n', PathsDict, '\n')

""" Get the DICOM directories and ROI Collection file paths: """
SrcDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
            ['experiments'][SrcExpLabel]['scans'][SrcScanId]\
            ['resources']['DICOM']['files']['Dir']

TrgDcmDir = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
            ['experiments'][TrgExpLabel]['scans'][TrgScanId]\
            ['resources']['DICOM']['files']['Dir']

SrcRoiColFname, SrcAsrId\
= GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, SrcExpLabel,
                        SrcRoiColMod, SrcRoiColName)

SrcRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                 ['experiments'][SrcExpLabel]['assessors'][SrcAsrId]\
                 ['resources'][SrcRoiColMod]['files'][SrcRoiColFname]\
                 ['Fpath']

if TrgRoiColName == None:
    TrgRoiColFpath = None
else:
    TrgRoiColFname, TrgAsrId\
    = GetFnameAndIdFromName(PathsDict, ProjId, SubjLabel, TrgExpLabel,
                            TrgRoiColMod, TrgRoiColName)

    TrgRoiColFpath = PathsDict['projects'][ProjId]['subjects'][SubjLabel]\
                     ['experiments'][TrgExpLabel]['assessors'][TrgAsrId]\
                     ['resources'][TrgRoiColMod]['files'][TrgRoiColFname]\
                     ['Fpath']


""" Import the 3D images. """
SrcIm = ImportImage(SrcDcmDir)
TrgIm = ImportImage(TrgDcmDir)


times.append(time.time())
            
""" Register SrcIm to TrgIm. """
if SelxOrSitk == 'Selx':
    RegIm, SelxImFiltOrSitkTx, TxParamMapDict\
    = RegisterImagesSelx(FixIm=TrgIm, MovIm=SrcIm, 
                         Transform=Transform,
                         MaxNumOfIters=MaxIters, 
                         LogToConsole=LogToConsole)

    """ Get the registration transform parameters: """
    TxParams = list(TxParamMapDict['TransformParameters'])

if SelxOrSitk == 'Sitk':
    RegIm, SelxImFiltOrSitkTx\
    = RegisterImagesSitk(FixIm=TrgIm, MovIm=SrcIm, 
                         Transform=Transform, InitMethod='moments',
                         FinalInterp='Linear',
                         FixFidsFpath='fixed_fiducials.txt', 
                         MovFidsFpath='moving_fiducials.txt',
                         FlipK=False, LogToConsole=LogToConsole)

    """ Get the registration transform parameters: """
    TxParams = [str(item) for item in SelxImFiltOrSitkTx.GetParameters()]

times.append(time.time())
Dtime = round(times[-1] - times[-2], 1)
msg = f'Took {Dtime} s to register the Source image to the Target '\
      + f'image using {Transform} transform.\n'
#ListOfTimings.append(msg)
if LogToConsole == False:
    print(f'*{msg}')